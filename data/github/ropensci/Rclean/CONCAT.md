<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- # ijtiff  <img src="man/figures/logo.png" height="140" align="right"> -->
<!-- Code status -->
[![Build
Status](https://travis-ci.org/ROpenSci/Rclean.svg?branch=master)](https://travis-ci.org/ROpenSci/Rclean)
[![Coverage
status](https://codecov.io/gh/ROpenSci/Rclean/branch/master/graph/badge.svg)](https://codecov.io/github/ROpenSci/Rclean?branch=master)

<!-- R status -->
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/Rclean)](https://cran.r-project.org/package=Rclean)
![RStudio CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/Rclean)
![RStudio CRAN monthly
downloads](http://cranlogs.r-pkg.org/badges/Rclean)
[![Rdocumentation](http://www.rdocumentation.org/badges/version/Rclean)](http://www.rdocumentation.org/packages/Rclean)

<!-- Dev status -->
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

<!-- Package Review -->
[![](https://badges.ropensci.org/327_status.svg)](https://github.com/ropensci/software-review/issues/327)
[![status](http://joss.theoj.org/papers/334d80d5508056dc6e7e17c6fd3ed5a6/status.svg)](http://joss.theoj.org/papers/334d80d5508056dc6e7e17c6fd3ed5a6)

<!-- Archiving -->
[![DOI](https://zenodo.org/badge/102645585.svg)](https://zenodo.org/badge/latestdoi/102645585)

Quick Start Guide
=================

-   [Rclean](https://github.com/MKLau/Rclean) was created to help
    scientists more *easily* write “cleaner” code.
-   Written with research scientists that are results oriented in mind,
    the package’s primary function provides a simple way to isolate the
    minimal code you need to produce a specific result, such as a
    statistical table or a figure. By focusing on specific results (aka.
    variables), large and/or complicated analytical scripts can be
    paired down to the essentials and easily re-factored to be more
    robust and easily shared.
-   Below, you'll find a brief introduction to get you started using
    the package. For more details, see `vignette("Rclean")`.

Install
=======

You can install the most up to date version easily with
[devtools](https://github.com/hadley/devtools):

    install.packages("devtools")
    devtools::install_github("MKLau/Rclean")

You will also likely need to install the
[RGraphViz](bioconductor.org/packages/release/bioc/html/Rgraphviz.html):


    install.packages("BiocManager")
    BiocManager::install("Rgraphviz")

Once installed, per usual R practice just load the *Rclean* package
with:

    library(Rclean)

Usage
=====

*Rclean* usage is simple. Have a script with code you want to clean
saved to disk. Then, just run the `clean` function with the path to the
script as the input. Here, we can use an example script that is included
with the package:

    script <- system.file("example", "simple_script.R", package = "Rclean")

Here's a quick look at the code:

    readLines(script)
    #>  [1] "## Make a data frame"                             
    #>  [2] "mat <- matrix(rnorm(400), nrow = 100)"            
    #>  [3] "dat <- as.data.frame(mat)"                        
    #>  [4] "dat[, \"V2\"] <- dat[, \"V2\"] + runif(nrow(dat))"
    #>  [5] "dat[, \"V5\"] <- gl(10, 10)"                      
    #>  [6] "## Conduct some analyses"                         
    #>  [7] "fit12 <- lm(V1 ~ V2, data = dat)"                 
    #>  [8] "fit13 <- lm(V1 ~ V3, data = dat)"                 
    #>  [9] "fit14 <- lm(V1 ~ V4, data = dat)"                 
    #> [10] "fit15.aov <- aov(V1 ~ V2 + V5, data = dat)"       
    #> [11] "## Summarize analyses"                            
    #> [12] "summary(fit15.aov)"                               
    #> [13] "tab.12 <- summary(fit12)"                         
    #> [14] "tab.13 <- summary(fit13)"                         
    #> [15] "tab.14 <- summary(fit14)"                         
    #> [16] "tab.15 <- append(fit15.aov, tab.14)"              
    #> [17] "## Conduct a calculation"                         
    #> [18] "dat <- 25 + 2"                                    
    #> [19] "dat[2] <- 10"                                     
    #> [20] "out <- dat * 2"

You can get a list of the variables found in an object with `get_vars`.

    get_vars(script)
    #>  [1] "mat"       "dat"       "fit12"     "fit13"     "fit14"     "fit15.aov"
    #>  [7] "tab.12"    "tab.13"    "tab.14"    "tab.15"    "out"

Sometimes for more complicated scripts, it can be helpful to see a
network graph showing the interdependencies of variables. `code_graph`
will produce a network diagram showing which lines of code produce or
use which variables:


    code_graph(script)

<img src="man/figures/README-unnamed-chunk-7-1.png" width="75%" />

Now, we can pick the result we want to focus on for cleaning:


    clean(script, "tab.15")
    #> mat <- matrix(rnorm(400), nrow = 100)
    #> dat <- as.data.frame(mat)
    #> dat[, "V2"] <- dat[, "V2"] + runif(nrow(dat))
    #> dat[, "V5"] <- gl(10, 10)
    #> fit14 <- lm(V1 ~ V4, data = dat)
    #> fit15.aov <- aov(V1 ~ V2 + V5, data = dat)
    #> tab.14 <- summary(fit14)
    #> tab.15 <- append(fit15.aov, tab.14)
    #> dat <- 25 + 2
    #> dat[2] <- 10

We can also select several variables at the same time:

    my.vars <- c("tab.12", "tab.15")
    clean(script, my.vars)
    #> mat <- matrix(rnorm(400), nrow = 100)
    #> dat <- as.data.frame(mat)
    #> dat[, "V2"] <- dat[, "V2"] + runif(nrow(dat))
    #> dat[, "V5"] <- gl(10, 10)
    #> fit12 <- lm(V1 ~ V2, data = dat)
    #> fit14 <- lm(V1 ~ V4, data = dat)
    #> fit15.aov <- aov(V1 ~ V2 + V5, data = dat)
    #> tab.12 <- summary(fit12)
    #> tab.14 <- summary(fit14)
    #> tab.15 <- append(fit15.aov, tab.14)
    #> dat <- 25 + 2
    #> dat[2] <- 10

While just taking a look at the simplified code can be very helpful, you
can also save the code for later use or sharing (e.g. creating a
reproducible example for getting help) with `keep`:

    my.code <- clean(script, my.vars)
    keep(my.code, file = "results_tables.R")

If you would like to copy your code to your clipboard, you can do that
by not specifying a file path. You can now paste the simplified as
needed.

    keep(my.code)

Contributing
============

This is an open-source project. If you would like to contribute to the
project, please check out [CONTRIBUTING.md](CONTRIBUTING.md).

Please note that the 'Rclean' project is released with a [Contributor
Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project,
you agree to abide by its terms.

![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)
# Rclean News
---

## version 1.1.7


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
(https://www.contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/.
Want to contribute to Rclean?
=============================

Great, thank you!

We foster open, respectful software development following the code of
conduct put forward by the R Community. Please take a moment to review
and adhere to these guidelines when contributing:

[R Community Code of Conduct](https://wiki.r-consortium.org/view/R_Consortium_and_the_R_Community_Code_of_Conduct)

# Issues

Just using the package and giving feedback is contributing. If
the package is not working as described, congratulations you've
discovered a bug! If you think this is the case, please submit an
issue to the github
[issue](https://github.com/MKLau/Rclean/issues) system. 

Include a reproducible example in the issue or a link to a
[gist](https://gist.github.com/), with the following:

1. Succinct title and description of expected results 
2. Code with library dependencies and representative data for input
3. Observed results description

Please check that the example runs in a fresh R session before
submitting. 

# Enhancements

If you haven't run into a bug, but would like the package to do
something that it currently doesn't, you're welcome to:

1. *Fork the repository and develop your own additions.* Take a look
   at how our code is structured (including the comments and help
   documentation). Also, develop appropriate tests to check that your
   code is correctly integrated and runs (see the [[tests]] directory).
2. *Submit an enhancement to the issue system.* We welcome new
   ideas. If you think there's a feature that would be useful for the
   community of users to have available in the package, but aren't
   able to develop it your self, please let us know.

# Retrospective Provenance

If you have used or are a developer who uses retrospective provenance
that is not currently handled by *Rclean*, please submit an
*enhancement* to the issue system. We'll be happy to help implement
it. 

# Thanks for contributing and keeping scientific software open!
## Test environments
* local OS X install, R 3.6.1
* ubuntu 16.04 (on travis-ci), R unstable r77518
* ubuntu 16.04 (on travis-ci), R 3.6.1
* ubuntu 16.04 (on travis-ci), R 3.5.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Comments

* This release fixes a previous violation of CRAN policy for writing
  to the user file system by using the clipr package. This is in
  keep.R, which is for the explicitly stated use of writing objects to
  disk or to the clipboard.
---
name: Bug report
about: Something went wrong with Rclean
title: ''
labels: bug
assignees: MKLau

---

## Check Here First

- [ ] Read and abide by `Rclean`'s
  [code of conduct](https://https://github.com/MKLau/Rclean/blob/master/CODE_OF_CONDUCT.md).
- [ ] Search for duplicates among the
  [existing issues](https://github.com/MKLau/Rclean/issues), both open
  and closed.
- [ ] Advanced users: mention the
  [SHA-1 hash](https://git-scm.com/book/en/v1/Getting-Started-Git-Basics#Git-Has-Integrity)
  of the
  [Git commit you install](https://github.com/MKLau/Rclean/commits/master).

## Description

Describe the bug clearly and concisely. 

## Reproducible example

Provide a minimal reproducible example with code and output that
demonstrates the problem. The `reprex()` function from the
[`reprex`](https://github.com/tidyverse/reprex) package is extremely
helpful for this.

To help us read your code, please try to follow the
[tidyverse style guide](https://style.tidyverse.org/). The
`style_text()` and `style_file()` functions from the
[`styler`](https://github.com/r-lib/styler) package make it easier.

## Expected result

What should have happened? Please be as specific as possible.

## Session info

End the reproducible example with a call to `sessionInfo()` in the
same session (e.g. `reprex(si = TRUE)`) and include the output.
---
title: 'Rclean: A Tool for Writing Cleaner, More Transparent Code'
tags:
  - R
  - reproducibility
  - transparency
  - code cleaning
  - data provenance
authors:
  - name: Matthew K. Lau
    orcid: 0000-0003-3758-2406
    affiliation: 1
  - name: Thomas F. J.-M. Pasquier
    orcid: 0000-0001-6876-1306
    affiliation: "2, 3" 
  - name: Margo Seltzer
    orcid: 0000-0002-2165-4658
    affiliation: "4"
affiliations:
 - name: Harvard Forest, Harvard University 
   index: 1
 - name: Department of Computer Science, University of Bristol 
   index: 2
 - name: School of Engineering and Applied Science, Harvard University
   index: 3
 - name: Department of Computer Science, University of British Columbia
   index: 4
date: 
bibliography: paper.bib
---


Introduction
============

The growth of programming in the sciences has been explosive in the last
decade. This has facilitated the rapid advancement of science through
the agile development of computational tools. However, concerns have
begun to surface about the reproducibility of scientific research in
general [@Baker2016] and the potential issues stemming from issues
with analytical software [@Stodden2018]. Specifically, there is a
growing recognition across disciplines that simply making data and
software "available" is not enough and that there is a need to improve
the transparency and stability of scientific software [@Pasquier2018].

At the core of the growth of scientific computation, the `R` statistical
programming language has grown exponentially to become one of the top
ten programming languages in use today. At its root R is a
*statistical* programming language. That is, it was designed for use in
analytical workflows, and the majority of the R community is focused on
producing code for idiosyncratic projects that are *results* oriented.
Also, R's design is intentionally at a level that abstracts many aspects
of programming that would otherwise act as a barrier to entry for many
users. This is good in that there are many people who use R
with little to no formal training in computer science or
software engineering, but these same users can also be frequently
frustrated by code that is fragile, buggy, and complicated enough to
quickly become obtuse even to the authors. The stability,
reproducibility, and re-use of scientific analyses in R would be improved
by refactoring, which is a common practice in software engineering
[@Martin2009CleanCraftsmanship]. From this perspective, tools that can
lower the time and energy required to refactor analytical scripts and
otherwise help to "clean" code, but abstracted enough to be easily
accessible, could have a significant impact on scientific
reproducibility across all disciplines [@Visser2015].

To provide support for easier refactoring in R, we have created
`Rclean`. The `Rclean` package provides tools to automatically reduce a
script to the parts that are specifically relevant to a research product
(e.g., a scientific report, academic talk, research article, etc.)
Although potentially useful to all R coders, it was designed to ease
refactoring for scientists who use R but do not have formal training in
software engineering.

Methods
=======

The goal of `Rclean` is to provide a set of tools that help someone
reduce and organize code based on results. More often then not, when
someone is writing an R script, the intent is to produce a set of
results, such as a statistical analysis, figure, table, etc. This set of
results is always a subset of a much larger set of possible ways to
explore a dataset, as there are many statistical approaches and tests,
let alone ways to create visualizations and other representations of
patterns in data. This commonly leads to lengthy, complicated scripts
from which researchers manually subset results, but likely never to be
refactored because of the difficulty in disentangling the code needed to
produce some results and not others. The 'Rclean' package uses an
automated technique based on data provenance to analyze existing scripts
and provide ways to identify and extract code to produce a desired
output.

Data Provenance
---------------

All of these processes rely on the generation of data provenance. The
term provenance means information about the origins of some object. Data
provenance is a formal representation of the execution of a
computational process (<https://www.w3.org/TR/prov-dm/>), to rigorously
determine the the unique computational pathway from inputs to results
[@Carata2014]. To avoid confusion, note that "data" in this context is
used in a broad sense to include all of the information generated during
computation, not just the data that are collected in a research project
that are used as input to an analysis. Having the formalized,
mathematically rigorous representation that data provenance provides
guarantees that the analyses that `Rclean` conducts are theoretically
sound. Most importantly, because the relationships defined by the
provenance can be represented as a graph, it is possible to apply
network search algorithms to determine the minimum and sufficient code
needed to generate the chosen result in the `clean` function.

There are multiple approaches to collecting data provenance, but
`Rclean` uses "prospective" provenance, which analyzes code and uses
language specific information to predict the relationship among
processes and data objects. `Rclean` relies on a library called
`CodeDepends` to gather the prospective provenance for each script. For
more information on the mechanics of the `CodeDepends` package, see
[@R-CodeDepends]. To get an idea of what data provenance is, take a
look at the `code_graph` function. The plot that it generates is a
graphical representation of the prospective provenance generated for
`Rclean` .

![Network diagram of the prospective data provenance generated for an
example script. Arrows indicate which lines of code (numbered) produced
which objects
(named).](paper_files/figure-markdown_strict/prov-graph-1.png)

In the future, it would also be useful to extend the existing framework
to support other provenance methods. One such possibility is
*retrospective provenance*, which tracks a computational process as it
is executing. Through this active, concurrent monitoring, retrospective
provenance can gather information that static prospective provenance
can't. Greater details of the computational process would enable other
features that could address some challenges, such as processing
information from comments, parsing control statements, and replicating
random processes. However, using retrospective provenance comes at a
cost. In order to gather it, the script needs to be executed. When
scripts are computationally intensive or contain bugs that stop
execution, then retrospective provenance can not be obtained for part or
all of the code. Some work has already been done in the direction of
implementing retrospective provenance for code cleaning in R (see
<http://end-to-end-provenance.github.io>).

Software Availability
---------------------

The software is currently hosted on GitHub, and we recommend using the
`devtools` library [@R-devtools] to install directly from the
repository (<https://github.com/ROpenSci/Rclean>). The package is
open-source and welcomes contributions. Please visit the repository page
to report issues, request features or provide other feedback.

Discussion
==========

We see promise in connecting `Rclean` with other clean code and
reproducibility tools. One example is the `reprex` package, which
provides a simple API for sharing reproducible examples [@R-reprex].
Another possibility is to help transition scripts to function, package
and workflow creation and refactoring via toolboxes like `drake`
[@R-drake]. `Rclean` could provide a reliable way to extract parts of
a larger script that would be piped to a simplified reproducible
example, in the case of `reprex`, or, since it can isolate the code from
inputs to one or more outputs, be used to extract all of the components
needed to write one or more functions that would be a part of a package
or workflow, as is the goal of `drake`. To conclude, we hope that
`Rclean` makes writing scientific software easier for the R community.
We look forward to feedback and help with extending its application,
particularly in the area of reproducibility. To get involved, report
bugs, suggest features, please visit the project page.

Acknowledgments
===============

This work was improved by discussions with ecologists at Harvard Forest
and through the helpful review provided by the ROpenSci community,
particularly Anna Krystalli, Will Landau, and Clemens Schmid. Much of the
work was funded by US National Science Foundation grant SSI-1450277 for
applications of End-to-End Data Provenance.

References
==========

<!-- Use overleaf papers.bib + knitr::write_bib -->
Title: Rclean: A Tool for Writing Cleaner, More Transparent Code

Submitting Author: Matthew K. Lau (@mklau)
Repository: https://github.com/provtools/rclean

---

-   Paste the full DESCRIPTION file inside a code block below:

```
Type: Package
Package: Rclean
Title: A Tool for Writing Cleaner, More Transparent Code
Version: 1.1.0
Date: 2019-04-24
Author: Matthew K. Lau
Maintainer: Matthew K. Lau <matthewklau@fas.harvard.edu>
Description: To create clearer, more concise code provides this
	     toolbox helps coders to isolate the essential parts of a script that
	     produces a chosen result, such as an object, tables and figures
	     written to disk and even warnings and errors. 
URL: https://github.com/ProvTools/Rclean
BugReports: https://github.com/ProvTools/Rclean/issues
License: GPL-3 | file LICENSE
Imports: igraph, jsonlite, formatR, CodeDepends
Suggests: roxygen2, testthat
RoxygenNote: 6.0.1

```


## Scope 

- Please indicate which category or categories from our [package fit policies](https://ropensci.github.io/dev_guide/policies.html#package-categories) this package falls under: (Please check an appropriate box below.:

	- [ ] data retrieval
	- [ ] data extraction
	- [ ] database access
	- [ ] data munging
	- [ ] data deposition
	- [X] reproducibility
	- [ ] geospatial data
	- [ ] text analysis
	

- Explain how the and why the package falls under these categories (briefly, 1-2 sentences).  Please note any areas you are unsure of:

In writing analytical scripts, software best practices are often a
lower priority than producing inferential results, leading to large,
complicated code bases that often need refactoring. The "code
cleaning" capabilities of the Rclean package provide a means to
rigorously identify the minimal code required to produce a given
result (e.g. object, table, plot, etc.), reducing the effort required
to create simpler, more transparent code that is easier to reproduce.

-   Who is the target audience and what are scientific applications of
    this package?
	
The target audience is domain scientists that have little to no formal
training in software engineering. Multiple studies on scientific
reproducibility have pointed to data and software availability as
limiting factors. This tool will provide an easy to use tool for
writing cleaner analytical code.

-   Are there other R packages that accomplish the same thing? If so,
    how does yours differ or meet [our criteria for
    best-in-category](https://ropensci.github.io/dev_guide/policies.html#overlap)?

There are other packages that analyze the syntax and structure of
code, such as lintr, formatr and cleanr. Rclean, as far as we are
aware, is the only package written for R that uses a data provenance
approach to construct the interdependencies of objects and functions
and then uses graph analytics to rigorously determine the desired
pathways to determine the minimal code-base needed to generate an
result.

-  Any other questions or issues we should be aware of?:

Not that I can think of at the moment.
# General Response

Thanks again to @wlandau and @nevrome for reviewing the package and
providing insightful comments. I have updated the package greatly
based on the various suggestions and questions. These updates are
documented in the package's github repo issue and change logs. I
detailed the changes as they pertain to the reviewers' comments point
by point below. For unchecked, please see the comments that follow. In
addition, I would like to point out that the package repo is now
hosted at https://gitub.com/MKLau/Rclean, which is the result of
shifts in how the project is now being managed. 

# Response to Reviewer Comments
## Documentation

- [x] @wlandau Please include a CODE_OF_CONDUCT.md file at the top
    level in your repo. You could generate one
    [here](https://github.com/ProvTools/Rclean/community) or call
    usethis::use_code_of_conduct() . On
    [this page](https://github.com/ProvTools/Rclean/community), you
    can even generate issue and pull request templates that remind
    people of the code of conduct.

    Added CODE_OF_CONDUCT.md at the top level via
    usethis::use_code_of_conduct().

- [x] @wlandau Consider renaming CONTRIBUTE.md to the more widely-used
     CONTRIBUTING.md .

    Changed to CONTRIBUTING.md.

- [x] @wlandau Do you have a longer, real-world script to show?
    simple_script.R gets the idea across, but if you show how to
    extract simple chunks from a confusing script with 1000+ lines,
    you could make Rclean really shine.

    I agree, and I added a longer, more realistic, script to the
    example scripts, see "examples/long_script.R". I used this in the
    vignette to demonstrate the utility of retrieving the code to
    produce a specific result when it is embedded in a longer,
    somewhat complicated script. I wasn't able to simply use an
    existing script from a real project because of issues with data
    dependencies. If the reviewer has a suggestion of an example, I
    would be happy to look at including it in a future update.

- [x] @wlandau At the top of the README, please consider reorganizing
    those bullet points into two paragraphs or subheadings: one
    introducing the existing problems in scripts and reproducibility,
    and then another focusing on how Rclean solves those
    problems. Maybe also elaborate on both points as well. I do not
    think this is the place to mention the provenance model or graphs,
    but you could definitely talk more about the survey you mentioned,
    the problem of messy code in science, and the fact that Rclean
    sifts through code to extract just what you need for a given
    result.

    I have edited and re-organized the content of the README to focus
    on the purpose of the package, installation and a quick start
    guide. The README refers users to the vignette for more information.

- [x] @wlandau Would you cite that survey from the third bullet in the
    README? If you cannot find the citation, feel free to remove the
    mention entirely and choose a different way to motivate the
    problem.

    > A recent survey of over 1500 scientists reported a crisis of
    reproducibility with “selective reporting” being the most cited
    contributing factor and 80% saying code availability is playing a
    role.

    This reference has been removed from the README.

- [x] @wlandau You might consider using a README.Rmd to generate the
    README.md so it is easier to refresh the docs when you change the
    package (example
    [here](https://github.com/tidyverse/dplyr/blob/master/README.Rmd)).

    Updated to use README.Rmd.

- [x] @wlandau I strongly agree with @nevrome to add examples to the
    help files for the exported functions you choose to keep (see
    below). It is difficult for me to understand how to create
    arguments for most of your exported functions.

    Exported functions now have simple examples to illustrate usage.

- [x] @wlandau The
    [vignette](https://github.com/ProvTools/Rclean/blob/master/vignettes/Rclean.rmd)
    repeats a lot of information from the README. Let's think about
    removing the redundant sections and refocusing the vignette on a
    deeper dive in a feature set not covered in the README. Perhaps
    retrospective provenance ?  See the next comment.

    README.md now functions as a quick start guide, and more detailed
    information has been moved to the vignette and JOSS paper.

- [x] @wlandau It is not intuitively clear from the docs why we need
    retrospective provenance. In the README or vignette, please
    demonstrate a situation where prospective provenance is not enough
    and retrospective provenance gives us what we seek.

    Thank you for this question, it has really pushed the development
    of a more focused package. I have removed the retrospective
    provenance functionality from the package, as I came to realize
    that it was unnecessarily slowing Rclean and complicating
    usage. The code to use retrospective provenance has been
    integrated into
    https://github.com/end-to-end-provenance/provClean, but it is now
    not used for Rclean. However, to explain the difference,
    retrospective provenance gets the lineage of a result by
    monitoring the execution of code, while prospective provenance
    infers similar information from the code prior to
    execution. Because of this difference, retrospective provenance
    can do things that prospective provenance can't, such as observe
    and record random numbers and which paths were taken through
    control statements. I have clarified this information in the
    vignette and similar information will be added to the JOSS
    article.

- [x] @wlandau Please consider citing a paper to make it clear that
    "retrospective" and "prospective" are widely used and accepted
    terms to describe provenance. A quick search reveals a couple
    possibilities (examples
    [here](https://ieeexplore.ieee.org/document/5557202) and
    [here](https://www.researchgate.net/publication/221212739_Provenance_and_scientific_workflows_Challenges_and_opportunities))?

    See the previous comment. Thanks for the references, I have added
    provenance explanations and other references to the vignette.


- [x] @nevrome Please add some more descriptive text and example code
    for clean , codeGraph and write.code . I would like to see
    codeGraph mentioned in the README/vignette.
    
    Examples are now present for exported functions and more
    expository text has been added.

- [x] @nevrome I would love to have an **example script** (maybe
     example/simple_script.R ) as a
     [package example dataset](http://r-pkgs.had.co.nz/data.html) that
     allows to explore the functions without the need to
     provide/prepare an own script. An example script will also allow
     you to add meaningful example code for the functions.

    Two example scripts, simple_script and long_script, are now
    available as R package data files. However, as clean() requires
    a file path, they cannot be passed directly into clean().

- [x] @nevrome **Retrospective Provenance** remained a mystery for
    me. The missing piece of the puzzle is how to generate the
    relevant json file. Maybe I just missed the crucial hint? I would
    prefer to see a more extensive explanation of this feature and how
    to use it in the README/vignette.

    See my explanation and comments above. 

- [x] @nevrome In the context of the Retrospective Provenance feature
    you use the ** option interface** to access the content of the
    json file. I do not understand why this rather unusual solution is
    necessary. Why can't the rd parameter simply get the relevant
    character vector? Instead of `options(prov.json =
    readLines("example/prov_micro.json")); clean(file =
    "example/micro.R", rp = TRUE) the interface could simply be
    clean(file = "example/micro.R", rp =
    readLines("example/prov_micro.json"))` or even just the path
    clean(file = "example/micro.R", rp = "example/prov_micro.json") .
    Why did you decide against this?

    Similarly here, see previous provenance comments.

- [x] @nevrome You could add a **help page** for the whole package as
    described [here](http://r-pkgs.had.co.nz/man.html#man-packages). I
    often use this page as an entry point for the function index.

    A whole package help page has been added. 

- [x] @nevrome This might be way beyond the purpose of this review,
    but as I personally consider these kind of things when I decide
    whether I want to use a package, I would like to add that I was
    not able to access neither your **homepage** http://mklau.info/
    nor the one the Github Orga ProvTools http://provtools.org/
    [2019-09-24].
    
    Thanks for checking these links. The provtools website is now
    defunct and materials have been moved to
    github.com/mklau/rclean. However, mklau.info is still supported
    and is working for me. Let me know if you happen to try accessing
    it again and have an issue.

## JOSS Co-Submission

@wlandau: missing all sections
@nevrome: missing short summary

- [ ] The **JOSS** paper manuscript seems to be in an unfinished state
    because several sections are mentioned with a keyword but not yet
    written. The text available so far is fine.

@mklau: The initial submission of the JOSS paper was based on an
example of a submitted article present on the JOSS
instructions. However, after this package review, there will be
substantial content added to the article. I am currently working on
these updates to the text.


## Functionality

- [x] @wlandau: installation and package guidelines missing.

   @mklau: This has been added to the README.

### API

- [x] @wlandau The clean() function tries to do three things: (1) list
    lines of execution, (2) list the names of possible results, and
    (3) get the code required to produce one of the results. I agree
    that clean() should do (3), but I think (1) and (2) are not really
    cleaning yet and should go in their own functions. Maybe
    list_lines() (which returns an integer vector) and list_results()?

    The clean function has been significantly refactored. Now that the
    retrospective provenance has been removed from the package, clean
    has been streamlined. It still returns a list of possible results,
    if no results are inputed, but the focus in now on getting the
    script and a variable or several variables to determine and return
    the minimal code. There is now an exported function get_vars to
    return a list of variables in a given script.


- [x] @wlandau Similarly, it is also worth considering an entirely new
    set of functions for retrospective provenance because they have
    different inputs and
    [your underlying code is entirely different](https://github.com/ProvTools/Rclean/blob/f6b7ffc3325638ba878d7bb00df28a630043920d/R/clean.R#L86-L182).

    Already discussed above, this feature has been removed from Rclean.

- [x] @wlandau In help(package = "Rclean") , I see several API
    functions that represent very technical steps and need esoteric
    data structures. Are there specific reasons why you exported these
    things to the user? I am not sure users will want to handle matrix
    representations of graphs or PROV-JSON-formatted provenance
    objects in memory. It seems like the easiest inputs to understand
    are R script files and JSON files. Similarly, the best outputs are
    human-readable: cleaned-up scripts and plotted code graph
    visualizations. Maybe even consider writing files directly from
    clean() .  Consistency and simplicity are extremely valuable.

    Again, the API has been substantially simplified. Specific to this
    comment, JSON is no longer required because the package no longer
    uses retrospective provenance, which stores provenance in
    PROV-JSON format.

- [x] @wlandau I am questioning the need for most of Rclean 's
    exported functions.  Not all your functions need to be exported,
    and not all of them need roxygen2 docstrings.

    Good suggestion. The API has been updated with a reduction in the
    exported functions, and significant refactoring has been done to
    the code base, particularly clean.R, which now has internal
    functions included in the clean.R file itself.

- [x] @wlandau Do we really need p.spine() and p.spine() ? If you add
    a function to get the igraph from a script, it seems like you can
    simply instruct users to call dfs() themselves.

    Per the previous comment, these functions (p.spine and get.spine)
    have been refactored and integrated into clean.R. In order to keep
    the workflow and API simple for users, manually conducting the
    depth first search is not described for the users. This could be
    done for some more advanced applications, but currently I can't
    see a purpose for exposing the user to this complexity.

- [x] @wlandau Should get.libs() , var.lineage() , and var.id() really
    be exported?  clean() needs them, but do you expect most users to
    ever invoke them manually?

    var.id() has been removed. var.lineage(), now var_lineage(), has
    been internalized. get.libs(), now get_libs(), is still a part of
    the API for users, as it allows users to quickly dig out the
    libraries used in code. 

- [x] @wlandau Do we really need a write.code() function?  readLines()
    already writes code to a file, and clipr::write_clip() writes to
    the clipboard. Optional: maybe clean() itself could have an
    argument to control where the output goes: either to the clipboard
    with write_clip() or to a file if a path is given.

    Thanks for these suggestions. write.code() has become keep(),
    which now uses the clipr package. Although it may seem like a
    small step, my sense in using Rclean is that abstracting
    writeLines(), which requires the creation of a file connection,
    and combining its functionality with clipr's write_clip(),
    substantially eases the management of cleaned code. 

- [x] @wlandau You explicitly list parse.info() , read.prov() , and
     parse.graph() as internal functions. Perhaps these should not be
     exported either.

     I agree here, however, as discussed above, this feature has been
     removed.

- [x] @wlandau The names of exported functions have different naming
     conventions:  codeGraph() (camel case) vs prov_json() (snake
     case) vs get.spine() (dot-separated). For the exported functions
     you choose to keep, please consider sticking to one of these
     conventions (I recommend snake case).  Depending on the user base
     you have right now, you may want to gracefully deprecate those
     functions as
     [described here](https://ropensci.org/technotes/2017/01/05/package-evolution/).

     Good suggestion. Switched to snake case throughout. 

### Installation and automated checks

- [x] @wlandau Re
    https://github.com/ropensci/software-review/issues/327#issuecomment-520062929,
    I also had some trouble installing Rclean . I needed to install
    Rgraphviz directly from Bioconductor. I suggest adding the
    following to the installation instructions in the README. At the
    very least, it is a useful stopgap until @annakrystalli's pull
    request to CodeDepends reaches CRAN.
    
    ```r
    install.packages("BiocManager")
    BiocManager::install(c("graph", "Rgraphviz"))
    ```

    This has been added to the README and vignette. 

- [x] @wlandau devtools::spell_check("Rclean") identifies several
    correctly spelled words. I recommend listing them in an
    inst/WORDLIST file
    ([example here](https://github.com/ropensci/drake/blob/master/inst/WORDLIST)).

    Thanks for this. I have added inst/wordlist.

- [x] @wlandau goddpractice::gp("Rclean") suggests a bunch of style
    fixes I agree with. The styler package can help you turn
    assignment =` into <-`.  Please let me know if you would like
    clarification on any of these notes.
    
    ```r
    It is good practice to
    
      ✖ write unit tests for all functions, and all package code
        in general. 86% of code lines are covered by test cases.
    
        R/clean.R:49:NA
        R/clean.R:72:NA
        R/clean.R:82:NA
        R/clean.R:91:NA
        R/clean.R:95:NA
        ... and 21 more lines
    
      ✖ use '<-' for assignment instead of '='. '<-' is the
        standard, and R users and developers are used it and it is easier
        to read your code for them if you use '<-'.
    
        R/var.id.R:41:44
        R/var.id.R:42:47
    
      ✖ avoid long code lines, it is bad for readability. Also,
        many people prefer editor windows that are about 80 characters
        wide. Try make your lines shorter than 80 characters
    
        R/clean.R:47:1
        R/clean.R:93:1
        R/get.spine.R:44:1
        R/parse.info.R:39:1
        R/parse.info.R:40:1
        ... and 1 more lines
    
      ✖ avoid sapply(), it is not type safe. It might return a
        vector, or a list, depending on the input data. Consider using
        vapply() instead.
    
        R/clean.R:150:29
        R/clean.R:151:55
        R/clean.R:169:34
    
      ✖ avoid 1:length(...), 1:nrow(...), 1:ncol(...),
        1:NROW(...) and 1:NCOL(...) expressions. They are error prone and
        result 1:0 if the expression on the right hand side is zero. Use
        seq_len() or seq_along() instead.
    
        R/get.libs.R:40:13
        R/get.libs.R:42:17
        R/get.libs.R:44:21
        R/var.id.R:39:13
        R/var.lineage.R:43:22
        ... and 4 more lines
    ```

    These have either been fixed or deleted from the
    package. goodpractice::gp() currently reports "Ah! Perfect
    package! Keep up the cat's pajamas work!".

### Miscellaneous

- [x] @wlandau I am a bit concerned about the scope of the built-in
    unit tests. It would be helpful to test at least one longer script
    to make sure the correct lines get captured in a serious
    workflow. Not only that, but when you actually run the longer
    script and various object-specific scripts generated by clean() ,
    do the output data objects agree?

    A longer script, that is arguably more realistic, has been added
    along with associated tests to check that minimized code for an
    entangled result can be correctly extracted.

- [x] @wlandau What is the reason you chose JSON to store
    retrospective provenance?  Have you considered a tidier format
    such as a long-form data frame? I like how
    [ pkgapi ](https://github.com/r-lib/pkgapi) represents its output.

    This functionality has been removed, but the PROV-JSON format was
    chosen to comply with standards used by the retrospective
    provenance community.

- [x] @wlandau The clean() function is ambitious, long, and deeply
    indented.  Please consider decomposing this monolithic function
    into smaller helper functions. No need to export the helper
    functions or write roxygen2 docstrings. A particularly good
    resource on this kind of refactoring is Jenny Bryan's talk on
    [code smells and feels](https://www.youtube.com/watch?v=7oyiPBjLAWY).

    As mentioned previously, clean() has been refactored and the API
    as a whole has been re-organized to make it more transparent and
    approachable.

- [x] @wlandau Are you open to using
    [ styler ](https://github.com/r-lib/styler) instead of formatR to
    reformat code? It is more modern and more actively maintained.

    Given the higher maintenance activity, Rclean now uses styler.


- [x] @nevrome When I tried to apply clean to some of my own scripts,
    I realized that all packages necessary to run this script have to
    be loaded in the current session. If this was not the case, then I
    got the error message **`Error in as.environment(pos): no item
    called "package:dplyr" on the search list`**.  I got this message
    10+ times for one of my more complex scripts before clean went
    through, until I loaded all the packages. This is a little bit
    tedious. Is there another way to implement package loading? Maybe
    the [automagic](https://github.com/cole-brokamp/automagic) package
    could help to automate that? If this is not possible or if you
    don't want to have an automatic solution, then I would at least
    ask you to catch the current error message and replace it with
    something more intelligible, e.g.  "Package xyz has to be loaded
    to analyze this script.".  

    Sorry you ran into this hurdle with the package. This should now
    no longer be an issue, as it was a product of how the get.libs()
    function found the packages that were being used in the script and
    the need for the retrospective provenance to run the code prior to
    analyzing it. Now, get_libs() no longer needs the packages to be
    loaded and retrospective provenance is no longer used. 

- [x] @nevrome **Internal functions** do not need to be exported (e.g.
    parse.info , read.prov ). Remove the @export tag for these
    functions and add the tags @keywords internal and @noRd .
    
    Agreed, the API has been updated and @noRd has been added where needed.

- [x] @nevrome The ** simple_script.R ** code appears many times in
    the package (in the directories test/ , inst/ , exec/ , example/
    ). I don't think you need this file so many times.

    **simple_script.R** now only occurs in inst/examples.

- [ ] @nevrome Although you not seem to work in RStudio I suggest to
    add a **Rstudio project file** in the root directory (+ the
    necessary additions to .Rbuildignore ). That simplifies
    contributing for RStudio users.

    Yes, I don't work in RStudio, except when I'm teaching. So, I'm
    not sure how to best set this up for using with the package
    dev. If someone who uses RStudio starts to contribute, I would be
    happy to have it added.

- [x] @nevrome In the **travis-ci** setup I suggest to treat warnings
    as errors ( warnings_are_errors: true ). I think this is usual for
    R packages, because no package with warnings can go to CRAN.

    This has been updated to treat warnings as errors.. 

