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
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
ezknitr - Avoid the typical working directory pain when using 'knitr'
=====================================================================

[![R Build Status](https://github.com/ropensci/ezknitr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ezknitr/actions)
[![CRAN
version](http://www.r-pkg.org/badges/version/ezknitr)](https://cran.r-project.org/package=ezknitr)
[![codecov.io](https://codecov.io/github/ropensci/ezknitr/coverage.svg?branch=master)](https://codecov.io/github/ropensci/ezknitr?branch=master)

> *Copyright 2016 [Dean Attali](http://deanattali.com). Licensed under
> the MIT license.*

[`knitr`](https://github.com/yihui/knitr) is a popular package for
generating dynamic reports in R using the concept of [literate
programming](http://www.literateprogramming.com/knuthweb.pdf). `ezknitr`
is an extension of `knitr` that adds flexibility in several ways and
solves a few issues that are commonly encountered with `knitr`.

One common source of frustration with `knitr` is that it assumes the
directory where the source file lives should be the working directory,
which is often not true. `ezknitr` addresses this problem by giving you
complete control over where all the inputs and outputs are, and adds
several other convenient features. The two main functions are `ezknit()`
and `ezspin()`, which are wrappers around `knitr`'s `knit()` and
`spin()`, used to make rendering markdown/HTML documents easier.

Table of contents
-----------------

-   [Installation](#installation)
-   [Overview](#overview)
-   [Motivation & simple use case](#motivation-simple-usecase)
-   [Use case: using one script to analyze multiple
    datasets](#usecase-advanced)
-   [Experiment with ezknitr](#experiment)
-   [spin() vs knit()](#spin-vs-knit)
-   [Using rmarkdown::render()](#using-render)

<h2 id="installation">
Installation
</h2>
`ezknitr` is available through both CRAN and GitHub.

To install the CRAN version:

    install.packages("ezknitr")

To install the latest developmental version from GitHub:

    install.packages("devtools")
    devtools::install_github("ropensci/ezknitr")

<h2 id="overview">
Overview
</h2>
If you have a very simple project with a flat directory structure, then
`knitr` works great. But even something as simple as trying to knit a
document that reads a file from a different directory or placing the
output rendered files in a different folder cannot be easily done with
`knitr`.

`ezknitr` improves basic `knitr` functionality in a few ways. You get to
decide:

-   What the working directory of the source file is
    -   Default is your current working directory, which often makes
        more sense than the `knitr` assumption that the working
        directory is wherever the input file is
-   Where the output files will go
    -   With `knitr`, all the rendered output files will be generated in
        the folder you're currently in
-   Where the figures generated in the markdown document will be stored
    -   `knitr` makes it cumbersome to change this directory
-   Any parameters to pass to the source file
    -   Useful if you want to run an identical source file multiple
        times with different parameters

<h2 id="motivation-simple-usecase">
Motivation & simple use case
</h2>
Assume you have an Rmarkdown file that reads a data file and produces a
short report while also generating a figure. Native `knit()` (or
`spin()` if you're starting with an R script instead of an Rmd file)
works great if you have a flat directory structure like this:

    - project/
      |- input.csv
      |- report.Rmd

But what happens if you have a slightly more complex structure? In a
real project, you rarely have everything just lying around in the same
folder. Here is an example of a more realistic initial directory
structure (assume the working directory/project root is set to
`project/`):

    - project/
      |- analysis/
         |- report.Rmd
      |- data/
         |- input.csv

Now if you want `knitr` to work, you'd have to ensure the path to
`input.csv` is relative to the `analysis/` directory because that's
where the Rmd file is. This is counter-intuitive because most people
expect to create paths relative to the *working directory/project root*
(`project/` in this case), but `knitr` will use the `analysis/` folder
as the working directory. Any code reading the input file needs to use
`../data/input.csv` instead of `data/input.csv`.

Other than being confusing, it also means that you cannot naively run
the Rmd code chunks manually because when you run the code in the
console, your working directory is not set to what `knitr` will use as
the working directory. More specifically, if you try to run the command
that reads the input file, your console will look in
`project/../data/input.csv` (which doesn't exist).

A similar problem arises when you want to create files in your report:
`knitr` will create the files relative to where the Rmd file is, rather
than relative to the project root.

Another problem with the flat directory structure is that you may want
to control where the resulting reports get generated. `knitr` will
create all the outputs in your working directory, and as far as I know
there is no way to control that.

`ezknitr` addresses these issues, and more. It provides wrappers to
`knit()` and `spin()` that allow you to set the working directory for
the input file, and also uses a more sensible default working directory:
the current working directory. `ezknitr` also lets you decide where the
output files and output figures will be generated, and uses a better
default path for the output files: the directory containing the input
file.

Assuming your working directory is currently set to the `project/`
directory, you could use the following `ezknitr` command to do what you
want:

    library(ezknitr)
    ezknit(file = "analysis/report.Rmd", out_dir = "reports", fig_dir = "myfigs")

    - project/
      |- analysis/
         |- report.Rmd
      |- data/
         |- input.csv
      |- reports/
         |- myfigs/
            |- fig1.png
         |- report.md
         |- report.HTML

We didn't explicitly have to set the working direcory, but you can use
the `wd` argument if you do require a different directory (for example,
if you are running this from some build script or from any arbitrary
directory). After running `ezknit()`, you can run `open_output_dir()` to
open the output directory in your file browser if you want to easily see
the resulting report. Getting a similar directory structure with `knitr`
is not simple, but with `ezknitr` it's trivial.

Note that `ezknitr` produces both a markdown and an HTML file for each
report (you can choose to discard them with the `keep_md` and
`keep_html` arguments).

<h2 id="usecase-advanced">
Use case: using one script to analyze multiple datasets
</h2>
As an example of a more complex realistic scenario where `ezknitr` would
be useful, imagine having multiple analysis scripts, with each one
needing to be run on multiple datasets. Being the organizer scientist
that you are, you want to be able to run each analysis on each dataset,
and keep the results neatly organized. I personally was involved in a
few projects requiring exactly this, and `ezknitr` was in fact born for
solving this exact issue. Assume you have the following files in your
project:

    - project/
      |- analysis/
         |- calculate.Rmd
         |- explore.Rmd
      |- data/
         |- human.dat
         |- mouse.dat

We can easily use `ezknitr` to run any of the analysis Rmarkdowns on any
of the datasets and assign the results to a unique output. Let's assume
that each analysis script expects there to be a variable named
`DATASET_NAME` inside the script that tells the script what data to
operate on. The following `ezknitr` code illustrates how to achieve the
desired output.

    library(ezknitr)
    ezknit(file = "analysis/explore.Rmd", out_dir = "reports/human",
            params = list("DATASET_NAME" = "human.dat"), keep_html = FALSE)
    ezknit(file = "analysis/explore.Rmd", out_dir = "reports/mouse",
            params = list("DATASET_NAME" = "mouse.dat"), keep_html = FALSE)
    ezknit(file = "analysis/calculate.Rmd", out_dir = "reports/mouse",
            params = list("DATASET_NAME" = "mouse.dat"), keep_html = FALSE)

    - project/
      |- analysis/
         |- calculate.Rmd
         |- explore.Rmd
      |- data/
         |- human.dat
         |- mouse.dat
      |- reports/
         |- human/
            |- explore.md
         |- mouse/
            |- calculate.md
            |- explore.md

Note that this example uses the `params = list()` argument, which lets
you pass variables to the input Rmarkdown. In this case, I use it to
tell the Rmarkdown what dataset to use, and the Rmarkdown assumes a
`DATASET_NAME` variable exists. This of course means that the analysis
script has an external dependency by having a variable that is not
defined inside of it. You can use the `set_default_params()` function
inside the Rmarkdown to ensure the variable uses a default value if none
was provided.

Also note that differentiating the species in the output could also have
been done using the `out_suffix` argument instead of the `out_dir`
argument. For example, using `out_suffix = "human"` would have resulted
in an ouput file named `explore-human.md`.

<h2 id="experiment">
Experiment with ezknitr
</h2>
After installing and loading the package (`library(ezknitr)`), you can
experiment with `ezknitr` using the `setup_ezknit_test()` or
`setup_ezspin_test()` functions to see their benefits. See
`?setup_ezknit_test` for more information.

<h2 id="spin-vs-knit">
spin() vs knit()
</h2>
`knit()` is the most popular and well-known function from `knitr`. It
lets you create a markdown document from an Rmarkdown file. You can
learn more about `knit()` on [the official knitr
page](http://yihui.name/knitr).

`spin()` is similar, but starts one step further back: it takes an R
script as input, creates an Rmarkdown document from the R script, and
then proceeds to create a markdown document from it. `spin()` can be
useful in situations where you develop a large R script and want to be
able to produce reports from it directly instead of having to copy
chunks into a separate Rmarkdown file. You can read more about why I
like `spin()` in the blog post ["knitr's best hidden gem:
spin"](http://deanattali.com/2015/03/24/knitrs-best-hidden-gem-spin/).

<h2 id="using-render">
Using rmarkdown::render()
</h2>
When the core of this package was developed, none of its functionality
was supported in any way by either `knitr` or `rmarkdown`. Over time,
`rmarkdown::render()` got some new features that are very similar to
features of `ezknitr`. Native support for parameters inside Rmarkdown
files using YAML is a big feature which makes the use of
`set_default_params()` and the `params` argument of `ezknitr` less
important. However, the core problem that `ezknitr` wants to solve is
the working directory issue, and this issue has yet to be addressed by
`rmarkdown` or `knitr`, which makes `ezknitr` still useful.

Please note that this project is released with a [Contributor Code of
Conduct](CONDUCT.md). By participating in this project you agree to
abide by its terms.

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
# ezknitr 0.6.1 2016-12-19

- GitHub repository ownership transferred from ropenscilabs to ropensci


# ezknitr 0.6 2016-09-15

- GitHub repository ownership transferred from daattali to rOpenSciLabs

# ezknitr 0.5 2016-08-21

- added JOSS paper
- revisions based on rOpenSci feedback, mostly modifications to README and other documentation

# ezknitr 0.4 2016-06-12

- added `...` argument to `ezspin` that passes any additional arguments to `knitr::spin` (#6)
- add param `move_intermediate_file` to `ezspin` that controls whether the intermediate Rmd file stays with the source R script or moves to the output folder (thanks @klmr)  

# ezknitr 0.3.0

2015-12-07

- add section about `rmarkdown::render` in vignette/readme

# ezknitr 0.2.0

2015-12-07

- add unit tests
- documentation additions, getting ready for CRAN release

# ezknitr 0.1.0

2015-12-05

- add `setup_ezknit_test` function
- add `open_output_dir` function

# ezknitr 0.0.0.9002

2015-11-06

- rename `ezrender` to `ezknitr`
- many internal changes and code improvements
- no longer need to delete an empty wrong figures folder because knitr fixed its bug #942 

# ezknitr 0.0.0.9001

2015-10-23

- `ezspin()` gains new arguments `keepRmd` (default: `FALSE`) and `keepMd` (default: `TRUE`) (#1)

# ezknitr 0.0.0.9000

- Import from `rsalad` package
# Round 1

## Test environments

* Windows 7, R 3.2.2 (local)
* ubuntu 12.04, R 3.2.2 (on travis-ci)
* ubuntu 14.04, R 3.2.1 (on my DigitalOcean droplet)

## Submission comments

2015-12-07

R CMD check has no ERRORs or WARNINGs, and 1 NOTE regarding an invalid URL pointing to the CRAN URL of this package. That's expected because this package is not yet on CRAN, but the URL will be valid once it is on CRAN.
  

## Reviewer comments

2015-12-08 Kurt Hornik

Thanks, on CRAN now.

Best

-k


---

# v0.4

# Round 1

## Submission comments

2016-06-12

R CMD check ran without any errors, warnings, or notes
  

## Reviewer comments

2016-06-12 Uwe Ligges

Thanks, on CRAN now.

# Version 0.5

### Test environments

* Windows 7, R 3.3.1 (local)
* ubuntu 12.04, R 3.3.1 (on travis-ci)
* ubuntu 14.04, R 3.2.1 (on my DigitalOcean droplet)

## Round 1

### Submission comments

2016-09-13

No errors, warnings, or notes

### Reviewer comments

2016-09-14 Kurt Hornik

Thanks, on CRAN now
---
title: "ezknitr: Avoid the Typical Working Directory Pain When Using 'knitr'"
tags:
  - reproducible research
  - R
authors:
 - name: Dean Attali
   orcid: 0000-0002-5645-3493
   affiliation: UBC
date: 7 August 2016
bibliography: paper.bib
---

# Summary

'knitr' [@R-knitr] is a popular R package [@R-base] that implements the concept of literate programming [@knuth1984literate] in R. It has been widely adopted in the R community as a means of generating dynamic reports. Despite its popularity, one common source of frustration with 'knitr' is that it makes important assumptions about the directory structure of the input files and does not offer much flexibility in deciding where the different output files are generated. The 'ezknitr' package [@ezknitr-archive] is an R package that addresses these problems and enhances the functionality of 'knitr'. 'ezknitr' provides the user with complete control over where all the inputs and outputs are, and adds several other convenient features to make the process of generating dynamic documents more streamlined and efficient.

# ReferencesTime: `r date()`

Where am I? My working directory is:

`r getwd()`


I like to assume the directory containing the 'R' and 'data' folders is the
working directory.  Let's try to read a file from the 'data' folder.

```{r }
dat <- scan("data/numbers.txt", quiet = TRUE)
```

How many numbers were read from the input?

```{r }
length(dat)
```

(If you see errors when using `knitr`, that's expected - it's because the
working directory was not what the script assumed.)

```{r setparam, echo = FALSE}
ezknitr::set_default_params(list(numPoints = 15))
```

Here's a plot of `r numPoints` random points:

```{r plot, echo = FALSE}
set.seed(100)
plot(rnorm(numPoints))
```

---
title: "Package ezknitr"
author: "Dean Attali"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Package ezknitr}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(tidy = FALSE, comment = "#>")
```

# ezknitr - Avoid the typical working directory pain when using 'knitr'

[![R Build Status](https://github.com/ropensci/ezknitr/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ezknitr/actions)
[![CRAN version](http://www.r-pkg.org/badges/version/ezknitr)](https://cran.r-project.org/package=ezknitr)

> *Copyright 2016 [Dean Attali](http://deanattali.com). Licensed under the MIT license.*

[`knitr`](https://github.com/yihui/knitr) is a popular package for generating dynamic reports in R using the concept of [literate programming](http://www.literateprogramming.com/knuthweb.pdf). `ezknitr` is an extension of `knitr` that adds flexibility in several ways and solves a few issues that are commonly encountered with `knitr`.

One common source of frustration with `knitr` is that it assumes the directory where the source file lives should be the working directory, which is often not true. `ezknitr` addresses this problem by giving you complete control over where all the inputs and outputs are, and adds several other convenient features. The two main functions are `ezknit()` and `ezspin()`, which are wrappers around `knitr`'s `knit()` and `spin()`, used to make rendering markdown/HTML documents easier. 

## Table of contents

- [Installation](#installation)
- [Overview](#overview)
- [Motivation & simple use case](#motivation-simple-usecase)
- [Use case: using one script to analyze multiple datasets](#usecase-advanced)
- [Experiment with ezknitr](#experiment)
- [spin() vs knit()](#spin-vs-knit)
- [Using rmarkdown::render()](#using-render)

<h2 id="installation">Installation</h2>

`ezknitr` is available through both CRAN and GitHub.

To install the CRAN version:

```
install.packages("ezknitr")
```

To install the latest developmental version from GitHub:

```
install.packages("devtools")
devtools::install_github("ropensci/ezknitr")
```

<h2 id="overview">Overview</h2>

If you have a very simple project with a flat directory structure, then `knitr` works great. But even something as simple as trying to knit a document that reads a file from a different directory or placing the output rendered files in a different folder cannot be easily done with `knitr`.

`ezknitr` improves basic `knitr` functionality in a few ways. You get to decide:

- What the working directory of the source file is
    - Default is your current working directory, which often makes more sense than the `knitr` assumption that the working directory is wherever the input file is
- Where the output files will go
    - With `knitr`, all the rendered output files will be generated in the folder you're currently in
- Where the figures generated in the markdown document will be stored
    - `knitr` makes it cumbersome to change this directory
- Any parameters to pass to the source file
    - Useful if you want to run an identical source file multiple times with different parameters

<h2 id="motivation-simple-usecase">Motivation & simple use case</h2>

Assume you have an Rmarkdown file that reads a data file and produces a short report while also generating a figure. Native `knit()` (or `spin()` if you're starting with an R script instead of an Rmd file) works great if you have a flat directory structure like this:

```
- project/
  |- input.csv
  |- report.Rmd
```

But what happens if you have a slightly more complex structure? In a real project, you rarely have everything just lying around in the same folder.  Here is an example of a more realistic initial directory structure (assume the working directory/project root is set to `project/`):

```
- project/
  |- analysis/
     |- report.Rmd
  |- data/
     |- input.csv
```

Now if you want `knitr` to work, you'd have to ensure the path to `input.csv` is relative to the `analysis/` directory because that's where the Rmd file is. This is counter-intuitive because most people expect to create paths relative to the *working directory/project root* (`project/` in this case), but `knitr` will use the `analysis/` folder as the working directory. Any code reading the input file needs to use `../data/input.csv` instead of `data/input.csv`.

Other than being confusing, it also means that you cannot naively run the Rmd code chunks manually because when you run the code in the console, your working directory is not set to what `knitr` will use as the working directory. More specifically, if you try to run the command that reads the input file, your console will look in `project/../data/input.csv` (which doesn't exist).

A similar problem arises when you want to create files in your report: `knitr` will create the files relative to where the Rmd file is, rather than relative to the project root.

Another problem with the flat directory structure is that you may want to control where the resulting reports get generated. `knitr` will create all the outputs in your working directory, and as far as I know there is no way to control that. 

`ezknitr` addresses these issues, and more. It provides wrappers to `knit()` and `spin()` that allow you to set the working directory for the input file, and also uses a more sensible default working directory: the current working directory. `ezknitr` also lets you decide where the output files and output figures will be generated, and uses a better default path for the output files: the directory containing the input file.

Assuming your working directory is currently set to the `project/` directory, you could use the following `ezknitr` command to do what you want:

```
library(ezknitr)
ezknit(file = "analysis/report.Rmd", out_dir = "reports", fig_dir = "myfigs")
```

```
- project/
  |- analysis/
     |- report.Rmd
  |- data/
     |- input.csv
  |- reports/
     |- myfigs/
        |- fig1.png
     |- report.md
     |- report.HTML
```

We didn't explicitly have to set the working direcory, but you can use the `wd` argument if you do require a different directory (for example, if you are running this from some build script or from any arbitrary directory). After running `ezknit()`, you can run `open_output_dir()` to open the output directory in your file browser if you want to easily see the resulting report. Getting a similar directory structure with `knitr` is not simple, but with `ezknitr` it's trivial.

Note that `ezknitr` produces both a markdown and an HTML file for each report (you can choose to discard them with the `keep_md` and `keep_html` arguments).

<h2 id="usecase-advanced">Use case: using one script to analyze multiple datasets</h2>

As an example of a more complex realistic scenario where `ezknitr` would be useful, imagine having multiple analysis scripts, with each one needing to be run on multiple datasets. Being the organizer scientist that you are, you want to be able to run each analysis on each dataset, and keep the results neatly organized. I personally was involved in a few projects requiring exactly this, and `ezknitr` was in fact born for solving this exact issue. Assume you have the following files in your project:

```
- project/
  |- analysis/
     |- calculate.Rmd
     |- explore.Rmd
  |- data/
     |- human.dat
     |- mouse.dat
```

We can easily use `ezknitr` to run any of the analysis Rmarkdowns on any of the datasets and assign the results to a unique output. Let's assume that each analysis script expects there to be a variable named `DATASET_NAME` inside the script that tells the script what data to operate on. The following `ezknitr` code illustrates how to achieve the desired output.

```
library(ezknitr)
ezknit(file = "analysis/explore.Rmd", out_dir = "reports/human",
        params = list("DATASET_NAME" = "human.dat"), keep_html = FALSE)
ezknit(file = "analysis/explore.Rmd", out_dir = "reports/mouse",
        params = list("DATASET_NAME" = "mouse.dat"), keep_html = FALSE)
ezknit(file = "analysis/calculate.Rmd", out_dir = "reports/mouse",
        params = list("DATASET_NAME" = "mouse.dat"), keep_html = FALSE)
```

```
- project/
  |- analysis/
     |- calculate.Rmd
     |- explore.Rmd
  |- data/
     |- human.dat
     |- mouse.dat
  |- reports/
     |- human/
        |- explore.md
     |- mouse/
        |- calculate.md
        |- explore.md
```

Note that this example uses the `params = list()` argument, which lets you pass variables to the input Rmarkdown.  In this case, I use it to tell the Rmarkdown what dataset to use, and the Rmarkdown assumes a `DATASET_NAME` variable exists. This of course means that the analysis script has an external dependency by having a variable that is not defined inside of it. You can use the `set_default_params()` function inside the Rmarkdown to ensure the variable uses a default value if none was provided.

Also note that differentiating the species in the output could also have been done using the `out_suffix` argument instead of the `out_dir` argument. For example, using `out_suffix = "human"` would have resulted in an ouput file named `explore-human.md`.

<h2 id="experiment">Experiment with ezknitr</h2>

After installing and loading the package (`library(ezknitr)`), you can experiment with `ezknitr` using the `setup_ezknit_test()` or `setup_ezspin_test()` functions to see their benefits. See `?setup_ezknit_test` for more information.

<h2 id="spin-vs-knit">spin() vs knit()</h2>

`knit()` is the most popular and well-known function from `knitr`. It lets you create a markdown document from an Rmarkdown file. You can learn more about `knit()` on [the official knitr page](http://yihui.name/knitr).

`spin()` is similar, but starts one step further back: it takes an R script as input, creates an Rmarkdown document from the R script, and then proceeds to create a markdown document from it. `spin()` can be useful in situations where you develop a large R script and want to be able to produce reports from it directly instead of having to copy chunks into a separate Rmarkdown file. You can read more about why I like `spin()` in the blog post ["knitr's best hidden gem: spin"](http://deanattali.com/2015/03/24/knitrs-best-hidden-gem-spin/).

<h2 id="using-render">Using rmarkdown::render()</h2>

When the core of this package was developed, none of its functionality was supported in any way by either `knitr` or `rmarkdown`. Over time, `rmarkdown::render()` got some new features that are very similar to features of `ezknitr`.  Native support for parameters inside Rmarkdown files using YAML is a big feature which makes the use of `set_default_params()` and the `params` argument of `ezknitr` less important.  However, the core problem that `ezknitr` wants to solve is the working directory issue, and this issue has yet to be addressed by `rmarkdown` or `knitr`, which makes `ezknitr` still useful.

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ezknitr.R
\docType{package}
\name{ezknitr}
\alias{ezknitr}
\alias{ezknitr-package}
\title{ezknitr}
\description{
Avoid the Typical Working Directory Pain When Using 'knitr'
}
\details{
\code{ezknitr} is an extension of \code{knitr} that adds flexibility in several
ways. One common source of frustration with \code{knitr} is that it assumes
the directory where the source file lives should be the working directory,
which is often not true. \code{ezknitr} addresses this problem by giving you
complete control over where all the inputs and outputs are, and adds several
other convenient features to make rendering markdown/HTML documents easier.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_test.R
\name{setup_ezspin_test}
\alias{setup_ezknit_test}
\alias{setup_ezspin_test}
\alias{setup_test}
\title{Set up a test directory to experiment with \code{ezspin} or \code{ezknit}}
\usage{
setup_ezspin_test()

setup_ezknit_test()
}
\value{
The path to the test directory.
}
\description{
Create a few directories that try to mimic a real
data-analysis project structure, and add a data file and a simple R script
(for \code{ezspin}) or Rmarkdown file (for \code{ezknit}).\cr\cr
After setting up these files and directories, you can run \code{ezknitr}
commands and their equivalent \code{knitr} commands to compare and see the
benefits of using \code{ezknitr}.\cr\cr
More specific instructions on how to interact with this test directory will
be printed to the console.
}
\examples{
\dontrun{
library(ezknitr)

# setup the test directory structures and run naive spin
setup_ezspin_test()
knitr::spin("ezknitr_test/R/ezspin_test.R")
file.remove(c("ezspin_test.md", "ezspin_test.html"))

# setup the test directory structures and run simple ezspin
setup_ezspin_test()
ezspin("R/ezspin_test.R", wd = "ezknitr_test")

# setup the test directory structures and run ezspin with more parameters
tmp <- setup_ezspin_test()
ezspin("R/ezspin_test.R", wd = "ezknitr_test",
        out_dir = "output", fig_dir = "coolplots")
unlink(tmp, recursive = TRUE, force = TRUE)
}
}
\seealso{
\code{\link[ezknitr]{ezspin}}
\code{\link[knitr]{spin}}
\code{\link[ezknitr]{ezknit}}
\code{\link[knitr]{knit}}
\code{\link[ezknitr]{open_output_dir}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/open_output_dir.R
\name{open_output_dir}
\alias{open_output_dir}
\title{Open the directory containing the output from the last ezknitr command}
\usage{
open_output_dir()
}
\description{
Call this function after running \link[ezknitr]{ezspin} or
\link[ezknitr]{ezknit} to open the resulting output directory in your
file browser. This is simply a convenience function so that if you want to
see the results you don't need to navigate to the appropriate folder manually.
}
\examples{
\dontrun{
library(ezknitr)
setup_ezspin_test()
ezspin("R/ezspin_test.R", wd = "ezknitr_test")
open_output_dir()
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{ezknitr_core}
\alias{ezknit}
\alias{ezknitr_core}
\alias{ezspin}
\title{Knit Rmd or spin R files without the typical pain of working directories}
\usage{
ezspin(file, wd, out_dir, fig_dir, out_suffix, params = list(),
  verbose = FALSE, chunk_opts = list(tidy = FALSE), keep_rmd = FALSE,
  keep_md = TRUE, keep_html = TRUE, move_intermediate_file = TRUE, ...)

ezknit(file, wd, out_dir, fig_dir, out_suffix, params = list(),
  verbose = FALSE, chunk_opts = list(tidy = FALSE), keep_md = TRUE,
  keep_html = TRUE)
}
\arguments{
\item{file}{The path to the input file (.Rmd file if using \code{ezknit} or 
.R script if using \code{ezspin}). If \code{wd} is provided, then this path is 
relative to \code{wd}.}

\item{wd}{The working directory to be used in the Rmd/R script. Defaults to
the current working directory (note that this is not the same behaviour as
\code{knitr}). See the 'Detailed Arguments' section for more details.}

\item{out_dir}{The output directory for the rendered markdown or HTML files
(if \code{wd} is provided, then this path is relative to \code{wd}).
Defaults to the directory containing the input file.}

\item{fig_dir}{The name (or path) of the directory where figures should
be generated. See the 'Detailed Arguments' section for more details.}

\item{out_suffix}{A suffix to add to the output files. Can be used to
differentiate outputs from runs with different parameters. The name of the
output files is the name of the input file appended by \code{out_suffix},
separated by a dash.}

\item{params}{A named list of parameters to be passed to use in the input
Rmd/R file. For example, if the script to execute assumes that there is a
variable named \code{DATASET_NAME}, then you can use
\code{params = list('DATASET_NAME' = 'oct10dat')}.}

\item{verbose}{If TRUE, then show the progress of knitting the document.}

\item{chunk_opts}{List of knitr chunk options to use. See 
\code{?knitr::opts_chunk} for a list of available chunk options.}

\item{keep_rmd, keep_md}{Should intermediate \code{Rmd} or \code{md} files be
kept (\code{TRUE}) or deleted (\code{FALSE})?}

\item{keep_html}{Should the final \code{html} file be kept (\code{TRUE})
or deleted (\code{FALSE})?}

\item{move_intermediate_file}{Should the intermediate \code{Rmd} file be
moved to the output folder (\code{TRUE}) or stay in the same folder as
the source \code{R} file (\code{FALSE})?}

\item{...}{Any extra parameters that should be passed to \code{knitr::spin}.}
}
\value{
The path to the output directory (invisibly).
}
\description{
\code{ezknitr} is an extension of \code{knitr} that adds flexibility in several
ways. One common source of frustration with \code{knitr} is that it assumes
the directory where the source file lives should be the working directory,
which is often not true. \code{ezknitr} addresses this problem by giving you
complete control over where all the inputs and outputs are, and adds several
other convenient features. The two main functions are \code{ezknit} and 
\code{ezspin}, which are wrappers around \code{knitr}'s \code{knit} and
\code{spin}, used to make rendering markdown/HTML documents easier.
}
\details{
If you have a very simple project with a flat directory structure, then
\code{knitr} works great. But even something as simple as trying to knit a 
document that reads a file from a different directory or placing the output 
rendered files in a different folder cannot be easily done with \code{knitr}.

\code{ezknitr} improves basic \code{knitr} functionality in a few ways.
You get to decide:
\itemize{
  \item What the working directory of the source file is
  \item Where the output files will go
  \item Where the figures used in the markdown will go
  \item Any parameters to pass to the source file
}
}
\section{Detailed Arguments}{

All paths given in the arguments can be either absolute or relative.

The \code{wd} argument is very important and is set to the current working
directory by default. The path of the input file and the path of the output
directory are both relative to \code{wd} (unless they are absolute paths).
Moreover, any code in the R script that reads or writes files will use
\code{wd} as the working directory.

The \code{fig_dir} argument is relative to the output directory, since the
figures accompanying a markdown file should be placed in the same
directory. It is recommended to either leave \code{fig_dir} as default or
set it to a different name but not to a different directory. Because of the
way \code{knitr} works, there are a few known minor issues if \code{fig_dir}
is set to a different directory.
}

\section{Difference between ezknit and ezspin}{

\code{ezknit} is a wrapper around \code{knitr::knit} while \code{ezspin}
is a wrapper around \code{ezspin}. The two functions are very similar.
\code{knit} is the more popular and well-known function. It is used
 to render a markdown/HTML document from an Rmarkdown source. 
\code{spin} takes an R script as its input, produces an
Rmarkdown document from the R script, and then calls \code{knit} on it.
}
\examples{
\dontrun{
   tmp <- setup_ezknit_test()
   ezknit("R/ezknit_test.Rmd", wd = "ezknitr_test")
   ezknit("R/ezknit_test.Rmd", wd = "ezknitr_test",
          out_dir = "output", fig_dir = "coolplots",
          params = list(numPoints = 50))
   open_output_dir()
   unlink(tmp, recursive = TRUE, force = TRUE)
 
   tmp <- setup_ezspin_test()
   ezspin("R/ezspin_test.R", wd = "ezknitr_test")
   ezspin("R/ezspin_test.R", wd = "ezknitr_test",
          out_dir = "output", fig_dir = "coolplots",
          params = list(numPoints = 50), keep_rmd = TRUE)
   open_output_dir()
   unlink(tmp, recursive = TRUE, force = TRUE)
}
}
\seealso{
\code{\link[ezknitr]{open_output_dir}}
\code{\link[ezknitr]{setup_ezknit_test}}
\code{\link[ezknitr]{setup_ezspin_test}}
\code{\link[ezknitr]{set_default_params}}
\code{\link[knitr]{knit}}
\code{\link[knitr]{spin}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_default_params.R
\name{set_default_params}
\alias{set_default_params}
\title{Set default parameters}
\usage{
set_default_params(params)
}
\arguments{
\item{params}{List of parameters.}
}
\description{
Create variables with the given values only if these variables do not currently
exist.
}
\details{
Sometimes it may be useful to define a variable only it hasn't been defined yet.
One example where this can be useful is when you have an Rmd script that 
uses some variables and you want to be able to use custom values for these
variables, but also give them a default value in the script in case they are
not set beforehand.
}
\examples{
exists("foo")
exists("bar")
foo <- 5
set_default_params(list(foo = 10, bar = 20))
print(foo)
print(bar)
}

