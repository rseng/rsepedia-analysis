---
title: 'Suppdata: Downloading Supplementary Data from Published Manuscripts'
authors:
- affiliation: 1
  name: William D Pearse
  orcid: 0000-0002-6241-3164
- affiliation: 2
  name: Scott A Chamberlain
  orcid: 0000-0003-1444-9135
date: "2 May 2017"
bibliography: paper.bib
tags:
- supplementary materials
- supplemental information
- open data
- meta-analysis
- DOI
affiliations:
- index: 1
  name: Department of Biology & Ecology Center, Utah State University, Logan, Utah,
    USA
- index: 2
  name: rOpenSci
---

# Summary

`suppdata` is an R [@R2018] package to provide easy, reproducible
access to supplemental materials within R. Thus `suppdata` facilitates
open, reproducible research workflows: scientists re-analyzing
published datasets can work with them as easily as if they were stored
on their own computer, and others can track their analysis workflow
painlessly.

For example, imagine you were conducting an analysis of the evolution
of body mass in mammals. Without `suppdata`, such an analysis would
require manually downloading body mass and phylogenetic data from
published manuscripts. This is time-consuming, difficult (if not
impossible) to make truly reproducible without re-distributing the
data, and hard to follow. With `suppdata`, such an analysis is
straightforward, reproducible, and the sources of the data
[@Fritz2009, @Jones2009] are clear because their DOIs are embedded
within the code:

```{R}
# Load phylogenetics packages
library(ape)
library(caper)
library(phytools)

# Load suppdata
library(suppdata)

# Load two published datasets
tree <- read.nexus(suppdata("10.1111/j.1461-0248.2009.01307.x", 1))[[1]]
traits <- read.delim(suppdata(
		"E090-184", "PanTHERIA_1-0_WR05_Aug2008.txt",
		"esa_archives"))

# Merge datasets
traits <- with(traits, data.frame(body.mass = log10(X5.1_AdultBodyMass_g),
             	       species=gsub(" ","_",MSW05_Binomial)))
c.data <- comparative.data(tree, traits, species)

# Calculate phylogenetic signal
phylosig(c.data$phy, c.data$data$body.mass)
```

The above example makes use of code from the packages `ape`
[@Paradis2004], `caper` [@Orme2013], and `phytools` [@Revell2012].

As `suppdata` was, originally, part of `fulltext` [@Chamberlain2018],
it is already being used in a number of research projects. One such
project is `natdb`, a package that builds a database of functional
traits from published sources. The software is currently available on
GitHub (https://github.com/ropensci/suppdata), and we plan to
distribute it through ROpenSci and CRAN.

# References<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/suppdata)](https://cran.r-project.org/package=suppdata)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://api.travis-ci.org/ropensci/suppdata.svg)](https://travis-ci.org/ropensci/suppdata)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/willpearse/suppdata?branch=master&svg=true)](https://ci.appveyor.com/project/willpearse/suppdata)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.00721/status.svg)](https://doi.org/10.21105/joss.00721)
[![](https://badges.ropensci.org/195_status.svg)](https://github.com/ropensci/onboarding/issues/195)
[![codecov](https://codecov.io/gh/ropensci/suppdata/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/suppdata)
<!-- badges: end -->

# Loading SUPPlementary DATA into R

William D. Pearse, Daniel Nuest, and Scott Chamberlain

## Overview

The aim of this package is to aid downloading data from published
papers. To download the supplementary data from a PLoS paper, for
example, you would simply type:

```{R}
library(suppdata)
suppdata("10.1371/journal.pone.0127900", 1)
```

...and this would download the first supplementary information (SI) from the paper.

This sort of thing is very useful if you're doing meta-analyses, or
just want to make sure that you know where all your data came from and
want a completely reproducible "audit trail" of what you've done.
It uses [`rcrossref`](https://cran.r-project.org/package=rcrossref) to lookup which journal the article is in.

## How to install and load the package
The version on CRAN is the most stable version. You can install and
load it like this:

```{R}
install.packages("suppdata")
library(suppdata)
```

If you want to load the development version, which probably contains
more features but is not always guaranteed to work, load the `master`
branch from this repository like this:

```{R}
library(devtools)
install_github("ropensci/suppdata")
library(suppdata)
```

This package depends on the packages `httr`, `xml2`,
`jsonlite`, and `rcrossref`.

## Supported publishers and repositories

- [bioRxiv](https://www.biorxiv.org/) (`biorxiv`)
- [Copernicus Publications](https://publications.copernicus.org/) (`copernicus`)
- [DRYAD](https://datadryad.org/) (`dryad`)
- [Ecological Society of Ameria - Ecological Archives](http://esapubs.org/archive/) (`esa_archives` and `esa_data_archives`)
- [Europe PMC](https://europepmc.org/) (`epmc`, multiple publishers from life-sciences upported including BMJ Journals, eLife, F1000Research, Wellcome Open Research, Gates Open Research)
- [figshare](https://figshare.com/) (`figshare`)
- [Journal of Statistical Software](https://www.jstatsoft.org/) (`jstatsoft`)
- [MDPI](https://www.mdpi.com/) (`mdpi`)
- [PeerJ](https://peerj.com/) (`peerj`)
- [PLOS | Public Library of Science](https://www.plos.org/) (`plos`)
- [Proceedings of the royal society Biology (RSBP)](https://rspb.royalsocietypublishing.org/) (`proceedings`)
- [Science](https://www.sciencemag.org/) (`science`)
- [Wiley](https://onlinelibrary.wiley.com/) (`wiley`)

See a list of potential sources at [#2](https://github.com/ropensci/suppdata/issues/2) - requests welcome!

## Contributing

[For more details on how to contribute to the package, check out the
guide in `CONTRIBUTING.MD`](CONTRIBUTING.md).

## A more detailed set of motivations for `suppdata`

`suppdata` is an R package to provide easy, reproducible
access to supplemental materials within R. Thus `suppdata` facilitates
open, reproducible research workflows: scientists re-analyzing
published datasets can work with them as easily as if they were stored
on their own computer, and others can track their analysis workflow
painlessly.

For example, imagine you were conducting an analysis of the evolution
of body mass in mammals. Without `suppdata`, such an analysis would
require manually downloading body mass and phylogenetic data from
published manuscripts. This is time-consuming, difficult (if not
impossible) to make truly reproducible without re-distributing the
data, and hard to follow. With `suppdata`, such an analysis is
straightforward, reproducible, and the sources of the data are clear
because their DOIs are embedded within the code:

```{R}
# Load phylogenetics packages
library(ape)
library(caper)
library(phytools)

# Load suppdata
library(suppdata)

# Load two published datasets
tree <- read.nexus(suppdata("10.1111/j.1461-0248.2009.01307.x", 1))[[1]]
traits <- read.delim(suppdata("E090-184", "PanTHERIA_1-0_WR05_Aug2008.txt", "esa_archives"))

# Merge datasets
traits <- with(traits, data.frame(body.mass = log10(X5.1_AdultBodyMass_g), species=gsub(" ","_",MSW05_Binomial)))
c.data <- comparative.data(tree, traits, species)

# Calculate phylogenetic signal
phylosig(c.data$phy, c.data$data$body.mass)
```

## A guided walk through `suppdata`

The aim of `suppdata` is to make it as easy as possible for you to write reproducible analysis scripts that make use of published data. So let's start with that first, simplest case: how to make use of published data in an analysis.

### Learning by example
Below is an example of an analysis run using `suppdata`. Read through it first, and then we'll go through what all the parts mean.

```{R}
# Load phylogenetics packages
library(ape)
library(caper)
library(phytools)

###############################
# LOAD TWO PUBLISHED DATASETS #
#       USING SUPPDATA        #
###############################
library(suppdata)
tree <- read.nexus(suppdata("10.1111/j.1461-0248.2009.01307.x", 1))[[1]]
traits <- read.delim(suppdata("E090-184", "PanTHERIA_1-0_WR05_Aug2008.txt", "esa_archives"))

# Merge datasets
traits <- with(traits, data.frame(body.mass = log10(X5.1_AdultBodyMass_g), species=gsub(" ","_",MSW05_Binomial)))
c.data <- comparative.data(tree, traits, species)

# Calculate phylogenetic signal
phylosig(c.data$phy, c.data$data$body.mass)
```

This short script loads some `R` packages focused on modelling the evolution of species' traits, then it gets to the "good stuff": using `suppdata`. First, we load the `suppdata` package using `library(suppdata)`. The next line uses a function called `read.nexus`, which loads something called a phylogeny (you might be familiar with this if you're a biologist). Normally, this function would take the location of a file on our hard-drive as a single argument, but now we're giving it the output from a call to the `suppdata` function.

`suppdata` is going to the website of the article whose DOI is _10.1111/j.1461-0248.2009.01307.x_ ([it's this paper by Fritz et al.](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1461-0248.2009.01307.x)), and then taking the first (`1`) supplement from that article. It saves that to a temporary location on your hard-drive, and then gives that location to `read.nexus`. _This works with any function that expects a file at a location on your hard-drive_. What particularly neat is that `suppdata` remembers that it has downloaded that file already (see below for more details), such that you only have to download something once per `R` session.

The second call to `suppdata`, which makes use of `read.delim`, shows two of the potential complexities of `suppdata`. First of all, because some journal publishers store their supplementary materials using numbers and others using specific file-names, `suppdata` takes either a number (like in the first example), or a name (like in the second example) depending on the journal publisher you're taking data from. If you look in the help file for `suppdata`, there is a table outlining those options. Sorry, you've just got to read up on it :-( Secondly, if you're an ecologist you might be familiar with the Ecological Society of America's data archives section. While they've moved over to a new way of storing data more recently, if you're hoping to load an older dataset from that journal you need to give the ESA data archive reference and specify that you're downloading from ESA (as in this example). If you're not an ecologist, don't worry about it, as this doesn't apply to you :D

That's it! You now know all you need to in order to use `suppdata`! The rest of the lines of code merge these datasets together, and then calculate something called _phylogenetic signal_ in these datasets. If you're an evolutionary biologist, those lines might be interesting to you. If you're not, then don't worry about them.

### Caching and saving to a specific directory

Sometimes, you will want to use `suppdata` to build a store of files on your hard-drive. If so, you should know that `suppdata` takes three optional arguments: `cache`, `dir`, and `save.name`. If you specify `cache=FALSE`, it will turn off `suppdata`'s caching of files: this will force it to download your data again. This is mostly useful if you somehow make `suppdata` foul itself up (maybe you hit control-c or stop during a download) and so `suppdata` has only half-downloaded a file, and so thinks it's cached something when it hasn't. If you get an error when using `suppdata`, this is a good thing to try setting.

`dir` specifies a directory where `suppdata` should store files, and `save.name` specifies the name that the file should be saved under when saved. This is useful if you want to make a folder on your computer that contains certain files you use a lot: `suppdata` will cache from this folder if you tell it to, and so you can build up a reproducible selection of data to use inbetween `R` sessions.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# suppdata 1.1-7
- Changed FigShare download URL to match recent change online

# suppdata 1.1-6
- Hotfixes for CRAN

# suppdata 1.1-5

- Added fall-back for zip downloads on Windows (thanks Alban Sagouis)
- New errors for moved download URLs
- Moved all help file examples back to \dontrun (at CRAN's request)
- Removed ESA data archives; site down at time of required submission

# suppdata 1.1-4

- FigShare website change --> no longer supports named SI

# suppdata 1.1-3

- Changes to available data on EPMC for a demonstration package
  required a change to the code on CRAN to pass its checks

# suppdata 1.1-0

- More thorough and numerous tests
- New journals added
- New co-author - Daniel Nüst!

# suppdata 1.0.0

- Initial release on CRAN
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at will.pearse@usu.edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [https://contributor-covenant.org/version/1/4][version]

[homepage]: https://contributor-covenant.org
[version]: https://contributor-covenant.org/version/1/4/
# Overview

Thank you for your interest in contributing to `suppdata`! The most important contribution you can make to `suppdata` is to add code to download data from another publisher's journals. There are five steps you have to go through to do that; I go through them in detail below, but briefly, they are:

1. Write a download function for a publisher
2. Modify the DOI lookup function so it knows about your download function
3. Modify the documentation for `suppdata` so the user knows what you've done
4. Write a brief unit test checking your function works
5. Submit a pull request

If you want to make a small change to `suppdata` (i.e., changing <= 5 lines of code) fork the repo, make the change, and then make a pull request with the suggestion. If you want to make a more sweeping change (i.e., > 5 lines of code) then _before writing any code_ make an issue and discuss it with @willpearse. The purpose of this is to make maximal use of everyone's time: small code changes are better off "just done" and then we can talk about it; larger changes require discussion before implementation. You're quite welcome to do whatever you wish with the code (within the boundaries of the license, of course), but please be aware that the maintainers of the package are not obligate to accept all pull requests. Of course, the ROpenSci maintainer rules apply, so we'll always be polite and we'll always let you know why we make any decision! :D

When making version changes, please follow the standards set by CRAN, so the next version after "1.2-9" would be "1.2-10". Package versions numbers are not decimal, so something like "1.2.9.5" won't pass CRAN's checks ([see the R extension guide](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#The-DESCRIPTION-file))

## (1) Write a download function for a publisher

All `suppdata` download functions start off with, at a minimum, something like:

```{R}
.suppdata.nameofpublisher <- function(doi, si, save.name=NA, dir=NA, cache=TRUE, ...)
```

...where `nameofpublisher` is replaced with the name of the publisher you're loading data from. This should then be followed with something like:

```{R}
#Argument handling
if(!is.numeric(si))
    stop("nameofpublisher download requires numeric SI info")
dir <- .tmpdir(dir)
save.name <- .save.name(doi, save.name, si)
```

Note the name of the section (`#Argument handling`), and how we're making sure the user gives us a supplement number if we need that, or using `is.character` if we're expecting some `character SI info`. We have a few internal functions that should make your life easier (...well, they should...); `.tmp.dir` will make a temporary directory to save files out to for you, and `.save.name` will generate a sensible name to save files out to as well. It's quite important that you use those functions, since they make the `cache` behaviour of the package work.

Next comes the hard work, where you get the data out of the publisher's website. [This is a lot of regular expression kung-foo](https://xkcd.com/208/); you may find the functions `.grep.url`, `.grep.text`, `.url.redir`, and all the other functions in [`utils.R`](https://github.com/willpearse/suppdata/blob/master/R/utils.R) useful. Please do have a poke around in there, and feel free to add any functions you think are missing (using the `.name.function.` convention for non-user-facing functions). If there's something you think would be useful to have, and you don't know how to write it, then just make an issue, tag me (@willpearse), and I'll see what I can do.

Finally, you need to return to the user a location of the file they want. Something like:

```{R}
destination <- file.path(dir, save.name)
return(.download(url, dir, save.name, cache))
```

...should suffice. Notice how we're using `R`'s `file.path` to make a sensible path on all distributions, and we're using the internal `.download` to download the file, and so guaranteeing that we'll obey all the `cache` instructions etc. We're also ensuring that the user will get sensible filename information ("oooh, this looks like a .csv file") as an attribute by using `.download`.

Save your function in [`journals.R`](https://github.com/willpearse/suppdata/blob/master/R/journals.R). There are plenty of examples in there if you get stuck. There is also a list of functions to be written sitting in the issues section on GitHub.

There is a 'hit list' of publishers that it would be great to write wrappers for up in the issues page - [click this link to see it](https://github.com/willpearse/suppdata/issues/2).

## (2) Modify the DOI lookup function so it knows about your download function

We're nearly there now, I promise! `suppdata` uses `rcrossref` to look up articles' publishers, so to hook your download function into the package you're going to have to figure out your journal's code. Take a paper's DOI that you know works, and run the following on it (replacing the DOI below with the DOI you're checking):

```{R}
library(rcrossref)
cr_works("10.1111/j.1461-0248.2009.01307.x")$data$member
```

...that number is your publisher's code. Modify the first `switch` statement in `suppdata.func` in [`utils.R`](https://github.com/willpearse/suppdata/blob/master/R/utils.R) to add your journal's number (as a character string) and then match it with the name of your download function. If that sounds complicated, once you open the function it will become obvious.

Next, modify the second (and last) `switch` in the same function to work if your publisher is known by name. This should match onto your function's name. So, for example, if your publishing company were called _Pearse Publishing_, and you'd called your function `.suppdata.pearse`, then you would add an entry like:

```{R}
"pearse" = .suppdata.pearse
```

Please remember that `case` statements are separated by commas; add a comma to the previous entry to keep the code syntactically correct!

## (3) Modify the documentation for `suppdata` so the user knows what you've done

Modify the `roxygen` documentation [here](https://github.com/willpearse/suppdata/blob/master/R/suppdata.R#L17) and [here](https://github.com/willpearse/suppdata/blob/master/man-roxygen/suppdata.R#L27) to give the user information about what your function expects (`numeric` vs. `character` SI information). Re-build the documentation when you're done by running something like:

```{R}
library(roxygen2)
roxygenize("path/to/suppdata")
```

...if you're an RStudio person there's a button for this in the "Build" tab ("More" > "Document"). If you're an emacs person like me, there are several and you probably have a strong opinion about which is best :p

## (4) Write a brief unit test checking your function works

Add tests to [`tests/testthat/`](https://github.com/willpearse/suppdata/blob/master/tests/testthat/) with a file of the format `test-<name of publisher>.R` (see existing tests to get started) to give the maintainers and the continuous integration services something to check that the publisher works.

## (5) Add publisher to list in README

Add the newly supported publisher to the alphabetically ordered list in `README.md`.

## (6) Submit a pull request

Commit your changes, then make a pull request to the `master` branch of `suppdata`. If you need help figuring out how to do that, [take a look at this website](https://help.github.com/articles/creating-a-pull-request/).

## (optional) Bask in glory

Thank you for helping make `suppdata` better!!!

## (side-note) History of the package

I always think a little history is useful when contributing to a package, so let me tell you how the code here came to be. Originally, this package was called `grabr` (this repo is the same repo, the name alone changed), which was then merged into [`fulltext`](https://github.com/ropensci/fulltext/). At that time, the structure was changed to match that of `fulltext`, but this code was then pulled back _out_ of `fulltext`. So if you ever find yourself wondering "why does it do that?", the answer is probably "because `fulltext` does".
