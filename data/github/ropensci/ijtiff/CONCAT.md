
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ijtiff <img src="man/figures/logo.png" align="right" height=140/>

This is a general purpose TIFF I/O utility for R. The [`tiff`
package](https://cran.r-project.org/package=tiff) already exists for
this purpose but `ijtiff` adds some functionality and overcomes some
bugs therein.

  - `ijtiff` can write TIFF files whose pixel values are real
    (floating-point) numbers; `tiff` cannot.
  - `ijtiff` can read and write *text images*; `tiff` cannot.
  - `tiff` struggles to interpret channel information and gives cryptic
    errors when reading TIFF files written by the *ImageJ* software;
    `ijtiff` works smoothly with these images.

The github repo of `ijtiff` is at <https://github.com/ropensci/ijtiff>.

## Installation

You can install the released version of `ijtiff` from
[CRAN](https://CRAN.R-project.org/package=ijtiff) with:

``` r
install.packages("ijtiff")
```

You can install the released version of `ijtiff` from
[GitHub](https://github.com/ropensci/ijtiff) with:

``` r
devtools::install_github("ropensci/ijtiff")
```

## How to use `ijtiff`

The [Reading and Writing
Images](https://docs.ropensci.org/ijtiff/articles/reading-and-writing-images.html)
article is probably all you need to know.

## More about `ijtiff`

  - [Text
    Images](https://docs.ropensci.org/ijtiff/articles/text-images.html)
    tells you more about what *text images* are and why you might ever
    use them.
  - [The *ImageJ*
    Problem](https://docs.ropensci.org/ijtiff/articles/the-imagej-problem.html)
    explains the problem that `tiff` has when reading TIFF files written
    by *ImageJ* and how `ijtiff` fixes this problem.
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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# ijtiff <img src="man/figures/logo.png" height="140" align="right">

[![R-CMD-check](https://github.com/ropensci/ijtiff/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ijtiff/actions)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ijtiff?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/ijtiff)
[![codecov](https://codecov.io/gh/ropensci/ijtiff/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/ijtiff)

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ijtiff)](https://cran.r-project.org/package=ijtiff)
![RStudio CRAN
downloads](http://cranlogs.r-pkg.org/badges/grand-total/ijtiff)
![RStudio CRAN monthly
downloads](http://cranlogs.r-pkg.org/badges/ijtiff)

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00633/status.svg)](https://doi.org/10.21105/joss.00633)

## Introduction

This is a general purpose TIFF I/O utility for R. The [`tiff`
package](https://cran.r-project.org/package=tiff) already exists for
this purpose but `ijtiff` adds some functionality and overcomes some
bugs therein.

-   `ijtiff` can write TIFF files whose pixel values are real
    (floating-point) numbers; `tiff` cannot.
-   `ijtiff` can read and write *text images*; `tiff` cannot.
-   `tiff` struggles to interpret channel information and gives cryptic
    errors when reading TIFF files written by the *ImageJ* software;
    `ijtiff` works smoothly with these images.

To learn about `ijtiff` and how to use it, visit the package website at
<https://docs.ropensci.org/ijtiff/>.

## Installation

### `libtiff`

`ijtiff` requires you to have the `libtiff` C library installed. To
install `libtiff`:

-   On **Debian Linux**, try `sudo apt-get install libtiff5-dev`, or if
    that fails, try  
    `sudo apt-get install libtiff4-dev`.
-   On **Fedora Linux**, try `sudo yum install libtiff5-dev`, or if that
    doesn’t work, try  
    `sudo yum install libtiff4-dev`.
-   On **Mac**, you need [Homebrew](https://brew.sh/). Then in the
    terminal, run `brew install libtiff`.
-   On **Windows**, no setup is required 😄.

### Installing the release version of the `ijtiff` R package

You can install `ijtiff` from CRAN (recommended) with:

``` r
install.packages("ijtiff")
```

### Installing the development version of the `ijtiff` R package

You can install the development version from GitHub with:

``` r
devtools::install_github("ropensci/ijtiff")
```

## Acknowledgement

This package uses a lot of code from the original `tiff` package by
Simon Urbanek.

## Contribution

Contributions to this package are welcome. The preferred method of
contribution is through a github pull request. Feel free to contact me
by creating an issue. Please note that this project is released with a
[Contributor Code of
Conduct](https://github.com/ropensci/ijtiff/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# `ijtiff` 2.2.7

## BUG FIXES
* Fix for new libtiff using C99's `<stdint.h>`.


# `ijtiff` 2.2.6

## BUG FIXES
* Suppress unhelpful warnings during configure when `pkg-config` doesn't find info for libtiff.
* Remove `LazyData` from `DESCRIPTION` (was causing CRAN note).


# `ijtiff` 2.2.5

## BUG FIXES
* Typo fix for `configure`. At one point there was a call of `pkg-configs` instead of `pkg-config`.
* Also now all compile flags from `pkg-config --libs` _and_ `pkg-config --libs --static` are used every time.


# `ijtiff` 2.2.4

## BUG FIXES
* Fix for `configure` error messages.


# `ijtiff` 2.2.3

## BUG FIXES
* Make `configure` more portable by using `sh` instead of `bash`.


# `ijtiff` 2.2.2

## BUG FIXES
* Insist on bug-fixed `strex` >= 1.4.


# `ijtiff` 2.2.1

## BUG FIXES
* Insist on `strex` >= 1.3.1 to avoid a garbage collection issue.


# `ijtiff` 2.2.0

## NEW FEATURES 
* The package now works on 32-bit Windows (thanks to PR #12 from Jeroen Ooms).

## BUG FIXES
* Fix tests by making use of `testthat::test_path()`.


# `ijtiff` 2.1.2

## BUG FIXES
* Fix some typos in the vignettes.


# `ijtiff` 2.1.1

## BUG FIXES
* Fix a PROTECTion error.


# `ijtiff` 2.1.0

## NEW FEATURES
* Add support for images with colormaps (also known as lookup tables (LUTs)).
* Add a print method for `ijtiff_img`s.


# `ijtiff` 2.0.5

## BUG FIXES
* Fix a test that failed due to breaking changes in `tibble`.


# `ijtiff` 2.0.4

## MINOR IMPROVEMENTS
* Include rOpenSci docs in `DESCRIPTION` as `URL`.

## BUG FIXES
* Sometimes `pkg-config` declares that `ijtiff` needs JBIG_KIT (compile flag `-ljbig`) at compile time. This is incorrect and it often causes users installation pain. This fix is a hack that removes this compile flag from the `pkg-config` output.


# `ijtiff` 2.0.3

## BUG FIXES
* `libjpeg` needs to be in `SystemRequirements`.


# `ijtiff` 2.0.2

## BUG FIXES
* For _ImageJ_-written images, if `n_slices` and `n_frames` are both specified, that should be OK if they're equal.


# `ijtiff` 2.0.1

## BUG FIXES
* Insist on latest, bug-fixed `filesstrings` 3.1.5.


# `ijtiff` 2.0.0

## BREAKING CHANGES
* `get_tiff_tags_reference()` is now `tif_tags_reference()`.
* `count_imgs()` is now `count_frames()`.

## NEW FEATURES
* It is now possible to read only certain frames of a TIFF image thanks to the `frames` argument of `read_tif()`.
* `read_tif()` and `read_tags()` now have the aliases `tif_read()` and `tags_read()` to comply with the rOpenSci `object_verb()` style.

## BUG FIXES
* Include `sys/types.h` for greater type compatibility.


# `ijtiff` 1.5.1

## BUG FIXES
* Require necessary version of `glue`.
* Fix dimension-related bug in `as_EBImage()`.
* Require latest (less-buggy) `filesstrings`.


# `ijtiff` 1.5.0

## NEW FEATURES
* Allow ZIP compression (which seems to be the best).

## BUG FIXES
* `write_txt_img()` was using decimal points for integers (e.g. 3.000 instead of just 3).


# `ijtiff` 1.4.2

## BUG FIXES
* Hacky fix for `configure` script to deal with lack of `-ljbig` on Solaris.
* Trim the package to below 5MB by compressing a few TIFF files.


# `ijtiff` 1.4.1

## NEW FEATURES
* The package is now lighter in appearance because it doesn't explicitly depend on `tibble`.

## BUG FIXES
* The configure script now allows for needing `--static` with `pkg-config`.


# `ijtiff` 1.4.0

## NEW FEATURES
* A `pkgdown` website.

## MINOR IMPROVEMENTS
* Better vignettes.
* Better error messages.


# `ijtiff` 1.3.0

## NEW FEATURES
* Conversion functions `linescan_to_stack()` and `stack_to_linescan()` useful for FCS data.


# `ijtiff` 1.2.0

## MINOR IMPROVEMENTS
* Improved the description of the package in DESCRIPTION, vignette and README.
* Added a hex sticker.
* Limited support for tiled images thanks to new author Kent Johnson.
* `write_tif()` is now slightly (<10%) faster.
* `write_tif()` messages are now more informative.


# `ijtiff` 1.1.0 

## NEW FEATURES
* `count_imgs()` counts the number of images in a TIFF file without reading the images themselves.
* `read_tags()` reads the tags from TIFF images without reading the images themselves.

## MINOR IMPROVEMENTS
* Now includes citation information.
* C code is more readable.
* `display()` is more flexible, accepting 3 and 4-dimensional arrays, just displaying the first frame from the first channel.


## `ijtiff` 1.0.0

#### PEER REVIEW
* The package is now peer reviewed by ROpenSci.


## `ijtiff` 0.3.0

#### MINOR IMPROVEMENTS
* Improve README and vignette with more tangible and fun example.

#### BUG FIXES
* Fix windows `libtiff` issues (thanks to Jeroen Ooms).
* Found some ImageJ-written TIFFs that weren't being read correctly and fixed that.
* Fix `protection stack overflow` error for TIFFs with many images.


## `ijtiff` 0.2.0
* First CRAN release.

#### MINOR IMPROVEMENTS
* Include handy shortcuts for 2- and 3-dimensional arrays.
* Messasges to inform the user about what kind of image is being read/written.


## `ijtiff` 0.1.0

* First github release.
# Test environments
* local OS X install, R 4.0.3
* ubuntu 16.04 (on travis-ci), R 4.0.3
* win-builder (devel and release)

### R CMD check results
0 ERRORs | 0 WARNINGs | 0 NOTEs

### Reverse dependencies
There are 3 reverse dependencies: `detrendr`, `nandb` and `autothresholdr`. This update does not break any of these. 
---
title: '`ijtiff`: An R package providing TIFF I/O for _ImageJ_ users'
tags:
  - R
  - file
  - string
authors:
  - name: Rory Nolan
    orcid: 0000-0002-5239-4043
    affiliation: 1
  - name: Sergi Padilla-Parra
    orcid: 0000-0002-8010-9481
    affiliation: 1, 2
affiliations:
  - name: Wellcome Centre Human Genetics, University of Oxford
    index: 1
  - name: Department of Structural Biology, University of Oxford
    index: 2
date: 26 February 2018
bibliography: paper.bib
nocite: | 
  @R, @RStudio, @checkmate, @magrittr, @filesstrings, @stringr, @readr, @purrr, @Rcpp, @fields, @grDevices, @knitr, @testthat, @rmarkdown, @covr, @devtools, @exampletestr, @BioFormats, @libtiff
---

# Summary
_ImageJ_ [@ImageJ] is the image viewing and processing GUI of choice for many in the fields of biology and microscopy. It is free and open-source. `ijtiff` is an R package which can correctly import TIFF files that were saved from _ImageJ_ and write TIFF files than can be correctly read by _ImageJ_. Due to the sometimes strange way that _ImageJ_ writes TIFF files, the original R `tiff` package [@tiff] may not correctly recognise their channel structure. 
`ijtiff` also goes beyond `tiff` in facilitating the writing of floating point (real-numbered) TIFF files from R. 

`ijtiff` reads TIFF pixel values in their native (usually integer) form, whereas `tiff` scales pixel values to the range [0, 1] by default. Hence and for other reasons, `ijtiff` should be viewed as a package with different capabilities and behaviours from the original `tiff` package, and not as an extension thereof. 

TIFF files are not always enough: they have maximum allowed values and their 32-bit floating point real-number representation can lack precision. For these extreme cases, `ijtiff` also supports text image I/O. Text images have no such limitations and are completely compatible with _ImageJ_. 

# References---
slug: "ijtiff"
title: Forcing Yourself to Make Your Life Easier
package_version: 1.0.0
authors:
  - name: Rory Nolan
    url: https://github.com/rorynolan
date: 2018-03-13
categories: blog
topicid: 
tags:
- r
- community
- software
- review
- onboarding
- package
- tiff
- data-import
- ImageJ
---

## The general struggle

Something that will make life easier in the long-run can be the most difficult thing to do _today_. For coders, prioritising the long term may involve an overhaul of current practice and the learning of a new skill. This can be painful for a number of reasons:

1.  We have to admit to ourselves that we’ve been doing things inefficiently (i.e. wasting time). This makes us feel stupid and fosters a sense of missed opportunity: we could’ve done something cool with the time we’d have saved (e.g. vacation).
2.  We’re fond of our existing methods, probably because we’re used to them and they’ve served us pretty well thus far.
3.  The learning we have to do may seem beneath us: “I’m an expert in R, I’ve even written my own package. Surely ‘Version Control for Beginners’ wasn’t intended for the likes of me.” This kind of thought permits us to dismiss good ideas.
4.  We’re tired. Overhauling workflows and learning new skills takes energy and concentration. Today, we’ll go through the motions just like we did yesterday. And hey, no one complained about the work we did yesterday.
5.  We feel like the task is actually beyond us. This is almost never true (so long as we’re good at asking for help).
6. In the context of a work environment, when we’re concentrated on learning, our apparent output drops (often to zero) until the learning period ends. A bad manager/supervisor who doesn’t appreciate the worth of taking the time to learn how to do things better may not forgive this short-term drop in apparent output. .

## One specific struggle: the _ImageJ_ TIFF problem

I’m doing a PhD in image analysis, working with a lot of microscopy images of cells, all of which are in the [TIFF](https://en.wikipedia.org/wiki/TIFF) format. [_ImageJ_](https://imagej.nih.gov/ij/) is a nice software for image viewing and processing. R is a nice software for image processing but is not as good as _ImageJ_ for viewing and playing around. For me, as an R enthusiast, it was ideal to do my image processing in R and my viewing in _ImageJ_.

On [CRAN](https://cran.r-project.org/), the [`tiff`](https://cran.r-project.org/package=tiff) and [`magick`](https://cran.r-project.org/package=magick) packages can both read TIFF files, and on [Bioconductor](https://www.bioconductor.org/), there’s [`EBImage`](https://bioconductor.org/packages/EBImage/). However, all of these packages sometimes struggle with TIFF files written from _ImageJ_ in that they wrongly perceive some images to have only 1 channel when in fact they have many (channels encode colour information: colour images have 3 channels (red, green and blue) whereas black and white images have one channel (grey)). Once the images are read into R, I am able to rejig them (with a combination of `aperm()` and `abind()`) into the format I want. However, the mistakes that the packages make vary, and thus the images require different rejigs. Every time I wanted to read an image into R, the process was _read, check, rejig, check_. This is a lot longer than just _read_. Nonetheless, this lengthy reading process still wasn’t that long in absolute terms: I could read in an image and have it in the format I wanted in about a minute. I processed thousands of images over the following two years (with various image analysis techniques). For each one, I went through _read, check, rejig, check_ . . .

In 2016, I attended Bioconductor’s [CSAMA](http://www-huber.embl.de/csama2018/) event in Italy. [Jenny Bryan](https://www.stat.ubc.ca/~jenny/) was there advocating the [tidyverse](https://www.tidyverse.org/), [happygitwithr](http://happygitwithr.com/) and various other good ideas. At that time, I wasn’t much into the R scene, so I’d never heard of Jenny nor of many of the things she mentioned, but I was very taken by her teaching (she grabbed my attention with the tip that `Alt`+`-` in RStudio gets you the assignment operator `<-`). The advice that resonated most was that _NOW is the time to do that thing that you know you ought to do eventually._ We all know that this is good advice, but that doesn’t mean we don’t need to be told it frequently. 

Jenny advised that every regular R user should write a _package_ containing the functions that they always use (rather than `source()`ing them from some random place atop every `.R` file). So I learned git, made a GitHub account and wrote my first R package. I decided to start easy, gathering together all of the functions that I use to rename and organise my files. This became (the terribly named) [`filesstrings`](https://cran.r-project.org/package=filesstrings) (now on CRAN). The most useful thing it can do is to extract numbers from strings with functions like `first_number()`. `filesstrings` used `Rcpp`, so I got comfortable writing R packages that use C++. Having completed `filesstrings`, I looked more seriously into writing a package that solved my _ImageJ_ TIFF issue. It was clear that I would require the `libtiff` C library, and thus would need to use R’s C interface, so I looked for documentation of it. The section in Hadley Wickham’s [_Advanced R_](http://adv-r.had.co.nz/C-interface.html) and his GitHub repo [r-internals](https://github.com/hadley/r-internals) are both good resources, but I still feel that this facet of R is under-documented. . As a result, I couldn’t really get to grips with R’s C interface (at least not within a week as I’d hoped), so I gave up and returned to _read, check, rejig, check_.

## Eventually

In the end, I did write `ijtiff`—the package I had been longing for—using `libtiff` and R’s C interface. What happened that pushed me over the edge to learn everything I needed to know about C and R’s interface to it and to slave away for six weeks with `clang`  errors and RStudio crashes and `valgrind`  in Linux virtual machines? Nothing in particular. I was doing yet another _read, check, rejig, check_  and I had had enough. I decided to take another crack at `ijtiff`. This time, with the help of what I’d learned during all of my previous failed attempts, I made it over the line.

So of course I was able to do it all along, it just took longer than I thought (it’s just reading a TIFF file into an array, right?). Most things take longer than we think, so that shouldn’t necessarily discourage us so immediately. Of course I should’ve done it right at the start of my PhD (rather than at the start of my 4th and final year), but for a combination of the silly reasons outlined above, I hesitated. _I hesitated for three years._ Anyway, it’s done now and I’m very happy with it. My feeling of delight that it’s now done easily trumps the feeling that I wasted time by waiting so long to do it. My R life is more efficient and less frustrating now.

## Immediate future

As mentioned previously, there is a lack of accessible information on R's C interface. The demand on a few experts like [Hadley Wickham](http://hadley.nz/), [Kevin Ushey](https://kevinushey.github.io/), [Jim Hester](https://github.com/jimhester) etc. to provide resources for the rest of us is huge. Although these people are very generous and well-appreciated, they only have a finite amount of time. As such, I want to contribute to the documentation on R's C interface. 

## Conclusions

1.  Find and heed all of [Jenny Bryan](https://www.stat.ubc.ca/~jenny/)’s advice.
2.  If there’s something that will make your life easier in the long-run but seems like a bit of a hurdle right now, do it. If necessary, take the time to explain to others around you why it’s a good idea and _do it_. If it takes longer than you thought to learn how to do it, don’t conclude from that that you can’t do it, just conclude that it’s not very easy and get all the help you can find.
3.  If it's a lack of skills that prevents you from doing what you wish you could do, take the time to learn those skills. They’ll make your life easier again and again in future.
4.  If you figure out how to do something that you wish was made easier in documentation, please—for the likes of me—contribute that new knowledge to existing documentation.

### Thanks

* Thanks to [ROpenSci](https://ropensci.org/) for their very helpful review of the `ijtiff` package. This type of review where the reviewers actively help you as well as objectively evaluating your work is a revelation. Getting the advice and help of people like [Jon Clayden](https://github.com/jonclayden) and [Jeroen Ooms](https://github.com/jeroen) was invaluable. Thanks also to editor [Scott Chamberlain](https://github.com/sckott) and to community manager [Stefanie Butland](https://github.com/stefaniebutland) who is helping me with this blog.
* Thanks to [Jenny Bryan](https://www.stat.ubc.ca/~jenny/) for great advice.
* Thanks to [Hadley Wickham](http://hadley.nz/), [Kevin Ushey](https://kevinushey.github.io/), [Jim Hester](https://github.com/jimhester) and others for making R usable for the rest of us.
* Thanks to my supervisor [Sergi Padilla-Parra](https://www.researchgate.net/profile/Sergi_Padilla-Parra2) and the other members of his group for having great ideas, believing in my ideas and making my PhD fun and stress free.
* Thanks to my lovely fiancée Naomi for her continued support throughout my PhD and for editing my lacklustre writing.# print method works

    Code
      read_tif(test_path("testthat-figs", "Rlogo-banana-red.tif"))
    Message <simpleMessage>
      Reading Rlogo-banana-red.tif: an 8-bit, 155x200 pixel image of
      unsigned integer type. Reading 1 channel and 2 frames . . .
       Done.
    Message <cliMessage>
      155x200 pixel ijtiff_img with 1 channel and 2 frames.
      Preview (top left of first channel of first frame):
    Output
           [,1] [,2] [,3] [,4] [,5] [,6]
      [1,]  255  255  255  255  255  255
      [2,]  255  255  255  255  255  255
      [3,]  255  255  255  255  255  255
      [4,]  255  255  255  255  255  255
      [5,]  255  255  255  255  255  255
      [6,]  255  255  255  255  255  255
      -- TIFF tags -------------------------------------------------------------------
    Message <cliMessage>
      * bits_per_sample: 8
      * samples_per_pixel: 1
      * sample_format: uint
      * planar_config: contiguous
      * rows_per_strip: 155
      * compression: none
      * description: ImageJ=1.51s images=2 slices=2 loop=false
      * color_space: palette
      * color map: matrix with 256 rows and 3 columns (red, green, blue)

<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/ijtiff/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/ropensci/ijtiff.git`
* Make sure to track progress upstream (i.e., on our version of `ijtiff` at `ropensci/ijtiff`) by doing `git remote add upstream https://github.com/ropensci/ijtiff.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/ijtiff`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
