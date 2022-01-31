
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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)
original_files <- dir()
```

# ijtiff  <img src="man/figures/logo.png" height="140" align="right">

[![R-CMD-check](https://github.com/ropensci/ijtiff/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/ijtiff/actions)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ijtiff?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/ijtiff)
[![codecov](https://codecov.io/gh/ropensci/ijtiff/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/ijtiff)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ijtiff)](https://cran.r-project.org/package=ijtiff)
![RStudio CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/ijtiff)
![RStudio CRAN monthly downloads](http://cranlogs.r-pkg.org/badges/ijtiff)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)

[![DOI](http://joss.theoj.org/papers/10.21105/joss.00633/status.svg)](https://doi.org/10.21105/joss.00633)


## Introduction

This is a general purpose TIFF I/O utility for R. The [`tiff` package](https://cran.r-project.org/package=tiff) already exists for this purpose but `ijtiff` adds some functionality and overcomes some bugs therein. 

* `ijtiff` can write TIFF files whose pixel values are real (floating-point) numbers; `tiff` cannot. 
* `ijtiff` can read and write _text images_; `tiff` cannot.
* `tiff` struggles to interpret channel information and gives cryptic errors when reading TIFF files written by the _ImageJ_ software; `ijtiff` works smoothly with these images.

To learn about `ijtiff` and how to use it, visit the package website at https://docs.ropensci.org/ijtiff/.


## Installation

### `libtiff`

`ijtiff` requires you to have the `libtiff` C library installed. To install `libtiff`:

* On __Debian Linux__, try `sudo apt-get install libtiff5-dev`, or if that fails, try  
`sudo apt-get install libtiff4-dev`.
* On __Fedora Linux__, try `sudo yum install libtiff5-dev`, or if that doesn't work, try  
`sudo yum install libtiff4-dev`.
* On __Mac__, you need [Homebrew](https://brew.sh/). Then in the terminal, run `brew install libtiff`.
* On __Windows__, no setup is required `r emo::ji("smile")`.


### Installing the release version of the `ijtiff` R package

You can install `ijtiff` from CRAN (recommended) with:
```{r CRAN-installation, eval=FALSE}
install.packages("ijtiff")
```


### Installing the development version of the `ijtiff` R package

You can install the development version from GitHub with:
```{r GitHub-installation, eval=FALSE}
devtools::install_github("ropensci/ijtiff")
```


## Acknowledgement
This package uses a lot of code from the original `tiff` package by Simon Urbanek.


## Contribution
Contributions to this package are welcome. The preferred method of contribution is through a github pull request. Feel free to contact me by creating an issue. Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/ijtiff/blob/master/CONDUCT.md). By participating in this project you agree to abide by its terms.

```{r cleanup, include = FALSE}
new_files <- setdiff(dir(), original_files)
file.remove(new_files)
```

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "The _ImageJ_ Problem"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{The ImageJ Problem}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(magrittr)
```

## Introduction

The _ImageJ_ software (https://imagej.nih.gov/ij/) is a widely-used image viewing and processing software, particularly popular in microscopy and life sciences. It supports the TIFF image format (and many others). It reads TIFF files perfectly, however it can sometimes write them in a peculiar way, meaning that when other softwares try to read TIFF files written by _ImageJ_, mistakes can be made. 

One goal of the `ijtiff` R package is to correctly import TIFF files that were saved from _ImageJ_.
 
### Frames and Channels in TIFF files

* In a volumetric image, _frames_ are typically the different z-slices. In a time-stack of images (i.e. a video), each frame represents a time-point.
* There is one _channel_ per colour. A conventional colour image is made up of 3 colour channels: red, green and blue. A grayscale (black and white) image has just one channel. It's possible to acquire two channels (e.g. red an blue but not green), five channels (e.g. infrared, red, green, blue and ultraviolet), or any number at all, but these cases are seen mostly in specialist imaging fields like microscopy.

### The Peculiarity of _ImageJ_ TIFF files

It is common to use `TIFFTAG_SAMPLESPERPIXEL` to record the number of channels in a TIFF image, however _ImageJ_ sometimes leaves `TIFFTAG_SAMPLESPERPIXEL` with a value of 1 and instead encodes the number of channels in `TIFFTAG_IMAGEDESCRIPTION` which might look something like  
`"ImageJ=1.51 images=16 channels=2 slices=8"`.

A conventional TIFF reader would miss this channel information (because it is in an unusual place). `ijtiff` does not miss it. We'll see an example below. 

_Note_: These peculiar _ImageJ_-written TIFF files are still bona fide TIFF files according to the TIFF specification. They just break with common conventions of encoding channel information.

## Reading _ImageJ_ TIFF files

```{r 2 channel path}
path_2ch_ij <- system.file("img", "Rlogo-banana-red_green.tif",
  package = "ijtiff"
)
```
`path_2ch_ij` is the path to a TIFF file which was made in _ImageJ_ from the R logo dancing banana GIF used in the README of Jeroen Ooms' `magick` package. The TIFF is a time-stack containing only the red and green channels of the first and third frames of the original GIF. Here's the full gif:

![](`r system.file("img", "Rlogo-banana.gif", package = "ijtiff")`)

Here are the red and green channels of the first and third frames of the TIFF:

```{r red and green banana, echo=FALSE, message=FALSE, dpi=300, fig.height=1, warning=FALSE, fig.width=2}
rgbanana_tif <- system.file("img", "Rlogo-banana-red_green.tif",
  package = "ijtiff"
) %>%
  ijtiff::read_tif()
d <- dim(rgbanana_tif)
reds <- purrr::map(seq_len(d[4]), ~ rgbanana_tif[, , 1, .]) %>%
  purrr::reduce(cbind)
greens <- purrr::map(seq_len(d[4]), ~ rgbanana_tif[, , 2, .]) %>%
  purrr::reduce(cbind)
to_display <- array(0, dim = c(2 * nrow(reds), ncol(reds), 3, 1))
to_display[seq_len(nrow(reds)), , 1, ] <- reds
to_display[seq_len(nrow(reds)) + nrow(reds), , 2, ] <- greens
ijtiff::display(to_display)
```

### The original `tiff` package

When we import it with the original `tiff` package:
```{r original tiff import}
img <- tiff::readTIFF(path_2ch_ij, all = TRUE)
str(img) # 10 images
img[[1]][100:105, 50:55, 1] # print a section of the first image in the series
```

* We just get a list of `r length(img)` frames, with wrong information about the `r dim(img[[1]][3])` channels (it looks like there are 3 channels per frame).
* The numbers in the image array(s) are (by default) normalized to the range [0, 1].

### The `ijtiff` package

When we import the same image with the `ijtiff` package:
```{r ijtiff import}
img <- ijtiff::read_tif(path_2ch_ij)
dim(img) # 2 channels, 2 frames
img[100:105, 50:55, 1, 1] # print a section of the first channel, first frame
```

* We see the image nicely represented as an array of `r dim(img[[1]][3])` channels of `r dim(img[[1]][4])` frames.
* The numbers in the image are integers, the same as would be seen if one opened the image with ImageJ.


## Note
The original `tiff` package reads several types of TIFFs correctly, including many that are saved from _ImageJ_. This is just an example of a TIFF type that it doesn't perform so well with.


## Advice for all _ImageJ_ users
Base _ImageJ_ (similar to the `tiff` R package) does not properly open some perfectly good TIFF files^[I think native _ImageJ_ only likes 1, 3 and 4-channel images and complains about the rest, but I'm not sure about this.] (including some TIFF files written by the `tiff` and `ijtiff` R packages).  Instead it often gives you the error message: _imagej can only open 8 and 16 bit/channel images_. These images in fact can be opened in _ImageJ_ using the wonderful _Bio-Formats_ plugin. 
---
title: "Reading and Writing Images"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Reading and Writing Images}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Reading TIFF files

Check out the following video:

![](`r system.file("img", "Rlogo-banana.gif", package = "ijtiff")`)

As you can see, it's a colour video of a banana dancing in front of the R logo. Hence, it has colour channel (red, green and blue) and frame (a video is comprised of several _frames_) information inside. I have this video saved in a TIFF file.

```{r dancing-banana-path}
path_dancing_banana <- system.file("img", "Rlogo-banana.tif",
                                   package = "ijtiff")
print(path_dancing_banana)
```

To read it in, you just need `read_tif()` and the path to the image.

```{r read-dancing-banana}
pacman::p_load(ijtiff, magrittr)
img_dancing_banana <- read_tif(path_dancing_banana)
```

Let's take a peek inside of `img_dancing_banana`.

```{r peek}
print(img_dancing_banana)
```

You can see it's a `r length(dim(img_dancing_banana))`-dimensional array. The last two dimensions are `r dplyr::nth(dim(img_dancing_banana), -2)` and `r dplyr::nth(dim(img_dancing_banana), -1)`; this is because these are the channel and frame slots respectively: the image has `r dplyr::nth(dim(img_dancing_banana), -2)` channels (red, green and blue) and `r dplyr::nth(dim(img_dancing_banana), -1)` frames. The first two dimensions tell us that the images in the video are `r dim(img_dancing_banana)[1]` pixels tall and `r dim(img_dancing_banana)[2]` pixels wide. The image object is of class `ijtiff_img`. This guarantees that it is a 4-dimensional array with this structure. The attributes of the `ijtiff_img` give information on the various TIFF tags that were part of the TIFF image. You can read more about various TIFF tags at https://www.awaresystems.be/imaging/tiff/tifftags.html. To read just the tags and not the image, use the `read_tags()` function.

Let's visualize the constituent parts of that 8-frame, colour TIFF.

```{r red-blue-green-banana, echo=FALSE, message=FALSE, out.width='100%', dpi=300, fig.height=0.9}
d <- dim(img_dancing_banana)
reds <- purrr::map(seq_len(d[4]), ~ img_dancing_banana[, , 1, .]) %>% 
  purrr::reduce(cbind)
greens <- purrr::map(seq_len(d[4]), ~ img_dancing_banana[, , 2, .]) %>% 
  purrr::reduce(cbind)
blues <- purrr::map(seq_len(d[4]), ~ img_dancing_banana[, , 3, .]) %>% 
  purrr::reduce(cbind)
to_display <- array(0, dim = c(3 * nrow(reds), ncol(reds), 3, 1))
to_display[seq_len(nrow(reds)), , 1, ] <- reds
to_display[seq_len(nrow(reds)) + nrow(reds), , 2, ] <- greens
to_display[seq_len(nrow(reds)) + 2 * nrow(reds), , 3, ] <- blues
display(to_display)
```

There you go: 8 frames in 3 colours.


### Reading only certain frames

It's possible to read only certain frames. This can be a massive time and memory saver when working with large images.

Suppose we only want frames 3, 5 and 7 from the image above.

```{r threefiveseven}
img_dancing_banana357 <- read_tif(path_dancing_banana, frames = c(3, 5, 7))
```

Let's visualize again.

```{r red-bblue-green-banana357, echo=FALSE, message=FALSE, out.width='100%', dpi=300, fig.height=0.9}
d <- dim(img_dancing_banana357)
reds <- purrr::map(seq_len(d[4]), ~ img_dancing_banana357[, , 1, .]) %>% 
  purrr::reduce(cbind)
greens <- purrr::map(seq_len(d[4]), ~ img_dancing_banana357[, , 2, .]) %>% 
  purrr::reduce(cbind)
blues <- purrr::map(seq_len(d[4]), ~ img_dancing_banana357[, , 3, .]) %>% 
  purrr::reduce(cbind)
to_display <- array(0, dim = c(3 * nrow(reds), ncol(reds), 3, 1))
to_display[seq_len(nrow(reds)), , 1, ] <- reds
to_display[seq_len(nrow(reds)) + nrow(reds), , 2, ] <- greens
to_display[seq_len(nrow(reds)) + 2 * nrow(reds), , 3, ] <- blues
display(to_display)
```

Just in case you're wondering, it's not currently possible to read only certain channels. 

### More examples

If you read an image with only one frame, the frame slot (4) will still be there:

```{r one-frame, dpi=300, fig.height=0.5, fig.width=0.5}
path_rlogo <- system.file("img", "Rlogo.tif", package = "ijtiff")
img_rlogo <- read_tif(path_rlogo) 
dim(img_rlogo)  # 4 channels, 1 frame
class(img_rlogo)
display(img_rlogo)
```

You can also have an image with only 1 channel:

```{r one-channel, dpi=300, fig.height=0.5, fig.width=0.5}
path_rlogo_grey <- system.file("img", "Rlogo-grey.tif", package = "ijtiff")
img_rlogo_grey <- read_tif(path_rlogo_grey)
dim(img_rlogo_grey)  # 1 channel, 1 frame
display(img_rlogo_grey)
```



## Writing TIFF files

To write an image, you need an object in the style of an `ijtiff_img` object (see `help("ijtiff_img", package = "ijtiff")`). The basic idea is to have your image in a 4-dimensional array with the structure `img[y, x, channel, frame]`. Then, to write this image to the location `path`, you just type `write_tif(img, path)`. 

```{r write-tif}
path <- tempfile(pattern = "dancing-banana", fileext = ".tif")
print(path)
write_tif(img_dancing_banana, path)
```



## Reading text images

Note: if you don't know what text images are, see `vignette("text-images", package = "ijtiff")`.

You may have a text image that you want to read (but realistically, you might never).

```{r read-txt-img}
path_txt_img <- system.file("img", "Rlogo-grey.txt", package = "ijtiff")
txt_img <- read_txt_img(path_txt_img)
```


## Writing text images

Writing a text image works as you'd expect.

```{r, write-txt-img}
write_txt_img(txt_img, path = tempfile(pattern = "txtimg", fileext = ".txt"))
```---
title: "Text Images"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Text Images}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## What are text images

Text images are just arrays of numbers stored in a tab-delimited file (https://en.wikipedia.org/wiki/Tab-separated_values), where each location in the file is a pixel and the value stored there is the pixel intensity. 


## What are they good for?

### You might need to read them

Many softwares (e.g. Microsoft excel) don't support the saving of arrays as TIFF files. Many of these such softwares do support the saving of arrays as (tab-separated) text files. So, whether you like it or not, you might come across an image that was saved as a text file. You might not (if you're lucky), but I have, so being able to read them is handy. Beware that `ijtiff::read_txt_img()` assumes a tab-separated file (so something else like a CSV file won't work). This is the type of text image that you can save from _ImageJ_.

### You might (once in a million years) want to write them

A 32-bit TIFF file can only hold values up to $2^{32} - 1$, that's approximately $4 \times 10^9$. For whatever reason, this might not be enough for you, what if you want to write a value of $10^{10}$ to an image? Then, you're out of luck with TIFF files (and most if not all other image formats), but a text image is your friend. Text images place no restriction on the values therein. They're awkward and inefficient, but they can get you out of a hole sometimes.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ijtiff.R
\docType{package}
\name{ijtiff}
\alias{ijtiff}
\alias{ijtiff-package}
\title{\code{ijtiff}: TIFF I/O for \emph{ImageJ} users}
\description{
This is a general purpose TIFF I/O utility for R. The \href{https://cran.r-project.org/package=tiff}{\code{tiff} package} already exists for this
purpose but \code{ijtiff} adds some functionality and overcomes some bugs therein.
}
\details{
\itemize{
\item \code{ijtiff} can write TIFF files whose pixel values are real (floating-point)
numbers; \code{tiff} cannot.
\item \code{ijtiff} can read and write \emph{text images}; \code{tiff}
cannot.
\item \code{tiff} struggles to interpret channel information and gives cryptic
errors when reading TIFF files written by the \emph{ImageJ} software; \code{ijtiff}
works smoothly with these images.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructors.R
\name{ijtiff_img}
\alias{ijtiff_img}
\alias{as_ijtiff_img}
\title{\code{ijtiff_img} class.}
\usage{
ijtiff_img(img, ...)

as_ijtiff_img(img, ...)
}
\arguments{
\item{img}{An array representing the image. \itemize{\item For a
single-plane, grayscale image, use a matrix \code{img[y, x]}. \item For a
multi-plane, grayscale image, use a 3-dimensional array \code{img[y, x, plane]}.
\item For a multi-channel, single-plane image, use a 4-dimensional array
with a redundant 4th slot \code{img[y, x, channel, ]} (see \link{ijtiff_img}
'Examples' for an example). \item For a multi-channel, multi-plane image,
use a 4-dimensional array \code{img[y, x, channel, plane]}.}}

\item{...}{Named arguments which are set as attributes.}
}
\value{
A 4 dimensional array representing an image, indexed by \code{img[y, x, channel, frame]}, with selected attributes.
}
\description{
A class for images which are read or to be written by the \code{ijtiff} package.
}
\examples{
img <- matrix(1:4, nrow = 2) # to be a single-channel, grayscale image
ijtiff_img(img, description = "single-channel, grayscale")
img <- array(seq_len(2^3), dim = rep(2, 3)) # 1 channel, 2 frame
ijtiff_img(img, description = "blah blah blah")
img <- array(seq_len(2^3), dim = c(2, 2, 2, 1)) #  2 channel, 1 frame
ijtiff_img(img, description = "blah blah")
img <- array(seq_len(2^4), dim = rep(2, 4)) # 2 channel, 2 frame
ijtiff_img(img, software = "R")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_tif}
\alias{read_tif}
\alias{tif_read}
\title{Read an image stored in the TIFF format}
\usage{
read_tif(path, frames = "all", list_safety = "error", msg = TRUE)

tif_read(path, frames = "all", list_safety = "error", msg = TRUE)
}
\arguments{
\item{path}{A string. The path to the tiff file to read.}

\item{frames}{Which frames do you want to read. Default all. To read the 2nd
and 7th frames, use \code{frames = c(2, 7)}.}

\item{list_safety}{A string. This is for type safety of this function. Since
returning a list is unlikely and probably unexpected, the default is to
error. You can instead opt to throw a warning (\code{list_safety = "warning"})
or to just return the list quietly (\code{list_safety = "none"}).}

\item{msg}{Print an informative message about the image being read?}
}
\value{
An object of class \link{ijtiff_img} or a list of \link{ijtiff_img}s.
}
\description{
Reads an image from a TIFF file/content into a numeric array or list.
}
\details{
TIFF files have the capability to store multiple images, each having multiple
channels. Typically, these multiple images represent the sequential frames in
a time-stack or z-stack of images and hence each of these images has the same
dimension. If this is the case, they are all read into a single 4-dimensional
array \code{img} where \code{img} is indexed as \code{img[y, x, channel, frame]} (where we
have \verb{y, x} to comply with the conventional \verb{row, col} indexing of a matrix -
it means that images displayed as arrays of numbers in the R console will
have the correct orientation). However, it is possible that the images in the
TIFF file have varying dimensions (most people have never seen this), in
which case they are read in as a list of images, where again each element of
the list is a 4-dimensional array \code{img}, indexed as \code{img[y, x, channel, frame]}.

A (somewhat random) set of TIFF tags are attributed to the read image. These
are IMAGEDEPTH, BITSPERSAMPLE, SAMPLESPERPIXEL, SAMPLEFORMAT, PLANARCONFIG,
COMPRESSION, THRESHHOLDING, XRESOLUTION, YRESOLUTION, RESOLUTIONUNIT, INDEXED
and ORIENTATION. More tags should be added in a subsequent version of this
package. You can read about TIFF tags at
https://www.awaresystems.be/imaging/tiff/tifftags.html.

TIFF images can have a wide range of internal representations, but only the
most common in image processing are supported (8-bit, 16-bit and 32-bit
integer and 32-bit float samples).
}
\note{
\itemize{ \item 12-bit TIFFs are not supported. \item There is no
standard for packing order for TIFFs beyond 8-bit so we assume big-endian
packing}.
}
\examples{
img <- read_tif(system.file("img", "Rlogo.tif", package = "ijtiff"))
}
\seealso{
\code{\link[=write_tif]{write_tif()}}
}
\author{
Simon Urbanek wrote most of this code for the 'tiff' package. Rory
Nolan lifted it from there and changed it around a bit for this 'ijtiff'
package. Credit should be directed towards Lord Urbanek.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{linescan-conversion}
\alias{linescan-conversion}
\alias{linescan_to_stack}
\alias{stack_to_linescan}
\title{Rejig linescan images.}
\usage{
linescan_to_stack(linescan_img)

stack_to_linescan(img)
}
\arguments{
\item{linescan_img}{A 4-dimensional array in which the time axis is the first
axis. Dimension 4 must be 1 i.e. \code{dim(linescan_img)[4] == 1}.}

\item{img}{A conventional \link{ijtiff_img}, to be turned into a linescan image.
Dimension 1 must be 1 i.e. \code{dim(img)[1] == 1}.}
}
\value{
The converted image, an object of class \link{ijtiff_img}.
}
\description{
\code{ijtiff} has the fourth dimension of an \link{ijtiff_img} as its time dimension.
However, some linescan images (images where a single line of pixels is
acquired over and over) have the time dimension as the y dimension, (to avoid
the need for an image stack). These functions allow one to convert this type
of image into a conventional \link{ijtiff_img} (with time in the fourth dimension)
and to convert back.
}
\examples{
linescan <- ijtiff_img(array(rep(1:4, each = 4), dim = c(4, 4, 1, 1)))
print(linescan)
stack <- linescan_to_stack(linescan)
print(stack)
linescan <- stack_to_linescan(stack)
print(linescan)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.R
\name{write_tif}
\alias{write_tif}
\alias{tif_write}
\title{Write images in TIFF format}
\usage{
write_tif(
  img,
  path,
  bits_per_sample = "auto",
  compression = "none",
  overwrite = FALSE,
  msg = TRUE
)

tif_write(
  img,
  path,
  bits_per_sample = "auto",
  compression = "none",
  overwrite = FALSE,
  msg = TRUE
)
}
\arguments{
\item{img}{An array representing the image. \itemize{\item For a
single-plane, grayscale image, use a matrix \code{img[y, x]}. \item For a
multi-plane, grayscale image, use a 3-dimensional array \code{img[y, x, plane]}.
\item For a multi-channel, single-plane image, use a 4-dimensional array
with a redundant 4th slot \code{img[y, x, channel, ]} (see \link{ijtiff_img}
'Examples' for an example). \item For a multi-channel, multi-plane image,
use a 4-dimensional array \code{img[y, x, channel, plane]}.}}

\item{path}{file name or a raw vector}

\item{bits_per_sample}{number of bits per sample (numeric scalar). Supported
values are 8, 16, and 32. The default \code{"auto"} automatically picks the
smallest workable value based on the maximum element in \code{img}. For example,
if the maximum element in \code{img} is 789, then 16-bit will be chosen because
789 is greater than 2 ^ 8 - 1 but less than or equal to 2 ^ 16 - 1.}

\item{compression}{A string, the desired compression algorithm. Must be one
of \code{"none"}, \code{"LZW"}, \code{"PackBits"}, \code{"RLE"}, \code{"JPEG"}, \code{"deflate"} or
\code{"Zip"}. If you want compression but don't know which one to go for, I
recommend \code{"Zip"}, it gives a large file size reduction and it's lossless.
Note that \code{"deflate"} and \code{"Zip"} are the same thing. Avoid using \code{"JPEG"}
compression in a TIFF file if you can; I've noticed it can be buggy.}

\item{overwrite}{If writing the image would overwrite a file, do you want to
proceed?}

\item{msg}{Print an informative message about the image being written?}
}
\value{
The input \code{img} (invisibly).
}
\description{
Write images into a TIFF file.
}
\examples{
img <- read_tif(system.file("img", "Rlogo.tif", package = "ijtiff"))
temp_dir <- tempdir()
write_tif(img, paste0(temp_dir, "/", "Rlogo"))
img <- matrix(1:4, nrow = 2)
write_tif(img, paste0(temp_dir, "/", "tiny2x2"))
list.files(temp_dir, pattern = "tif$")
}
\seealso{
\code{\link[=read_tif]{read_tif()}}
}
\author{
Simon Urbanek wrote most of this code for the 'tiff' package. Rory
Nolan lifted it from there and changed it around a bit for this 'ijtiff'
package. Credit should be directed towards Lord Urbanek.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructors.R
\name{as_EBImage}
\alias{as_EBImage}
\title{Convert an \link{ijtiff_img} to an \link[EBImage:Image]{EBImage::Image}.}
\usage{
as_EBImage(img, colormode = NULL, scale = TRUE, force = TRUE)
}
\arguments{
\item{img}{An \link{ijtiff_img} object (or something coercible to one).}

\item{colormode}{A numeric or a character string containing the color mode
which can be either \code{"Grayscale"} or \code{"Color"}. If not specified, a guess
is made. See 'Details'.}

\item{scale}{Scale values in an integer image to the range \verb{[0, 1]}? Has no
effect on floating-point images.}

\item{force}{This function is designed to take \link{ijtiff_img}s as input. To
force any old array through this function, use \code{force = TRUE}, but take
care to check that the result is what you'd like it to be.}
}
\value{
An \link[EBImage:Image]{EBImage::Image}.
}
\description{
This is for interoperability with the the \code{EBImage} package.
}
\details{
The guess for the \code{colormode} is made as follows: * If \code{img} has an attribute
\code{color_space} with value \code{"RGB"}, then \code{colormode} is set to \code{"Color"}. *
Else if \code{img} has 3 or 4 channels, then \code{colormode} is set to \code{"Color"}. *
Else \code{colormode} is set to "Grayscale".
}
\examples{
if (rlang::is_installed("EBImage")) {
  img <- read_tif(system.file("img", "Rlogo.tif", package = "ijtiff"))
  str(img)
  str(as_EBImage(img))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{tif_tags_reference}
\alias{tif_tags_reference}
\title{TIFF tag reference.}
\source{
\url{https://www.awaresystems.be}
}
\usage{
tif_tags_reference()
}
\description{
A dataset containing the information on all known baseline and extended TIFF
tags.
}
\details{
A data frame with 96 rows and 10 variables: \describe{
\item{code_dec}{decimal numeric code of the TIFF tag}
\item{code_hex}{hexadecimal numeric code of the TIFF tag} \item{name}{the
name of the TIFF tag} \item{short_description}{a short description of the
TIFF tag} \item{tag_type}{the type of TIFF tag: either "baseline" or
"extended"} \item{url}{the URL of the TIFF tag at
\url{https://www.awaresystems.be}} \item{libtiff_name}{the TIFF tag name in
the libtiff C library} \item{c_type}{the C type of the TIFF tag data in
libtiff} \item{count}{the number of elements in the TIFF tag data}
\item{default}{the default value of the data held in the TIFF tag} }
}
\examples{
tif_tags_reference()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_tags}
\alias{read_tags}
\alias{tags_read}
\title{Read TIFF tag information without actually reading the image array.}
\usage{
read_tags(path, frames = 1)

tags_read(path, frames = 1)
}
\arguments{
\item{path}{A string. The path to the tiff file to read.}

\item{frames}{Which frames do you want to read tags from. Default first frame
only. To read from the 2nd and 7th frames, use \code{frames = c(2, 7)}, to read
from all frames, use \code{frames = "all"}.}
}
\value{
A list of lists.
}
\description{
TIFF files contain metadata about images in their \emph{TIFF tags}. This function
is for reading this information without reading the actual image.
}
\examples{
read_tags(system.file("img", "Rlogo.tif", package = "ijtiff"))
read_tags(system.file("img", "Rlogo-banana.tif", package = "ijtiff"),
  frames = c(2, 4)
)
}
\seealso{
\code{\link[=read_tif]{read_tif()}}
}
\author{
Simon Urbanek, Kent Johnson, Rory Nolan.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/text-io.R
\name{text-image-io}
\alias{text-image-io}
\alias{write_txt_img}
\alias{read_txt_img}
\alias{txt_img_write}
\alias{txt_img_read}
\title{Read/write an image array to/from disk as text file(s).}
\usage{
write_txt_img(img, path, rds = FALSE, msg = TRUE)

read_txt_img(path, msg = TRUE)

txt_img_write(img, path, rds = FALSE, msg = TRUE)

txt_img_read(path, msg = TRUE)
}
\arguments{
\item{img}{An image, represented by a 4-dimensional array, like an
\link{ijtiff_img}.}

\item{path}{The name of the input/output output file(s), \emph{without} a
file extension.}

\item{rds}{In addition to writing a text file, save the image as an RDS (a
single R object) file?}

\item{msg}{Print an informative message about the image being read?}
}
\description{
Write images (arrays) as tab-separated \code{.txt} files on disk. Each
channel-frame pair gets its own file.
}
\examples{
img <- read_tif(system.file("img", "Rlogo.tif", package = "ijtiff"))
tmptxt <- tempfile(pattern = "img", fileext = ".txt")
write_txt_img(img, tmptxt)
tmptxt_ch1_path <- paste0(strex::str_before_last_dot(tmptxt), "_ch1.txt")
print(tmptxt_ch1_path)
txt_img <- read_txt_img(tmptxt_ch1_path)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.R
\name{print.ijtiff_img}
\alias{print.ijtiff_img}
\title{Print method for an \code{ijtiff_img}.}
\usage{
\method{print}{ijtiff_img}(x, ...)
}
\arguments{
\item{x}{An object of class \link{ijtiff_img}.}

\item{...}{Not currently used.}
}
\value{
The input (invisibly).
}
\description{
Print method for an \code{ijtiff_img}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{display}
\alias{display}
\title{Basic image display.}
\usage{
display(img, method = NULL, basic = FALSE, normalize = TRUE)
}
\arguments{
\item{img}{An \link{ijtiff_img} object.}

\item{method}{The way of displaying images. Defaults to "browser" when R is
used interactively, and to "raster" otherwise. The default behavior can be
overridden by setting options("EBImage.display"). This has no effect when
\code{basic = TRUE}.}

\item{basic}{Force the basic (non-\code{EBImage}) display.}

\item{normalize}{Normalize the image before displaying (for better contrast)?
This only has an effect if the EBImage functionality is used. The basic
display always normalizes.}
}
\description{
Display an image that has been read in by \code{\link[=read_tif]{read_tif()}} as it would look in
'ImageJ'. This function is really just \code{\link[EBImage:display]{EBImage::display()}} on the inside. If
you do not have \code{EBImage} installed, a more basic display is offered.
}
\examples{
if (requireNamespace("EBImage")) {
  img <- read_tif(system.file("img", "Rlogo.tif", package = "ijtiff"))
  display(img)
  display(img[, , 1, 1]) # first (red) channel, first frame
  display(img[, , 2, ]) # second (green) channel, first frame
  display(img[, , 3, ]) # third (blue) channel, first frame
  display(img, basic = TRUE) # displays first (red) channel, first frame
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{count_frames}
\alias{count_frames}
\alias{frames_count}
\title{Count the number of frames in a TIFF file.}
\usage{
count_frames(path)

frames_count(path)
}
\arguments{
\item{path}{A string. The path to the tiff file to read.}
}
\value{
A number, the number of frames in the TIFF file. This has an
attribute \code{n_dirs} which holds the true number of directories in the TIFF
file, making no allowance for the way ImageJ may write TIFF files.
}
\description{
TIFF files can hold many frames. Often this is sensible, e.g. each frame
could be a time-point in a video or a slice of a z-stack.
}
\details{
For those familiar with TIFF files, this function counts the number of
directories in a TIFF file. There is an adjustment made for some
ImageJ-written TIFF files.
}
\examples{
count_frames(system.file("img", "Rlogo.tif", package = "ijtiff"))
}
