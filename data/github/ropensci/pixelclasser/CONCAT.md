# pixelclasser

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci
peer-review](https://badges.ropensci.org/406_status.svg)](https://github.com/ropensci/software-review/issues/406)
<!-- badges: end -->

This package contains a set of tools to classify the pixels of digital
images into colour categories arbitrarily defined by the user. It is a
simple version of the multivariate technique known as Support Vector
Machine, adapted to this particular use.

The procedure is simple. A digital image in JPEG or TIFF format is
imported into R. The original image contains three colour variables (or
bands): \(R\), \(G\), and \(B\). The first step is to transform them
into proportions (\(r\), \(g\) and \(b\)), which simplifies the problem
into a bivariate one. The pixels of the test images can then be
represented in the plane defined by two of the variables (the user
judges which two are more convenient by trial and error) and, hopefully,
they would form separate clusters (pixel categories). The user then
traces straight lines (classification rules) that enclose the pixel
clusters. Using the mathematical expression for these rules and the
values of the transformed variables, each pixel can be classified in one
category. This produces a set of logical matrices (incidence matrices)
indicating which pixels belong to each category, stored in appropriate R
objects. These can be submitted to posterior analysis or used to create
a new version of the original image showing the category of each pixel.

`pixelclasser` contains functions to visualize the pixels of the images
and the rules created by the user, to create the rules and to store them
in objects that can be passed to function `classify_pixels()` for the
analysis of the image, and functions to import and export the original
and the classified images.

## Installation

You can install the last version from the rOpenSci repository in GitHub
using packages `remotes` or `devtools`, which install `remotes`

``` r
remotes::install_github("ropensci/pixelclasser", build_vignettes = TRUE)
devtools::install_github("ropensci/pixelclasser", build_vignettes = TRUE)
```

## Using pixelclasser

The manual with the description of each function and use examples is the
file `/doc/pixelclasser_1.0.0.pdf` (see the link to source code on the
right), but its contents can be found in the Reference section of this
website.

An example session is described in the vignette included in the package,
which can be accessed after installation in the usual way:

``` r
vignette("pixelclasser")
```

It also can be accessed in the section Get started in the top menu of
this page.

# Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# Contributing to pixelclasser

<!-- This CONTRIBUTING.md is adapted from https://gist.github.com/peterdesmet/e90a1b0dc17af6c12daf6e8b2f044e7c -->

If you are reading this document, you have probably been using `pixelclasser` in your work, so it is a pleasure to receive your comments, suggestions and bug reports. 

[repo]: https://github.com/ropensci/pixelclasser
[email]: mailto:carlos.real@usc.es

## How you can contribute

The purpose of `pixelclasser` is to classify RGB images using a simplified form of the multivariate technique known as Support Vector Machine. The functions that it contains are simple and few, and are designed to be integrated in your workflow as one of the initial steps in the analysis of your images.

Because it is a small piece of code, I hope that the number of bugs would be correspondingly small. Being simple, creating your own scripts using pixelclasser should be simple as well. Having said this, it is also true that bugs always creep in, and good ideas and improvements are always out there, as demonstrated the peer-review process of the package. So if you wish to report a bug or propose some improvement that you consider interesting or necessary, there are two ways:
* Open an issue in the code repository of pixelclaser: https://github.com/ropensci/pixelclasser
* Send me an e-mail to mailto:carlos.real@usc.es, if you are not used to gitHub or prefer a private conversation.

Making me know that you used `pixelclasser` in your research and how badly (or well) things gone, is another way to help to improve the package. Also, if you are creating your scripts (or package) that use `pixelclasser`, do not hesitate to ask me if you have some question or problem.

## Code of conduct

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.
