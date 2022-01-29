---
title: 'TidyTensor: Utilities for multidimensional arrays as named hierarchical structures'
tags:
  - R
  - deep learning
  - array
  - tensor
  - tidy
authors:
  - name: Shawn T. O'Neil
    orcid: 0000-0001-6220-7080
    affiliation: 1
affiliations:
 - name: Center for Health Artificial Intelligence, University of Colorado Anschutz Medical Campus
   index: 1
date: 13 July 2021
bibliography: paper.bib
---

# Summary

Deep learning applications commonly employ the use of *tensors* (multidimensional arrays) for data storage and manipulation. In addition to being space efficient tensors provide a common framework for, the engine of deep neural network training [@baydin2018automatic]. TidyTensor is an R package for inspecting and manipulating tensors. It provides an improved `print()` function for summarizing structure, named tensors, conversion to data frames, and high-level manipulation functions. Designed to complement `keras` package (or other packages that utilize native R arrays), functionality is layered on top of base R types.


# Statement of need

While R supports multidimensional arrays natively, there are important differences between these types and tensors as typically used in deep learning applications. First, R assumes array data are organized in a column-major order as opposed to the row-major order typically used in Python, rendering the default `print()` and other inspection methods far less useful. TidyTensor additionally provides support for "named" tensors, allowing researchers to use semantically-relevant names for working with and manipulating tensors via a familiar tidy-friendly interface. While other packages exist for named arrays in R, TidyTensor leverages a hierarchical interpretation of tensor data that makes it easy to manipulate and investigate tensor data of various kinds. Early versions of TidyTensor supported a 3-week researcher-oriented workshop "Deep Learning for Life Scientists" at Oregon State University's Center for Genome Research and Biocomputing, and these features supported graduate students and faculty new to the field in applying deep learning techniques to datasets from their own work.


# Related Work

Other packages implement named tensors, though without the "tidy" features of TidyTensor or the simplification of asserting nesting/set-of-sets semantics. In Python, the [xarray](http://xarray.pydata.org/en/stable/) package provides named tensors and a wide variety of manipulation and related features [@hoyer2017xarray], while in R the [garray](https://cran.r-project.org/web/packages/garray/garray.pdf) package provides features similar to TidyTensor [@garray], including implementing names in the same way allowing for easy conversion to and from TidyTensors. The [stars](https://r-spatial.github.io/stars/) R package also provides named arrays with a focus on spatiotemporal data [@stars]. TidyTensor utilizes [abind](https://cran.r-project.org/web/packages/abind/index.html) internally for array manipulation [@abind].

TidyTensor was designed primarily as a companion to the [keras](https://keras.rstudio.com/) and [tensorflow](https://tensorflow.rstudio.com/) ports for R which utilize native R types [@rkeras; @rtensorflow]. By contrast, the recently released [torch](https://torch.mlverse.org/) port for R does not use native R types for tensors [@rtorch]; using TidyTensor with this package would require costly conversion to-and-from Torch-native types. 


# References
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4968727.svg)](https://doi.org/10.5281/zenodo.4968727)
[![status](https://joss.theoj.org/papers/07ef2e53d083c0eea30c0d08eef0f1cb/status.svg)](https://joss.theoj.org/papers/07ef2e53d083c0eea30c0d08eef0f1cb)
[![R-CMD-check](https://github.com/oneilsh/tidytensor/workflows/R-CMD-check/badge.svg)](https://github.com/oneilsh/tidytensor/actions)
[![codecov](https://codecov.io/gh/oneilsh/tidytensor/branch/master/graph/badge.svg?token=GWMT57CGDK)](https://codecov.io/gh/oneilsh/tidytensor) 


<br />
<img src="man/figures/tidytensor_transparent.png" height=200px/> 

*If you are reading this on [GitHub](https://github.com/oneilsh/tidytensor) be sure to check out the full [documentation](https://oneilsh.github.io/tidytensor/) on GitHub pages.*

TidyTensor is an R package for inspecting and manipulating tensors (multidimensional arrays). It provides an improved `print()` for summarizing structure, named tensors, conversion to data frames, and high-level manipulation functions. Designed to complement the  `keras` package for deep learning in R, functionality is layered on top of base R types.

TidyTensor was inspired by a workshop I taught in deep learning with R, and a desire to explain and explore tensors in a more intuitive way.  


<br />
<br />

## Installation

A simple `devtools::install_github("oneilsh/tidytensor")` will do it. If you don't have `devtools`, first grab it with `install.packages("devtools")`.

If you use TidyTensor, let us know in an [issue](https://github.com/oneilsh/tidytensor/issues/new)!

<br />
<br />

## Examples and Usage

Here we provide just two basic examples of how TidyTensor can help illuminate data used for deep learning. See the [Getting Started](https://oneilsh.github.io/tidytensor/articles/tidytensor.html) vignette for more examples of visualizing tensor structure, filtering and data augmentation, producing train/test splits, and other handy features.


Consider the `CIFAR10` dataset distributed with the `keras` library:

```
## library(keras)
cifar10_raw <- dataset_cifar10()$train$x 
```

TidyTensor can be used to plot the contained image data with the help of other `tidyverse` packages:

```
## library(tidytensor)
## library(ggplot2)
## library(tidyr)

cifar10_raw %>%
  as.tidytensor() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  subset(image = 1:4) %>%
  as.data.frame() %>%
  spread(channel, value) %>%
  ggplot() +
    geom_tile(aes(x = col, y = row, fill = rgb(R, G, B, maxColorValue = 255))) +
    facet_wrap( ~ image) +
    coord_equal() +
    scale_y_reverse() +
    scale_fill_identity()
```

<img src="man/figures/cifar_color.png" style="width: 60%; align: left" />

For a second example, we can start with a `keras` model and generate a function mapping inputs to internal feature map representations:

```
vgg_model <- application_vgg16(include_top = FALSE, input_shape = c(32, 32, 3))

input <- vgg_model$input
output <- get_layer(vgg_model, name = "block1_conv2")$output

# input shape (N, 32, 32, 3)
# output shape (N, 32, 32, 64) tensor, where last rank are feature maps
compute_featuremaps <- k_function(input, output)
```

TidyTensor can then be used to feed a given set of input tensors through the function and visualize the resulting convolutional feature maps:

```
(cifar10_raw / 255) %>%
  as.tidytensor() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  subset(image = 1:4) %>%
  compute_featuremaps() %>% 
  as.tidytensor() %>%
  set_ranknames(image, row, col, featuremap) %>%
  subset(featuremap = 1:6) %>%
  as.data.frame() %>%
  ggplot() +
    geom_tile(aes(x = col, y = row, fill = value)) +
    facet_grid(image ~ featuremap) +
    coord_equal() +
    scale_y_reverse() 
```

<img src="man/figures/feature_maps.png" style="width: 80%; align: left" />

<br />
<br />


## Contributing

Pull requests welcome! Please see the [`CONTRIBUTING.md`](https://github.com/oneilsh/tidytensor/blob/master/CONTRIBUTING.md) file for details. 
We are especially interested in additional methods for visualizing or summarizing the structure and content of TidyTensors. 

<br />
<br />

## Changelog

* v1.0.0: Minor documentation improvements, version 1.0 minted for JOSS publication
* v0.9.1: preparation for JOSS submission, many documentation improvements, removing `allow_huge` option from `as.data.frame.tidytensor()`
* v0.9.0: refactor `bottom` paramter of `print()` to `base_rank`
* v0.8.2: minor bugfixes, new combine_ranks() function
* v0.8.1: add [] and [] <- functionality
* v0.8: first version on github





# Contributing to TidyTensor

In general, the Tidyverse [contributing guide](https://www.tidyverse.org/contribute/) is a good place to start, and we welcome issues, documentation, and pull requests. 

When contributing to this repository, please first discuss the change you wish to make via issue. Please note we have a code of conduct (below), please follow it in all your interactions with the project.

<br />
<br />

## Specific Needs

Although TidyTensor's `print()` function is the original *raison d'être* for the package, there is still a lot of room for improvement (up to and including a full reimplementation). We also welcome additional suggestions for inspecting, browsing, or visualizing the structure and content of large tensors, for example via plots or interactive Shiny dashboards. 

Assigning values to TidyTensors is possible (with e.g. `t[, , , "R"] <- 0` to zero out the red channel in a channels-last image tensor), but the implementation is significantly slower than for base-R arrays where `[<-` is dispatched to a primitive function. Improvements or suggestions here are welcome.

Aside from these improvements to base functionality, additional functions for 'tensorifying' common data types would be beneficial. For example, early versions of TidyTensor included functionality to produce tensor-encoded versions of [fasta files](https://en.wikipedia.org/wiki/FASTA_format) commonly used in bioinformations applications--remnants of these efforts can be found in the source. Because these require dependencies (on e.g. `Rsamtools` and `GenomicRanges`) and are field-specific, these may be better placed in a separate package with a dependency on TidyTensor. 

<br />
<br />

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement via email to shawn@tislab.org.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
