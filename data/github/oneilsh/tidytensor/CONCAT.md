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
---
title: "Getting Started"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteEngine{knitr::knitr}
  %\VignetteIndexEntry{Getting Started}
  %\usepackage[UTF-8]{inputenc}
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Terminology

*Tensors* as used in machine and deep learning are multi-dimensional array structures. In addition to being space efficient they provide a common framework for [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation), the engine of deep neural network training. 

Some packages, notably the R interfaces to [tensorflow](https://tensorflow.rstudio.com/) and [keras](https://keras.rstudio.com/), use native R structures for tensor storage, in the form of vectors (for 1-dimentional arrays), matrices (2 dimensions), and higher dimensional arrays. 

Tensors usually represent "set of" relationships. For example, a vector like `point <- c(20.2, 19.1, 24.8)` might represent a 1-dimensional tensor storing a point in space as X, Y, and Z coordinate values. A matrix of the form `points <- matrix(runif(10*3), ncol = 3, nrow = 10)` would represent 10 such points, with `points[6, ]` accessing the sixth. However, tensors are not often referred to by their "dimensionality" to avoid confusing array dimensions with data dimensions; rather, we'd call our point (a 1-d array) a "rank 1 tensor" and our matrix of points (a 2-d array) a "rank 2 tensor". 

On the other hand, it is common to refer to specific "dimensions" of an array as "axes" (Python's `numpy` refers to axes of multidimensional arrays for example). TidyTensor eschews this as well in favor of rank-based nomenclature; rather than referring to the first and second axes of `points` (slicing data by point and spatial coordinate respectively), we'll refer to the first and second ranks. Tensors also have a *shape*, which is conveniently returned by R's `dim()` function. The `point` tensor has shape `(3)` (parentheses are used to denote a shape vector) while `points` has shape `(10, 3)`. Lastly, "margins" is frequently used to refer to individual ranks/axes in R functions such as `apply()`. 

<br />
<br />

## Printing to check structure

While simple rank 1, 2, or 3 tensors are easy enough to understand and work with, things become more complex for higher-rank tensors. However, higher-rank tensors often encode "set of" relationships that we can exploit to quickly navigate data. 

Take for example a rank 2 tensor with shape `(28, 28)`, which might represent a small grayscale image. We can generate one with some fake data, and we run the result through `as.tidytensor()` which converts it to an array type (for consistency) and adds the `tidytensor` class to provide a custom printout:

```{r}
library(tidytensor)
library(magrittr)

image <- matrix(runif(28 * 28), nrow = 28, ncol = 28) %>% as.tidytensor()
image
```

Three of these stored in a rank 3 tensor of shape `(3, 28, 28)` (as an R `array()`) could represent an image with R, G, and B channel values. The `tt()` function is a shortcut for `as.tidytensor()`. 

```{r}
image <- array(runif(3*28*28), dim = c(3, 28, 28)) %>% tt()
image
```

The printout above is showing only the first of the 3 channel tensors, with `# ...` on the last line indicating more exist inside the first rank; we can see more by calling `print()` explicitly and adding a larger `max_per_level` argument:

```{r}
print(image, max_per_level = 2)
```

Similarly, a set of 24 of these might represent one second of a color video clip (24 RGB color images at a rate of 24 frames/second).

```{r}
video <- array(runif(24*3*28*28), dim = c(24, 3, 28, 28)) %>% tt() 
print(video, max_per_level = 2)
```

With `max_per_level = 2`, this printout is showing the first 2 channels of the first 2 frames of the video. 

We can continue in this way: a set of 256 of these might represent a training batch of 256 one-second color videos. 

```{r}
videos <- array(runif(256*24*3*28*28), dim = c(256, 24, 3, 28, 28)) %>% tt()
videos
```

This view provides a quick summary of the structure of a tensor, much better than R's default which assumes we're using [column-major ordering](https://www.r-bloggers.com/2013/07/column-major-confusion/) for high-dimensional arrays (output scrolls):

```{r, attr.output='style="max-height: 205px;"'}
unclass(videos) # remove tidytensor class to convert back to basic R array
```

The above tensor stores data in a "channels first" representation, where the channels rank comes before the image row and column ranks. It's not uncommon to use a "channels last" representation, in this case `print()` will notice the last rank has a size of 3 and assume we're printing an image in channels-last configuration. 

To help make the printout fit nicely we supply some extra parameters to `print()`. In the result, we can quickly see that each image is represented by a 28x28 matrix of length-3 R,G,B values (as opposed to three 28x28 matrices as above). 

```{r}
videos <- array(runif(256*24*28*28*3), dim = c(256, 24, 28, 28, 3)) %>% tt()
print(videos, max_rows = 4, max_cols = 3, signif_digits = 2, max_per_level = 2)
```

Not all tensors store image data, for example a 1 second audio clip sampled at 44.1 kHz might be stored as a rank 1 tensor of shape `(44100)`, ten of these would have shape `(10, 44100)`.

```{r}
samples <- array(runif(10*44100), dim = c(10, 44100)) %>% tt()
samples
```

Notice that the "base" of this tensor is displayed as a set of rank 1 tensors; in our channels-first representation we saw the base as a rank 2 tensor, and in the channels-last representation we saw rank 3 tensors (each displayed a grid of rank 1 tensors). These can be explicitly chosen with the `base_rank` option, which is possibly better illustrated with a small tensor of shape `(10, 3, 4, 2)`.

```{r}
data <- array(runif(10*3*4*2), dim = c(10, 3, 4, 2)) %>% tt()
print(data, base_rank = 1)
print(data, base_rank = 2)
print(data, base_rank = 3)
```

The default for `base_rank` is `NULL`, in which case the base-rank will be determined automatically by the shape of the tensor; if it looks like the tensor contains images in a channels-last configuration (last rank is of size 3 or 1) then `base_rank` will be set to `3`, if it looks like an image in channels-first configuration (3rd-last-rank is of size 3, or the last two ranks are of the same size) it will be set to `2`, otherwise it will be set to `1`.

<br />
<br />

## Named Ranks

TidyTensors support `ranknames()` in addition to traditional `dimnames()` (implemented as `names(dimnames())`) for meaningful rank annotation. Let's consider the CIFAR10 dataset distributed with `keras`. 

```{r}
library(keras)
cifar10_raw <- dataset_cifar10()$train$x 

cifar10_raw %>% tt()
```

This appears to be a set of 50000 32x32 RGB images in a channels-last organization. We can set ranknames with either `ranknames(t) <-` syntax or the tidy-friendly `set_ranknames()`, and the `.dots` parameter (accepted by most TidyTensor functions) can be used to provide ranknames as a character vector.

```{r}
# Base-R style
images <- tt(cifar10_raw)
ranknames(images) <- c("image", "row", "col", "channel")

# OR
images <- cifar10_raw %>% 
  tt() %>%
  set_ranknames(image, row, col, channel)

# OR
images <- cifar10_raw %>%
  tt() %>%
  set_ranknames(.dots = c("image", "row", "col", "channel"))
```

Named ranks make a variety of operations more explicit. For example to set dimension names for a single rank we can use `set_dimnames_for_rank()` rather than base-R `dimnames()`.

```{r}
cifar10_raw %>%
  tt() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  set_dimnames_for_rank(row, .dots = paste("row", 1:32)) %>%
  print(show_dimnames = T, max_cols = 5)
```
Dimension names are currently only shown for the last ranks displayed in the grid when using `show_dimnames`, and ranks without dimnames will show index numbers instead. In the above display the value `59` is accessible at index `[1, "row 1", 1, "R"]`.

<br />
<br />

## Subsetting, Data Frame Conversion

Named ranks are also useful when converting a tensor to a data frame. We’ll convert just the first four images, because tensors in data frame representation are significantly larger (approximately number-of-ranks times as large). Rather than selecting the first four images with `cifar10_raw[1:4, , , ]`, we'll demonstrate the `subset()` function for TidyTensors which supports a variety of syntax (see the `subset.tidytensor()` help doc examples):

```{r}
cifar10_raw %>%
  tt() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  set_dimnames_for_rank(row, .dots = paste("row", 1:32)) %>%
  subset(image = 1:4) %>%
  as.data.frame() %>%
  head()
```

Non-named tensors get generic column names. In the case of ranks with dimnames set as in `channel` and `row` above, factors are created with the level ordering determined by the ordering of the dimension names. This provides a nice interface for visualizing image tensor data with `ggplot2`.

```{r}
library(ggplot2)

cifar10_raw %>%
  tt() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  # set_dimnames_for_rank(row, .dots = paste("row", 1:32)) %>%    # no need for this for plotting
  subset(image = 1:4) %>%
  as.data.frame() %>%
  ggplot() +
    geom_tile(aes(x = col, y = row, fill = value)) +
    facet_grid(channel ~ image) +
    coord_equal()
```

It’s a little hard to make out, but these images are upside-down, because image data are typically encoded with an inverted y-axis, so next time we’ll add a `scale_y_reverse()` as well. To get fancy, we can use `tidyr::spread()` to create individual `R`, `G`, and `B` columns, combined with `rgb()` and `scale_fill_identity()` to merge the channels into color images.

```{r}
library(tidyr)

cifar10_raw %>%
  tt() %>%
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

These techniques work nicely for model investigation, for example in plotting internal feature maps produced by deep models. For an example, we'll start by importing a predefined model and creating a function that maps input tensors to feature maps using the `keras` API.

```{r}
vgg_model <- application_vgg16(include_top = FALSE, input_shape = c(32, 32, 3))

input <- vgg_model$input
output <- get_layer(vgg_model, name = "block1_conv2")$output

# input shape (N, 32, 32, 3)
# output shape (N, 32, 32, 64) tensor, where last rank are feature maps
compute_featuremaps <- k_function(input, output)
```

To visualize the feature maps we generate an output tensor, name it, convert it to data frame, select only the first six featuremaps with `subset()` to keep it reasonable, and then plot the values. Since this model assumes input data are in a 0-1 range but these image values are in a 0-255 range, we'll first divide the dataset by 255 (see the section on Applying below for examples of other normalization strategies).

```{r}
(cifar10_raw / 255) %>%
  tt() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B) %>%
  subset(image = 1:4) %>%
  compute_featuremaps() %>% 
  tt() %>%                           # compute_featuremaps() doesn't return a TidyTensor
  set_ranknames(image, row, col, featuremap) %>%
  subset(featuremap = 1:6) %>%
  as.data.frame() %>%
  ggplot() +
    geom_tile(aes(x = col, y = row, fill = value)) +
    facet_grid(image ~ featuremap) +
    coord_equal() +
    scale_y_reverse() 
```


<br />
<br />

## Permuting, Stitching, and Binding

A few functions are included for transforming or working with tensors. First up is `permute()` which allows us to re-order the ranks of a TidyTensor, for example to convert to a channels-first representation.

```{r}
cifar10 <- cifar10_raw %>%
  tt() %>%
  set_ranknames(image, row, col, channel) %>%
  set_dimnames_for_rank(channel, R, G, B)

cifar10

cifar10 %>% permute(image, channel, row, col)
```
The `stitch()` and `bind()` functions help us combine multiple tensors into one. To start with let's use `subset()` to get a few different subsets of our images.

```{r}
seta <- subset(cifar10, image = 1:4)      # shape (4, 32, 32, 3)
setb <- subset(cifar10, image = 5:8)      # shape (4, 32, 32, 3)
setc <- subset(cifar10, image = 9:12)     # shape (4, 32, 32, 3)
setd <- subset(cifar10, image = 101:200)  # shape (100, 32, 32, 3)
```

The `stitch()` function concatenates multiple tensors of the same shape *except* for the first rank, into a single tensor representing a concatenated set. (`stitch()` can also take a list of such tensors, which is why it's named `stitch()` and not `c()`.)

```{r}
merged <- stitch(seta, setb, setc, setd)

print(merged, max_cols = 5)
```
`bind()`, on the other hand, collects multiple tensors of the *exact* same shape into a new tensor with one higher rank. 

```{r}
bound <- bind(seta, setb, setc, new_rank_name = "set")
bound
```
<br />
<br />

## Partitioning, List Conversion, and Shuffling

The `partition()` function works as an inverse of `stitch()`, partitioning a tensor into a list of tensors of relative proportions along the first rank, which may be handy for generating train/validate/test splits (particularly when combined with `shuffle()`, below). The result is a list of TidyTensors.

```{r}
split_images <- partition(cifar10, c(0.1, 0.1, 0.8))   # given proportions will be normalized to sum to 1.0
split_images
```

Converting to a list with `as.list()` produces a list of *each* tensor from the first rank (with the first rank dropped). Note that the list names indicate the rank name, entry index, and shape of the contained tensor.

```{r}
cifar10 %>%
  subset(image = 1:3) %>%
  as.list() %>% 
  str()
```

By default `as.list()` splits along the first rank. It’s possible to split at ranks lower than the first, in which case all of the sub-tensors at that level become list elements. (This function can take while if the resulting list is very large, which is quite possible with large tensors.) In this example we first go to a channels-first representation, and then split out each channel of each list into a list element. 

```{r}
cifar10 %>%
  subset(image = 1:3) %>%
  permute(image, channel, row, col) %>%
  as.list(rank = channel)
```

If `flatten = FALSE`, a nested list is returned (compatible with for example the `data.tree` package, note that this output scrolls).

```{r, attr.output='style="max-height: 205px;"'}
cifar10 %>%
  subset(image = 1:3) %>%
  permute(image, channel, row, col) %>%
  as.list(rank = channel, flatten = FALSE) %>%
  str()
```

As a reminder with all this list-of-tensors making, `bind()` and `stitch()` can each take a list of tensors as opposed to one-per-parameter, so they can be used on results of `partition()` and `as.list()`. Here's a more involved example of grabbing every image where the green channel is brightest (from the first 100 to keep the build time for this vignette small): we first deconstruct the tensor into a list, then use `purrr::keep()` to filter them, and since each element is also a TidyTensor we can use `subset()` to select by channel dimension and run the resulting tensors though `mean()` before returning `TRUE` or `FALSE`. The resulting filtered list is then run through `bind()` re-specifying the new binding rankname as `"image"` (since `as.list()` strips off the first ranks for the contained tensors). 


```{r}
library(purrr)

cifar10 %>%
  subset(image = 1:100) %>%
  as.list() %>%
  keep(function(image_tensor) {
    mean_green <- image_tensor %>% subset(channel = "G") %>% mean()
    mean_other <- image_tensor %>% subset(channel = c("R", "B")) %>% mean()
    mean_green > mean_other
  }) %>% 
  bind(new_rank_name = "image")
```
Because the tensor is deconstructed into a list and then reconstructed, this isn't a very efficient way to solve this problem. We'll see a faster alternative below using `tt_apply()`. 

The `shuffle()` function simply shuffles a tensor along the first rank, for example `cifar10 %>% shuffle()` would return the same data, in the same shape, but entries in the first rank will be randomly permuted. A `seed` parameter can be used to set the shuffling random seed for repeatbility; in conjunction with `partition()` this can be used to create matched train/test splits. To illustrate this, we'll start by getting the corresponding target values for the `cifar10` images. First we check the structure:

```{r}
cifar10_targets_raw <- dataset_cifar10()$train$y
cifar10_targets_raw %>% tt()
```

Then we'll create a named version:

```{r}
cifar10_targets <- cifar10_targets_raw %>% tt() %>% set_ranknames(image, target)
cifar10_targets
```

Now we can similarly shuffle both datasets before running them through `partition()`. 

```{r}
cifar10_targets_split <- cifar10_targets %>%
  shuffle(seed = 42) %>%
  partition(c(0.2, 0.8))
  
print(cifar10_targets_split)

cifar10_split <- cifar10 %>%
  shuffle(seed = 42) %>%
  partition(c(0.2, 0.8))

print(cifar10_split, max_cols = 3)
```

<br />
<br />

## Applying, Combining Ranks


Finally we have `tt_apply()`, TidyTensor’s alternative to base-R `apply()` (we don't override `apply()` since TidyTensors are also native R types). While `apply()` allows applying a function over arbitrary ranks of an array with the `MARGIN` parameter, it not only removes `class` and other attributes from the returned output, it *deconstructs* each call’s return value into a vector, before prepending it to the remaining ranks. If `t` has shape `(500, 3, 26, 26)` and `normalize()` is a function that returns an array of the same shape as its input, `apply(t, MARGIN = c(1, 2), normalize)` doesn’t return a tensor of the same shape as `t`, but rather shape `c(26 * 26, 500, 3)`. The `tt_apply()` function uses `apply()` under the hood, but puts all the pieces back together nicely and is rankname-aware for both function inputs and outputs.

However, `tt_apply()` is limited compared to `apply()` in one respect, it only applies over the first N ranks, to stick with the sets-of-sets metaphor for tensors. If `t` has ranknames `c("video", "image", "channel", "row", "col")`, we can `tt_apply()` over every video, or every image (within each video), or every channel (within each within each image inside of each video), and so on. Should a different grouping be desired, `permute()` can help.

Consider the following `normalize()` function:

```{r}
# scales to [0, 1]
normalize <- function(t) {
  t <- t - min(t)
  t <- t / max(t)
  return(t)
}
```

For pre-processing, rather than scale the entire dataset to a [0,1] range, perhaps we’d like to do so for each channel of each image, maximizing the dynamic range utilized (that this might not actually be a useful strategy in a deep-learning context). If given a TidyTensor, this function normalizes values to a [0, 1] range and returns it with TidyTensor attributes intact. (Some R functions strip class and other attributes for their output, `apply()` and `scale()` are notable examples. Calling these on a TidyTensor will not result in a TidyTensor being returned.)

```{r}
images_channel_normalized <- cifar10 %>%
  subset(image = 1:100) %>%
  permute(image, channel, row, col) %>%       # permute to channels-first represtation
  tt_apply(channel, normalize) %>%            # apply normalize() to each channel
  permute(image, row, col, channel)           # permute back to channels-last

print(images_channel_normalized, max_cols = 3, signif_digits = 2)
```
Notice that here we’ve organized the data channel-first so we can do our per-image, per-channel normalization, then permuted back to the original. Even though `tt_apply()` is built on base-R `apply()`, applying can take some time as a function call is being executed for each entry asked for and there can be quite a few in a large tensor.

The function call results needn’t match their input shape: `tt_apply()` just needs returned result to be tensors (vectors, matrices, arrays, or TidyTensors) all of the same shape. TidyTensors are usually numeric, but they can also be logical or character; this allows us to use `tt_apply()` to filter to only images with a green channel brighter than user using an indexing strategy.



```{r}
selection <- cifar10 %>%
  subset(image = 1:100) %>%
  tt_apply(image, function(image_tensor) {
    mean_green <- image_tensor %>% subset(channel = "G") %>% mean()
    mean_other <- image_tensor %>% subset(channel = c("R", "B")) %>% mean()
    mean_green > mean_other
  })

selection

cifar10 %>% 
  subset(image = 1:100) %>%
  subset(image = selection)
```
Applying a function can also be useful for [data augmentation](https://en.wikipedia.org/wiki/Data_augmentation), where we want to generate
adjusted versions of input examples. Here's a function which duplicates a tensor 3 times adding Gaussian noise to each, returning the result as a new batch of tensors with a `replicate` rank to distinguish them. (This is not a very useful data augmentation method, but just meant to illustrate the technique.)

```{r}
augment <- function(t) {
  len <- length(t)             # length() returns number of elements regardless of shape
  t2 <- t + rnorm(n = len)     # default of 0 mean, 1 variance 
  t3 <- t + rnorm(n = len)
  t4 <- t + rnorm(n = len)
  return(bind(t, t2, t3, t4, new_rank_name = "replicate"))
}

```

Now we can apply this to each image of our normalized data:

```{r}
cifar10_augmented <- cifar10 %>% 
  subset(image = 1:100) %>%
  tt_apply(image, augment)

print(cifar10_augmented, max_cols = 4, max_per_level = 2)
```

However, we likely don't want to explictly model the `replicate` rank when training on this data; instead of a set of images each each containing a number of replicates, we'd prefer to just consider the replicas as individual images. The `combine_ranks()` function helps to accomplish exactly this, and can only be applied to consecutive ranks (again, `permute()` can be useful to get ranks in the right order).

```{r}
cifar10_augmented_combined <- cifar10_augmented %>% combine_ranks(image, replicate)

print(cifar10_augmented_combined, max_cols = 4, max_per_level = 2)
```
By default, the combined new rank name is built from the existing rank names being combined, we could alternatively specify the new name explicitly with e.g. `combine_ranks(image, replicate, new_rank_name = "image")`.



<!-- <div class="fold o"> -->
<!-- ```{r} -->
<!-- samples <- array(runif(10*44100), dim = c(10, 44100)) %>% tt() -->
<!-- samples -->
<!-- ``` -->
<!-- </div> -->
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranknames_various.R
\name{set_dimnames_for_rank}
\alias{set_dimnames_for_rank}
\title{Set dimnames() via a standard function call, for a particular rank.}
\usage{
set_dimnames_for_rank(x, rank, ..., .dots = NULL)
}
\arguments{
\item{x}{input tidytensor to set dimnames on.}

\item{rank}{rank to set the dimnames on.}

\item{...}{dimnames to assign (quoted or unquoted).}

\item{.dots}{character vector of dimnames to assign (quoted or unquoted).}
}
\value{
a tidytensor with dimnames set.
}
\description{
Sets the dimensions names for a particular rank, without requiring dimnames for the other ranks.
}
\details{
If all dimnames are unset, they will be set to NA for the other ranks, otherwise they will be left alone.
}
\examples{
t <- as.tidytensor(array(1:(3 * 2), dim = c(3, 2)))
t <- set_dimnames_for_rank(t, 2, valset1, valset2)
t <- set_dimnames_for_rank(t, 2, .dots = c("valset1", "valset2"))
print(t)
}
\seealso{
\code{\link{ranknames<-}}, \code{\link{dimnames}}, \code{\link{set_dimnames}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_ranks.R
\name{combine_ranks}
\alias{combine_ranks}
\title{Combine multiple ranks of a tensor into a single rank}
\usage{
combine_ranks(x, ..., new_rank_name = NULL, .dots = NULL)
}
\arguments{
\item{x}{the tidytensor to combine ranks for.}

\item{...}{ranknames or integers to combine (quoted or unquoted).}

\item{new_rank_name}{Name to give the newly combined rank; by default the new rank name is constructed from the names of the combined ranks.}

\item{.dots}{character or integer vector of ranknames.}
}
\value{
a new tidytensor.
}
\description{
Combine multiple ranks of a tensor into a single rank, for example for use in data augmentation.
}
\details{
If all ranks being combined have dimension names, the dimension names of the newly produced rank will be combinations of those specified.

It is only possible to combine consecutive ranks; use \code{permute()} to first organize ranks.
}
\examples{
# shape [5, 20, 26, 26] for 5 batches of 20 26x26 "images"
t <- as.tidytensor(array(rnorm(5 * 20 * 26 * 26), dim = c(5, 20, 26, 26)))
ranknames(t) <- c("batch", "image", "row", "col")

# given an image tidytensor (26x26), return a set of replicates with noise added
make_noisy_images <- function(t2) {
  res <- bind(t2,
              t2 + rnorm(length(t2)),
              t2 + rnorm(length(t2)),
              t2 + rnorm(length(t2)), new_rank_name = "replicate")
}

# augment the original data by replacing each image with a set of
# noisy replicates
t <- tt_apply(t, image, make_noisy_images)

# now t is shape (5, 20, 4, 26, 26)
# with ranknames (batch, image, replicate, row, col)
# let's set some dimension names

# setting to "1", "2", "3", ...
t <- set_dimnames_for_rank(t, image, .dots = 1:20)

# setting to "original", "rep1", "rep2", "rep3"
t <- set_dimnames_for_rank(t, replicate, original, rep1, rep2, rep3)

# to make it compatible with the original shape we
# combine images and replicates
t2 <- combine_ranks(t, image, replicate)

print(t2)

# since the combined ranks both have dimension names, the newly
# created rank does as well and we can verify contents
# here we see that the second batch, image 3, replicate 2 is indeed the same
print(t[2, "3", "rep2", , ])
print(t2[2, "3_rep2", , ])
}
\seealso{
\code{\link{permute}}, \code{\link{bind}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tt_apply.tidytensor.R
\name{tt_apply}
\alias{tt_apply}
\title{Apply a function over lower ranks of a tidytensor}
\usage{
tt_apply(x, rank = 1, FUN, flatten = FALSE, drop_final_1 = TRUE, ...)
}
\arguments{
\item{x}{the tidytensor to apply over.}

\item{rank}{an indicator of the rank to apply over (see details).}

\item{FUN}{the function to apply}

\item{flatten}{whether to preserve the higher-rank structure, or collapse into a single rank (see description).}

\item{drop_final_1}{If FUN returns a rank-0 tensor (length-1 vector), should it be collapsed? E.g. if final shape is (10, 10, 1), adjusts shape to (10, 10)}

\item{...}{additional arguments passed to FUN.}
}
\value{
a new tidytensor.
}
\description{
Applies a function over the lower ranks of a tidytensor, collecting
the results into a tidytensor. For example, if \code{FUN} is a function that takes a tidytensor
of shape [26, 26] and returns a tidytensor of shape [13, 13], then we could apply \code{FUN}
on a tidytensor of shape [3, 100, 26, 26] starting at rank 2 to get back one with shape [3, 100, 13, 13].
If \code{flatten = TRUE}, the higher ranks are collapsed to produce shape [300, 26, 26]

Ranknames are respected for both inputs and return values.
}
\details{
The \code{rank} argument should specify a single rank to apply over;
if \code{ranknames(t) <- c("sample", "rows", "cols", "channels")} then \code{rank = 2}, \code{rank = "rows"},
and \code{rank = c(FALSE, TRUE, FALSE, FALSE)} all indicate that \code{FUN} will be called on tidytensors
with ranknames \code{c("rows", "cols", "channels")}.
}
\examples{
# shape [20, 26, 26]
t <- as.tidytensor(array(rnorm(20 * 26 * 26), dim = c(20, 26, 26)))
ranknames(t) <- c("sample", "row", "col")
print(t)

# compute the deviation from median for each sample
dev_median <- function(t) {
  return(t - median(t))
}

median_deviations <- tt_apply(t, sample, dev_median)
print(median_deviations)
}
\seealso{
\code{\link{permute}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bind.tidytensor.R
\name{bind}
\alias{bind}
\title{Bind two or more tidytensors to create a new one with a new rank.}
\usage{
bind(..., new_rank_name = NULL)
}
\arguments{
\item{...}{one or more tidytensors, or a single list of them, to bind}

\item{new_rank_name}{a name (length-1 character vector) for the newly created rank.}
}
\value{
a new tidytensor.
}
\description{
Given multiple tidytensors, or a list of tidytensors, binds them together to create a tidytensor
of higher rank. For example, bind(x, y, z) where x, y, and z have shape [2, 3, 5] returns a new tidytensor
of shape [3, 2, 3, 5].
}
\details{
All input tidytensors must have the same shape. It's also possible to set a new rankname for the
newly created dimension; if ranknames were prevously unset lower ranknames are set to NA. If the input ranknames
conflict, only those of the first input tidytensor will be used, and a warning will be generated.
}
\examples{
# Three tidytensors of the same shape
t1 <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
t2 <- as.tidytensor(array(10 * 1:(3 * 4 * 5), dim = c(3, 4, 5)))
t3 <- as.tidytensor(array(100 * 1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t1) <- c("sample", "row", "col")
ranknames(t2) <- c("sample", "row", "col")
ranknames(t3) <- c("sample", "row", "col")
t4 <- bind(t1, t2, t3, new_rank_name = "batch")
print(t4)
}
\seealso{
\code{\link{ranknames}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranknames_various.R
\name{set_ranknames}
\alias{set_ranknames}
\title{Assign ranknames to a tidytensor via a standard function call.}
\usage{
set_ranknames(x, ..., .dots = NULL)
}
\arguments{
\item{x}{input tidytensor to set ranknames on.}

\item{...}{new ranknames to assign (quoted or unquoted).}

\item{.dots}{character vector of new ranknames to assign.}
}
\value{
a tidytensor with ranknames set.
}
\description{
A tidytensor t may have ranknames(t); this is a character vector of the same length as dim(t)
for future use. Note that ranknames(t) is independent of names(t) or dimnames(t); we are not naming elements,
or the dimension names for each rank, but rank names themselves.
Like names() and dimnames(), unset ranknames() are NULL.
}
\details{
Ranknames for a tidytensor t are stored as the names() attribute of dimnames(t). If dimnames(t) happens
to be null, before setting ranknames() we create valid dimnames() filled with NA values. The tidytensor
package also provides a specialized dimnames() which preserves ranknames when setting dimnames().
}
\examples{
t <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
t <- set_ranknames(t, sample, row, col)
t <- set_ranknames(t, .dots = c("sample", "row", "col"))
print(t)
}
\seealso{
\code{\link{ranknames<-}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subset.R
\name{subset.tidytensor}
\alias{subset.tidytensor}
\title{Subset dimensions of a tidytensor}
\usage{
\method{subset}{tidytensor}(x, ..., drop = TRUE)
}
\arguments{
\item{x}{the tidytensor to apply over.}

\item{...}{named or unnamed parameters specifying subsetting criteria (see examples)}

\item{drop}{whether to drop ranks with size 1 (similar to \code{x[..., drop = TRUE]})}
}
\value{
a tidytensor
}
\description{
A functional form of e.g. \code{tensor[1:10, 3, ]}, supporting selecting by ranknames, usage with %>%, and
indexing when the rank is unknown.
}
\details{
Subsetting a tidytensor with \code{subset()} as opposed to \code{[]} allows for subsetting even when the number of ranks of the input is unknown; see examples.
}
\examples{
# shape [100, 26, 26, 3]
t <- as.tidytensor(array(rnorm(100 * 26 * 26 * 3), dim = c(100, 26, 26, 3)))
ranknames(t) <- c("sample", "row", "col", "channel")
t <- set_dimnames_for_rank(t, channel, R, G, B)
print(t)

t2 <- subset(t, row = 1:10, sample = 27, drop = FALSE)
print(t2)

# same thing, but without named ranks (not a good idea to mixes named
# subsetting and non-named subsetting)
t2 <- subset(t, 27, 1:10, drop = FALSE)
print(t2)

# equiv of t3[1:20, , , c("G", "R", "B")] # since the last rank has dimnames()
# note the re-ordering of channel levels
t3 <- subset(t, sample = 1:20, channel = c("G", "R", "B"), drop = FALSE)
print(t3)

}
\seealso{
\code{\link{shuffle}}, \code{\link{permute}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.data.frame.tidytensor.R
\name{as.data.frame.tidytensor}
\alias{as.data.frame.tidytensor}
\title{Convert a tidytensor to a data.frame representation.}
\usage{
\method{as.data.frame}{tidytensor}(x, row.names = NULL, optional = FALSE, ...)
}
\arguments{
\item{x}{input to convert to a data.frame}

\item{row.names}{NULL (default) or character vector giving the new row names for the data frame (included for method compatibility with base \code{as.data.frame}).}

\item{optional}{Ignored (included for method compatibility with base \code{as.data.frame})}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\value{
a data.frame
}
\description{
Given a tidytensor, returns a data.frame, with each rank of the tensor being represented by a column.
Produces an error if the resulting data.frame would have more than 10 million entries and \code{allow_huge = FALSE}.
}
\details{
Note that this produces a row for each value in the tensor, and a column for each rank; data.frames are a much less
efficient representation, but can be useful for e.g. visualization purposes. This method thus produces an error if
the resulting data.frame would have more than 10 million entries and \code{allow_huge = FALSE} is set (default is \code{TRUE}).
If dimnames() are set (naming each dimension withina rank), then the columns will be factors, rather than integer indices.

If the tidytensor ranks are not named, columns will be named \code{index_1}, \code{index_2}, etc., otherwise they will be
set to ranknames.
Tensor values will be in a column named \code{value}.
}
\examples{
# From an array (representing e.g. 30 26x26 images (30 sets of 26 rows of 26 pixels))
a <- array(rnorm(30 * 26 * 26), dim = c(30, 26, 26))
t <- as.tidytensor(a)
ranknames(t) <- c("sample", "row", "pixel")
df <- as.data.frame(t)
print(head(df))

# Example with named dimensions:
dimnames(t)[[1]] <- paste("sample", 1:30, sep = "_")
dimnames(t)[[2]] <- paste("row", 1:26, sep = "_")
dimnames(t)[[3]] <- paste("pixel", 1:26, sep = "_")
# or with a list:
dimnames(t) <- list(paste("sample", 1:30, sep = "_"),
                    paste("row", 1:26, sep = "_"),
                    paste("pixel", 1:26, sep = "_"))

print(head(as.data.frame(t)))
}
\seealso{
\code{\link{ranknames}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.tidytensor.R
\name{as.tidytensor}
\alias{as.tidytensor}
\title{Convert a vector, matrix, or array to a tidytensor type.}
\usage{
as.tidytensor(x, ...)
}
\arguments{
\item{x}{input to convert to a tidytensor.}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\value{
a new tidytensor.
}
\description{
Given a vector, matrix, or array, returns a tidytensor.
If given a vector, converts to a 1-d array supporting dim(), matrices are left as matrices,
and in all cases the class 'tidytensor' is added.
}
\details{
Matrices are synonymous with 2-d arrays, so these are left as is. Vectors are converted
to 1-d arrays so that they can support dim().
}
\examples{
# From an array (representing e.g. 30 26x26 images (30 sets of 26 rows of 26 pixels))
a <- array(rnorm(30 * 26 * 26), dim = c(30, 26, 26))
t <- as.tidytensor(a)
ranknames(t) <- c("sample", "row", "pixel")
print(t)

# From a matrix (representing e.g. a 26x26 image (26 rows of 26 pixels))
m <- matrix(rnorm(26 * 26), nrow = 26, ncol = 26)
t <- as.tidytensor(m)
ranknames(t) <- c("row", "pixel")
print(t)

# From a vector (representing e.g. 26 pixel values)
v <- rnorm(26)
t <- as.tidytensor(v)
ranknames(t) <- c("pixel")
print(t)
}
\seealso{
\code{\link{tt}}, \code{\link{ranknames}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.tidytensor.R
\name{print.tidytensor}
\alias{print.tidytensor}
\title{Print a tidytensor.}
\usage{
\method{print}{tidytensor}(
  x,
  show_dimnames = FALSE,
  max_per_level = 1,
  base_rank = NULL,
  max_rows = 6,
  max_cols = 6,
  max_depth = 3,
  signif_digits = 3,
  indent = 0,
  ...
)
}
\arguments{
\item{x}{a tidytensor to summarize.}

\item{show_dimnames}{show the dimension names, if present, or dimension indices if not in base-level prints.}

\item{max_per_level}{only show this many sub-tensors per level.}

\item{base_rank}{either NULL, 1, 2, or 3 - specifies whether the inner/bottom-most tensors should be represented as rank 1, 2, or 3 in a grid (NULL for autodetect based on tensor shape, see details).}

\item{max_rows}{limit the base-level prints to include this many rows (also applies to 1d prints).}

\item{max_cols}{limit the base-level prints to include this many columns.}

\item{max_depth}{in 3d representation, limit the base-level prints to include this many entries of the last rank.}

\item{signif_digits}{number of significant digits to print for numeric tensors.}

\item{indent}{indent the printout by this much (used internally).}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\description{
Prints a summary of a tidytensor as a nested hierarchy of tensors of lower rank.
}
\details{
The \code{base_rank} argument specifies whether the lowest ranks of the tensor (displayed as a grid) should be shown as rank 2 tensors, rank 3 tensors, or rank 1 tensors; the default of \code{NULL} will
select 3 if the last rank is of size 3 or 1 (assuming an image and a "channels-last" convention), 2 if the 3rd-to-last rank is length 3 or 1 (assuming an image
and a "channels-first" convention) or if there are only two ranks or if the last two ranks are equal (assuming an image channel of some kind), and otherwise will default to 1.

\code{max_per_level} indicates how many replicates
}
\examples{
t <- as.tidytensor(array(1:(2 * 3 * 4 * 5), dim = c(2, 3, 4, 5)))
ranknames(t) <- c("samples", "batches", "rows", "cols")
print(t, base_rank = 2)

t <- as.tidytensor(array(1:(2 * 3 * 40 * 50 * 3), dim = c(2, 3, 40, 50, 3)))
ranknames(t) <- c("sample", "batch", "row", "pixel", "channel")
print(t, max_rows = 6, max_cols = 6, max_depth = 3, show_dimnames = TRUE, base_rank = 3)

}
\seealso{
\code{\link{print.tidytensor}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranknames_various.R
\name{set_dimnames}
\alias{set_dimnames}
\title{Set dimnames() via a standard function call.}
\usage{
set_dimnames(x, newnames, ...)
}
\arguments{
\item{x}{input tidytensor to set dimnames on.}

\item{newnames}{list of dimnames to assign.}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\value{
a tidytensor with dimnames set.
}
\description{
Since tidytensors are arrays, they support dimnames(). The usuall syntax dimnames(x) <- works;
this function provides a Magritte-compatible regular function, set_dimnames(x, newnames) which returns a new tidytensor.
}
\details{
Setting dimnames with set_dimnames() preserves any ranknames present.
}
\examples{
t <- as.tidytensor(array(1:(3 * 2), dim = c(3, 2)))
t <- set_dimnames(t, list(c("sample1", "sample2", "sample3"), c("valset1", "valset2")))
print(t)

# We can also assign ranknames:
ranknames(t) <- c("sample", "valset")
print(t)

}
\seealso{
\code{\link{ranknames<-}}, \code{\link{dimnames}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.list.tidytensor.R
\name{as.list.tidytensor}
\alias{as.list.tidytensor}
\title{Convert a tidytensor into a nested list of tensors.}
\usage{
\method{as.list}{tidytensor}(x, rank = 1, flatten = TRUE, state = NULL, ...)
}
\arguments{
\item{x}{the tidytensor to convert.}

\item{rank}{an indicator of the rank defining the contained tensors.}

\item{flatten}{whether to return a nested list (\code{FALSE}) or a flattened list (\code{TRUE}).}

\item{state}{an internally used parameter for tracking build state-do not set manually.}

\item{...}{additional arguments passed to methods (unusued).}
}
\value{
a list.
}
\description{
Convert a tidytensor into a nested list of tensors, nested down to level specified in \code{rank}.
If \code{flatten = TRUE}, returns a flattens the structure to a list of tensors (not nested).
}
\details{
The \code{state} parameter is for internal use, and needn't be set during normal usage.
}
\examples{
# Three tidytensors of the same shape
t1 <- as.tidytensor(array(100 * 1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t1) <- c("sample", "row", "col")
l1 <- as.list(t1)
str(l1)
}
\seealso{
\code{\link{as.data.frame.tidytensor}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranknames_various.R
\name{ranknames}
\alias{ranknames}
\title{Get ranknames of a tidytensor.}
\usage{
ranknames(x, ...)
}
\arguments{
\item{x}{input tidytensor to get ranknames for.}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\value{
character vector of the same length as dim(x), or NULL if unset.
}
\description{
A tidytensor t may have ranknames(t); this is a character vector of the same length as dim(t)
for future use. Note that ranknames(t) is independent of names(t) or dimnames(t); we are not naming elements,
or the dimension names for each rank, but rank names themselves.
Like names() and dimnames(), unset ranknames() are NULL.
}
\details{
Ranknames for a tidytensor t are stored as the names() attribute of dimnames(t). If dimnames(t) happens
to be null, before setting ranknames() we create valid dimnames() filled with NA values. The tidytensor
package also provides a specialized dimnames() which preserves ranknames when setting dimnames().
}
\examples{
t <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t) <- c("sample", "row", "col")
print(ranknames(t))
}
\seealso{
\code{\link{set_ranknames}}, \code{\link{ranknames<-}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permute.tidytensor.R
\name{permute}
\alias{permute}
\title{Permute the ranks of a tensor}
\usage{
permute(tensor, ..., .dots = NULL)
}
\arguments{
\item{tensor}{the tidytensor permute.}

\item{...}{ranknames or integers to permute by (quoted or unquoted).}

\item{.dots}{character or integer vector to permute by.}
}
\value{
a new tidytensor.
}
\description{
Permute the ranks of a tensor, for example to convert between "channels first" and "channels last" representations.

Ranknames are respected for both inputs and return values.
}
\details{
The \code{rank} parameter may be an integer numeric vector (for permuting by index), or character vector (for permuting by rankname).
}
\examples{
# shape [20, 26, 26]
t <- as.tidytensor(array(rnorm(20 * 26 * 26), dim = c(20, 26, 26)))
ranknames(t) <- c("sample", "row", "col")
print(t)

t2 <- permute(t, col, sample, row)
t2 <- permute(t, 3, 1, 2)
t2 <- permute(t, .dots = c(3, 1, 2))
t2 <- permute(t, .dots = c("col", "sample", "row"))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shuffle.R
\name{shuffle}
\alias{shuffle}
\title{Shuffle a tidytensor in the first rank.}
\usage{
shuffle(t, seed = NULL)
}
\arguments{
\item{t}{the tidytensor to apply over.}

\item{seed}{random seed to be used for shuffling.}
}
\value{
a tidytensor of the same shape.
}
\description{
Shuffle's the entries in the first rank of a tensor. For example, if
\code{x} has shape (3, 5, 5), it may be indexed as \code{x[c(2, 3, 1), , ]}.
It's possible to set a custom seed for repeatable shuffling (amongst tensors with
the same size in the first rank).
}
\details{
Since tidytensor consider tensors as representing hierarchical "set of" relationships,
shuffling in any rank other than the first would permute lower entities across set boundaries
in higher ranks. For example, in a set of color images of shape (500, 28, 28, 3), shuffling the last rank
would re-order the channels, but identically for all the images. See \code{\link{tt_apply}} for applying functions
(such as shuffle) over lower ranks of a tensor.
}
\examples{
# shape [100, 26, 26]
t <- as.tidytensor(array(rnorm(100 * 26 * 26), dim = c(100, 26, 26)))
ranknames(t) <- c("sample", "row", "col")
print(t)

t <- shuffle(t, seed = 42)

}
\seealso{
\code{\link{tt_apply}}, \code{\link{permute}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/partition.R
\name{partition}
\alias{partition}
\title{Partition a tidytensor into a list of smaller tidytensors of the same rank}
\usage{
partition(x, sizes = c(0.5, 0.5))
}
\arguments{
\item{x}{the tidytensor to apply over.}

\item{sizes}{relative sizes of partitions}
}
\value{
a list of tidytensors.
}
\description{
Partitions a tensor into pieces of sizes relative to \code{sizes}; e.g. a
tensor with shape (24, 50, 50, 3) partitioned with \code{partition(sizes = c(0.5, 0.5))}
results in a list of two tensors of shape (12, 50, 50, 3).

Ranknames are respected for both inputs and return values.
}
\details{
Entries in \code{sizes} are treated as relative, so \code{sizes = c(2, 1, 1)}
is equivalent to \code{sizes = c(0.5, 0.25, 0.25)}. Non-integer parition boundaries are
rounded down, and this may result in entries with shape (0, ...), but only when
the size of the first rank is smaller than the number of partitions requested.
}
\examples{
# shape [100, 26, 26]
t <- as.tidytensor(array(rnorm(100 * 26 * 26), dim = c(100, 26, 26)))
ranknames(t) <- c("sample", "row", "col")
print(t)

partitions <- partition(t, c(0.2, 0.8))
print(partitions)
}
\seealso{
\code{\link{c}}, \code{\link{permute}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stitch.tidytensor.R
\name{stitch}
\alias{stitch}
\title{Concatenate two or more tidytensors to create a new one with the same rank}
\usage{
stitch(...)
}
\arguments{
\item{...}{a number of tidytensors of the same shape, or a single list of them.}
}
\value{
a new tidytensor.
}
\description{
Given multiple tidytensors of the same shape except the first rank, concatenates them together to create a tidytensor
of the same shape, but larger in the first. For example, c(x, y, z) where x and have shape [2, 3, 5] and z has shape [10, 3, 5] returns a new tidytensor
of shape [14, 3, 5].
}
\details{
All input tidytensors must have the same shape except for the first rank. If the input ranknames
conflict, only those of the first input tidytensor will be used, and a warning will be generated.
}
\examples{
# Three tidytensors of the same shape
t1 <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
t2 <- as.tidytensor(array(10 * 1:(3 * 4 * 5), dim = c(3, 4, 5)))
t3 <- as.tidytensor(array(100 * 1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t1) <- c("sample", "row", "col")
ranknames(t2) <- c("sample", "row", "col")
ranknames(t3) <- c("sample", "row", "col")
t4 <- stitch(t1, t2, t3)
print(t4)

list_example <- list(t1, t2, t3)
t5 <- stitch(list_example)
print(t5)
}
\seealso{
\code{\link{bind}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.tidytensor.R
\name{tt}
\alias{tt}
\title{Convert a vector, matrix, or array to a tidytensor type.}
\usage{
tt(x, ...)
}
\arguments{
\item{x}{input to convert to a tidytensor.}

\item{...}{additional arguments to be passed to or from methods (ignored).}
}
\value{
a new tidytensor.
}
\description{
\code{tt()} is a convenience shorthand for \code{as.tidytensor()}. Given a vector, matrix, or array, returns a tidytensor.
If given a vector, converts to a 1-d array supporting \code{dim()}, matrices are left as matrices,
and in all cases the class 'tidytensor' is added.
}
\details{
Matrices are synonymous with 2-d arrays, so these are left as is. Vectors are converted
to 1-d arrays so that they can support \code{dim()}.
}
\examples{
# From an array (representing e.g. 30 26x26 images (30 sets of 26 rows of 26 pixels))
a <- array(rnorm(30 * 26 * 26), dim = c(30, 26, 26))
t <- tt(a)
ranknames(t) <- c("sample", "row", "pixel")
print(t)

# From a matrix (representing e.g. a 26x26 image (26 rows of 26 pixels)) using \%>\%
library(magrittr)
t <- matrix(rnorm(26 * 26), nrow = 26, ncol = 26) \%>\% tt()
ranknames(t) <- c("row", "pixel")
print(t)

# From a vector (representing e.g. 26 pixel values)
v <- rnorm(26)
t <- tt(rnorm(26))
ranknames(t) <- c("pixel")
print(t)
}
\seealso{
\code{\link{print.tidytensor}}, \code{\link{ranknames}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ranknames_various.R
\name{ranknames<-}
\alias{ranknames<-}
\title{Assign ranknames to a tidytensor.}
\usage{
ranknames(x) <- value
}
\arguments{
\item{x}{input tidytensor to set ranknames on.}

\item{value}{what to store in ranknames(x).}
}
\description{
A tidytensor t may have ranknames(t); this is a character vector of the same length as dim(t)
for future use. Note that ranknames(t) is independent of names(t) or dimnames(t); we are not naming elements,
or the dimension names for each rank, but rank names themselves.
Like names() and dimnames(), unset ranknames() are NULL.
}
\details{
Ranknames for a tidytensor t are stored as the names() attribute of dimnames(t). If dimnames(t) happens
to be null, before setting ranknames() we create valid dimnames() filled with NA values. The tidytensor
package also provides a specialized dimnames() which preserves ranknames when setting dimnames().
}
\examples{
t <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t) <- c("sample", "row", "col")
print(t)

# works like names():
t <- as.tidytensor(array(1:(3 * 4 * 5), dim = c(3, 4, 5)))
ranknames(t) <- c("sample", "row", "col")
print(ranknames(t))
ranknames(t)[3] <- "pixel"
print(t)
}
\seealso{
\code{\link{set_ranknames}}, \code{\link{dimnames<-}}
}
