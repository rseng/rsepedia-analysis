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
---
title: pixelclasser
author: Carlos Real
date: January 24, 2021
output:
  md_document:
    variant: gfm
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# pixelclasser

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![rOpenSci peer-review](https://badges.ropensci.org/406_status.svg)](https://github.com/ropensci/software-review/issues/406)
<!-- badges: end -->

This package contains a set of tools to classify the pixels of digital images
into colour categories arbitrarily defined by the user. It is a simple version
of the multivariate technique known as Support Vector Machine, adapted to this
particular use.

The procedure is simple. A digital image in JPEG or TIFF format is imported into
R. The original image contains three colour variables (or bands): $R$, $G$, and
$B$. The first step is to transform them into proportions ($r$, $g$ and $b$),
which simplifies the problem into a bivariate one. The pixels of the test images
can then be represented in the plane defined by two of the variables (the user
judges which two are more convenient by trial and error) and, hopefully, they
would form separate clusters (pixel categories). The user then traces straight
lines (classification rules) that enclose the pixel clusters. Using the
mathematical expression for these rules and the values of the transformed
variables, each pixel can be classified in one category.  This produces a set of
logical matrices (incidence matrices) indicating which pixels belong to each
category, stored in appropriate R objects. These can be submitted to posterior
analysis or used to create a new version of the original image showing the
category of each pixel.

`pixelclasser` contains functions to visualize the pixels of the images and the
rules created by the user, to create the rules and to store them in objects that
can be passed to function `classify_pixels()` for the analysis of the image, and
functions to import and export the original and the classified images.

## Installation

You can install the last version from the rOpenSci repository in  GitHub using
packages `remotes` or `devtools`, which install `remotes`

```{r eval=FALSE}
remotes::install_github("ropensci/pixelclasser", build_vignettes = TRUE)
devtools::install_github("ropensci/pixelclasser", build_vignettes = TRUE)
```

## Using pixelclasser

The manual with the description of each function and use examples is the file
`/doc/pixelclasser_1.0.0.pdf` (see the link to source code on the right), but
its contents can be found in the Reference section of this website.

An example session is described in the vignette included in the package, which
can be accessed after installation in the usual way:

```{r eval=FALSE}
vignette("pixelclasser")
```
It also can be accessed in the section Get started in the top menu of this page.

# Code of conduct

Please note that this package is released with a [Contributor
Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.
---
title: "A pixelclasser sample session"
output: [html_document]
fig.height: 5
fig.width: 5
vignette: >
  %\VignetteIndexEntry{A pixelclasser sample session}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction #

This package contains a set of tools to classify the pixels of digital images into colour categories arbitrarily defined by the user. It contains functions to

* visualize the distribution of the pixel colours in the images,
* define classification rules
* classify the pixels and to store this information in R objects,
* save these as image files.

It is a simple version of the multivariate technique known as Support Vector Machine (Cortes and Vapnik, 1995; Bennet and Campbell, 2000), adapted to this particular use. A manuscript describing the package `pixelclasser` and its use in real research has been submitted to Methods in Ecology and Evolution (Real et al., 2021). It also describes the procedure in more detail than the next paragraphs. 

### The procedure ###

The basic steps of the procedure are the following:

* One or more digital images in JPEG or TIFF format is imported into R. The categories to identify are represented in this set (the test set).
* The values of the three three colour variables (or bands) that compose each image (*R*, *G*, and *B*) are transformed into proportions (*r*, *g* and *b*).
* The pixels of the image are plotted in the plane defined by two of the transformed variables (the user can select them arbitrarily) and, hopefully, they would form separate clusters (pixel categories).
* The user then traces straight lines that separate the pixel clusters. Using the mathematical expression for these rules and the *rgb* values, each pixel can be tested for membership in each category (see below).
* Recording the results of the tests as 1 or 0 (pass/fail), an incidence matrix is build for that rule. This is the result of the procedure, which can be submitted to posterior analysis or used to create a new version of the original image showing the category of each pixel.

The second step simplifies the problem because it makes one of the variables dependent on the other two (as *r + g + b* = 1). Moreover, the transformation eliminates colour variations due to differences in illumination.

The expressions for classification rules are the same as the expression for a straight line but using one of the comparison operators $<$, $\leq$, $>$ or $\geq$. For example: $r \geq a g +c$, being $a$ and $c$ the slope and intercept of the line, and $r$ and $g$ the colour variables selected for the classification. A single line can produce two classification rules.

### Using several rules per category ###

When there are more than two categories, or when the cluster of points has a complex shape, a single rule is not enough. In these cases the procedure has additional steps:

* several rules are defined for each category,
* incidence matrices are created for each rule,
* the incidence matrices are combined with the `&` operator to obtain the category incidence matrix.

The last step is equivalent to estimate the union of the incidence matrices, i e $\mathbf{M} = \mathbf{M}_{1} \cap \mathbf{M}_{2} \cap \ldots \cap \mathbf{M}_{p}$, being *p* the number of rules.

### Concave category shapes ###

A caveat of the method is that the rules must delimit a convex polygon to combine the individual rule results successfully (in a convex polygon, a line joining any two internal points is contained in the polygon). Not all clusters have convex shape. In these cases, the cluster must be divided in convex sub-polygons (subcategories) for which rules are defined as before. The incidence matrices of the subcategories are combined using the `|` operator, i.e. $\mathbf{M} = \mathbf{M}_{1} \cup \mathbf{M}_{2} \cup \ldots \cup \mathbf{M}_{s}$, being *s* the number of subcategories. Note that any polygon, convex or not, can be subdivided in triangles and, as triangles are convex polygons, it is always possible to solve this problem. Note that the goal is to obtain a minimal set of convex polygons, not a complete triangulation. The example presented below is one of such cases.

# The session #

What follows is a sample session illustrating both the method and the use of the package functions. It uses an example image and a test set created by cutting small areas out of the example image. It is not a good test set, see below, but it is enough to show how the method works, and its problems.

### Loading the functions ###

The package is loaded in the usual way:

```{r}
library(pixelclasser)
```
  
### Image loading and transforming ###

These are the images included in the package as examples. The goal is to identify the pixels corresponding to dead, oak and ivy leaves that compose the image. The small images are fragments of the main image and are the test set. In an actual case, more than one image per category should be used to represent the whole variation of the category:

```{r echo=FALSE, fig.align='center', out.width = "50%"}
knitr::include_graphics('../inst/extdata/ExampleImages.png')
```

As the images are included in the package as external (non R) data, they are loaded with the following code:

```{r}
ivy_oak_rgb <- read_image(system.file("extdata", "IvyOak400x300.JPG", package = "pixelclasser"))
test_ivy_rgb <- read_image(system.file("extdata", "TestIvy.JPG", package = "pixelclasser"))
test_oak_rgb <- read_image(system.file("extdata", "TestOak.JPG", package = "pixelclasser"))
test_dead_rgb <- read_image(system.file("extdata", "TestDeadLeaves.JPG", package = "pixelclasser"))
```

The function `read_image()`  performs the first step of the procedure. It stores the image as an array of *rgb* values, which are the proportion of each colour variable (i.e. *R /(R+G+B)*, and so on). This uses functions from packages `jpeg` or `tiff`, and uses the extension in the file name to identify which one to use.

### Pixel distributions in *rgb* space ###

Before plotting pixels and lines, it is convenient to define a set of colours to use throughout the session:
  
```{r}
transparent_black <- "#00000008"
brown <- "#c86432ff"
yellow <- "#ffcd0eff"
blue <- "#5536ffff"
green <- "#559800ff"
```

The next step is to visualize the distribution of the pixels in *rgb* space, but only two variables are needed. Any pair of variables would do, but a particular combination might produce a better display of the clusters. It is a matter of try the three possible combinations to select the most convenient.

Plotting the pixels is a two-step procedure: a void plot is drawn first and then the pixels are added to the plot (the use of a transparent black colour, `#00000008`, creates a "density plot" effect):
  
```{r, out.width = "50%", fig.align="center", out.width = "50%"}
plot_rgb_plane("r", "b", main = "Image: ivy and oak")
plot_pixels(ivy_oak_rgb, "r", "b", col = transparent_black)
```

The coloured lines are an aid to interpret the graph: no pixels could be found outside the blue lines, and the red lines converging in the barycentre of the triangle *(r, g, b)* = (1/3, 1/3, 1/3), define the areas where a colour is dominant. Note that graphical parameters (`main` in this example) can be passed to the function to change the final appearance of the graph. All the auxiliary lines can be deleted, as in the following example, which uses different colour variables to create the graph.

```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("r", "g", plot_limits = F, plot_guides = F, plot_grid = F)
plot_pixels(ivy_oak_rgb, "r", "g", col = transparent_black)
```

There are two clear pixel clusters and a small, but noticeable, quantity of pixels in between. Also visible are linear patterns that are artefacts created because the *RGB* data are discrete variables (eight bit in the most common cases). These are more appreciable in the following graphs, which are restricted to the area occupied by the pixels. In the following examples *g* and *b* will be used as variables *x* and *y* for plotting and pixel classification. 

### Adding the pixels of the test images ###

The following code plots the pixels of the example image on the *gb* plane and then adds the pixels of the test images, using arbitrary colours. Here, the graphic parameters `xlim` and `ylim` were used to limit the extent of the plot to the area occupied by the pixels:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33))
plot_pixels(ivy_oak_rgb, "g", "b", col = transparent_black)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)
plot_pixels(test_dead_rgb, "g", "b", col = brown)
```

The plot shows that the clusters of pixels in the `ivy_oak_rgb` image correspond to dead leaves (on the left), and oak and ivy (on the right).

The small areas taken as test images were not representative of the whole pixel set in the image, as they do not cover the same area as the black pixels. This is not a surprise given that a single sample was collected for each type of pixel.

Warning: plotting several million points in an R graph is an slow process. Be patient or use images as small as possible. Using a nice smartphone with a petapixel camera sensor to capture images is good for artistic purposes, but not always for efficient scientific work.

### Defining the rules ###

Defining the rules that classify the pixels is a matter of tracing straight lines to separate the clusters. In this example, a single line more or less equidistant to both clusters should suffice to separate them. The intermediate points will be arbitrarily ascribed to one category.

The rules are defined by setting the name of the rule, the colour variables to use, the coordinates of two points in the plane and a comparison operator. The exact placement of the line is an arbitrary decision, as the method does not include any mechanism to place it automatically.

There are two methods to create the rule. The first uses the function `create_rule()`, which receives a list with the coordinates of two points defining a line in the selected subspace. In the following example, the points with coordinates (*g*, *b*) = (0.345, 1/3), and *(g,b)* = (0.40, 0.10) defined the position of the first line, and were selected by trial and error. The adequate operator must be included in the rule definition:
  
```{r}
rule_01 <- define_rule("rule_01", "g", "b", list(c(0.345, 1/3), c(0.40, 0.10)), "<")
rule_02 <- define_rule("rule_02", "g", "b", list(c(0.345, 1/3), c(0.40, 0.10)), ">=")
```

Both rules are described by the same line but use different comparison operator. `rule_01` includes the pixels at the left (under) of the line and `rule_02` those at the right (over) and on the line, i.e. the dead leaves and the fresh leaves, respectively. Each line can generate two rules, but beware: if `>` and `<` define the rules, then the points on the line will not belong to any category, and if `>=` and `<=` are used, the points on the line will belong to the two categories simultaneously. The function that classifies the pixels can identify the second type of error, but if there are legitimate unclassified points, the first type can pass unnoticed.

The second method uses `place_rule()`, which is a wrapper for `graphics:locator()` that allows the user to select the two points by clicking in the rgb plot with the mouse:

```
rp01 <- place_rule("g", "b")
rp02 <- place_rule("g", "b", "v")
rule_07 <- define_rule("rule_07", "g", "b", rp02, ">=")
```
The function returns an object of class `rule_points` that can then be passed to `define_rule()` in the parameter `rule_points`. The second example produces a vertical line ("h" for horizontal lines), which would be difficult to produce by hand. To make the code self-contained, `create_rule()` is used in this vignette, but using `place_rule()` is the easiest way to define the rules. It is even easier to place the call to `place_rule()` in the call to `create_rule()` to avoid creating the intermediate `rule_points` object:

```
rule_07 <- define_rule("rule_07", "g", "b", place_rule("g", "b"), "<")
```

Note that both `define_rule()` and the `rule_points` object must use the same colour variables as axis. `define_rule()` throws an error if this condition does not hold.

The rule objects store the values passed as parameters, the parameters of the equation of the line (*a* and *c*), and a textual representation of the equation which will be evaluated by the classification function. To check the correctness of the rules, the line can be added to the plot:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33))
plot_pixels(ivy_oak_rgb, "g", "b", col = transparent_black)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)
plot_pixels(test_dead_rgb, "g", "b", col = brown)
plot_rule(rule_01, lty = 2, col = brown)
```

In order to classify the fresh leaves into ivy and oak categories, more rules are needed. The pixels of the oak test image were plotted and used to define additional rules:

```{r}
rule_03 <- define_rule("rule_03","g", "b", list(c(0.35, 0.30), c(0.565, 0.10)), "<")
rule_04 <- define_rule("rule_04","g", "b", list(c(0.35, 0.25), c(0.5, 0.25)), "<")
```

Here is the plot of the pixels and the rules. Line type and colour were set using the graphical parameters `lty` and `col` (see `graphics::par)`) `...` argument of `plot_rule()`:

```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33), plot_limits = F, plot_guides = F)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_rule(rule_01, lty = 2, col = green)
plot_rule(rule_03, lty = 2, col = green)
plot_rule(rule_04, lty = 2, col = green)
```

The ivy pixels are now plotted to check whether the rules can identify them. Labels to identify the lines and their associated rules were added to the plot, using the parameter `shift` to place them conveniently:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2,0.6), ylim=c(0.1,0.33), plot_limits = F, plot_guides = F)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)

plot_rule(rule_02, lty = 1, col = green)
label_rule(rule_02, label = expression('L'[1]*' (R'[1]*',R'[2]*')'), shift = c(0.035, -0.004), col = green)

plot_rule(rule_03, lty = 1, col = green)
label_rule(rule_03, label = expression('L'[2]*' (R'[3]*',R'[5]*')'), shift = c(0.20, -0.15), col = green)

plot_rule(rule_04, lty = 1, col = green)
label_rule(rule_04, label = expression('L'[3]*' (R'[4]*',R'[6]*')'), shift = c(0.19, 0.0), col = green)
```

The graph shows two problems: a) part of the ivy pixels are inside the area delimited by the oak rules, i.e. both categories overlap. As a consequence, some ivy pixels will be miss-classified. b) The shape of the ivy cluster is not convex.

To solve the second problem, two subcategories must be defined as explained before. The first is delimited by *L~1~* and *L~3~*, and the second by *L~2~* and *L~3~*. To do this, two new rules are needed:
  
```{r}
rule_05 <- define_rule("rule_05", "g", "b", list(c(0.35, 0.30), c(0.565, 0.16)), ">=")
rule_06 <- define_rule("rule_06", "g", "b", list(c(0.35, 0.25), c(0.5, 0.25)), ">=")
```

There is an intentional error in the coordinates of the second point of `rule_05`, which are not the same as in `rule_03`. It is left here to show later how the internal checks of the classification function allow to detect it.

Note that because no points can be found outside the blue triangle, its borders are implicit rules that close the polygons defined by the explicit rules, but that do not need to be created.

### Creating the classifier objects ###

After the rules have been defined, they must be included in classifier objects which will be used later by `classify_pixels()`. This function receives a list of objects of class `pixel_cat`, each containing the information needed to identify the pixels belonging to a particular category. These objects contain a list of objects `pixel_subcat`, each containing a list of objects of class `pixel_rule`. This is a nested structure which always has three levels (rule, subcategory, and category) even when no subcategories would be needed for the classification. In these simple cases, a subcategory object containing the rules is internally added to the category object. This consistency in the structure of the objects simplifies the code of `classify_pixels()`.

Creating the classifiers is simple once the rules have been defined. The following code defines a class classifier that can identify the dead leaves:
  
```{r}
cat_dead_leaves <- define_cat("dead_leaves", blue, rule_01)
```

`define_cat()` needs a label for the category, a colour to identify the pixels if an image file is generated, and a list of rules that define the category. Here, the list contains a single rule. This is a simple case where no subcategories are needed and a list of rules suffice to classify the pixels. See below for a more complex case. The corresponding classifier for the living leaves is:

```{r}
cat_living_leaves <- define_cat("living_leaves", yellow, rule_02)
```

In these examples `define_cat()` detects that the the list contains only rules, not subcategories, and wraps them into an object of type `pixel_subcat`.

A classifier object for oak pixels needs three rules:

```{r}
cat_oak_leaves <- define_cat("oak_leaves", green, rule_02, rule_03, rule_04)
```

Finally, the classifier for ivy pixels is the most complex, as it is composed of two subcategory objects that must be defined explicitly and then included in the class classifier:

```{r}
subcat_ivy01 <- define_subcat("ivy01", rule_02, rule_06)
subcat_ivy02 <- define_subcat("ivy02", rule_04, rule_05)

cat_ivy_leaves <- define_cat("ivy_leaves", yellow, subcat_ivy01, subcat_ivy02)
```

Note that rules and subcategories cannot be mixed in the list, so sometimes a subcategory object containing a single rule should be created by the user before creating the category object. `define_cat()` checks for the type of the objects in the list and complains if they are not adequate.

### Classifying the pixels ###

Function `classify_pixels()` uses a list of categories to classify the pixels. As a preliminary example, the example image will be classified in dead and living leaves. The parameters are the object to classify and the list of category objects:
  
```{r}
dead_live_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_living_leaves)
```

Note that a category named `unclassified` is automatically added to the classes defined by the user. With the rule set used in this example, the `unclassified` class must contain zero pixels. The function outputs counts of pixels in each classes, which are useful to verify the consistency of the rules. The function also detects duplicate pixels, i e those counted in more than one class (the sum of the pixels in each class is larger than the total number of pixels). If the consistency of the rules has been verified, these messages can be suppressed by `verbose = FALSE` in the function call.

The result can be saved as a JPEG (or TIFF) file:

```{r, eval=FALSE}
save_classif_image(dead_live_classified, "DeadLiveClassified.JPG", quality = 1)
```

The type of the file is automatically selected from the file name (only `JPEG` or `TIFF` files allowed). Note the use of the `quality` parameter, which is passed to the underlying function, to set the quality of the JPEG file produced to its maximum value.

The final classification includes the three categories:

```{r}
ivy_oak_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_ivy_leaves, cat_oak_leaves)
```

The function informs that several points were left unclassified. This is the consequence of the error in the definition of `rule_05` noted above. If the error is corrected, and the image classified again:

```{r}
rule_05 <- define_rule("rule_05", "g", "b", list(c(0.35, 0.30), c(0.565, 0.10)), ">=")
subcat_ivy02 <- define_subcat("ivy02", rule_04, rule_05)
cat_ivy_leaves <- define_cat("ivy_leaves", yellow, subcat_ivy01, subcat_ivy02)
ivy_oak_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_ivy_leaves, cat_oak_leaves)
```

the result is correct and it can be saved, as a `TIFF` file in this case:

```{r, eval = FALSE}
save_classif_image(ivy_oak_classified, "IvyOakClassified.TIFF")
```

The following figure shows the original image and the results of the two classifications.

```{r echo=FALSE, fig.align='center'}
knitr::include_graphics('../inst/extdata/ClassifResults.png')
```

Dead and fresh leaves were correctly differentiated in the first classification. The second classification was accurate for dead and oak pixels but, as expected, part of the ivy pixels were miss-classified as oak pixels because of the overlap between these two categories.

# References #

Bennet, K. P. and C. Campbell (2000). Support vector machines: hype or Halleluiah. SIGKDD Explorations 2, 1–11.

Cortes, C. and V. Vapnik (1995). Support-vector networks. Machine Learning 20, 273–297.

Real, C., Cruz de Carvalho, R., García-Seoane, R., Branquinho, C. and Varela, Z. (2021). A simplified support vector machine method for image classification. Submitted to Methods in Ecology and Evolution.---
title: "A pixelclasser sample session"
output: [pdf_document]
fig.height: 5
fig.width: 5
vignette: >
  %\VignetteIndexEntry{A pixelclasser sample session}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction #

This package contains a set of tools to classify the pixels of digital images into colour categories arbitrarily defined by the user. It contains functions to

* visualize the distribution of the pixel colours in the images,
* define classification rules
* classify the pixels and to store this information in R objects,
* save these as image files.

It is a simple version of the multivariate technique known as Support Vector Machine (Cortes and Vapnik, 1995; Bennet and Campbell, 2000), adapted to this particular use. A manuscript describing the package `pixelclasser` and its use in real research has been submitted to Methods in Ecology and Evolution (Real et al., 2021). It also describes the procedure in more detail than the next paragraphs. 

### The procedure ###

The basic steps of the procedure are the following:

* One or more digital images in JPEG or TIFF format is imported into R. The categories to identify are represented in this set (the test set).
* The values of the three three colour variables (or bands) that compose each image (*R*, *G*, and *B*) are transformed into proportions (*r*, *g* and *b*).
* The pixels of the image are plotted in the plane defined by two of the transformed variables (the user can select them arbitrarily) and, hopefully, they would form separate clusters (pixel categories).
* The user then traces straight lines that separate the pixel clusters. Using the mathematical expression for these rules and the *rgb* values, each pixel can be tested for membership in each category (see below).
* Recording the results of the tests as 1 or 0 (pass/fail), an incidence matrix is build for that rule. This is the result of the procedure, which can be submitted to posterior analysis or used to create a new version of the original image showing the category of each pixel.

The second step simplifies the problem because it makes one of the variables dependent on the other two (as *r + g + b* = 1). Moreover, the transformation eliminates colour variations due to differences in illumination.

The expressions for classification rules are the same as the expression for a straight line but using one of the comparison operators $<$, $\leq$, $>$ or $\geq$. For example: $r \geq a g +c$, being $a$ and $c$ the slope and intercept of the line, and $r$ and $g$ the colour variables selected for the classification. A single line can produce two classification rules.

### Using several rules per category ###

When there are more than two categories, or when the cluster of points has a complex shape, a single rule is not enough. In these cases the procedure has additional steps:

* several rules are defined for each category,
* incidence matrices are created for each rule,
* the incidence matrices are combined with the `&` operator to obtain the category incidence matrix.

The last step is equivalent to estimate the union of the incidence matrices, i e $\mathbf{M} = \mathbf{M}_{1} \cap \mathbf{M}_{2} \cap \ldots \cap \mathbf{M}_{p}$, being *p* the number of rules.

### Concave category shapes ###

A caveat of the method is that the rules must delimit a convex polygon to combine the individual rule results successfully (in a convex polygon, a line joining any two internal points is contained in the polygon). Not all clusters have convex shape. In these cases, the cluster must be divided in convex sub-polygons (subcategories) for which rules are defined as before. The incidence matrices of the subcategories are combined using the `|` operator, i.e. $\mathbf{M} = \mathbf{M}_{1} \cup \mathbf{M}_{2} \cup \ldots \cup \mathbf{M}_{s}$, being *s* the number of subcategories. Note that any polygon, convex or not, can be subdivided in triangles and, as triangles are convex polygons, it is always possible to solve this problem. Note that the goal is to obtain a minimal set of convex polygons, not a complete triangulation. The example presented below is one of such cases.

# The session #

What follows is a sample session illustrating both the method and the use of the package functions. It uses an example image and a test set created by cutting small areas out of the example image. It is not a good test set, see below, but it is enough to show how the method works, and its problems.

### Loading the functions ###

The package is loaded in the usual way:

```{r}
library(pixelclasser)
```
  
### Image loading and transforming ###

These are the images included in the package as examples. The goal is to identify the pixels corresponding to dead, oak and ivy leaves that compose the image. The small images are fragments of the main image and are the test set. In an actual case, more than one image per category should be used to represent the whole variation of the category:

```{r echo=FALSE, fig.align='center', out.width = "50%"}
knitr::include_graphics('../inst/extdata/ExampleImages.png')
```

As the images are included in the package as external (non R) data, they are loaded with the following code:

```{r}
ivy_oak_rgb <- read_image(system.file("extdata", "IvyOak400x300.JPG", package = "pixelclasser"))
test_ivy_rgb <- read_image(system.file("extdata", "TestIvy.JPG", package = "pixelclasser"))
test_oak_rgb <- read_image(system.file("extdata", "TestOak.JPG", package = "pixelclasser"))
test_dead_rgb <- read_image(system.file("extdata", "TestDeadLeaves.JPG", package = "pixelclasser"))
```

The function `read_image()`  performs the first step of the procedure. It stores the image as an array of *rgb* values, which are the proportion of each colour variable (i.e. *R /(R+G+B)*, and so on). This uses functions from packages `jpeg` or `tiff`, and uses the extension in the file name to identify which one to use.

### Pixel distributions in *rgb* space ###

Before plotting pixels and lines, it is convenient to define a set of colours to use throughout the session:
  
```{r}
transparent_black <- "#00000008"
brown <- "#c86432ff"
yellow <- "#ffcd0eff"
blue <- "#5536ffff"
green <- "#559800ff"
```

The next step is to visualize the distribution of the pixels in *rgb* space, but only two variables are needed. Any pair of variables would do, but a particular combination might produce a better display of the clusters. It is a matter of try the three possible combinations to select the most convenient.

Plotting the pixels is a two-step procedure: a void plot is drawn first and then the pixels are added to the plot (the use of a transparent black colour, `#00000008`, creates a "density plot" effect):
  
```{r, out.width = "50%", fig.align="center", out.width = "50%"}
plot_rgb_plane("r", "b", main = "Image: ivy and oak")
plot_pixels(ivy_oak_rgb, "r", "b", col = transparent_black)
```

The coloured lines are an aid to interpret the graph: no pixels could be found outside the blue lines, and the red lines converging in the barycentre of the triangle *(r, g, b)* = (1/3, 1/3, 1/3), define the areas where a colour is dominant. Note that graphical parameters (`main` in this example) can be passed to the function to change the final appearance of the graph. All the auxiliary lines can be deleted, as in the following example, which uses different colour variables to create the graph.

```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("r", "g", plot_limits = F, plot_guides = F, plot_grid = F)
plot_pixels(ivy_oak_rgb, "r", "g", col = transparent_black)
```

There are two clear pixel clusters and a small, but noticeable, quantity of pixels in between. Also visible are linear patterns that are artefacts created because the *RGB* data are discrete variables (eight bit in the most common cases). These are more appreciable in the following graphs, which are restricted to the area occupied by the pixels. In the following examples *g* and *b* will be used as variables *x* and *y* for plotting and pixel classification. 

### Adding the pixels of the test images ###

The following code plots the pixels of the example image on the *gb* plane and then adds the pixels of the test images, using arbitrary colours. Here, the graphic parameters `xlim` and `ylim` were used to limit the extent of the plot to the area occupied by the pixels:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33))
plot_pixels(ivy_oak_rgb, "g", "b", col = transparent_black)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)
plot_pixels(test_dead_rgb, "g", "b", col = brown)
```

The plot shows that the clusters of pixels in the `ivy_oak_rgb` image correspond to dead leaves (on the left), and oak and ivy (on the right).

The small areas taken as test images were not representative of the whole pixel set in the image, as they do not cover the same area as the black pixels. This is not a surprise given that a single sample was collected for each type of pixel.

Warning: plotting several million points in an R graph is an slow process. Be patient or use images as small as possible. Using a nice smartphone with a petapixel camera sensor to capture images is good for artistic purposes, but not always for efficient scientific work.

### Defining the rules ###

Defining the rules that classify the pixels is a matter of tracing straight lines to separate the clusters. In this example, a single line more or less equidistant to both clusters should suffice to separate them. The intermediate points will be arbitrarily ascribed to one category.

The rules are defined by setting the name of the rule, the colour variables to use, the coordinates of two points in the plane and a comparison operator. The exact placement of the line is an arbitrary decision, as the method does not include any mechanism to place it automatically.

There are two methods to create the rule. The first uses the function `create_rule()`, which receives a list with the coordinates of two points defining a line in the selected subspace. In the following example, the points with coordinates (*g*, *b*) = (0.345, 1/3), and *(g,b)* = (0.40, 0.10) defined the position of the first line, and were selected by trial and error. The adequate operator must be included in the rule definition:
  
```{r}
rule_01 <- define_rule("rule_01", "g", "b", list(c(0.345, 1/3), c(0.40, 0.10)), "<")
rule_02 <- define_rule("rule_02", "g", "b", list(c(0.345, 1/3), c(0.40, 0.10)), ">=")
```

Both rules are described by the same line but use different comparison operator. `rule_01` includes the pixels at the left (under) of the line and `rule_02` those at the right (over) and on the line, i.e. the dead leaves and the fresh leaves, respectively. Each line can generate two rules, but beware: if `>` and `<` define the rules, then the points on the line will not belong to any category, and if `>=` and `<=` are used, the points on the line will belong to the two categories simultaneously. The function that classifies the pixels can identify the second type of error, but if there are legitimate unclassified points, the first type can pass unnoticed.

The second method uses `place_rule()`, which is a wrapper for `graphics:locator()` that allows the user to select the two points by clicking in the rgb plot with the mouse:

```
rp01 <- place_rule("g", "b")
rp02 <- place_rule("g", "b", "v")
rule_07 <- define_rule("rule_07", "g", "b", rp02, ">=")
```
The function returns an object of class `rule_points` that can then be passed to `define_rule()` in the parameter `rule_points`. The second example produces a vertical line ("h" for horizontal lines), which would be difficult to produce by hand. To make the code self-contained, `create_rule()` is used in this vignette, but using `place_rule()` is the easiest way to define the rules. It is even easier to place the call to `place_rule()` in the call to `create_rule()` to avoid creating the intermediate `rule_points` object:

```
rule_07 <- define_rule("rule_07", "g", "b", place_rule("g", "b"), "<")
```

Note that both `define_rule()` and the `rule_points` object must use the same colour variables as axis. `define_rule()` throws an error if this condition does not hold.

The rule objects store the values passed as parameters, the parameters of the equation of the line (*a* and *c*), and a textual representation of the equation which will be evaluated by the classification function. To check the correctness of the rules, the line can be added to the plot:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33))
plot_pixels(ivy_oak_rgb, "g", "b", col = transparent_black)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)
plot_pixels(test_dead_rgb, "g", "b", col = brown)
plot_rule(rule_01, lty = 2, col = brown)
```

In order to classify the fresh leaves into ivy and oak categories, more rules are needed. The pixels of the oak test image were plotted and used to define additional rules:

```{r}
rule_03 <- define_rule("rule_03","g", "b", list(c(0.35, 0.30), c(0.565, 0.10)), "<")
rule_04 <- define_rule("rule_04","g", "b", list(c(0.35, 0.25), c(0.5, 0.25)), "<")
```

Here is the plot of the pixels and the rules. Line type and colour were set using the graphical parameters `lty` and `col` (see `graphics::par)`) `...` argument of `plot_rule()`:

```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2, 0.6), ylim = c(0.1, 0.33), plot_limits = F, plot_guides = F)
plot_pixels(test_oak_rgb, "g", "b", col = green)
plot_rule(rule_01, lty = 2, col = green)
plot_rule(rule_03, lty = 2, col = green)
plot_rule(rule_04, lty = 2, col = green)
```

The ivy pixels are now plotted to check whether the rules can identify them. Labels to identify the lines and their associated rules were added to the plot, using the parameter `shift` to place them conveniently:
  
```{r, fig.align="center", out.width = "50%"}
plot_rgb_plane("g", "b", xlim = c(0.2,0.6), ylim=c(0.1,0.33), plot_limits = F, plot_guides = F)
plot_pixels(test_ivy_rgb, "g", "b", col = blue)

plot_rule(rule_02, lty = 1, col = green)
label_rule(rule_02, label = expression('L'[1]*' (R'[1]*',R'[2]*')'), shift = c(0.035, -0.004), col = green)

plot_rule(rule_03, lty = 1, col = green)
label_rule(rule_03, label = expression('L'[2]*' (R'[3]*',R'[5]*')'), shift = c(0.20, -0.15), col = green)

plot_rule(rule_04, lty = 1, col = green)
label_rule(rule_04, label = expression('L'[3]*' (R'[4]*',R'[6]*')'), shift = c(0.19, 0.0), col = green)
```

The graph shows two problems: a) part of the ivy pixels are inside the area delimited by the oak rules, i.e. both categories overlap. As a consequence, some ivy pixels will be miss-classified. b) The shape of the ivy cluster is not convex.

To solve the second problem, two subcategories must be defined as explained before. The first is delimited by *L~1~* and *L~3~*, and the second by *L~2~* and *L~3~*. To do this, two new rules are needed:
  
```{r}
rule_05 <- define_rule("rule_05", "g", "b", list(c(0.35, 0.30), c(0.565, 0.16)), ">=")
rule_06 <- define_rule("rule_06", "g", "b", list(c(0.35, 0.25), c(0.5, 0.25)), ">=")
```

There is an intentional error in the coordinates of the second point of `rule_05`, which are not the same as in `rule_03`. It is left here to show later how the internal checks of the classification function allow to detect it.

Note that because no points can be found outside the blue triangle, its borders are implicit rules that close the polygons defined by the explicit rules, but that do not need to be created.

### Creating the classifier objects ###

After the rules have been defined, they must be included in classifier objects which will be used later by `classify_pixels()`. This function receives a list of objects of class `pixel_cat`, each containing the information needed to identify the pixels belonging to a particular category. These objects contain a list of objects `pixel_subcat`, each containing a list of objects of class `pixel_rule`. This is a nested structure which always has three levels (rule, subcategory, and category) even when no subcategories would be needed for the classification. In these simple cases, a subcategory object containing the rules is internally added to the category object. This consistency in the structure of the objects simplifies the code of `classify_pixels()`.

Creating the classifiers is simple once the rules have been defined. The following code defines a class classifier that can identify the dead leaves:
  
```{r}
cat_dead_leaves <- define_cat("dead_leaves", blue, rule_01)
```

`define_cat()` needs a label for the category, a colour to identify the pixels if an image file is generated, and a list of rules that define the category. Here, the list contains a single rule. This is a simple case where no subcategories are needed and a list of rules suffice to classify the pixels. See below for a more complex case. The corresponding classifier for the living leaves is:

```{r}
cat_living_leaves <- define_cat("living_leaves", yellow, rule_02)
```

In these examples `define_cat()` detects that the the list contains only rules, not subcategories, and wraps them into an object of type `pixel_subcat`.

A classifier object for oak pixels needs three rules:

```{r}
cat_oak_leaves <- define_cat("oak_leaves", green, rule_02, rule_03, rule_04)
```

Finally, the classifier for ivy pixels is the most complex, as it is composed of two subcategory objects that must be defined explicitly and then included in the class classifier:

```{r}
subcat_ivy01 <- define_subcat("ivy01", rule_02, rule_06)
subcat_ivy02 <- define_subcat("ivy02", rule_04, rule_05)

cat_ivy_leaves <- define_cat("ivy_leaves", yellow, subcat_ivy01, subcat_ivy02)
```

Note that rules and subcategories cannot be mixed in the list, so sometimes a subcategory object containing a single rule should be created by the user before creating the category object. `define_cat()` checks for the type of the objects in the list and complains if they are not adequate.

### Classifying the pixels ###

Function `classify_pixels()` uses a list of categories to classify the pixels. As a preliminary example, the example image will be classified in dead and living leaves. The parameters are the object to classify and the list of category objects:
  
```{r}
dead_live_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_living_leaves)
```

Note that a category named `unclassified` is automatically added to the classes defined by the user. With the rule set used in this example, the `unclassified` class must contain zero pixels. The function outputs counts of pixels in each classes, which are useful to verify the consistency of the rules. The function also detects duplicate pixels, i e those counted in more than one class (the sum of the pixels in each class is larger than the total number of pixels). If the consistency of the rules has been verified, these messages can be suppressed by `verbose = FALSE` in the function call.

The result can be saved as a JPEG (or TIFF) file:

```{r, eval=FALSE}
save_classif_image(dead_live_classified, "DeadLiveClassified.JPG", quality = 1)
```

The type of the file is automatically selected from the file name (only `JPEG` or `TIFF` files allowed). Note the use of the `quality` parameter, which is passed to the underlying function, to set the quality of the JPEG file produced to its maximum value.

The final classification includes the three categories:

```{r}
ivy_oak_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_ivy_leaves, cat_oak_leaves)
```

The function informs that several points were left unclassified. This is the consequence of the error in the definition of `rule_05` noted above. If the error is corrected, and the image classified again:

```{r}
rule_05 <- define_rule("rule_05", "g", "b", list(c(0.35, 0.30), c(0.565, 0.10)), ">=")
subcat_ivy02 <- define_subcat("ivy02", rule_04, rule_05)
cat_ivy_leaves <- define_cat("ivy_leaves", yellow, subcat_ivy01, subcat_ivy02)
ivy_oak_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves, cat_ivy_leaves, cat_oak_leaves)
```

the result is correct and it can be saved, as a `TIFF` file in this case:

```{r, eval = FALSE}
save_classif_image(ivy_oak_classified, "IvyOakClassified.TIFF")
```

The following figure shows the original image and the results of the two classifications.

```{r echo=FALSE, fig.align='center'}
knitr::include_graphics('../inst/extdata/ClassifResults.png')
```

Dead and fresh leaves were correctly differentiated in the first classification. The second classification was accurate for dead and oak pixels but, as expected, part of the ivy pixels were miss-classified as oak pixels because of the overlap between these two categories.

# References #

Bennet, K. P. and C. Campbell (2000). Support vector machines: hype or Halleluiah. SIGKDD Explorations 2, 1–11.

Cortes, C. and V. Vapnik (1995). Support-vector networks. Machine Learning 20, 273–297.

Real, C., Cruz de Carvalho, R., García-Seoane, R., Branquinho, C. and Varela, Z. (2021). A simplified support vector machine method for image classification. Submitted to Methods in Ecology and Evolution.% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_clasif_image.R
\name{save_classif_image}
\alias{save_classif_image}
\title{Saves a classified image in TIFF or JPEG format}
\usage{
save_classif_image(classified_image, file_name, ...)
}
\arguments{
\item{classified_image}{an object of class \code{classified_image}.}

\item{file_name}{a character string with the name of the output file,
including the extension.}

\item{...}{further parameters to pass to functions \code{writeJPG} and
\code{writeTIFF}. If void, the default values of these functions are used.}
}
\value{
It does not return anything, only creates the file.
}
\description{
Creates an image file in TIFF or JPEG format from an array of class
\code{classified_image}.
}
\details{
The type of the output file (jpeg or tiff) is selected from the
  extension included in the file name. It must be one of \code{("jpg", "JPG",
  "jpeg", "JPEG", "tif", "TIF", "tiff", "TIFF")}.

  Note that the default value for jpg quality is 0.7. For maximal quality set
  \code{quality = 1} using the \dots argument. Such adjustments are not
  needed with \code{tiff} files, as this is a lossless format.
}
\examples{
\dontrun{

# Saving an hypothetical image. Note the use of quality to set the
# maximum quality level in the JPEG file
save_classif_image(image01_class, "./myimages/image01_classified.jpg",
                   quality = 1)
}

}
\seealso{
\code{\link{classify_pixels}}

  For more information about the options for file formatting see see the help
  pages of \code{\link[jpeg]{readJPEG}} and \code{\link[tiff]{readTIFF}}
  functions in packages \code{jpeg} and \code{tiff}, respectively.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_pixels.R
\name{classify_pixels}
\alias{classify_pixels}
\title{Classifies the pixels of an image}
\usage{
classify_pixels(image_prop, ..., unclassed_colour = "black", verbose = TRUE)
}
\arguments{
\item{image_prop}{an array containing the image. It must be an
object produced with function \code{read_image()}.}

\item{\dots}{a list of objects of class \code{pixel_cat} containing the
classification rules.}

\item{unclassed_colour}{a character string setting the colour to be assigned
to unclassified pixels. Defaults to "black".}

\item{verbose}{a logical value. When TRUE (default) the function prints some
statistics about the classification.}
}
\value{
Returns an object of class \code{classified_image}, which is a list
  containing nested lists. Each first-level element corresponds to one of the
  pixel categories and its name is the category name. They contains the
  second-level list, which have the following elements:
\itemize{
  \item \code{colour}: a matrix defining a colour to paint the pixels in the
  classified image. Inherited from the \code{pixel_class} object defining the
  class.
  \item \code{incid_mat}: a logical matrix where \code{TRUE} values indicate
  that the pixel belongs to this pixel category.
}
}
\description{
Classifies the pixels represented in an object of class
\code{transformed_image} using the rules contained in a list of objects of
class \code{pixel_cat}.
}
\details{
This function generates a set of incidence matrices indicating
  whether a pixel belongs to a pixel category or not. An additional matrix
  identifies the pixels that do not belong to the defined categories, i e
  unclassed pixels. Depending on how the rules were defined, it can be void
  or contain pixels, but it is always present and named \code{unclassified}.

  To create the incidence matrices for each category, a matrix for each rule
  is created and then combined with the matrices of the other using the
  \code{and} operator.

  When a set of subcategories is used, the procedure is the same for each
  subcategory and then the matrices of the subcategories are combined again,
  this time using the \code{or} operator. See the help for
  \code{define_subcat} for more details.

  \code{unclassed_colour} can be specified in any form understood by
  \code{grDevices::col2grb}.
}
\examples{

# The series of steps to classify a image supplied in the package

yellow <- "#ffcd0eff"
blue <- "#5536ffff"

ivy_oak_rgb <- read_image(system.file("extdata", "IvyOak400x300.JPG",
                          package = "pixelclasser"))

rule_01 <- define_rule("rule_01", "g", "b",
                       list(c(0.345, 1/3), c(0.40, 0.10)), comp_op = "<")
rule_02 <- define_rule("rule_02", "g", "b",
                       list(c(0.345, 1/3), c(0.40, 0.10)), comp_op = ">=")

cat_dead_leaves <- define_cat("dead_leaves", blue, rule_01)
cat_living_leaves <- define_cat("living_leaves", yellow, rule_02)

ivy_oak_classified <- classify_pixels(ivy_oak_rgb, cat_dead_leaves,
                        cat_living_leaves)

}
\seealso{
\code{\link{define_cat}}, \code{\link[grDevices]{col2rgb}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_subcat.R
\name{define_subcat}
\alias{define_subcat}
\title{Creates a subcategory object}
\usage{
define_subcat(subcat_name, ...)
}
\arguments{
\item{subcat_name}{a character string containing the name of the subcategory.}

\item{\dots}{a list of objects of class \code{pixel_rule}.}
}
\value{
An object of class \code{pixel_subcat}, which is a list with the
  following elements:
  \itemize{
  \item \code{name} a character string containing the name of the
  subcategory.
  \item \code{rules_list} a list of  \code{pixel_rule} objects.
  }
}
\description{
Creates an object of class \code{pixel_subcat} from a list of objects of
class \code{pixel_rule}.
}
\details{
When the shape of the cluster of pixels belonging to one category is
  not convex, the rules become inconsistent (the lines cross in awkward ways)
  and the classification produced is erroneous. To solve this problem, the
  complete set of rules is divided into several subsets (subcategories) that
  break the original non-convex shape into a set of convex polygons. Note
  that any polygon can be divided in a number of triangles, so this problem
  always has solution. However, in many cases (such as the one presented in
  the pixelclasser vignette) a complete triangulation is not needed.
  
  Internally, \code{classify_pixels()} classifies the points belonging to
  each subcategory and then joins the incidence matrices using the \code{or}
  operator, to create the matrix for the whole category.
}
\examples{
rule01 <- define_rule("R01", "g", "b",
                      list(c(0.35, 0.30), c(0.45, 0.10)), ">=")
rule02 <- define_rule("R02", "g", "b",
                      list(c(0.35, 0.253), c(0.45, 0.253)), ">=")

subcat01 <- define_subcat("Subcat_01", rule01, rule02)

}
\seealso{
\code{\link{define_rule}}, \code{\link{define_cat}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pixels.R
\name{plot_pixels}
\alias{plot_pixels}
\title{Plot the pixels of a transformed image}
\usage{
plot_pixels(image_rgb, x_axis, y_axis, ...)
}
\arguments{
\item{image_rgb}{an object produced by \code{read_image()}.}

\item{x_axis}{a character string indicating which colour variable use as x.}

\item{y_axis}{a character string indicating which colour variable use as y.}

\item{\dots}{additional graphical parameters, mainly to set the colour
(\code{col}) of the points.}
}
\value{
The function does not return any value.
}
\description{
This function is a wrapper for function \code{points()} in package
\code{graphics} for plotting the pixels of a transformed rgb image on the
triangular diagram previously created by \code{plot_rgb_plane()}.
}
\details{
It is advantageous to specify a colour such as \code{"#00000005"}
  which is black but almost transparent. In this way a kind of density plot
  is created because the clustering of points creates areas of darker colour.
  Note that a colour without specific transparency information defaults to an
  opaque colour, so \code{"#000000"} is the same as \code{"#000000ff"}. The
  colours can be specified in any form understandable by
  \code{grDevices::col2rgb}, but the hexadecimal string allows setting the
  colour transparency and it is the preferred style. Note also that the
  points are plotted using pch = ".", as any other symbol would clutter the
  graph.
  
  Warning: plotting several million points in an R graph is a slow process. 
  Be patient or reduce the size of the images as much as possible.
  Having a nice smartphone with a petapixel camera sensor is good for
  artistic purposes, but not always for efficient scientific work.
}
\examples{

# Plotting the pixels of the example image included in this package
ivy_oak_rgb <- read_image(system.file("extdata", "IvyOak400x300.JPG",
                                       package = "pixelclasser"))
plot_rgb_plane("g", "b")
plot_pixels(ivy_oak_rgb, "g", "b", col = "#00000005")

}
\seealso{
\code{\link{plot_rgb_plane}}, \code{\link[grDevices]{col2rgb}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_rule.R
\name{plot_rule}
\alias{plot_rule}
\title{Plots the line that defines a rule}
\usage{
plot_rule(rule, label = "", ...)
}
\arguments{
\item{rule}{an object of class \code{pixel_rule} produced by
\code{define_rule()}.}

\item{label}{a string to label the line. It is attached at the coordinates of
the second point used to define the line.}

\item{\dots}{additional graphical parameters passed to the underlying
\code{lines()} function, for example to define the line colour or dashing
style. They are also used for the line label.}
}
\value{
The function does not return any value.
}
\description{
This function draws the line that defines a rule on the plot created by
\code{plot_rgb_plane()}.
}
\details{
The function uses the information stored in the \code{pixel_rule
  object} to plot the line.

  Use the \dots to set the colour and other characteristics of the line. Use
  any character string understood by \code{col2rgb()}.
  
  Labels can be added to the rule using \code{label_rule()}.
}
\examples{
rule_01 <- define_rule("rule_01", "g", "b",
                      list(c(0.345, 1/3), c(0.40, 0.10)), "<")

plot_rgb_plane("g", "b")
plot_rule(rule_01, col = "green")

}
\seealso{
\code{\link{plot_rgb_plane}}, \code{\link{define_rule}},
  \code{\link{label_rule}} \code{\link[grDevices]{col2rgb}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_cat.R
\name{define_cat}
\alias{define_cat}
\title{Creates a category object}
\usage{
define_cat(cat_name, cat_colour, ...)
}
\arguments{
\item{cat_name}{a character string containing the name of the category.}

\item{cat_colour}{a character string defining the colour to paint the pixels
with when creating a classified picture.}

\item{\dots}{a list of \code{pixel_subcat} objects, or \code{pixel_rule}
objects in case that subcategories are not needed. A mixed list of
\code{pixel_rule} and \code{pixel_subcat} objects is not allowed.}
}
\value{
A list of class \code{pixel_cat} with the following elements:
  \itemize{ \item \code{name}: a character string containing the name of the
  pixel category. \item \code{colour}: a character string describing the
  colour of the pixels of the category in the classified images. \item
  \code{subcats}: a list containing the subcategories. Their names are the
  names of the elements of the list. }
}
\description{
Creates an object of class \code{pixel_cat}, which contains a list of objects
of class \code{pixel_subcat}.
}
\details{
The function receives a list of objects of class \code{pixel_subcat}
  and creates a list of class \code{pixel_cat} with them. However, for cases
  that does not need subcategories, i e that only need a set of rules,need a
  single set of rules, these can be passed to the function, which creates an
  internal subcategory object to contain them. See the examples below.

  Note that it is an error to pass a mixture of \code{pixel_rule} and
  \code{pixel_subcat} objects.

  \code{colour} can be specified in any form understood by
  \code{grDevices::col2grb}.
}
\examples{
# The rules are not consistent, they are only useful as examples
rule01 <- define_rule("R01", "g", "b",
                      list(c(0.35, 0.30), c(0.45, 0.10)), ">=")
rule02 <- define_rule("R02", "g", "b",
                      list(c(0.35, 0.253), c(0.45, 0.253)), ">=")
rule03 <- define_rule("R03", "g", "b",
                      list(c(0.35, 0.29), c(0.49, 0.178)), ">=")
rule04 <- define_rule("R04", "g", "b",
                      list(c(0.35, 0.253), c(0.45, 0.253)), "<")

subcat01 <- define_subcat("Subcat01", rule01, rule02)
subcat02 <- define_subcat("Subcat02", rule03, rule04)

cat01 <- define_cat("Cat01", "#ffae2a", subcat01, subcat02)

# A single category defined by a set of rules, not subcategories
cat02 <- define_cat("Cat02", "#00ae2a", rule01, rule02, rule03)

}
\seealso{
\code{\link{define_rule}}, \code{\link{define_subcat}},
  \code{\link[grDevices]{col2rgb}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pixelclasser-package.R
\docType{package}
\name{pixelclasser}
\alias{pixelclasser}
\title{pixelclasser: Functions to classify pixels by colour}
\description{
\code{pixelclasser} contains functions to classify the pixels of an image
file (in format jpeg or tiff) by its colour. It uses a simple form of the
technique known as Support Vector Machine, adapted to this particular
problem. The original colour variables (\code{R, G, B}) are transformed into
colour proportions (\code{r, g, b}), and the resulting two dimensional plane,
defined by any convenient pair of the transformed variables is divided in
several subsets (categories) by one or more straight lines (rules) selected
by the user. Finally, the pixels belonging to each category are identified
using the rules, and a classified image can be created and saved.
}
\details{
To classify the pixels of an image, a series of steps must be done
in the following order, using the functions shown in parenthesis:
\itemize{
\item import the image into an R array of transformed (\code{rgb}) data
(\code{read_image()}).
\item plot the pixels of the image on the plane of two transformed variables
that shows the categories of pixels most clearly (\code{plot_rgb_plane()},
\code{plot_pixels}).
\item trace lines between the pixel clusters and use them to create
classification rules (\code{place_rule()}, \code{define_rule},
\code{plot_rule()}).
\item combine the rules to define categories. Sometimes the rules are
combined into subcategories and these into categories (\code{define_cat()}, 
\code{define_subcat()}).
\item use the categories to classify the pixels (\code{classify_pixels()}).
\item save the results of the classification as an image, if needed
(\code{save_clasif_image()}).
}

These steps are explained in depth in the vignette included in the package.
}
\author{
Carlos Real (carlos.real@usc.es)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_colours.R
\name{transform_colours}
\alias{transform_colours}
\title{Transforms RGB values into proportions (rgb values)}
\usage{
transform_colours(image_array)
}
\arguments{
\item{image_array}{an array of class \code{image_array} created by function
\code{read_image()}.}
}
\value{
Returns an array of class \code{transformed_image} containing the
  proportions of each colour variable in the pixels of the image. The third
  dimension of the array is named "bands" and its elements are labelled as
  "r", "g" and "b", respectively.
}
\description{
This function transforms an array of RGB absolute values into a similar array
containing the proportion of each band (= colour variable): r g and b.
}
\details{
The proportions are calculated as \code{r} = \code{R} / (\code{R + G
  + B}), and so on. It is used by function read_image().
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define_rule.R
\name{define_rule}
\alias{define_rule}
\title{Creates a rule object}
\usage{
define_rule(rule_name, x_axis, y_axis, rule_points, comp_op)
}
\arguments{
\item{rule_name}{a character string containing the name of the rule.}

\item{x_axis}{a character string selecting the colour variable used as x
axis, one of \code{"r"}, \code{"g"} or \code{"b"}.}

\item{y_axis}{a character string selecting the colour variable used as y
axis, one of \code{"r"}, \code{"g"} or \code{"b"}.}

\item{rule_points}{either an object of  of class \code{"rule_points"} created
with function \code{place_rule()}, or a list containing the coordinates of
two points defining the line.}

\item{comp_op}{a character string containing one of the comparison operators
\code{">", ">=", "<", "<="}.}
}
\value{
A list of class \code{pixel_rule} containing the following elements:
  \itemize{
  \item \code{rule_name}: a character string containing the rule name.
  \item \code{rule_text}: a character string containing the mathematical
  expression of the rule.
  \item \code{comp_op}: a character string containing the comparison operator
  used in the rule.
  \item \code{a}: a numerical vector containing the parameter \code{a}
  (slope) of the line.
  \item \code{c}: a numerical vector containing the parameter \code{c}
  (intercept) of the line.
  \item \code{x_axis}: a character string containing the colour variable
  selected as \code{x} axis.
  \item \code{y_axis}: a character string containing the colour variable
  selected as \code{y} axis.
  \item \code{first_point}: a numerical vector containing the coordinates of
  the first point used to estimate the line equation.
  \item \code{second_point}: a numerical vector containing the coordinates of
  the second point.
  }
}
\description{
Creates an object of class \code{pixel_rule} from a line in \code{rgb} space,
defined by the user, and a relational operator.
}
\details{
This function estimates the slope (\code{a}) and intercept
  (\code{c}) of the line \code{y = ax + c} using the coordinates of two
  points on the line. \code{x} and \code{y} are two colour variables selected
  by the user (\code{r}, \code{g}, or \code{b}). The line divides the plane
  in two subsets and the comparison operator selects the subset that contains
  the points (pixels) of interest.

  When a list of two points is passed in \code{rule_points}, it is internally
  converted into an an object of class \code{rule_points}.

  The pair of points used to define the line are not constrained to belong to
  the area occupied by the pixels, but they are used by \code{plot_rule()} as
  the start and end of the plotted line. Therefore, the extremes of the line
  can be selected in the most convenient way, provided that the line divides
  correctly the categories. Convenience means that the line should seem nice
  in the plot, if this matters.
  
  Because the variables were transformed into proportions, the pixel are
  always inside the triangle defined by the points \code{(0, 0), (1, 0), (0,
  1)}. So, the sides of this triangle can be considered as implicit rules
  which do not need to be created. In this way, a single line creates two
  polygons by cutting the triangle in two. The implicit rules can reduce the
  number of rules to create in most cases.
}
\examples{
# Creating the line by passing the coordinates of two points on the line:
rule01 <- define_rule("rule01", "g", "b",
                      list(c(0.35, 0.30), c(0.45, 0.10)),">")

# A vertical line as a rule; note that the equation is simplified
rule02 <- define_rule("rule02", "g", "b",
                      list(c(0.35, 0.30), c(0.35, 0.00)), ">")
\dontrun{
# Creating the rule by passing an object of type rule_point:
rule_points01 <- place_rule("g", "b")
rule03 <- define_rule("rule03", "g", "b", rule_points01,">")

# Note that the creation of the intermediate object can be avoided:
rule04 <- define_rule("rule04", "g", "b", place_rule("g", "b"),">")
}
}
\seealso{
\code{\link{define_subcat}}, \code{\link{define_cat}},
  \code{\link{plot_rule}}, \code{\link{plot_rgb_plane}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_image.R
\name{read_image}
\alias{read_image}
\title{Imports a jpg or tiff file.}
\usage{
read_image(file_name)
}
\arguments{
\item{file_name}{A character string containing the name of the image file.}
}
\value{
Returns an array of dimensions \code{r x c x 3} and class
  \code{transformed_image}, being \code{r} and \code{c} the number of rows
  and columns in the image. The last dimension corresponds to the \code{R},
  \code{G} and \code{B} variables (= bands) that define the colours of the
  pixels. The values in the array are the proportions
  of each colour (\code{r, g, b}), i.e. \code{r} = \code{R} / (\code{R + G +
  B}), and so on.
}
\description{
Imports an image file (in JPEG or TIFF format) into an array, and converts
the original \code{R}, \code{G} and \code{B} values in the file into
proportions (\code{r}, \code{g} and \code{b} variables).
}
\details{
This function calls the functions \code{readJPEG()} and
  \code{readTIFF()} in packages \code{jpeg} and \code{tiff} to import the
  data into an R array. Then it transforms the data into proportions
}
\examples{

# An example that loads the example file included in the package
ivy_oak_rgb <- read_image(system.file("extdata", "IvyOak400x300.JPG", 
                                       package = "pixelclasser"))

}
\seealso{
For more information about jpeg and tiff file formats, see the help
  pages of \code{\link[jpeg]{readJPEG}} and \code{\link[tiff]{readTIFF}}
  functions in packages \code{jpeg} and \code{tiff}, respectively.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_rgb_plane.R
\name{plot_rgb_plane}
\alias{plot_rgb_plane}
\title{Plots a triangular plot to be filled with pixels and rules}
\usage{
plot_rgb_plane(
  x_axis,
  y_axis,
  plot_limits = TRUE,
  plot_guides = TRUE,
  plot_grid = TRUE,
  ...
)
}
\arguments{
\item{x_axis}{a character string indicating which colour variable use as x.}

\item{y_axis}{a character string indicating which colour variable use as y.}

\item{plot_limits}{a logical value. When TRUE (default) the limits of the
area where the pixels can be found are plotted.}

\item{plot_guides}{a logical value. When TRUE (default) the limits of the
area where one variable dominates are plotted.}

\item{plot_grid}{a logical value; if TRUE (default) a grid is added.}

\item{\dots}{allows passing graphical parameters to the plotting functions.}
}
\value{
The function does not return any value.
}
\description{
Plots a plane of the two variables selected by the user (\code{r}, \code{g}
or \code{b}) and, to serve as visual references, the lines limiting the
triangular area that can contain pixels (in blue) and the areas where one of
the colour variables has the larger proportion values (in red). Points
representing the pixels of a transformed image and lines representing the
rules can be later added to the plot using functions \code{plot_pixels()} and
\code{plot_rule()}.
}
\details{
Graphical parameters can be passed to the function to modify the
  appearance of the plot. Intended for passing \code{xlim} and \code{ylim}
  values to plot the part of the graph where the points are concentrated.
  
  Because the variables were transformed into proportions, the pixel are
  always inside the triangle defined by the points \code{(0, 0), (1, 0), (0,
  1)}. This triangle is plotted in blue. The point where all three variables
  have the same value is \code{(1/3, 1/3)}. The lines joining this point with
  the centers of the triangle sides divide the areas where one of the three
  variables has higher proportions than the other two. These lines are
  plotted as visual aids, so they can be deleted at will.
}
\examples{
# Simplest call
plot_rgb_plane("g", "b")

# Plane without the red lines
plot_rgb_plane("g", "b", plot_guides = FALSE)

# Restricting the plane area to show
plot_rgb_plane("g", "b", xlim = c(0.2, 0.5), ylim = c(0.0, 0.33))

}
\seealso{
\code{\link{plot_pixels}}, \code{\link{plot_rule}},
  \code{\link{define_rule}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/label_rule.R
\name{label_rule}
\alias{label_rule}
\title{Adds a label to the rule}
\usage{
label_rule(rule, label = "", shift = c(0, 0), ...)
}
\arguments{
\item{rule}{an object of class \code{pixel_rule}.}

\item{label}{a string to label the line. It is attached at the coordinates of
the start (first point) of the line.}

\item{shift}{a numeric vector to set the displacement of the label from the
start of the line. Expressed in graph units, defaults to c(0, 0).}

\item{\dots}{additional graphical parameters passed to the underlying
\code{graphics::text()} function.}
}
\value{
The function does not return any value.
}
\description{
This function adds a label to the line that represents a rule on a plot
created by \code{plot_rgb_plane()}.
}
\details{
The function uses the information stored in the pixel_rule object to
  plot the label at the start of the line. The \code{shift} values, expressed
  in plot coordinates, are added to the coordinates of that point to place
  the label elsewhere. Note that \dots can be used to pass values for the
  \code{adj} parameter to the underlying \code{graphics::text()} function,
  which also modifies the position of the label.
  
  Use a character string understood by \code{grDevices::col2rgb()} to set
  the colour of the label.
}
\examples{
rule_01 <- define_rule("rule_01", "g", "b",
                       list(c(0.1, 0.8), c(0.40, 0.10)), "<")
plot_rgb_plane("g", "b")

# The rule is represented as a green line
plot_rule(rule_01, col = "green")

# And the label is added in three different positions by passing col and adj
# to the underlying function
label_rule(rule_01, label = expression('R'[1]*''), shift = c(0,0),
           col = 'black', adj = 1.5)
label_rule(rule_01, label = expression('R'[1]*''), shift = c(0.2, -0.4),
           col = 'blue', adj = 0)
label_rule(rule_01, label = expression('R'[1]*''), shift = c(0.3, -0.7),
           col = 'black', adj = -0.5)

}
\seealso{
\code{\link{plot_rgb_plane}}, \code{\link{define_rule}},
  \code{\link[grDevices]{col2rgb}}, \code{\link[graphics]{text}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/place_rule.R
\name{place_rule}
\alias{place_rule}
\title{Places a line on the rgb plot}
\usage{
place_rule(x_axis, y_axis, line_type = "f")
}
\arguments{
\item{x_axis}{a character string indicating the colour variable that
corresponds to the x axis, one of \code{"r"}, \code{"g"} or \code{"b"}.}

\item{y_axis}{a character string indicating the colour variable that
corresponds to the y axis, one of \code{"r"}, \code{"g"} or \code{"b"}.}

\item{line_type}{a character string indicating that the line is vertical
\code{"v"}, horizontal \code{"h"} or free (\code{"f"}, the default).}
}
\value{
A list of class \code{rule_points} containing the following elements:
\itemize{
\item \code{x_axis}: a character string containing the colour variable
  selected as \code{x} axis.
\item \code{y_axis}: a character string containing the colour variable
  selected as \code{y} axis.
\item \code{first_point}: coordinates of the start point of the line.
\item \code{second_point}: coordinates of the end point of the line.
}
}
\description{
A wrapper function for \code{graphics::locator} that makes the creation
of rules easier.
}
\details{
This function calls \code{graphics::locator} allowing to select two
  points, plots the line joining these points and returns a list
  containing their coordinates. The coordinates are rearranged to
  pass them to \code{define_rule()}.

  True horizontal and vertical lines are difficult to create by hand. In
  these cases, specifying \code{"vertical"} or \code{"horizontal"} (partial
  match allowed, i e "h") will copy the appropriate coordinate value from the
  first point to the second. Note that this is done after \code{locator()}
  returns, so the plot will show the line joining the original points, not
  the corrected ones. Use \code{plot_rule()} to see corrected line.
}
\examples{
\dontrun{
plot_rgb_plane("r", "g")
line01 <- place_rule("r", "g")          # A "free" line
line02 <- place_rule("r", "g", "h")     # A horizontal line
}

}
\seealso{
\code{\link[graphics]{locator}}, \code{\link{define_rule}},
  \code{\link{plot_rule}}, \code{\link{plot_rgb_plane}}
}
