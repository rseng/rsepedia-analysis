# rOpenSci: The *magick* package <img src="hexlogo.png" align="right" height="134.5" />

> Advanced Image-Processing in R

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/magick?branch=master&svg=true)](https://ci.appveyor.com/project/jeroen/magick)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/magick)](https://cran.r-project.org/package=magick)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/magick)](https://cran.r-project.org/package=magick)
[![R-CMD-check](https://github.com/ropensci/magick/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/magick/actions)
<!-- badges: end -->

Bindings to ImageMagick: the most comprehensive open-source image
processing library available. Supports many common formats (png, jpeg, tiff,
pdf, etc) and manipulations (rotate, scale, crop, trim, flip, blur, etc).
All operations are vectorized via the Magick++ STL meaning they operate either
on a single frame or a series of frames for working with layers, collages,
or animation. In RStudio images are automatically previewed when printed to
the console, resulting in an interactive editing environment.

## Documentation

About the R package:

  - [Getting started: The magick package: Advanced Image-Processing in R](https://docs.ropensci.org/magick/articles/intro.html)
  - [rOpenSci Community Call (recording)](https://vimeo.com/channels/rocommunitycalls/180799058)

About the underlying library:

 - [Magick++ Tutorial](https://www.imagemagick.org/Magick++/tutorial/Magick++_tutorial.pdf)
 - [Magick++ STL Documentation](https://www.imagemagick.org/Magick++/STL.html)

## Hello World

**Run examples in RStudio** to see live previews of the images! If you do not use RStudio, use `image_browse` to open images. On Linux you can also use `image_display` to get an X11 preview.

```r
library(magick)
frink <- image_read("https://jeroen.github.io/images/frink.png")
image_trim(frink)
image_scale(frink, "200x200")
image_flip(frink)
image_rotate(frink, 45) ## <-- result of this is shown
image_negate(frink)
frink %>% 
  image_background("green") %>% 
  image_flatten() %>%
  image_border("red", "10x10")
```

```r
image_rotate(frink, 45) %>% image_write("man/figures/frink-rotated.png")
```

![](man/figures/frink-rotated.png)

Effects

```r
image_oilpaint(frink)
image_implode(frink)
image_charcoal(frink) ## <-- result of this is shown
image_blur(frink)
image_edge(frink)
```

```r
image_charcoal(frink) %>% image_write("man/figures/frink-charcoal.png")
```

![](man/figures/frink-charcoal.png)

Create GIF animation:

```r
# Download images
oldlogo <- image_read("https://developer.r-project.org/Logo/Rlogo-2.png")
newlogo <- image_read("https://jeroen.github.io/images/Rlogo-old.png")
logos <- c(oldlogo, newlogo)
logos <- image_scale(logos, "400x400")

# Create GIF
(animation1 <- image_animate(logos))
image_write(animation1, "man/figures/anim1.gif")

# Morph effect  <-- result of this is shown
(animation2 <- image_animate(image_morph(logos, frames = 10)))
image_write(animation2, "man/figures/anim2.gif")
```

![](man/figures/anim2.gif)

Read GIF animation frames. See the [rotating earth example GIF](https://upload.wikimedia.org/wikipedia/commons/2/2c/Rotating_earth_%28large%29.gif).

```r
earth <- image_read("https://upload.wikimedia.org/wikipedia/commons/2/2c/Rotating_earth_%28large%29.gif")
length(earth)
earth[1]
earth[1:3]
earth1 <- rev(image_flip(earth)) ## How Austrialans see earth
image_write(earth1, "man/figures/earth1.gif") ## <-- result of this is shown
```

![](man/figures/earth1.gif)

R logo with dancing banana

```r
logo <- image_read("https://www.r-project.org/logo/Rlogo.png")
banana <- image_read("https://jeroen.github.io/images/banana.gif")
front <- image_scale(banana, "300")
background <- image_scale(logo, "400")
frames <- lapply(as.list(front), function(x) image_flatten(c(background, x)))
image_write(image_animate(image_join(frames)), "man/figures/Rlogo-banana.gif")
```

![](man/figures/Rlogo-banana.gif)

## Use magick in Shiny Apps

This demo application shows how to use magick with shiny: https://github.com/jeroen/shinymagick

## Installation

Binary packages for __macOS__ or __Windows__ can be installed directly from CRAN:

```r
install.packages("magick")
```

Installation from source on Linux or OSX requires the imagemagick [`Magick++`](https://www.imagemagick.org/Magick++/Documentation.html) library. On __Debian or Ubuntu__ install [libmagick++-dev](https://packages.debian.org/testing/libmagick++-dev):

```
sudo apt-get install -y libmagick++-dev
```

If you are on __Ubuntu__ 14.04 (trusty) or 16.04 (xenial) you can get a more recent backport from the ppa:

```
sudo add-apt-repository -y ppa:cran/imagemagick
sudo apt-get update
sudo apt-get install -y libmagick++-dev 
```

On __Fedora__,  __CentOS or RHEL__ we need [ImageMagick-c++-devel](https://src.fedoraproject.org/rpms/ImageMagick). However on CentOS the system version of ImageMagick is quite old. More recent versions are available from the [ImageMagick downloads](https://www.imagemagick.org/download/linux/CentOS/x86_64/) website.

```
sudo yum install ImageMagick-c++-devel
````

On __macOS__ use [imagemagick@6](https://github.com/Homebrew/homebrew-core/blob/master/Formula/imagemagick@6.rb) from Homebrew.

```
brew install imagemagick@6
```

The unversioned homebrew formula`imagemagick` can also be used, however it has some unsolved OpenMP problems. 

There is also a fork of imagemagick called graphicsmagick, but this doesn't work for this package.
---
title: "The magick package: Advanced Image-Processing in R"
date: "`r Sys.Date()`"
output:
  html_document:
    fig_caption: false
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
    toc_depth: 3
vignette: >
  %\VignetteIndexEntry{The magick package: Advanced Image-Processing in R}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE}
dev.off <- function(){
  invisible(grDevices::dev.off())
}
```

The [magick](https://cran.r-project.org/package=magick) package provide a modern and simple toolkit for image processing in R. It wraps the [ImageMagick STL](https://www.imagemagick.org/Magick++/STL.html) which is the most comprehensive open-source image processing library available today.

The ImageMagick library has an overwhelming amount of functionality. Magick exposes a decent subset of it, but it is impossible to document everything in detail. This article introduces some basic concepts and examples to get started.

## Installing `magick`

On Windows or macOS the package is most easily installed via CRAN.

```r
install.packages("magick")
```

The binary CRAN packages work out of the box and have most important features enabled.
Use `magick_config` to see which features and formats are supported by your version of ImageMagick.


```{r}
library(magick)
str(magick::magick_config())
```


### Build from source

On Linux you need to install the ImageMagick++ library: on Debian/Ubuntu this is called [libmagick++-dev](https://packages.debian.org/testing/libmagick++-dev):

```
sudo apt-get install libmagick++-dev
```

On Fedora or CentOS/RHEL we need [ImageMagick-c++-devel](https://src.fedoraproject.org/rpms/ImageMagick):

```
sudo yum install ImageMagick-c++-devel
```

To install from source on macOS you need `imagemagick@6` from homebrew.

```
brew install imagemagick@6
```

Unfortunately the current `imagemagick@6` configuration on homebrew disables a bunch of features, including librsvg and fontconfig. Therefore the quality of fonts and svg rendering might be suboptimal. The is not a problem for the CRAN binary package.

## Image IO

What makes magick so magical is that it automatically converts and renders all common image formats. ImageMagick supports dozens of formats and automatically detects the type. Use `magick::magick_config()` to list the formats that your version of ImageMagick supports.

### Read and write

Images can be read directly from a file path, URL, or raw vector with image data with `image_read`. The `image_info` function shows some meta data about the image, similar to the imagemagick `identify` command line utility.

```{r, eval = require(rsvg, quietly = TRUE)}
library(magick)
tiger <- image_read_svg('http://jeroen.github.io/images/tiger.svg', width = 350)
print(tiger)
```

We use `image_write` to export an image in any format to a file on disk, or in memory if `path = NULL`.

```r
# Render svg to png bitmap
image_write(tiger, path = "tiger.png", format = "png")
```

If `path` is a filename, `image_write` returns `path` on success such that the result can be piped into function taking a file path.

### Converting formats

Magick keeps the image in memory in its original format. Specify the `format` parameter `image_write` to convert to another format. You can also internally convert the image to another format earlier, before applying transformations. This can be useful if your original format is lossy.

```{r, eval = require(rsvg, quietly = TRUE)}
tiger_png <- image_convert(tiger, "png")
image_info(tiger_png)
```

Note that size is currently 0 because ImageMagick is lazy (in the good sense) and does not render until it has to.

### Preview

IDE's with a built-in web browser (such as RStudio) automatically display magick images in the viewer. This results in a neat interactive image editing environment.

<img id="rstudioimg" alt="rstudio">
<script>
//this is a hack to prevent pandoc 'self_contained' from embedding this image
//in future version of pandoc we can use the image and set 'data-external=1' instead
window.onload = function(){
  document.getElementById("rstudioimg").src = "https://jeroen.github.io/images/magick-rstudio.png";
}
</script>

Alternatively, on Linux you can use `image_display` to preview the image in an X11 window. Finally `image_browse` opens the image in your system's default application for a given type.

```r
# X11 only
image_display(tiger)

# System dependent
image_browse(tiger)
```

Another method is converting the image to a raster object and plot it on R's graphics display. However this is very slow and only useful in combination with other plotting functionality. See [#raster](#raster-images) below.

## Transformations

The best way to get a sense of available transformations is walk through the examples in the `?transformations` help page in RStudio. Below a few examples to get a sense of what is possible.

### Cut and edit

Several of the transformation functions take an `geometry` parameter which requires a special syntax of the form `AxB+C+D` where each element is optional. Some examples:

  - `image_crop(image, "100x150+50")`: *crop out `width:100px` and `height:150px` starting `+50px` from the left*
  - `image_scale(image, "200")`: *resize proportionally to width: `200px`*
  - `image_scale(image, "x200")`: *resize proportionally to height: `200px`*
  - `image_fill(image, "blue", "+100+200")`: *flood fill with blue starting at the point at `x:100, y:200`*
  - `image_border(frink, "red", "20x10")`: *adds a border of 20px left+right and 10px top+bottom*

The full syntax is specified in the [Magick::Geometry](http://www.imagemagick.org/Magick++/Geometry.html) documentation.

```{r}
# Example image
frink <- image_read("https://jeroen.github.io/images/frink.png")
```

```{r}

print(frink)

# Add 20px left/right and 10px top/bottom
image_border(image_background(frink, "hotpink"), "#000080", "20x10")

# Trim margins
image_trim(frink)

# Passport pica
image_crop(frink, "100x150+50")

# Resize
image_scale(frink, "300") # width: 300px
image_scale(frink, "x300") # height: 300px

# Rotate or mirror
image_rotate(frink, 45)
image_flip(frink)
image_flop(frink)

# Brightness, Saturation, Hue
image_modulate(frink, brightness = 80, saturation = 120, hue = 90)

# Paint the shirt orange
image_fill(frink, "orange", point = "+100+200", fuzz = 20)
```

With `image_fill` we can flood fill starting at pixel `point`. The `fuzz` parameter allows for the fill to cross for adjacent pixels with similarish colors. Its value must be between 0 and 256^2 specifying the max geometric distance between colors to be considered equal. Here we give professor frink an orange shirt for the World Cup.

### Filters and effects

ImageMagick also has a bunch of standard effects that are worth checking out.

```{r}
# Add randomness
image_blur(frink, 10, 5)
image_noise(frink)

# Silly filters
image_charcoal(frink)
image_oilpaint(frink)
image_negate(frink)
```

### Kernel convolution

The `image_convolve()` function applies a [kernel](https://en.wikipedia.org/wiki/Kernel_(image_processing)) over the image. Kernel convolution means that each pixel value is recalculated using the weighted neighborhood sum defined in the kernel matrix. For example lets look at this simple kernel:


```{r}
kern <- matrix(0, ncol = 3, nrow = 3)
kern[1, 2] <- 0.25
kern[2, c(1, 3)] <- 0.25
kern[3, 2] <- 0.25
kern
```

This kernel changes each pixel to the mean of its horizontal and vertical neighboring pixels, which results in a slight blurring effect in the right-hand image below:

```{r}
img <- image_resize(logo, "300x300")
img_blurred <- image_convolve(img, kern)
image_append(c(img, img_blurred))
```

Or use any of the [standard kernels](https://legacy.imagemagick.org/Usage/convolve/)

```{r}
img %>% image_convolve('Sobel') %>% image_negate()
img %>% image_convolve('DoG:0,0,2') %>% image_negate()
```

### Text annotation

Finally it can be useful to print some text on top of images:

```{r}
# Add some text
image_annotate(frink, "I like R!", size = 70, gravity = "southwest", color = "green")

# Customize text
image_annotate(frink, "CONFIDENTIAL", size = 30, color = "red", boxcolor = "pink",
  degrees = 60, location = "+50+100")

# Fonts may require ImageMagick has fontconfig
image_annotate(frink, "The quick brown fox", font = 'Times', size = 30)
```

Fonts that are supported on most platforms include `"sans"`, `"mono"`, `"serif"`, `"Times"`, `"Helvetica"`, `"Trebuchet"`, `"Georgia"`, `"Palatino"`or `"Comic Sans"`.

### Combining with pipes

Each of the image transformation functions returns a **modified copy** of the original image. It does not affect the original image.

```{r}
frink <- image_read("https://jeroen.github.io/images/frink.png")
frink2 <- image_scale(frink, "100")
image_info(frink)
image_info(frink2)
```

Hence to combine transformations you need to chain them:

```{r}
test <- image_rotate(frink, 90)
test <- image_background(test, "blue", flatten = TRUE)
test <- image_border(test, "red", "10x10")
test <- image_annotate(test, "This is how we combine transformations", color = "white", size = 30)
print(test)
```

Using `magrittr` pipe syntax makes it a bit more readable

```{r}
image_read("https://jeroen.github.io/images/frink.png") %>%
  image_rotate(270) %>%
  image_background("blue", flatten = TRUE) %>%
  image_border("red", "10x10") %>%
  image_annotate("The same thing with pipes", color = "white", size = 30)
```


## Image Vectors

The examples above concern single images. However all functions in magick have been vectorized to support working with layers, compositions or animation.

The standard base methods `[` `[[`, `c()` and `length()` are used to manipulate vectors of images which can then be treated as layers or frames. 

```{r}
# Download earth gif and make it a bit smaller for vignette
earth <- image_read("https://jeroen.github.io/images/earth.gif") %>%
  image_scale("200x") %>%
  image_quantize(128)

length(earth)
earth
head(image_info(earth))

rev(earth) %>% 
  image_flip() %>% 
  image_annotate("meanwhile in Australia", size = 20, color = "white")
```

### Layers

We can stack layers on top of each other as we would in Photoshop:

```{r}
bigdata <- image_read('https://jeroen.github.io/images/bigdata.jpg')
frink <- image_read("https://jeroen.github.io/images/frink.png")
logo <- image_read("https://jeroen.github.io/images/Rlogo.png")
img <- c(bigdata, logo, frink)
img <- image_scale(img, "300x300")
image_info(img)
```

A mosaic prints images on top of one another, expanding the output canvas such that that everything fits:

```{r}
image_mosaic(img)
```

Flattening combines the layers into a single image which has the size of the first image:

```{r}
image_flatten(img)
```

Flattening and mosaic allow for specifying alternative [composite operators](https://www.imagemagick.org/Magick++/Enumerations.html#CompositeOperator):

```{r}
image_flatten(img, 'Add')
image_flatten(img, 'Modulate')
image_flatten(img, 'Minus')
```

### Combining

Appending means simply putting the frames next to each other:

```{r}
image_append(image_scale(img, "x200"))
```

Use `stack = TRUE` to position them on top of each other:

```{r}
image_append(image_scale(img, "100"), stack = TRUE)
```

Composing allows for combining two images on a specific position:

```{r}
bigdatafrink <- image_scale(image_rotate(image_background(frink, "none"), 300), "x200")
image_composite(image_scale(bigdata, "x400"), bigdatafrink, offset = "+180+100")
```

### Pages

When reading a PDF document, each page becomes an element of the vector. Note that PDF gets rendered while reading so you need to specify the density immediately.

```{r, eval = require(pdftools, quietly = TRUE)}
manual <- image_read_pdf('https://cloud.r-project.org/web/packages/magick/magick.pdf', density = 72)
image_info(manual)
manual[1]
```

### Animation

Instead of treating vector elements as layers, we can also make them frames in an animation!

```{r}
image_animate(image_scale(img, "200x200"), fps = 1, dispose = "previous")
```

Morphing creates a sequence of `n` images that gradually morph one image into another. It makes animations 

```{r}
newlogo <- image_scale(image_read("https://jeroen.github.io/images/Rlogo.png"))
oldlogo <- image_scale(image_read("https://jeroen.github.io/images/Rlogo-old.png"))
image_resize(c(oldlogo, newlogo), '200x150!') %>%
  image_background('white') %>%
  image_morph() %>%
  image_animate(optimize = TRUE)
```


If you read in an existing GIF or Video file, each frame becomes a layer:

```{r}
# Foreground image
banana <- image_read("https://jeroen.github.io/images/banana.gif")
banana <- image_scale(banana, "150")
image_info(banana)
```

Manipulate the individual frames and put them back into an animation:

```{r}
# Background image
background <- image_background(image_scale(logo, "200"), "white", flatten = TRUE)

# Combine and flatten frames
frames <- image_composite(background, banana, offset = "+70+30")

# Turn frames into animation
animation <- image_animate(frames, fps = 10, optimize = TRUE)
print(animation)
```

Animations can be saved as GIF of MPEG files:

```r
image_write(animation, "Rlogo-banana.gif")
```

## Drawing and Graphics

A relatively recent addition to the package is a native R graphics device which produces a magick image object. This can either be used like a regular device for making plots, or alternatively to open a device which draws onto an existing image using pixel coordinates.

### Graphics device

The `image_graph()` function opens a new graphics device similar to e.g. `png()` or `x11()`. It returns an image object to which the plot(s) will be written. Each "page" in the plotting device will become a frame in the image object.

```{r}
# Produce image using graphics device
fig <- image_graph(width = 400, height = 400, res = 96)
ggplot2::qplot(mpg, wt, data = mtcars, colour = cyl)
dev.off()
```

We can easily post-process the figure using regular image operations.

```{r}
# Combine
out <- image_composite(fig, frink, offset = "+70+30")
print(out)
```

### Drawing device

Another way to use the graphics device is to draw on top of an exiting image using pixel coordinates.

```{r}
# Or paint over an existing image
img <- image_draw(frink)
rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
abline(h = 300, col = 'blue', lwd = '10', lty = "dotted")
text(30, 250, "Hoiven-Glaven", family = "monospace", cex = 4, srt = 90)
palette(rainbow(11, end = 0.9))
symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
  bg = 1:11, inches = FALSE, add = TRUE)
dev.off()
```

```{r}
print(img)
```

By default `image_draw()` sets all margins to 0 and uses graphics coordinates to match image size in pixels (width x height) where (0,0) is the top left corner. Note that this means the y axis increases from top to bottom which is the opposite of typical graphics coordinates. You can override all this by passing custom `xlim`, `ylim` or `mar` values to `image_draw`.

### Animated Graphics

The graphics device supports multiple frames which makes it easy to create animated graphics. The code below shows how you would implement the example from the very cool [gganimate](https://gganimate.com/) package using the magick graphics device.


```{r}
library(gapminder)
library(ggplot2)
img <- image_graph(600, 340, res = 96)
datalist <- split(gapminder, gapminder$year)
out <- lapply(datalist, function(data){
  p <- ggplot(data, aes(gdpPercap, lifeExp, size = pop, color = continent)) +
    scale_size("population", limits = range(gapminder$pop)) + geom_point() + ylim(20, 90) + 
    scale_x_log10(limits = range(gapminder$gdpPercap)) + ggtitle(data$year) + theme_classic()
  print(p)
})
dev.off()
animation <- image_animate(img, fps = 2, optimize = TRUE)
print(animation)
```

To write it to a file you would simply do:

```r
image_write(animation, "gapminder.gif")
```

## Raster Images

Magick images can also be converted to raster objects for use with R's graphics device. Thereby we can combine it with other graphics tools. However do note that R's graphics device is very slow and has a very different coordinate system which reduces the quality of the image.

### Base R rasters

Base R has an `as.raster` format which converts the image to a vector of strings. The paper [Raster Images in R Graphics](https://journal.r-project.org/archive/2011-1/RJournal_2011-1_Murrell.pdf) by Paul Murrell gives a nice overview.

```{r}
plot(as.raster(frink))
```

```{r, fig.width=7, fig.height=5}
# Print over another graphic
plot(cars)
rasterImage(frink, 21, 0, 25, 80)
```

### The `grid` package

The `grid` package makes it easier to overlay a raster on the graphics device without having to adjust for the x/y coordinates of the plot.

```{r, fig.width=5, fig.height=3}
library(ggplot2)
library(grid)
qplot(speed, dist, data = cars, geom = c("point", "smooth"))
grid.raster(frink)
```

## OCR text extraction

A recent addition to the package is to extract text from images using OCR. This requires the tesseract package:

```{r eval=FALSE}
install.packages("tesseract")
```

```{r, eval = require(tesseract, quietly = TRUE)}
img <- image_read("http://jeroen.github.io/images/testocr.png")
print(img)

# Extract text
cat(image_ocr(img))
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{_index_}
\alias{_index_}
\alias{magick}
\alias{magick-package}
\alias{imagemagick}
\title{Magick Image Processing}
\description{
The \code{magick} package for graphics and image processing in R. Important resources:
\itemize{
\item \href{https://docs.ropensci.org/magick/articles/intro.html}{R introduction vignette}: getting started
\item \href{https://www.imagemagick.org/Magick++/Image++.html}{Magick++ API} and
\href{https://www.imagemagick.org/Magick++/STL.html}{Magick++ STL} detailed descriptions of methods and parameters
}
}
\details{
Documentation is split into the following pages:
\itemize{
\item \link{analysis} - metrics and calculations: \code{compare}, \code{fft}
\item \link{animation} - manipulate or combine multiple frames: \code{animate},
\code{morph}, \code{mosaic}, \code{montage}, \code{average}, \code{append}, \code{apply}
\item \link{attributes} - image properties: \code{comment}, \code{info}
\item \link{color} - contrast, brightness, colors: \code{modulate}, \code{quantize}, \code{map}, \code{transparent},
\code{background}, \code{colorize}, \code{contrast}, \code{normalize}, \code{enhance}, \code{equalize}, \code{median}
\item \link{composite} - advanced joining: \code{composite}, \code{border}, \code{frame}
\item \link{device} - creating graphics and drawing on images
\item \link{editing} - basic image IO: \code{read}, \code{write}, \code{convert}, \code{join}, \code{display}, \code{brose}
\item \link{effects} - fun effects: \code{despecle}, \code{reducenoise}, \code{noise}, \code{blur}, \code{charcoal},
\code{edge}, \code{oilpaint}, \code{emboss}, \code{implode}, \code{negate}
\item \link{geometry} - specify points, areas and sizes using geometry syntax
\item \link{ocr} - extract text from image using \link[tesseract:tesseract]{tesseract} package
\item \link{options} - list option types and values supported in your version of ImageMagick
\item \link{painting} - flood fill and annotating text
\item \link{transform} - shape operations: \code{trim}, \code{chop}, \code{rotate}, \code{resize}, \code{scale}, \code{sample}
\code{crop}, \code{flip}, \code{flop}, \code{deskew}, \code{page}
}
}
\seealso{
Other image: 
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometry.R
\name{geometry}
\alias{geometry}
\alias{geometry_point}
\alias{geometry_area}
\alias{geometry_size_pixels}
\alias{geometry_size_percent}
\title{Geometry Helpers}
\usage{
geometry_point(x, y)

geometry_area(width = NULL, height = NULL, x_off = 0, y_off = 0)

geometry_size_pixels(width = NULL, height = NULL, preserve_aspect = TRUE)

geometry_size_percent(width = 100, height = NULL)
}
\arguments{
\item{x}{left offset in pixels}

\item{y}{top offset in pixels}

\item{width}{in pixels}

\item{height}{in pixels}

\item{x_off}{offset in pixels on x axis}

\item{y_off}{offset in pixels on y axis}

\item{preserve_aspect}{if FALSE, resize to width and height exactly, loosing original
aspect ratio. Only one of \code{percent} and \code{preserve_aspect} may be \code{TRUE}.}
}
\description{
ImageMagick uses a handy geometry syntax to specify coordinates and shapes
for use in image transformations. You can either specify these manually as
strings or use the helper functions below.
}
\details{
See \href{http://www.imagemagick.org/Magick++/Geometry.html}{ImageMagick Manual}
for details about the syntax specification.
Examples of \code{geometry} strings:
\itemize{
\item \strong{\code{"500x300"}}       -- \emph{Resize image keeping aspect ratio, such that width does not exceed 500 and the height does not exceed 300.}
\item \strong{\code{"500x300!"}}      -- \emph{Resize image to 500 by 300, ignoring aspect ratio}
\item \strong{\code{"500x"}}          -- \emph{Resize width to 500 keep aspect ratio}
\item \strong{\code{"x300"}}          -- \emph{Resize height to 300 keep aspect ratio}
\item \strong{\code{"50\%x20\%"}} -- \emph{Resize width to 50 percent and height to 20 percent of original}
\item \strong{\code{"500x300+10+20"}} -- \emph{Crop image to 500 by 300 at position 10,20}
}
}
\examples{
# Specify a point
logo <- image_read("logo:")
image_annotate(logo, "Some text", location = geometry_point(100, 200), size = 24)

# Specify image area
image_crop(logo, geometry_area(300, 300), repage = FALSE)
image_crop(logo, geometry_area(300, 300, 100, 100), repage = FALSE)

# Specify image size
image_resize(logo, geometry_size_pixels(300))
image_resize(logo, geometry_size_pixels(height = 300))
image_resize(logo, geometry_size_pixels(300, 300, preserve_aspect = FALSE))

# resize relative to current size
image_resize(logo, geometry_size_percent(50))
image_resize(logo, geometry_size_percent(50, 20))

}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ocr.R
\name{ocr}
\alias{ocr}
\alias{image_ocr}
\alias{image_ocr_data}
\title{Image Text OCR}
\usage{
image_ocr(image, language = "eng", HOCR = FALSE, ...)

image_ocr_data(image, language = "eng", ...)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{language}{passed to \link[tesseract:tesseract]{tesseract}. To install additional languages see
instructions in \link[tesseract:tessdata]{tesseract_download()}.}

\item{HOCR}{if \code{TRUE} return results as HOCR xml instead of plain text}

\item{...}{additional parameters passed to \link[tesseract:tesseract]{tesseract}}
}
\description{
Extract text from an image using the \link[tesseract:tesseract]{tesseract} package.
}
\details{
To use this function you need to tesseract first:\preformatted{  install.packages("tesseract")
}

Best results are obtained if you set the correct language in \link[tesseract:tesseract]{tesseract}.
To install additional languages see instructions in \link[tesseract:tessdata]{tesseract_download()}.
}
\examples{
if(require("tesseract")){
img <- image_read("http://jeroen.github.io/images/testocr.png")
image_ocr(img)
image_ocr_data(img)
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paint.R
\name{painting}
\alias{painting}
\alias{image_fill}
\alias{image_annotate}
\title{Image Painting}
\usage{
image_fill(image, color, point = "+1+1", fuzz = 0, refcolor = NULL)

image_annotate(
  image,
  text,
  gravity = "northwest",
  location = "+0+0",
  degrees = 0,
  size = 10,
  font = "",
  style = "normal",
  weight = 400,
  kerning = 0,
  decoration = NULL,
  color = NULL,
  strokecolor = NULL,
  boxcolor = NULL
)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{color}{a valid \href{https://www.imagemagick.org/Magick++/Color.html}{color string} such as
\code{"navyblue"} or \code{"#000080"}. Use \code{"none"} for transparency.}

\item{point}{a geometry_point string indicating the starting point of the flood-fill}

\item{fuzz}{relative color distance (value between 0 and 100) to be considered similar
in the filling algorithm}

\item{refcolor}{if set, \code{fuzz} color distance will be measured against this color,
not the color of the starting \code{point}. Any color (within \code{fuzz} color distance of
the given \code{refcolor}), connected to starting point will be replaced with the \code{color}.
If the pixel at the starting point does not itself match the given \code{refcolor}
(according to \code{fuzz}) then no action will be taken.}

\item{text}{character vector of length equal to 'image' or length 1}

\item{gravity}{string with
\href{https://www.imagemagick.org/Magick++/Enumerations.html#GravityType}{gravity}
value from \link{gravity_types}.}

\item{location}{geometry string with location relative to \code{gravity}}

\item{degrees}{rotates text around center point}

\item{size}{font-size in pixels}

\item{font}{string with font family such as \code{"sans"}, \code{"mono"}, \code{"serif"},
\code{"Times"}, \code{"Helvetica"}, \code{"Trebuchet"}, \code{"Georgia"}, \code{"Palatino"} or \code{"Comic Sans"}.}

\item{style}{value of \link{style_types} for example \code{"italic"}}

\item{weight}{thickness of the font, 400 is normal and 700 is bold.}

\item{kerning}{increases or decreases whitespace between letters}

\item{decoration}{value of \link{decoration_types} for example \code{"underline"}}

\item{strokecolor}{a \href{https://www.imagemagick.org/Magick++/Color.html}{color string}
adds a stroke (border around the text)}

\item{boxcolor}{a \href{https://www.imagemagick.org/Magick++/Color.html}{color string}
for background color that annotation text is rendered on.}
}
\description{
The \code{\link[=image_fill]{image_fill()}} function performs flood-fill by painting starting point and all
neighboring pixels of approximately the same color. Annotate prints some text on
the image.
}
\details{
Note that more sophisticated drawing mechanisms are available via the graphics
device using \link{image_draw}.

Setting a font, weight, style only works if your imagemagick is compiled
with fontconfig support.
}
\examples{
logo <- image_read("logo:")
logo <- image_background(logo, 'white')
image_fill(logo, "pink", point = "+450+400")
image_fill(logo, "pink", point = "+450+400", fuzz = 25)
# Add some text to an image
image_annotate(logo, "This is a test")
image_annotate(logo, "CONFIDENTIAL", size = 50, color = "red", boxcolor = "pink",
 degrees = 30, location = "+100+100")

# Setting fonts requires fontconfig support (and that you have the font)
image_annotate(logo, "The quick brown fox", font = "monospace", size = 50)
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edges.R
\name{edges}
\alias{edges}
\alias{image_edge}
\alias{image_canny}
\alias{image_hough_draw}
\alias{image_hough_txt}
\title{Edge / Line Detection}
\usage{
image_edge(image, radius = 1)

image_canny(image, geometry = "0x1+10\%+30\%")

image_hough_draw(
  image,
  geometry = NULL,
  color = "red",
  bg = "transparent",
  size = 3,
  overlay = FALSE
)

image_hough_txt(image, geometry = NULL, format = c("mvg", "svg"))
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{radius}{edge size in pixels}

\item{geometry}{geometry string, see details.}

\item{color}{a valid \href{https://www.imagemagick.org/Magick++/Color.html}{color string} such as
\code{"navyblue"} or \code{"#000080"}. Use \code{"none"} for transparency.}

\item{bg}{background color}

\item{size}{size in points to draw the line}

\item{overlay}{composite the drawing atop the input image. Only for \code{bg = 'transparent'}.}

\item{format}{output format of the text, either \code{svg} or \code{mvg}}
}
\description{
Best results are obtained by finding edges with \code{\link[=image_canny]{image_canny()}} and
then performing Hough-line detection on the edge image.
}
\details{
For Hough-line detection, the geometry format is \verb{\{W\}x\{H\}+\{threshold\}}
defining the size and threshold of the filter used to find 'peaks' in
the intermediate search image. For canny edge detection the format is
\verb{\{radius\}x\{sigma\}+\{lower\%\}+\{upper\%\}}. More details and examples are
available at the \href{https://legacy.imagemagick.org/Usage/transform/}{imagemagick website}.
}
\examples{
if(magick_config()$version > "6.8.9"){
shape <- demo_image("shape_rectangle.gif")
rectangle <- image_canny(shape)
rectangle \%>\% image_hough_draw('5x5+20')
rectangle \%>\% image_hough_txt(format = 'svg') \%>\% cat()
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effects.R
\name{effects}
\alias{effects}
\alias{image_despeckle}
\alias{image_reducenoise}
\alias{image_noise}
\alias{image_blur}
\alias{image_motion_blur}
\alias{image_charcoal}
\alias{image_oilpaint}
\alias{image_emboss}
\alias{image_implode}
\alias{image_negate}
\title{Image Effects}
\usage{
image_despeckle(image, times = 1L)

image_reducenoise(image, radius = 1L)

image_noise(image, noisetype = "gaussian")

image_blur(image, radius = 1, sigma = 0.5)

image_motion_blur(image, radius = 1, sigma = 0.5, angle = 0)

image_charcoal(image, radius = 1, sigma = 0.5)

image_oilpaint(image, radius = 1)

image_emboss(image, radius = 1, sigma = 0.5)

image_implode(image, factor = 0.5)

image_negate(image)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{times}{number of times to repeat the despeckle operation}

\item{radius}{radius, in pixels, for various transformations}

\item{noisetype}{string with a
\href{https://www.imagemagick.org/Magick++/Enumerations.html#NoiseType}{noisetype} value
from \link{noise_types}.}

\item{sigma}{the standard deviation of the Laplacian, in pixels.}

\item{angle}{angle, in degrees, for various transformations}

\item{factor}{image implode factor (special effect)}
}
\description{
High level effects applied to an entire image.
These are mostly just for fun.
}
\examples{
logo <- image_read("logo:")
image_despeckle(logo)
image_reducenoise(logo)
image_noise(logo)
image_blur(logo, 10, 10)
image_motion_blur(logo, 10, 10, 45)
image_charcoal(logo)
image_oilpaint(logo, radius = 3)
image_emboss(logo)
image_implode(logo)
image_negate(logo)
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coder.R
\name{coder_info}
\alias{coder_info}
\alias{magick_config}
\alias{magick_set_seed}
\title{Magick Configuration}
\usage{
coder_info(format)

magick_config()

magick_set_seed(seed)
}
\arguments{
\item{format}{image format such as \code{png}, \code{tiff} or \code{pdf}.}

\item{seed}{integer with seed value to use}
}
\description{
ImageMagick can be configured to support various additional tool and formats
via external libraries. These functions show which features ImageMagick supports
on your system.
}
\details{
Note that \code{coder_info} raises an error for unsupported formats.
}
\examples{
coder_info("png")
coder_info("jpg")
coder_info("pdf")
coder_info("tiff")
coder_info("gif")
# Reproduce random image
magick_set_seed(123)
image_blank(200,200, pseudo_image = "plasma:fractal")
}
\references{
\url{https://www.imagemagick.org/Magick++/CoderInfo.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot2.R
\name{image_ggplot}
\alias{image_ggplot}
\title{Image to ggplot}
\usage{
image_ggplot(image, interpolate = FALSE)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{interpolate}{passed to \link[ggplot2:annotation_raster]{ggplot2::annotation_raster}}
}
\description{
Create a ggplot with axes set to pixel coordinates and plot the raster image
on it using \link[ggplot2:annotation_raster]{ggplot2::annotation_raster}. See examples for how to plot an image
onto an existing ggplot.
}
\examples{
# Plot with base R
plot(logo)

# Plot image with ggplot2
library(ggplot2)
myplot <- image_ggplot(logo)
myplot + ggtitle("Test plot")

# Show that coordinates are reversed:
myplot + theme_classic()

# Or add to plot as annotation
image <- image_fill(logo, 'none')
raster <- as.raster(image)
myplot <- qplot(mpg, wt, data = mtcars)
myplot + annotation_raster(raster, 25, 35, 3, 5)

# Or overplot image using grid
library(grid)
qplot(speed, dist, data = cars, geom = c("point", "smooth"))
grid.raster(image)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defines.R
\name{defines}
\alias{defines}
\alias{image_set_defines}
\title{Set encoder defines}
\usage{
image_set_defines(image, defines)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{defines}{a named character vector with extra options to control reading.
These are the \verb{-define key\{=value\}} settings in the \href{http://www.imagemagick.org/script/command-line-options.php#define}{command line tool}.
Use an empty string for value-less defines, and NA to unset a define.}
}
\description{
So called 'defines' are properties that are passed along to external
filters and libraries. Usually defines are used in \link{image_read} or
\link{image_write} to control the image encoder/decoder, but you can also
set these manually on the image object.
}
\details{
The defines values must be a character string, where the names contain
the defines keys. Each name must be of the format "enc:key" where the
first part is the encoder or filter to which the key is passed. For
example \code{"png:...."} defines can control the encoding and decoding of
png images.

The \link{image_set_defines} function does not make a copy of the image, so
the defined values remain in the image object until they are overwritten
or unset.
}
\examples{
# Write an image
x <- image_read("https://jeroen.github.io/images/frink.png")
image_write(x, "frink.png")

# Pass some properties to PNG encoder
defines <- c("png:compression-filter" = "1", "png:compression-level" = "0")
image_set_defines(x, defines)
image_write(x, "frink-uncompressed.png")

# Unset properties
defines[1:2] = NA
image_set_defines(x, defines)
image_write(x, "frink-final.png")

# Compare size and cleanup
file.info(c("frink.png", "frink-uncompressed.png", "frink-final.png"))
unlink(c("frink.png", "frink-uncompressed.png", "frink-final.png"))
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attributes.R
\name{attributes}
\alias{attributes}
\alias{image_comment}
\alias{image_info}
\title{Image Attributes}
\usage{
image_comment(image, comment = NULL)

image_info(image)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{comment}{string to set an image comment}
}
\description{
Attributes are properties of the image that might be present on some images
and might affect image manipulation methods.
}
\details{
Each attribute can be get and set with the same function. The \code{\link[=image_info]{image_info()}}
function returns a data frame with some commonly used attributes.
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation.R
\name{segmentation}
\alias{segmentation}
\alias{image_connect}
\alias{image_split}
\alias{image_fuzzycmeans}
\title{Image Segmentation}
\usage{
image_connect(image, connectivity = 4)

image_split(image, keep_color = TRUE)

image_fuzzycmeans(image, min_pixels = 1, smoothing = 1.5)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{connectivity}{number neighbor colors which are considered part of a unique object}

\item{keep_color}{if TRUE the output images retain the color of the input pixel.
If FALSE all matching pixels are set black to retain only the image mask.}

\item{min_pixels}{the minimum number of pixels contained in a hexahedra before it can be considered valid (expressed as a percentage)}

\item{smoothing}{the smoothing threshold which eliminates noise in the second derivative of the histogram (higher values gives smoother second derivative)}
}
\description{
Basic image segmentation like connected components labelling, blob extraction and fuzzy c-means
}
\details{
\itemize{
\item \link{image_connect} Connect adjacent pixels with the same pixel intensities to do blob extraction
\item \link{image_split} Splits the image according to pixel intensities
\item \link{image_fuzzycmeans} Fuzzy c-means segmentation of the histogram of color components
}

\link{image_connect} performs blob extraction by scanning the image, pixel-by-pixel from top-left
to bottom-right where regions of adjacent pixels which share the same set of intensity values
get combined.
}
\examples{
# Split an image by color
img <- image_quantize(logo, 4)
layers <- image_split(img)
layers

# This returns the original image
image_flatten(layers)

# From the IM website
objects <- image_convert(demo_image("objects.gif"), colorspace = "Gray")
objects

# Split image in blobs of connected pixel levels
if(magick_config()$version > "6.9.0"){
objects \%>\%
  image_connect(connectivity = 4) \%>\%
  image_split()

# Fuzzy c-means
image_fuzzycmeans(logo)

logo \%>\%
  image_convert(colorspace = "HCL") \%>\%
  image_fuzzycmeans(smoothing = 5)
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{analysis}
\alias{analysis}
\alias{image_compare}
\alias{image_compare_dist}
\alias{image_fft}
\title{Image Analysis}
\usage{
image_compare(image, reference_image, metric = "", fuzz = 0)

image_compare_dist(image, reference_image, metric = "", fuzz = 0)

image_fft(image)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{reference_image}{another image to compare to}

\item{metric}{string with a \href{http://www.imagemagick.org/script/command-line-options.php#metric}{metric}
from \link[=metric_types]{metric_types()} such as \code{"AE"} or \code{"phash"}}

\item{fuzz}{relative color distance (value between 0 and 100) to be considered similar
in the filling algorithm}
}
\description{
Functions for image calculations and analysis. This part of the package needs more work.
}
\details{
For details see \href{https://www.imagemagick.org/Magick++/Image++.html}{Image++}
documentation. Short descriptions:
\itemize{
\item \link{image_compare} calculates a metric by comparing image with a reference image.
\item \link{image_fft} returns Discrete Fourier Transform (DFT) of the image as a
magnitude / phase image pair. I wish I knew what this means.
}

Here \code{image_compare()} is vectorized over the first argument and returns the diff image
with the calculated distortion value as an attribute.
}
\examples{
out1 <- image_blur(logo, 3)
out2 <- image_oilpaint(logo, 3)
input <- c(logo, out1, out2, logo)
if(magick_config()$version >= "6.8.7"){
  diff_img <- image_compare(input, logo, metric = "AE")
  attributes(diff_img)
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/composite.R
\name{composite}
\alias{composite}
\alias{image_composite}
\alias{image_border}
\alias{image_frame}
\alias{image_shadow_mask}
\alias{image_shadow}
\alias{image_shade}
\title{Image Composite}
\usage{
image_composite(
  image,
  composite_image,
  operator = "atop",
  offset = "+0+0",
  gravity = "northwest",
  compose_args = ""
)

image_border(image, color = "lightgray", geometry = "10x10", operator = "copy")

image_frame(image, color = "lightgray", geometry = "25x25+6+6")

image_shadow_mask(image, geometry = "50x10+30+30")

image_shadow(
  image,
  color = "black",
  bg = "white",
  geometry = "50x10+30+30",
  operator = "atop",
  offset = "+20+20"
)

image_shade(image, azimuth = 30, elevation = 30, color = FALSE)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{composite_image}{composition image}

\item{operator}{string with a
\href{https://www.imagemagick.org/Magick++/Enumerations.html#CompositeOperator}{composite operator}
from \link[=compose_types]{compose_types()}}

\item{offset}{string with either a \link[=gravity_types]{gravity_type} or a \link{geometry_point}
to set position of top image.}

\item{gravity}{string with
\href{https://www.imagemagick.org/Magick++/Enumerations.html#GravityType}{gravity}
value from \link{gravity_types}.}

\item{compose_args}{additional arguments needed for some composite operations}

\item{color}{Set to true to shade the red, green, and blue components of the image.}

\item{geometry}{a \href{https://www.imagemagick.org/Magick++/Geometry.html}{geometry string}
to set height and width of the border, e.g. \code{"10x8"}. In addition \link{image_frame} allows
for adding shadow by setting an offset e.g. \code{"20x10+7+2"}.}

\item{bg}{background color}

\item{azimuth}{position of light source}

\item{elevation}{position of light source}
}
\description{
Similar to the ImageMagick \code{composite} utility: compose an image on top of another one using a
\href{https://www.imagemagick.org/Magick++/Enumerations.html#CompositeOperator}{CompositeOperator}.
}
\details{
The \code{image_composite} function is vectorized over both image arguments: if the first image has
\code{n} frames and the second \code{m} frames, the output image will contain \code{n} * \code{m} frames.

The \link{image_border} function creates a slightly larger solid color frame and then composes the
original frame on top. The \link{image_frame} function is similar but has an additional feature to
create a shadow effect on the border (which is really ugly).
}
\examples{
# Compose images using one of many operators
imlogo <- image_scale(image_read("logo:"), "x275")
rlogo <- image_read("https://jeroen.github.io/images/Rlogo-old.png")

# Standard is atop
image_composite(imlogo, rlogo)

# Same as 'blend 50' in the command line
image_composite(imlogo, rlogo, operator = "blend", compose_args="50")

# Offset can be geometry or gravity
image_composite(logo, rose, offset = "+100+100")
image_composite(logo, rose, gravity = "East")

# Add a border frame around the image
image_border(imlogo, "red", "10x10")
image_frame(imlogo)
image_shadow(imlogo)
image_shade(imlogo)
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform.R
\name{transform}
\alias{transform}
\alias{image_trim}
\alias{image_chop}
\alias{image_rotate}
\alias{image_resize}
\alias{image_scale}
\alias{image_sample}
\alias{image_crop}
\alias{image_extent}
\alias{image_flip}
\alias{image_flop}
\alias{image_deskew}
\alias{image_deskew_angle}
\alias{image_page}
\alias{image_repage}
\alias{image_orient}
\alias{image_shear}
\alias{image_distort}
\title{Image Transform}
\usage{
image_trim(image, fuzz = 0)

image_chop(image, geometry)

image_rotate(image, degrees)

image_resize(image, geometry = NULL, filter = NULL)

image_scale(image, geometry = NULL)

image_sample(image, geometry = NULL)

image_crop(image, geometry = NULL, gravity = NULL, repage = TRUE)

image_extent(image, geometry, gravity = "center", color = "none")

image_flip(image)

image_flop(image)

image_deskew(image, threshold = 40)

image_deskew_angle(image, threshold = 40)

image_page(image, pagesize = NULL, density = NULL)

image_repage(image)

image_orient(image, orientation = NULL)

image_shear(image, geometry = "10x10", color = "none")

image_distort(image, distortion = "perspective", coordinates, bestfit = FALSE)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{fuzz}{relative color distance (value between 0 and 100) to be considered similar
in the filling algorithm}

\item{geometry}{a \link{geometry} string specifying area (for cropping) or size (for resizing).}

\item{degrees}{value between 0 and 360 for how many degrees to rotate}

\item{filter}{string with \href{https://www.imagemagick.org/Magick++/Enumerations.html#FilterTypes}{filter}
type from: \link{filter_types}}

\item{gravity}{string with
\href{https://www.imagemagick.org/Magick++/Enumerations.html#GravityType}{gravity}
value from \link{gravity_types}.}

\item{repage}{resize the canvas to the cropped area}

\item{color}{a valid \href{https://www.imagemagick.org/Magick++/Color.html}{color string} such as
\code{"navyblue"} or \code{"#000080"}. Use \code{"none"} for transparency.}

\item{threshold}{straightens an image. A threshold of 40 works for most images.}

\item{pagesize}{geometry string with preferred size and location of an image canvas}

\item{density}{geometry string with vertical and horizontal resolution in pixels of
the image. Specifies an image density when decoding a Postscript or PDF.}

\item{orientation}{string to set image orientation one of the \link{orientation_types}.
If \code{NULL} it applies auto-orientation which tries to infer the correct orientation
from the Exif data.}

\item{distortion}{string to set image orientation one of the \link{distort_types}.}

\item{coordinates}{numeric vector (typically of length 12) with distortion coordinates}

\item{bestfit}{if set to \code{TRUE} the size of the output image can be different from input}
}
\description{
Basic transformations like rotate, resize, crop and flip. The \link{geometry} syntax
is used to specify sizes and areas.
}
\details{
For details see \href{https://www.imagemagick.org/Magick++/STL.html}{Magick++ STL}
documentation. Short descriptions:
\itemize{
\item \link{image_trim} removes edges that are the background color from the image.
\item \link{image_chop} removes vertical or horizontal subregion of image.
\item \link{image_crop} cuts out a subregion of original image
\item \link{image_rotate} rotates and increases size of canvas to fit rotated image.
\item \link{image_deskew} auto rotate to correct skewed images
\item \link{image_resize} resizes using custom \href{https://www.imagemagick.org/Magick++/Enumerations.html#FilterTypes}{filterType}
\item \link{image_scale} and \link{image_sample} resize using simple ratio and pixel sampling algorithm.
\item \link{image_flip} and \link{image_flop} invert image vertically and horizontally
}

The most powerful resize function is \link{image_resize} which allows for setting
a custom resize filter. Output of \link{image_scale} is similar to \code{image_resize(img, filter = "point")}.

For resize operations it holds that if no \code{geometry} is specified, all frames
are rescaled to match the top frame.
}
\examples{
logo <- image_read("logo:")
logo <- image_scale(logo, "400")
image_trim(logo)
image_chop(logo, "100x20")
image_rotate(logo, 45)
# Small image
rose <- image_convert(image_read("rose:"), "png")

# Resize to 400 width or height:
image_resize(rose, "400x")
image_resize(rose, "x400")

# Resize keeping ratio
image_resize(rose, "400x400")

# Resize, force size losing ratio
image_resize(rose, "400x400!")

# Different filters
image_resize(rose, "400x", filter = "Triangle")
image_resize(rose, "400x", filter = "Point")
# simple pixel resize
image_scale(rose, "400x")
image_sample(rose, "400x")
image_crop(logo, "400x400+200+200")
image_extent(rose, '200x200', color = 'pink')
image_flip(logo)
image_flop(logo)
skewed <- image_rotate(logo, 5)
deskewed <- image_deskew(skewed)
attr(deskewed, 'angle')
if(magick_config()$version > "6.8.6")
  image_orient(logo)
image_shear(logo, "10x10")
building <- demo_image('building.jpg')
image_distort(building, 'perspective', c(7,40,4,30,4,124,4,123,85,122,100,123,85,2,100,30))
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoview.R
\name{autoviewer}
\alias{autoviewer}
\alias{autoviewer_enable}
\alias{autoviewer_disable}
\title{RStudio Graphics AutoViewer}
\usage{
autoviewer_enable()

autoviewer_disable()
}
\description{
This enables a \link{addTaskCallback} that automatically updates the viewer after
the state of a magick graphics device has changed. This is enabled by default in
RStudio.
}
\examples{
# Only has effect in RStudio (or other GUI with a viewer):
autoviewer_enable()

img <- magick::image_graph()
plot(1)
abline(0, 1, col = "blue", lwd = 2, lty = "solid")
abline(0.1, 1, col = "red", lwd = 3, lty = "dotted")

autoviewer_disable()
abline(0.2, 1, col = "green", lwd = 4, lty = "twodash")
abline(0.3, 1, col = "black", lwd = 5, lty = "dotdash")

autoviewer_enable()
abline(0.4, 1, col = "purple", lwd = 6, lty = "dashed")
abline(0.5, 1, col = "yellow", lwd = 7, lty = "longdash")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/threshold.R
\name{thresholding}
\alias{thresholding}
\alias{image_threshold}
\alias{image_level}
\alias{image_lat}
\title{Image thresholding}
\usage{
image_threshold(
  image,
  type = c("black", "white"),
  threshold = "50\%",
  channel = NULL
)

image_level(
  image,
  black_point = 0,
  white_point = 100,
  mid_point = 1,
  channel = NULL
)

image_lat(image, geometry = "10x10+5\%")
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{type}{type of thresholding, either one of lat, black or white (see details below)}

\item{threshold}{pixel intensity threshold percentage for black or white thresholding}

\item{channel}{a value of \code{\link[=channel_types]{channel_types()}} specifying which channel(s) to set}

\item{black_point}{value between 0 and 100, the darkest color in the image}

\item{white_point}{value between 0 and 100, the lightest color in the image}

\item{mid_point}{value between 0 and 10 used for gamma correction}

\item{geometry}{pixel window plus offset for LAT algorithm}
}
\description{
Thresholding an image can be used for simple and straightforward image segmentation.
The function \code{\link[=image_threshold]{image_threshold()}} allows to do black and white thresholding whereas
\code{\link[=image_lat]{image_lat()}} performs local adaptive thresholding.
}
\details{
\itemize{
\item \code{image_threshold(type = "black")}: Forces all pixels below the threshold into black while leaving all pixels
at or above the threshold unchanged
\item \code{image_threshold(type = "white")}: Forces all pixels above the threshold into white while leaving all pixels
at or below the threshold unchanged
\item \code{image_lat()}: Local Adaptive Thresholding. Looks in a box (width x height) around the
pixel neighborhood if the pixel value is bigger than the average minus an offset.
}
}
\examples{
test <- image_convert(logo, colorspace = "Gray")
image_threshold(test, type = "black", threshold = "50\%")
image_threshold(test, type = "white", threshold = "50\%")

# Turn image into BW
test \%>\%
  image_threshold(type = "white", threshold = "50\%") \%>\%
  image_threshold(type = "black", threshold = "50\%")

# adaptive thresholding
image_lat(test, geometry = '10x10+5\%')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/animation.R
\name{animation}
\alias{animation}
\alias{image_animate}
\alias{image_coalesce}
\alias{image_morph}
\alias{image_mosaic}
\alias{image_flatten}
\alias{image_average}
\alias{image_append}
\alias{image_apply}
\alias{image_montage}
\title{Image Frames and Animation}
\usage{
image_animate(
  image,
  fps = 10,
  delay = NULL,
  loop = 0,
  dispose = c("background", "previous", "none"),
  optimize = FALSE
)

image_coalesce(image)

image_morph(image, frames = 8)

image_mosaic(image, operator = NULL)

image_flatten(image, operator = NULL)

image_average(image)

image_append(image, stack = FALSE)

image_apply(image, FUN, ...)

image_montage(
  image,
  geometry = NULL,
  tile = NULL,
  gravity = "Center",
  bg = "white",
  shadow = FALSE
)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{fps}{frames per second. Ignored if \code{delay} is not \code{NULL}.}

\item{delay}{delay after each frame, in 1/100 seconds.
Must be length 1, or number of frames. If specified, then \code{fps} is ignored.}

\item{loop}{how many times to repeat the animation. Default is infinite.}

\item{dispose}{a frame \href{https://legacy.imagemagick.org/Usage/anim_basics/#dispose}{disposal method}
from \link[=dispose_types]{dispose_types()}}

\item{optimize}{optimize the \code{gif} animation by storing only the differences
between frames. Input images must be exactly the same size.}

\item{frames}{number of frames to use in output animation}

\item{operator}{string with a
\href{https://www.imagemagick.org/Magick++/Enumerations.html#CompositeOperator}{composite operator}
from \link[=compose_types]{compose_types()}}

\item{stack}{place images top-to-bottom (TRUE) or left-to-right (FALSE)}

\item{FUN}{a function to be called on each frame in the image}

\item{...}{additional parameters for \code{FUN}}

\item{geometry}{a \link{geometry} string that defines the size the individual
thumbnail images, and the spacing between them.}

\item{tile}{a \link{geometry} string for example "4x5 with limits on how the
tiled images are to be laid out on the final result.}

\item{gravity}{a gravity direction, if the image is smaller than the frame, where
in the frame is the image to be placed.}

\item{bg}{a background color string}

\item{shadow}{enable shadows between images}
}
\description{
Operations to manipulate or combine multiple frames of an image. Details below.
}
\details{
For details see \href{https://www.imagemagick.org/Magick++/STL.html}{Magick++ STL}
documentation. Short descriptions:
\itemize{
\item \link{image_animate} coalesces frames by playing the sequence and converting to \code{gif} format.
\item \link{image_morph} expands number of frames by interpolating intermediate frames to blend
into each other when played as an animation.
\item \link{image_mosaic} inlays images to form a single coherent picture.
\item \link{image_montage} creates a composite image by combining frames.
\item \link{image_flatten} merges frames as layers into a single frame using a given operator.
\item \link{image_average} averages frames into single frame.
\item \link{image_append} stack images left-to-right (default) or top-to-bottom.
\item \link{image_apply} applies a function to each frame
}

The \link{image_apply} function calls an image function to each frame and joins
results back into a single image. Because most operations are already vectorized
this is often not needed. Note that \code{FUN()} should return an image. To apply other
kinds of functions to image frames simply use \link{lapply}, \link{vapply}, etc.
}
\examples{
# Combine images
logo <- image_read("https://jeroen.github.io/images/Rlogo.png")
oldlogo <- image_read("https://jeroen.github.io/images/Rlogo-old.png")

# Create morphing animation
both <- image_scale(c(oldlogo, logo), "400")
image_average(image_crop(both))
image_animate(image_morph(both, 10))

# Create thumbnails from GIF
banana <- image_read("https://jeroen.github.io/images/banana.gif")
length(banana)
image_average(banana)
image_flatten(banana)
image_append(banana)
image_append(banana, stack = TRUE)

# Append images together
wizard <- image_read("wizard:")
image_append(image_scale(c(image_append(banana[c(1,3)], stack = TRUE), wizard)))

image_composite(banana, image_scale(logo, "300"))

# Break down and combine frames
front <- image_scale(banana, "300")
background <- image_background(image_scale(logo, "400"), 'white')
frames <- image_apply(front, function(x){image_composite(background, x, offset = "+70+30")})
image_animate(frames, fps = 10)
# Simple 4x3 montage
input <- rep(logo, 12)
image_montage(input, geometry = 'x100+10+10', tile = '4x3', bg = 'pink', shadow = TRUE)

# With varying frame size
input <- c(wizard, wizard, logo, logo)
image_montage(input, tile = '2x2', bg = 'pink', gravity = 'southwest')
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit.R
\name{editing}
\alias{editing}
\alias{image_read}
\alias{image_read_svg}
\alias{image_read_pdf}
\alias{image_read_video}
\alias{image_write}
\alias{image_convert}
\alias{image_data}
\alias{image_raster}
\alias{image_display}
\alias{image_browse}
\alias{image_strip}
\alias{image_blank}
\alias{image_destroy}
\alias{image_join}
\alias{image_attributes}
\alias{image_get_artifact}
\alias{demo_image}
\title{Image Editing}
\usage{
image_read(
  path,
  density = NULL,
  depth = NULL,
  strip = FALSE,
  coalesce = TRUE,
  defines = NULL
)

image_read_svg(path, width = NULL, height = NULL)

image_read_pdf(path, pages = NULL, density = 300, password = "")

image_read_video(path, fps = 1, format = "png")

image_write(
  image,
  path = NULL,
  format = NULL,
  quality = NULL,
  depth = NULL,
  density = NULL,
  comment = NULL,
  flatten = FALSE,
  defines = NULL,
  compression = NULL
)

image_convert(
  image,
  format = NULL,
  type = NULL,
  colorspace = NULL,
  depth = NULL,
  antialias = NULL,
  matte = NULL,
  interlace = NULL
)

image_data(image, channels = NULL, frame = 1)

image_raster(image, frame = 1, tidy = TRUE)

image_display(image, animate = TRUE)

image_browse(image, browser = getOption("browser"))

image_strip(image)

image_blank(width, height, color = "none", pseudo_image = "", defines = NULL)

image_destroy(image)

image_join(...)

image_attributes(image)

image_get_artifact(image, artifact = "")

demo_image(path)
}
\arguments{
\item{path}{a file, url, or raster object or bitmap array}

\item{density}{resolution to render pdf or svg}

\item{depth}{color depth (either 8 or 16)}

\item{strip}{drop image comments and metadata}

\item{coalesce}{automatically \code{\link[=image_coalesce]{image_coalesce()}} gif images}

\item{defines}{a named character vector with extra options to control reading.
These are the \verb{-define key\{=value\}} settings in the \href{http://www.imagemagick.org/script/command-line-options.php#define}{command line tool}.
Use an empty string for value-less defines, and NA to unset a define.}

\item{width}{in pixels}

\item{height}{in pixels}

\item{pages}{integer vector with page numbers. Defaults to all pages.}

\item{password}{user \link[pdftools:pdf_render_page]{password} to open protected pdf files}

\item{fps}{how many images to capture per second of video. Set to
\code{NULL} to get all frames from the input video.}

\item{format}{output format such as \code{"png"}, \code{"jpeg"}, \code{"gif"}, \code{"rgb"} or \code{"rgba"}.}

\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{quality}{number between 0 and 100 for jpeg quality. Defaults to 75.}

\item{comment}{text string added to the image metadata for supported formats}

\item{flatten}{should image be flattened before writing? This also replaces
transparency with background color.}

\item{compression}{a string with compression type from \link{compress_types}}

\item{type}{string with \href{https://www.imagemagick.org/Magick++/Enumerations.html#ImageType}{imagetype}
value from \link{image_types} for example \code{grayscale} to convert into black/white}

\item{colorspace}{string with a \href{https://www.imagemagick.org/Magick++/Enumerations.html#ColorspaceType}{\code{colorspace}}
from \link{colorspace_types} for example \code{"gray"}, \code{"rgb"} or \code{"cmyk"}}

\item{antialias}{enable anti-aliasing for text and strokes}

\item{matte}{set to \code{TRUE} or \code{FALSE} to enable or disable transparency}

\item{interlace}{string with \href{https://www.imagemagick.org/Magick++/Enumerations.html#InterlaceType}{interlace}}

\item{channels}{string with image channel(s) for example \code{"rgb"}, \code{"rgba"},
\code{"cmyk"},\code{"gray"}, or \code{"ycbcr"}. Default is either \code{"gray"}, \code{"rgb"} or \code{"rgba"}
depending on the image}

\item{frame}{integer setting which frame to extract from the image}

\item{tidy}{converts raster data to long form for use with \link[ggplot2:geom_tile]{geom_raster}.
If \code{FALSE} output is the same as \code{as.raster()}.}

\item{animate}{support animations in the X11 display}

\item{browser}{argument passed to \link[utils:browseURL]{browseURL}}

\item{color}{a valid \href{https://www.imagemagick.org/Magick++/Color.html}{color string} such as
\code{"navyblue"} or \code{"#000080"}. Use \code{"none"} for transparency.}

\item{pseudo_image}{string with \href{http://www.imagemagick.org/script/formats.php#pseudo}{pseudo image}
specification for example \code{"radial-gradient:purple-yellow"}}

\item{...}{several images or lists of images to be combined}

\item{artifact}{string with name of the artifact to extract, see the
\link{image_deskew} for an example.}
}
\description{
Read, write and join or combine images. All image functions are vectorized, meaning
they operate either on a single frame or a series of frames (e.g. a collage, video,
or animation). Besides paths and URLs, \code{\link[=image_read]{image_read()}} supports commonly used bitmap
and raster object types.
}
\details{
All standard base vector methods such as \link{[}, \link{[[}, \code{\link[=c]{c()}}, \code{\link[=as.list]{as.list()}},
\code{\link[=as.raster]{as.raster()}}, \code{\link[=rev]{rev()}}, \code{\link[=length]{length()}}, and \code{\link[=print]{print()}} can be used to work with magick
image objects. Use the standard \code{img[i]} syntax to extract a subset of the frames
from an image. The \code{img[[i]]} method is an alias for \code{\link[=image_data]{image_data()}} which extracts
a single frame as a raw bitmap matrix with pixel values.

For reading svg or pdf it is recommended to use \code{image_read_svg()} and \code{image_read_pdf()}
if the \link[rsvg:rsvg]{rsvg} and \link[pdftools:pdf_render_page]{pdftools} R packages are available.
These functions provide more rendering options (including rendering of literal svg) and
better quality than built-in svg/pdf rendering delegates from imagemagick itself.

X11 is required for \code{image_display()} which is only works on some platforms. A more
portable method is \code{image_browse()} which opens the image in a browser. RStudio has
an embedded viewer that does this automatically which is quite nice.

Image objects are automatically released by the garbage collector when they are no longer
reachable. Because the GC only runs once in a while, you can also call \code{image_destroy()}
explicitly to release the memory immediately. This is usually only needed if you create
a lot of images in a short period of time, and you might run out of memory.
}
\examples{
# Download image from the web
frink <- image_read("https://jeroen.github.io/images/frink.png")
worldcup_frink <- image_fill(frink, "orange", "+100+200", 20)
image_write(worldcup_frink, "output.png")

# extract raw bitmap array
bitmap <- frink[[1]]

# replace pixels with #FF69B4 ('hot pink') and convert back to image
bitmap[,50:100, 50:100] <- as.raw(c(0xff, 0x69, 0xb4, 0xff))
image_read(bitmap)

# Plot to graphics device via legacy raster format
raster <- as.raster(frink)
par(ask=FALSE)
plot(raster)

# Read bitmap arrays from from other image packages
curl::curl_download("https://jeroen.github.io/images/example.webp", "example.webp")
if(require(webp)) image_read(webp::read_webp("example.webp"))
unlink(c("example.webp", "output.png"))
if(require(rsvg)){
tiger <- image_read_svg("http://jeroen.github.io/images/tiger.svg")
svgtxt <- '<?xml version="1.0" encoding="UTF-8"?>
<svg width="400" height="400" viewBox="0 0 400 400" fill="none">
 <circle fill="steelblue" cx="200" cy="200" r="100" />
 <circle fill="yellow" cx="200" cy="200" r="90" />
</svg>'
circles <- image_read_svg(svgtxt)
}
if(require(pdftools))
image_read_pdf(file.path(R.home('doc'), 'NEWS.pdf'), pages = 1, density = 100)
# create a solid canvas
image_blank(600, 400, "green")
image_blank(600, 400, pseudo_image = "radial-gradient:purple-yellow")
image_blank(200, 200, pseudo_image = "gradient:#3498db-#db3a34",
  defines = c('gradient:direction' = 'east'))
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/video.R
\name{video}
\alias{video}
\alias{image_write_video}
\alias{image_write_gif}
\title{Write Video}
\usage{
image_write_video(image, path = NULL, framerate = 10, ...)

image_write_gif(image, path = NULL, delay = 1/10, ...)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{path}{filename of the output gif or video. This is also the return value.}

\item{framerate}{frames per second, passed to \link[av:encoding]{av_encode_video}}

\item{...}{additional parameters passed to \link[av:encoding]{av_encode_video} and
\link[gifski:gifski]{gifski}.}

\item{delay}{duration of each frame in seconds (inverse of framerate)}
}
\description{
High quality video / gif exporter based on external packages \link[gifski:gifski]{gifski}
and \link[av:encoding]{av}.
}
\details{
This requires an image with multiple frames. The GIF exporter accomplishes the same
thing as \link{image_animate} but much faster and with better quality.
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fx.R
\name{fx}
\alias{fx}
\alias{image_fx}
\alias{image_fx_sequence}
\title{Image FX}
\usage{
image_fx(image, expression = "p", channel = NULL)

image_fx_sequence(image, expression = "p")
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{expression}{string with an \href{https://www.imagemagick.org/script/fx.php}{fx expression}}

\item{channel}{a value of \code{\link[=channel_types]{channel_types()}} specifying which channel(s) to set}
}
\description{
Apply a custom an \href{https://www.imagemagick.org/script/fx.php}{fx expression} to the image.
}
\details{
There are two different interfaces. The \link{image_fx} function simply applies
the same fx to each frame in the input image. The \link{image_fx_sequence} function
on the other hand treats the entire input vector as a sequence, allowing you
to apply an expression with multiple input images. See examples.
}
\examples{
# Show image_fx() expression
img <- image_convert(logo, colorspace = "Gray")
gradient_x <- image_convolve(img, kernel = "Prewitt")
gradient_y <- image_convolve(img, kernel = "Prewitt:90")
gradient <- c(image_fx(gradient_x, expression = "p^2"),
                image_fx(gradient_y, expression = "p^2"))
gradient <- image_flatten(gradient, operator = "Plus")
#gradient <- image_fx(gradient, expression = "sqrt(p)")
gradient

\donttest{
image_fx(img, expression = "pow(p, 0.5)")
image_fx(img, expression = "rand()")
}
# Use multiple source images
\donttest{
input <- c(logo, image_flop(logo))
image_fx_sequence(input, "(u+v)/2")
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EBImage.R
\name{as_EBImage}
\alias{as_EBImage}
\title{Convert to EBImage}
\usage{
as_EBImage(image)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}
}
\description{
Convert a Magick image to \href{https://bioconductor.org/packages/release/bioc/html/EBImage.html}{EBImage}
class. Note that EBImage only supports multi-frame images in greyscale.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/morphology.R
\name{morphology}
\alias{morphology}
\alias{image_morphology}
\alias{image_convolve}
\title{Morphology}
\usage{
image_morphology(
  image,
  method = "convolve",
  kernel = "Gaussian",
  iterations = 1,
  opts = list()
)

image_convolve(
  image,
  kernel = "Gaussian",
  iterations = 1,
  scaling = NULL,
  bias = NULL
)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{method}{a string with a valid method from \code{\link[=morphology_types]{morphology_types()}}}

\item{kernel}{either a square matrix or a string. The string can either be a
parameterized \link[=kernel_types]{kerneltype} such as: \code{"DoG:0,0,2"} or \code{"Diamond"}
or it can contain a custom matrix (see examples)}

\item{iterations}{number of iterations}

\item{opts}{a named list or character vector with custom attributes}

\item{scaling}{string with kernel scaling. The special flag \code{"!"} automatically scales to full
dynamic range, for example: \code{"50\%!"}}

\item{bias}{output bias string, for example \code{"50\%"}}
}
\description{
Apply a morphology method. This is a very flexible function which
can be used to apply any morphology method with custom parameters.
See \href{https://legacy.imagemagick.org/Usage/morphology/}{imagemagick website}
for examples.
}
\examples{
#example from IM website:
if(magick_config()$version > "6.8.8"){
pixel <- image_blank(1, 1, 'white') \%>\% image_border('black', '5x5')

# See the effect of Dilate method
pixel \%>\% image_scale('800\%')
pixel \%>\% image_morphology('Dilate', "Diamond") \%>\% image_scale('800\%')

# These produce the same output:
pixel \%>\% image_morphology('Dilate', "Diamond", iter = 3) \%>\% image_scale('800\%')
pixel \%>\% image_morphology('Dilate', "Diamond:3") \%>\% image_scale('800\%')

# Plus example
pixel \%>\% image_morphology('Dilate', "Plus", iterations = 2) \%>\% image_scale('800\%')

# Rose examples
rose \%>\% image_morphology('ErodeI', 'Octagon', iter = 3)
rose \%>\% image_morphology('DilateI', 'Octagon', iter = 3)
rose \%>\% image_morphology('OpenI', 'Octagon', iter = 3)
rose \%>\% image_morphology('CloseI', 'Octagon', iter = 3)

# Edge detection
man <- demo_image('man.gif')
man \%>\% image_morphology('EdgeIn', 'Octagon')
man \%>\% image_morphology('EdgeOut', 'Octagon')
man \%>\% image_morphology('Edge', 'Octagon')

# Octagonal Convex Hull
 man \%>\%
   image_morphology('Close', 'Diamond') \%>\%
   image_morphology('Thicken', 'ConvexHull', iterations = 1)

# Thinning down to a Skeleton
man \%>\% image_morphology('Thinning', 'Skeleton', iterations = 1)

# Specify custom kernel matrix usingn a string:
img <- demo_image("test_mag.gif")
i <- image_convolve(img, kernel = '4x5:
       0 -1  0  0
      -1 +1 -1  0
      -1 +1 -1  0
      -1 +1 +1 -1
       0 -1 -1  0 ', bias = "50\%")
}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{magrittr}{\code{\link[magrittr:pipe]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/device.R
\name{device}
\alias{device}
\alias{image_graph}
\alias{image_device}
\alias{image_draw}
\alias{image_capture}
\title{Magick Graphics Device}
\usage{
image_graph(
  width = 800,
  height = 600,
  bg = "white",
  pointsize = 12,
  res = 72,
  clip = TRUE,
  antialias = TRUE
)

image_draw(image, pointsize = 12, res = 72, antialias = TRUE, ...)

image_capture()
}
\arguments{
\item{width}{in pixels}

\item{height}{in pixels}

\item{bg}{background color}

\item{pointsize}{size of fonts}

\item{res}{resolution in pixels}

\item{clip}{enable clipping in the device. Because clipping can slow things down
a lot, you can disable it if you don't need it.}

\item{antialias}{TRUE/FALSE: enables anti-aliasing for text and strokes}

\item{image}{an existing image on which to start drawing}

\item{...}{additional device parameters passed to \link{plot.window} such as
\code{xlim}, \code{ylim}, or \code{mar}.}
}
\description{
Graphics device that produces a Magick image. Can either be used like a regular
device for making plots, or alternatively via \code{image_draw} to open a device
which draws onto an existing image using pixel coordinates. The latter is vectorized,
i.e. drawing operations are applied to each frame in the image.
}
\details{
The device is a relatively recent feature of the package. It should support all
operations but there might still be small inaccuracies. Also it is a bit slower than
some of the other devices, in particular for rendering text and clipping. Hopefully
this can be optimized in the next version.

By default \code{image_draw} sets all margins to 0 and uses graphics coordinates to
match image size in pixels (width x height) where \code{(0,0)} is the top left corner.
Note that this means the y axis increases from top to bottom which is the opposite
of typical graphics coordinates.  You can override all this by passing custom
\code{xlim}, \code{ylim} or \code{mar} values to \code{image_draw}.

The \code{image_capture} function returns the current device as an image. This only
works if the current device is a magick device or supports \link{dev.capture}.
}
\examples{
# Regular image
frink <- image_read("https://jeroen.github.io/images/frink.png")

# Produce image using graphics device
fig <- image_graph(res = 96)
ggplot2::qplot(mpg, wt, data = mtcars, colour = cyl)
dev.off()

# Combine
out <- image_composite(fig, frink, offset = "+70+30")
print(out)

# Or paint over an existing image
img <- image_draw(frink)
rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
abline(h = 300, col = 'blue', lwd = '10', lty = "dotted")
text(10, 250, "Hoiven-Glaven", family = "monospace", cex = 4, srt = 90)
palette(rainbow(11, end = 0.9))
symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
  bg = 1:11, inches = FALSE, add = TRUE)
dev.off()
print(img)

# Vectorized example with custom coordinates
earth <- image_read("https://jeroen.github.io/images/earth.gif")
img <- image_draw(earth, xlim = c(0,1), ylim = c(0,1))
rect(.1, .1, .9, .9, border = "red", lty = "dashed", lwd = 5)
text(.5, .9, "Our planet", cex = 3, col = "white")
dev.off()
print(img)
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wizard.R
\docType{data}
\name{wizard}
\alias{wizard}
\alias{logo}
\alias{rose}
\alias{granite}
\title{Example Images}
\format{
An object of class \code{magick-image} of length 1.
}
\usage{
logo
}
\description{
Example images included with ImageMagick:
}
\details{
\itemize{
\item \code{logo}: ImageMagick Logo, 640x480
\item \code{wizard}: ImageMagick Wizard, 480x640
\item \code{rose} : Picture of a rose, 70x46
\item \code{granite} : Granite texture pattern, 128x128
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/options.R
\name{options}
\alias{options}
\alias{magick_options}
\alias{option_types}
\alias{filter_types}
\alias{metric_types}
\alias{dispose_types}
\alias{compose_types}
\alias{colorspace_types}
\alias{channel_types}
\alias{image_types}
\alias{kernel_types}
\alias{noise_types}
\alias{gravity_types}
\alias{orientation_types}
\alias{morphology_types}
\alias{style_types}
\alias{decoration_types}
\alias{compress_types}
\alias{distort_types}
\title{Magick Options}
\usage{
magick_options()

option_types()

filter_types()

metric_types()

dispose_types()

compose_types()

colorspace_types()

channel_types()

image_types()

kernel_types()

noise_types()

gravity_types()

orientation_types()

morphology_types()

style_types()

decoration_types()

compress_types()

distort_types()
}
\description{
List option types and values supported in your version of ImageMagick. For
descriptions see
\href{https://www.imagemagick.org/Magick++/Enumerations.html}{ImageMagick Enumerations}.
}
\references{
ImageMagick Manual: \href{https://www.imagemagick.org/Magick++/Enumerations.html}{Enumerations}
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{color}},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/color.R
\name{color}
\alias{color}
\alias{image_modulate}
\alias{image_quantize}
\alias{image_map}
\alias{image_ordered_dither}
\alias{image_channel}
\alias{image_separate}
\alias{image_combine}
\alias{image_transparent}
\alias{image_background}
\alias{image_colorize}
\alias{image_contrast}
\alias{image_normalize}
\alias{image_enhance}
\alias{image_equalize}
\alias{image_median}
\title{Image Color}
\usage{
image_modulate(image, brightness = 100, saturation = 100, hue = 100)

image_quantize(
  image,
  max = 256,
  colorspace = "rgb",
  dither = TRUE,
  treedepth = NULL
)

image_map(image, map, dither = FALSE)

image_ordered_dither(image, threshold_map)

image_channel(image, channel = "lightness")

image_separate(image, channel = "default")

image_combine(image, colorspace = "sRGB", channel = "default")

image_transparent(image, color, fuzz = 0)

image_background(image, color, flatten = TRUE)

image_colorize(image, opacity, color)

image_contrast(image, sharpen = 1)

image_normalize(image)

image_enhance(image)

image_equalize(image)

image_median(image, radius = 1)
}
\arguments{
\item{image}{magick image object returned by \code{\link[=image_read]{image_read()}} or \code{\link[=image_graph]{image_graph()}}}

\item{brightness}{modulation of brightness as percentage of the current value (100 for no change)}

\item{saturation}{modulation of saturation as percentage of the current value (100 for no change)}

\item{hue}{modulation of hue is an absolute rotation of -180 degrees to +180 degrees from the
current position corresponding to an argument range of 0 to 200 (100 for no change)}

\item{max}{preferred number of colors in the image. The actual number of colors in the image may
be less than your request, but never more.}

\item{colorspace}{string with a \href{https://www.imagemagick.org/Magick++/Enumerations.html#ColorspaceType}{\code{colorspace}}
from \link{colorspace_types} for example \code{"gray"}, \code{"rgb"} or \code{"cmyk"}}

\item{dither}{a boolean (defaults to \code{TRUE}) specifying whether to apply Floyd/Steinberg error
diffusion to the image: averages intensities of several neighboring pixels}

\item{treedepth}{depth of the quantization color classification tree. Values of 0 or 1 allow
selection of the optimal tree depth for the color reduction algorithm. Values between 2 and 8
may be used to manually adjust the tree depth.}

\item{map}{reference image to map colors from}

\item{threshold_map}{A string giving the dithering pattern to use. See
\href{https://legacy.imagemagick.org/Usage/option_link.cgi?ordered-dither}{the ImageMagick documentation}
for possible values}

\item{channel}{a string with a
\href{https://www.imagemagick.org/Magick++/Enumerations.html#ChannelType}{channel} from
\link{channel_types} for example \code{"alpha"} or \code{"hue"} or \code{"cyan"}}

\item{color}{a valid \href{https://www.imagemagick.org/Magick++/Color.html}{color string} such as
\code{"navyblue"} or \code{"#000080"}. Use \code{"none"} for transparency.}

\item{fuzz}{relative color distance (value between 0 and 100) to be considered similar
in the filling algorithm}

\item{flatten}{should image be flattened before writing? This also replaces
transparency with background color.}

\item{opacity}{percentage of opacity used for coloring}

\item{sharpen}{enhance intensity differences in image}

\item{radius}{replace each pixel with the median color in a circular neighborhood}
}
\description{
Functions to adjust contrast, brightness, colors of the image. Details below.
}
\details{
For details see \href{https://www.imagemagick.org/Magick++/STL.html}{Magick++ STL}
documentation. Short descriptions:
\itemize{
\item \link{image_modulate} adjusts brightness, saturation and hue of image relative to current.
\item \link{image_quantize} reduces number of unique colors in the image.
\item \link{image_ordered_dither} reduces number of unique colors using a dithering threshold map.
\item \link{image_map} replaces colors of image with the closest color from a reference image.
\item \link{image_channel} extracts a single channel from an image and returns as grayscale.
\item \link{image_transparent} sets pixels approximately matching given color to transparent.
\item \link{image_background} sets background color. When image is flattened, transparent pixels get background color.
\item \link{image_colorize} overlays a solid color frame using specified opacity.
\item \link{image_contrast} enhances intensity differences in image
\item \link{image_normalize} increases contrast by normalizing the pixel values to span the full range of colors
\item \link{image_enhance} tries to minimize noise
\item \link{image_equalize} equalizes using histogram equalization
\item \link{image_median} replaces each pixel with the median color in a circular neighborhood
}

Note that
colors are also determined by image properties
\href{https://www.imagemagick.org/Magick++/Enumerations.html#ImageType}{imagetype} and
\href{https://www.imagemagick.org/Magick++/Enumerations.html#ColorspaceType}{colorspace}
which can be modified via \code{\link[=image_convert]{image_convert()}}.
}
\examples{
# manually adjust colors
logo <- image_read("logo:")
image_modulate(logo, brightness = 200)
image_modulate(logo, saturation = 150)
image_modulate(logo, hue = 200)

# Reduce image to 10 different colors using various spaces
image_quantize(logo, max = 10, colorspace = 'gray')
image_quantize(logo, max = 10, colorspace = 'rgb')
image_quantize(logo, max = 10, colorspace = 'cmyk')

image_ordered_dither(logo, 'o8x8')
# Change background color
translogo <- image_transparent(logo, 'white')
image_background(translogo, "pink", flatten = TRUE)

# Compare to flood-fill method:
image_fill(logo, "pink", fuzz = 20)

# Other color tweaks
image_colorize(logo, 50, "red")
image_contrast(logo)
image_normalize(logo)
image_enhance(logo)
image_equalize(logo)
image_median(logo)

# Alternate way to convert into black-white
image_convert(logo, type = 'grayscale')
}
\seealso{
Other image: 
\code{\link{_index_}},
\code{\link{analysis}},
\code{\link{animation}},
\code{\link{attributes}()},
\code{\link{composite}},
\code{\link{defines}},
\code{\link{device}},
\code{\link{edges}},
\code{\link{editing}},
\code{\link{effects}()},
\code{\link{fx}},
\code{\link{geometry}},
\code{\link{morphology}},
\code{\link{ocr}},
\code{\link{options}()},
\code{\link{painting}},
\code{\link{segmentation}},
\code{\link{transform}()},
\code{\link{video}}
}
\concept{image}
