tidyHeatmap
================

<!-- badges: start -->

[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02472/status.svg)](https://doi.org/10.21105/joss.02472)
<!-- badges: end -->

## Citation

Mangiola et al., (2020). tidyHeatmap: an R package for modular heatmap
production based on tidy principles. Journal of Open Source Software,
5(52), 2472, <https://doi.org/10.21105/joss.02472>

Please have a look also to

-   [tidygate](https://github.com/stemangiola/tidygate/) for adding
    custom gate information to your tibble
-   [tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/)
    for tidy manipulation of Seurat objects
-   [tidyseurat](https://stemangiola.github.io/tidyseurat/) for tidy
    manipulation of Seurat objects
-   [tidybulk](https://stemangiola.github.io/tidybulk/) for tidy
    high-level data analysis and manipulation
-   [tidySummarizedExperiment](https://stemangiola.github.io/tidySummarizedExperiment/)
    for heatmaps produced with tidy principles

website:
[stemangiola.github.io/tidyHeatmap](https://stemangiola.github.io/tidyHeatmap/)

`tidyHeatmap` is a package that introduces tidy principles to the
creation of information-rich heatmaps. This package uses
[ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
as graphical engine.

**Advantages:**

-   Modular annotation with just specifying column names
-   Custom grouping of rows is easy to specify providing a grouped tbl.
    For example `df |> group_by(...)`
-   Labels size adjusted by row and column total number
-   Default use of Brewer and Viridis palettes

## Functions/utilities available

| Function           | Description                                 |
|--------------------|---------------------------------------------|
| `heatmap`          | Plots base heatmap                          |
| `add_tile`         | Adds tile annotation to the heatmap         |
| `add_point`        | Adds point annotation to the heatmap        |
| `add_bar`          | Adds bar annotation to the heatmap          |
| `add_line`         | Adds line annotation to the heatmap         |
| `layer_point`      | Adds layer of symbols on top of the heatmap |
| `layer_square`     | Adds layer of symbols on top of the heatmap |
| `layer_diamond`    | Adds layer of symbols on top of the heatmap |
| `layer_arrow_up`   | Adds layer of symbols on top of the heatmap |
| `layer_arrow_down` | Add layer of symbols on top of the heatmap  |
| `split_rows`       | Splits the rows based on the dendogram      |
| `split_columns`    | Splits the columns based on the dendogram   |
| `save_pdf`         | Saves the PDF of the heatmap                |

## Installation

To install the most up-to-date version

``` r
devtools::install_github("stemangiola/tidyHeatmap")
```

To install the most stable version (however please keep in mind that
this package is under a maturing lifecycle stage)

``` r
install.packages("tidyHeatmap")
```

## Contribution

If you want to contribute to the software, report issues or problems
with the software or seek support please open an issue
[here](https://github.com/stemangiola/tidyHeatmap/issues)

## Input data frame

The heatmaps visualise a multi-element, multi-feature dataset, annotated
with independent variables. Each observation is a element-feature pair
(e.g., person-physical characteristics).

| element         | feature         | value     | independent_variables |
|-----------------|-----------------|-----------|-----------------------|
| `chr` or `fctr` | `chr` or `fctr` | `numeric` | …                     |

Let’s transform the mtcars dataset into a tidy
“element-feature-independent variables” data frame. Where the
independent variables in this case are ‘hp’ and ‘vs’.

``` r
mtcars_tidy <- 
    mtcars |> 
    as_tibble(rownames="Car name") |> 
    
    # Scale
    mutate_at(vars(-`Car name`, -hp, -vs), scale) |>
    
    # tidyfy
    pivot_longer(cols = -c(`Car name`, hp, vs), names_to = "Property", values_to = "Value")

mtcars_tidy
```

    ## # A tibble: 288 × 5
    ##    `Car name`       hp    vs Property Value[,1]
    ##    <chr>         <dbl> <dbl> <chr>        <dbl>
    ##  1 Mazda RX4       110     0 mpg          0.151
    ##  2 Mazda RX4       110     0 cyl         -0.105
    ##  3 Mazda RX4       110     0 disp        -0.571
    ##  4 Mazda RX4       110     0 drat         0.568
    ##  5 Mazda RX4       110     0 wt          -0.610
    ##  6 Mazda RX4       110     0 qsec        -0.777
    ##  7 Mazda RX4       110     0 am           1.19 
    ##  8 Mazda RX4       110     0 gear         0.424
    ##  9 Mazda RX4       110     0 carb         0.735
    ## 10 Mazda RX4 Wag   110     0 mpg          0.151
    ## # … with 278 more rows

## Plotting

For plotting, you simply pipe the input data frame into heatmap,
specifying:

-   The rows, cols relative column names (mandatory)
-   The value column name (mandatory)
-   The annotations column name(s)

mtcars

``` r
mtcars_heatmap <- 
    mtcars_tidy |> 
    heatmap(`Car name`, Property, Value ) |>
    add_tile(hp)
```

    ## tidyHeatmap says: (once per session) from release 1.2.3 the grouping labels have white background by default. To add color for one-ay grouping specify palette_grouping = list(c("red", "blue"))

``` r
mtcars_heatmap
```

![](man/fragments/figures/unnamed-chunk-7-1.png)<!-- -->

## Saving

``` r
mtcars_heatmap |> save_pdf("mtcars_heatmap.pdf")
```

## Grouping and splitting

We can easily group the data (one group per dimension maximum, at the
moment only the vertical dimension is supported) with dplyr, and the
heatmap will be grouped accordingly

``` r
# Make up more groupings
mtcars_tidy_groupings = 
    mtcars_tidy |>
    mutate(property_group = if_else(Property %in% c("cyl", "disp"), "Engine", "Other"))

mtcars_tidy_groupings |> 
    group_by(vs, property_group) |>
    heatmap(`Car name`, Property, Value ) |>
    add_tile(hp)
```

![](man/fragments/figures/unnamed-chunk-9-1.png)<!-- -->

We can provide colour palettes to groupings

``` r
mtcars_tidy_groupings |> 
    group_by(vs, property_group) |>
    heatmap(
        `Car name`, Property, Value ,
        palette_grouping = list(
            
            # For first grouping (vs)
            c("#66C2A5", "#FC8D62"), 
            
            # For second grouping (property_group)
            c("#b58b4c", "#74a6aa")
        )
    ) |>
    add_tile(hp)
```

![](man/fragments/figures/unnamed-chunk-10-1.png)<!-- -->

We can split based on the cladogram

``` r
mtcars_tidy |> 
    heatmap(`Car name`, Property, Value ) |>
    split_rows(2) |>
    split_columns(2)
```

![](man/fragments/figures/unnamed-chunk-11-1.png)<!-- -->

We can split on kmean clustering (using ComplexHeatmap options, it is
stochastic)

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, Property, Value ,
        row_km = 2,
        column_km = 2
    ) 
```

![](man/fragments/figures/unnamed-chunk-12-1.png)<!-- -->

## Custom palettes

We can easily use custom palette, using strings, hexadecimal color
character vector,

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, 
        Property, 
        Value,
        palette_value = c("red", "white", "blue")
    )
```

![](man/fragments/figures/unnamed-chunk-13-1.png)<!-- -->

A better-looking blue-to-red palette

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, 
        Property, 
        Value,
        palette_value = circlize::colorRamp2(
            seq(-2, 2, length.out = 11), 
            RColorBrewer::brewer.pal(11, "RdBu")
        )
    )
```

![](man/fragments/figures/unnamed-chunk-14-1.png)<!-- -->

Or a grid::colorRamp2 function for higher flexibility

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, 
        Property, 
        Value,
        palette_value = circlize::colorRamp2(c(-2, -1, 0, 1, 2), viridis::magma(5))
    )
```

![](man/fragments/figures/unnamed-chunk-15-1.png)<!-- -->

## Multiple groupings and annotations

``` r
tidyHeatmap::pasilla |>
    group_by(location, type) |>
    heatmap(
        .column = sample,
        .row = symbol,
        .value = `count normalised adjusted`
    ) |>
    add_tile(condition) |>
    add_tile(activation)
```

![](man/fragments/figures/unnamed-chunk-16-1.png)<!-- -->

Remove legends, adding aesthetics to annotations in a modular fashion,
using `ComplexHeatmap` arguments

``` r
tidyHeatmap::pasilla |>
    group_by(location, type) |>
    heatmap(
        .column = sample,
        .row = symbol,
        .value = `count normalised adjusted`,
        show_heatmap_legend = FALSE
    ) |>
    add_tile(condition, show_legend = FALSE) |>
    add_tile(activation, show_legend = FALSE)
```

![](man/fragments/figures/unnamed-chunk-17-1.png)<!-- -->

## Annotation types

“tile”, “point”, “bar” and “line” are available

``` r
# Create some more data points
pasilla_plus <- 
    tidyHeatmap::pasilla |>
    dplyr::mutate(act = activation) |> 
    tidyr::nest(data = -sample) |>
    dplyr::mutate(size = rnorm(n(), 4,0.5)) |>
    dplyr::mutate(age = runif(n(), 50, 200)) |>
    tidyr::unnest(data) 

# Plot
pasilla_plus |>
    heatmap(
        .column = sample,
        .row = symbol,
        .value = `count normalised adjusted`
    ) |>
    add_tile(condition) |>
    add_point(activation) |>
    add_tile(act) |>
    add_bar(size) |>
    add_line(age)
```

![](man/fragments/figures/unnamed-chunk-18-1.png)<!-- -->

## Annotation size

We can customise annotation sizes using the `grid::unit()`, and the size
of their names using in-built `ComplexHeatmap` arguments

``` r
pasilla_plus |>
    heatmap(
        .column = sample,
        .row = symbol,
        .value = `count normalised adjusted`
    ) |>
    add_tile(condition, size = unit(0.3, "cm"), annotation_name_gp= gpar(fontsize = 8)) |>
    add_point(activation, size = unit(0.3, "cm"),   annotation_name_gp= gpar(fontsize = 8)) |>
    add_tile(act, size = unit(0.3, "cm"),   annotation_name_gp= gpar(fontsize = 8)) |>
    add_bar(size, size = unit(0.3, "cm"),   annotation_name_gp= gpar(fontsize = 8)) |>
    add_line(age, size = unit(0.3, "cm"),   annotation_name_gp= gpar(fontsize = 8))
```

![](man/fragments/figures/unnamed-chunk-19-1.png)<!-- -->

# Layer symbol

Add a layer on top of the heatmap

``` r
tidyHeatmap::pasilla |>
    
    # filter
    filter(symbol %in% head(unique(tidyHeatmap::pasilla$symbol), n = 10)) |>
    
    heatmap(
        .column = sample,
        .row = symbol,
        .value = `count normalised adjusted`
    ) |> 
    layer_point(
        `count normalised adjusted log` > 6 & sample == "untreated3" 
    )
```

![](man/fragments/figures/unnamed-chunk-20-1.png)<!-- -->

# ComplexHeatmap further styling

## Add cell borders

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, Property, Value, 
        rect_gp = grid::gpar(col = "#161616", lwd = 0.5)
    ) 
```

![](man/fragments/figures/unnamed-chunk-21-1.png)<!-- -->

## Drop row clustering

``` r
mtcars_tidy |> 
    heatmap(
        `Car name`, Property, Value, 
        cluster_rows = FALSE
    ) 
```

![](man/fragments/figures/unnamed-chunk-22-1.png)<!-- -->

## Reorder rows elements

``` r
library(forcats)
```

    ## Warning: package 'forcats' was built under R version 4.1.2

``` r
mtcars_tidy |> 
    mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
    heatmap(
        `Car name`, Property, Value, 
        cluster_rows = FALSE
    ) 
```

![](man/fragments/figures/unnamed-chunk-23-1.png)<!-- -->

## Size of dendrograms

``` r
mtcars_tidy |> 
    mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
    heatmap(
        `Car name`, Property, Value, 
        column_dend_height = unit(0.2, "cm"), 
        row_dend_width = unit(0.2, "cm")
    ) 
```

![](man/fragments/figures/unnamed-chunk-24-1.png)<!-- -->

## Size of rows/columns titles and names

``` r
mtcars_tidy |> 
    mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
    heatmap(
        `Car name`, Property, Value, 
        row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7),
        column_title_gp = gpar(fontsize = 7),
        row_title_gp = gpar(fontsize = 7)
    ) 
```

![](man/fragments/figures/unnamed-chunk-25-1.png)<!-- -->

## Using patchwork to integrate heatmaps

``` r
library(ggplot2)
library(patchwork)

p_heatmap =
    mtcars_tidy |> 
    heatmap(
        `Car name`, Property, Value, 
            show_heatmap_legend = FALSE,
        row_names_gp = gpar(fontsize = 7)
    ) 

p_ggplot = tibble(value = 1:10) %>% ggplot(aes(value)) + geom_density()

wrap_heatmap(p_heatmap) + 
    p_ggplot +
    wrap_heatmap(p_heatmap) + 
    plot_layout(width = c(1, 0.3, 1))
```

![](man/fragments/figures/unnamed-chunk-26-1.png)<!-- -->
---
title: 'tidyHeatmap: an R package for modular heatmap production based on tidy principles'
tags:
- R
- tidyverse
- tidy
- heatmap
- data visualization
date: "07 July 2020"
output:
  pdf_document: default
  word_document: default
  html_document:
    df_print: paged
authors:
- name: Stefano Mangiola^[Corresponding author]
  orcid: 0000-0001-7474-836X
  affiliation: 1, 2
- name: Anthony T. Papenfuss^[Corresponding author]
  orcid: 0000-0002-1102-8506
  affiliation: 1, 2, 3, 4
bibliography: paper.bib
affiliations:
- name: Bioinformatics Division, The Walter and Eliza Hall Institute of Medical Research,
    Parkville, Victoria, Australia
  index: 1
- name: Department of Medical Biology, University of Melbourne, Melbourne, Victoria,
    Australia.
  index: 2
- name: Peter MacCallum Cancer Centre, Melbourne, VIC 3000, Australia.
  index: 3
- name: School of Mathematics and Statistics, University of Melbourne, Melbourne,
    VIC 3010, Australia.
  index: 4
---

# Background
The heatmap is a powerful tool for visualising multi-dimensional data, where individual values can be organised in a two-dimensional matrix and their values expressed as colours. Rows and columns can be ordered according to their reciprocal similarity using hierarchical clustering, and dendrograms can be added to the plot to facilitate the interpretation. Row- and column-wise visual annotations, such as coloured tiles, can also be included. Within the R environment, several packages have been developed to produce heatmaps. The simplest and most readily available tool, `heatmap`, is provided within the `stats` package [@RCoreTeam:2013] and offers basic heatmaps with simple tile annotations. The versatile package `ggplot2` can be also used to produce basic heatmaps [@Hadley:2016]. More powerful software exists for producing fully annotated and/or multi-panel heatmaps, such as `Pheatmap` [@Kolde2012-tu], `superheat` [@superheatmap] and `ComplexHeatmap` [@Gu2016-cd]. The versatility of these packages comes at the cost of adding complexity to the user interface, characterised by many parameters and annotation functions that introduce a steep learning curve to produce complex, clear, and good-looking graphics.

Recently, efforts have been made toward the harmonisation of data frame structures and data analysis workflows using the concept of tidiness. Tidy data frames are characterised by having a specific structure where each variable is a column and each observation is a row. They provide ease of manipulation, modelling, and visualisation. The `tidyverse` is a suite of R libraries that define the standard for tidy data and APIs [@Hadley:2019]. The unique correspondence between quantities and annotations, characteristic of tidy data frames, allows complex operations to be performed from simple user inputs, such as a list of column names. 

# Statement of need
Considering (i) the utility and complexity of creating information-rich heatmaps, and (ii) the opportunity of increased coding efficiency and robustness offered by the tidy paradigm, a bridge between the two is very much needed. Recently, many tools for data science have been implemented according with tidy principles, this package aims to fill the gap for one of the most used data explorations tools.

# Tidy paradigm for visualisation
`tidyHeatmap` is a R package that introduces tidy principles to the creation of information-rich heatmaps. It is available in the CRAN R repository. This package currently uses <br>`ComplexHeatmap` as its graphical engine; however, due to its modular design it can be readily expanded to interface other engines. The command-line user interface is organised into (i) a main plotting utility; (ii) annotation layer utilities; and (iii) file-output utilities. The input is a tidy data frame with element (*e.g.*, person), feature (*e.g.*, physical characteristics) and value columns, with additional columns for independent variables for either elements (*e.g.*, number of sport medals) or features (*e.g.*, macroscopic or molecular characteristics). In this data structure, each observation is an element-feature pair. 

| element         | feature         | value     | annotation | group |
| --------------- | --------------- | --------- | ---------- | ----- |
| `chr` or `fctr` | `chr` or `fctr` | `numeric` | …          | …     |

The input data frame streams along the utility path using the pipe operator from <br>`magrittr`, allowing high modularity. The main utility allows the user to plot a base heatmap with dendrograms. The annotation utilities allow to serially add tile, point, bar and/or line annotation boxes to the side of the heatmap. The orientation of the annotations (row- or column-wise) is inferred based on the input data frame. The file-output utility allows the user to save vector or bitmap images directly from the R object, in the style of `ggplot2`. Row- or column-wise clusters can be defined effortlessly by applying the `group_by` function from `dplyr` [@Hadley:2020] to the input data frame. Data transformation and row/column scaling is done internally. Together, this leads to a decrease of coding burden of 3 and 5 folds for lines and characters respectively compared to `ComplexHeatmap` (*e.g.*, for \autoref{fig:example}). Besides offering a modular and user-friendly interface, `tidyHeatmap` applies publication-ready aesthetics such as `viridis` [@Garnier:2018] and `brewer` [@Neuwirth:2014] colour palettes and automatic sizing of row and column labels to avoid overlapping (\autoref{fig:example}). 

![Heatmap of the pasilla dataset including grouping and multiple annotations. Some annotation data was simulated for visualisation purposes. \label{fig:example}](paper_tables_and_figures_files/figure-gfm/example_figure-1.png)

The code interface consists of modular functions linked through the pipe operator. Custom colour palettes can be used by passing an array of colours or a colour function (*e.g.*, circlize [@Zuguang:2014]) to the palette argument of the annotation utilities.

```r
my_heatmap = 

	# Grouping
	input_df %>%
	group_by(pathway) %>%
		
	# Plotting
	heatmap(feature, element, value) %>%
    
	# Annotation
	add_tile(condition) %>%
	add_tile(act) %>%
	add_point(activation) %>%
	add_bar(size) %>%
	add_line(age)

# Saving
my_heatmap %>% save_pdf("my_file.pdf")
```

# Conclusions
In order to perform complex tasks, the use of disjointed data structures demands time consuming and bug-prone information matching. Joint, tidy data frames decrease the cost/benefit ratio for the user, automating a large part on the data manipulation. `tidyHeatmap` introduces a modular paradigm for specifying information-rich heatmaps, just requiring column names as input. Due to its intuitive user interface and its advanced default aesthetic features, `tidyHeatmap` is ideal for the quick production of publication-ready heatmaps. This software is designed for modular expandability. Future directions include the incorporation of more static and interactive heatmap visualisation engines.

# Acknowledgements

We acknowledge contributions all the Papenfuss Lab for feedback, and Maria Doyle for constant support.

# References

## Minimal input data frame

| element         | feature         | value     | annotation | group |
| --------------- | --------------- | --------- | ---------- | ----- |
| `chr` or `fctr` | `chr` or `fctr` | `numeric` | …          | …     |

``` r
library(tidyverse)
```

    ## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.1     ✓ purrr   0.3.4
    ## ✓ tibble  3.0.1     ✓ dplyr   1.0.0
    ## ✓ tidyr   1.1.0     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(tidyHeatmap)
```

    ## 
    ## Attaching package: 'tidyHeatmap'

    ## The following object is masked from 'package:stats':
    ## 
    ##     heatmap

``` r
pasilla_plus = 
    tidyHeatmap::pasilla %>%
        dplyr::mutate(act = activation) %>% 
        tidyr::nest(data = -sample) %>%
        dplyr::mutate(size = rnorm(n(), 4,0.5)) %>%
        dplyr::mutate(age = runif(n(), 50, 200)) %>%
        tidyr::unnest(data) %>%
        dplyr::rename(count = `count normalised adjusted`) %>%
        dplyr::mutate(pathway = if_else(location == "Secretory", "cluster 1", "cluster 2"))
```

``` r
pasilla_plus %>%
  group_by(pathway) %>%
  heatmap( symbol, sample, count ) %>%
    add_tile(condition) %>%
    add_point(activation) %>%
    add_tile(act) %>%
    add_bar(size) %>%
    add_line(age)
```

    ## Adding missing grouping variables: `pathway`
    ## Adding missing grouping variables: `pathway`
    ## Adding missing grouping variables: `pathway`
    ## Adding missing grouping variables: `pathway`

![](paper_tables_and_figures_files/figure-gfm/example_figure-1.png)<!-- -->
---
title: "tidyHeatmap"
output: github_document
---

```{r echo=FALSE}
knitr::opts_chunk$set( fig.path = "man/fragments/figures/")
```

<!-- badges: start -->
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02472/status.svg)](https://doi.org/10.21105/joss.02472)
<!-- badges: end -->

```{r child="man/fragments/intro.Rmd"}
```

---
output: github_document
---


## Minimal input data frame

element | feature | value | annotation | group 
------------ | ------------- | ------------- | ------------- | -------------
`chr` or `fctr` | `chr` or `fctr` | `numeric` | ... | ...

```{r}
library(tidyverse)
library(tidyHeatmap)
pasilla_plus = 
    tidyHeatmap::pasilla %>%
        dplyr::mutate(act = activation) %>% 
        tidyr::nest(data = -sample) %>%
        dplyr::mutate(size = rnorm(n(), 4,0.5)) %>%
        dplyr::mutate(age = runif(n(), 50, 200)) %>%
        tidyr::unnest(data) %>%
        dplyr::rename(count = `count normalised adjusted`) %>%
        dplyr::mutate(pathway = if_else(location == "Secretory", "cluster 1", "cluster 2"))
```

```{r example_figure}
pasilla_plus %>%
  group_by(pathway) %>%
  heatmap( symbol, sample, count ) %>%
    add_tile(condition) %>%
    add_point(activation) %>%
    add_tile(act) %>%
    add_bar(size) %>%
    add_line(age)
```
---
title: "Overview of the tidyHeatmap package"
author: "Stefano Mangiola"
date: "`r Sys.Date()`"
package: tidyHeatmap
output:
  html_vignette:
    toc_float: true
vignette: >
  %\VignetteIndexEntry{Overview of the tidyHeatmap package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

<!-- badges: start -->
[![Lifecycle:maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02472/status.svg)](https://doi.org/10.21105/joss.02472)
<!-- badges: end -->

```{r child="../man/fragments/intro.Rmd"}
```

# Session Info

```{r}
sessionInfo()
```

## Citation
Mangiola et al., (2020). tidyHeatmap: an R package for modular heatmap production based on tidy principles. Journal of Open Source Software, 5(52), 2472, https://doi.org/10.21105/joss.02472

Please have a look also to 

- [tidygate](https://github.com/stemangiola/tidygate/) for adding custom gate information to your tibble 
- [tidySingleCellExperiment](https://stemangiola.github.io/tidySingleCellExperiment/) for tidy manipulation of Seurat objects
- [tidyseurat](https://stemangiola.github.io/tidyseurat/) for tidy manipulation of Seurat objects
- [tidybulk](https://stemangiola.github.io/tidybulk/) for tidy high-level data analysis and manipulation 
- [tidySummarizedExperiment](https://stemangiola.github.io/tidySummarizedExperiment/) for heatmaps produced with tidy principles

website: [stemangiola.github.io/tidyHeatmap](https://stemangiola.github.io/tidyHeatmap/)

`tidyHeatmap` is a package that introduces tidy principles to the creation of information-rich heatmaps. 
This package uses [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html) as graphical engine.

**Advantages:**

  - Modular annotation with just specifying column names
  - Custom grouping of rows is easy to specify providing a grouped tbl.
    For example `df |> group_by(...)`
  - Labels size adjusted by row and column total number
  - Default use of Brewer and Viridis palettes

## Functions/utilities available

Function | Description
------------ | -------------
`heatmap` | Plots base heatmap
`add_tile` | Adds tile annotation to the heatmap
`add_point` | Adds point annotation to the heatmap
`add_bar` | Adds bar annotation to the heatmap
`add_line` | Adds line annotation to the heatmap
`layer_point` | Adds layer of symbols on top of the heatmap
`layer_square` | Adds layer of symbols on top of the heatmap
`layer_diamond` | Adds layer of symbols on top of the heatmap
`layer_arrow_up` | Adds layer of symbols on top of the heatmap
`layer_arrow_down` | Add layer of symbols on top of the heatmap
`split_rows` | Splits the rows based on the dendogram
`split_columns` | Splits the columns based on the dendogram
`save_pdf` | Saves the PDF of the heatmap

## Installation

To install the most up-to-date version

```{r, eval=FALSE}

devtools::install_github("stemangiola/tidyHeatmap")


```


To install the most stable version (however please keep in mind that this package is under a maturing lifecycle stage)

```{r, eval=FALSE}

install.packages("tidyHeatmap")

```

## Contribution

If you want to contribute to the software, report issues or problems with the software or seek support please open an issue [here](https://github.com/stemangiola/tidyHeatmap/issues)

## Input data frame

The heatmaps visualise a multi-element, multi-feature dataset, annotated with independent variables. Each observation is a element-feature pair (e.g., person-physical characteristics).

element | feature | value | independent_variables
------------ | ------------- | ------------- | -------------
`chr` or `fctr` | `chr` or `fctr` | `numeric` | ...

Let's transform the mtcars dataset into a tidy "element-feature-independent variables" data frame. Where the independent variables in this case are 'hp' and 'vs'.

```{r, echo=FALSE, include=FALSE}
library(dplyr)
library(tidyr)
library(tidyHeatmap)
library(grid)
```

```{r}
mtcars_tidy <- 
	mtcars |> 
	as_tibble(rownames="Car name") |> 
	
	# Scale
	mutate_at(vars(-`Car name`, -hp, -vs), scale) |>
	
	# tidyfy
	pivot_longer(cols = -c(`Car name`, hp, vs), names_to = "Property", values_to = "Value")

mtcars_tidy
```


## Plotting

For plotting, you simply pipe the input data frame into heatmap, specifying:

- The rows, cols relative column names (mandatory)
- The value column name (mandatory)
- The annotations column name(s)

mtcars
```{r}
mtcars_heatmap <- 
	mtcars_tidy |> 
	heatmap(`Car name`, Property, Value	) |>
	add_tile(hp)

mtcars_heatmap
```

## Saving

```{r eval=F}
mtcars_heatmap |> save_pdf("mtcars_heatmap.pdf")
```

## Grouping and splitting

We can easily group the data (one group per dimension maximum, at the moment only the vertical dimension is supported) with dplyr, and the heatmap will be grouped accordingly

```{r}
# Make up more groupings
mtcars_tidy_groupings = 
	mtcars_tidy |>
	mutate(property_group = if_else(Property %in% c("cyl", "disp"), "Engine", "Other"))

mtcars_tidy_groupings |> 
	group_by(vs, property_group) |>
	heatmap(`Car name`, Property, Value	) |>
	add_tile(hp)
```

We can provide colour palettes to groupings

```{r}
mtcars_tidy_groupings |> 
	group_by(vs, property_group) |>
	heatmap(
		`Car name`, Property, Value	,
		palette_grouping = list(
			
			# For first grouping (vs)
			c("#66C2A5", "#FC8D62"), 
			
			# For second grouping (property_group)
			c("#b58b4c", "#74a6aa")
		)
	) |>
	add_tile(hp)

```

We can split based on the cladogram

```{r}
mtcars_tidy |> 
	heatmap(`Car name`, Property, Value	) |>
	split_rows(2) |>
	split_columns(2)
```

We can split on kmean clustering (using ComplexHeatmap options, it is stochastic)

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, Property, Value	,
		row_km = 2,
		column_km = 2
	) 
```



## Custom palettes

We can easily use custom palette, using strings, hexadecimal color character vector, 

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, 
		Property, 
		Value,
		palette_value = c("red", "white", "blue")
	)
```

A better-looking blue-to-red palette

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, 
		Property, 
		Value,
		palette_value = circlize::colorRamp2(
			seq(-2, 2, length.out = 11), 
			RColorBrewer::brewer.pal(11, "RdBu")
		)
	)

```

Or a grid::colorRamp2 function for higher flexibility

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, 
		Property, 
		Value,
		palette_value = circlize::colorRamp2(c(-2, -1, 0, 1, 2), viridis::magma(5))
	)
```

## Multiple groupings and annotations

```{r}
tidyHeatmap::pasilla |>
	group_by(location, type) |>
	heatmap(
		.column = sample,
		.row = symbol,
		.value = `count normalised adjusted`
	) |>
	add_tile(condition) |>
	add_tile(activation)
```

Remove legends, adding aesthetics to annotations in a modular fashion, using `ComplexHeatmap` arguments

```{r}

tidyHeatmap::pasilla |>
	group_by(location, type) |>
	heatmap(
		.column = sample,
		.row = symbol,
		.value = `count normalised adjusted`,
		show_heatmap_legend = FALSE
	) |>
	add_tile(condition, show_legend = FALSE) |>
	add_tile(activation, show_legend = FALSE)
```

## Annotation types

"tile", "point", "bar" and "line" are available

```{r}
# Create some more data points
pasilla_plus <- 
	tidyHeatmap::pasilla |>
	dplyr::mutate(act = activation) |> 
	tidyr::nest(data = -sample) |>
	dplyr::mutate(size = rnorm(n(), 4,0.5)) |>
	dplyr::mutate(age = runif(n(), 50, 200)) |>
	tidyr::unnest(data) 

# Plot
pasilla_plus |>
	heatmap(
		.column = sample,
		.row = symbol,
		.value = `count normalised adjusted`
	) |>
	add_tile(condition) |>
	add_point(activation) |>
	add_tile(act) |>
	add_bar(size) |>
	add_line(age)
```

## Annotation size 

We can customise annotation sizes using the `grid::unit()`, and the size of their names using in-built `ComplexHeatmap` arguments

```{r}
pasilla_plus |>
	heatmap(
		.column = sample,
		.row = symbol,
		.value = `count normalised adjusted`
	) |>
	add_tile(condition, size = unit(0.3, "cm"),	annotation_name_gp= gpar(fontsize = 8)) |>
	add_point(activation, size = unit(0.3, "cm"),	annotation_name_gp= gpar(fontsize = 8)) |>
	add_tile(act, size = unit(0.3, "cm"),	annotation_name_gp= gpar(fontsize = 8)) |>
	add_bar(size, size = unit(0.3, "cm"),	annotation_name_gp= gpar(fontsize = 8)) |>
	add_line(age, size = unit(0.3, "cm"),	annotation_name_gp= gpar(fontsize = 8))
```

# Layer symbol

Add a layer on top of the heatmap

```{r}
tidyHeatmap::pasilla |>
	
	# filter
	filter(symbol %in% head(unique(tidyHeatmap::pasilla$symbol), n = 10)) |>
	
	heatmap(
		.column = sample,
		.row = symbol,
		.value = `count normalised adjusted`
	) |> 
	layer_point(
		`count normalised adjusted log` > 6 & sample == "untreated3" 
	)
```

# ComplexHeatmap further styling

## Add cell borders

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, Property, Value, 
		rect_gp = grid::gpar(col = "#161616", lwd = 0.5)
	) 
```

## Drop row clustering

```{r}
mtcars_tidy |> 
	heatmap(
		`Car name`, Property, Value, 
		cluster_rows = FALSE
	) 
```

## Reorder rows elements

```{r}
library(forcats)
mtcars_tidy |> 
	mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
	heatmap(
		`Car name`, Property, Value, 
		cluster_rows = FALSE
	) 
```

## Size of dendrograms

```{r}
mtcars_tidy |> 
	mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
	heatmap(
		`Car name`, Property, Value, 
		column_dend_height = unit(0.2, "cm"), 
		row_dend_width = unit(0.2, "cm")
	) 
```

## Size of rows/columns titles and names

```{r}
mtcars_tidy |> 
	mutate(`Car name` = fct_reorder(`Car name`, `Car name`, .desc = TRUE)) %>% 
	heatmap(
		`Car name`, Property, Value, 
		row_names_gp = gpar(fontsize = 7),
		column_names_gp = gpar(fontsize = 7),
		column_title_gp = gpar(fontsize = 7),
		row_title_gp = gpar(fontsize = 7)
	) 
```

## Using patchwork to integrate heatmaps

```{r}
library(ggplot2)
library(patchwork)

p_heatmap =
	mtcars_tidy |> 
	heatmap(
		`Car name`, Property, Value, 
			show_heatmap_legend = FALSE,
		row_names_gp = gpar(fontsize = 7)
	) 

p_ggplot = tibble(value = 1:10) %>% ggplot(aes(value)) + geom_density()

wrap_heatmap(p_heatmap) + 
	p_ggplot +
	wrap_heatmap(p_heatmap) + 
	plot_layout(width = c(1, 0.3, 1))

```
---
title: "Overview of the tidyHeatmap package"
author: "Stefano Mangiola"
date: "`r Sys.Date()`"
package: tidyHeatmap
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Overview of the tidyHeatmap package}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


```{r, echo=FALSE, include=FALSE, }
library(knitr)
knitr::opts_chunk$set(cache = TRUE, warning = FALSE, message = FALSE, cache.lazy = FALSE)

library(magrittr)
library(dplyr)
library(tidyr)
library(tidyHeatmap)
library(purrr)

```

Tidy heatmap. This package is a tidy wrapper of the package [ComplexHeatmap](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html). The goal of this package is to interface tidy data frames with this powerful tool.

**Advantages:**

- Modular annotation with just specifying column names
- Custom grouping of rows is easy to specify providing a grouped tbl. For example `df %>% group_by(...)`
- Labels size adjusted by row and column total number
- Default use of Brewer and Viridis palettes

## Functions/utilities available

Function | Description
------------ | -------------
`heatmap` | Plot base heatmap
`add_tile` | Add tile annotation to the heatmap
`add_point` | Add point annotation to the heatmap
`add_bar` | Add bar annotation to the heatmap
`add_line` | Add line annotation to the heatmap
`save_pdf` | Save the PDF of the heatmap

## Installation

To install the most up-to-date version

```{r, eval=FALSE}

devtools::install_github("stemangiola/tidyHeatmap")

```

To install the most stable version (however please keep in mind that this package is under a maturing lifecycle stage)

```{r, eval=FALSE}

install.packages("tidyHeatmap")

```

## Contribution

If you want to contribute to the software, report issues or problems with the software or seek support please open an issue [here](https://github.com/stemangiola/tidyHeatmap/issues)

## Input data frame

The heatmaps visualise a multi-element, multi-feature dataset, annotated with independent variables. Each observation is a element-feature pair (e.g., person-physical characteristics).

element | feature | value | independent_variables
------------ | ------------- | ------------- | -------------
`chr` or `fctr` | `chr` or `fctr` | `numeric` | ...

Let's transform the mtcars dataset into a tidy "element-feature-independent variables" data frame. Where the independent variables in this case are 'hp' and 'vs'.

```{r, echo=FALSE, include=FALSE}
library(dplyr)
library(tidyr)
library(tidyHeatmap)
```

```{r}
mtcars_tidy = 
	mtcars %>% 
	as_tibble(rownames="Car name") %>% 
	
	# Scale
	mutate_at(vars(-`Car name`, -hp, -vs), scale) %>%
	
	# tidyfy
	gather(Property, Value, -`Car name`, -hp, -vs)

mtcars_tidy
```


## Plot

For plotting, you simply pipe the input data frame into heatmap, specifying:

- The rows, cols relative column names (mandatory)
- The value column name (mandatory)
- The annotations column name(s)

mtcars
```{r}
mtcars_heatmap = 
	mtcars_tidy %>% 
		heatmap(`Car name`, Property, Value	) %>%
		add_tile(hp)

mtcars_heatmap
```

## Save

```{r eval=F}
mtcars_heatmap %>% save_pdf("mtcars_heatmap.pdf")
```

## Grouping

We can easily group the data (one group per dimension maximum, at the moment only the vertical dimension is supported) with dplyr, and the heatmap will be grouped accordingly

```{r}
mtcars_tidy %>% 
	group_by(vs) %>%
	heatmap(`Car name`, Property, Value	) %>%
	add_tile(hp)
```

## Custom palettes

We can easily use custom palette, using strings, hexadecimal color character vector, 

```{r}
mtcars_tidy %>% 
	heatmap(
		`Car name`, 
		Property, 
		Value,
		palette_value = c("red", "white", "blue")
	)
```

Or a grid::colorRamp2 function for higher flexibility

```{r}
mtcars_tidy %>% 
	heatmap(
		`Car name`, 
		Property, 
		Value,
		palette_value = circlize::colorRamp2(c(-2, -1, 0, 1, 2), viridis::magma(5))
	)
```

## Multiple groupings and annotations

```{r}
tidyHeatmap::pasilla %>%
	group_by(location, type) %>%
	heatmap(
			.column = sample,
			.row = symbol,
			.value = `count normalised adjusted`
		) %>%
	add_tile(condition) %>%
	add_tile(activation)
```

## Annotation types

**This feature requires >= 0.99.20 version**

"tile" (default), "point", "bar" and "line" are available

```{r}
# Create some more data points
pasilla_plus = 
	tidyHeatmap::pasilla %>%
		dplyr::mutate(act = activation) %>% 
		tidyr::nest(data = -sample) %>%
		dplyr::mutate(size = rnorm(n(), 4,0.5)) %>%
		dplyr::mutate(age = runif(n(), 50, 200)) %>%
		tidyr::unnest(data) 

# Plot
pasilla_plus %>%
		heatmap(
			.column = sample,
			.row = symbol,
			.value = `count normalised adjusted`
		) %>%
	add_tile(condition) %>%
	add_point(activation) %>%
	add_tile(act) %>%
	add_bar(size) %>%
	add_line(age)
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{input_heatmap}
\alias{input_heatmap}
\title{input_heatmap}
\usage{
input_heatmap(
  .data,
  .horizontal,
  .vertical,
  .abundance,
  transform = NULL,
  .scale = "row",
  palette_value = c("#440154FF", "#21908CFF", "#fefada"),
  palette_grouping = list(),
  ...
)
}
\arguments{
\item{.data}{A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |}

\item{.horizontal}{The name of the column horizontally presented in the heatmap}

\item{.vertical}{The name of the column vertically presented in the heatmap}

\item{.abundance}{The name of the transcript/gene abundance column}

\item{transform}{A function, used to transform .value, for example log1p}

\item{.scale}{A character string. Possible values are c(\"none\", \"row\", \"column\", \"both\")}

\item{palette_value}{A character vector, or a function for higher customisation (colorRamp2). This is the palette that will be used as gradient for abundance. If palette_value is a vector of hexadecimal colours, it should have 3 values. If you want more customisation, you can pass to palette_value a function, that is derived as for example `colorRamp2(c(-2, 0, 2), palette_value)`}

\item{palette_grouping}{A list of character vectors. This is the list of palettes that will be used for grouping}

\item{...}{Further arguments to be passed to ComplexHeatmap::Heatmap}
}
\value{
A `ComplexHeatmap` object
}
\description{
input_heatmap() takes a tbl object and easily produces a ComplexHeatmap plot, with integration with tibble and dplyr frameworks.
}
\details{
To be added.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/patchwork.R
\docType{methods}
\name{wrap_heatmap}
\alias{wrap_heatmap}
\alias{wrap_heatmap,InputHeatmap-method}
\title{Wrap tidyHeatmap (ComplexHeatmap) in a patchwork-compliant patch}
\usage{
wrap_heatmap(
  panel = NULL,
  plot = NULL,
  full = NULL,
  clip = TRUE,
  ignore_tag = FALSE
)

\S4method{wrap_heatmap}{InputHeatmap}(
  panel = NULL,
  plot = NULL,
  full = NULL,
  clip = TRUE,
  ignore_tag = FALSE
)
}
\arguments{
\item{panel, plot, full}{A grob, ggplot, patchwork, formula, raster, or
nativeRaster object to add to the respective area.}

\item{clip}{Should the grobs be clipped if expanding outside its area}

\item{ignore_tag}{Should tags be ignored for this patch. This is relevant
when using automatic tagging of plots and the content of the patch does not
qualify for a tag.}
}
\value{
A wrapped_patch object

A wrapped_patch object
}
\description{
In order to add tidyHeatmap (ComplexHeatmap) element to a patchwork they can be
converted to a compliant representation using the `wrap_heatmap()` function.
This allows you to position either grobs, ggplot objects, patchwork
objects, or even base graphics (if passed as a formula) in either the full
area, the full plotting area (anything between and
including the axis label), or the panel area (only the actual area where data
is drawn).
}
\examples{


tidyHeatmap::N52 |>
tidyHeatmap::heatmap(
 .row = symbol_ct,
 .column = UBR,
 .value = `read count normalised log`,
) |> 
wrap_heatmap()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{scale_robust}
\alias{scale_robust}
\title{Scale counts in a robust way against sd == 0}
\usage{
scale_robust(y)
}
\arguments{
\item{y}{A numerical array}
}
\value{
A scaled and centred numerical array
}
\description{
Scale counts in a robust way against sd == 0
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{scale_design}
\alias{scale_design}
\title{Scale design matrix}
\usage{
scale_design(df, .formula)
}
\arguments{
\item{df}{A tibble}

\item{.formula}{a formula}
}
\value{
A tibble
}
\description{
Scale design matrix
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{annot_to_list}
\alias{annot_to_list}
\title{annot_to_list}
\usage{
annot_to_list(.data)
}
\arguments{
\item{.data}{A data frame}
}
\value{
A list
}
\description{
annot_to_list
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{select_closest_pairs}
\alias{select_closest_pairs}
\title{Sub function of remove_redundancy_elements_though_reduced_dimensions}
\usage{
select_closest_pairs(df)
}
\arguments{
\item{df}{A tibble}
}
\value{
A tibble with pairs to drop
}
\description{
Sub function of remove_redundancy_elements_though_reduced_dimensions
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{layer_point}
\alias{layer_point}
\alias{layer_point,InputHeatmap-method}
\title{Adds a layers of symbols above the heatmap tiles to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
layer_point(.data, ...)

\S4method{layer_point}{InputHeatmap}(.data, ...)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{...}{Expressions that return a logical value, and are defined in terms of the variables in .data. If multiple expressions are included, they are combined with the & operator. Only rows for which all conditions evaluate to TRUE are kept.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
layer_point() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% layer_point()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{add_line}
\alias{add_line}
\alias{add_line,InputHeatmap-method}
\title{Adds a line annotation layer to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
add_line(.data, .column, palette = NULL, size = NULL, ...)

\S4method{add_line}{InputHeatmap}(.data, .column, palette = NULL, size = NULL, ...)
}
\arguments{
\item{.data}{A `tbl_df` formatted as | <ELEMENT> | <FEATURE> | <VALUE> | <...> |}

\item{.column}{Vector of quotes}

\item{palette}{A character vector of colors  This is the list of palettes that will be used for horizontal and vertical discrete annotations. The discrete classification of annotations depends on the column type of your input tibble (e.g., character and factor).}

\item{size}{A grid::unit object, e.g. unit(2, "cm"). This is the height or width of the annotation depending on the orientation.}

\item{...}{The arguments that will be passed to top_annotation or left_annotation of the ComplexHeatmap container}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
add_line() from a `InputHeatmap` object, adds a line annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% add_line()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_elements_features}
\alias{get_elements_features}
\title{Get column names either from user or from attributes}
\usage{
get_elements_features(.data, .element, .feature, of_samples = TRUE)
}
\arguments{
\item{.data}{A tibble}

\item{.element}{A character name of the sample column}

\item{.feature}{A character name of the transcript/gene column}

\item{of_samples}{A boolean}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{save_pdf}
\alias{save_pdf}
\title{Save plot on PDF file}
\usage{
save_pdf(
  .heatmap,
  filename,
  width = NULL,
  height = NULL,
  units = c("in", "cm", "mm")
)
}
\arguments{
\item{.heatmap}{A `Heatmap`}

\item{filename}{A character string. The name of the output file/path}

\item{width}{A `double`. Plot width}

\item{height}{A `double`. Plot height}

\item{units}{A character string. units ("in", "cm", or "mm")}
}
\value{
NA
}
\description{
save_pdf() takes as input a Heatmap from ComplexHeatmap and save it to PDF file
}
\details{
\lifecycle{maturing}

It simply save an `Heatmap` to a PDF file use pdf() function in the back end
}
\examples{


library(dplyr)
	tidyHeatmap::heatmap(
  dplyr::group_by(tidyHeatmap::pasilla,		location, type),
  .column = sample,
  .row = symbol,
  .value = `count normalised adjusted`,
 ) \%>\%
 save_pdf(tempfile())


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_x_y_annotation_columns}
\alias{get_x_y_annotation_columns}
\title{get_x_y_annotation_columns}
\usage{
get_x_y_annotation_columns(.data, .column, .row, .abundance)
}
\arguments{
\item{.data}{A `tbl` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |}

\item{.column}{The name of the column horizontally presented in the heatmap}

\item{.row}{The name of the column vertically presented in the heatmap}

\item{.abundance}{The name of the transcript/gene abundance column}
}
\value{
A list
}
\description{
get_x_y_annotation_columns
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_sample_transcript_counts}
\alias{get_sample_transcript_counts}
\title{Get column names either from user or from attributes}
\usage{
get_sample_transcript_counts(.data, .sample, .transcript, .abundance)
}
\arguments{
\item{.data}{A tibble}

\item{.sample}{A character name of the sample column}

\item{.transcript}{A character name of the transcript/gene column}

\item{.abundance}{A character name of the read count column}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{check_if_wrong_input}
\alias{check_if_wrong_input}
\title{Check whether there are NA counts}
\usage{
check_if_wrong_input(.data, list_input, expected_type)
}
\arguments{
\item{.data}{A tibble of read counts}

\item{list_input}{A list}

\item{expected_type}{A character string}
}
\value{
A tbl
}
\description{
Check whether there are NA counts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pasilla}
\alias{pasilla}
\title{Example data set Pasilla}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 504 rows and 8 columns.
}
\usage{
pasilla
}
\description{
Example data set Pasilla
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{check_if_duplicated_genes}
\alias{check_if_duplicated_genes}
\title{Check whether there are duplicated genes/transcripts}
\usage{
check_if_duplicated_genes(.data, .sample, .transcript, .abundance)
}
\arguments{
\item{.data}{A tibble of read counts}

\item{.sample}{A character name of the sample column}

\item{.transcript}{A character name of the transcript/gene column}

\item{.abundance}{A character name of the read count column}
}
\value{
A tbl
}
\description{
Check whether there are duplicated genes/transcripts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{layer_square}
\alias{layer_square}
\alias{layer_square,InputHeatmap-method}
\title{Adds a layers of symbols above the heatmap tiles to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
layer_square(.data, ...)

\S4method{layer_square}{InputHeatmap}(.data, ...)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{...}{Expressions that return a logical value, and are defined in terms of the variables in .data. If multiple expressions are included, they are combined with the & operator. Only rows for which all conditions evaluate to TRUE are kept.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
layer_square() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% layer_square()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{split_rows}
\alias{split_rows}
\alias{split_rows,InputHeatmap-method}
\alias{split_columns}
\alias{split_columns,InputHeatmap-method}
\title{Split the heatmap row-wise depending on the biggest branches in the cladogram.}
\usage{
split_rows(.data, number_of_groups)

\S4method{split_rows}{InputHeatmap}(.data, number_of_groups)

split_columns(.data, number_of_groups)

\S4method{split_columns}{InputHeatmap}(.data, number_of_groups)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{number_of_groups}{An integer. The number of groups to split the cladogram into.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
split_rows() from a `InputHeatmap` object, split the row cladogram.

split_columns() from a `InputHeatmap` object, split the column cladogram.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.

\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% split_rows(2)


library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% split_columns(2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ifelse_pipe}
\alias{ifelse_pipe}
\title{This is a generalisation of ifelse that accepts an object and return an objects}
\usage{
ifelse_pipe(.x, .p, .f1, .f2 = NULL)
}
\arguments{
\item{.x}{A tibble}

\item{.p}{A boolean}

\item{.f1}{A function}

\item{.f2}{A function}
}
\value{
A tibble
}
\description{
This is a generalisation of ifelse that accepts an object and return an objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{add_class}
\alias{add_class}
\title{Add class to abject}
\usage{
add_class(var, name)
}
\arguments{
\item{var}{A tibble}

\item{name}{A character name of the attribute}
}
\value{
A tibble with an additional attribute
}
\description{
Add class to abject
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{drop_class}
\alias{drop_class}
\title{Remove class to abject}
\usage{
drop_class(var, name)
}
\arguments{
\item{var}{A tibble}

\item{name}{A character name of the class}
}
\value{
A tibble with an additional attribute
}
\description{
Remove class to abject
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{prepend}
\alias{prepend}
\title{From rlang deprecated}
\usage{
prepend(x, values, before = 1)
}
\arguments{
\item{x}{An array}

\item{values}{An array}

\item{before}{A boolean}
}
\value{
An array
}
\description{
From rlang deprecated
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{save_pdf,InputHeatmap-method}
\alias{save_pdf,InputHeatmap-method}
\title{save_pdf}
\usage{
\S4method{save_pdf}{InputHeatmap}(
  .heatmap,
  filename,
  width = NULL,
  height = NULL,
  units = c("in", "cm", "mm")
)
}
\arguments{
\item{.heatmap}{A `Heatmap`}

\item{filename}{A character string. The name of the output file/path}

\item{width}{A `double`. Plot width}

\item{height}{A `double`. Plot height}

\item{units}{A character string. units ("in", "cm", or "mm")}
}
\description{
save_pdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{check_if_counts_is_na}
\alias{check_if_counts_is_na}
\title{Check whether there are NA counts}
\usage{
check_if_counts_is_na(.data, .abundance)
}
\arguments{
\item{.data}{A tibble of read counts}

\item{.abundance}{A character name of the read count column}
}
\value{
A tbl
}
\description{
Check whether there are NA counts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ifelse2_pipe}
\alias{ifelse2_pipe}
\title{This is a generalisation of ifelse that accepts an object and return an objects}
\usage{
ifelse2_pipe(.x, .p1, .p2, .f1, .f2, .f3 = NULL)
}
\arguments{
\item{.x}{A tibble}

\item{.p1}{A boolean}

\item{.p2}{ELSE IF condition}

\item{.f1}{A function}

\item{.f2}{A function}

\item{.f3}{A function}
}
\value{
A tibble
}
\description{
This is a generalisation of ifelse that accepts an object and return an objects
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{as_matrix}
\alias{as_matrix}
\title{Get matrix from tibble}
\usage{
as_matrix(tbl, rownames = NULL, do_check = TRUE)
}
\arguments{
\item{tbl}{A tibble}

\item{rownames}{A character string of the rownames}

\item{do_check}{A boolean}
}
\value{
A matrix
}
\description{
Get matrix from tibble
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_elements_features_abundance}
\alias{get_elements_features_abundance}
\title{Get column names either from user or from attributes}
\usage{
get_elements_features_abundance(
  .data,
  .element,
  .feature,
  .abundance,
  of_samples = TRUE
)
}
\arguments{
\item{.data}{A tibble}

\item{.element}{A character name of the sample column}

\item{.feature}{A character name of the transcript/gene column}

\item{.abundance}{A character name of the read count column}

\item{of_samples}{A boolean}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_sample_counts}
\alias{get_sample_counts}
\title{Get column names either from user or from attributes}
\usage{
get_sample_counts(.data, .sample, .abundance)
}
\arguments{
\item{.data}{A tibble}

\item{.sample}{A character name of the sample column}

\item{.abundance}{A character name of the read count column}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{add_point}
\alias{add_point}
\alias{add_point,InputHeatmap-method}
\title{Adds a point annotation layer to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
add_point(.data, .column, palette = NULL, size = NULL, ...)

\S4method{add_point}{InputHeatmap}(.data, .column, palette = NULL, size = NULL, ...)
}
\arguments{
\item{.data}{A `tbl_df` formatted as | <ELEMENT> | <FEATURE> | <VALUE> | <...> |}

\item{.column}{Vector of quotes}

\item{palette}{A character vector of colors  This is the list of palettes that will be used for horizontal and vertical discrete annotations. The discrete classification of annotations depends on the column type of your input tibble (e.g., character and factor).}

\item{size}{A grid::unit object, e.g. unit(2, "cm"). This is the height or width of the annotation depending on the orientation.}

\item{...}{The arguments that will be passed to top_annotation or left_annotation of the ComplexHeatmap container}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
add_point() from a `InputHeatmap` object, adds a point annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% add_point()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{add_annotation}
\alias{add_annotation}
\title{add_annotation}
\usage{
add_annotation(
  my_input_heatmap,
  annotation,
  type = rep("tile", length(quo_names(annotation))),
  palette_discrete = list(),
  palette_continuous = list(),
  size = NULL,
  ...
)
}
\arguments{
\item{my_input_heatmap}{A `InputHeatmap` formatted as | <SAMPLE> | <TRANSCRIPT> | <COUNT> | <...> |}

\item{annotation}{Vector of quotes}

\item{type}{A character vector of the set c(\"tile\", \"point\", \"bar\", \"line\")}

\item{palette_discrete}{A list of character vectors. This is the list of palettes that will be used for horizontal and vertical discrete annotations. The discrete classification of annotations depends on the column type of your input tibble (e.g., character and factor).}

\item{palette_continuous}{A list of character vectors. This is the list of palettes that will be used for horizontal and vertical continuous annotations. The continuous classification of annotations depends on the column type of your input tibble (e.g., integer, numerical, double).}

\item{size}{A grid::unit object, e.g. unit(2, "cm"). This is the height or width of the annotation depending on the orientation.}

\item{...}{The arguments that will be passed to top_annotation or left_annotation of the ComplexHeatmap container}
}
\value{
A `ComplexHeatmap` object
}
\description{
add_annotation() takes a tbl object and easily produces a ComplexHeatmap plot, with integration with tibble and dplyr frameworks.
}
\details{
To be added.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{N52}
\alias{N52}
\title{Example data set N52}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 520 rows and 15 columns.
}
\usage{
N52
}
\description{
Example data set N52
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{add_bar}
\alias{add_bar}
\alias{add_bar,InputHeatmap-method}
\title{Adds a bar annotation layer to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
add_bar(.data, .column, palette = NULL, size = NULL, ...)

\S4method{add_bar}{InputHeatmap}(.data, .column, palette = NULL, size = NULL, ...)
}
\arguments{
\item{.data}{A `tbl_df` formatted as | <ELEMENT> | <FEATURE> | <VALUE> | <...> |}

\item{.column}{Vector of quotes}

\item{palette}{A character vector of colors  This is the list of palettes that will be used for horizontal and vertical discrete annotations. The discrete classification of annotations depends on the column type of your input tibble (e.g., character and factor).}

\item{size}{A grid::unit object, e.g. unit(2, "cm"). This is the height or width of the annotation depending on the orientation.}

\item{...}{The arguments that will be passed to top_annotation or left_annotation of the ComplexHeatmap container}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
add_bar() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% add_bar()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{heatmap}
\alias{heatmap}
\alias{heatmap_}
\alias{heatmap,tbl-method}
\alias{heatmap,tbl_df-method}
\title{Creates a  `InputHeatmap` object from `tbl_df` on evaluation creates a `ComplexHeatmap`}
\usage{
heatmap(
  .data,
  .row,
  .column,
  .value,
  transform = NULL,
  .scale = "row",
  palette_value = c("#440154FF", "#21908CFF", "#fefada"),
  palette_grouping = list(),
  annotation = NULL,
  type = rep("tile", length(quo_names(annotation))),
  palette_discrete = list(),
  palette_continuous = list(),
  ...
)

heatmap_(
  .data,
  .row,
  .column,
  .value,
  transform = NULL,
  .scale = "row",
  palette_value = c("#440154FF", "#21908CFF", "#fefada"),
  palette_grouping = list(),
  annotation = NULL,
  type = rep("tile", length(quo_names(annotation))),
  palette_discrete = list(),
  palette_continuous = list(),
  ...
)

\S4method{heatmap}{tbl}(
  .data,
  .row,
  .column,
  .value,
  transform = NULL,
  .scale = "row",
  palette_value = c("#440154FF", "#21908CFF", "#fefada"),
  palette_grouping = list(),
  annotation = NULL,
  type = rep("tile", length(quo_names(annotation))),
  palette_discrete = list(),
  palette_continuous = list(),
  ...
)

\S4method{heatmap}{tbl_df}(
  .data,
  .row,
  .column,
  .value,
  transform = NULL,
  .scale = "row",
  palette_value = c("#440154FF", "#21908CFF", "#fefada"),
  palette_grouping = list(),
  annotation = NULL,
  type = rep("tile", length(quo_names(annotation))),
  palette_discrete = list(),
  palette_continuous = list(),
  ...
)
}
\arguments{
\item{.data}{A `tbl_df` formatted as | <ELEMENT> | <FEATURE> | <VALUE> | <...> |}

\item{.row}{The name of the column vertically presented in the heatmap}

\item{.column}{The name of the column horizontally presented in the heatmap}

\item{.value}{The name of the column for the value of the element/feature pair}

\item{transform}{A function, used to transform .value row-wise (e.g., transform = log1p)}

\item{.scale}{A character string. Possible values are c(\"none\", \"row\", \"column\", \"both\")}

\item{palette_value}{A character vector This is the palette that will be used as gradient for .value. For example c("red", "white", "blue"). For higher flexibility you can use circlize::colorRamp2\(c\(-2, -1, 0, 1, 2\), viridis::magma\(5\)\)}

\item{palette_grouping}{A list of character vectors. This is the list of palettes that will be used for grouping. For example list(RColorBrewer::brewer.pal(8, "Accent")) or list(c("#B3E2CD", "#FDCDAC", "#CBD5E8")) or list(c("black", "red"))}

\item{annotation}{DEPRECATED. please use the annotation functions add_* function \(\* one of tile, point, bar, line  \).}

\item{type}{DEPRECATED. please use the annotation functions add_* function \(\* one of tile, point, bar, line  \).}

\item{palette_discrete}{DEPRECATED. please use the annotation functions add_* function \(\* one of tile, point, bar, line  \).}

\item{palette_continuous}{DEPRECATED. please use the annotation functions add_* function \(\* one of tile, point, bar, line  \).}

\item{...}{The arguments that will be passed to the Heatmap function of ComplexHeatmap backend}
}
\value{
A `InputHeatmap` objects that gets evaluated to a `ComplexHeatmap` object

A `InputHeatmap` object

A `InputHeatmap` object

A `InputHeatmap` object
}
\description{
heatmap() takes a tbl object and easily produces a ComplexHeatmap plot, with integration with tibble and dplyr frameworks.
}
\details{
\lifecycle{maturing}

This function takes a tbl as an input and creates a `ComplexHeatmap` plot. The information is stored in a `InputHeatmap` object that is updated along the pipe statement, for example adding annotation layers.
}
\examples{

library(dplyr)

tidyHeatmap::N52 \%>\%
group_by( `Cell type`) \%>\%
tidyHeatmap::heatmap(
 .row = symbol_ct,
 .column = UBR,
 .value = `read count normalised log`,
)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{layer_diamond}
\alias{layer_diamond}
\alias{layer_diamond,InputHeatmap-method}
\title{Adds a layers of symbols above the heatmap tiles to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
layer_diamond(.data, ...)

\S4method{layer_diamond}{InputHeatmap}(.data, ...)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{...}{Expressions that return a logical value, and are defined in terms of the variables in .data. If multiple expressions are included, they are combined with the & operator. Only rows for which all conditions evaluate to TRUE are kept.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
layer_diamond() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% layer_diamond()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{add_tile}
\alias{add_tile}
\alias{add_tile,InputHeatmap-method}
\title{Adds a tile annotation layer to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
add_tile(.data, .column, palette = NULL, size = NULL, ...)

\S4method{add_tile}{InputHeatmap}(.data, .column, palette = NULL, size = NULL, ...)
}
\arguments{
\item{.data}{A `tbl_df` formatted as | <ELEMENT> | <FEATURE> | <VALUE> | <...> |}

\item{.column}{Vector of quotes}

\item{palette}{A character vector of colors  This is the list of palettes that will be used for horizontal and vertical discrete annotations. The discrete classification of annotations depends on the column type of your input tibble (e.g., character and factor).}

\item{size}{A grid::unit object, e.g. unit(2, "cm"). This is the height or width of the annotation depending on the orientation.}

\item{...}{The arguments that will be passed to top_annotation or left_annotation of the ComplexHeatmap container}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
add_tile() from a `InputHeatmap` object, adds a tile annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% add_tile(CAPRA_TOTAL)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{add_attr}
\alias{add_attr}
\title{Add attribute to abject}
\usage{
add_attr(var, attribute, name)
}
\arguments{
\item{var}{A tibble}

\item{attribute}{An object}

\item{name}{A character name of the attribute}
}
\value{
A tibble with an additional attribute
}
\description{
Add attribute to abject
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_sample_transcript}
\alias{get_sample_transcript}
\title{Get column names either from user or from attributes}
\usage{
get_sample_transcript(.data, .sample, .transcript)
}
\arguments{
\item{.data}{A tibble}

\item{.sample}{A character name of the sample column}

\item{.transcript}{A character name of the transcript/gene column}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_abundance_norm_if_exists}
\alias{get_abundance_norm_if_exists}
\title{Get column names either from user or from attributes}
\usage{
get_abundance_norm_if_exists(.data, .abundance)
}
\arguments{
\item{.data}{A tibble}

\item{.abundance}{A character name of the abundance column}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{save_pdf,Heatmap-method}
\alias{save_pdf,Heatmap-method}
\title{save_pdf}
\usage{
\S4method{save_pdf}{Heatmap}(
  .heatmap,
  filename,
  width = NULL,
  height = NULL,
  units = c("in", "cm", "mm")
)
}
\arguments{
\item{.heatmap}{A `Heatmap`}

\item{filename}{A character string. The name of the output file/path}

\item{width}{A `double`. Plot width}

\item{height}{A `double`. Plot height}

\item{units}{A character string. units ("in", "cm", or "mm")}
}
\description{
save_pdf
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{get_elements}
\alias{get_elements}
\title{Get column names either from user or from attributes}
\usage{
get_elements(.data, .element, of_samples = TRUE)
}
\arguments{
\item{.data}{A tibble}

\item{.element}{A character name of the sample column}

\item{of_samples}{A boolean}
}
\value{
A list of column enquo or error
}
\description{
Get column names either from user or from attributes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{quo_names}
\alias{quo_names}
\title{Convert array of quosure (e.g. c(col_a, col_b)) into character vector}
\usage{
quo_names(v)
}
\arguments{
\item{v}{A array of quosures (e.g. c(col_a, col_b))}
}
\value{
A character vector
}
\description{
Convert array of quosure (e.g. c(col_a, col_b)) into character vector
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{parse_formula}
\alias{parse_formula}
\title{.formula parser}
\usage{
parse_formula(fm)
}
\arguments{
\item{fm}{a formula}
}
\value{
A character vector
}
\description{
.formula parser
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{layer_arrow_down}
\alias{layer_arrow_down}
\alias{layer_arrow_down,InputHeatmap-method}
\title{Adds a layers of symbols above the heatmap tiles to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
layer_arrow_down(.data, ...)

\S4method{layer_arrow_down}{InputHeatmap}(.data, ...)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{...}{Expressions that return a logical value, and are defined in terms of the variables in .data. If multiple expressions are included, they are combined with the & operator. Only rows for which all conditions evaluate to TRUE are kept.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
layer_arrow_down() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% layer_arrow_down()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\docType{methods}
\name{layer_arrow_up}
\alias{layer_arrow_up}
\alias{layer_arrow_up,InputHeatmap-method}
\title{Adds a layers of symbols above the heatmap tiles to a `InputHeatmap`, that on evaluation creates a `ComplexHeatmap`}
\usage{
layer_arrow_up(.data, ...)

\S4method{layer_arrow_up}{InputHeatmap}(.data, ...)
}
\arguments{
\item{.data}{A `InputHeatmap`}

\item{...}{Expressions that return a logical value, and are defined in terms of the variables in .data. If multiple expressions are included, they are combined with the & operator. Only rows for which all conditions evaluate to TRUE are kept.}
}
\value{
A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`

A `InputHeatmap` object that gets evaluated to a `ComplexHeatmap`
}
\description{
layer_arrow_up() from a `InputHeatmap` object, adds a bar annotation layer.
}
\details{
\lifecycle{maturing}

It uses `ComplexHeatmap` as visualisation tool.
}
\examples{

library(dplyr)

hm = 
  tidyHeatmap::N52 \%>\%
  tidyHeatmap::heatmap(
    .row = symbol_ct,
    .column = UBR,
    .value = `read count normalised log`
)

hm \%>\% layer_arrow_up()


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{error_if_log_transformed}
\alias{error_if_log_transformed}
\title{Check whether a numeric vector has been log transformed}
\usage{
error_if_log_transformed(x, .abundance)
}
\arguments{
\item{x}{A numeric vector}

\item{.abundance}{A character name of the transcript/gene abundance column}
}
\value{
NA
}
\description{
Check whether a numeric vector has been log transformed
}
