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
