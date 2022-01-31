
<!-- README.md is generated from README.Rmd. Please edit that file -->

# microViz <a href='https://david-barnett.github.io/microViz/index.html'><img src="man/figures/logo.png" align="right" height="180" width="156"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/david-barnett/microViz/workflows/R-CMD-check/badge.svg)](https://github.com/david-barnett/microViz/actions)
[![codecov](https://codecov.io/gh/david-barnett/microViz/branch/main/graph/badge.svg?token=C1EoVkhnxA)](https://codecov.io/gh/david-barnett/microViz)
[![Docker Cloud Build
Status](https://img.shields.io/docker/cloud/build/barnettdavid/microviz-rocker-verse)](https://hub.docker.com/r/barnettdavid/microviz-rocker-verse)
[![status](https://joss.theoj.org/papers/4547b492f224a26d96938ada81fee3fa/status.svg)](https://joss.theoj.org/papers/4547b492f224a26d96938ada81fee3fa)
[![DOI](https://zenodo.org/badge/307119750.svg)](https://zenodo.org/badge/latestdoi/307119750)

<!-- badges: end -->

## Overview

:package: `microViz` is an R package for analysis and visualization of
microbiome sequencing data.

:hammer: `microViz` functions are intended to be easy to use and
flexible.

:microscope: `microViz` extends and complements popular microbial
ecology packages like `phyloseq`, `vegan`, & `microbiome`.

## Learn more

:paperclip: This website is the best place for documentation and
examples: <https://david-barnett.github.io/microViz/>

-   [**This ReadMe**](https://david-barnett.github.io/microViz/) shows a
    few example analyses

-   **The
    [Reference](https://david-barnett.github.io/microViz/reference/index.html)
    page** lists all functions and links to help pages and examples

-   **The
    [Changelog](https://david-barnett.github.io/microViz/news/index.html)**
    describes important changes in new microViz package versions

-   **The Articles pages** give tutorials and further examples

    -   [Working with phyloseq
        objects](https://david-barnett.github.io/microViz/articles/web-only/phyloseq.html)

    -   [Fixing your taxa table with
        tax_fix](https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html)

    -   [Creating ordination
        plots](https://david-barnett.github.io/microViz/articles/web-only/ordination.html)
        (e.g. PCA or PCoA)

    -   [Interactive ordination plots with
        ord_explore](https://david-barnett.github.io/microViz/articles/web-only/ordination-interactive.html)

    -   [Visualising taxonomic compositions with
        comp_barplot](https://david-barnett.github.io/microViz/articles/web-only/compositions.html)

    -   [Heatmaps of microbiome composition and
        correlation](https://david-barnett.github.io/microViz/articles/web-only/heatmaps.html)

    -   [Modelling and plotting individual taxon associations with
        taxatrees](https://david-barnett.github.io/microViz/articles/web-only/modelling-taxa.html)

    -   More coming soon(ish)! Post on [GitHub
        discussions](https://github.com/david-barnett/microViz/discussions)
        if you have questions/requests

## Installation

You can install the latest available microViz package version using the
following instructions.

``` r
# Installing from github requires the devtools package
install.packages("devtools") 

# To install the latest "released" version of this package
devtools::install_github("david-barnett/microViz@0.9.0") # check 0.9.0 is the latest release

# To install the very latest version:
devtools::install_github("david-barnett/microViz")
# If you encounter a bug please try the latest version & let me know if the bug persists!

# If the Bioconductor dependencies don't automatically install you can install
# them yourself like this:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"))
```

:computer: **Windows users** - will need to have RTools installed so
that your computer can build this package (follow instructions here:
<http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/>)

**:apple: macOS** **users** - might need to install
[xquartz](https://www.xquartz.org/) to make the heatmaps work (to do
this with homebrew, run the following command in your mac’s Terminal:
`brew install --cask xquartz`

:package: I highly recommend using
[renv](https://rstudio.github.io/renv/index.html) for managing your R
package installations across multiple projects.

:whale: Alternatively, for Docker users an image with the main branch
installed is available at:
<https://hub.docker.com/r/barnettdavid/microviz-rocker-verse>

:date: microViz is tested to work with R version 4 on Windows, MacOS,
and Ubuntu 18 and 20. R version 3.6.\* should probably work, but I don’t
formally test this.

## Interactive ordination exploration

``` r
library(microViz)
```

microViz provides a Shiny app for an easy way to start exploring your
microbiome data: all you need is a phyloseq object.

``` r
# example data from corncob package
pseq <- corncob::ibd_phylo %>%
  tax_fix() %>%
  phyloseq_validate()
```

``` r
ord_explore(pseq) # gif generated with microViz version 0.7.4 (plays at 1.75x speed)
```

![](vignettes/web-only/images/20210429_ord_explore_x175.gif)

## Example analyses

``` r
library(phyloseq)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
```

``` r
# get some example data
data("dietswap", package = "microbiome")

# create a couple of numerical variables to use as constraints or conditions
dietswap <- dietswap %>%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  )
# add a couple of missing values to show how microViz handles missing data
sample_data(dietswap)$african[c(3, 4)] <- NA
```

### Looking at your data

You have quite a few samples in your phyloseq object, and would like to
visualise their compositions. Perhaps these example data differ
participant nationality?

``` r
dietswap %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  facet_wrap("nationality", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
#> Registered S3 method overwritten by 'seriation':
#>   method         from 
#>   reorder.hclust vegan
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
htmp <- dietswap %>%
  ps_mutate(nationality = as.character(nationality)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(dietswap, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(palette = viridis::turbo(11)),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Nationality = anno_sample_cat(
        var = "nationality", col = c(AAM = "grey35", AFR = "grey85"),
        box_col = NA, legend_title = "Nationality", size = grid::unit(4, "mm")
      )
    )
  )

ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

### Example ordination plot workflow

Ordination methods can also help you to visualise if overall microbial
ecosystem composition differs markedly between groups, e.g. BMI.

Here is one option as an example:

1.  Filter out rare taxa (e.g. remove Genera not present in at least 10%
    of samples) - use `tax_filter()`
2.  Aggregate the taxa into bacterial families (for example) - use
    `tax_agg()`
3.  Transform the microbial data with the centre-log-ratio
    transformation - use `tax_transform()`
4.  Perform PCA with the clr-transformed features (equivalent to
    aitchison distance PCoA) - use `ord_calc()`
5.  Plot the first 2 axes of this PCA ordination, colouring samples by
    group and adding taxon loading arrows to visualise which taxa
    generally differ across your samples - use `ord_plot()`
6.  Customise the theme of the ggplot as you like and/or add features
    like ellipses or annotations

``` r
# perform ordination
unconstrained_aitchison_pca <- dietswap %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>%
  tax_agg("Family") %>%
  tax_transform("clr") %>%
  ord_calc()
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.
# ord_calc will automatically infer you want a "PCA" here
# specify explicitly with method = "PCA", or you can pick another method

# create plot
pca_plot <- unconstrained_aitchison_pca %>%
  ord_plot(
    plot_taxa = 1:6, colour = "bmi_group", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 0.5),
    auto_caption = 8
  )

# customise plot
customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group), size = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 0.5, clip = "off") # makes rotated labels align correctly

# show plot
customised_plot
```

<img src="man/figures/README-ordination-plot-1.png" width="100%" />

### PERMANOVA

You visualised your ordinated data in the plot above. Below you can see
how to perform a PERMANOVA to test the significance of BMI’s association
with overall microbial composition. This example uses the Family-level
aitchison distance to correspond with the plot above.

``` r
# calculate distances
aitchison_dists <- dietswap %>%
  tax_filter(min_prevalence = 0.1) %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("aitchison")
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>% 
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99, # you should use at least 999!
    variables = "bmi_group"
  )
#> 2021-12-02 19:22:03 - Starting PERMANOVA with 99 perms with 1 processes
#> 2021-12-02 19:22:03 - Finished PERMANOVA

# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()
#>            Df  SumOfSqs         R2        F Pr(>F)
#> bmi_group   2  104.0678 0.04177157 4.773379   0.01
#> Residual  219 2387.2862 0.95822843       NA     NA
#> Total     221 2491.3540 1.00000000       NA     NA

# view the info stored about the distance calculation
info_get(aitchison_perm)
#> ps_extra info:
#> tax_agg = Family 
#> tax_transform = identity 
#> tax_scale = NA 
#> distMethod = aitchison 
#> ordMethod = NA 
#> constraints = NA 
#> conditions = NA
```

### Constrained partial ordination

You could visualise the effect of the (numeric/logical) variables in
your permanova directly using the `ord_plot` function with constraints
(and conditions).

``` r
perm2 <- aitchison_dists %>% 
  dist_permanova(variables = c("weight", "african", "sex"), seed = 321)
#> Dropping samples with missings: 2
#> 2021-12-02 19:22:03 - Starting PERMANOVA with 999 perms with 1 processes
#> 2021-12-02 19:22:04 - Finished PERMANOVA
```

We’ll visualise the effect of nationality and bodyweight on sample
composition, after first removing the effect of sex.

``` r
perm2 %>%
  ord_calc(constraints = c("weight", "african"), conditions = "female") %>%
  ord_plot(
    colour = "nationality", size = 2.5, alpha = 0.35,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(size = 1.5, colour = "grey15"),
    constraint_lab_style = constraint_lab_style(
      max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = nationality), size = 0.2) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed(ratio = 0.35, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
#> 
#> Centering (mean) and scaling (sd) the constraints and conditions:
#>  weight
#>  african
#>  female
```

<img src="man/figures/README-constrained-ord-plot-1.png" width="100%" />

### Correlation Heatmaps

microViz heatmaps are powered by `ComplexHeatmap` and annotated with
taxa prevalence and/or abundance.

``` r
# set up the data with numerical variables and filter to top taxa
psq <- dietswap %>%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  ) %>%
  tax_filter(
    tax_level = "Genus", min_prevalence = 1 / 10, min_sample_abundance = 1 / 10
  ) %>%
  tax_transform("identity", rank = "Genus")
#> Proportional min_prevalence given: 0.1 --> min 23/222 samples.

# randomly select 30 taxa from the 50 most abundant taxa (just for an example)
set.seed(123)
taxa <- sample(tax_top(psq, n = 50), size = 30)
# actually draw the heatmap
cor_heatmap(
  data = psq, taxa = taxa,
  taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
  tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(undetected = 50),
    Log2 = anno_tax_box(undetected = 50, trans = "log2", zero_replace = 1)
  )
)
```

<img src="man/figures/README-heatmap-1.png" width="100%" />

## Citation

:innocent: If you find any part of microViz useful to your work, please
consider citing the JOSS article:

Barnett et al., (2021). microViz: an R package for microbiome data
visualization and statistics. Journal of Open Source Software, 6(63),
3201, <https://doi.org/10.21105/joss.03201>

## Contributing

Bug reports, questions, suggestions for new features, and other
contributions are all welcome. Feel free to create a [GitHub
Issue](https://github.com/david-barnett/microViz/issues) or write on the
[Discussions](https://github.com/david-barnett/microViz/discussions)
page. Alternatively you could also contact me (David) on Twitter
[@\_david_barnett\_](https://twitter.com/_david_barnett_) .

This project is released with a [Contributor Code of
Conduct](https://david-barnett.github.io/microViz/CODE_OF_CONDUCT.html)
and by participating in this project you agree to abide by its terms.

## Session info

``` r
sessionInfo()
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.6 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
#>  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#>  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] ggplot2_3.3.5   dplyr_1.0.7     phyloseq_1.36.0 microViz_0.9.0 
#> [5] devtools_2.4.2  usethis_2.1.3   pkgdown_2.0.0  
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rtsne_0.15             colorspace_2.0-2       rjson_0.2.20          
#>   [4] ellipsis_0.3.2         rprojroot_2.0.2        circlize_0.4.13       
#>   [7] markdown_1.1           XVector_0.32.0         GlobalOptions_0.1.2   
#>  [10] fs_1.5.1               gridtext_0.1.4         ggtext_0.1.1          
#>  [13] clue_0.3-60            rstudioapi_0.13        farver_2.1.0          
#>  [16] remotes_2.4.1          fansi_0.5.0            xml2_1.3.2            
#>  [19] codetools_0.2-18       splines_4.1.2          doParallel_1.0.16     
#>  [22] cachem_1.0.6           knitr_1.36             pkgload_1.2.3         
#>  [25] ade4_1.7-18            jsonlite_1.7.2         Cairo_1.5-12.2        
#>  [28] cluster_2.1.2          png_0.1-7              compiler_4.1.2        
#>  [31] assertthat_0.2.1       Matrix_1.3-4           fastmap_1.1.0         
#>  [34] cli_3.1.0              htmltools_0.5.2        prettyunits_1.1.1     
#>  [37] tools_4.1.2            igraph_1.2.9           gtable_0.3.0          
#>  [40] glue_1.5.1             GenomeInfoDbData_1.2.7 reshape2_1.4.4        
#>  [43] Rcpp_1.0.7             Biobase_2.52.0         vctrs_0.3.8           
#>  [46] Biostrings_2.60.2      rhdf5filters_1.4.0     multtest_2.48.0       
#>  [49] ape_5.5                nlme_3.1-153           iterators_1.0.13      
#>  [52] xfun_0.28              stringr_1.4.0          ps_1.6.0              
#>  [55] testthat_3.1.0         lifecycle_1.0.1        zlibbioc_1.38.0       
#>  [58] MASS_7.3-54            scales_1.1.1           TSP_1.1-11            
#>  [61] parallel_4.1.2         biomformat_1.20.0      rhdf5_2.36.0          
#>  [64] RColorBrewer_1.1-2     ComplexHeatmap_2.8.0   yaml_2.2.1            
#>  [67] memoise_2.0.0          gridExtra_2.3          stringi_1.7.6         
#>  [70] highr_0.9              S4Vectors_0.30.2       desc_1.4.0            
#>  [73] foreach_1.5.1          permute_0.9-5          seriation_1.3.1       
#>  [76] BiocGenerics_0.38.0    pkgbuild_1.2.0         shape_1.4.6           
#>  [79] GenomeInfoDb_1.28.4    rlang_0.4.12           pkgconfig_2.0.3       
#>  [82] bitops_1.0-7           matrixStats_0.61.0     evaluate_0.14         
#>  [85] lattice_0.20-45        purrr_0.3.4            Rhdf5lib_1.14.2       
#>  [88] labeling_0.4.2         processx_3.5.2         tidyselect_1.1.1      
#>  [91] plyr_1.8.6             magrittr_2.0.1         R6_2.5.1              
#>  [94] magick_2.7.3           IRanges_2.26.0         generics_0.1.1        
#>  [97] DBI_1.1.1              pillar_1.6.4           withr_2.4.3           
#> [100] mgcv_1.8-38            survival_3.2-13        RCurl_1.98-1.5        
#> [103] tibble_3.1.6           corncob_0.2.0          crayon_1.4.2          
#> [106] utf8_1.2.2             microbiome_1.14.0      rmarkdown_2.11        
#> [109] GetoptLong_1.0.5       viridis_0.6.2          grid_4.1.2            
#> [112] data.table_1.14.2      callr_3.7.0            vegan_2.5-7           
#> [115] digest_0.6.29          tidyr_1.1.4            stats4_4.1.2          
#> [118] munsell_0.5.0          registry_0.5-1         viridisLite_0.4.0     
#> [121] sessioninfo_1.2.1
```
# microViz 0.9.0

## Heatmaps major changes
__`comp_heatmap` and `cor_heatmap` and their helpers are largely rewritten.__ 
For guidance, see the new website article on heatmaps.

- `taxAnnotation` and `varAnnotation` annotation helper functions replace the deprecated `tax_anno` and `var_anno` functions
- `sampleAnnotation` is added for coordinating sample annotations on `comp_heatmap`
- various `anno_*` helpers for each of the annotation coordinating functions above
- arguments for 

# microViz 0.8.2

## Features
- `taxatree_plots` can now show symbols indicating multiple levels of statistical significance
- `taxatree_plots`: it is now easier to use other colour palettes with their default luminance

# microViz 0.8.1

## Features
- `ord_explore` can perform "binary" transformations, unlocking interactive use of Binary Jaccard etc.
- `tax_transform` gains `add` argument to simply add a constant value to all otu_table values before transformation (as an alternative to `zero_replace`)
- `tax_scale` gains `keep_counts` argument for consistency with `tax_transform`

## Fixes
- `tax_model` now only requires `corncob` to be installed when actually used
- `ord_explore` barplot numerical inputs are now debounced to prevent lag from repeated redrawing 

# microViz 0.8.0 - "autumn leaves"

## Breaking changes

### Trees
__The taxatree_* family of functions are largely rewritten.__ For guidance, see the new website article on statistical modelling of taxa.

- `taxatree_models` now attaches resulting list to ps_extra
- `taxatree_models2stats` must be run on the output of `taxatree_models` before using `taxatree_plots`
- `taxatree_plots` has different arguments and can now be directly labelled with `taxatree_plot_labels` when `taxatree_label` is run first to identify which taxa to label.
- `taxatree_plotkey` has different arguments, with more flexible labelling conditions and a smarter label positioning approach.
- `tax_model` and `taxatree_models` now use "lm" type by default, instead of "bbdml", as `corncob` is only a suggested dependency.

### Ordination
- `ord_plot` default labels now have `alpha` = 1 for both taxa and constraints (previously 0.8)
- `ord_plot` "auto"matic loading/constraint vector length scalar adjustment improvement: now uses both axes

### Barplots
- `comp_barplot` now uses bray-curtis by default for sample ordering (instead of aitchison) as this generally looks better
- `comp_barplot` now expects palette argument colours in first-to-last order, which is more intuitive than the previous reverse order
- `distinct_palette` now adds "lightgrey" to end by default

### Other
- `tax_filter`'s `is_counts` argument replaced by `use_counts`, allowing it to filter ps_extra objects using stored count data (i.e. after `tax_transform`).
- `comp_heatmap` can no longer transform data internally, but accepts data already transformed with `tax_transform` and uses stored count data in the ps_extra for any taxa annotations
- `tax_anno` heatmap annotation default style slightly changed.
- `tax_names2rank` replaces deprecated `tax_names2tt`

## Features 
- `ord_plot` arrow labels can now be rotated with the help of `tax_lab_style()` and `constraint_lab_style()` 
- `ps_calc_dominant` function added, for conveniently identifying the dominant taxon in each phyloseq sample
- `distinct_palette` gains "kelly" and "greenArmytage" palettes and helpfully adds "lightgrey" to the end by default for convenient use as the palette argument to `comp_barplot`
- `tax_transform` can now chain multiple transformations together and records these transformations in the ps_extra output
- `tax_mutate` function added, for modifying the tax_table rank variables with `dplyr` `mutate` syntax
- `tax_sort` can now sort ps_extra objects
- heatmap annotation helper `tax_anno` no longer requires 'column' or 'row' to be specified in advance
- `prev`, a low level helper function for calculating taxon prevalence now exported
- Various phyloseq accessor functions now work with ps_extra objects e.g. taxa_names, sample_names

## Fixes
- `cor_heatmap` and `comp_heatmap` now respect column seriation arguments when different to row seriation
- `comp_barplot` now actually orders grouped samples by similarity AFTER splitting into groups, as documented

# microViz 0.7.10

Updated citation information for JOSS publication. No other changes.

# microViz 0.7.9

Release accompanying JOSS manuscript acceptance. 

- Includes a fix (hopefully temporary) for incorrect barplot legend in ord_explore caused by bug introduced by ggplot2 version 3.3.4 noted at https://github.com/tidyverse/ggplot2/issues/4511

# microViz 0.7.8

- `ord_plot` gains `vec_*` helper functions for generating lists for styling taxa and constraint vectors/arrows (`vec_constraint`, `vec_tax_all` and `vec_tax_sel`)

# microViz 0.7.7

- `ord_explore` shapes selection bug fixed by limiting to 5 shapes returned by new `scale_shape_girafe_filled` function

# microViz 0.7.6

## Features
- `stat_chulls` added, for drawing convex hulls on ord_plot and ord_explore ordinations
- `add_paths` added, for drawing geom_paths connecting a subset of samples (over time) on an ordination plot

## Fixes
- `phyloseq_validate` no longer checks "unique" tax_table rank column for NAs or nchar<4 
(avoiding warnings in ord_explore caused by short taxa_names)

# microViz 0.7.5

- `ord_explore` can now draw stat_ellipse or taxa loading vectors
- `tax_agg` error messages now include personalised suggested tax_fix code

# microViz 0.7.4

- `ord_explore` Shiny app GUI can now also be used to interactively generate ordination plots, and to generate `ord_plot` code
- `ord_plot` bug fix - can now plot any dimension
- Removed deprecated `tax_fill_unknowns` function

# microViz 0.7.3

- `phyloseq_validate` verbose = FALSE is actually silent now.

# microViz 0.7.2

- `ord_explore` now compatible with Shiny version >=1.5.0
- `tax_fix` now sends messages about fixing completely anonymous rows, instead of warnings

# microViz 0.7.1

Allows `ps_seriate`, `ps_arrange`, `ps_reorder`, `ps_mutate`, and `ps_select` to work directly with `ps_extra` objects, as this can be helpful when quickly exploring / printing aggregated data, as in the new "Working with phyloseq objects" tutorial.

# microViz 0.7.0 - "fickle fixes"

## Breaking changes
- `tax_fix` replaces the deprecated `tax_fill_unknowns`, `tax_fix` has all the same arguments except the 'levels' argument, which was removed

## Features
- `tax_fix_interactive` Shiny app will help you clean up your taxonomy table with `tax_fix`

# microViz 0.6.1

## Breaking changes
- `ord_plot_iris` and `ord_explore` no longer take ps argument of untransformed counts, because (by default) `tax_transform` now keeps the untransformed counts otu_table in the ps_extra object

## Features
- `ord_explore` now allows much better control over selection of points (using `ggiraph` functionality)
- `ord_plot` now has interactive option with `ggiraph` package
- `ord_plot_iris` gains ord_plot argument, allowing a simple pairing of iris plot and ordination to be made more easily
- `comp_barplot` (and by extension `ord_plot_iris`) can now be made interactive in a simple fashion, using ggiraph for hover/tooltip interaction with taxa

# microViz 0.6.0 - "open sesame"

This is the first public release version of microViz. It is still under active development, so pay attention to the following:

- Minor version changes (e.g. 0.5.* to 0.6.0) will signal that breaking changes have been made, i.e. installing the new version may break some previously working code. Breaking changes will be listed in this document.
- Patch versions (e.g. 0.5.1 to 0.5.2) will be used to release new features and bug fixes, which should not break existing code. Please let me know if it does anyway!

## Breaking changes
* `tax_agg` argument agg_level renamed to rank. `tax_agg` returns taxa in different order than before (and now different order from, but same aggregation as, `microbiome::aggregate_taxa()`). tax_agg now checks if taxa cannot be uniquely identified at the specified rank level. (now also about twice as fast)
* `tax_fill_unknowns` x arg renamed to ps. Also now stops when unknown values are detected to the left of known values in the tax_table, as this should always be wrong/need fixing. Also, by default it now searches a larger list of probably unknown/uninformative tax_table values e.g. "k__NA", "p__Unknown" will now be replaced

## Features
* `comp_barplot` gets merge_other argument, which, when FALSE, shows full diversity of samples by outlining individual taxa within the (grey) "other" category!
* `tax_sort` for sorting taxa in tax_table and otu_table by several name or abundance options (deletes phy_tree if present!)
* `tax_transform` can take a rank argument, to perform aggregation (internally using tax_agg) and transformation with one function, and record the results. **This is now the recommended usage!**
* `tax_transform` new transformation = "binary" can convert to presence/absence data (used by `tax_sort` for by = "prev")
* `tax_top` for flexibly returning top n taxa, at chosen rank, with ordering by `tax_sort`

# microViz 0.5.0 - "hot maps"

## Breaking changes
* `cor_heatmap` and `comp_heatmap` argument changed: 'taxa_which' replaced with 'taxa_side' for easier control over where taxa annotations are placed (default behaviour stays the same)

## Features
* Optionally annotate `cor_heatmap` with variable distributions using `var_anno` and its helpers: `anno_var_box` and `anno_var_hist`
* `cor_heatmap` gets 'var_fun' argument for transforming variables before correlating
* `phyloseq_validate` checks for zero taxa, which can happen after filtering samples
* `tax_filter` gets undetected arg (greater than), as optional alternative to prev_detection_threshold (greater than or equal to)

## Bug fixes
* heatmaps should handle NAs better (`viz_heatmap` internal function fix)
* `heat_palette` can set range arg manually now without errors

## Other
* minor versions from 0.5 will now get a memorable name, probably referring to features added since the last minor version

# microViz 0.4.3

## Features
* `cor_heatmap` for microbe-metadata correlation heatmaps
* `comp_heatmap` for visualising taxonomic composition across samples (ordered/clustered)
* `ord_calc` can now guess which method the user wants by default (by checking for presence of distance matrix and constraints)
* `ord_plot` auto_caption size can now be set, and it also now exposes the `coord_*` args: expand and clip 
* `ord_plot_iris` now handles multiple rings of anno_binary annotations, and anno_binary position is now closer when no anno_colour is set
* `ord_explore` shiny app menu styling is a little cleaner (but still needs some love)
* `tax_scale` for applying `base::scale()` to phyloseq otu_table
* `tax_name` for easily setting informative unique phyloseq taxa_names

## Bug fixes
* `dist_permanova`'s obsolete return arg removed
* `tax_filter` gets explicit handling of compositional input
* `tax_fill_unknowns` gets better handling of fully unclassified taxa
* `taxatree_nodes` now checks for cross-rank name duplications, which would mess up the taxatree graph structure

# microViz 0.4.2

## Bug fixes
* `ord_explore` now works with ps_extra classes
* `tax_fill_unknowns` gives informative error on data with only one taxonomic rank
* `tax_filter` ps argument default removed

# microViz 0.4.1

## Bug fixes:
* `ord_calc` now correctly returns subsetted phyloseq in ps_extra obj after constrained ordination
* `taxatree_plots` and `taxatree_plotkey` now have reactive automatic **minimum** node and edge sizes that depend on the set maxes

# microViz 0.4.0

## Breaking changes
* `dist_permanova` replaces `permanova` for naming consistency and guiding user
* `comp_barplot` replaces `plot_comp_bar` in anticipation of (a) future heatmap function(s) named comp_heatmap or similar

## Features

New **"ps_extra"** class (S3) conveniently stores phyloseq object alongside any calculated distance matrix, ordination and permanova models, as well as records of `tax_agg` and `tax_transform` calls. "ps_extra" class objects have a pretty and compact print method, a simple list structure, and convenient accessor functions to return each component: `ps_get`, `info_get`, `dist_get`, `ord_get`, `perm_get`, `bdisp_get`.

# microViz 0.3.2

## Features
* Taxon modelling updates: 
    * `tax_model` and `taxatree_models` can handle linear modeling e.g. on compositional (TS-Scaled) data
    * `taxatree_plots` has more sensible defaults (automatic symmetrical colour limits and variable selection based on model type)

# microViz 0.3.1

## Breaking changes
* `ord_plot_iris` the `data` arg is replaced with `ord` and conditionally optional `ps` arg for when data in `ord` have been transformed
* `permanova` always uses adonis2 now, so that arg is **removed**, and replaced with `by` argument to set sums of squares choice

## Features
* `ord_plot` gets a `center` argument to center expand the plot limits to center around zero (useful when pairing with `ord_plot_iris`)
* `ord_explore` can now also display ordinations that don't use distances like PCA and RDA (as well as PCoA of course)
* `ord_explore` gains a `ps` arg (for untransformed version) and other tweaks to facilitate using transformed data in `ord`

## Other fixes
* `ord_plot_iris` annotation args are now NULL by default. 
* `microbiome::aggregate_top_taxa` is copied directly into `microViz` to avoid its coming deprecation and the current warning.

# microViz 0.2.0

## Features

* `taxatree_models` and `taxatree_plots` bring tree-based visualisations of statistical models

## Breaking changes

* `tax_model` replaces `tax_model_corncob` as a function that will eventually be generalised to use other model types 
* `taxatree_plots` and `taxatree_plotkey` replaces the rudimentary `taxatree_plot`

# microViz 0.1.2

## Features

* `ord_explore` allows interactive exploration of the composition of ordinated samples using `Shiny`

# microViz 0.1.1

## Features

* `ord_plot_iris` for pca-ordered circular compositional barplots to pair with `ord_plot` output
* `ps_select` for easily selecting variables within phyloseq sample data
* `ord_plot` exposes scaling argument

# microViz 0.1.0

## Major renaming

* `calc_dist` --> `dist_calc`
* `beta_disper` --> `dist_bdisp`
* `ordin8` --> `ord_calc`
* `plot_ordin8` --> `ord_plot`
* `model_tax_corncob` --> `tax_model_corncob`

## Other breaking changes

* `plot_comp_bar` 
    * groups argument renamed to group_by for consistency with new facet_by argument
    * drop_unused_vars replaced with keep_all_vars which defaults to TRUE

## Features

* `ps_seriate` for ordering phyloseq samples by similarity


# microViz 0.0.6

## Features

* `ps_filter` allows filtering of `phyloseq` samples by values of variables in `sample_data`  (wrapper around `dplyr`'s `filter` function)

## Breaking changes

* `ps_mutate` no longer needs `.across` argument to use `dplyr::across()`

# microViz 0.0.5

## Features

* `ps_dedupe` bug fix method = "readcount" now correctly keeps all samples that "first" or "last" methods do

# microViz 0.0.4

## Features

* `ps_arrange` allows reordering of `phyloseq` samples by values of variables in `sample_data` or `otu_table` (wrapper around `dplyr`'s `arrange` function)
* `ps_otu2samdat` add taxon abundance variables from `phyloseq` `otu_table` to `sample_data` for use in plotting etc

# microViz 0.0.3

## Features

* `ps_mutate` allowing easy and piped modification of `phyloseq` `sample_data` (wrapper around `dplyr`'s `mutate` function)
* `ps_join` allows easy and piped addition of a dataframe of data to `phyloseq` `sample_data` (wrapper around `dplyr`'s `*_join` functions)

# microViz 0.0.2

## Breaking changes

* `plot_ordin8` arguments changed to allow easier sizing and styling of all vectors and labels
* `plot_ordin8` default styling of taxon and constraint vectors and labels is changed: background vectors are now semi-transparent and dashed lines are not used any more by default (but can be set)

## Other

* `ordin8` and `plot_ordin8` get basic support for CCA and NMDS methods finally
* `plot_comp_bar` gets `order_samples_with_all_taxa` and `tax_transform_for_ordering`
* Documentation of manual sample ordering across `plot_comp_bar` groups added to Visualising Compositions article on website.
* Documentation of experimental polar coordinates and PCA_angle sorting added in new article called PCA-sorted polar composition plots, on website.

# microViz 0.0.1

## Main changes

* `phyloseq_validate` now fixes otu_tables stored as integers and messages user about suspicious or NA tax_table entries and this happens as part of `plot_comp_bar` and `plot_ordin8`
* `plot_comp_bar` and `plot_ordin8` gain `taxon_renamer` argument to allow you to customise the taxon names on these plots
* `plot_comp_bar` can now handle missing values in the grouping variable by converting NAs to "NA"s

* Some functions renamed for naming consistency: 
    - `prepend_ranks` -> `tax_prepend_ranks` 
    - `tax_model_corncob` -> `model_tax_corncob`
    - `corncob_models_to_var_stats` -> `modelsmodels2stats_corncob`
    - `tax_tree_` `nodes`/`edges`/`plot` -> `taxatree_` ...

## Other

* Added a `NEWS.md` file to track changes to the package.

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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# phyloseq_validate warns about removing all zero taxa

    Code
      phyloseq_validate(ps = dietswap, remove_undetected = TRUE, verbose = TRUE)
    Warning <simpleWarning>
      Some taxa_sums were zero, removing the following taxa:
      	Aerococcus 
      	Aneurinibacillus 
      	Asteroleplasma et rel. 
      	Clostridium felsineum et rel. 
      	Clostridium thermocellum et rel. 
      	Methylobacterium 
      	Micrococcaceae
      This may be caused by using `subset_samples()`.
      Try using `ps_filter()` instead, with .keep_all_taxa = FALSE.
      Otherwise, to avoid this warning, try filtering out taxa summing to zero with `tax_filter()`.
      If you have already transformed and/or scaled your taxa, e.g. with a log transformation or scale,
      seeing this warning is possible, but very unlikely and possibly a bug. Please report this.
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 123 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 123 taxa by 3 taxonomic ranks ]

# phyloseq_validate fixes missing sam_data with message

    Code
      phyloseq_validate(ps = dietswap, verbose = TRUE)
    Message <simpleMessage>
      Note: Replacing missing sample_data with a dataframe of only sample_names.
      Try `ps <- phyloseq_validate(ps, verbose = FALSE)` to avoid this message
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 1 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

# phyloseq_validate fixes missing tax_table with message

    Code
      phyloseq_validate(soilrep, verbose = TRUE)
    Message <simpleMessage>
      Note: Replacing missing tax_table with a 1-column table of only taxa_names.
      Try `ps <- phyloseq_validate(ps, verbose = FALSE)` to avoid this message
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 16825 taxa and 56 samples ]
      sample_data() Sample Data:       [ 56 samples by 4 sample variables ]
      tax_table()   Taxonomy Table:    [ 16825 taxa by 1 taxonomic ranks ]

# expected names order for comparison hasn't changed

    Code
      cat(df[["names"]])
    Output
      Prevotella melaninogenica et rel. Oscillospira guillermondii et rel. Bacteroides vulgatus et rel. Clostridium cellulosi et rel. Prevotella oralis et rel. Faecalibacterium prausnitzii et rel. Sporobacter termitidis et rel. Clostridium symbiosum et rel. Allistipes et rel. Clostridium orbiscindens et rel. Subdoligranulum variable at rel. Ruminococcus obeum et rel. Butyrivibrio crossotus et rel. Bacteroides fragilis et rel. Akkermansia Bacteroides ovatus et rel. Parabacteroides distasonis et rel. Dorea formicigenerans et rel. Bacteroides uniformis et rel. Dialister Bryantella formatexigens et rel. Uncultured Clostridiales I Coprococcus eutactus et rel. Clostridium leptum et rel. Clostridium sphenoides et rel. Escherichia coli et rel. Streptococcus bovis et rel. Uncultured Clostridiales II Bifidobacterium Anaerotruncus colihominis et rel. Lachnospira pectinoschiza et rel. Anaerostipes caccae et rel. Ruminococcus callidus et rel. Bacteroides splachnicus et rel. Ruminococcus bromii et rel. Prevotella tannerae et rel. Lachnobacillus bovis et rel. Eubacterium rectale et rel. Mitsuokella multiacida et rel. Outgrouping clostridium cluster XIVa Clostridium nexile et rel. Uncultured Mollicutes Streptococcus mitis et rel. Megasphaera elsdenii et rel. Bacteroides plebeius et rel. Sutterella wadsworthia et rel. Oxalobacter formigenes et rel. Tannerella et rel. Phascolarctobacterium faecium et rel. Papillibacter cinnamivorans et rel. Ruminococcus gnavus et rel. Clostridium colinum et rel. Clostridium difficile et rel. Eubacterium biforme et rel. Bacteroides stercoris et rel. Clostridium (sensu stricto) Roseburia intestinalis et rel. Eubacterium hallii et rel. Enterobacter aerogenes et rel. Eubacterium ventriosum et rel. Lactobacillus plantarum et rel. Anaerovorax odorimutans et rel. Klebisiella pneumoniae et rel. Collinsella Enterococcus Ruminococcus lactaris et rel. Clostridium stercorarium et rel. Veillonella Lactobacillus gasseri et rel. Alcaligenes faecalis et rel. Streptococcus intermedius et rel. Fusobacteria Proteus et rel. Eggerthella lenta et rel. Desulfovibrio et rel. Bulleidia moorei et rel. Campylobacter Peptococcus niger et rel. Eubacterium siraeum et rel. Weissella et rel. Coprobacillus catenaformis et rel. Oceanospirillum Vibrio Yersinia et rel. Eubacterium cylindroides et rel. Lactobacillus salivarius et rel. Clostridium ramosum et rel. Peptostreptococcus micros et rel. Helicobacter Lactococcus Bilophila et rel. Eubacterium limosum et rel. Propionibacterium Megamonas hypermegale et rel. Bacillus Corynebacterium Prevotella ruminicola et rel. Bacteroides intestinalis et rel. Uncultured Bacteroidetes Xanthomonadaceae Catenibacterium mitsuokai et rel. Burkholderia Haemophilus Serratia Leminorella Moraxellaceae Uncultured Selenomonadaceae Lactobacillus catenaformis et rel. Atopobium Aquabacterium Brachyspira Granulicatella Actinomycetaceae Staphylococcus Wissella et rel. Pseudomonas Anaerofustis Uncultured Chroococcales Gemella Aeromonas Anaerobiospirillum Novosphingobium Peptostreptococcus anaerobius et rel. Aerococcus Aneurinibacillus Asteroleplasma et rel. Clostridium felsineum et rel. Clostridium thermocellum et rel. Methylobacterium Micrococcaceae

# cor_heatmap doesn't change: 

    Code
      p@matrix_param[names(p@matrix_param) != "cell_fun"]
    Output
      $row_km
      [1] 1
      
      $row_km_repeats
      [1] 1
      
      $row_gap
      [1] 1mm
      
      $column_km
      [1] 1
      
      $column_km_repeats
      [1] 1
      
      $column_gap
      [1] 1mm
      
      $jitter
      [1] FALSE
      
      $gp
      $col
      [1] "white"
      
      $lwd
      [1] 0.5
      
      $lineheight
      [1] 0.9
      
      
      $border
      [1] NA
      
      $border_gp
      $col
      [1] "black"
      
      
      $width
      [1] 5null
      
      $height
      [1] 30null
      

---

    Code
      p@matrix_color_mapping
    Output
      Continuous color mapping:
      name: pearson 
      default breaks:
      [1] -1.0 -0.5  0.0  0.5  1.0
      
      colors:
      [1] "#4A6FE3FF" "#6B82E1FF" "#E2E2E2FF" "#DA5F7DFF" "#D33F6AFF"
      

---

    Code
      p@right_annotation@anno_list
    Output
      $prev
      A single annotation with anno_barplot() function
        name: prev 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      
      $abund
      A single annotation with anno_boxplot() function
        name: abund 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      

---

    Code
      str(p@column_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 5 members at h = 2.23]
        |--leaf "african" 
        `--[dendrogram w/ 2 branches and 4 members at h = 0.946]
           |--[dendrogram w/ 2 branches and 2 members at h = 0.415]
           |  |--leaf "timepoint" 
           |  `--leaf "timepoint.within.group" 
           `--[dendrogram w/ 2 branches and 2 members at h = 0.738]
              |--leaf "weight" 
              `--leaf "female" 

---

    Code
      str(p@row_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 30 members at h = 1.62]
        |--[dendrogram w/ 2 branches and 6 members at h = 0.301]
        |  |--[dendrogram w/ 2 branches and 2 members at h = 0.0866]
        |  |  |--leaf "Allistipes et rel." 
        |  |  `--leaf "Bacteroides vulgatus et rel." 
        |  `--[dendrogram w/ 2 branches and 4 members at h = 0.223]
        |     |--[dendrogram w/ 2 branches and 3 members at h = 0.139]
        |     |  |--[dendrogram w/ 2 branches and 2 members at h = 0.119]
        |     |  |  |--leaf "Bacteroides plebeius et rel." 
        |     |  |  `--leaf "Tannerella et rel." 
        |     |  `--leaf "Parabacteroides distasonis et rel." 
        |     `--leaf "Bacteroides uniformis et rel." 
        `--[dendrogram w/ 2 branches and 24 members at h = 1.11]
           |--[dendrogram w/ 2 branches and 17 members at h = 0.716]
           |  |--[dendrogram w/ 2 branches and 11 members at h = 0.547]
           |  |  |--[dendrogram w/ 2 branches and 7 members at h = 0.39]
           |  |  |  |--[dendrogram w/ 2 branches and 3 members at h = 0.309]
           |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 0.165]
           |  |  |  |  |  |--leaf "Prevotella tannerae et rel." 
           |  |  |  |  |  `--leaf "Anaerostipes caccae et rel." 
           |  |  |  |  `--leaf "Bryantella formatexigens et rel." 
           |  |  |  `--[dendrogram w/ 2 branches and 4 members at h = 0.213]
           |  |  |     |--[dendrogram w/ 2 branches and 2 members at h = 0.0816]
           |  |  |     |  |--leaf "Subdoligranulum variable at rel." 
           |  |  |     |  `--leaf "Bifidobacterium" 
           |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 0.132]
           |  |  |        |--leaf "Clostridium symbiosum et rel." 
           |  |  |        `--leaf "Clostridium sphenoides et rel." 
           |  |  `--[dendrogram w/ 2 branches and 4 members at h = 0.352]
           |  |     |--leaf "Bacteroides fragilis et rel." 
           |  |     `--[dendrogram w/ 2 branches and 3 members at h = 0.249]
           |  |        |--[dendrogram w/ 2 branches and 2 members at h = 0.214]
           |  |        |  |--leaf "Lachnospira pectinoschiza et rel." 
           |  |        |  `--leaf "Lachnobacillus bovis et rel." 
           |  |        `--leaf "Akkermansia" 
           |  `--[dendrogram w/ 2 branches and 6 members at h = 0.295]
           |     |--[dendrogram w/ 2 branches and 2 members at h = 0.119]
           |     |  |--leaf "Escherichia coli et rel." 
           |     |  `--leaf "Clostridium cellulosi et rel." 
           |     `--[dendrogram w/ 2 branches and 4 members at h = 0.23]
           |        |--[dendrogram w/ 2 branches and 3 members at h = 0.119]
           |        |  |--leaf "Phascolarctobacterium faecium et rel." 
           |        |  `--[dendrogram w/ 2 branches and 2 members at h = 0.0523]
           |        |     |--leaf "Streptococcus mitis et rel." 
           |        |     `--leaf "Streptococcus bovis et rel." 
           |        `--leaf "Ruminococcus obeum et rel." 
           `--[dendrogram w/ 2 branches and 7 members at h = 0.422]
              |--[dendrogram w/ 2 branches and 6 members at h = 0.402]
              |  |--[dendrogram w/ 2 branches and 4 members at h = 0.278]
              |  |  |--[dendrogram w/ 2 branches and 3 members at h = 0.228]
              |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 0.109]
              |  |  |  |  |--leaf "Clostridium nexile et rel." 
              |  |  |  |  `--leaf "Outgrouping clostridium cluster XIVa" 
              |  |  |  `--leaf "Clostridium orbiscindens et rel." 
              |  |  `--leaf "Uncultured Mollicutes" 
              |  `--[dendrogram w/ 2 branches and 2 members at h = 0.221]
              |     |--leaf "Sporobacter termitidis et rel." 
              |     `--leaf "Uncultured Clostridiales II" 
              `--leaf "Prevotella oralis et rel." 

# cor_heatmap with var_anno doesn't change: 

    Code
      v@matrix_param[names(v@matrix_param) != "cell_fun"]
    Output
      $row_km
      [1] 1
      
      $row_km_repeats
      [1] 1
      
      $row_gap
      [1] 1mm
      
      $column_km
      [1] 1
      
      $column_km_repeats
      [1] 1
      
      $column_gap
      [1] 1mm
      
      $jitter
      [1] FALSE
      
      $gp
      $col
      [1] "white"
      
      $lwd
      [1] 0.5
      
      $lineheight
      [1] 0.9
      
      
      $border
      [1] NA
      
      $border_gp
      $col
      [1] "black"
      
      
      $width
      [1] 5null
      
      $height
      [1] 30null
      

---

    Code
      v@matrix_color_mapping
    Output
      Continuous color mapping:
      name: pearson 
      default breaks:
      [1] -1.0 -0.5  0.0  0.5  1.0
      
      colors:
      [1] "#4A6FE3FF" "#6B82E1FF" "#E2E2E2FF" "#DA5F7DFF" "#D33F6AFF"
      

---

    Code
      v@right_annotation@anno_list
    Output
      $prev
      A single annotation with anno_barplot() function
        name: prev 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      
      $abund
      A single annotation with anno_boxplot() function
        name: abund 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      

---

    Code
      v@top_annotation@anno_list
    Output
      $x
      A single annotation with anno_histogram() function
        name: x 
        position: column 
        no legend
        items: 5 
        width: 1npc 
        height: 10mm 
        this object is  subsetable
        3.56915555555556mm extension on the left 
        2.56915555555556mm extension on the right 
      
      $`log10(x+1)`
      A single annotation with anno_boxplot() function
        name: log10(x+1) 
        position: column 
        no legend
        items: 5 
        width: 1npc 
        height: 20mm 
        this object is  subsetable
        5.92288888888889mm extension on the left 
        15.0377333333333mm extension on the right 
      

---

    Code
      str(v@column_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 5 members at h = 2.23]
        |--leaf "african" 
        `--[dendrogram w/ 2 branches and 4 members at h = 0.946]
           |--[dendrogram w/ 2 branches and 2 members at h = 0.415]
           |  |--leaf "timepoint" 
           |  `--leaf "timepoint.within.group" 
           `--[dendrogram w/ 2 branches and 2 members at h = 0.738]
              |--leaf "weight" 
              `--leaf "female" 

---

    Code
      str(v@row_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 30 members at h = 1.62]
        |--[dendrogram w/ 2 branches and 6 members at h = 0.301]
        |  |--[dendrogram w/ 2 branches and 2 members at h = 0.0866]
        |  |  |--leaf "Allistipes et rel." 
        |  |  `--leaf "Bacteroides vulgatus et rel." 
        |  `--[dendrogram w/ 2 branches and 4 members at h = 0.223]
        |     |--[dendrogram w/ 2 branches and 3 members at h = 0.139]
        |     |  |--[dendrogram w/ 2 branches and 2 members at h = 0.119]
        |     |  |  |--leaf "Bacteroides plebeius et rel." 
        |     |  |  `--leaf "Tannerella et rel." 
        |     |  `--leaf "Parabacteroides distasonis et rel." 
        |     `--leaf "Bacteroides uniformis et rel." 
        `--[dendrogram w/ 2 branches and 24 members at h = 1.11]
           |--[dendrogram w/ 2 branches and 17 members at h = 0.716]
           |  |--[dendrogram w/ 2 branches and 11 members at h = 0.547]
           |  |  |--[dendrogram w/ 2 branches and 7 members at h = 0.39]
           |  |  |  |--[dendrogram w/ 2 branches and 3 members at h = 0.309]
           |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 0.165]
           |  |  |  |  |  |--leaf "Prevotella tannerae et rel." 
           |  |  |  |  |  `--leaf "Anaerostipes caccae et rel." 
           |  |  |  |  `--leaf "Bryantella formatexigens et rel." 
           |  |  |  `--[dendrogram w/ 2 branches and 4 members at h = 0.213]
           |  |  |     |--[dendrogram w/ 2 branches and 2 members at h = 0.0816]
           |  |  |     |  |--leaf "Subdoligranulum variable at rel." 
           |  |  |     |  `--leaf "Bifidobacterium" 
           |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 0.132]
           |  |  |        |--leaf "Clostridium symbiosum et rel." 
           |  |  |        `--leaf "Clostridium sphenoides et rel." 
           |  |  `--[dendrogram w/ 2 branches and 4 members at h = 0.352]
           |  |     |--leaf "Bacteroides fragilis et rel." 
           |  |     `--[dendrogram w/ 2 branches and 3 members at h = 0.249]
           |  |        |--[dendrogram w/ 2 branches and 2 members at h = 0.214]
           |  |        |  |--leaf "Lachnospira pectinoschiza et rel." 
           |  |        |  `--leaf "Lachnobacillus bovis et rel." 
           |  |        `--leaf "Akkermansia" 
           |  `--[dendrogram w/ 2 branches and 6 members at h = 0.295]
           |     |--[dendrogram w/ 2 branches and 2 members at h = 0.119]
           |     |  |--leaf "Escherichia coli et rel." 
           |     |  `--leaf "Clostridium cellulosi et rel." 
           |     `--[dendrogram w/ 2 branches and 4 members at h = 0.23]
           |        |--[dendrogram w/ 2 branches and 3 members at h = 0.119]
           |        |  |--leaf "Phascolarctobacterium faecium et rel." 
           |        |  `--[dendrogram w/ 2 branches and 2 members at h = 0.0523]
           |        |     |--leaf "Streptococcus mitis et rel." 
           |        |     `--leaf "Streptococcus bovis et rel." 
           |        `--leaf "Ruminococcus obeum et rel." 
           `--[dendrogram w/ 2 branches and 7 members at h = 0.422]
              |--[dendrogram w/ 2 branches and 6 members at h = 0.402]
              |  |--[dendrogram w/ 2 branches and 4 members at h = 0.278]
              |  |  |--[dendrogram w/ 2 branches and 3 members at h = 0.228]
              |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 0.109]
              |  |  |  |  |--leaf "Clostridium nexile et rel." 
              |  |  |  |  `--leaf "Outgrouping clostridium cluster XIVa" 
              |  |  |  `--leaf "Clostridium orbiscindens et rel." 
              |  |  `--leaf "Uncultured Mollicutes" 
              |  `--[dendrogram w/ 2 branches and 2 members at h = 0.221]
              |     |--leaf "Sporobacter termitidis et rel." 
              |     `--leaf "Uncultured Clostridiales II" 
              `--leaf "Prevotella oralis et rel." 

# comp_heatmap doesn't change: 

    Code
      p@matrix_param
    Output
      $row_km
      [1] 1
      
      $row_km_repeats
      [1] 1
      
      $row_gap
      [1] 1mm
      
      $column_km
      [1] 1
      
      $column_km_repeats
      [1] 1
      
      $column_gap
      [1] 1mm
      
      $jitter
      [1] FALSE
      
      $gp
      $col
      [1] "white"
      
      $lwd
      [1] 0.5
      
      $lineheight
      [1] 0.9
      
      
      $border
      [1] NA
      
      $border_gp
      $col
      [1] "black"
      
      
      $width
      [1] 222null
      
      $height
      [1] 30null
      

---

    Code
      p@matrix_color_mapping
    Output
      Continuous color mapping:
      name: Abd. 
      default breaks:
      [1] -2  0  2  4  6  8
      
      colors:
      [1] "#FDF5EBFF" "#EF9B66FF" "#DE2F52FF" "#910062FF" "#40043FFF" "#070707FF"
      

---

    Code
      p@right_annotation@anno_list
    Output
      $prev
      A single annotation with anno_barplot() function
        name: prev 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      
      $abund
      A single annotation with anno_boxplot() function
        name: abund 
        position: row 
        no legend
        items: 30 
        width: 15mm 
        height: 1npc 
        this object is  subsetable
        8.94924444444444mm extension on the bottom 
      

---

    Code
      str(p@column_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 222 members at h = 61]
        |--[dendrogram w/ 2 branches and 94 members at h = 22.7]
        |  |--[dendrogram w/ 2 branches and 50 members at h = 17.1]
        |  |  |--[dendrogram w/ 2 branches and 32 members at h = 11.8]
        |  |  |  |--[dendrogram w/ 2 branches and 14 members at h = 9.27]
        |  |  |  |  |--[dendrogram w/ 2 branches and 5 members at h = 6.76]
        |  |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 5.37]
        |  |  |  |  |  |  |--leaf "Sample-35" 
        |  |  |  |  |  |  `--leaf "Sample-132" 
        |  |  |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 4.54]
        |  |  |  |  |     |--leaf "Sample-112" 
        |  |  |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 3.12]
        |  |  |  |  |        |--leaf "Sample-118" 
        |  |  |  |  |        `--leaf "Sample-69" 
        |  |  |  |  `--[dendrogram w/ 2 branches and 9 members at h = 7.26]
        |  |  |  |     |--[dendrogram w/ 2 branches and 6 members at h = 4.96]
        |  |  |  |     |  |--[dendrogram w/ 2 branches and 3 members at h = 3.18]
        |  |  |  |     |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.72]
        |  |  |  |     |  |  |  |--leaf "Sample-220" 
        |  |  |  |     |  |  |  `--leaf "Sample-31" 
        |  |  |  |     |  |  `--leaf "Sample-15" 
        |  |  |  |     |  `--[dendrogram w/ 2 branches and 3 members at h = 4.02]
        |  |  |  |     |     |--[dendrogram w/ 2 branches and 2 members at h = 2.98]
        |  |  |  |     |     |  |--leaf "Sample-5" 
        |  |  |  |     |     |  `--leaf "Sample-208" 
        |  |  |  |     |     `--leaf "Sample-54" 
        |  |  |  |     `--[dendrogram w/ 2 branches and 3 members at h = 3.59]
        |  |  |  |        |--[dendrogram w/ 2 branches and 2 members at h = 1.92]
        |  |  |  |        |  |--leaf "Sample-53" 
        |  |  |  |        |  `--leaf "Sample-102" 
        |  |  |  |        `--leaf "Sample-79" 
        |  |  |  `--[dendrogram w/ 2 branches and 18 members at h = 7.47]
        |  |  |     |--[dendrogram w/ 2 branches and 11 members at h = 6.28]
        |  |  |     |  |--[dendrogram w/ 2 branches and 2 members at h = 2]
        |  |  |     |  |  |--leaf "Sample-193" 
        |  |  |     |  |  `--leaf "Sample-190" 
        |  |  |     |  `--[dendrogram w/ 2 branches and 9 members at h = 5.39]
        |  |  |     |     |--[dendrogram w/ 2 branches and 7 members at h = 4.69]
        |  |  |     |     |  |--[dendrogram w/ 2 branches and 5 members at h = 3.79]
        |  |  |     |     |  |  |--[dendrogram w/ 2 branches and 3 members at h = 3.27]
        |  |  |     |     |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.77]
        |  |  |     |     |  |  |  |  |--leaf "Sample-187" 
        |  |  |     |     |  |  |  |  `--leaf "Sample-107" 
        |  |  |     |     |  |  |  `--leaf "Sample-94" 
        |  |  |     |     |  |  `--[dendrogram w/ 2 branches and 2 members at h = 2.53]
        |  |  |     |     |  |     |--leaf "Sample-168" 
        |  |  |     |     |  |     `--leaf "Sample-182" 
        |  |  |     |     |  `--[dendrogram w/ 2 branches and 2 members at h = 3.19]
        |  |  |     |     |     |--leaf "Sample-161" 
        |  |  |     |     |     `--leaf "Sample-212" 
        |  |  |     |     `--[dendrogram w/ 2 branches and 2 members at h = 3.22]
        |  |  |     |        |--leaf "Sample-207" 
        |  |  |     |        `--leaf "Sample-61" 
        |  |  |     `--[dendrogram w/ 2 branches and 7 members at h = 6.61]
        |  |  |        |--[dendrogram w/ 2 branches and 3 members at h = 3.89]
        |  |  |        |  |--[dendrogram w/ 2 branches and 2 members at h = 2.87]
        |  |  |        |  |  |--leaf "Sample-209" 
        |  |  |        |  |  `--leaf "Sample-210" 
        |  |  |        |  `--leaf "Sample-213" 
        |  |  |        `--[dendrogram w/ 2 branches and 4 members at h = 4.65]
        |  |  |           |--[dendrogram w/ 2 branches and 3 members at h = 3.33]
        |  |  |           |  |--leaf "Sample-50" 
        |  |  |           |  `--[dendrogram w/ 2 branches and 2 members at h = 1.91]
        |  |  |           |     |--leaf "Sample-89" 
        |  |  |           |     `--leaf "Sample-104" 
        |  |  |           `--leaf "Sample-56" 
        |  |  `--[dendrogram w/ 2 branches and 18 members at h = 11.5]
        |  |     |--[dendrogram w/ 2 branches and 8 members at h = 6.82]
        |  |     |  |--[dendrogram w/ 2 branches and 5 members at h = 4.69]
        |  |     |  |  |--[dendrogram w/ 2 branches and 2 members at h = 4.2]
        |  |     |  |  |  |--leaf "Sample-55" 
        |  |     |  |  |  `--leaf "Sample-95" 
        |  |     |  |  `--[dendrogram w/ 2 branches and 3 members at h = 3.34]
        |  |     |  |     |--leaf "Sample-81" 
        |  |     |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.83]
        |  |     |  |        |--leaf "Sample-88" 
        |  |     |  |        `--leaf "Sample-101" 
        |  |     |  `--[dendrogram w/ 2 branches and 3 members at h = 4.18]
        |  |     |     |--[dendrogram w/ 2 branches and 2 members at h = 2.6]
        |  |     |     |  |--leaf "Sample-60" 
        |  |     |     |  `--leaf "Sample-77" 
        |  |     |     `--leaf "Sample-66" 
        |  |     `--[dendrogram w/ 2 branches and 10 members at h = 8.16]
        |  |        |--[dendrogram w/ 2 branches and 8 members at h = 7.88]
        |  |        |  |--[dendrogram w/ 2 branches and 5 members at h = 5.95]
        |  |        |  |  |--[dendrogram w/ 2 branches and 2 members at h = 3.64]
        |  |        |  |  |  |--leaf "Sample-67" 
        |  |        |  |  |  `--leaf "Sample-86" 
        |  |        |  |  `--[dendrogram w/ 2 branches and 3 members at h = 4.61]
        |  |        |  |     |--[dendrogram w/ 2 branches and 2 members at h = 4.26]
        |  |        |  |     |  |--leaf "Sample-6" 
        |  |        |  |     |  `--leaf "Sample-59" 
        |  |        |  |     `--leaf "Sample-93" 
        |  |        |  `--[dendrogram w/ 2 branches and 3 members at h = 4.85]
        |  |        |     |--leaf "Sample-52" 
        |  |        |     `--[dendrogram w/ 2 branches and 2 members at h = 2.52]
        |  |        |        |--leaf "Sample-103" 
        |  |        |        `--leaf "Sample-99" 
        |  |        `--[dendrogram w/ 2 branches and 2 members at h = 3.7]
        |  |           |--leaf "Sample-105" 
        |  |           `--leaf "Sample-9" 
        |  `--[dendrogram w/ 2 branches and 44 members at h = 18.2]
        |     |--[dendrogram w/ 2 branches and 10 members at h = 13.9]
        |     |  |--[dendrogram w/ 2 branches and 5 members at h = 8.2]
        |     |  |  |--[dendrogram w/ 2 branches and 3 members at h = 5.28]
        |     |  |  |  |--leaf "Sample-14" 
        |     |  |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.57]
        |     |  |  |     |--leaf "Sample-30" 
        |     |  |  |     `--leaf "Sample-22" 
        |     |  |  `--[dendrogram w/ 2 branches and 2 members at h = 4.73]
        |     |  |     |--leaf "Sample-4" 
        |     |  |     `--leaf "Sample-40" 
        |     |  `--[dendrogram w/ 2 branches and 5 members at h = 6.67]
        |     |     |--[dendrogram w/ 2 branches and 2 members at h = 2.06]
        |     |     |  |--leaf "Sample-121" 
        |     |     |  `--leaf "Sample-126" 
        |     |     `--[dendrogram w/ 2 branches and 3 members at h = 4.37]
        |     |        |--[dendrogram w/ 2 branches and 2 members at h = 2.83]
        |     |        |  |--leaf "Sample-137" 
        |     |        |  `--leaf "Sample-147" 
        |     |        `--leaf "Sample-49" 
        |     `--[dendrogram w/ 2 branches and 34 members at h = 15.4]
        |        |--[dendrogram w/ 2 branches and 19 members at h = 8.97]
        |        |  |--[dendrogram w/ 2 branches and 6 members at h = 6.78]
        |        |  |  |--[dendrogram w/ 2 branches and 4 members at h = 4.22]
        |        |  |  |  |--leaf "Sample-16" 
        |        |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 2.56]
        |        |  |  |     |--leaf "Sample-218" 
        |        |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.23]
        |        |  |  |        |--leaf "Sample-29" 
        |        |  |  |        `--leaf "Sample-21" 
        |        |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.84]
        |        |  |     |--leaf "Sample-7" 
        |        |  |     `--leaf "Sample-2" 
        |        |  `--[dendrogram w/ 2 branches and 13 members at h = 8.36]
        |        |     |--[dendrogram w/ 2 branches and 7 members at h = 6.49]
        |        |     |  |--[dendrogram w/ 2 branches and 4 members at h = 5.35]
        |        |     |  |  |--[dendrogram w/ 2 branches and 2 members at h = 3.82]
        |        |     |  |  |  |--leaf "Sample-25" 
        |        |     |  |  |  `--leaf "Sample-222" 
        |        |     |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.61]
        |        |     |  |     |--leaf "Sample-3" 
        |        |     |  |     `--leaf "Sample-214" 
        |        |     |  `--[dendrogram w/ 2 branches and 3 members at h = 3.87]
        |        |     |     |--[dendrogram w/ 2 branches and 2 members at h = 2.85]
        |        |     |     |  |--leaf "Sample-8" 
        |        |     |     |  `--leaf "Sample-20" 
        |        |     |     `--leaf "Sample-23" 
        |        |     `--[dendrogram w/ 2 branches and 6 members at h = 5.93]
        |        |        |--[dendrogram w/ 2 branches and 4 members at h = 4.12]
        |        |        |  |--[dendrogram w/ 2 branches and 3 members at h = 3.08]
        |        |        |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.7]
        |        |        |  |  |  |--leaf "Sample-217" 
        |        |        |  |  |  `--leaf "Sample-12" 
        |        |        |  |  `--leaf "Sample-28" 
        |        |        |  `--leaf "Sample-13" 
        |        |        `--[dendrogram w/ 2 branches and 2 members at h = 4.71]
        |        |           |--leaf "Sample-11" 
        |        |           `--leaf "Sample-10" 
        |        `--[dendrogram w/ 2 branches and 15 members at h = 9.37]
        |           |--[dendrogram w/ 2 branches and 11 members at h = 7.06]
        |           |  |--[dendrogram w/ 2 branches and 6 members at h = 5.67]
        |           |  |  |--[dendrogram w/ 2 branches and 2 members at h = 4.95]
        |           |  |  |  |--leaf "Sample-76" 
        |           |  |  |  `--leaf "Sample-83" 
        |           |  |  `--[dendrogram w/ 2 branches and 4 members at h = 3.53]
        |           |  |     |--[dendrogram w/ 2 branches and 2 members at h = 3.15]
        |           |  |     |  |--leaf "Sample-63" 
        |           |  |     |  `--leaf "Sample-64" 
        |           |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.66]
        |           |  |        |--leaf "Sample-90" 
        |           |  |        `--leaf "Sample-57" 
        |           |  `--[dendrogram w/ 2 branches and 5 members at h = 4.2]
        |           |     |--[dendrogram w/ 2 branches and 2 members at h = 2.3]
        |           |     |  |--leaf "Sample-58" 
        |           |     |  `--leaf "Sample-91" 
        |           |     `--[dendrogram w/ 2 branches and 3 members at h = 3.69]
        |           |        |--leaf "Sample-24" 
        |           |        `--[dendrogram w/ 2 branches and 2 members at h = 3.22]
        |           |           |--leaf "Sample-65" 
        |           |           `--leaf "Sample-85" 
        |           `--[dendrogram w/ 2 branches and 4 members at h = 6.22]
        |              |--[dendrogram w/ 2 branches and 2 members at h = 2.31]
        |              |  |--leaf "Sample-92" 
        |              |  `--leaf "Sample-98" 
        |              `--[dendrogram w/ 2 branches and 2 members at h = 3.18]
        |                 |--leaf "Sample-96" 
        |                 `--leaf "Sample-97" 
        `--[dendrogram w/ 2 branches and 128 members at h = 33.3]
           |--[dendrogram w/ 2 branches and 63 members at h = 24]
           |  |--[dendrogram w/ 2 branches and 38 members at h = 17.1]
           |  |  |--[dendrogram w/ 2 branches and 9 members at h = 7.82]
           |  |  |  |--[dendrogram w/ 2 branches and 6 members at h = 5.62]
           |  |  |  |  |--[dendrogram w/ 2 branches and 3 members at h = 3.55]
           |  |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.17]
           |  |  |  |  |  |  |--leaf "Sample-75" 
           |  |  |  |  |  |  `--leaf "Sample-84" 
           |  |  |  |  |  `--leaf "Sample-100" 
           |  |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 2]
           |  |  |  |     |--leaf "Sample-18" 
           |  |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 1.63]
           |  |  |  |        |--leaf "Sample-26" 
           |  |  |  |        `--leaf "Sample-215" 
           |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 4.21]
           |  |  |     |--leaf "Sample-221" 
           |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 1.86]
           |  |  |        |--leaf "Sample-216" 
           |  |  |        `--leaf "Sample-27" 
           |  |  `--[dendrogram w/ 2 branches and 29 members at h = 13.9]
           |  |     |--[dendrogram w/ 2 branches and 11 members at h = 8.61]
           |  |     |  |--[dendrogram w/ 2 branches and 7 members at h = 7.38]
           |  |     |  |  |--[dendrogram w/ 2 branches and 4 members at h = 5.27]
           |  |     |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 3.15]
           |  |     |  |  |  |  |--leaf "Sample-80" 
           |  |     |  |  |  |  `--leaf "Sample-87" 
           |  |     |  |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.87]
           |  |     |  |  |     |--leaf "Sample-51" 
           |  |     |  |  |     `--leaf "Sample-108" 
           |  |     |  |  `--[dendrogram w/ 2 branches and 3 members at h = 3.48]
           |  |     |  |     |--leaf "Sample-133" 
           |  |     |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.55]
           |  |     |  |        |--leaf "Sample-143" 
           |  |     |  |        `--leaf "Sample-44" 
           |  |     |  `--[dendrogram w/ 2 branches and 4 members at h = 4.56]
           |  |     |     |--[dendrogram w/ 2 branches and 2 members at h = 4.07]
           |  |     |     |  |--leaf "Sample-34" 
           |  |     |     |  `--leaf "Sample-111" 
           |  |     |     `--[dendrogram w/ 2 branches and 2 members at h = 2.99]
           |  |     |        |--leaf "Sample-142" 
           |  |     |        `--leaf "Sample-43" 
           |  |     `--[dendrogram w/ 2 branches and 18 members at h = 10.5]
           |  |        |--[dendrogram w/ 2 branches and 8 members at h = 8.35]
           |  |        |  |--[dendrogram w/ 2 branches and 3 members at h = 4.4]
           |  |        |  |  |--leaf "Sample-149" 
           |  |        |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.29]
           |  |        |  |     |--leaf "Sample-129" 
           |  |        |  |     `--leaf "Sample-123" 
           |  |        |  `--[dendrogram w/ 2 branches and 5 members at h = 6.68]
           |  |        |     |--[dendrogram w/ 2 branches and 3 members at h = 5.33]
           |  |        |     |  |--[dendrogram w/ 2 branches and 2 members at h = 4.23]
           |  |        |     |  |  |--leaf "Sample-163" 
           |  |        |     |  |  `--leaf "Sample-17" 
           |  |        |     |  `--leaf "Sample-165" 
           |  |        |     `--[dendrogram w/ 2 branches and 2 members at h = 4.04]
           |  |        |        |--leaf "Sample-115" 
           |  |        |        `--leaf "Sample-136" 
           |  |        `--[dendrogram w/ 2 branches and 10 members at h = 8.68]
           |  |           |--[dendrogram w/ 2 branches and 4 members at h = 5.97]
           |  |           |  |--[dendrogram w/ 2 branches and 2 members at h = 4.39]
           |  |           |  |  |--leaf "Sample-146" 
           |  |           |  |  `--leaf "Sample-162" 
           |  |           |  `--[dendrogram w/ 2 branches and 2 members at h = 5.26]
           |  |           |     |--leaf "Sample-19" 
           |  |           |     `--leaf "Sample-156" 
           |  |           `--[dendrogram w/ 2 branches and 6 members at h = 6.07]
           |  |              |--[dendrogram w/ 2 branches and 2 members at h = 3.86]
           |  |              |  |--leaf "Sample-211" 
           |  |              |  `--leaf "Sample-219" 
           |  |              `--[dendrogram w/ 2 branches and 4 members at h = 4.63]
           |  |                 |--leaf "Sample-169" 
           |  |                 `--[dendrogram w/ 2 branches and 3 members at h = 2.77]
           |  |                    |--leaf "Sample-174" 
           |  |                    `--[dendrogram w/ 2 branches and 2 members at h = 2.53]
           |  |                       |--leaf "Sample-178" 
           |  |                       `--leaf "Sample-184" 
           |  `--[dendrogram w/ 2 branches and 25 members at h = 21.3]
           |     |--[dendrogram w/ 2 branches and 17 members at h = 11.5]
           |     |  |--[dendrogram w/ 2 branches and 8 members at h = 9.71]
           |     |  |  |--[dendrogram w/ 2 branches and 5 members at h = 6.95]
           |     |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 4.13]
           |     |  |  |  |  |--leaf "Sample-78" 
           |     |  |  |  |  `--leaf "Sample-106" 
           |     |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 4.16]
           |     |  |  |     |--leaf "Sample-62" 
           |     |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 3.62]
           |     |  |  |        |--leaf "Sample-82" 
           |     |  |  |        `--leaf "Sample-74" 
           |     |  |  `--[dendrogram w/ 2 branches and 3 members at h = 5.25]
           |     |  |     |--[dendrogram w/ 2 branches and 2 members at h = 3.38]
           |     |  |     |  |--leaf "Sample-196" 
           |     |  |     |  `--leaf "Sample-198" 
           |     |  |     `--leaf "Sample-204" 
           |     |  `--[dendrogram w/ 2 branches and 9 members at h = 7.6]
           |     |     |--[dendrogram w/ 2 branches and 5 members at h = 5.8]
           |     |     |  |--[dendrogram w/ 2 branches and 3 members at h = 3.88]
           |     |     |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.06]
           |     |     |  |  |  |--leaf "Sample-202" 
           |     |     |  |  |  `--leaf "Sample-200" 
           |     |     |  |  `--leaf "Sample-206" 
           |     |     |  `--[dendrogram w/ 2 branches and 2 members at h = 4.53]
           |     |     |     |--leaf "Sample-125" 
           |     |     |     `--leaf "Sample-189" 
           |     |     `--[dendrogram w/ 2 branches and 4 members at h = 5.26]
           |     |        |--[dendrogram w/ 2 branches and 3 members at h = 4.13]
           |     |        |  |--[dendrogram w/ 2 branches and 2 members at h = 3.12]
           |     |        |  |  |--leaf "Sample-192" 
           |     |        |  |  `--leaf "Sample-159" 
           |     |        |  `--leaf "Sample-181" 
           |     |        `--leaf "Sample-38" 
           |     `--[dendrogram w/ 2 branches and 8 members at h = 8.76]
           |        |--[dendrogram w/ 2 branches and 2 members at h = 2.78]
           |        |  |--leaf "Sample-47" 
           |        |  `--leaf "Sample-46" 
           |        `--[dendrogram w/ 2 branches and 6 members at h = 6.91]
           |           |--[dendrogram w/ 2 branches and 4 members at h = 5.32]
           |           |  |--[dendrogram w/ 2 branches and 2 members at h = 3.39]
           |           |  |  |--leaf "Sample-145" 
           |           |  |  `--leaf "Sample-135" 
           |           |  `--[dendrogram w/ 2 branches and 2 members at h = 4.2]
           |           |     |--leaf "Sample-114" 
           |           |     `--leaf "Sample-37" 
           |           `--[dendrogram w/ 2 branches and 2 members at h = 1.5]
           |              |--leaf "Sample-124" 
           |              `--leaf "Sample-120" 
           `--[dendrogram w/ 2 branches and 65 members at h = 18.8]
              |--[dendrogram w/ 2 branches and 51 members at h = 16.7]
              |  |--[dendrogram w/ 2 branches and 13 members at h = 11.6]
              |  |  |--[dendrogram w/ 2 branches and 7 members at h = 11.4]
              |  |  |  |--[dendrogram w/ 2 branches and 4 members at h = 5.66]
              |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 4.24]
              |  |  |  |  |  |--leaf "Sample-127" 
              |  |  |  |  |  `--leaf "Sample-45" 
              |  |  |  |  `--[dendrogram w/ 2 branches and 2 members at h = 3.79]
              |  |  |  |     |--leaf "Sample-134" 
              |  |  |  |     `--leaf "Sample-144" 
              |  |  |  `--[dendrogram w/ 2 branches and 3 members at h = 2.83]
              |  |  |     |--leaf "Sample-140" 
              |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.18]
              |  |  |        |--leaf "Sample-41" 
              |  |  |        `--leaf "Sample-130" 
              |  |  `--[dendrogram w/ 2 branches and 6 members at h = 5.41]
              |  |     |--[dendrogram w/ 2 branches and 3 members at h = 3.55]
              |  |     |  |--leaf "Sample-119" 
              |  |     |  `--[dendrogram w/ 2 branches and 2 members at h = 1.68]
              |  |     |     |--leaf "Sample-131" 
              |  |     |     `--leaf "Sample-141" 
              |  |     `--[dendrogram w/ 2 branches and 3 members at h = 4.24]
              |  |        |--[dendrogram w/ 2 branches and 2 members at h = 2.89]
              |  |        |  |--leaf "Sample-33" 
              |  |        |  `--leaf "Sample-110" 
              |  |        `--leaf "Sample-42" 
              |  `--[dendrogram w/ 2 branches and 38 members at h = 13.8]
              |     |--[dendrogram w/ 2 branches and 18 members at h = 12.5]
              |     |  |--[dendrogram w/ 2 branches and 12 members at h = 9.53]
              |     |  |  |--[dendrogram w/ 2 branches and 8 members at h = 7.75]
              |     |  |  |  |--[dendrogram w/ 2 branches and 3 members at h = 4.43]
              |     |  |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 2.9]
              |     |  |  |  |  |  |--leaf "Sample-48" 
              |     |  |  |  |  |  `--leaf "Sample-148" 
              |     |  |  |  |  `--leaf "Sample-1" 
              |     |  |  |  `--[dendrogram w/ 2 branches and 5 members at h = 6.59]
              |     |  |  |     |--[dendrogram w/ 2 branches and 3 members at h = 5.02]
              |     |  |  |     |  |--leaf "Sample-39" 
              |     |  |  |     |  `--[dendrogram w/ 2 branches and 2 members at h = 0.84]
              |     |  |  |     |     |--leaf "Sample-154" 
              |     |  |  |     |     `--leaf "Sample-116" 
              |     |  |  |     `--[dendrogram w/ 2 branches and 2 members at h = 3.48]
              |     |  |  |        |--leaf "Sample-36" 
              |     |  |  |        `--leaf "Sample-113" 
              |     |  |  `--[dendrogram w/ 2 branches and 4 members at h = 5.04]
              |     |  |     |--[dendrogram w/ 2 branches and 2 members at h = 2.61]
              |     |  |     |  |--leaf "Sample-177" 
              |     |  |     |  `--leaf "Sample-173" 
              |     |  |     `--[dendrogram w/ 2 branches and 2 members at h = 2.78]
              |     |  |        |--leaf "Sample-73" 
              |     |  |        `--leaf "Sample-155" 
              |     |  `--[dendrogram w/ 2 branches and 6 members at h = 4.78]
              |     |     |--[dendrogram w/ 2 branches and 3 members at h = 3.48]
              |     |     |  |--[dendrogram w/ 2 branches and 2 members at h = 2.64]
              |     |     |  |  |--leaf "Sample-157" 
              |     |     |  |  `--leaf "Sample-164" 
              |     |     |  `--leaf "Sample-179" 
              |     |     `--[dendrogram w/ 2 branches and 3 members at h = 3.03]
              |     |        |--[dendrogram w/ 2 branches and 2 members at h = 1.74]
              |     |        |  |--leaf "Sample-171" 
              |     |        |  `--leaf "Sample-185" 
              |     |        `--leaf "Sample-175" 
              |     `--[dendrogram w/ 2 branches and 20 members at h = 9.49]
              |        |--[dendrogram w/ 2 branches and 14 members at h = 7.23]
              |        |  |--[dendrogram w/ 2 branches and 10 members at h = 6.85]
              |        |  |  |--[dendrogram w/ 2 branches and 4 members at h = 6.57]
              |        |  |  |  |--[dendrogram w/ 2 branches and 2 members at h = 5.1]
              |        |  |  |  |  |--leaf "Sample-183" 
              |        |  |  |  |  `--leaf "Sample-167" 
              |        |  |  |  `--[dendrogram w/ 2 branches and 2 members at h = 4.44]
              |        |  |  |     |--leaf "Sample-32" 
              |        |  |  |     `--leaf "Sample-109" 
              |        |  |  `--[dendrogram w/ 2 branches and 6 members at h = 5.64]
              |        |  |     |--[dendrogram w/ 2 branches and 3 members at h = 3.15]
              |        |  |     |  |--[dendrogram w/ 2 branches and 2 members at h = 2.53]
              |        |  |     |  |  |--leaf "Sample-195" 
              |        |  |     |  |  `--leaf "Sample-197" 
              |        |  |     |  `--leaf "Sample-203" 
              |        |  |     `--[dendrogram w/ 2 branches and 3 members at h = 1.71]
              |        |  |        |--[dendrogram w/ 2 branches and 2 members at h = 1.4]
              |        |  |        |  |--leaf "Sample-199" 
              |        |  |        |  `--leaf "Sample-205" 
              |        |  |        `--leaf "Sample-201" 
              |        |  `--[dendrogram w/ 2 branches and 4 members at h = 3.85]
              |        |     |--leaf "Sample-160" 
              |        |     `--[dendrogram w/ 2 branches and 3 members at h = 2.41]
              |        |        |--leaf "Sample-188" 
              |        |        `--[dendrogram w/ 2 branches and 2 members at h = 1.1]
              |        |           |--leaf "Sample-191" 
              |        |           `--leaf "Sample-194" 
              |        `--[dendrogram w/ 2 branches and 6 members at h = 5.37]
              |           |--[dendrogram w/ 2 branches and 3 members at h = 2.5]
              |           |  |--leaf "Sample-72" 
              |           |  `--[dendrogram w/ 2 branches and 2 members at h = 1.38]
              |           |     |--leaf "Sample-153" 
              |           |     `--leaf "Sample-152" 
              |           `--[dendrogram w/ 2 branches and 3 members at h = 4.28]
              |              |--[dendrogram w/ 2 branches and 2 members at h = 2.61]
              |              |  |--leaf "Sample-150" 
              |              |  `--leaf "Sample-71" 
              |              `--leaf "Sample-151" 
              `--[dendrogram w/ 2 branches and 14 members at h = 12.5]
                 |--[dendrogram w/ 2 branches and 8 members at h = 7.08]
                 |  |--[dendrogram w/ 2 branches and 6 members at h = 5.35]
                 |  |  |--[dendrogram w/ 2 branches and 2 members at h = 1.41]
                 |  |  |  |--leaf "Sample-128" 
                 |  |  |  `--leaf "Sample-122" 
                 |  |  `--[dendrogram w/ 2 branches and 4 members at h = 3.09]
                 |  |     |--[dendrogram w/ 2 branches and 3 members at h = 2.79]
                 |  |     |  |--[dendrogram w/ 2 branches and 2 members at h = 1.94]
                 |  |     |  |  |--leaf "Sample-138" 
                 |  |     |  |  `--leaf "Sample-139" 
                 |  |     |  `--leaf "Sample-70" 
                 |  |     `--leaf "Sample-68" 
                 |  `--[dendrogram w/ 2 branches and 2 members at h = 5.12]
                 |     |--leaf "Sample-117" 
                 |     `--leaf "Sample-166" 
                 `--[dendrogram w/ 2 branches and 6 members at h = 6.27]
                    |--[dendrogram w/ 2 branches and 4 members at h = 4.93]
                    |  |--leaf "Sample-186" 
                    |  `--[dendrogram w/ 2 branches and 3 members at h = 3.18]
                    |     |--[dendrogram w/ 2 branches and 2 members at h = 2.18]
                    |     |  |--leaf "Sample-172" 
                    |     |  `--leaf "Sample-176" 
                    |     `--leaf "Sample-170" 
                    `--[dendrogram w/ 2 branches and 2 members at h = 3.85]
                       |--leaf "Sample-158" 
                       `--leaf "Sample-180" 

---

    Code
      str(p@row_dend_param$obj)
    Output
      --[dendrogram w/ 2 branches and 30 members at h = 91.4]
        |--[dendrogram w/ 2 branches and 8 members at h = 40.2]
        |  |--[dendrogram w/ 2 branches and 7 members at h = 38.2]
        |  |  |--leaf "Prevotella oralis et rel." 
        |  |  `--[dendrogram w/ 2 branches and 6 members at h = 27.8]
        |  |     |--leaf "Clostridium cellulosi et rel." 
        |  |     `--[dendrogram w/ 2 branches and 5 members at h = 16.8]
        |  |        |--leaf "Sporobacter termitidis et rel." 
        |  |        `--[dendrogram w/ 2 branches and 4 members at h = 14]
        |  |           |--[dendrogram w/ 2 branches and 2 members at h = 11.6]
        |  |           |  |--leaf "Clostridium orbiscindens et rel." 
        |  |           |  `--leaf "Ruminococcus obeum et rel." 
        |  |           `--[dendrogram w/ 2 branches and 2 members at h = 13.7]
        |  |              |--leaf "Subdoligranulum variable at rel." 
        |  |              `--leaf "Clostridium symbiosum et rel." 
        |  `--leaf "Bacteroides vulgatus et rel." 
        `--[dendrogram w/ 2 branches and 22 members at h = 49.2]
           |--[dendrogram w/ 2 branches and 7 members at h = 31.6]
           |  |--[dendrogram w/ 2 branches and 3 members at h = 15.6]
           |  |  |--[dendrogram w/ 2 branches and 2 members at h = 14.8]
           |  |  |  |--leaf "Allistipes et rel." 
           |  |  |  `--leaf "Parabacteroides distasonis et rel." 
           |  |  `--leaf "Bacteroides fragilis et rel." 
           |  `--[dendrogram w/ 2 branches and 4 members at h = 22.8]
           |     |--leaf "Bacteroides uniformis et rel." 
           |     `--[dendrogram w/ 2 branches and 3 members at h = 12.2]
           |        |--leaf "Prevotella tannerae et rel." 
           |        `--[dendrogram w/ 2 branches and 2 members at h = 8.81]
           |           |--leaf "Tannerella et rel." 
           |           `--leaf "Bacteroides plebeius et rel." 
           `--[dendrogram w/ 2 branches and 15 members at h = 39.7]
              |--[dendrogram w/ 2 branches and 5 members at h = 29.5]
              |  |--[dendrogram w/ 2 branches and 2 members at h = 10.2]
              |  |  |--leaf "Streptococcus bovis et rel." 
              |  |  `--leaf "Streptococcus mitis et rel." 
              |  `--[dendrogram w/ 2 branches and 3 members at h = 24.6]
              |     |--[dendrogram w/ 2 branches and 2 members at h = 20.7]
              |     |  |--leaf "Phascolarctobacterium faecium et rel." 
              |     |  `--leaf "Escherichia coli et rel." 
              |     `--leaf "Uncultured Mollicutes" 
              `--[dendrogram w/ 2 branches and 10 members at h = 27.4]
                 |--[dendrogram w/ 2 branches and 7 members at h = 19.2]
                 |  |--[dendrogram w/ 2 branches and 4 members at h = 15.7]
                 |  |  |--leaf "Lachnobacillus bovis et rel." 
                 |  |  `--[dendrogram w/ 2 branches and 3 members at h = 12.3]
                 |  |     |--[dendrogram w/ 2 branches and 2 members at h = 9.93]
                 |  |     |  |--leaf "Clostridium nexile et rel." 
                 |  |     |  `--leaf "Outgrouping clostridium cluster XIVa" 
                 |  |     `--leaf "Lachnospira pectinoschiza et rel." 
                 |  `--[dendrogram w/ 2 branches and 3 members at h = 14.7]
                 |     |--[dendrogram w/ 2 branches and 2 members at h = 10.8]
                 |     |  |--leaf "Clostridium sphenoides et rel." 
                 |     |  `--leaf "Bryantella formatexigens et rel." 
                 |     `--leaf "Anaerostipes caccae et rel." 
                 `--[dendrogram w/ 2 branches and 3 members at h = 20.1]
                    |--[dendrogram w/ 2 branches and 2 members at h = 17.6]
                    |  |--leaf "Bifidobacterium" 
                    |  `--leaf "Uncultured Clostridiales II" 
                    `--leaf "Akkermansia" 

# ord_explore_init stays the same

    Code
      ord_explore_init(dietswap)
    Output
      $data
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 9 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 4 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = unique tax_transform = identity
      
      $info
      $info$rank
      [1] "unique"
      
      $info$trans
      [1] "identity"
      
      $info$scale
      [1] NA
      
      $info$dist
      [1] "none"
      
      $info$ord
      [1] "auto"
      
      $info$constraints
      NULL
      
      $info$conditions
      NULL
      
      $info$isCon
      [1] FALSE
      
      
      $vars
      $vars$all
      [1] "subject"                "sex"                    "nationality"           
      [4] "group"                  "sample"                 "timepoint"             
      [7] "timepoint.within.group" "bmi_group"              "SAMPLE"                
      
      $vars$num
      [1] "timepoint"              "timepoint.within.group"
      
      $vars$cat
      [1] "subject"     "sex"         "nationality" "group"       "sample"     
      [6] "bmi_group"   "SAMPLE"     
      
      $vars$shapeSafe
      [1] "sex"                    "nationality"            "group"                 
      [4] "timepoint.within.group" "bmi_group"             
      
      
      $ranks
      [1] "Phylum" "Family" "Genus"  "unique"
      
      $warn
      [1] FALSE
      

---

    Code
      ord_explore_init(ord)
    Output
      $data
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 11 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 4 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = clr
      
      ordination of class: rda cca 
      rda(formula = OTU ~ weight + female, data = data)
      constraints: weight+female
      
      $counts OTU Table: [ 130 taxa and 222 samples ]
      
      $info
      $info$rank
      [1] "Genus"
      
      $info$trans
      [1] "clr"
      
      $info$scale
      [1] NA
      
      $info$dist
      [1] "none"
      
      $info$ord
      [1] "RDA"
      
      $info$constraints
      [1] "weight" "female"
      
      $info$conditions
      NULL
      
      $info$isCon
      [1] TRUE
      
      
      $vars
      $vars$all
       [1] "subject"                "sex"                    "nationality"           
       [4] "group"                  "sample"                 "timepoint"             
       [7] "timepoint.within.group" "bmi_group"              "weight"                
      [10] "female"                 "SAMPLE"                
      
      $vars$num
      [1] "timepoint"              "timepoint.within.group" "weight"                
      [4] "female"                
      
      $vars$cat
      [1] "subject"     "sex"         "nationality" "group"       "sample"     
      [6] "bmi_group"   "SAMPLE"     
      
      $vars$shapeSafe
      [1] "sex"                    "nationality"            "group"                 
      [4] "timepoint.within.group" "bmi_group"              "weight"                
      [7] "female"                
      
      
      $ranks
      [1] "Phylum" "Family" "Genus"  "unique"
      
      $warn
      [1] FALSE
      

---

    Code
      ord_explore_init(esophagus)
    Message <simpleMessage>
      Note: Replacing missing sample_data with a dataframe of only sample_names.
      Try `ps <- phyloseq_validate(ps, verbose = FALSE)` to avoid this message
      Note: Replacing missing tax_table with a 1-column table of only taxa_names.
      Try `ps <- phyloseq_validate(ps, verbose = FALSE)` to avoid this message
    Output
      $data
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 58 taxa and 3 samples ]
      sample_data() Sample Data:       [ 3 samples by 1 sample variables ]
      tax_table()   Taxonomy Table:    [ 58 taxa by 1 taxonomic ranks ]
      phy_tree()    Phylogenetic Tree: [ 58 tips and 57 internal nodes ]
      
      ps_extra info:
      tax_agg = unique tax_transform = identity
      
      $info
      $info$rank
      [1] "unique"
      
      $info$trans
      [1] "identity"
      
      $info$scale
      [1] NA
      
      $info$dist
      [1] "none"
      
      $info$ord
      [1] "auto"
      
      $info$constraints
      NULL
      
      $info$conditions
      NULL
      
      $info$isCon
      [1] FALSE
      
      
      $vars
      $vars$all
      [1] "SAMPLE"
      
      $vars$num
      character(0)
      
      $vars$cat
      NULL
      
      $vars$shapeSafe
      NULL
      
      
      $ranks
      [1] "unique"
      
      $warn
      [1] FALSE
      

# dist_choices helper works

    Code
      dist_choices(dietswap, type = "tree")
    Output
      named character(0)

---

    Code
      dist_choices(esophagus, type = "tree")
    Output
             gunifrac: Generalised UniFrac, alpha=0.5 
                                           "gunifrac" 
                           wunifrac: weighted UniFrac 
                                           "wunifrac" 
                          unifrac: unweighted UniFrac 
                                            "unifrac" 
      va-wunifrac: variance-adjusted weighted UniFrac 
                                        "va-wunifrac" 
                                                dpcoa 
                                              "dpcoa" 

# ord_choices helper works

    Code
      cat(type)
    Output
      all
    Code
      cat(ord_choices(type))
    Output
      auto PCA PCoA RDA CAP CCA NMDS

---

    Code
      cat(type)
    Output
      constrained
    Code
      cat(ord_choices(type))
    Output
      auto RDA CAP CCA

---

    Code
      cat(type)
    Output
      unconstrained
    Code
      cat(ord_choices(type))
    Output
      auto PCA PCoA NMDS

---

    Code
      cat(type)
    Output
      dist
    Code
      cat(ord_choices(type))
    Output
      auto PCoA CAP NMDS

---

    Code
      cat(type)
    Output
      noDist
    Code
      cat(ord_choices(type))
    Output
      auto PCA RDA CCA

---

    Code
      cat(type)
    Output
      constrained dist
    Code
      cat(ord_choices(type))
    Output
      auto CAP

---

    Code
      cat(type)
    Output
      unconstrained dist
    Code
      cat(ord_choices(type))
    Output
      auto PCoA NMDS

---

    Code
      cat(type)
    Output
      constrained noDist
    Code
      cat(ord_choices(type))
    Output
      auto RDA CCA

---

    Code
      cat(type)
    Output
      unconstrained noDist
    Code
      cat(ord_choices(type))
    Output
      auto PCA

# ord_code helper works

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      FALSE 	0.5 	 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        colour = "v", fill = "v",
        shape = "var", alpha = 0.5,
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      FALSE 	0.5 	test1 test2 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        constraints = c("test1", "test2"),
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        colour = "v", fill = "v",
        shape = "var", alpha = 0.5,
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      FALSE 	aVariable 	 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        colour = "v", fill = "v",
        shape = "var", alpha = "aVariable",
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      FALSE 	aVariable 	test1 test2 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        constraints = c("test1", "test2"),
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        colour = "v", fill = "v",
        shape = "var", alpha = "aVariable",
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      1 2 3 4 5 6 	0.5 	 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        plot_taxa = 1:6,
        colour = "v", fill = "v",
        shape = "var", alpha = 0.5,
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      1 2 3 4 5 6 	0.5 	test1 test2 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        constraints = c("test1", "test2"),
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        plot_taxa = 1:6,
        colour = "v", fill = "v",
        shape = "var", alpha = 0.5,
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      1 2 3 4 5 6 	aVariable 	 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        plot_taxa = 1:6,
        colour = "v", fill = "v",
        shape = "var", alpha = "aVariable",
        size = 1
       )

---

    Code
      for (x in list(p, a, c)) cat(x, "\t")
    Output
      1 2 3 4 5 6 	aVariable 	test1 test2 	
    Code
      ord_code(rank = "Genus", trans = "identity", dist = "none", ord = "RDA", const = c,
        conds = NULL, x = 1, y = 2, colour = "v", fill = "v", shape = "var", alpha = a,
        size = 1, plot_taxa = p, ellipses = FALSE, chulls = FALSE, paths = NULL)
    Output
      your_phyloseq %>%
       tax_transform(rank = "Genus", trans = "identity") %>%
       ord_calc(
        constraints = c("test1", "test2"),
        method = "RDA"
       ) %>% 
       ord_plot(
        axes = c(1, 2),
        plot_taxa = 1:6,
        colour = "v", fill = "v",
        shape = "var", alpha = "aVariable",
        size = 1
       )

# ord_code_dist helper works

    Code
      cat(ord_code_dist("aitchison"))
    Output
       dist_calc(dist = "aitchison") %>%

---

    Code
      cat(ord_code_dist("none"))

# ord_code_stat and paths helpers work

    Code
      cat(ord_code_stat(ellipses = TRUE, chulls = FALSE, colour = "aVar"))
    Output
       ) +
       ggplot2::stat_ellipse(
        ggplot2::aes(colour = aVar)
       )

---

    Code
      cat(ord_code_stat(ellipses = FALSE, chulls = FALSE, colour = "aVar"))
    Output
       )

---

    Code
      cat(ord_code_stat(ellipses = FALSE, chulls = TRUE, colour = "aVar"))
    Output
       ) +
       stat_chull(
        ggplot2::aes(colour = aVar)
       )

---

    Code
      cat(ord_code_paths(paths = list(colour = "aVar", id_var = "bVar", id_values = letters[
        1:4], all_vars = "aVar")))
    Output
       ) %>%
       add_paths(
        id_var = "bVar", 
        id_values = c("a", "b", "c", "d"),
        mapping = ggplot2::aes(colour = aVar)
       )

---

    Code
      cat(ord_code_paths(paths = list(colour = "aVar", id_var = "bVar", id_values = letters[
        1:4], all_vars = c("otherVar", "anotherVar"))))
    Output
       ) %>%
       add_paths(
        id_var = "bVar", 
        id_values = c("a", "b", "c", "d"),
        colour = "aVar"
       )

# ord_build works

    Code
      ord_build(data = dietswap, rank = "Genus", trans = "identity", dist = "bray",
        method = "PCoA", constraints = NULL, conditions = NULL)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = identity
      
      bray distance matrix of size 222 
      0.7639533 0.7851213 0.6680796 0.7699252 0.80507 ...
      
      ordination of class: capscale rda cca 
      capscale(formula = distance ~ 1, data = data)
      

---

    Code
      ord_build(data = dietswap, rank = "Genus", trans = "clr", dist = NA, method = "auto",
        constraints = NULL, conditions = NULL)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = clr
      
      ordination of class: rda cca 
      rda(formula = OTU ~ 1, data = data)
      
      
      $counts OTU Table: [ 130 taxa and 222 samples ]

# ord_explore_palet_fun works

    Code
      ord_explore_palet_fun(dietswap, "Genus")
    Output
         Prevotella melaninogenica et rel.   Oscillospira guillermondii et rel. 
                                 "#A6CEE3"                            "#1F78B4" 
              Bacteroides vulgatus et rel.        Clostridium cellulosi et rel. 
                                 "#B2DF8A"                            "#33A02C" 
                 Prevotella oralis et rel. Faecalibacterium prausnitzii et rel. 
                                 "#FB9A99"                            "#E31A1C" 
            Sporobacter termitidis et rel.        Clostridium symbiosum et rel. 
                                 "#FDBF6F"                            "#FF7F00" 
                        Allistipes et rel.     Clostridium orbiscindens et rel. 
                                 "#CAB2D6"                            "#6A3D9A" 
          Subdoligranulum variable at rel.           Ruminococcus obeum et rel. 
                                 "#FFFF99"                            "#B15928" 
            Butyrivibrio crossotus et rel.         Bacteroides fragilis et rel. 
                                 "#1ff8ff"                            "#1B9E77" 
                               Akkermansia           Bacteroides ovatus et rel. 
                                 "#D95F02"                            "#7570B3" 
        Parabacteroides distasonis et rel.        Dorea formicigenerans et rel. 
                                 "#E7298A"                            "#66A61E" 
             Bacteroides uniformis et rel.                            Dialister 
                                 "#E6AB02"                            "#A6761D" 
          Bryantella formatexigens et rel.           Uncultured Clostridiales I 
                                 "#666666"                            "#4b6a53" 
              Coprococcus eutactus et rel.           Clostridium leptum et rel. 
                                 "#b249d5"                            "#7edc45" 
            Clostridium sphenoides et rel.             Escherichia coli et rel. 
                                 "#5c47b8"                            "#cfd251" 
               Streptococcus bovis et rel.          Uncultured Clostridiales II 
                                 "#ff69b4"                            "#69c86c" 
                           Bifidobacterium    Anaerotruncus colihominis et rel. 
                                 "#cd3e50"                            "#83d5af" 
         Lachnospira pectinoschiza et rel.          Anaerostipes caccae et rel. 
                                 "#da6130"                            "#5e79b2" 
             Ruminococcus callidus et rel.      Bacteroides splachnicus et rel. 
                                 "#c29545"                            "#532a5a" 
               Ruminococcus bromii et rel.          Prevotella tannerae et rel. 
                                 "#5f7b35"                            "#c497cf" 
              Lachnobacillus bovis et rel.          Eubacterium rectale et rel. 
                                 "#773a27"                            "#7cb9cb" 
            Mitsuokella multiacida et rel. Outgrouping clostridium cluster XIVa 
                                 "#594e50"                            "#d3c4a8" 
                Clostridium nexile et rel.                                other 
                                 "#c17e7f"                             "grey90" 

---

    Code
      ord_explore_palet_fun(ps = dietswap, tax_level = "Family", top_by = median,
        other = "colourz")
    Output
                  Bacteroidetes    Clostridium cluster IV  Clostridium cluster XIVa 
                      "#A6CEE3"                 "#1F78B4"                 "#B2DF8A" 
                 Proteobacteria    Clostridium cluster IX                   Bacilli 
                      "#33A02C"                 "#FB9A99"                 "#E31A1C" 
       Uncultured Clostridiales            Actinobacteria           Verrucomicrobia 
                      "#FDBF6F"                 "#FF7F00"                 "#CAB2D6" 
         Clostridium cluster XI     Clostridium cluster I   Clostridium cluster XVI 
                      "#6A3D9A"                 "#FFFF99"                 "#B15928" 
          Uncultured Mollicutes Clostridium cluster XVIII              Fusobacteria 
                      "#1ff8ff"                 "#1B9E77"                 "#D95F02" 
        Clostridium cluster III  Clostridium cluster XIII    Clostridium cluster XV 
                      "#7570B3"                 "#E7298A"                 "#66A61E" 
       Clostridium cluster XVII            Asteroleplasma              Spirochaetes 
                      "#E6AB02"                 "#A6761D"                 "#666666" 
                  Cyanobacteria                     other 
                      "#4b6a53"                 "colourz" 

# full join stays the same

    Code
      ps_join(x = x, y = y, match_sample_names = "ID_var", type = j)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 553 taxa and 280 samples ]
      sample_data() Sample Data:       [ 280 samples by 12 sample variables ]
      tax_table()   Taxonomy Table:    [ 553 taxa by 1 taxonomic ranks ]

# inner join stays the same

    Code
      ps_join(x = x, y = y, match_sample_names = "ID_var", type = j)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 533 taxa and 100 samples ]
      sample_data() Sample Data:       [ 100 samples by 12 sample variables ]
      tax_table()   Taxonomy Table:    [ 533 taxa by 1 taxonomic ranks ]

# left join stays the same

    Code
      ps_join(x = x, y = y, match_sample_names = "ID_var", type = j)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 553 taxa and 280 samples ]
      sample_data() Sample Data:       [ 280 samples by 12 sample variables ]
      tax_table()   Taxonomy Table:    [ 553 taxa by 1 taxonomic ranks ]

# anti join stays the same

    Code
      ps_join(x = x, y = y, match_sample_names = "ID_var", type = j)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 549 taxa and 180 samples ]
      sample_data() Sample Data:       [ 180 samples by 10 sample variables ]
      tax_table()   Taxonomy Table:    [ 549 taxa by 1 taxonomic ranks ]

# semi join stays the same

    Code
      ps_join(x = x, y = y, match_sample_names = "ID_var", type = j)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 533 taxa and 100 samples ]
      sample_data() Sample Data:       [ 100 samples by 10 sample variables ]
      tax_table()   Taxonomy Table:    [ 533 taxa by 1 taxonomic ranks ]

# tax_fix defaults allow agg: enterotype Genus

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 553 taxa and 280 samples ]
      sample_data() Sample Data:       [ 280 samples by 9 sample variables ]
      tax_table()   Taxonomy Table:    [ 553 taxa by 1 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = NA

# tax_fix defaults allow agg: modified_dietswap Phylum

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 12 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 12 taxa by 1 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Phylum tax_transform = NA

# tax_fix defaults allow agg: modified_dietswap Family

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 28 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 28 taxa by 2 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Family tax_transform = NA

# tax_fix defaults allow agg: modified_dietswap Genus

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 110 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 110 taxa by 3 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Kingdom

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 1 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 1 taxa by 1 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Kingdom tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Phylum

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 9 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 9 taxa by 2 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Phylum tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Class

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 19 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 19 taxa by 3 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Class tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Order

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 32 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 32 taxa by 4 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Order tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Family

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 59 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 59 taxa by 5 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Family tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Genus

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 178 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 178 taxa by 6 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = NA

# tax_fix defaults allow agg: ibd_phylo Species

    Code
      tax_agg(ps = fixed, rank = r)
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 178 taxa and 91 samples ]
      sample_data() Sample Data:       [ 91 samples by 15 sample variables ]
      tax_table()   Taxonomy Table:    [ 178 taxa by 7 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Species tax_transform = NA

# list of plots generated by group_by works

    Code
      print(p)
    Output
      [1] "overweight"

---

    Code
      print(p)
    Output
      [1] "lean"

---

    Code
      print(p)
    Output
      [1] "obese"

# tt_add_topN_var helper works as expected

    Code
      print(head(x = tt, 6), width = 80)
    Output
                                   Phylum            Family           
      Actinomycetaceae             "Actinobacteria"  "Actinobacteria" 
      Aerococcus                   "Firmicutes"      "Bacilli"        
      Aeromonas                    "Proteobacteria"  "Proteobacteria" 
      Akkermansia                  "Verrucomicrobia" "Verrucomicrobia"
      Alcaligenes faecalis et rel. "Proteobacteria"  "Proteobacteria" 
      Allistipes et rel.           "Bacteroidetes"   "Bacteroidetes"  
                                   Genus                          test              
      Actinomycetaceae             "Actinomycetaceae"             "Actinomycetaceae"
      Aerococcus                   "Aerococcus"                   "Aerococcus"      
      Aeromonas                    "Aeromonas"                    "Aeromonas"       
      Akkermansia                  "Akkermansia"                  "Akkermansia"     
      Alcaligenes faecalis et rel. "Alcaligenes faecalis et rel." "things"          
      Allistipes et rel.           "Allistipes et rel."           "things"          
      attr(,"class")
      [1] "taxonomyTable"
      attr(,"class")attr(,"package")
      [1] "phyloseq"

# microbiome's dietswap data hasn't changed

    Code
      dietswap
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]

# microbiome::aggregate_taxa output hasn't changed: Phylum

    Code
      microbiome::aggregate_taxa(x = dietswap, level = level)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 8 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 8 taxa by 2 taxonomic ranks ]

# microViz::tax_agg output hasn't changed: Phylum

    Code
      ps_get(tax_agg(ps = dietswap, level, sort_by = "name", add_unique = TRUE))
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 8 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 8 taxa by 2 taxonomic ranks ]

# microbiome::aggregate_taxa output hasn't changed: Family

    Code
      microbiome::aggregate_taxa(x = dietswap, level = level)
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 22 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 22 taxa by 3 taxonomic ranks ]

# microViz::tax_agg output hasn't changed: Family

    Code
      ps_get(tax_agg(ps = dietswap, level, sort_by = "name", add_unique = TRUE))
    Output
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 22 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 22 taxa by 3 taxonomic ranks ]

# tax_fix error prompt looks right

    Code
      cat(taxFixPrompt())
    Output
      
      
      To fix the problem, try:
        `yourData %>% tax_fix()`
      
      Try tax_fix_interactive() to find and fix further problems

---

    Code
      cat(taxFixPrompt(unknowns = c("anUnknown", "another")))
    Output
      
      
      To fix the problem, try:
        `yourData %>% tax_fix(unknowns = c("anUnknown", "another"))`
      
      Try tax_fix_interactive() to find and fix further problems

# tax_agg errors on NAs, '', or convergent values

    Code
      tax_agg(dietswap, rank = "Genus")
    Error <simpleError>
      NAs in tax_table at rank: Genus
      
      To fix the problem, try:
        `yourData %>% tax_fix()`
      
      Try tax_fix_interactive() to find and fix further problems

---

    Code
      tax_agg(dietswap, rank = "Genus")
    Error <simpleError>
      zero-length name(s) in tax_table at rank: Genus
      
      To fix the problem, try:
        `yourData %>% tax_fix()`
      
      Try tax_fix_interactive() to find and fix further problems

---

    Code
      tax_agg(dietswap, rank = "Genus")
    Message <simpleMessage>
      Problematic Genus values detected in tax_table:
      g__
      -
      taxa_name / Phylum / Family / Genus
      Aeromonas / Proteobacteria / Proteobacteria / g__
      Akkermansia / Verrucomicrobia / Verrucomicrobia / g__
      Allistipes et rel. / Bacteroidetes / Bacteroidetes / g__
      Anaerofustis / Firmicutes / Clostridium cluster XV / g__
      Anaerostipes caccae et rel. / Firmicutes / Clostridium cluster XIVa / g__
      Anaerotruncus colihominis et rel. / Firmicutes / Clostridium cluster IV / g__
      -
    Error <simpleError>
      Taxa not unique at rank: Genus
      See last messages for convergent taxa rows.
      
      To fix the problem, try:
        `yourData %>% tax_fix(unknowns = c("g__"))`
      
      Try tax_fix_interactive() to find and fix further problems

# constrained rda plot gives correct positions

    Code
      cat(p$data[1:50, 1, drop = TRUE])
    Output
      -2.130476 0.5894143 0.3160391 1.08751 -0.05355557 3.503399 1.203936 0.6234953 5.592124 2.525964 2.730401 1.39254 3.615096 2.609643 0.7825542 0.6269883 -0.3570368 -0.335859 0.5819399 0.4438977 1.663886 0.149256 0.3422652 0.2049955 0.7025216 -0.3234688 0.3380817 2.476972 1.913976 1.312951 1.080631 -1.103206 -1.392823 -0.2848868 1.558326 -1.85179 0.5686295 -2.050395 -1.256443 1.342032 -1.156612 -1.018807 -1.394098 0.1725057 -0.04134601 -0.6274757 -1.368677 -2.246543 -1.530365 1.691596

---

    Code
      cat(p$data[1:50, 2, drop = TRUE])
    Output
      2.618072 -2.748471 -3.076232 -1.148784 0.3788259 -0.1015143 -1.35111 -1.367386 -0.5596091 -3.083707 -1.623868 -1.807981 -0.5441718 -0.5529239 0.9604805 -1.28105 -1.098181 -2.594557 -1.21884 -1.861085 -1.052584 -0.4884243 -0.5190012 -2.966743 -2.431879 -2.63563 -1.925107 -1.227659 -0.7474585 -2.331083 -0.1704811 3.099588 3.292661 -1.575772 -1.654633 2.536926 1.863032 -0.03288809 1.431886 -0.3492331 2.171615 2.694277 -1.276792 -1.185101 3.757507 0.6055251 0.3254076 1.146757 -1.454654 -1.136692

---

    Code
      cat(p$layers[[2]]$data[, 1, drop = TRUE])
    Output
      -0.3565107 -0.1062509 0.1547571 -0.8509554 -1.899592 -0.2096152 -0.2829297 -1.779616 -0.4937304 0.5586688 0.3846523 0.4645738 -0.059117 -2.16059 -1.011784 -1.82159 -0.9692845 -1.742553 0.6740945 -1.020301 -2.346078 -0.2790361 -0.1204818 -0.5271442 -1.78475 0.6517035 -0.5044418 -0.9840494 -0.07047137 2.96672 0.4707 1.574522 -2.596042 1.221985 0.1057324 -0.902472 1.130822 -0.1781154 -0.159024 0.1011578 -0.02875207 0.8307857 -0.6201586 -0.3972224 -0.04517522 1.492099 0.771006 0.3775464 0.1834584 0.6782735 0.4648295 3.962965 2.859382 0.5720216 0.08039252 0.527929 -0.6067356 -0.6899561 -1.121057 -0.8502372 0.1333341 -0.4010931 0.5326941 1.674685 -0.08324529 1.282892 -0.337443 -1.879871 -0.2650652 0.1850139 0.1382349 0.4058195 0.1380115 -0.2692422 0.1113067 4.776038 3.61217 0.8894821 -0.2152379 -0.414358 -1.159676 -0.8853212 -1.144903 -0.4437925 -1.527541 -1.223474 -0.2949115 -0.3068263 1.469021 -0.6934613 -0.2475891 0.02910698 -0.2801657 -0.2688476 0.2432086 -0.5492048 -0.3459145 0.5510748 -0.03510124 0.02694244 -0.1714998 1.310444 2.871598 -1.903205 -0.09571269 1.979713 2.077833 3.246218 -1.032558 -0.8310774 -0.8383837 -1.479313 -0.6247188 -1.976837 -0.4051784 -0.7227111 1.405419 0.634562 0.2666264 0.3054786 0.00385781 -0.7068847 0.07462064

---

    Code
      cat(p$layers[[3]]$data[, 1, drop = TRUE])
    Output
      4.776038 3.962965 2.96672

---

    Code
      cat(p$layers[[3]]$data[, 2, drop = TRUE])
    Output
      -0.09814972 0.2533772 -2.395719

# partialed bray CAP plot gives correct positions

    Code
      cat(p2$data[1:50, 1, drop = TRUE])
    Output
      -0.9433384 -0.4360683 0.01614725 -1.363398 1.665281 1.65784 1.362464 1.129672 1.881359 0.3532599 1.631433 1.388797 1.250819 1.72303 0.8478204 -0.8442538 -0.7064204 -1.10628 -0.2730612 0.4129244 -0.7327066 -0.6039317 1.781876 -2.437726 -0.7128981 -0.5591445 -0.5503607 2.121992 -0.7871148 -0.2522774 0.6067446 0.08712734 2.221565 -1.241554 0.6175632 1.628865 2.537898 -0.03497828 -1.09659 3.99456 -0.07894219 2.536319 -0.8418586 -0.2491924 0.2117898 -0.05593509 -1.513947 -1.795448 -2.313529 -0.2527132

---

    Code
      cat(p2$data[1:50, 2, drop = TRUE])
    Output
      -0.777159 0.319758 0.6808594 -0.3605004 0.7990729 0.9501213 -0.2204894 0.7256214 -0.2381873 0.3182067 0.7976874 0.6710309 -0.007768279 -0.3096304 0.9760564 0.2008539 0.3406989 -0.5329119 -0.3638514 0.7576985 -0.1883829 -0.3308651 0.7528426 0.05331558 0.06812347 -0.5354913 0.0230663 0.8174722 -0.331778 -0.3277457 0.8594061 -0.8416075 -0.7947944 -0.4042679 0.1839848 -0.7063301 -0.6911026 0.6860141 -0.6901189 -0.3483738 -0.7896026 -0.690313 -0.05402992 0.6581042 -0.7226321 0.4378195 0.06440707 -0.8558818 0.4833695 0.7706338

---

    Code
      cat(p2$layers[[2]]$data[, 1, drop = TRUE])
    Output
      2.686247

# aitchison plot hasn't changed

    Code
      cat(abs(p3$data[1:50, 1, drop = TRUE]))
    Output
      1.106491 0.8693114 0.7769138 0.1722429 0.5494207 1.342983 0.6422256 0.6158025 1.381771 1.076831 1.047452 0.9241626 1.112544 0.6735143 0.9536031 0.7115146 0.4646272 0.09576095 0.1244906 0.7623212 0.6518461 0.09280437 0.7242988 0.5598014 0.8549104 0.2034272 0.2821326 1.169105 0.3978237 0.5269234 0.8514695 1.310387 1.287183 0.0458201 0.5668595 1.1262 0.6357289 0.044228 0.9423226 0.0271642 1.315038 1.145424 0.1641273 0.2225042 0.8148921 0.1271763 0.3235937 0.8928878 0.0088699 0.9016179

---

    Code
      cat(abs(p3$data[1:50, 2, drop = TRUE]))
    Output
      0.5787119 1.110654 1.004009 1.55686 0.476133 0.3927158 1.640435 0.4383189 0.3547506 1.092631 0.5076862 1.161595 1.458424 2.499746 0.4842615 0.7047841 0.3034085 1.131942 1.505218 0.7052604 1.343786 1.922793 0.1882586 0.7064074 1.607724 1.576125 0.8710285 0.9390623 1.286664 2.279908 0.4152985 0.2621039 0.2757177 0.1157476 0.2378094 0.2285831 0.1264438 0.9902263 0.1533971 1.291777 0.08253041 0.191776 0.2938356 0.6350377 1.308276 1.170076 1.064223 0.4547523 0.06588292 0.6282335

---

    Code
      p3$layers
    Output
      [[1]]
      mapping: colour = ~bmi_group 
      geom_point: na.rm = FALSE
      stat_identity: na.rm = FALSE
      position_identity 
      

# clr PCA plot hasn't changed

    Code
      cat(p4$data[1:50, 1, drop = TRUE])
    Output
      -1.106491 0.8693114 0.7769138 0.1722429 0.5494207 1.342983 0.6422256 0.6158025 1.381771 1.076831 1.047452 0.9241626 1.112544 0.6735143 0.9536031 0.7115146 0.4646272 -0.09576095 0.1244906 0.7623212 0.6518461 0.09280437 0.7242988 0.5598014 0.8549104 -0.2034272 0.2821326 1.169105 0.3978237 0.5269234 0.8514695 -1.310387 -1.287183 -0.0458201 0.5668595 -1.1262 -0.6357289 -0.044228 -0.9423226 -0.0271642 -1.315038 -1.145424 -0.1641273 0.2225042 -0.8148921 -0.1271763 -0.3235937 -0.8928878 0.0088699 0.9016179

---

    Code
      cat(p4$data[1:50, 2, drop = TRUE])
    Output
      -0.5787119 1.110654 1.004009 1.55686 -0.476133 -0.3927158 1.640435 0.4383189 0.3547506 1.092631 0.5076862 1.161595 1.458424 2.499746 -0.4842615 0.7047841 0.3034085 1.131942 1.505218 0.7052604 1.343786 1.922793 0.1882586 0.7064074 1.607724 1.576125 0.8710285 0.9390623 1.286664 2.279908 -0.4152985 -0.2621039 0.2757177 0.1157476 -0.2378094 0.2285831 -0.1264438 -0.9902263 0.1533971 1.291777 -0.08253041 0.191776 -0.2938356 -0.6350377 -1.308276 -1.170076 -1.064223 -0.4547523 0.06588292 -0.6282335

---

    Code
      p4$layers
    Output
      [[1]]
      mapping: colour = ~bmi_group 
      geom_point: na.rm = FALSE
      stat_identity: na.rm = FALSE
      position_identity 
      

# tax_transform 'maaslin2-default' chaining works

    Code
      ord
    Output
      ps_extra object - a list with phyloseq and extras:
      
      phyloseq-class experiment-level object
      otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
      sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
      tax_table()   Taxonomy Table:    [ 130 taxa by 3 taxonomic ranks ]
      
      ps_extra info:
      tax_agg = Genus tax_transform = compositional&log2
      
      ordination of class: rda cca 
      rda(formula = OTU ~ 1, data = data)
      
      
      $counts OTU Table: [ 130 taxa and 222 samples ]

# unifrac distances work

    Code
      dist_get(dist_calc(esophagus, dist = d))
    Output
                B         C
      C 0.4404284          
      D 0.4332325 0.4969773

---

    Code
      dist_get(dist_calc(esophagus, dist = d))
    Output
                B         C
      C 0.5175550          
      D 0.5182284 0.5422394

---

    Code
      dist_get(dist_calc(esophagus, dist = d))
    Output
                B         C
      C 0.2035424          
      D 0.2603371 0.2477016

---

    Code
      dist_get(dist_calc(esophagus, dist = d))
    Output
                B         C
      C 0.3605292          
      D 0.3594336 0.3610123

---
title: 'microViz: an R package for microbiome data visualization and statistics'
authors:
- name: David J.M. Barnett
  orcid: 0000-0003-1961-7206
  affiliation: "1, 2"
- name: Ilja C.W. Arts
  orcid: 0000-0001-6462-6692
  affiliation: 1
- name: John Penders
  orcid: 0000-0001-9146-5919
  affiliation: "2, 3, 4"
date: "7th April 2021"
affiliations:
- name: Maastricht Centre for Systems Biology (MaCSBio)
  index: 1
- name: Maastricht University Medical Center+, Department of Medical Microbiology
  index: 2
- name: Maastricht University, Public Health Research Institute (CAPHRI)
  index: 3
- name: Maastricht University, School of Nutrition and Translational Research in Metabolism (NUTRIM)
  index: 4
bibliography: microViz.bib
tags:
- R
- microbiome
- visualization
---

# Summary

microViz is an R package for the statistical analysis and visualization of microbiota data. This package extends the functionality of popular microbial ecosystem data analysis R packages, including phyloseq [@mcmurdie_phyloseq_2013], vegan [@oksanen_vegan_2020] and microbiome [@lahti_microbiome_2012]. microViz provides a selection of powerful additions to the toolbox of researchers already familiar with phyloseq and microbiome, as well as assisting researchers with less R programming experience to independently explore and analyse their data and to generate publication-ready figures.

The tools offered by microViz include:

-   A Shiny app [@shiny] for interactive exploration of microbiota data within R, pairing ordination plots with abundance bar charts
-   Easy to use functions for generating publication-ready ordination plots with ggplot2 [@wickham_ggplot2_2016], accommodating constrained and partial ordination, bi-plots and tri-plots, and automatic captioning designed to promote methodological transparency and reproducibility
-   A novel visualization approach pairing ordination plots with circular bar charts (iris plots) for comprehensive, intuitive and compact visualization of the similarity and composition of hundreds of microbial ecosystems
-   Correlation and composition heatmaps for microbiome data annotated with plots showing each taxon's prevalence and/or abundance
-   A compact cladogram visualization approach for intuitive comparison of numerous microbe-metadata associations derived from (multivariable) statistical models (taxonomic association trees)

# Statement of need

Modern microbiome research typically involves the use of next-generation-sequencing methods to profile the relative abundance of hundreds of microbial taxa across tens, hundreds or thousands of samples. Alongside increasing sample sizes, the amount of relevant metadata collected is growing, particularly in human cohort studies. These trends all increase the size and complexity of the resulting dataset, which makes its exploration, statistical analysis, and presentation increasingly challenging.

Ordination methods like Principle Coordinates/Components Analyses (PCoA/PCA) are a staple method in microbiome research. The vegan R package implements the majority of distance and ordination calculation methods and microViz makes available two further dissimilarity measures: Generalized UniFrac from the GUniFrac package, and the Aitchison distance. The former provides a balanced intermediate between unweighted and weighted UniFrac [@chen_2012], and the latter is a distance measure designed for use with compositional data, such as sequencing read counts [@gloor_2017].

The phyloseq R package provides an interface for producing ordination plots, with the ggplot2 R package. microViz streamlines the computation and presentation of ordination methods including the constrained analyses: redundancy analysis (RDA), distance-based RDA, partial RDA, and canonical correspondence analysis (CCA). microViz can generate highly customizable ggplot2 bi-plots and tri-plots, showing labelled arrows for microbial loadings and constraint variables when applicable. Furthermore, these figures are captioned automatically, by default. The captions are intended to promote better reporting of ordination methods in published research, where too often insufficient information is given to reproduce the ordination plot. To provide the automated captioning, microViz implements a simple S3 list class, ps_extra, for provenance tracking, by storing distance matrices and ordination objects alongside the phyloseq object they were created from, as well as relevant taxonomic aggregation and transformation information.

Moreover, microViz provides a Shiny app interface [@shiny] that allows the user to interactively create and explore ordination plots directly from phyloseq objects. The Shiny app generates code that can be copy-pasted into a script to reproduce the interactively designed ordination plot. The user can click and drag on the interactive ordination plot to select samples and directly examine their taxonomic compositions on a customizable stacked bar chart with a clear colour scheme.

Alternatively, for a comprehensive and intuitive static presentation of both sample variation patterns and underlying microbial composition, microViz provides an easy approach to pair ordination plots with attractive circular bar charts (iris plots) by ordering the bar chart in accordance with the rotational position of samples around the origin point on the ordination plot, e.g. \autoref{fig:one}. Bar charts do have limitations when visualizing highly diverse samples, such as the adult gut microbiome, at a detailed taxonomic level. This is why microViz also offers an enhanced heatmap visualization approach, pairing an ordered heatmap of (transformed and/or scaled) microbial abundances with compact plots showing each taxon's overall prevalence and/or abundance distribution. The same annotation can easily be added to metadata-to-microbe correlation heatmaps.

microViz provides a flexible wrapper around methods for the statistical modelling of microbial abundances, including e.g. beta-binomial regression models from the corncob R package [@martin_modeling_2020], and compositional linear regression. To visualize metadata-to-microbiome associations derived from more complicated statistical models, microViz offers a visualization approach that combines multiple annotated cladograms to comprehensively and compactly display patterns of microbial associations with multiple covariates from the same multivariable statistical model. These "taxonomic association trees" facilitate direct comparison of the direction, strength and significance of microbial associations between covariates and across multiple taxonomic ranks. This visualisation also provides an intuitive reminder of the balancing act inherent in compositional data analysis: if one clade/branch goes up, others must go down. Other packages in R, such as ggtree [@yu_ggtree_2017] or metacoder [@foster_metacoder_2017], can be used to make annotated cladograms similar to the microViz taxonomic association tree visualizations, but the microViz style has a few advantages for the purpose of reporting multivariable model results: Firstly, microViz cladogram generating functions are directly paired with functions to compute the statistical model results for all taxa in a phyloseq. Secondly, the tree layouts are more compact, by default, for displaying multiple trees for easy comparison.

Finally, beyond the main visualization functionality, microViz provides a suite of tools for working easily with phyloseq objects including wrapper functions that bring approaches from the popular dplyr package to phyloseq, to help the researcher easily filter, select, join, mutate and arrange phyloseq sample data. All microViz functions are designed to work with magrittr's pipe operator (%\>%), to chain successive functions together and improve code readability [@magrittr]. Lastly, for user convenience, microViz documentation and tutorials are hosted online via a pkgdown [@pkgdown] website on Github Pages, with extensive examples of code and output generated with example datasets.

![Simple example of a microViz figure pairing an ordination plot of microbial samples (left) with an "iris plot" (right): a circular stacked barchart showing the microbial compositions of samples ordered in accordance with the ordination plot. This figure is created with a subset of the "dietswap" dataset available within the microbiome R package. The ordination plot is a PCA bi-plot created using centered-log-ratio transformed species-like HITChip microbial features. The dark grey filled points on both plots indicate samples where the participant's nationality is AFR. AFR = African; AAM = African American. \label{fig:one}](fig1.jpg)

# Acknowledgements

This work was completed as part of a project jointly funded by the Dutch Research Council (NWO), AVEBE, FrieslandCampina and NuScience, as coordinated by the Carbohydrate Competence Center (CCC-CarboBiotics; www.cccresearch.nl).

# References
---
output: github_document
editor_options: 
  markdown: 
    wrap: sentence
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", dev = "CairoPNG",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# microViz <a href='https://david-barnett.github.io/microViz/index.html'><img src="man/figures/logo.png" align="right" height="180" width="156"/></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/david-barnett/microViz/workflows/R-CMD-check/badge.svg)](https://github.com/david-barnett/microViz/actions) [![codecov](https://codecov.io/gh/david-barnett/microViz/branch/main/graph/badge.svg?token=C1EoVkhnxA)](https://codecov.io/gh/david-barnett/microViz) [![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/barnettdavid/microviz-rocker-verse)](https://hub.docker.com/r/barnettdavid/microviz-rocker-verse) [![status](https://joss.theoj.org/papers/4547b492f224a26d96938ada81fee3fa/status.svg)](https://joss.theoj.org/papers/4547b492f224a26d96938ada81fee3fa) [![DOI](https://zenodo.org/badge/307119750.svg)](https://zenodo.org/badge/latestdoi/307119750)

<!-- badges: end -->

## Overview

:package: `microViz` is an R package for analysis and visualization of microbiome sequencing data.

:hammer: `microViz` functions are intended to be easy to use and flexible.

:microscope: `microViz` extends and complements popular microbial ecology packages like `phyloseq`, `vegan`, & `microbiome`.

## Learn more

:paperclip: This website is the best place for documentation and examples: <https://david-barnett.github.io/microViz/>

-   [**This ReadMe**](https://david-barnett.github.io/microViz/) shows a few example analyses

-   **The [Reference](https://david-barnett.github.io/microViz/reference/index.html) page** lists all functions and links to help pages and examples

-   **The [Changelog](https://david-barnett.github.io/microViz/news/index.html)** describes important changes in new microViz package versions

-   **The Articles pages** give tutorials and further examples

    -   [Working with phyloseq objects](https://david-barnett.github.io/microViz/articles/web-only/phyloseq.html)

    -   [Fixing your taxa table with tax_fix](https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html)

    -   [Creating ordination plots](https://david-barnett.github.io/microViz/articles/web-only/ordination.html) (e.g. PCA or PCoA)

    -   [Interactive ordination plots with ord_explore](https://david-barnett.github.io/microViz/articles/web-only/ordination-interactive.html)

    -   [Visualising taxonomic compositions with comp_barplot](https://david-barnett.github.io/microViz/articles/web-only/compositions.html)
    
    -   [Heatmaps of microbiome composition and correlation](https://david-barnett.github.io/microViz/articles/web-only/heatmaps.html)

    -   [Modelling and plotting individual taxon associations with taxatrees](https://david-barnett.github.io/microViz/articles/web-only/modelling-taxa.html)

    -   More coming soon(ish)!
        Post on [GitHub discussions](https://github.com/david-barnett/microViz/discussions) if you have questions/requests

## Installation

You can install the latest available microViz package version using the following instructions.

``` r
# Installing from github requires the devtools package
install.packages("devtools") 

# To install the latest "released" version of this package
devtools::install_github("david-barnett/microViz@0.9.0") # check 0.9.0 is the latest release

# To install the very latest version:
devtools::install_github("david-barnett/microViz")
# If you encounter a bug please try the latest version & let me know if the bug persists!

# If the Bioconductor dependencies don't automatically install you can install
# them yourself like this:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") 
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"))
```

:computer: **Windows users** - will need to have RTools installed so that your computer can build this package (follow instructions here: <http://jtleek.com/modules/01_DataScientistToolbox/02_10_rtools/>)

**:apple: macOS** **users** - might need to install [xquartz](https://www.xquartz.org/) to make the heatmaps work (to do this with homebrew, run the following command in your mac's Terminal: `brew install --cask xquartz`

:package: I highly recommend using [renv](https://rstudio.github.io/renv/index.html) for managing your R package installations across multiple projects.

:whale: Alternatively, for Docker users an image with the main branch installed is available at: <https://hub.docker.com/r/barnettdavid/microviz-rocker-verse>

:date: microViz is tested to work with R version 4 on Windows, MacOS, and Ubuntu 18 and 20.
R version 3.6.\* should probably work, but I don't formally test this.

## Interactive ordination exploration

```{r load, message=FALSE}
library(microViz)
```

microViz provides a Shiny app for an easy way to start exploring your microbiome data: all you need is a phyloseq object.

```{r example ord_explore}
# example data from corncob package
pseq <- corncob::ibd_phylo %>%
  tax_fix() %>%
  phyloseq_validate()
```

``` {.r}
ord_explore(pseq) # gif generated with microViz version 0.7.4 (plays at 1.75x speed)
```

![](vignettes/web-only/images/20210429_ord_explore_x175.gif)

## Example analyses

```{r packages}
library(phyloseq)
library(dplyr)
library(ggplot2)
```

```{r data setup}
# get some example data
data("dietswap", package = "microbiome")

# create a couple of numerical variables to use as constraints or conditions
dietswap <- dietswap %>%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  )
# add a couple of missing values to show how microViz handles missing data
sample_data(dietswap)$african[c(3, 4)] <- NA
```

### Looking at your data

You have quite a few samples in your phyloseq object, and would like to visualise their compositions.
Perhaps these example data differ participant nationality?

```{r fig.height=5, fig.width=6, dpi=120}
dietswap %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 15, other_name = "Other",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    palette = distinct_palette(n = 15, add = "grey90"),
    merge_other = FALSE, bar_outline_colour = "darkgrey"
  ) +
  coord_flip() +
  facet_wrap("nationality", nrow = 1, scales = "free") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
```


```{r, fig.width=5.5, fig.height=3.5, dpi=120}
htmp <- dietswap %>%
  ps_mutate(nationality = as.character(nationality)) %>%
  tax_transform("log2", add = 1, chain = TRUE) %>%
  comp_heatmap(
    taxa = tax_top(dietswap, n = 30), grid_col = NA, name = "Log2p",
    taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
    colors = heat_palette(palette = viridis::turbo(11)),
    row_names_side = "left", row_dend_side = "right", sample_side = "bottom",
    sample_anno = sampleAnnotation(
      Nationality = anno_sample_cat(
        var = "nationality", col = c(AAM = "grey35", AFR = "grey85"),
        box_col = NA, legend_title = "Nationality", size = grid::unit(4, "mm")
      )
    )
  )

ComplexHeatmap::draw(
  object = htmp, annotation_legend_list = attr(htmp, "AnnoLegends"),
  merge_legends = TRUE
)
```


### Example ordination plot workflow

Ordination methods can also help you to visualise if overall microbial ecosystem composition differs markedly between groups, e.g. BMI.

Here is one option as an example:

1.  Filter out rare taxa (e.g. remove Genera not present in at least 10% of samples) - use `tax_filter()`
2.  Aggregate the taxa into bacterial families (for example) - use `tax_agg()`
3.  Transform the microbial data with the centre-log-ratio transformation - use `tax_transform()`
4.  Perform PCA with the clr-transformed features (equivalent to aitchison distance PCoA) - use `ord_calc()`
5.  Plot the first 2 axes of this PCA ordination, colouring samples by group and adding taxon loading arrows to visualise which taxa generally differ across your samples - use `ord_plot()`
6.  Customise the theme of the ggplot as you like and/or add features like ellipses or annotations

```{r ordination-plot, dpi=120}
# perform ordination
unconstrained_aitchison_pca <- dietswap %>%
  tax_filter(min_prevalence = 0.1, tax_level = "Genus") %>%
  tax_agg("Family") %>%
  tax_transform("clr") %>%
  ord_calc()
# ord_calc will automatically infer you want a "PCA" here
# specify explicitly with method = "PCA", or you can pick another method

# create plot
pca_plot <- unconstrained_aitchison_pca %>%
  ord_plot(
    plot_taxa = 1:6, colour = "bmi_group", size = 1.5,
    tax_vec_length = 0.325,
    tax_lab_style = tax_lab_style(max_angle = 90, aspect_ratio = 0.5),
    auto_caption = 8
  )

# customise plot
customised_plot <- pca_plot +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group), size = 0.3) +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "bottom") +
  coord_fixed(ratio = 0.5, clip = "off") # makes rotated labels align correctly

# show plot
customised_plot
```

### PERMANOVA

You visualised your ordinated data in the plot above.
Below you can see how to perform a PERMANOVA to test the significance of BMI's association with overall microbial composition.
This example uses the Family-level aitchison distance to correspond with the plot above.

```{r permanova}
# calculate distances
aitchison_dists <- dietswap %>%
  tax_filter(min_prevalence = 0.1) %>%
  tax_transform("identity", rank = "Family") %>%
  dist_calc("aitchison")

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_perm <- aitchison_dists %>% 
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 99, # you should use at least 999!
    variables = "bmi_group"
  )

# view the permanova results
perm_get(aitchison_perm) %>% as.data.frame()

# view the info stored about the distance calculation
info_get(aitchison_perm)
```

### Constrained partial ordination

You could visualise the effect of the (numeric/logical) variables in your permanova directly using the `ord_plot` function with constraints (and conditions). 

```{r constrained-ord}
perm2 <- aitchison_dists %>% 
  dist_permanova(variables = c("weight", "african", "sex"), seed = 321)
```

We'll visualise the effect of nationality and bodyweight on sample composition, after first removing the effect of sex.

```{r constrained-ord-plot, fig.height=6, fig.width=6, dpi=120}
perm2 %>%
  ord_calc(constraints = c("weight", "african"), conditions = "female") %>%
  ord_plot(
    colour = "nationality", size = 2.5, alpha = 0.35,
    auto_caption = 7,
    constraint_vec_length = 1,
    constraint_vec_style = vec_constraint(size = 1.5, colour = "grey15"),
    constraint_lab_style = constraint_lab_style(
      max_angle = 90, size = 3, aspect_ratio = 0.35, colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = nationality), size = 0.2) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed(ratio = 0.35, clip = "off") +
  theme(legend.position = c(0.9, 0.1), legend.background = element_rect())
```

### Correlation Heatmaps

microViz heatmaps are powered by `ComplexHeatmap` and annotated with taxa prevalence and/or abundance.

```{r heatmap, dpi=120, fig.width=7, fig.height=5.5}
# set up the data with numerical variables and filter to top taxa
psq <- dietswap %>%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  ) %>%
  tax_filter(
    tax_level = "Genus", min_prevalence = 1 / 10, min_sample_abundance = 1 / 10
  ) %>%
  tax_transform("identity", rank = "Genus")

# randomly select 30 taxa from the 50 most abundant taxa (just for an example)
set.seed(123)
taxa <- sample(tax_top(psq, n = 50), size = 30)
# actually draw the heatmap
cor_heatmap(
  data = psq, taxa = taxa,
  taxon_renamer = function(x) stringr::str_remove(x, " [ae]t rel."),
  tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(undetected = 50),
    Log2 = anno_tax_box(undetected = 50, trans = "log2", zero_replace = 1)
  )
)
```

## Citation

:innocent: If you find any part of microViz useful to your work, please consider citing the JOSS article:

Barnett et al., (2021).
microViz: an R package for microbiome data visualization and statistics.
Journal of Open Source Software, 6(63), 3201, <https://doi.org/10.21105/joss.03201>

## Contributing

Bug reports, questions, suggestions for new features, and other contributions are all welcome.
Feel free to create a [GitHub Issue](https://github.com/david-barnett/microViz/issues) or write on the [Discussions](https://github.com/david-barnett/microViz/discussions) page.
Alternatively you could also contact me (David) on Twitter [\@\_david_barnett\_](https://twitter.com/_david_barnett_) .

This project is released with a [Contributor Code of Conduct](https://david-barnett.github.io/microViz/CODE_OF_CONDUCT.html) and by participating in this project you agree to abide by its terms.

## Session info

```{r session}
sessionInfo()
```
---
title: "This page has moved"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


See the homepage Readme here: https://david-barnett.github.io/microViz/ 

The brief example analyses with the atlas1006 dataset are now here:
https://david-barnett.github.io/microViz/articles/web-only/atlas1006.html

More detailed tutorials on various microViz features are listed here:
https://david-barnett.github.io/microViz/articles/


---
title: "Statistical modelling of individual taxa"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article will give an example of statistical modelling of the
abundances of individual taxa, and visual presentation of the results
using microViz *taxonomic association tree* plots.

## Setup

```{r setup}
library(microViz)
library(corncob)
library(dplyr)
library(ggplot2)
```

First we'll get some OTU abundance data from inflammatory bowel disease
patients and controls from the corncob package.

```{r}
data("ibd_phylo")
ibd_phylo
```

We'll keep only the Ulcerative Colitis and Healthy Control samples, to
simplify the analyses for this example. We'll also remove the Species
rank information, as most OTUs in this dataset are not assigned to a
species. We'll also use `tax_fix` to fill any gaps where the Genus is
unknown, with the family name or whatever higher rank classification is
known.

```{r}
phylo <- ibd_phylo %>% 
  ps_filter(DiseaseState %in% c("UC", "nonIBD")) %>% 
  tax_mutate(Species = NULL) %>% 
  tax_fix()
phylo 
```

Let's have a quick look at the sample data using the `skimr` package.

```{r}
phylo %>% 
  samdat_tbl() %>% 
  dplyr::mutate(across(where(is.character), as.factor)) %>% 
  skimr::skim()
```

Let's make some sample data variables that are easier to use and compare
in the statistical modelling ahead. We will convert dichotomous
categorical variables into similar binary variables (values: 1 for true,
or 0 for false). We will also scale and center the numeric variable for
age.

```{r}
phylo <- phylo %>% 
  ps_mutate(
    UC = ifelse(DiseaseState == "UC", yes = 1, no = 0),
    female = ifelse(gender == "female", yes = 1, no = 0),
    antibiotics = ifelse(abx == "abx", yes = 1, no = 0),
    steroids = ifelse(steroids == "steroids", yes = 1, no = 0),
    age_scaled = scale(age, center = TRUE, scale = TRUE)
  )
```

## Modelling one taxon at a time

### TSS log2 linear regression

We will start by creating a linear regression model for one genus,
*Bacteroides*. We will transform the count data by first making it
proportions, and then taking the logarithm, with base 2. This is what
MaAsLin2 does by default (except they call the compositional
transformation "Total Sum Scaling (TSS)"). We will replace the zeros
with half the minimum observed abundance proportion (of any taxon)
before log2 transformation.

```{r}
parabacteroides_lm <- phylo %>% 
  tax_fix() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>% 
  tax_model(
    type = "lm", rank = "Genus", taxa = "Parabacteroides", 
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  ) 
parabacteroides_lm$Parabacteroides
```

```{r}
summary(parabacteroides_lm$Parabacteroides)
```

This model suggests *Parabacteroides* abundances are significantly lower
in Ulcerative Colitis patients than controls, on average.

### Plotting TSS log2 data

Let's boxplot the transformed data to see what this *Parabacteroides*
association looks like as a crude group difference.

```{r, fig.width=4, dpi=120}
plot_data <- phylo %>% 
  tax_fix() %>% 
  tax_transform("compositional", rank = "Genus") %>% 
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>% 
  ps_get() %>% 
  ps_otu2samdat("Parabacteroides") %>% # adds Parabacteroides as sample data!
  samdat_tbl()

ggplot(plot_data, aes(x = DiseaseState, y = Parabacteroides)) + 
  geom_boxplot(width = 0.5, colour = "grey35") +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  scale_y_continuous(
    breaks = log2(1/2^(0:13)), 
    labels = function(x) paste0(100 * round(2^x, digits = 5), "%"), 
    limits = c(log2(0.00005), log2(0.25))
  ) + 
  theme_bw()
```

### Beta binomial regression

You can also use other regression modelling functions that take a
formula. For example the beta binomial modelling provided in the corncob
package. This approach models both abundance and dispersion, and
directly uses untransformed counts. By default, microViz's `tax_model()`
will use the same formula for both abundance and dispersion modelling,
but you can override this by setting the `phi.formula` argument
yourself. See `vignette("corncob-intro", package = "corncob")` for more
info on these models.

```{r}
parabacteroides_bb <- phylo %>% 
  tax_fix() %>% 
  tax_model(
    type = corncob::bbdml, rank = "Genus", taxa = "Parabacteroides", 
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  ) 
parabacteroides_bb$Parabacteroides
```

## Model all the taxa!

Now we will fit a similar model for almost every taxon at every rank.
The code for `taxatree_models` is quite similar to `tax_model`. However,
you might need to run `tax_prepend_ranks` to ensure that each taxon at
each rank is always unique. As an example of the problem, Actinobacteria
is the name of both a Phylum and a Class!

```{r}
lm_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  # it makes sense to perform the compositional transformation BEFORE filtering
  tax_transform("compositional", rank = "Genus", keep_counts = TRUE) %>% 
  tax_filter(min_prevalence = 0.1, undetected = 0, use_counts = TRUE) %>% 
  tax_transform(
    trans = "log2", chain = TRUE, zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = NULL, # uses every rank available except the first
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  )
```

Why filter the taxa? *It's less likely that we are interested in rare
taxa, and models of rare taxon abundances are more likely to be
unreliable. Reducing the the number of taxa modelled also makes
visualising the results easier!*

```{r}
lm_models
```

### Getting stats from the models

Next we will get a data.frame containing the regression coefficient
estimates, test statistics and corresponding p values from all these
regression models. The function `taxatree_models2stats()` can do this
for any type of model that has a `broom::tidy()` method, as well as for
beta binomial regression models calculated with the `corncob` package
`bbdml()` function.

```{r}
lm_stats <- taxatree_models2stats(lm_models)
lm_stats
```

```{r}
lm_stats$taxatree_stats
```

### Adjusting p values

Using the `taxatree_stats_p_adjust()` function, you can correct for
multiple testing / control the false discovery rate or family-wise error
rate.

Instead of applying these adjustment methods across all 88 taxa models
at all ranks, the default behaviour is to control the family-wise error
rate per rank.

```{r}
lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)
# notice the new variable
lm_stats$taxatree_stats
```

## Plot all the taxatree_stats!

`taxatree_plots()` allows you to plot statistics (e.g. point estimates
and significance) from all of the taxa models onto a tree layout. The
taxon models are organised by rank, radiating out from the central root
node from e.g. Phyla around the center to Genera in the outermost ring.

`taxatree_plots()` itself returns a list of plots, which you can arrange
into one figure with the
[`patchwork`](https://patchwork.data-imaginist.com/) package for example
(and/or
[`cowplot`](https://wilkelab.org/cowplot/articles/plot_grid.html)).

```{r, fig.width=7, fig.height=10}
lm_stats %>% taxatree_plots(
  node_size_range = c(1, 3), var_renamer = toupper
) %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )
```

### Taxatree Key

But how do we know which taxa are which nodes? We can create a labelled
grey tree with `taxatree_plotkey`. This labels taxa based on certain
conditions.

```{r, fig.height=5, fig.width=6}
set.seed(123) # label position 
key <- taxatree_plotkey(
  data = lm_stats, 
  taxon_renamer = function(x) stringr::str_remove(x, "[PFG]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Phylum" | rank == "Genus" & prevalence > 0.25, 
  !grepl("Kingdom", taxon)
) + 
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))
# all phyla are labelled, and all genera with a prevalence of over 0.2
# except for any taxa whose names (partly) match "Kingdom" 
# (i.e. an unclassified taxon)
key
```

#### Key + Trees

Let's put the key and some of the trees together in one `patchwork`
figure. Getting the sizing right on these combined plots can be very
tricky! Pay attention to the absolute height and width of the plot
output.

`gridExtra::grid.arrange()` or `cowplot::plot_grid()` are alternatives
you can also try. `cowplot::get_legend()` can be particularly useful.

```{r, fig.width=13, fig.height=5.5, dpi=120}
trees <- lm_stats %>% 
  taxatree_plots(node_size_range = c(1, 2.25)) %>%
  .[1:4] %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )

panel <- patchwork::wrap_plots(key, trees, nrow = 1, widths = 8:7)
set.seed(111)
panel
```

You could save the plot with `ggsave()` like this.

```{r eval=FALSE}
set.seed(111)
ggsave("test.png", panel, width = 13, height = 5.5, dpi = 120, device = "png")
```

#### Alternative label styling

You can change the default styling of the labels by first suppressing
the automatic drawing of labels with `.draw_label = FALSE` in
`taxatree_plotkey()` and then adding your own custom-style labels with
`taxatree_plot_labels()`. Here we will draw some yellow labels.

```{r, fig.height=5, fig.width=7}
taxatree_plotkey(
  data = lm_stats, .draw_label = FALSE,
  rank %in% c("Phylum", "Family") & !grepl("Bacteria", taxon),
  prevalence > 0.2 | rank == "Phylum"
) %>% 
  taxatree_plot_labels(
    taxon_renamer = function(x) stringr::str_remove(x, "[PFGO]: "),
    # default fun is ggrepel::geom_text_repel
    fun = ggrepel::geom_label_repel, 
    # arguments modifying label style
    size = 2.5, alpha = 0.5, colour = "black", fill = "gold1", 
    label.size = 0.05, label.r = unit(0.05, "lines"), 
    label.padding = unit(0.15, "lines"), segment.size = 0.5,
    # arguments affecting label positioning
    box.padding = 0.05, x_nudge = 0.4, y_nudge = 0.05,
    hjust = 0.5, seed = 123
  )
```

### Directly labelling taxa

You can directly label taxatree_plots too, but it is better to only do
this for a few taxa. You must run `taxatree_label()` first to create a
"label" indicator variable.

```{r, fig.width=6.5, fig.height=4}
lm_stats %>% 
  taxatree_label(
    rank == "Genus", p.value < 0.05 | prevalence > 0.5, estimate > 0
  ) %>% 
  taxatree_plots() %>% 
  .[[1]] %>% # show only the first plot
  taxatree_plot_labels(
    taxon_renamer = function(x) stringr::str_remove(x, "G: "),
    fun = ggrepel::geom_label_repel, x_nudge = 0.7, hjust = 0.5, size = 2
  ) +
  # the lines below allow expansion and make some space for labels
  coord_fixed(clip = "off", expand = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.3))
```

### Changing color palette

Choosing another color palette is easy, just name any diverging palette
from colorspace hcl diverging palettes. See your options below.

```{r, fig.width=3, fig.height=2.5}
colorspace::hcl_palettes(type = "diverging", plot = TRUE, n = 11)
```

By default, the colour scale is transformed by taking the square root of
the absolute values. But you can change this to "identity" to have no
palette transformation.

By default the range of the data is used to set symmetrical limits on
the colour scale, which are the same for all plots in the list produced.
You can set alternative limits. If some data lie outside these limits,
their values will be "squished" to match the max or min limit value.

For finer control of the palette luminance range, you can set custom
values for l1 and l2, e.g. if the extremes are too bright or dark. This
is done by default for the Green-Brown palette.

```{r, fig.width=6.5, fig.height=4}
lm_stats %>% 
  taxatree_label(
    rank == "Genus", p.value < 0.05 | prevalence > 0.5, estimate > 0
  ) %>% 
  taxatree_plots(
    colour_lims = c(-20, 20), colour_trans = "identity",
    palette = "Blue-Red 3", l2 = 90
  ) %>% 
  .[[1]] %>% # show only the first plot
  taxatree_plot_labels(
    taxon_renamer = function(x) stringr::str_remove(x, "G: "),
    fun = ggrepel::geom_label_repel, x_nudge = 0.7, hjust = 0.5, size = 2
  ) +
  # the lines below allow expansion and make some space for labels
  coord_fixed(clip = "off", expand = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.3)) 
```

Palettes like "Berlin" that go through a black midpoint would probably
only make sense with a darker background!

```{r, fig.width=6.5, fig.height=4}
lm_stats %>% 
  taxatree_label(
    rank == "Genus", p.value < 0.05 | prevalence > 0.5, estimate > 0
  ) %>% 
  taxatree_plots(
    palette = "Berlin", colour_lims = c(-20, 20)
  ) %>% 
  .[[1]] %>% # show only the first plot
  taxatree_plot_labels(
    taxon_renamer = function(x) stringr::str_remove(x, "G: "),
    fun = ggrepel::geom_label_repel, 
    x_nudge = 0.7, xlim = c(-1.7, 1.5), hjust = 0.5, size = 2
  ) +
  # the lines below allow expansion and make some space for labels
  coord_fixed(clip = "off", expand = TRUE) +
  scale_x_continuous(expand = expansion(mult = 0.3)) +
  theme(
    text = element_text(colour = "white"),
    plot.background = element_rect(fill = "grey30"), 
    plot.title = element_text(size = 20, colour = "white")
  )
```

### Sorting taxa nodes

If you like, you can sort the nodes by sorting the taxa in the ps_extra
object.

```{r, fig.width=5, fig.height=4}
lm_stats %>% 
  tax_sort(by = "prev", at = "Genus") %>%
  taxatree_plots() %>% 
  .[[1]] # show only the first plot
```

You can chain multiple `tax_sort()` calls together to fine-tune the
order of the nodes on the tree to your own preference.

```{r, fig.width=5, fig.height=4}
lm_stats %>% 
  tax_sort(by = "prev", at = "Family") %>%
  tax_sort(by = "name", at = "Phylum") %>%
  tax_sort(by = "rev") %>%
  taxatree_plots() %>% 
  .[[1]] # show only the first plot
```

### Plotting adjusted p values

Remember we made adjusted p values earlier? Let's plot those instead.
Just to show how it's done, we'll also change the symbol used to
identify the significant sites to a cross, and we'll also relax the
significance threshold to 0.1.

It looks like only the disease state (having ulcerative colitis) shows
any significant associations after this FDR correction.

```{r, fig.width=7, fig.height=6}
lm_stats %>% 
  taxatree_plots(
    sig_stat = "p.adj.BH.rank", sig_threshold = 0.1, 
    sig_shape = "cross", sig_size = 1.5,
    node_size_range = c(1, 3), var_renamer = toupper
) %>% 
  .[1:4] %>% # keep only first 4 plots
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )
```

### Plotting multiple significance markers

You can also plot multiple significance markers. You must start with the
strictest threshold. Here we will plot the FDR corrected significance
markers at for p.adj \< 0.05 (as thick white crosses) and then also
unadjusted significance markers for p \< 0.05 (as outlined white
circles).

```{r, fig.width=7, fig.height=6}
lm_stats %>% 
  taxatree_plots(
    sig_stat = c("p.adj.BH.rank", "p.value"), sig_threshold = 0.05, 
    sig_shape = c("cross", "circle filled"), sig_colour = "white",
    sig_size = c(1.5, 1), sig_stroke = c(1, 0.25),
    node_size_range = c(1, 3), var_renamer = toupper
) %>% 
  .[1:4] %>% # keep only first 4 plots
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )
```

## Beta-binomial regression example

The corncob package provides beta-binomial regression models. See the
paper [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7514055/), and
the helpful package vignette:
`vignette("corncob-intro", package = "corncob")`.

We will filter the taxa more strictly (by a higher prevalence threshold)
before this type of modelling. We do not need to transform the data, as
this approach uses counts.

```{r}
bb_models <- phylo %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  tax_filter(min_prevalence = 0.3) %>% 
  taxatree_models(
    type = corncob::bbdml, 
    ranks = c("Phylum", "Class", "Order", "Family"),
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  )
bb_models
```

When extracting stats from corncob beta-binomial models, you need to
specify which parameter estimate you want, "mu" for differential
abundance, or "phi" for differential variability or overdispersion.

```{r}
bb_stats <- taxatree_models2stats(bb_models, param = "mu")
bb_stats
bb_stats$taxatree_stats
```

```{r, fig.width=5, fig.height=5}
bb_stats %>% taxatree_plots(
  node_size_range = c(1, 4), colour_trans = "identity"
) %>% 
  # keep only first 4 plots
  .[1:4] %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
```

### Alternative layouts

You do not need to make circular tree plots if you don't want to!

```{r, fig.width=5, fig.height=5}
alt_trees <- bb_stats %>% 
  taxatree_plots(
    node_size_range = c(1, 4), circular = FALSE, colour_trans = "identity"
  ) %>% 
  # keep only first 4 plots
  .[1:4] %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) & # & is used by patchwork to modify multiple ggplots (instead of +)
  coord_flip() & 
  scale_y_reverse()

alt_trees
```

Let's add the key for this layout and label it manually with
`taxatree_plot_labels()`.

```{r, fig.width=5, fig.height=8}
alt_tree_key <- taxatree_plotkey(
  data = bb_stats, circular = FALSE, .draw_label = FALSE,
  rank == "Family"
) %>% 
  taxatree_plot_labels(
    circular = FALSE, 
    taxon_renamer = function(x) stringr::str_remove(x, "[PFG]: "),
    hjust = 0.5, force = 0, nudge_y = 2, size = 3
  )  + 
  coord_flip() + 
  scale_y_reverse(expand = expansion(mult = c(0.1, 1))) +
  theme(plot.title = element_text(size = 14))

patchwork::wrap_plots(alt_tree_key, alt_trees, nrow = 2, heights = 1:2)
```

You don't have to use a regular tree!

Alternative layouts from the igraph package are possible, such as the
Kamada and Kawai spring algorithm ("kk") or Fruchterman and Reingold
force-directed algorithm ("fr"). You must set a layout_seed number for
these layouts to ensure they are always the same.

```{r, fig.width=5, fig.height=5}
bb_stats %>% taxatree_plots(
  node_size_range = c(1, 4), 
  colour_trans = "identity", layout = "kk", layout_seed = 321
) %>% 
  # keep only first 4 plots
  .[1:4] %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )
```

# Session info

```{r}
devtools::session_info()
```
---
title: "microViz annotated heatmaps"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article will show you how to plot annotated correlation and
microbial composition heatmaps with microViz.

## Setup

```{r setup}
library(dplyr)
library(phyloseq)
library(microViz)
```

First we'll get some OTU abundance data from inflammatory bowel disease
patients and controls from the corncob package.

```{r}
data("ibd_phylo", package = "corncob")
ibd_phylo
```

Remove the mostly unclassified species-level data, drop the rare taxa
and fix the taxonomy of the rest. Also drop patients with unclassified
IBD.

```{r}
psq <- ibd_phylo %>% 
  tax_mutate(Species = NULL) %>% 
  tax_filter(min_prevalence = 5) %>% 
  tax_fix() %>% 
  ps_filter(DiseaseState != "IBDundef")
psq
```

## Microbiome heatmaps

Visualise the microbial composition of your samples.

The samples and taxa are sorted by similarity. (*By default this uses
hierarchical clustering with optimal leaf ordering, using euclidean
distances on the transformed data*).

In this example we use a "compositional" transformation, so the Class
abundances are shown as proportions of each sample.

```{r, fig.height=2.5, fig.width=7, dpi=120}
psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap()
```

You can easily swap to a symmetrical colour palette for transformations
like "clr" or "standardize". This is the default symmetrical palette but
you can pick from many.

```{r, fig.height=4, fig.width=7, dpi=120}
psq %>% 
  tax_transform("clr", rank = "Family") %>% 
  comp_heatmap(colors = heat_palette(sym = TRUE), name = "CLR")
```

### Annotating taxa

```{r, fig.height=2.5, fig.width=7, dpi=120}
psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap(tax_anno = taxAnnotation(
    Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
  ))
```

### Legend positioning

Positioning the heatmap legend at the bottom is possible. You can assign
the heatmap to a name and then call `ComplexHeatmap`'s `draw` function.

```{r, fig.height=2.5, fig.width=7, dpi=120}
heat <- psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))),
    heatmap_legend_param = list(
      at = 0:5 / 5,
      direction = "horizontal", title_position = "leftcenter", 
      legend_width = grid::unit(4, "cm"), grid_height = grid::unit(5, "mm")
      )
  )

ComplexHeatmap::draw(
  object = heat, heatmap_legend_side = "bottom", 
  adjust_annotation_extension = FALSE
  )
```


### Annotating samples

#### Group membership

2 different methods for annotating each sample's values of categorical
metadata are possible.

-   `anno_sample()` cannot have borders around each cell, but
    automatically adds a legend.

-   `anno_sample_cat()` can have cell borders, but requires an extra
    step to draw a legend

```{r, fig.height=3.5, fig.width=7, dpi=120}
cols <- distinct_palette(n = 3, add = NA)
names(cols) <- unique(samdat_tbl(psq)$DiseaseState)

psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    sample_anno = sampleAnnotation(
      State1 = anno_sample("DiseaseState"),
      col = list(State1 = cols), border = FALSE,
      State2 = anno_sample_cat("DiseaseState", col = cols)
    )
  ) 
```

Let's try drawing equivalent categorical annotations by two methods.
Both methods can draw annotations with borders and no individual boxes.
This style suits heatmaps with no gridlines (i.e. `grid_col = NA`).

In the example below we have suppressed row ordering with
`cluster_rows = FALSE`, and added spaces between taxa by splitting at
every row with `row_split = 1:11`, which are both
`ComplexHeatmap::Heatmap()` arguments.

```{r, fig.height=4, fig.width=7, dpi=120}
psqC <- psq %>% tax_transform("compositional", rank = "Class") 

htmp <- psqC %>% 
  comp_heatmap(
    grid_col = NA,
    cluster_rows = FALSE, row_title = NULL,
    row_split = seq_len(ntaxa(ps_get(psqC))),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.9, size = grid::unit(1, "cm"), border = F)
    ),
    sample_anno = sampleAnnotation(
      # method one
      State1 = anno_sample("DiseaseState"), 
      col = list(State1 = cols), border = TRUE,
      # method two
      State2 = anno_sample_cat(
        var = "DiseaseState", col = cols, box_col = NA, border_col = "black",
        legend_title = "State2"
      )
    )
  ) 
htmp %>% ComplexHeatmap::draw(
  annotation_legend_list = attr(htmp, "AnnoLegends")
)
```

You can also manually draw a legend with the convenience function 
`anno_cat_legend()`. 


```{r, fig.height=0.5, fig.width=4, dpi=120}
grid::grid.newpage()
anno_cat_legend(
  col = c("a level" = "red", "another level" = "blue", c = "white"), 
  border = "black", gap = grid::unit(2, "cm"), ncol = 3
) %>% 
ComplexHeatmap::draw()
```


### Numbering cells

If you have fewer samples (and taxa) you might like to label the cells
with their values. By default, the raw counts are shown.

```{r, fig.height=3, fig.width=7, dpi=120}
psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap(samples = 1:15, numbers = heat_numbers())
```

You can easily change to showing the same values as the colours by
setting `numbers_use_counts = FALSE`, and you can/should change the
number of decimals shown too.

```{r, fig.height=3, fig.width=7, dpi=120}
psq %>% 
  tax_transform("compositional", rank = "Class") %>% 
  comp_heatmap(
    samples = 1:15, numbers_use_counts = FALSE,
    numbers = heat_numbers(decimals = 2)
  )
```

The numbers can any transformation of counts, irrespective of what
transformations were used for the colours, or seriation.

```{r, fig.height=3, fig.width=7, dpi=120}
psq %>% 
  tax_transform("binary", undetected = 0, rank = "Class") %>% 
  comp_heatmap(
    samples = 1:15, numbers_use_counts = TRUE, numbers_trans = "compositional", 
    numbers = heat_numbers(decimals = 2, col = "white"), 
    show_heatmap_legend = FALSE
  )
```

To demonstrate that coloration, numbering and seriation can all use
different transformations of the original count data, the example below
we specifies seriating the taxa and samples using the same numerical
values used for the numbers transformation, not the colours, which are
just presence/absence!

```{r, fig.height=3, fig.width=7, dpi=120}
psq %>% 
  tax_transform("binary", undetected = 0, rank = "Class") %>% 
  comp_heatmap(
    samples = 1:15, 
    sample_ser_counts = TRUE, sample_ser_trans = "compositional",
    tax_ser_counts = TRUE, tax_ser_trans = "compositional",
    numbers_use_counts = TRUE, numbers_trans = "compositional", 
    numbers = heat_numbers(decimals = 2, col = "white"), 
    show_heatmap_legend = FALSE
  )
```

## Correlation heatmaps

Correlation heatmaps can be a nice way to quickly assess patterns of
associations between numerical variables in your dataset, such as
microbial abundances and other metadata.

Let's make some fake numeric variables to exemplify this.

```{r}
set.seed(111) # ensures making same random variables every time!
psq <- psq %>% 
  ps_arrange(ibd) %>% 
  ps_mutate(
    var1 = rnorm(nsamples(psq), mean = 10, sd = 3),
    var2 = c(
      rnorm(nsamples(psq)*0.75, mean = 4, sd = 2),
      rnorm(1 + nsamples(psq)/4, mean = 9, sd = 3)
      ),
    var3 = runif(nsamples(psq), 2, 10),
    var4 = rnorm(nsamples(psq), mean = 100 + nsamples(psq):0, sd = 20) / 20,
    var5 = rnorm(nsamples(psq), mean = 5, sd = 2),
    var6 = rnbinom(nsamples(psq), size = 1:75/10, mu = 5)
  )
```

### Calculating correlations

By default, the `cor_heatmap` function will correlate all taxa to all
numerical sample data, using pearson correlation method.

```{r, fig.width=7, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(vars = c("var1", "var2", "var3", "var4", "var5", "var6"))
```

It's easy to change to a different method, i.e. spearman's rank
correlation or kendall's tau, which will be reflected in the legend
title. We will also specify to use only the 15 most abundant taxa, by
maximum count, just to make these tutorial figures a little more
compact!

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), cor = "spearman"
  )
```

*Older versions of microViz `cor_heatmap` had a `tax_transform`
argument. But for flexibility, you must now transform your taxa
**before** passing the ps_extra object to `cor_heatmap`.*

Here we have transformed our taxa with the "clr" or centered-log-ratio
transformation prior to correlating. *Notice that the annotations stay
on the same scale, as by default the annotation functions extract the
stored counts data from the ps_extra input, not the transformed data.*

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  tax_transform("clr", zero_replace = "halfmin") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6)
  )
```

Let's transform and scale the taxon abundances before correlating.

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  tax_transform("clr", zero_replace = "halfmin") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6)
  )
```

### Taxon annotations

As seen in the previous plots taxa are annotated by default with
prevalence and relative abundance.

You can transform the taxa for the abundance annotation. The `trans` and
`zero_replace` arguments are sent to `tax_transform()`.

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), 
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1), 
      Log10. = anno_tax_box(trans = "log10", zero_replace = "halfmin")
    )
  )
```

You can do multiple transformations and or scaling by supplying a
function, that takes a ps_extra or phyloseq object, transforms it, and
returns it.

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), 
    tax_anno = taxAnnotation(
      Log2 = anno_tax_density(
        joyplot_scale = 2, gp = grid::gpar(fill = "black", alpha = 0.2),
        trans = "log2", zero_replace = 1
      ),
      `prop Log2` = anno_tax_density(
        joyplot_scale = 1.5, gp = grid::gpar(fill = "black", alpha = 0.2),
        trans = function(ps) {
          ps %>% 
            tax_transform("compositional", zero_replace = 1) %>% 
            tax_transform("log2", chain = TRUE) 
        }
      )
    )
  )
```

Note that by default the relative abundance is shown only for samples
where the taxon is detected! You can include values for all samples for
all taxa by setting `only_detected = FALSE`.

Let's try this with a heatmap-style density plot annotation. We'll
replace zeroes with ones for an interpretable minimum value on the plot.

We'll compare it side-by-side with the default setting of showing only
distribution of values above the detection threshold.

For zero-inflated microbiome data, showing prevalence and "abundance
when detected" often seems like a more informative annotation.

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), 
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(size = grid::unit(10, "mm"), ylim = 0:1), 
      All = anno_tax_density(
        size = grid::unit(20, "mm"),
        trans = "log10", zero_replace = 1, 
        heatmap_colors = viridisLite::turbo(n = 15),
        type = "heatmap", only_detected = FALSE
      ),
      Default = anno_tax_density(
        size = grid::unit(20, "mm"),
        trans = "log10", zero_replace = 1, 
        heatmap_colors = viridisLite::turbo(n = 15),
        type = "heatmap", only_detected = TRUE
      )
    )
  )
```

#### Sorting

By default, rows and columns are sorted using hierarchical clustering
with optimal leaf ordering `"OLO_ward"`. You can use any valid method
from the `seriation` package. You can suppress ordering by using
`seriation_method = "Identity"`. By default this also suppresses column
ordering, so you can set `seriation_method_col = OLO_ward` to keep
ordering.

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  tax_sort(by = prev, at = "Family") %>% 
  cor_heatmap(
    seriation_method = "Identity",
    seriation_method_col = "OLO_ward",
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), 
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1), 
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )
```

#### Taxa annotation side

You can easily put the taxa annotations on another of the heatmap with
e.g. `taxa_side = "left"`

```{r, fig.width=7, fig.height=3, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  tax_sort(by = prev, at = "Family") %>% 
  cor_heatmap(
    seriation_method = "Identity",
    seriation_method_col = "OLO_ward",
    taxa_side = "left",
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1), 
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )
```

Or on the top or bottom is also possible, this will rotate the heatmap.
Remember to swap the seriation method arguments around!

```{r, fig.width=5, fig.height=4, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  tax_sort(by = prev, at = "Family") %>% 
  cor_heatmap(
    seriation_method_col = "Identity", #swapped!
    seriation_method = "OLO_ward", #swapped!
    taxa_side = "top",
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1), 
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )
```

### Variable annotation

As well as annotating the taxa, you can also annotate the variables.

```{r, fig.width=6, fig.height=6, dpi=120}
psq %>% 
  tax_agg("Family") %>% 
  cor_heatmap(
    taxa = tax_top(psq, 15, by = max, rank = "Family"),
    vars = paste0("var", 1:6), 
    var_anno = varAnnotation(
      Value = anno_var_box(),
      Zscore = anno_var_density(fun = scale, type = "violin")
    )
  )
```

## Other stuff

Complicated stuff demonstrated down here, not necessarily useful.

### Custom breaks and seriation

Two approaches to custom colour scale breaks. The first way is better, because
the colour scale is interpolated through the default 11 colours, instead of only 5.

Transform data and customise only labels.

```{r, fig.height=2.5, fig.width=8, dpi=120}
psq %>%
  tax_transform("compositional", rank = "Class") %>%
  tax_transform("log10", zero_replace = "halfmin", chain = TRUE) %>% 
  comp_heatmap(
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    heatmap_legend_param = list(
      labels = rev(c("100%", " 10%", "  1%", " 0.1%", "0.01%"))
    )
  )
```


This alternative way might be helpful in some cases, maybe... 
It demonstrates that custom breaks can be set in `heat_palette()`.

```{r, fig.height=2.5, fig.width=8, dpi=120}
# seriation transform
serTrans <- function(x) {
  tax_transform(x, trans = "log10", zero_replace = "halfmin", chain = TRUE)
}

psq %>%
  tax_transform("compositional", rank = "Class") %>%
  comp_heatmap(
    sample_ser_trans = serTrans, tax_ser_trans = serTrans,
    colors = heat_palette(breaks = c(0.0001, 0.001, 0.01, 0.1, 1), rev = T),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(bar_width = 0.3, size = grid::unit(1, "cm"))
    ),
    heatmap_legend_param = list(at = c(0.0001, 0.001, 0.01, 0.1, 1), break_dist = 1)
  )
```

## Session info

```{r}
devtools::session_info()
```

---
title: "Interactive Ordination Plot"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This article shows you how to quickly get started with interactive exploration of your data/ ordination plot.

```{r setup}
library(phyloseq)
library(microViz)
```

## Example

Example data loaded from the corncob package. All you need is a valid phyloseq object, and to run `tax_fix` to [ensure the tax_table doesn't contain problematic names](https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html).

```{r example}
pseq <- corncob::ibd_phylo %>% tax_fix() %>% phyloseq_validate()
```

The gif animation below shows the result of running `ord_explore`, the animation starts immediately after interactively selecting "Genus" level aggregation, "clr" transformation, and the "PCA" ordination method from the "Edit" menu.

``` {.R}
ord_explore(pseq)
```

![](images/20210429_ord_explore_x175.gif)

## Another (old) example

Get example dataset from the phyloseq package and clean up the taxa just a little.

```{r get data}
data("enterotype", package = "phyloseq")
taxa_names(enterotype)[1] <- "Unclassified" # replace strange "-1" name
ps <- tax_fix(enterotype) # remove NA taxa and similar
```

Create simple Bray-Curtis PCoA to explore interactively.

```{r create ordination}
ord1 <- ps %>%
  tax_transform("identity", rank = "Genus") %>% 
  dist_calc("bray") %>% # bray curtis
  ord_calc() # automagically picks PCoA
```

Start interactive Shiny app. Note that the gif animation shown is from an outdated version of microViz. More recent versions of `ord_explore` allow editing the ordination shown, and generating `ord_plot` code.

``` {.r}
ord_explore(data = ord1, auto_caption = NA)
```

![Note: GIF is sped up x2: redrawing plots is not instantaneous, but pretty quick unless your dataset has many 1000s of samples. Date of creation: 29/03/2021](images/20210329_ord_explore_dev_x2.gif)

# Session info

```{r session info}
devtools::session_info()
```
---
title: "Videos"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

<div style="padding:56.25% 0 0 0;position:relative;"><iframe src="https://player.vimeo.com/video/559627248?badge=0&amp;autopause=0&amp;player_id=0&amp;app_id=58479" frameborder="0" allow="autoplay; fullscreen; picture-in-picture" allowfullscreen style="position:absolute;top:0;left:0;width:100%;height:100%;" title="microViz 10-minute intro/demo"></iframe></div><script src="https://player.vimeo.com/api/player.js"></script>

This demo shows how to use microViz to get started with easy exploration of your data with a Shiny app for creating interactive and reproducible ordination plots.

This demo was made for the Microbiome Data Congress 2021 and BioSB congress 2021.
---
title: "Ordination plots"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Ordination plots}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Intro

This article will show you how to create and customise ordination plots, like [PCA](#pca) and [RDA](#rda), with microViz.

For an even quicker start with ordinating your phyloseq data, check out the `ord_explore` [Shiny app](https://david-barnett.github.io/microViz/articles/web-only/ordination-interactive.html), which allows you to create ordination plots with a point and click interface, and generates `ord_plot` code for you to copy. Return to this article to understand more about creating and customising your ordination plotting script.

```{r setup}
library(phyloseq)
library(ggplot2)
library(microViz)
knitr::opts_chunk$set(fig.width = 7, fig.height = 6)
```

We will use example data from stool samples from an inflammatory bowel disease (IBD) study, borrowed from the great `corncob` package. See the article about [working with phyloseq objects](https://david-barnett.github.io/microViz/articles/web-only/phyloseq.html) if you want to get started with your own data, or just to learn more about manipulating phyloseq objects with microViz.

```{r get data}
ibd <- corncob::ibd_phylo
ibd
```

First we fix any uninformative tax_table entries and check the phyloseq object (see the article on [fixing your tax_table](https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html) for more info).

```{r fix it}
ibd <- tax_fix(ibd) # try tax_fix_interactive if you have problems with your own data
ibd <- phyloseq_validate(ibd, remove_undetected = TRUE)
```

## Motivation

Ordination plots are a great way to see any clustering or other patterns of microbiota (dis)similarity in (many) samples. Ordinations like PCA or PCoA show the largest patterns of variation in your data, and constrained ordination techniques like RDA or CCA can show you microbial variation that could be explained by other variables in your sample_data (but interpret constrained ordinations with care, and ideally test for the statistical significance of any hypothesised associations using a method like PERMANOVA, see `dist_permanova()`).

Check out the GUide to STatistical Analysis in Microbial Ecology (GUSTA ME)[website](https://sites.google.com/site/mb3gustame/) for a gentle theoretical introduction to [PCA](https://sites.google.com/site/mb3gustame/indirect-gradient-analysis/pca), [PCoA](https://sites.google.com/site/mb3gustame/dissimilarity-based-methods/principal-coordinates-analysis), [RDA](https://sites.google.com/site/mb3gustame/constrained-analyses/rda), [CCA](https://sites.google.com/site/mb3gustame/constrained-analyses/cca) and more.

Ordination plots can also be paired with barplots for greater insight into microbial compositions, e.g. see `ord_plot_iris()` and the `ord_explore()` interactive [Shiny app](https://david-barnett.github.io/microViz/articles/web-only/ordination-interactive.html).

## Prepare your microbes

When creating an ordination plot, you first need to prepare the microbiota variables.

-   Decide at which taxonomic rank to aggregate your data, e.g. "Genus"

-   Consider transforming the microbial counts, e.g. using the "clr" (centred log ratio) transformation, which is often recommended for compositional data ([like sequencing data](https://doi.org/10.3389/fmicb.2017.02224 "Gloor 2017"))

```{r clr}
ibd %>% 
  tax_transform(trans = "clr", rank = "Genus")
```

-   Some methods, such as [PCoA](#pcoa), require a matrix of pairwise distances between samples, which you can easily calculate with `dist_calc()`. Normally you should NOT transform your data when using a distance-based method, but it is useful to record an "identity" transformation anyway, to make it clear you have not transformed your data.

```{r identity}
ibd %>% 
  tax_transform(trans = "identity", rank = "Genus") %>% 
  dist_calc("bray") # bray curtis distance
```

-   Some dissimilarity measures, such as unweighted UniFrac, do not consider the abundance of each taxon when calculating dissimilarity, and so may be (overly) sensitive to differences in rare/low-abundance taxa. So you might want to filter out very rare taxa, with `tax_filter()` before using `dist_calc(ps, dist = "unifrac")`. Distances that are (implicitly) abundance weighted, including Generalised UniFrac, Bray-Curtis and Aitchison distance, should be less sensitive to rare taxa / filtering threshold choices.

### ps_extra

Notice that the objects created above are of class "ps_extra". This is just a simple S3 list object that holds your phyloseq object and additional stuff created from this phyloseq object, such as a distance matrix, as well as info on any transformation and aggregation applied to your taxa. microViz uses this to automatically create plot captions, to help you and your collaborators remember how you made each plot! You can access the phyloseq object, distance matrix and other parts of a ps_extra object with `ps_get()`, `dist_get()`, and friends.

## PCA {#pca}

Principle ***Components*** Analysis is an unconstrained method that does not use a distance matrix. PCA directly uses the (transformed) microbial variables, so you do not need `dist_calc()`. `ord_calc` performs the ordination (adding it to the ps_extra object) and `ord_plot()` creates the ggplot2 scatterplot (which you can customise like other ggplot objects).

Each point is a sample, and samples that appear closer together are typically more similar to each other than samples which are further apart. So by colouring the points by IBD status you can see that the microbiota from people with IBD is often, but not always, highly distinct from people without IBD.

```{r}
ibd %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2")
```

One benefit of not using a distance matrix, is that you can plot taxa "loadings" onto your PCA axes, using the plot_taxa argument. microViz plots all of the taxa loading vectors in light grey, and you choose how many of the vectors to label, starting with the longest arrows (alternatively you can name particular taxa to label).

The relative length of each loading vector indicates its contribution to each PCA axis shown, and allows you to roughly estimate which samples will contain more of that taxon e.g. samples on the left of the plot below, will typically contain more *Escherichia*/*Shigella* than samples on the right, and this taxon contributes heavily to the PC1 axis.

```{r}
ibd %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "ibd", shape = "DiseaseState", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "Dark2")
```

microViz also allows you directly visualize the sample compositions on a circular barplot or "iris plot" (named because it looks kinda like an eyeball) alongside the PCA plot. The samples on the iris plot are automatically arranged by their rotational position around the center/origin of the PCA plot.

```{r, fig.height=8, fig.width=7}
ibd %>% 
  tax_transform("clr", rank = "Genus") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot_iris(tax_level = "Genus", ord_plot = "above", anno_colour = "ibd")
```

Here we created the ordination plot as a quick accompaniment to the circular barchart, but it is more flexible to create and customise the ordination plot and iris plot separately, and then pair them afterwards with patchwork. See the `ord_plot_iris` [docs](https://david-barnett.github.io/microViz/reference/ord_plot_iris.html) for examples.

## PCoA {#pcoa}

Principle ***Co-ordinates*** Analysis is also an unconstrained method, but it does require a distance matrix. In an ecological context, a distance (or more generally a "dissimilarity") measure indicates how different a pair of (microbial) ecosystems are. This can be calculated in many ways.

### Aitchison distance

The [Euclidean distance](https://en.wikipedia.org/wiki/Euclidean_distance) is similar to the distance we humans are familiar with in the physical world. The results of a PCA is practically equivalent to a PCoA with Euclidean distances. The Aitchison distance is a dissimilarity measure calculated as the Euclidean distance between observations (samples) after performing a centred log ratio ("clr") transformation. That is why the Aitchison distance PCoA, below, looks the same as the PCA we made earlier. However, we cannot use plot_taxa, as the taxa loadings are only available for PCA (and related methods like RDA).

```{r}
ibd %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2")
```

Note that PCoA is also known as MDS, for (metric) Multi-Dimensional Scaling, hence the axes names.

### Ecological dissimilarities

Over the years, ecologists have invented numerous ways of quantifying dissimilarity between pairs of ecosystems. One ubiquitous example is the [Bray-Curtis](https://en.wikipedia.org/wiki/Bray–Curtis_dissimilarity) dissimilarity measure, shown below.

```{r}
ibd %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2")
```

Beyond Bray-Curtis, microViz `dist_calc()` can also help you calculate all the other ecological distances listed in `phyloseq::distanceMethodList` such the Jensen-Shannon Divergence, `"jsd"`, or Jaccard dissimilarity `"jaccard"`. Beware that if you want a binary dissimilarity measure from `vegan::vegdist()` (i.e. only using presence/absence info, and noting all the caveats about sensitivity to low abundance taxa) you will need to pass `binary = TRUE`, as below.

```{r}
ibd %>% 
  tax_transform("identity", rank = "Genus") %>% 
  dist_calc(dist = "jaccard", binary = TRUE) %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2")
```

### UniFrac distances

If you have a phylogenetic tree available, and [attached to your phyloseq object](https://david-barnett.github.io/microViz/articles/web-only/phyloseq.html#getting-your-data-into-phyloseq). You can calculate dissimilarities from the [UniFrac family](https://en.wikipedia.org/wiki/UniFrac) of methods, which take into account the phylogenetic relatedness of the taxa / sequences in your samples when calculating dissimilarity. Un-weighted UniFrac, `dist_calc(dist = "unifrac")`, does not consider the relative abundance of taxa, only their presence (detection) or absence, which can make it (overly) sensitive to rare taxa, sequencing artefacts, and abundance filtering choices. Conversely, weighted UniFrac, `"wunifrac"`, does put (perhaps too much) more importance on highly abundant taxa, when determining dissimilarities. The Generalised UniFrac, `"gunifrac"`, finds a balance between these two extremes, and by adjusting the `gunifrac_alpha` argument of `dist_calc()`, you can tune this balance to your liking (although the 0.5 default should be fine!).

Below is a Generalised UniFrac example using a different, and tiny, example dataset from the phyloseq package that has a phylogenetic tree.

You should **not** aggregate taxa before using a phylogenetic distance measure, but you can and probably should register the choice not to transform or aggregate, as below.

```{r}
data("esophagus", package = "phyloseq")
esophagus %>%
  phyloseq_validate(verbose = FALSE) %>% 
  tax_transform("identity", rank = "unique") %>% 
  dist_calc("gunifrac", gunifrac_alpha = 0.5)

```

## Further dimensions

You can show other dimensions / axes of an ordination than just the first two, by setting the axes argument. You can judge from the variation explained by each successive axis (on a scree plot) whether this is worthwhile information to show, e.g. in the example below, it could be interesting to also show the 3rd axis, but not any others.

```{r, fig.width=5, fig.height=3}
ibd %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_get() %>% 
  phyloseq::plot_scree() + theme(axis.text.x = element_text(size = 6))
```

Let us view the 1st and 3rd axes.

```{r}
ibd %>% 
  tax_transform("identity", rank = "Genus") %>% # don't transform!
  dist_calc("bray") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(axes = c(1, 3), color = "ibd", shape = "DiseaseState", size = 2) +
  scale_colour_brewer(palette = "Dark2") 
```

## RDA {#rda}

Redundancy analysis is a constrained ordination method. It displays the microbial variation that can also be explained by selected constraint variables.

*Behind the scenes, a linear regression model is created for each microbial abundance variable (using the constraints as the explanatory variables) and a PCA is performed using the fitted values of the microbial abundances.*

Starting from the same phyloseq object `ibd` the code below first creates a couple of binary (0/1) numeric variables to use a constraint variables. This is easy enough with ps_mutate.

Then we aggregate and transform our taxa, and like PCA we skip the dist_calc step.

```{r, fig.height=7, fig.width=7}
ibd %>%
  ps_mutate(
    IBD = as.numeric(ibd == "ibd"),
    Female = as.numeric(gender == "female"),
    Abx. = as.numeric(abx == "abx")
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("IBD", "Female", "Abx."),
    # method = "RDA", # Note: you can specify RDA explicitly, and it is good practice to do so, but microViz can guess automatically that you want an RDA here (helpful if you don't remember the name?)
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "DiseaseState", size = 2, alpha = 0.5, shape = "active",
    plot_taxa = 1:8
  )
```

### Customising your ordination plot

The plot above looks okay by default, but it is fairly easy to tweak ord_plot further to get the style just how you want it. The code below has comments to explain which part makes which changes to the plot.

```{r, fig.width=7, fig.height=8}
# first we make a function that replaces any unwanted "_" in our taxa labels with spaces
library(stringr) 
renamer <- function(x) str_replace(x, pattern = "_", replacement = " ")

ibd %>%
  ps_mutate(
    IBD = as.numeric(ibd == "ibd"),
    Female = as.numeric(gender == "female"),
    Abx. = as.numeric(abx == "abx")
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("IBD", "Female", "Abx."),
    method = "RDA",
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "DiseaseState", size = 2, alpha = 0.5, shape = "active",
    auto_caption = NA, # remove the helpful automatic caption
    plot_taxa = 1:8, taxon_renamer = renamer, # renamer is the function we made earlier
    tax_vec_length = 5, # this value is a scalar multiplier for the biplot score vectors
    tax_lab_length = 6, # this multiplier moves the labels, independently of the arrowheads
    tax_lab_style = tax_lab_style(size = 1.8, alpha = 0.5), # create a list of options to tweak the taxa labels' default style
    constraint_vec_length = 3, # this adjusts the length of the constraint arrows, and the labels track these lengths by default
    constraint_vec_style = vec_constraint(size = 1.5, alpha = 0.5), # this styles the constraint arrows
    constraint_lab_style = constraint_lab_style(size = 3) # this styles the constraint labels
  ) +
  # the functions below are from ggplot2:
  # You can pick a different colour scale, such as a color_brewer palette
  scale_colour_brewer(palette = "Set1") +
  # You can set any scale's values manually, such as the shapes used
  scale_shape_manual(values = c(
    active = "circle", mild = "circle cross",
    inactive = "circle open", control = "square open"
  )) +
  # this is how you add a title and subtitle
  ggtitle(
    label = "[Insert your exciting interpretations here?]",
    subtitle = "RDA with clr-transformed genera: constraints in red, taxa in black"
  ) +
  # and this is how you make your own caption
  labs(caption = "91 samples, 178 genera. Type 2 scaling.") +
  # this is one way to set the aspect ratio of the plot
  coord_fixed(ratio = 1, clip = "off")
```

#### Custom labels

`tax_lab_style()` is a helper function that gives you some options for 
customising the look of the taxa loading labels, including, in this example, 
using rotated, bold and italic text for the taxa names.

`constraint_lab_style()` is a similar helper function for customising the 
constraint labels. When rotating labels (not text) the ggtext package must be
installed.

```{r, fig.width=7, fig.height=6}
ibd %>%
  ps_mutate(
    IBD = as.numeric(ibd == "ibd"),
    Female = as.numeric(gender == "female"),
    Abx. = as.numeric(abx == "abx")
  ) %>%
  tax_transform("clr", rank = "Genus") %>%
  ord_calc(
    constraints = c("IBD", "Female", "Abx."),
    method = "RDA",
    scale_cc = FALSE # doesn't make a difference
  ) %>%
  ord_plot(
    colour = "DiseaseState", size = 2, alpha = 0.5, shape = "active",
    auto_caption = NA, 
    plot_taxa = 1:8, taxon_renamer = renamer,
    # with rotated labels, it is nicer to keep lab_length closer to vec_length
    tax_vec_length = 4.5, tax_lab_length = 4.6, 
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, fontface = "bold.italic"
    ), 
    constraint_vec_style = vec_constraint(size = 1.5, alpha = 0.5), 
    constraint_vec_length = 3, constraint_lab_length = 3.3,
    constraint_lab_style = constraint_lab_style(
      alpha = 0.8, size = 3, max_angle = 90, perpendicular = TRUE
    ) 
  ) +
  # SETTING A FIXED RATIO IDENTICAL TO THE aspect_ratio ARGUMENT IN
  # tax_lab_style() IS NECESSARY FOR THE ANGLES OF TEXT TO ALIGN WITH ARROWS!
  coord_fixed(ratio = 1, clip = "off", xlim = c(-6, 6)) +
  # The scales below are set the same as in the previous customisation:
  scale_colour_brewer(palette = "Set1") +
  scale_shape_manual(values = c(
    active = "circle", mild = "circle cross",
    inactive = "circle open", control = "square open"
  )) 
```


### Partial ordinations

Tutorial coming soon, for now see `ord_plot()` for examples.

## Session info

```{r session}
devtools::session_info()
```
---
title: "Fixing your tax_table"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Quick start / TLDR

Use `tax_fix()` on your phyloseq data with default arguments to repair most `tax_table` problems (missing or uninformative values). If you still encounter errors using e.g. `tax_agg`, try using the Shiny app `tax_fix_interactive()` to help you generate `tax_fix` code that will fix your particular `tax_table` problems.

------------------------------------------------------------------------

```{r setup}
library(phyloseq)
suppressPackageStartupMessages(library(microViz))
```

## Intro

This article will explain some of the common problems that can occur in your phyloseq `tax_table`, and that might cause problems for e.g. `tax_agg`. You can fix these problems with the help of `tax_fix` and `tax_fix_interactive`.

## Fixing problems

Let's look at some example data from the corncob package:

```{r example data}
pseq <- corncob::ibd_phylo
pseq
```

The Species rank appears to be blank for many entries. This is a problem you may well encounter in your data: unique sequences or OTUs often cannot be annotated at lower taxonomic ranks.

```{r example tt look}
tax_table(pseq)[40:54, 4:7] # highest 3 ranks not shown, to save space
```

If we would try to aggregate at Genus or Family rank level, we discover that blank values at these ranks prevent taxonomic aggregation. This is because, for example, it looks like OTU.43 and OTU.54 share the same (empty) Genus name, "", despite being different at a higher rank, Family.

```{r tax_agg fail}
# tax_agg(pseq, rank = "Family") # this fails, and sends (helpful) messages about taxa problems
```

So we should run `tax_fix` first, which will fix **most** problems with default settings, allowing taxa to be aggregated successfully (at any rank). If you still have errors when using tax_agg after tax_fix, carefully read the error and accompanying messages. Often you can copy suggested tax_fix code from the tax_agg error. You should generally also have a look around your tax_table for other uninformative values, using [tax_fix_interactive](https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html#interactive-solutions).

```{r tax_agg success}
pseq %>% tax_fix() %>% tax_agg(rank = "Family")
```

### What does tax_fix do?

`tax_fix` searches all the ranks of the phyloseq object `tax_table` for:

-   short values, like "g\_\_", "", " ", etc. (any with fewer characters than min_length)

-   common, longer but uninformative values like "unknown" (see full list at `?tax_fix`)

-   NAs

`tax_fix` replaces these values with the next higher taxonomic rank, e.g. an "unknown" Family within the Order Clostridiales will be renamed "Clostridiales Order", as seen below.

```{r tax_fix look}
pseq %>% tax_fix(min_length = 4) %>% 
  tax_agg("Family") %>% 
  ps_get() %>% # gets the phyloseq from the ps_extra list
  tax_table() %>% 
  .[1:8, 3:5] # removes the first 2 ranks and shows only first 8 rows for nice printing
```

### Interactive solutions

You can use `tax_fix_interactive()` to explore your data's `tax_table` visually, and interactively find and fix problematic entries. You can then copy your automagically personalised `tax_fix` code from `tax_fix_interactive`'s output, to paste into your script. Below is a screen capture video of `tax_fix` in action, using some other artificially mangled example data (see details at `?tax_fix_interactive()`).

``` {.r .R}
tax_fix_interactive(example_data)
```

![](images/20210411_tax_fix_interactive_x2.gif)

### Other possible problems

1.  **Completely unclassified taxa** (aka taxa where all values in their `tax_table` row are either too short or listed in unknowns argument) will be replaced at all ranks with their unique row name by default (or alternatively with a generic name of "unclassified [highest rank]", which is useful if you want to aggregate all the unclassified sequences together with `tax_agg()`)
2.  **Unclassified taxa that also have short / unknown row names**, e.g. the unclassified taxon called "-1" in the example "enterotype" dataset from phyloseq. If something like this happens in your data, rename the taxa manually, (e.g. `taxa_names(enterotype)[1] <- "unclassified taxon"` or give them all completely different names with `tax_name()`.
3.  T**axa with the same `tax_table` entry repeated across multiple ranks**: This is a problem for functions like `taxatree_plots()`, which need distinct entries at each rank to build the tree structure for plotting. This might happen after you `tax_fix` data with problem 1. of this list, or in data from e.g. microarray methods like HITchip. The solution is to use `tax_prepend_ranks()` (after `tax_fix`) to add the first character of the rank to all tax_table entries (you will also need set the `tax_fix` argument suffix_rank = "current").
4.  **Informative but duplicated `tax_table` entries:** e.g. you don't want to delete/replace a genus name completely, but it is shared by two families and thus blocking `tax_agg`. The solution is to rename (one of) these values manually to make them distinct. `tax_table(yourPhyloseq)["targetTaxonName", "targetRank"] <- "newBetterGenusName"`
5.  **Really long taxa_names():** e.g. you have DNA sequences as names. See `tax_name()` for an easy way to rename all your taxa.

## Abundance filtering as a solution

Sequences that are unclassified at fairly high ranks e.g. Class are often very low abundance (or possibly represent sequencing errors/chimeras), if you are using data from an environment that is typically well represented in reference databases. So if you are struggling with what to do with unclassified taxa, consider if you can just remove them first using `tax_filter()` (perhaps using fairly relaxed filtering criteria like min_prevalence of 2 samples, or min_total_abundance of 1000 reads, and keeping the tax_level argument as NA, so that no aggregation is attempted before filtering).

## Alternatives

`microbiome::aggregate_taxa()` also solves some `tax_table` problems, e.g. where multiple distinct genera converge again to the same species name like "" or "s\_\_", it will make unique taxa names by pasting together **all** of the rank names. However this can produce some very long names, which need to be manually shortened before use in plots. Plus, it doesn't replace names like "s\_\_" if they only occur once. Moreover, when creating ordination plots with microViz, only `tax_agg()` will record the aggregation level for provenance tracking and automated plot captioning.

## Session info

```{r session info}
devtools::session_info()
```
---
title: "Visualising compositions"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Visualising compositions}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
editor_options: 
  markdown: 
    wrap: 72
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r, message=FALSE}
options(width = 100)
library(microViz)
library(phyloseq)
library(ggplot2)
library(patchwork) # for arranging groups of plots
knitr::opts_chunk$set(
  fig.height = 6,
  fig.width = 9
)
```

```{r example data}
# get example phyloseq data from corncob package and tidy up
pseq <- corncob::ibd_phylo %>% 
  tax_filter(min_prevalence = 2) %>% 
  tax_fix() %>% 
  phyloseq_validate()
```

The `comp_barplot` function allows you to visualise the taxonomic
compositions of your microbiome samples in a flexible, scalable,
group-able, and visually appealing way.

## Quick example barplot

Visualise the top Genera across all the female samples from this
inflammatory bowel disease study dataset. The order of the samples is
automatically set by their "bray"-curtis dissimilarity.

By default, the top 8 taxa are shown. These taxa are chosen by their
total count abundance across all plotted samples.

```{r simple}
pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(tax_level = "Genus") + 
  coord_flip() # horizontal bars are often more readable
```

## Customising this barplot

The output of comp_barplot can be customised in several ways. See the
comment alongside each argument for an explanation of its effect.

```{r}
pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(
    tax_level = "Genus",
    label = "DiseaseState", # name an alternative variable to label axes
    n_taxa = 15, # give more taxa unique colours
    merge_other = FALSE, # split the "other" category to display alpha diversity
    bar_width = 0.7, # reduce the bar width to 70% of one row
    bar_outline_colour = "grey5" # is the default (use NA to remove outlines)
  ) + 
  coord_flip()
```

**Other notes:**

-   Dissimilarity is calculated using only the visibly distinct taxa, to
    optimise sorting for visual similarity. You can change this by
    setting `order_with_all_taxa = TRUE`, to always use all taxa for
    similarity sorting.
-   The colour palette is important, to allow (adjacent) taxa to be
    distinguished. The palette microViz uses is generated by the
    `distinct_palette` function, which starts with the Paired and Dark2
    palettes from ColorBrewer and continues with further distinct
    colours generated at <http://medialab.github.io/iwanthue/> (all
    colors, soft k-means).

## Averages, faceting or grouping?

### Averaging compositions

Sometimes, to compare microbial compositions across groups, average
compositions are presented. However that "group-averaging" approach
hides a lot of within-group variation, as well as any imbalance in group
sizes.

```{r}
pseq %>%
  ps_select(age, DiseaseState) %>% # avoids lots of phyloseq::merge_samples warnings
  ps_filter(DiseaseState != "IBDundef") %>% 
  phyloseq::merge_samples(group = "DiseaseState") %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 12,
    bar_width = 0.8
  ) +
  coord_flip() + labs(x = NULL, y = NULL)
```

### Faceting

Faceting is where you show each group on a small subplot.

In the plot below can you see that at minority of UC samples have a high
abundance of Escherichia/Shigella or Streptococcus. The merged bars
above might have misled you into thinking all UC samples had somewhat
increased abundances of these taxa.

```{r}
pseq %>%
  ps_filter(DiseaseState != "IBDundef") %>% # only one sample in this group
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    DiseaseState = factor(
      DiseaseState, levels = c("UC", "nonIBD", "CD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    bar_outline_colour = NA, facet_by = "DiseaseState"
  ) +
  coord_flip() 
```

Instead of using the `facet_by` argument in `comp_barplot` you can have
more control over faceting by doing it yourself afterwards. You can use
facet_grid to create row facets.

```{r, fig.height=8, fig.width=6}
pseq %>%
  ps_filter(DiseaseState != "IBDundef") %>% # only one sample in this group
  # convert DiseaseState into ordered factor to control order of facets
  ps_mutate(
    DiseaseState = factor(
      DiseaseState, levels = c("UC", "CD", "nonIBD"), ordered = TRUE)
  ) %>% 
  comp_barplot(
    tax_level = "Genus", n_taxa = 15,
    sample_order = "bray", bar_outline_colour = NA,
  ) +
  facet_grid(
    rows = vars(DiseaseState), 
    scales = "free", space = "free" # these options are critically important!
  ) +
  coord_flip() + 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
```

### Grouping

For even greater control than faceting, `comp_barplot` allows you to
generate separate ggplot objects for each group, whilst maintaining the
same taxa colour scheme.

You can assemble these plots into one figure with, for example, the
`patchwork` package, or keep them separate.

*Note that the ordering of the samples may differ between facet and
group_by approaches. In the group_by method, the ordering of the samples
by similarity is done separately for each group, whereas in the facet_by
method, similarity-based ordering is done with all samples and then the
samples are separated by facet afterwards.*

```{r}
plot_list <- pseq %>% 
  ps_filter(DiseaseState != "IBDundef") %>% 
  comp_barplot(
  n_taxa = 15, tax_level = "Genus", group_by = "DiseaseState"
)

# Plot them side by side with the patchwork package.
patch <- patchwork::wrap_plots(plot_list, nrow = 1, guides = "collect")
patch & coord_flip() # make all plots horizontal (note: use & instead of +)
```

Notice how you can theme all plots with the `&` operator.

See <https://patchwork.data-imaginist.com/index.html> for more examples
of arranging multiple plots.

```{r}
patch &
  coord_flip() & labs(x = NULL, y = NULL) &
  theme(
    axis.text.y = element_text(size = 5),
    legend.text = element_text(size = 8)
  ) &
  plot_annotation(
    title = "Microbial composition across disease groups",
    caption = "Caption: patchwork is a great package!",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )
```

## Sorting the barplot

### Sorting by similarity

Sorting the samples on compositional barplots by similarity can make
patterns in the data much easier to see. Check out this unsorted version
of the first barplot in this article.

```{r}
pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(tax_level = "Genus", sample_order = "default") + 
  coord_flip() +
  ggtitle("Unsorted barcharts are hard to read!")
```

You can play with the dissimilarity measure (set in `sample_order`
argument) and `seriate_method` if you like, but the defaults (Bray
Curtis and OLO Ward) seem to work pretty well most of the time.

When sorting samples by similarity, the default is to treat the "other"
taxa as one group, i.e. when `merge_other = TRUE` and
`order_with_all_taxa = FALSE`.

If you set `order_with_all_taxa = TRUE`, samples are sorted BEFORE
merging taxa. The resulting sample order is then the same as when
`merge_other = FALSE`.

```{r}
pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(tax_level = "Genus") + 
  coord_flip() +
  ggtitle("Samples sorted AFTER merging 'other' taxa")

pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(tax_level = "Genus", order_with_all_taxa = TRUE) + 
  coord_flip() +
  ggtitle("Samples sorted BEFORE merging 'other' taxa")


pseq %>%
  ps_filter(gender == "female") %>%
  comp_barplot(tax_level = "Genus", merge_other = FALSE) + 
  coord_flip() +
  ggtitle("'other' taxa not merged")

```

### Sort by ordination

Coming soon.

### Sort by 1 taxon

To study the distribution of a single taxonomic group across your
samples, you can use `ps_arrange` (with the `.target` argument set to
"otu_table") and the 'default' sample_order setting in `comp_barplot`.

```{r}
pseq %>% 
   tax_agg("Phylum") %>% 
   tax_transform("compositional") %>% 
   ps_arrange(desc(Firmicutes), .target = "otu_table") %>% 
   comp_barplot(tax_level = "Phylum", sample_order = "default") +
  coord_flip()
```

### Sorting by metadata

Sorting across timepoint groups in another dataset, dietswap.

```{r diet only}
data("dietswap", package = "microbiome")
ps <- dietswap %>% ps_filter(group == "DI")

# subset to participants/"subjects" with samples for both timepoints
timepoint_groups <- sample_data(ps) %>%
  data.frame() %>%
  split.data.frame(f = .$timepoint.within.group)
have_both_timepoints <- intersect(timepoint_groups$`1`$subject, timepoint_groups$`2`$subject)
ps <- ps %>% subset_samples(subject %in% have_both_timepoints)

# define grouping variable for plotting
ps <- ps %>% ps_mutate(
  plot_groups = interaction(nationality, timepoint.within.group)
)
```

Grouped by both timepoint and another grouping factor, nationality in
this example.

```{r, fig.height=7.5, fig.width=10}
times_list <- ps %>%
  ps_arrange(timepoint.within.group, nationality, desc(subject)) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 10, sample_order = "default",
    group_by = "plot_groups", bar_width = 0.7, label = "subject"
  )

times <- wrap_plots(
  times_list[c("AAM.1", "AAM.2", "AFR.1", "AFR.2")],
  byrow = TRUE, ncol = 2, guides = "collect", heights = c(5, 4)
) &
  coord_flip() &
  theme(
    text = element_text(size = 10),
    axis.title = element_blank()
  )

times
```

Same grouping, now showing diversity of taxa within other, with
`merge_other = FALSE`

```{r, fig.height=7.5, fig.width=10}
times_list <- ps %>%
  ps_arrange(timepoint.within.group, nationality, desc(subject)) %>%
  comp_barplot(
    tax_level = "Genus", n_taxa = 10, sample_order = "default",
    merge_other = FALSE, bar_outline_colour = "grey10",
    group_by = "plot_groups", bar_width = 0.7, label = "subject"
  )

times <- wrap_plots(
  times_list[c("AAM.1", "AAM.2", "AFR.1", "AFR.2")],
  byrow = TRUE, ncol = 2, guides = "collect", heights = c(5, 4)
) &
  coord_flip() &
  theme(
    text = element_text(size = 10),
    axis.title = element_blank()
  )

times
```

```{r}
devtools::session_info()
```
---
title: "Working with phyloseq objects"
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Intro

This tutorial will show you how microViz makes working with phyloseq objects easier.

```{r setup}
library(dplyr)
library(phyloseq)
library(microViz)
```

## Getting your data into phyloseq

phyloseq objects are probably the most commonly used data format for working with microbiome data in R.

The creator of phyloseq, [Paul J. McMurdie](https://github.com/joey711), explains the structure of phyloseq objects and how to construct them on the [phyloseq website](https://joey711.github.io/phyloseq/import-data.html#phyloseq-ize_data_already_in_r).

-   [biom format](http://biom-format.org/) files can be imported to phyloseq with the [`import_biom`](https://joey711.github.io/phyloseq/import-data.html#the_import_family_of_functions) function. Such biom files are generated (or can be) from many processing tools including QIIME [1](http://qiime.org/tutorials/index.html) / [2](https://docs.qiime2.org/2021.8/tutorials/exporting/), [MetaPhlAn](http://segatalab.cibio.unitn.it/tools/metaphlan/), and [NG-Tax](https://doi.org/10.3389/fgene.2019.01366).

-   [mothur](https://mothur.org/) output is also directly supported via phyloseq's [`import_mothur`](https://joey711.github.io/phyloseq/import-data.html#_import_mothur) function.

-   [DADA2](https://benjjneb.github.io/dada2/) output can converted to phyloseq according to these [DADA2 handoff instructions](https://benjjneb.github.io/dada2/tutorial.html#bonus-handoff-to-phyloseq)

Don't worry too much about getting all of your sample metadata into your biom file or phyloseq object at the start, as `ps_join()` makes it easy to add sample data later.

```{r example data}
# This tutorial will just use some example data 
# It is already available as a phyloseq object from the corncob package
example_ps <- corncob::ibd_phylo
example_ps
```

### Validating your phyloseq

phyloseq checks that your sample and taxa names are consistent across the different slots of the phyloseq object. microViz provides `phyloseq_validate()` to check for and fix other possible problems with your phyloseq that might cause problems in later analyses. It is recommended to run this at the start of your analyses, and fix any problems identified.

```{r validate ps}
example_ps <- phyloseq_validate(example_ps, remove_undetected = TRUE)
example_ps <- tax_fix(example_ps)

```

## microViz and phyloseq overview

Once you have a valid phyloseq object, microViz provides several helpful functions for manipulating that object. The names and syntax of some functions will be familiar to users of dplyr, and reading [the dplyr help pages](https://dplyr.tidyverse.org/) may be useful for getting the most out of these functions.

-   Add a data.frame of metadata to the sample_data slot with `ps_join()`

-   Compute new sample_data variables with `ps_mutate()`

-   Subset or reorder variables in sample_data with `ps_select()`

-   Subset samples based on sample_data variables with `ps_filter()`

-   Reorder samples (can be useful for examining sample data or plotting) with:

    -   `ps_reorder()` : manually set sample order

    -   `ps_arrange()` : order samples using sample_data variables

    -   `ps_seriate()` : order samples according to microbiome similarity

-   Remove duplicated/repeated samples with `ps_dedupe()`

-   Remove samples with missing values in sample_data with `ps_drop_incomplete()`

Check out the examples section on each function's [reference](https://david-barnett.github.io/microViz/reference/index.html#section-manipulating-sample-data) page for extensive usage examples.

## Example sample manipulation

Lets look at the sample data already in our example phyloseq.

```{r example data tbl}
# 91 samples, with 15 sample_data variables
example_ps 
# return sample data as a tibble for pretty printing
samdat_tbl(example_ps)
```

Maybe you want to only select participants who have IBD (not controls). You can do that by filtering samples based on the values of the sample data variable: ibd

```{r select ibd}
example_ps %>% ps_filter(ibd == "ibd") # it is essential to use `==` , not just one `=`
# notice that taxa that no longer appear in the remaining 67 samples have been removed!
```

More complicated filtering rules can be applied. Let's say you want female IBD patients with "mild" or "severe" activity, who are at least 13 years old.

```{r select ibd 2}
partial_ps <- example_ps %>% 
  ps_filter(
    gender == "female",
    activity %in% c("mild", "severe"),
    age >= 13
  )
partial_ps
```

Let's have a look at the sample data of these participants. We will also arrange the samples grouped by disease and in descending age order, and select only a few interesting variables to show.

```{r arrange and select}
partial_ps %>% 
  ps_arrange(DiseaseState, desc(age)) %>% 
  ps_select(DiseaseState, age, matches("activ"), abx) %>% # selection order is respected
  samdat_tbl() # this adds the .sample_name variable
```

You can also sort sample by microbiome similarity with `ps_seriate()`.

```{r}
partial_ps %>% 
  tax_agg("Genus") %>% 
  ps_seriate(dist = "bray", method = "OLO_ward") %>% # these are the defaults
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10)
# note that comp_barplot with sample_order = "bray" will run 
# this ps_seriate call internally, so you don't have to!

```

You can also arrange samples by abundance of one of more microbes using `ps_arrange()` with `.target = "otu_table"`. Arranging by taxon can only be done at the current taxonomic rank, so we will aggregate to Genus level first.

```{r}
# Arranging by decreasing Bacteroides abundance
partial_ps %>% 
  tax_agg("Genus") %>% 
  ps_arrange(desc(Bacteroides), .target = "otu_table") %>% 
  otu_get() %>% # get the otu table
  .[, 1:6] # show only a subset of the otu_table

# Plot samples' compositions in this order
partial_ps %>% 
  tax_agg("Genus") %>% 
  ps_arrange(desc(Bacteroides), .target = "otu_table") %>% 
  comp_barplot(tax_level = "Genus", sample_order = "default", n_taxa = 10)
# Notice this is sorted by bacteroides counts 
# (this doesn't quite match relative abundance % due to sequencing depth variation)
```

## Other notes:

### `ps_filter()` vs. `phyloseq::subset_samples()`

As well as filtering your samples, `ps_filter()` might also modify the otu_table and tax_table of the phyloseq object (unlike `phyloseq::subset_samples()`, which never does this).

Why does it do this?\
If you remove many samples from your dataset, often your phyloseq object will be left with taxa that never occur in any of the remaining samples (i.e. total counts of zero). `ps_filter()` removes those absent taxa by default.

If you don't want this, you can set the `.keep_all_taxa` argument to `TRUE` in `ps_filter`.

## Technical log

```{r session info}
devtools::session_info()
```
---
title: "Example analyses with atlas1006 data"
---

This document shows an example `microViz` analysis workflow using
example data `atlas1006` from the `microbiome` package.

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(width = 100)
```

## Load necessary packages

```{r load packages, message=FALSE}
library(dplyr)
library(ggplot2)
library(patchwork)
library(phyloseq)
library(microbiome)
library(microViz)
```

## Get example data

```{r data setup}
# get some example data
data("atlas1006", package = "microbiome")

# create a couple of numerical variables (often more useful than character strings)
ps <- atlas1006 %>%
  ps_mutate(
    weight = recode(
      .x = bmi_group, morbidobese = 5, severeobese = 4, 
      obese = 3, overweight = 2, lean = 1, underweight = 0
    ),
    lean = if_else(weight < 2, true = 1, false = 0, missing = 0),
    female = if_else(sex == "female", true = 1, false = 0),
    extract_P = if_else(DNA_extraction_method == "p", true = 1, false = 0)
  ) %>%
  # only look at the baseline time point if multiple samples available
  # and drop samples with no DNA extraction method info
  ps_filter(time == 0, !is.na(DNA_extraction_method)) %>%
  # remove the sample data variables about time
  ps_select(-time)

# add a couple of missing values to show how microViz handles missing data
sample_data(ps)$female[c(3, 4)] <- NA

# look at phyloseq object description
ps
```

```{r data prep}
# this example data has a slightly odd tax_table because it comes from HITChip data (instead of sequencing)
# the taxonomic assignment is done differently, so we need to ensure taxon names are not repeated across columns
# it can otherwise be used the same as sequencing data for this example
tax_table(ps) %>% head(3)
ps <- ps %>% tax_prepend_ranks()
tax_table(ps) %>% head(3)
# replace any uninformative tax_table values
ps <- ps %>% tax_fix()

# look at the effect of removing rare Genera, e.g. how many Genera are present in at least 5% of samples?
ps %>% tax_filter(min_prevalence = 5 / 100, tax_level = "Genus")
# we will use this udring other analyses, but not overwrite the original phyloseq object as the
# unfiltered set of taxa would be required if we were performing e.g. alpha diversity analyses
```

## Exploring your data

### Ordination plots

Ordination methods can help you to visualize variation in overall
microbial ecosystem composition, and look at whether it might differ
markedly between groups, e.g. weight.

Here is one option to try first:

1.  Filter out rare taxa (e.g. remove Genera not present in at least 5%
    of samples) - use `tax_filter()`
2.  Aggregate the taxa into bacterial families (for example) - use
    `tax_agg()`
3.  Transform the microbial data with the centered-log-ratio
    transformation - use `tax_transform()`
4.  Perform PCA with the CLR-transformed features (equivalent to
    Aitchison distance PCoA) - use `ord_calc()`
5.  Plot the first 2 axes of this PCA ordination, colouring samples by
    group and adding taxon loading arrows to visualise which taxa
    generally differ across your samples - use `ord_plot()`
6.  Customise the theme of the ggplot as you like and/or add features
    like ellipses or annotations

```{r calculate clr pca}
clr_pca <- ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>% # aggregate taxa at Genus level
  tax_transform("clr") %>% # centered log ratio transformation
  ord_calc(method = "PCA") # Note: PCA is RDA without constraints (& ord_calc uses an RDA method to perform PCA)
```

```{r print clr pca}
clr_pca
```

```{r}
clr_pca %>%
  ord_plot(colour = "weight", shape = "DNA_extraction_method", alpha = 0.7, size = 1.5) +
  scale_colour_viridis_c(option = "inferno", direction = -1, end = 0.8) +
  scale_shape_manual(
    values = c(o = "square open", r = "triangle open", p = "circle"),
    name = "DNA\nextraction\nmethod"
  )
```

### Taxonomic compositions?

Using a PCA ordination allows us reliably draw biplots, showing which
taxa are associated with this major variation in sample microbiota
composition as represented on axes PC1 and PC2.

```{r clr pca with tax loading}
pca <- clr_pca %>%
  ord_plot(
    colour = "weight", shape = "DNA_extraction_method", alpha = 0.7, size = 1,
    plot_taxa = 12:1, tax_vec_length = 0.3, 
    taxon_renamer = function(x) stringr::str_remove_all(x, "^G: | [ae]t rel."),
    center = TRUE,
    tax_lab_style = tax_lab_style(
      type = "label", max_angle = 90, fontface = "bold", 
      alpha = 0.8, size = 2
    )
  ) +
  scale_colour_viridis_c(option = "inferno", direction = -1, end = 0.8) +
  scale_shape_manual(
    values = c(o = "square open", r = "triangle open", p = "circle"),
    name = "DNA\nextraction\nmethod"
  ) +
  # essential for correct label angles
  coord_fixed(clip = "off", xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5))
```

```{r, fig.height=6, fig.width=8}
pca
```

What about the actual compositions of these 800 samples?

```{r, fig.width=8, dpi=120}
iris <- clr_pca %>%
  ord_plot_iris(
    tax_level = "Genus", n_taxa = 15,
    anno_binary = "extract_P",
    anno_binary_style = list(y = 1.1, size = 0.5, alpha = 0.3, shape = 1),
    taxon_renamer = function(x) stringr::str_remove_all(x, "^G: | [ae]t rel.")
  )

iris
```

You can arrange both plots together with patchwork package.

```{r, fig.height=10, fig.width=8, dpi=120}
pca / iris 
```

## Testing hypotheses

### PERMANOVA

What variables is the overall microbial community variation associated
with?

```{r permanova}
ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  dist_calc("aitchison") %>%
  dist_permanova(
    variables = c("DNA_extraction_method", "weight", "female"),
    n_perms = 99, # this is a low number of permutations for speed, you should set more e.g. 9999
    seed = 12345, complete_cases = TRUE, verbose = "max"
  )
```

### Visualising significant predictors?

(Partial) Constrained ordination can be useful to show microbial
variation explained by predictors significant in PERMANOVA analyses.

```{r constraining prep}
# constraints need to be on the same or similar scales for comparability
# so make binary variables and scale the weight variable
ps <- ps %>%
  ps_mutate(
    wt_scaled = c(scale(weight, center = TRUE, scale = TRUE)),
    P = if_else(DNA_extraction_method == "p", 1, 0),
    R = if_else(DNA_extraction_method == "r", 1, 0),
    O = if_else(DNA_extraction_method == "o", 1, 0)
  )
```

```{r constrained ord}
constr_ord <- ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  tax_transform("clr") %>%
  ord_calc(
    method = "RDA", constraints = c("female", "wt_scaled", "P", "R", "O")
  )
```

```{r print constrained ord}
constr_ord
```

```{r}
ord_plot(constr_ord, plot_taxa = 10:1, colour = "DNA_extraction_method", shape = 1) +
  scale_color_brewer(palette = "Set2", name = "DNA\nextraction\nmethod")
```

As DNA extraction method dominates the plot above, we could try
"partialing out" the variation explained by DNA extraction method, prior
to constraining on the other factors of interest.

```{r partial constrained ord}
ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  tax_agg("Genus") %>%
  tax_transform("clr") %>%
  ord_calc(
    method = "RDA", conditions = c("P", "R", "O"),
    constraints = c("female", "wt_scaled")
  ) %>%
  ord_plot(plot_taxa = 10:1, colour = "DNA_extraction_method", shape = 1) +
  scale_color_brewer(palette = "Set2", name = "DNA\nextraction\nmethod")
```

### Taxon models

What are the effects of these factors on individual taxa? Let's use beta
binomial regression models to find out. We will skip fitting models for
genera here, only for speed of generating this example. You should
include all ranks (which is the default) for the best use of
taxatree_plots.

See the taxon modelling article for a more comprehensive look at the 
taxatree_* family of functions.

```{r taxon modelling}
library(corncob)
tt_models <- ps %>%
  tax_filter(min_prevalence = 5 / 100, tax_level = "Genus") %>%
  taxatree_models(
    ranks = c("Phylum", "Family"), 
    variables = c("female", "wt_scaled", "P", "R"),
    type = "bbdml", verbose = "max"
  )
tt_stats <- taxatree_models2stats(tt_models, param = "mu")
```

Visualize the results compactly on a microViz taxonomic association
tree.

```{r, fig.height=6}
tt_stats %>% 
  taxatree_plotkey(
    node_size_range = c(2, 8), size_stat = list(mean = mean),
    rank == "Phylum", 
    taxon_renamer = function(x) stringr::str_remove_all(x, "^.: | [ae]t rel.")
  )
```

```{r taxatree plots, fig.height=7}
tt_stats %>% 
  taxatree_plots(
    node_size_range = c(1, 5), size_stat = list(mean = mean)
  ) %>%
    .[c("P", "R", "female", "wt_scaled")] %>% 
  wrap_plots(., ncol = 2, guides = "collect")
```

#### Example interpretation (illustrative only):

-   DNA extraction methods P and R are associated with significantly
    higher levels of Actinobacteria and lower Bacteroides, relative to
    the reference of extraction method O.
-   There are associations between weight and various taxa, but these
    are not as strong as the associations with extraction method.

## Disclaimer

This document is intended only to be an example of the kind of analyses
and visualization you can do with microViz. The analysis of the
atlas1006 data is not intended to be considered theoretically sound or
biologically interpretable.

## Session info

```{r}
devtools::session_info()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_fix_interactive.R
\name{tax_fix_interactive}
\alias{tax_fix_interactive}
\title{Shiny app to help you use tax_fix}
\usage{
tax_fix_interactive(data, app_options = list(launch.browser = TRUE))
}
\arguments{
\item{data}{a phyloseq object}

\item{app_options}{options list passed to shinyApp()}
}
\value{
nothing
}
\description{
Try this app if you get errors with \code{tax_fix()} that are tricky to work past,
or suggestions to use \code{tax_fix()} that you don't understand.

The app shows you the tax_table of your data (searchable) with unknown values highlighted.

It allows you to interactively modify minimum allowed length and to select
further values to be defined as unknown.

It will show you the correct \code{tax_fix()} code to copy paste into your script
to reproduce the interactive filtering.
}
\examples{
library(dplyr)
library(phyloseq)

# create some problem-filled example tax_table data
data(dietswap, package = "microbiome")
ps <- dietswap
# create unknowns to test filling
tt <- tax_table(ps)
ntax <- ntaxa(ps)
set.seed(123)
g <- sample(1:ntax, 30)
f <- sample(g, 10)
p <- sample(f, 3)
tt[g, 3] <- "g__"
tt[f, 2] <- "f__"
tt[p, 1] <- "p__"
tt[sample(1:ntax, 10), 3] <- "unknown"
# create a row with only NAs
tt[1, ] <- NA
tax_table(ps) <- tax_table(tt)

# function takes a phyloseq and shows code for how to fix the tax_table
# tax_fix_interactive(data = ps)
}
\seealso{
\code{\link{tax_fix}} for the non-interactive function to use in your scripts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-taxAnnotation.R
\name{anno_tax_density}
\alias{anno_tax_density}
\title{Helper to specify heatmap annotation for showing taxa abundance density plot}
\usage{
anno_tax_density(
  undetected = 0,
  only_detected = TRUE,
  trans = "log10p",
  zero_replace = 0,
  use_counts = TRUE,
  size = grid::unit(30, "mm"),
  type = c("lines", "violin", "heatmap"),
  xlim = NULL,
  heatmap_colors = c("white", "forestgreen"),
  joyplot_scale = 1.5,
  border = TRUE,
  gp = grid::gpar(fill = "lightgrey"),
  axis = TRUE,
  ...,
  data = NULL,
  taxa = NULL,
  which = NULL
)
}
\arguments{
\item{undetected}{the value above which taxa are classed as detected/present in a sample}

\item{only_detected}{only plot values for samples where the taxon abundance is > undetected}

\item{trans}{name of transformation suitable for tax_transform,
or a function calling tax_transform, and/or tax_scale,
(a function must take a phyloseq or ps_extra, and return one)}

\item{zero_replace}{zero_replace value for for tax_transform, ignored if trans is a function}

\item{use_counts}{try to retrieve counts from data object?}

\item{size}{width or height as a grid unit object}

\item{type}{Type of graphics to represent density distribution. "lines" for normal density plot; "violine" for violin plot and "heatmap" for heatmap visualization of density distribution.}

\item{xlim}{Range on x-axis.}

\item{heatmap_colors}{A vector of colors for interpolating density values.}

\item{joyplot_scale}{Relative height of density distribution. A value higher than 1 increases the height of the density distribution and the plot will represented as so-called "joyplot".}

\item{border}{Wether draw borders of the annotation region?}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_density]{ComplexHeatmap::anno_density}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{taxa}{OPTIONAL selection vector of taxa (names, numbers or logical),
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to taxAnnotation(),
which itself is used by cor_heatmap and comp_heatmap as tax_anno argument.
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- tax_filter(ibd_phylo, min_prevalence = 5)
psq <- tax_mutate(psq, Species = NULL)
psq <- tax_fix(psq)
psq <- tax_agg(psq, rank = "Family")
taxa <- tax_top(psq, n = 15, rank = "Family")
# makes a function that takes data, taxa and which (at minimum)
fun <- anno_tax_density()
# manually specify the density plot function by giving it data etc.
heatmapAnnoFunction <- fun(data = psq, which = "column", taxa = taxa)

# draw the density plot without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)

# let's change some style options and specify the data up front
grid.newpage()
pushViewport(vp)
draw(anno_tax_density(
  data = psq, taxa = taxa, which = "row",
  gp = grid::gpar(fill = "red"), border = FALSE
))

# heatmap type, with alternative transformation and axis_param
grid.newpage()
pushViewport(vp)
draw(anno_tax_density(
  data = psq, taxa = taxa, which = "row", type = "heatmap",
  trans = "log2", zero_replace = "halfmin", axis_param = list(labels_rot = 0)
))

grid.newpage()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-sampleAnnotation.R
\name{sampleAnnotation}
\alias{sampleAnnotation}
\title{Helper to specify a HeatmapAnnotation for variables in cor_heatmap}
\usage{
sampleAnnotation(
  ...,
  name,
  annotation_legend_param = list(),
  show_legend = TRUE,
  gp = grid::gpar(col = NA),
  border = FALSE,
  gap = grid::unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_label = NULL,
  annotation_name_gp = grid::gpar(),
  annotation_name_offset = NULL,
  annotation_name_rot = NULL,
  annotation_name_align = FALSE,
  annotation_name_side = "auto",
  .data = NULL,
  .samples = NULL,
  .side = NULL
)
}
\arguments{
\item{...}{Name-value pairs where the names correspond to annotation names and values
are the output of sample annotation functions
such as anno_sample(), or manually specified AnnotationFunction objects}

\item{name}{Name of the heatmap annotation, optional.}

\item{annotation_legend_param}{A list which contains parameters for annotation legends. See \code{\link[ComplexHeatmap]{color_mapping_legend,ColorMapping-method}} for all possible options.}

\item{show_legend}{Whether show annotation legends. The value can be one single value or a vector.}

\item{gp}{Graphic parameters for simple annotations (with \code{fill} parameter ignored).}

\item{border}{border of single annotations.}

\item{gap}{Gap between annotations. It can be a single value or a vector of \code{\link[grid]{unit}} objects.}

\item{show_annotation_name}{Whether show annotation names? For column annotation, annotation names are drawn either on the left or the right, and for row annotations, names are draw either on top or at the bottom. The value can be a vector.}

\item{annotation_label}{Labels for the annotations. By default it is the same as individual annotation names.}

\item{annotation_name_gp}{Graphic parameters for anntation names. Graphic paramters can be vectors.}

\item{annotation_name_offset}{Offset to the annotation names, a \code{\link[grid]{unit}} object. The value can be a vector.}

\item{annotation_name_rot}{Rotation of the annotation names. The value can be a vector.}

\item{annotation_name_align}{Whether to align the annotation names.}

\item{annotation_name_side}{Side of the annotation names.}

\item{.data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{.samples}{OPTIONAL selection vector of sample names,
only set this if providing .data argument to override default}

\item{.side}{OPTIONAL string, indicating the side for the variable annotations:
only set this to override default}
}
\value{
HeatmapAnnotation object
}
\description{
Helper to specify a HeatmapAnnotation for variables in cor_heatmap
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- tax_filter(ibd_phylo, min_prevalence = 5)
psq <- tax_mutate(psq, Species = NULL)
psq <- tax_fix(psq)
psq <- tax_agg(psq, rank = "Family")
taxa <- tax_top(psq, n = 15, rank = "Family")
samples <- phyloseq::sample_names(psq)

set.seed(42) # random colours used in first example
# sampleAnnotation returns a function that takes data, samples, and which
fun <- sampleAnnotation(
  gap = grid::unit(2.5, "mm"),
  Dis1 = anno_sample(var = "DiseaseState"),
  IBD = anno_sample_cat(var = "ibd"),
  Dis2 = anno_sample_cat(var = "DiseaseState", col = 1:4)
)

# manually specify the sample annotation function by giving it data etc.
heatmapAnnoFunction <- fun(.data = psq, .side = "top", .samples = samples)

# draw the annotation without a heatmap, you will never normally do this!
grid.newpage()
vp <- viewport(width = 0.65, height = 0.75)
pushViewport(vp)
draw(heatmapAnnoFunction)
pushViewport(viewport(x = 0.7, y = 0.6))
draw(attr(heatmapAnnoFunction, "Legends"))
}
\seealso{
\code{\link[=taxAnnotation]{taxAnnotation()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_labels.R
\name{textAngleCalc}
\alias{textAngleCalc}
\alias{textHjustCalc}
\title{Helpers for ord_plot label adjustments}
\usage{
textAngleCalc(xvec, yvec, max = 90, ratio = 1, perpendicular = FALSE)

textHjustCalc(xvec, adjust = TRUE)
}
\arguments{
\item{xvec}{numeric vector of values used for x axis}

\item{yvec}{numeric vector of values used for y axis}

\item{max}{maximum absolute numeric value of angle in degrees to return
(for rotating text/labels)}

\item{ratio}{adjustment for aspect ratio of plot when setting a fixed coordinate aspect
ratio with coord_fixed (advised)}

\item{adjust}{logical, apply hjust or not (FALSE means return only 0.5)}
}
\value{
numeric vector representing either angles to rotate geom_text
labels, or hjust values
}
\description{
Consider moving these functions to tax_lab_style() man page/.R file.

See functions section.
}
\section{Functions}{
\itemize{
\item \code{textAngleCalc}: Calculate rotation of text labels for ordination plot

\item \code{textHjustCalc}: Calculate hjust of text labels for ordination plot
}}

\examples{
library(ggplot2)
library(dplyr)

# create basic ggplot for labelling

df <- mtcars \%>\% mutate(across(everything(), scale))

p <- ggplot(df, aes(mpg, hp, label = rownames(df))) +
  geom_segment(xend = 0, yend = 0, color = "lightgrey") +
  annotate(x = 0, y = 0, geom = "point", size = 4) +
  theme_minimal()

p

# calculate new variable within aes mapping non-standard evaluation
p +
  geom_text(size = 2.5, mapping = aes(angle = textAngleCalc(mpg, hp))) +
  coord_fixed(ratio = 1)

# equivalent: calculate variable outside aes by referring to dataframe again
p +
  geom_text(size = 2.5, angle = textAngleCalc(df$mpg, df$hp)) +
  coord_fixed(ratio = 1)

# fixing aspect ratio is important
# see how angles may be incorrect otherwise
p +
  geom_text(size = 2.5, mapping = aes(angle = textAngleCalc(mpg, hp)))

# ratio argument allows matching angles with alternative aspect ratio
p +
  geom_text(size = 2.5, angle = textAngleCalc(df$mpg, df$hp, ratio = .5)) +
  coord_fixed(ratio = .5)

p +
  geom_text(size = 2.5, angle = textAngleCalc(df$mpg, df$hp, ratio = 1.5)) +
  coord_fixed(ratio = 1.5)

# perpendicular argument makes text perpendicular instead of parallel
p +
  geom_text(
    check_overlap = TRUE, size = 2.5,
    angle = textAngleCalc(df$mpg, df$hp, perpendicular = TRUE, ratio = 1.5)
  ) +
  coord_fixed(ratio = 1.5, clip = "off")

# max angle limits extreme text angles
p +
  geom_text(
    size = 2.5, check_overlap = TRUE,
    angle = textAngleCalc(df$mpg, df$hp, ratio = 2, max = 10),
    hjust = textHjustCalc(xvec = df$mpg, adjust = TRUE)
  ) +
  coord_fixed(ratio = 2, clip = "off")
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_agg.R
\name{tax_agg}
\alias{tax_agg}
\title{Aggregate taxa and create ps_extra}
\usage{
tax_agg(
  ps,
  rank = NA,
  sort_by = NA,
  top_N = NA,
  force = FALSE,
  add_unique = FALSE
)
}
\arguments{
\item{ps}{phyloseq object}

\item{rank}{NA (for tax_names level) or name of valid taxonomic rank (try phyloseq::rank_names(ps)) or "unique"}

\item{sort_by}{if not NA, how should the taxa be sorted, uses tax_sort(), takes same options as \code{by} arg}

\item{top_N}{NA does nothing, but if top_N is a number, it creates an extra tax_table column called top,
which is the same as the unique column for the first top_N number of taxa, and "other" otherwise.}

\item{force}{If TRUE, this forces aggregation at chosen rank to occur regardless of if the output will be sensible!
This avoids the "Taxa not unique at rank: ..." error, but may allow very inappropriate aggregation to occur.
Do not use force = TRUE unless you know why you are doing this, and what the result will be.
If you are getting an error with force = FALSE, it is almost certainly better to examine the tax_table and fix the problem.
force = TRUE is similar to microbiome::aggregate_taxa,
which also does not check that the taxa are uniquely defined by only the aggregation level.}

\item{add_unique}{if TRUE, adds a rank named unique, identical to the rownames after aggregation}
}
\value{
ps_extra list object including phyloseq and tax_agg rank info
}
\description{
\code{tax_agg} sums the abundances of the phyloseq taxa at the given rank.
It records the tax_agg rank argument in the info of the ps_extra object output.
This ps_extra object tracks aggregation, and any further transformations and scaling,
to help you keep track of what you have done with your phyloseq object and automatically caption ordination plots.

Instead of tax_agg, consider using \code{tax_transform()} with a rank argument instead, to both aggregate and transform the taxa.
This is also useful when you want to aggregate but not transform the taxa,
and yet still log the "identity" transformation in ps_extra for captioning your ordination plots.
e.g. \code{tax_transform(rank = "Genus", trans = "identity")}

tax_agg allows you to pass NA or "unique" to the rank argument which will NOT aggregate the taxa.
If you use rank = "unique" or add_unique = TRUE, it will add a new rank called unique, identical to the taxa_names (after any aggregation)

Be aware: you should not use the top_N argument yourself without good reason.
top_N provides a feature inspired by the deprecated microbiome function aggregate_top_taxa
which is primarily useful for decluttering compositional barplots.
microViz comp_barplot (and ord_plot_iris) already run tax_agg with a top_N argument for you, so you should not.
The tax_table produced when using top_N is otherwise INVALID FOR MOST OTHER ANALYSES.
}
\details{
This function is inspired by \code{microbiome::aggregate_taxa}.
However if \code{microbiome::aggregate_taxa} is used, microViz cannot track this aggregation.

Comparing aggregate_taxa and tax_agg:

Except for the ordering of taxa, and the addition of a "unique" rank being optional,
the resulting phyloseq objects are identical for aggregating a phyloseq with no ambiguous taxa.
Taxa are ambiguous when the tax_table converges at a lower rank after branching,
such as if two different genera share the same species (e.g. "s__").
\code{microbiome::aggregate_taxa} handles ambiguous taxa by creating a "unique" rank with all
of the taxonomic rank info pasted together into one, often very long, name.
\code{tax_agg} throws an error, and directs the user to \code{tax_fix()} to fix the ambiguous taxa before aggregation,
which should then result in (much) shorter unique names at the aggregation rank.
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")

tax_agg(ps = dietswap, "Phylum") \%>\%
  ps_get() \%>\%
  tax_table()
tax_agg(ps = dietswap, "Family") \%>\%
  ps_get() \%>\%
  tax_table()

# create some missing values
tax_table(dietswap)[3:7, "Genus"] <- "g__"

# this will produce an error, instructing the user to use tax_fix
# tax_agg(ps = dietswap, "Genus")

# this will then work:
dietswap \%>\%
  tax_fix() \%>\%
  tax_agg("Genus")

# you can replace unknown values with `tax_fix()`
# which will fix most problems, like the common "g__" and "s__"
# but default tax_fix settings won't catch this long unknown
tax_table(dietswap)[13:17, "Family"] <- "some_unknown_family"
dietswap \%>\%
  tax_fix(unknowns = "some_unknown_family") \%>\%
  tax_agg("Family")

# try tax_fix_interactive() to help you find and fix all the uninformative
# and converging values in your taxonomy table.

# the code below won't aggregate taxa,
# but just adds a new rank called unique, equal to taxa_names
tax_agg(ps = dietswap, rank = NA, add_unique = TRUE)
identical(tax_agg(dietswap, NA, add_unique = TRUE), tax_agg(dietswap, "unique")) # TRUE
}
\seealso{
\code{\link{tax_fix}}

\code{\link{tax_fix_interactive}}

\code{\link{tax_transform}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_models.R
\name{taxatree_models}
\alias{taxatree_models}
\title{Statistical modelling for individual taxa across multiple ranks}
\usage{
taxatree_models(
  ps,
  ranks = NULL,
  type = "lm",
  variables = NULL,
  formula = NULL,
  verbose = "rank",
  ...
)
}
\arguments{
\item{ps}{phyloseq object or ps_extra}

\item{ranks}{vector of rank names at which to aggregate taxa for modelling}

\item{type}{name of regression modelling function, or the function itself}

\item{variables}{vector of variable names, to be used as model formula right hand side
(ignored if formula given)}

\item{formula}{(alternative to variables arg) right hand side of a formula,
as a formula object or character value}

\item{verbose}{message about progress}

\item{...}{extra arguments are passed directly to modelling function}
}
\description{
\code{taxatree_models} runs \code{tax_model} on every taxon at multiple taxonomic
ranks (you choose which ranks with the plural \code{ranks} argument).
It returns the results as a named nested list of models
attached to a ps_extra object.
One list per rank, one model per taxon at each rank.

The result can then be used with \code{taxatree_models2stats} to extract a
dataframe of statistics for use with \code{taxatree_plots}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_funs.R
\name{taxatree_funs}
\alias{taxatree_funs}
\alias{taxatree_nodes}
\alias{taxatree_edges}
\title{Create node and edge dataframes for taxatree_plots}
\usage{
taxatree_nodes(
  ps,
  fun = list(sum = sum),
  ranks = "all",
  .sort = NULL,
  .use_counts = TRUE
)

taxatree_edges(nodes_df)
}
\arguments{
\item{ps}{phyloseq object or ps_extra}

\item{fun}{function to calculate for each taxon/node}

\item{ranks}{selection of taxonomic ranks to make nodes for ("all", or names)}

\item{.sort}{sort nodes by "ascending" or "descending" values of fun function result}

\item{.use_counts}{use count data if available (instead of transformed data)}

\item{nodes_df}{dataframe output from taxatree_nodes}
}
\description{
Mostly you will not have to use these functions directly:
instead call \code{taxatree_plots} with the output of \code{taxatree_stats}
\itemize{
\item \code{taxatree_nodes} creates taxon nodes and calculates a summary statistic
about each taxon (given by \code{fun}). Takes a ps_extra or phyloseq object.
\item \code{taxatree_edges} uses the output of \code{taxatree_nodes} to create a
dataframe of edges.
}
}
\details{
\code{taxatree_nodes} makes nodes for taxa at all ranks or for a list of
consecutive ranks (plus a root rank if tree is not rooted).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distinct_palette.R
\name{distinct_palette}
\alias{distinct_palette}
\title{Colour palettes suitable for 20+ categories}
\usage{
distinct_palette(n = NA, pal = "brewerPlus", add = "lightgrey")
}
\arguments{
\item{n}{number of colours to return}

\item{pal}{palette name, one of "brewerPlus", "kelly", "greenArmytage"}

\item{add}{colour to append to end of palette, as colour n+1,
lightgrey by default for the use as "other" taxa in comp_barplot,
or NA for no additional colour.}
}
\value{
vector of colours
}
\description{
Available palettes (max colors) are "brewerPlus" (41), "kelly" (20) and "greenArmytage" (25).
\itemize{
\item "brewerPlus" is an arbitrary expansion of the "Paired" and "Dark2"
colorbrewer palettes. The philosophy behind this expansion was to ensure
that similar colours are far apart, and the earlier colours are attractive.
\item "kelly" is based on the 22-colour palette developed by Kenneth Kelly but
with white and black starting colours removed. This palette is ordered
such that the first colours are most distinct.
\item "greenArmytage" is based on a 26-colour palette proposed by Paul
Green-Armytage, with black removed. This palette is not ordered by maximum
contrast.
}
}
\details{
Hex color codes for 'kelly' and 'greenArmytage' palettes are copied and
slightly modified from the Polychrome R package:
i.e. Polychrome::kelly.colors() and Polychrome::green.armytage.colors()

Please consider also citing Coombes 2019
\doi{10.18637/jss.v090.c01}
if you use either of these palettes.

See the Polychrome reference manual for more information:
\url{https://CRAN.R-project.org/package=Polychrome}
}
\examples{
brewerPlus <- distinct_palette()
scales::show_col(brewerPlus)

kelly <- distinct_palette(pal = "kelly")
scales::show_col(kelly)

greenArmytage <- distinct_palette(pal = "greenArmytage")
scales::show_col(greenArmytage)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_top.R
\name{tax_top}
\alias{tax_top}
\title{Print names of "top" n taxa}
\usage{
tax_top(data, n = 10, by = sum, rank = "unique", ...)
}
\arguments{
\item{data}{phyloseq object or ps_extra}

\item{n}{how many taxa names to return, or NA for all
(can return fewer than n values, if there are fewer to return)}

\item{by}{how to sort taxa (see \code{?tax_sort()}),
defaults to \code{sum} which sorts by total abundance across all samples}

\item{rank}{taxonomic rank to aggregate at before calculating
("unique" = no aggregation)}

\item{...}{extra optional args passed to tax_sort}
}
\value{
vector of taxa names at chosen rank
}
\description{
Simple wrapper function that:
\enumerate{
\item optionally aggregates taxa at \code{rank}
\item sorts the aggregated taxa according to \code{by}
\item returns the top \code{n} number of taxa names
}
}
\examples{
data("dietswap", package = "microbiome")
tax_top(dietswap)
tax_top(dietswap, n = 4, by = "prev", rank = "Phylum", undetected = 30)
}
\seealso{
\code{\link{tax_agg}} for more info on taxonomic aggregation

\code{\link{tax_sort}} for more info on sorting taxa
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_model.R
\name{tax_model}
\alias{tax_model}
\title{Statistical modelling for individual taxa in a phyloseq}
\usage{
tax_model(
  ps,
  rank,
  type = "lm",
  variables = NULL,
  formula = NULL,
  taxa = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{ps}{phyloseq object}

\item{rank}{name of taxonomic rank to aggregate to and model taxa at}

\item{type}{name of modelling function to use, or the function itself}

\item{variables}{vector of variable names to use in statistical model as right hand side (ignored if formula given)}

\item{formula}{(alternative to variables arg) right hand side of a formula, as a formula object or character value}

\item{taxa}{taxa to model (named, numbered, logical selection, or defaulting to all if NULL)}

\item{verbose}{message about progress and any taxa name modifications}

\item{...}{extra args passed directly to modelling function}
}
\value{
named list of model objects or list of lists
}
\description{
\code{tax_model} provides a simple framework to statistically model the abundance
of individual taxa in your data.
You can choose which type of statistical model you want to fit, and you can
choose at which rank and (optionally) which specific taxa to fit statistical models for.
\code{tax_model} takes a phyloseq and returns a list of statistical models, one model for each taxon.
The same independent variables are used for all models,
as specified in \code{variables} or \code{formula} argument (latter takes precedence).

\code{taxatree_models} runs \code{tax_model} on every taxon at multiple taxonomic ranks
(you choose which ranks with the plural \code{ranks} argument),
and returns the results as a named nested list designed for use with \code{taxatree_plots}.
One list per rank, one model per taxon at each rank.

\code{type = "bbdml"} will run beta binomial regression model(s) using the \code{corncob} package.
For bbdml the same formula/variables is/are used for modelling both the
abundance and dispersion parameters.
}
\details{
\code{tax_model} and \code{taxatree_models} can use parallel processing with the \code{future} package.
This can speed up analysis if you have many taxa to model.
Run a line like this beforehand: \code{future::plan(future::multisession, workers = 3)}
}
\examples{
library(corncob)
library(dplyr)

data(dietswap, package = "microbiome")
ps <- dietswap

# create some binary variables for easy visualisation
ps <- ps \%>\% ps_mutate(
  female = if_else(sex == "female", 1, 0, NaN),
  overweight = if_else(bmi_group == "overweight", 1, 0, NaN),
  obese = if_else(bmi_group == "obese", 1, 0, NaN)
)

# This example HITChip dataset has some taxa with the same name for phylum and family...
# We can fix problems like this with the tax_prepend_ranks function
ps <- tax_prepend_ranks(ps)

# filter out rare taxa (it is often difficult to fit multivariable models to rare taxa)
ps <- ps \%>\% tax_filter(min_prevalence = 0.1, min_total_abundance = 10000)

# specify variables used for modelling
VARS <- c("female", "overweight", "obese")

# Model first 3 genera using all VARS as predictors (just for a quick test)
models <- tax_model(ps, type = "bbdml", rank = "Genus", taxa = 1:3, variables = VARS)

# Alternative method using formula arg instead of variables to produce identical results
models2 <- tax_model(
  ps = ps, rank = "Genus", type = "bbdml",
  taxa = 1:3, formula = ~ female + overweight + obese
)
all.equal(models, models2) # should be TRUE

# Model only one genus, NOTE the modified name,
# which was returned by tax_prepend_ranks defaults
models3 <- ps \%>\%
  tax_model(
    rank = "Genus", type = "bbdml",
    taxa = "G: Bacteroides fragilis et rel.", variables = VARS
  )

# Model all taxa at multiple taxonomic ranks (ranks 1 and 2)
# using only female variable as predictor
models4 <- taxatree_models(
  ps = ps, type = "bbdml", ranks = 1:2, formula = ~female, verbose = FALSE
)

# modelling proportions with simple linear regression is also possible via type = lm
# and transforming the taxa to compositional first
models_lm <- ps \%>\%
  tax_transform("compositional") \%>\%
  tax_model(rank = "Genus", taxa = 1:3, variables = VARS, type = "lm")
}
\seealso{
\code{\link{taxatree_plots}} for how to plot the output of \code{taxatree_models}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_sort.R
\name{tax_reorder}
\alias{tax_reorder}
\title{Reorder taxa in phyloseq object using vector of names}
\usage{
tax_reorder(ps, tax_order, tree_warn = TRUE)
}
\arguments{
\item{ps}{phyloseq object}

\item{tax_order}{names or current numerical indices of taxa
in desired order and same length as taxa_names(ps)}

\item{tree_warn}{If phylogenetic tree is present in phyloseq phy_tree slot, taxa cannot be reordered.
Default behaviour of tax_sort is to remove the phylogenetic tree and warn about this.
tree_warn = FALSE will suppress the warning message, but still remove the tree!}
}
\value{
phyloseq object (always without phy_tree)
}
\description{
Reorder taxa in phyloseq object using vector of names
}
\examples{
data("dietswap", package = "microbiome")
new_order <- c(
  "Fusobacteria", "Cyanobacteria", "Verrucomicrobia", "Spirochaetes",
  "Actinobacteria", "Firmicutes", "Proteobacteria", "Bacteroidetes"
)
tax_agg(dietswap, rank = "Phylum")[["ps"]] \%>\%
  phyloseq::taxa_names()
tax_agg(dietswap, rank = "Phylum")[["ps"]] \%>\%
  microViz:::tax_reorder(tax_order = new_order) \%>\%
  phyloseq::taxa_names()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_melt.R
\name{ps_melt}
\alias{ps_melt}
\title{Melt phyloseq data object into large data.frame (tibble)}
\usage{
ps_melt(ps)
}
\arguments{
\item{ps}{(Required). An \code{\link{otu_table-class}} or
\code{\link{phyloseq-class}}. Function most useful for phyloseq-class.}
}
\value{
A \code{\link{tibble}} class data frame.
}
\description{
The ps_melt function is a specialized melt function for melting phyloseq objects
(instances of the phyloseq class), usually for producing graphics
with \code{\link[ggplot2]{ggplot}2}.
The naming conventions used in downstream phyloseq graphics functions
have reserved the following variable names that should not be used
as the names of \code{\link{sample_variables}}
or taxonomic \code{\link{rank_names}}.
These reserved names are \code{c("Sample", "Abundance", "OTU")}.
Also, you should not have identical names for
sample variables and taxonomic ranks.
That is, the intersection of the output of the following two functions
\code{\link{sample_variables}}, \code{\link{rank_names}}
should be an empty vector
(e.g. \code{intersect(sample_variables(ps), rank_names(ps))}).
All of these potential name collisions are checked-for
and renamed automatically with a warning.
However, if you (re)name your variables accordingly ahead of time,
it will reduce confusion and eliminate the warnings.

NOTE: Code and documentation copied (and very slightly modified) from an old version of speedyseq by Michael McLaren.
speedyseq reimplements \code{psmelt} from \code{phyloseq} to use functions from the \code{tidyr}
and \code{dplyr} packages. The name in microViz is changed to \code{ps_melt} for consistency with other functions.
}
\details{
Note that ``melted'' phyloseq data is stored much less efficiently, and so
RAM storage issues could arise with a smaller dataset (smaller number of
samples/OTUs/variables) than one might otherwise expect.  For common sizes
of graphics-ready datasets, however, this should not be a problem.  Because
the number of OTU entries has a large effect on the RAM requirement, methods
to reduce the number of separate OTU entries -- for instance by
agglomerating OTUs based on phylogenetic distance using
\code{\link{tip_glom}} -- can help alleviate RAM usage problems.  This
function is made user-accessible for flexibility, but is also used
extensively by plot functions in phyloseq.
}
\examples{
library(ggplot2)
library(phyloseq)
data("GlobalPatterns")
gp_ch <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
mdf <- ps_melt(gp_ch)
mdf2 <- psmelt(gp_ch) # slower
# same dataframe, except with somewhat different row orders
dplyr::all_equal(tibble::as_tibble(mdf), mdf2, convert = TRUE) # TRUE
nrow(mdf2)
ncol(mdf)
colnames(mdf)
head(rownames(mdf))
p <- ggplot(mdf, aes(x = SampleType, y = Abundance, fill = Genus))
p <- p + geom_bar(color = "black", stat = "identity", position = "stack")
# This example plot doesn't make any sense
print(p + coord_flip())
# TODO replace this...
}
\seealso{
\code{\link{psmelt}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-sampleAnnotation.R
\name{anno_cat_legend}
\alias{anno_cat_legend}
\title{Convenience function for generating levels of convenience.}
\usage{
anno_cat_legend(col, x = NULL, renamer = identity, title = "", ...)
}
\arguments{
\item{col}{vector of colors, named by all levels of data (e.g. x) or not named}

\item{x}{optional: vector of data to pair with unnamed col or check against named col}

\item{renamer}{function applied to generate labels: from names(col) or levels of x}

\item{title}{title of legend}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:Legend]{ComplexHeatmap::Legend}}
  \describe{
    \item{\code{labels}}{Labels corresponding to \code{at}. If it is not specified, the values of \code{at} are taken as labels.}
    \item{\code{nrow}}{For legend which is represented as grids, \code{nrow} controls number of rows of the grids if the grids are arranged into multiple rows.}
    \item{\code{ncol}}{Similar as \code{nrow}, \code{ncol} controls number of columns of the grids if the grids are arranged into multiple columns. Note at a same time only one of \code{nrow} and \code{ncol} can be specified.}
    \item{\code{by_row}}{Are the legend grids arranged by rows or by columns?}
    \item{\code{grid_height}}{The height of legend grid. It can also control the height of the continuous legend if it is horizontal.}
    \item{\code{grid_width}}{The width of legend grid. It can also control the width of the continuous legend if it is vertical.}
    \item{\code{gap}}{If legend grids are put into multiple rows or columns, this controls the gap between neighbouring rows or columns, measured as a \code{\link[grid]{unit}} object.}
    \item{\code{labels_gp}}{Graphic parameters for labels.}
    \item{\code{labels_rot}}{Text rotation for labels. It should only be used for horizontal continuous legend.}
    \item{\code{border}}{Color of legend grid borders. It also works for the ticks in the continuous legend.}
    \item{\code{type}}{Type of legends. The value can be one of \code{grid}, \code{points}, \code{lines} and \code{boxplot}.}
    \item{\code{direction}}{Direction of the legend, vertical or horizontal?}
    \item{\code{title_position}}{Position of title relative to the legend. \code{topleft}, \code{topcenter}, \code{leftcenter-rot} and \code{lefttop-rot} are only for vertical legend and \code{leftcenter}, \code{lefttop} are only for  horizontal legend.}
    \item{\code{title_gap}}{Gap between title and the legend body.}
  }}
}
\value{
a ComplexHeatmap Legend class object
}
\description{
Convenience function for generating levels of convenience.
}
\examples{
grid::grid.newpage()
ComplexHeatmap::draw(
  anno_cat_legend(
    col = c("ibd" = "blue", "nonibd" = "grey90"),
    renamer = toupper, title = "Hi there, I'm a title"
  )
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_permanova.R
\name{dist_permanova}
\alias{dist_permanova}
\title{Calculate PERMANOVA after dist_calc()}
\usage{
dist_permanova(
  data,
  variables = NULL,
  interactions = NULL,
  complete_cases = TRUE,
  n_processes = 1,
  n_perms = 999,
  seed = NULL,
  by = "margin",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{ps_extra output from dist_calc()}

\item{variables}{character vector of variables to include in model or character representation of the right-hand side of a formula, e.g "varA + varB + varA:varB"}

\item{interactions}{optional argument to define any interactions between variables, written in the style of e.g. "var_a * var_b"}

\item{complete_cases}{if TRUE, drops observations if they contain missing values (otherwise stops if missings are detected)}

\item{n_processes}{how many parallel processes to use? (on windows this uses parallel::makePSOCKcluster())}

\item{n_perms}{how many permutations? e.g. 9999. Less is faster but more is better!}

\item{seed}{set a random number generator seed to ensure you get the same results each run}

\item{by}{passed to vegan::adonis2() \code{by} argument: what type of sums of squares to calculate? "margin" or "terms"}

\item{verbose}{sends messages about progress if TRUE}

\item{...}{additional arguments are passed directly to vegan::adonis2() (e.g. strata, add, sqrt.dist etc.)}
}
\value{
ps_extra list containing permanova results and (filtered) input objects
}
\description{
\code{dist_permanova} runs a Permutational Multivariate ANOVA (aka Non-parametric MANOVA).
This is a way to test for the statistical significance of (independent)
associations between variables in your phyloseq::sample_data(),
and a microbiota distance matrix you have already calculated with dist_calc().

This function is a wrapper around vegan's \code{adonis2()} function. See \code{?vegan::adonis2()} for more insight.

You can also read this excellent book chapter on PERMANOVA by Marti Anderson:
\doi{10.1002/9781118445112.stat07841}

Or this NPMANOVA page on GUSTA ME:
\url{https://sites.google.com/site/mb3gustame/hypothesis-tests/manova/npmanova}
}
\details{
The variables argument will be collapsed into one string (if length > 1) by pasting together, separated by "+".
Any interaction terms described in the interactions argument will be pasted onto the end of the pasted variables argument.
Alternatively, you can supply the complete right hand side of the formula yourself e.g variables = "varA + varB + varC\*varD"

Watch out, if any of your variable names contain characters that would normally separate variables in a formula then
you should rename the offending variable (e.g. avoid any of "+" "\*" "|" or ":" ) otherwise permanova will split that variable into pieces.
}
\examples{
data("dietswap", package = "microbiome")

# add some missings to demonstrate automated removal
phyloseq::sample_data(dietswap)$sex[3:6] <- NA

# compute distance
testDist <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  tax_transform("identity") \%>\%
  dist_calc("bray")

PERM <- testDist \%>\%
  dist_permanova(
    seed = 1,
    variables = c("sex", "bmi_group"),
    n_processes = 1,
    n_perms = 99 # only 99 perms used in examples for speed (use 9999+!)
  )
PERM
str(PERM, max.level = 1)

# try permanova with interaction terms
PERM2 <- testDist \%>\%
  dist_permanova(
    seed = 1,
    variables = "nationality + sex * bmi_group",
    n_processes = 1, n_perms = 99
  )
PERM2$permanova

# specify the same model in alternative way
PERM3 <- testDist \%>\%
  dist_permanova(
    seed = 1,
    variables = c("nationality", "sex", "bmi_group"),
    interactions = "sex * bmi_group",
    n_processes = 1, n_perms = 99
  )
PERM3$permanova

identical(PERM3, PERM2) # TRUE

# take same distance matrix used for the permanova and plot an ordination
PERM2 \%>\%
  ord_calc(method = "PCoA") \%>\%
  ord_plot(color = "bmi_group")
# this trick ensures any samples dropped from the permanova
# for having missing values in the covariates are NOT included
# in the corresponding ordination plot
}
\seealso{
\code{\link{dist_calc}} for calculating the required distance matrix input

\code{\link{ord_plot}} with constraints as a way to visualise the microbial associations of significant predictors

\code{vegan::\link[vegan]{adonis2}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_join.R
\name{ps_join}
\alias{ps_join}
\title{Join a dataframe to phyloseq sample data}
\usage{
ps_join(
  x,
  y,
  by = NULL,
  match_sample_names = NULL,
  keep_sample_name_col = TRUE,
  sample_name_natural_join = FALSE,
  type = "left",
  .keep_all_taxa = FALSE
)
}
\arguments{
\item{x}{phyloseq (or dataframe)}

\item{y}{dataframe (or phyloseq for e.g. type = "right")}

\item{by}{A character vector of variables to join by (col must be present in both x and y or paired via a named vector like c("xname" = "yname", etc.))}

\item{match_sample_names}{match against the phyloseq sample_names by naming a variable in the additional dataframe (this is in addition to any variables named in by)}

\item{keep_sample_name_col}{should the column named in match_sample_names be kept in the returned phyloseq's sample_data? (only relevant if match_sample_names is not NULL)}

\item{sample_name_natural_join}{if TRUE, use sample_name AND all shared colnames to match rows (only relevant if match_sample_names is not NULL, this arg takes precedence over anything also entered in \code{by} arg)}

\item{type}{name of type of join e.g. "left", "right", "inner", "semi" (see dplyr help pages)}

\item{.keep_all_taxa}{if FALSE (the default), remove taxa which are no longer present in the dataset after filtering}
}
\value{
phyloseq with modified sample_data (and possibly filtered)
}
\description{
You can use most types of join from the dplyr::*_join function family, including e.g. "inner", "left", "semi", "anti" (see details below).
Defaults to type = "left" which calls left_join(), this supports x as a phyloseq and y as a dataframe.
Most of the time you'll want "left" (adds variables with no sample filtering), or "inner" (adds variables and filters samples).
This function simply:
\enumerate{
\item extracts the sample_data from the phyloseq as a dataframe
\item performs the chosen type of join (with the given arguments)
\item filters the phyloseq if type = inner, semi or anti
\item reattaches the modified sample_data to the phyloseq and returns the phyloseq
}
}
\details{
\strong{Mutating joins}, which will add columns from a dataframe to phyloseq sample data, matching rows based on the key columns named in the \code{by} argument:
\itemize{
\item "inner": includes all rows in present in both x and y.
\item "left": includes all rows in x. (so x must be the phyloseq)
\item "right": includes all rows in y. (so y must be the phyloseq)
\item "full": includes all rows present in x or y. (will likely NOT work, as additional rows cannot be added to sample_data!)
}

If a row in x matches multiple rows in y (based on variables named in the \code{by} argument),
all the rows in y will be added once for each matching row in x.
This will cause this function to fail, as additional rows cannot be added to the phyloseq sample_data!

\strong{Filtering joins} filter rows from x based on the presence or absence of matches in y:
\itemize{
\item "semi": return all rows from x with a match in y.
\item "anti": return all rows from x without a match in y.
}
}
\examples{
library(phyloseq)
data("enterotype", package = "phyloseq")

x <- enterotype
y <- data.frame(
  ID_var = sample_names(enterotype)[c(1:50, 101:150)],
  SeqTech = sample_data(enterotype)[c(1:50, 101:150), "SeqTech"],
  arbitrary_info = rep(c("A", "B"), 50)
)

# simply match the new data to samples that exist in x, as default is a left_join
# where some sample names of x are expected to match variable ID_var in dataframe y
out1A <- ps_join(x = x, y = y, match_sample_names = "ID_var")
out1A
sample_data(out1A)[1:6, ]


# use sample_name and all shared variables to join
# (a natural join is not a type of join per se,
# but it indicates that all shared variables should be used for matching)
out1B <- ps_join(
  x = x, y = y, match_sample_names = "ID_var",
  sample_name_natural_join = TRUE, keep_sample_name_col = FALSE
)
out1B
sample_data(out1B)[1:6, ]

# if you only want to keep phyloseq samples that exist in the new data, try an inner join
# this will add the new variables AND filter the phyloseq
# this example matches sample names to ID_var and by matching the shared SeqTech variable
out1C <- ps_join(x = x, y = y, type = "inner", by = "SeqTech", match_sample_names = "ID_var")
out1C
sample_data(out1C)[1:6, ]

# the id variable is named Sample_ID in x and ID_var in y
# semi_join is only a filtering join (doesn't add new variables but just filters samples in x)
out2A <- ps_join(x = x, y = y, by = c("Sample_ID" = "ID_var"), type = "semi")
out2A
sample_data(out2A)[1:6, ]

# anti_join is another type of filtering join
out2B <- ps_join(x = x, y = y, by = c("Sample_ID" = "ID_var"), type = "anti")
out2B
sample_data(out2B)[1:6, ]

# semi and anti joins keep opposite sets of samples
intersect(sample_names(out2A), sample_names(out2B))

# you can mix and match named and unnamed values in the `by` vector
# inner is like a combination of left join and semi join
out3 <- ps_join(x = x, y = y, by = c("Sample_ID" = "ID_var", "SeqTech"), type = "inner")
out3
sample_data(out3)[1:6, ]
}
\seealso{
\code{\link{ps_mutate}} for computing new variables from existing sample data

\code{\link{ps_select}} for selecting only some sample_data variables

\url{https://www.garrickadenbuie.com/project/tidyexplain/} for an animated introduction to joining dataframes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_mutate.R
\name{ps_mutate}
\alias{ps_mutate}
\title{Modify or compute new sample_data in phyloseq object}
\usage{
ps_mutate(ps, ..., .target)
}
\arguments{
\item{ps}{phyloseq object with sample data}

\item{...}{passed straight to dplyr::mutate (see examples and dplyr::mutate help)}

\item{.target}{DEPRECATED. See tax_mutate for manipulation of tax_table}
}
\value{
phyloseq object with modified sample_data
}
\description{
Add or compute new phyloseq sample_data variables.
Uses \code{dplyr::mutate()} syntax.
}
\examples{
library(phyloseq)
library(dplyr)
data("enterotype")

sample_data(enterotype)[1:10, ]

months_in_year <- 12
ps <- enterotype \%>\%
  ps_mutate(
    age_months = Age * months_in_year,
    IDs_match = as.character(Sample_ID) == as.character(SampleID),
    placeholder = "Word"
  )

sample_data(ps)[1:10, ]

# Using the dplyr::across functionality is also possible
ps <- ps \%>\%
  ps_mutate(
    dplyr::across(where(is.factor), toupper),
    another_var = TRUE,
    SeqTech = NULL # deletes SeqTech variable
  )

head(sample_data(ps))
}
\seealso{
\code{\link{tax_mutate}} for manipulation of tax_table variables

\code{\link[dplyr]{mutate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{heat_numbers}
\alias{heat_numbers}
\title{Aesthetic settings for drawing numbers on heatmap tiles}
\usage{
heat_numbers(
  decimals = 0,
  fontsize = 7,
  col = "darkgrey",
  fontface = "bold",
  fmt = NULL,
  ...
)
}
\arguments{
\item{decimals}{number of decimal places to print}

\item{fontsize}{fontsize specification,}

\item{col}{colour of font}

\item{fontface}{plain, bold, italic}

\item{fmt}{NULL or number print format, see ?sprintf, overrides decimals arg if set}

\item{...}{passed to grid::gpar() for grid.text}
}
\value{
list
}
\description{
Works with comp_heatmap() and cor_heatmap().
See the help for those functions.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_labels.R
\name{Ordination-labels}
\alias{Ordination-labels}
\alias{tax_lab_style}
\alias{constraint_lab_style}
\title{Create list for ord_plot() *_lab_style arguments}
\usage{
tax_lab_style(
  type = "label",
  max_angle = 0,
  perpendicular = FALSE,
  aspect_ratio = 1,
  justify = "auto",
  size = 2,
  alpha = 1,
  colour = "black",
  ...
)

constraint_lab_style(
  type = "label",
  max_angle = 0,
  perpendicular = FALSE,
  aspect_ratio = 1,
  justify = "auto",
  size = 2.5,
  alpha = 1,
  colour = "brown",
  ...
)
}
\arguments{
\item{type}{'label', 'text' or 'richtext'
('richtext' also used if 'label' type are rotated, when max_angle > 0)}

\item{max_angle}{maximum angle of rotation to allow to match vector angle
(requires ggtext package to rotate "label" type)}

\item{perpendicular}{if TRUE, sets rotated labels perpendicular to desired angle, not parallel}

\item{aspect_ratio}{aspect ratio of plot (y/x) must also be used in coord_fixed() ratio argument
(must be set when rotated labels are used, to ensure match to arrow angles)}

\item{justify}{"center", "side", or "auto"?
Should the text/label align with the arrows at the text center or text sides
(uses hjust, if 'auto', picks based on whether max_angle is greater than 0)}

\item{size}{fixed size of text or label}

\item{alpha}{fixed alpha of text or label}

\item{colour}{fixed colour of text or label}

\item{...}{further named arguments passed to geom_text, geom_label or geom_richtext}
}
\value{
named list
}
\description{
Customise taxa and constraint labels on your ordination plots.
Choose 'text' or 'label' type, rotate and/or justify the text/labels
and set aesthetic appearances using \code{tax_lab_style()} or
\code{constraint_lab_style()}.
}
\examples{
# These examples show styling of taxa labels with tax_lab_style().
# The same options are available for constraint labels in constrained
# ordinations. constraint_lab_style() just has different default settings.

library(ggplot2)

# get example inflammatory bowel disease stool dataset from corncob package
data("ibd_phylo", package = "corncob")

# filter out rare taxa and clean up names etc
ibd <- ibd_phylo \%>\%
  tax_filter(min_prevalence = 3) \%>\%
  tax_fix() \%>\%
  phyloseq_validate()

# calculate a centered-log-ratio transformed PCA ordination
ibd_ord <- ibd \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc("PCA")

# basic plot with default label style
ibd_ord \%>\% ord_plot(color = "ibd", plot_taxa = 1:10)

# Rotating labels: requires the ggtext package #
# A fixed coordinate ratio must be set to ensure label rotation
# matches the vectors. It is also helpful to set the vector and label length
# multipliers manually for a good look. Rotated labels are justified to the
# 'sides' automatically by tax_lab_style() with justify = 'auto'
ibd_ord \%>\%
  ord_plot(
    color = "ibd", plot_taxa = 1:7,
    tax_vec_length = 1.3, tax_lab_length = 1.3,
    tax_lab_style = tax_lab_style(max_angle = 90)
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

# You can use text instead of labels
# - a bold fontface helps text to stand out
# - see ?ggplot2::geom_text for all settings available
ibd_ord \%>\%
  ord_plot(
    color = "ibd", plot_taxa = 1:7,
    tax_vec_length = 1.3, tax_lab_length = 1.4,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, size = 2.5, fontface = "bold.italic"
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

# With text you can prevent overlaps with check_overlap = TRUE
ibd_ord \%>\%
  ord_plot(
    color = "ibd", plot_taxa = 1:12,
    tax_vec_length = 1.3, tax_lab_length = 1.4,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, size = 3, fontface = "bold.italic",
      check_overlap = TRUE
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

# With labels, you can reduce the padding and line weight to free space
# but check_overlap is not available
# see ?ggtext::geom_richtext for more possibilities
ibd_ord \%>\%
  ord_plot(
    color = "ibd", plot_taxa = 1:7,
    tax_vec_length = 1.3, tax_lab_length = 1.35,
    tax_lab_style = tax_lab_style(
      max_angle = 90, fontface = "italic", size = 2.5, fill = "grey95",
      label.size = 0.1, # width outline
      label.padding = unit(0.1, "lines"),
      label.r = unit(0, "lines") # reduces rounding of corners to radius 0
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

# Perpendicular angled labels/text are possible
ibd_ord \%>\%
  ord_plot(
    color = "ibd", plot_taxa = 1:12,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, perpendicular = TRUE, size = 3,
      check_overlap = TRUE
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))


# You can limit and/or attenuate the angle of rotation by:
#  - setting a lower max_angle
#  - decreasing the aspect_ratio in the tax_lab_style call
ibd_ord \%>\%
  ord_plot(
    shape = "circle", color = "ibd", plot_taxa = 1:7,
    tax_vec_length = 1.3, tax_lab_length = 1.3,
    tax_lab_style = tax_lab_style(
      max_angle = 10, size = 2, label.size = 0.1,
      label.padding = unit(0.1, "lines"), label.r = unit(0, "lines")
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

ibd_ord \%>\%
  ord_plot(
    shape = "circle", color = "ibd", plot_taxa = 1:7,
    tax_vec_length = 1.3, tax_lab_length = 1.3,
    tax_lab_style = tax_lab_style(
      max_angle = 90, size = 2, label.size = 0.1, aspect_ratio = 0.5,
      label.padding = unit(0.1, "lines"), label.r = unit(0, "lines")
    )
  ) +
  coord_fixed(ratio = 1, clip = "off", xlim = c(-3.5, 3.5))

# another example with some extras #
ibd_ord \%>\%
  ord_plot(
    shape = "circle filled", fill = "ibd",
    plot_taxa = 1:10,
    taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "),
    tax_vec_length = 2, tax_lab_length = 2.1,
    tax_lab_style = tax_lab_style(
      type = "text", max_angle = 90, size = 2.5,
      fontface = "bold.italic", check_overlap = TRUE
    )
  ) +
  coord_fixed(1, clip = "off", xlim = c(-5, 5)) +
  theme(legend.position = c(0.8, 0.2), legend.background = element_rect()) +
  stat_chull(mapping = aes(colour = ibd, fill = ibd), alpha = 0.1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_sort_ord.R
\name{ps_sort_ord}
\alias{ps_sort_ord}
\alias{ord_order_samples}
\title{Sort phyloseq samples by ordination axes scores}
\usage{
ps_sort_ord(ps, ord, axes = 1:2, scaling = 2)

ord_order_samples(ord, axes = 1:2, scaling = 2)
}
\arguments{
\item{ps}{phyloseq object to be sorted}

\item{ord}{ps_extra with ordination object}

\item{axes}{which axes to use for sorting? numerical vector of length 1 or 2}

\item{scaling}{Type 2, or type 1 scaling. For more info, see \url{https://sites.google.com/site/mb3gustame/constrained-analyses/rda}.
Either "species" or "site" scores are scaled by (proportional) eigenvalues, and the other set of scores is left unscaled (from ?vegan::scores.cca)}
}
\value{
\code{ps_sort_ord} returns a phyloseq

\code{ord_order_samples} returns a character vector
}
\description{
\code{ps_sort_ord} reorders samples in a phyloseq object based on their relative
position on 1 or 2 ordination axes.

\code{ord_order_samples} gets the sample_names in order from the ordination
contained in a ps_extra list. This is used internally by \code{ps_sort_ord}

If 2 axes given, the samples are sorted by anticlockwise rotation around the
selected ordination axes, starting on first axis given, upper right quadrant.
(This is used by ord_plot_iris.)

If 1 axis is given, samples are sorted by increasing value order along this axis.
This could be used to arrange samples on a rectangular barplot in order of
appearance along a parallel axis of a paired ordination plot.
}
\examples{
# attach other necessary packages
library(ggplot2)

# example data
ibd <- corncob::ibd_phylo \%>\%
  tax_filter(min_prevalence = 2) \%>\%
  tax_fix() \%>\%
  phyloseq_validate()

# create numeric variables for constrained ordination
ibd <- ibd \%>\%
  ps_mutate(
    ibd = as.numeric(ibd == "ibd"),
    steroids = as.numeric(steroids == "steroids"),
    abx = as.numeric(abx == "abx"),
    female = as.numeric(gender == "female"),
    # and make a shorter ID variable
    id = stringr::str_remove_all(sample, "^[0]{1,2}|-[A-Z]+")
  )

# create an ordination
ordi <- ibd \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc()

ord_order_samples(ordi, axes = 1) \%>\% head(8)
ps_sort_ord(ibd, ordi, axes = 1) \%>\%
  phyloseq::sample_names() \%>\%
  head(8)

p1 <- ord_plot(ordi, colour = "grey90", plot_taxa = 1:8, tax_vec_length = 1) +
  geom_text(aes(label = id), size = 2.5, colour = "red")

b1 <- ibd \%>\%
  ps_sort_ord(ord = ordi, axes = 1) \%>\%
  comp_barplot(
    tax_level = "Genus", n_taxa = 12, label = "id",
    order_taxa = ord_order_taxa(ordi, axes = 1),
    sample_order = "default",
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

patchwork::wrap_plots(p1, b1, ncol = 1)

# constrained ordination example (and match vertical axis) #
cordi <- ibd \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc(
    constraints = c("steroids", "abx", "ibd"), conditions = "female",
    scale_cc = FALSE
  )

cordi \%>\% ord_plot(plot_taxa = 1:6, axes = 2:1)
}
\seealso{
\itemize{
\item These functions were created to support ordering of samples on \code{ord_plot_iris}
\item \code{tax_sort_ord} for ordering taxa in phyloseq by ordination
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_prepend_ranks.R
\name{tax_prepend_ranks}
\alias{tax_prepend_ranks}
\title{Add rank prefixes to phyloseq tax_table values}
\usage{
tax_prepend_ranks(ps, sep = ": ", nchar = 1)
}
\arguments{
\item{ps}{phyloseq object}

\item{sep}{characters to paste in between rank initial and taxon name}

\item{nchar}{number of characters to use from start of rank_names}
}
\value{
phyloseq
}
\description{
Prepend the start of rank names to each taxon at each rank
(useful particularly in case of duplicated taxa names across ranks, e.g. dietswap dataset)
}
\examples{
data("dietswap", package = "microbiome")
phyloseq::tax_table(dietswap) \%>\% head()
dietswap \%>\%
  tax_prepend_ranks() \%>\%
  phyloseq::tax_table() \%>\%
  head()
}
\seealso{
\code{\link{tax_fix}} for fixing other tax_table problems
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{comp_heatmap}
\alias{comp_heatmap}
\title{Draw heatmap of microbiome composition across samples}
\usage{
comp_heatmap(
  data,
  taxa = NA,
  taxa_side = "right",
  tax_anno = NULL,
  taxon_renamer = identity,
  samples = NA,
  sample_side = adjacent_side(taxa_side),
  sample_anno = NULL,
  sample_names_show = FALSE,
  colors = heat_palette(palette = "Rocket", rev = TRUE),
  numbers = NULL,
  sample_seriation = "OLO_ward",
  sample_ser_dist = "euclidean",
  sample_ser_counts = !sample_ser_dist \%in\% c("euclidean", "maximum", "manhattan",
    "canberra", "binary"),
  sample_ser_trans = NULL,
  tax_seriation = "OLO_ward",
  tax_ser_dist = "euclidean",
  tax_ser_counts = FALSE,
  tax_ser_trans = NULL,
  numbers_trans = NULL,
  numbers_zero_replace = 0,
  numbers_use_counts = TRUE,
  grid_col = "white",
  grid_lwd = 0.5,
  name = "Abd.",
  anno_tax = NULL,
  ...
)
}
\arguments{
\item{data}{phyloseq or phyloseq extra}

\item{taxa}{list of taxa to include, or NA for all}

\item{taxa_side}{"top"/"right"/"bottom"/"left": controls heatmap orientation and where any
annotations specified in tax_anno are placed}

\item{tax_anno}{NULL or annotation function for taxa: taxAnnotation() output.}

\item{taxon_renamer}{function to rename taxa before plotting}

\item{samples}{list of samples to include on plot}

\item{sample_side}{which side to show any sample annotation on, must be adjacent to taxa_side}

\item{sample_anno}{NULL or annotation function for samples: sampleAnnotation() output.}

\item{sample_names_show}{show sample names? (you can control side and rotation of names with
other ComplexHeatmap::Heatmap arguments)}

\item{colors}{output of heat_palette() to set heatmap fill color scheme}

\item{numbers}{output of heat_numbers() to draw numbers on heatmap cells}

\item{sample_seriation}{name of method used to order the samples (from seriation::seriate)}

\item{sample_ser_dist}{name of distance to use with sample_seriation method (if needed)}

\item{sample_ser_counts}{insist on using count data for sample seriation?}

\item{sample_ser_trans}{function for transformation of data used for sample seriation
(such as a call to \code{tax_transform()})}

\item{tax_seriation}{name of method used to order the taxa (from seriation::seriate)}

\item{tax_ser_dist}{name of distance to use with tax_seriation method (if needed)}

\item{tax_ser_counts}{insist on using count data for taxa seriation?}

\item{tax_ser_trans}{function for transformation of data used for taxa seriation
(such as a call to \code{tax_transform()})}

\item{numbers_trans}{name of tax_transform transformation, or a function
for transformation of data used for drawing numbers on cells}

\item{numbers_zero_replace}{zero replacement method used if named transformation given to number_trans}

\item{numbers_use_counts}{insist on using count data for number drawing?
(if TRUE, any numbers_trans transformation would be applied to count data)}

\item{grid_col}{colour of gridlines, or NA for none}

\item{grid_lwd}{width of gridlines}

\item{name}{used as legend title (colourbar)}

\item{anno_tax}{DEPRECATED:
optional annotation of taxa distributions: tax_anno() list output,
or a pre-made ComplexHeatmap HeatmapAnnotation}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap}}
  \describe{
    \item{\code{row_dend_side}}{Should the row dendrogram be put on the left or right of the heatmap?}
    \item{\code{row_dend_width}}{Width of the row dendrogram, should be a \code{\link[grid]{unit}} object.}
    \item{\code{show_row_dend}}{Whether show row dendrogram?}
    \item{\code{row_dend_gp}}{Graphic parameters for the dendrogram segments. If users already provide a \code{\link[stats]{dendrogram}} object with edges rendered, this argument will be ignored.}
    \item{\code{show_row_names}}{Whether show row names.}
    \item{\code{row_names_gp}}{Graphic parameters for row names.}
    \item{\code{row_names_rot}}{Rotation of row names.}
    \item{\code{row_names_centered}}{Should row names put centered?}
  }}
}
\description{
Heatmap made with \code{ComplexHeatmap::Heatmap()},
with optional annotation of taxa prevalence/abundance,
and/or other sample data.

Transform your data with \code{tax_transform()} prior to plotting
(and/or scale with \code{tax_scale()}).

See the heatmaps vignette for more examples of use.

Plotting "compositional" data can give an idea of the dominant taxa in each sample.
Plotting some form of log or clr transformed (or scaled) microbial features can highlight
other patterns.

The data will be ordered via your selected seriation methods and distances
on either the transformed data (default) or the original count data
(or with any other transformation).

Any cell numbers printed can be transformed independently of the colour scheme,
and do not affect ordering.
}
\examples{
library(dplyr)
data("dietswap", package = "microbiome")
# create a couple of numerical variables to use
psq <- dietswap \%>\%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  )
psq <- tax_filter(psq, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq <- tax_agg(psq, "Genus")

# randomly select 20 taxa from the 40 top taxa, and 40 random samples

set.seed(123)
taxa <- sample(tax_top(psq, n = 40), size = 20)
samples <- sample(1:122, size = 40)

comp_heatmap(data = psq, taxa = taxa, samples = samples)

# transforming taxon abundances #

# NOTE: if you plan on transforming taxa (e.g. to compositional data or clr)
# but only want to plot a subset of the taxa (e.g. most abundant)
# you should NOT subset the original phyloseq before transformation!
# Instead, choose the subset of taxa plotted with:

# Note 2, choose a symmetrical palette for clr-transformed data
psq \%>\%
  tax_transform("clr", zero_replace = "halfmin") \%>\%
  comp_heatmap(
    taxa = taxa, samples = samples, colors = heat_palette(sym = TRUE)
  )

# Almost all the taxa have high values (>> 0) because they are a highly
# abundant subset taken after clr transformation was calculated on all taxa

# See how just taking the first 30 taxa from the dataset gives more balance
psq \%>\%
  tax_transform("clr", zero_replace = "halfmin") \%>\%
  comp_heatmap(
    taxa = 1:30, samples = samples, colors = heat_palette(sym = TRUE)
  )

# annotating taxa #

# Notes:
# - Unlike cor_heatmap, taxa are not annotated by default
# - Detection threshold set to 50 as HITchip example data seems to have background noise

comp_heatmap(
  data = psq, taxa = taxa, samples = samples,
  tax_anno = taxAnnotation(Prev = anno_tax_prev(undetected = 50))
)

# annotating samples #

htmp <- psq \%>\%
  tax_transform("clr", zero_replace = "halfmin") \%>\%
  comp_heatmap(
    taxa = taxa, samples = samples, colors = heat_palette(sym = TRUE),
    sample_anno = sampleAnnotation(
      Nation. = anno_sample_cat("nationality", legend_title = "Nation.")
    )
  )
htmp

# legends from `anno_sample_cat()` are stored as an attribute of the Heatmap
ComplexHeatmap::draw(
  object = htmp,
  annotation_legend_list = attr(htmp, "AnnoLegends"), merge_legends = TRUE
)
}
\seealso{
\code{\link[=cor_heatmap]{cor_heatmap()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_plot.R
\name{ord_plot}
\alias{ord_plot}
\title{Customisable ggplot2 ordination plot}
\usage{
ord_plot(
  data,
  axes = 1:2,
  plot_taxa = FALSE,
  tax_vec_length = NA,
  tax_vec_style_all = vec_tax_all(),
  tax_vec_style_sel = vec_tax_sel(),
  tax_lab_length = tax_vec_length * 1.1,
  tax_lab_style = list(),
  taxon_renamer = function(x) identity(x),
  constraint_vec_length = NA,
  constraint_vec_style = vec_constraint(),
  constraint_lab_length = constraint_vec_length * 1.1,
  constraint_lab_style = list(),
  var_renamer = function(x) identity(x),
  plot_samples = TRUE,
  scaling = 2,
  auto_caption = 8,
  center = FALSE,
  clip = "off",
  expand = !center,
  interactive = FALSE,
  ...
)
}
\arguments{
\item{data}{ps_extra list object, output from ord_calc}

\item{axes}{which axes to plot: numerical vector of length 2, e.g. 1:2 or c(3,5)}

\item{plot_taxa}{if ord_calc method was "PCA/RDA" draw the taxa loading vectors (see details)}

\item{tax_vec_length}{taxon arrow vector scale multiplier.
NA = auto-scaling, or provide a numeric multiplier yourself.}

\item{tax_vec_style_all}{list of named aesthetic attributes for all (background) taxon vectors}

\item{tax_vec_style_sel}{list of named aesthetic attributes for taxon vectors for the taxa selected by plot_taxa}

\item{tax_lab_length}{scale multiplier for label distance/position for any selected taxa}

\item{tax_lab_style}{list of style options for the taxon labels, see tax_lab_style() function.}

\item{taxon_renamer}{function that takes any plotted taxon names and returns modified names for labels}

\item{constraint_vec_length}{constraint arrow vector scale multiplier.
NA = auto-scaling, or provide a numeric multiplier yourself.}

\item{constraint_vec_style}{list of aesthetics/arguments (colour, alpha etc) for the constraint vectors}

\item{constraint_lab_length}{label distance/position for any constraints
(relative to default position which is proportional to correlations with each axis)}

\item{constraint_lab_style}{list of aesthetics/arguments (colour, size etc) for the constraint labels}

\item{var_renamer}{function to rename constraining variables for plotting their labels}

\item{plot_samples}{if TRUE, plot sample points with geom_point}

\item{scaling}{Type 2, or type 1 scaling. For more info, see \url{https://sites.google.com/site/mb3gustame/constrained-analyses/rda}.
Either "species" or "site" scores are scaled by (proportional) eigenvalues, and the other set of scores is left unscaled (from ?vegan::scores.cca)}

\item{auto_caption}{size of caption with info about the ordination, NA for none}

\item{center}{expand plot limits to center around origin point (0,0)}

\item{clip}{clipping of labels that extend outside plot limits?}

\item{expand}{expand plot limits a little bit further than data range?}

\item{interactive}{creates plot suitable for use with ggiraph (used in ord_explore)}

\item{...}{pass aesthetics arguments for sample points,
drawn with geom_point using aes_string}
}
\value{
ggplot
}
\description{
Draw ordination plot. Utilises results of \code{\link{ord_calc}}.
\itemize{
\item For an extensive tutorial see \href{https://david-barnett.github.io/microViz/articles/web-only/ordination.html}{the ordination vignette}.
\item For interpretation see the the relevant pages on PCA, PCoA, RDA, or CCA on the GUide to STatistical Analysis in Microbial Ecology (GUSTA ME) website: \url{https://sites.google.com/site/mb3gustame/}
}
}
\details{
How to specify the plot_taxa argument (when using PCA, CCA or RDA):
\itemize{
\item FALSE --> plot no taxa vectors or labels
\item integer vector e.g. 1:3 --> plot labels for top 3 taxa (by longest line length)
\item single numeric value e.g. 0.75 --> plot labels for taxa with line length > 0.75
\item character vector e.g. c('g__Bacteroides', 'g__Veillonella') --> plot labels for the exactly named taxa
}
}
\examples{
library(ggplot2)
data("dietswap", package = "microbiome")

# create a couple of numerical variables to use as constraints or conditions
dietswap <- dietswap \%>\%
  ps_mutate(
    weight = dplyr::recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = dplyr::if_else(sex == "female", true = 1, false = 0)
  )

# unconstrained PCA ordination
unconstrained_aitchison_pca <- dietswap \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc() # method = "auto" --> picks PCA as no constraints or distances

unconstrained_aitchison_pca \%>\%
  ord_plot(colour = "bmi_group", plot_taxa = 1:5) +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group))

# you can generate an interactive version of the plot by specifying
# interactive = TRUE, and passing a variable name to another argument
# called `data_id` which is required for interactive point selection
interactive_plot <- unconstrained_aitchison_pca \%>\%
  ord_plot(
    colour = "bmi_group", plot_taxa = 1:5,
    interactive = TRUE, data_id = "sample"
  )

# to start the html viewer, and allow selecting points, we must use a
# ggiraph function called girafe and set some options and css
ggiraph::girafe(
  ggobj = interactive_plot,
  options = list(
    ggiraph::opts_selection(
      css = ggiraph::girafe_css(
        css = "fill:orange;stroke:black;",
        point = "stroke-width:1.5px"
      ),
      type = "multiple", # this activates lasso selection (click top-right)
      only_shiny = FALSE # allows interactive plot outside of shiny app
    )
  )
)


# remove effect of weight with conditions arg
# scaling weight with scale_cc is not necessary as only 1 condition is used
dietswap \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc(conditions = "weight", scale_cc = FALSE) \%>\%
  ord_plot(colour = "bmi_group") +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group))

# alternatively, constrain variation on weight and female
constrained_aitchison_rda <- dietswap \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc(constraints = c("weight", "female")) # constraints --> RDA

constrained_aitchison_rda \%>\%
  ord_plot(colour = "bmi_group", constraint_vec_length = 2) +
  stat_ellipse(aes(linetype = bmi_group, colour = bmi_group))

# ggplot allows additional customisation of the resulting plot
p <- constrained_aitchison_rda \%>\%
  ord_plot(colour = "bmi_group", plot_taxa = 1:3) +
  lims(x = c(-5, 6), y = c(-5, 5)) +
  scale_colour_brewer(palette = "Set1")

p + stat_ellipse(aes(linetype = bmi_group, colour = bmi_group))
p + stat_density2d(aes(colour = bmi_group))

# you can rename the taxa on the labels with any function that
# takes and modifies a character vector
constrained_aitchison_rda \%>\%
  ord_plot(
    colour = "bmi_group",
    plot_taxa = 1:3,
    taxon_renamer = function(x) stringr::str_extract(x, "^.")
  ) +
  lims(x = c(-5, 6), y = c(-5, 5)) +
  scale_colour_brewer(palette = "Set1")

# You can plot PCoA and constrained PCoA plots too.
# You don't typically need/want to use transformed taxa variables for PCoA
# But it is good practice to call tax_transform("identity") so that
# the automatic caption can record that no transformation was applied
dietswap \%>\%
  tax_agg("Genus") \%>\%
  tax_transform("identity") \%>\%
  # so caption can record (lack of) transform
  dist_calc("bray") \%>\%
  # bray curtis
  ord_calc() \%>\%
  # guesses you want an unconstrained PCoA
  ord_plot(colour = "bmi_group")

# it is possible to facet these plots
# (although I'm not sure it makes sense to)
# but only unconstrained ordination plots and with plot_taxa = FALSE
unconstrained_aitchison_pca \%>\%
  ord_plot(color = "sex", auto_caption = NA) +
  facet_wrap("sex") +
  theme(line = element_blank()) +
  stat_density2d(aes(colour = sex)) +
  guides(colour = FALSE)

unconstrained_aitchison_pca \%>\%
  ord_plot(color = "bmi_group", plot_samples = FALSE, auto_caption = NA) +
  facet_wrap("sex") +
  theme(line = element_blank(), axis.text = element_blank()) +
  stat_density2d_filled(show.legend = FALSE) +
  geom_point(size = 1, shape = 21, colour = "black", fill = "white")
}
\seealso{
\code{\link{tax_lab_style}} / \code{\link{tax_lab_style}} for styling labels

\code{\link{ord_explore}} for interactive ordination plots

\code{\link{ord_calc}} for calculating an ordination to plot with ord_plot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_plots.R
\name{taxatree_plots}
\alias{taxatree_plots}
\title{Plot statistical model results for all taxa on a taxonomic tree}
\usage{
taxatree_plots(
  data,
  colour_stat = "estimate",
  palette = "Green-Brown",
  reverse_palette = FALSE,
  colour_lims = NULL,
  colour_oob = scales::oob_squish,
  colour_trans = "abs_sqrt",
  size_stat = list(prevalence = prev),
  node_size_range = c(1, 4),
  edge_width_range = node_size_range * 0.8,
  size_guide = "legend",
  size_trans = "identity",
  sig_stat = "p.value",
  sig_threshold = 0.05,
  sig_shape = "circle filled",
  sig_size = 0.75,
  sig_stroke = 0.75,
  sig_colour = "white",
  edge_alpha = 0.7,
  vars = "term",
  var_renamer = identity,
  title_size = 10,
  layout = "tree",
  layout_seed = NA,
  circular = identical(layout, "tree"),
  node_sort = NULL,
  add_circles = isTRUE(circular),
  drop_ranks = TRUE,
  l1 = if (palette == "Green-Brown") 10 else NULL,
  l2 = if (palette == "Green-Brown") 85 else NULL,
  colour_na = "grey35"
)
}
\arguments{
\item{data}{ps_extra with taxatree_stats, e.g. output of \code{taxatree_models2stats()}}

\item{colour_stat}{name of variable to scale colour/fill of nodes and edges}

\item{palette}{diverging hcl colour palette name from \code{colorspace::hcl_palettes("diverging")}}

\item{reverse_palette}{reverse direction of colour palette?}

\item{colour_lims}{limits of colour and fill scale, NULL will infer lims from range of all data}

\item{colour_oob}{scales function to handle colour_stat values outside of colour_lims
(default simply squishes "out of bounds" values into the given range)}

\item{colour_trans}{name of transformation for colour scale:
default is "abs_sqrt", the square-root of absolute values,
but you can use the name of any transformer from the \code{scales} package,
such as "identity" or "exp"}

\item{size_stat}{named list of length 1, giving function calculated for each taxon,
to determine the size of nodes (and edges). Name used as size legend title.}

\item{node_size_range}{min and max node sizes, decrease to avoid node overlap}

\item{edge_width_range}{min and max edge widths}

\item{size_guide}{guide for node sizes, try "none", "legend" or ggplot2::guide_legend()}

\item{size_trans}{transformation for size scale
you can use (the name of) any transformer from the scales package,
such as "identity", "log1p", or "sqrt"}

\item{sig_stat}{name of variable indicating statistical significance}

\item{sig_threshold}{value of sig_stat variable indicating statistical significance (below this)}

\item{sig_shape}{fixed shape for significance marker}

\item{sig_size}{fixed size for significance marker}

\item{sig_stroke}{fixed stroke width for significance marker}

\item{sig_colour}{fixed colour for significance marker (used as fill for filled shapes)}

\item{edge_alpha}{fixed alpha value for edges}

\item{vars}{name of column indicating terms in models (one plot made per term)}

\item{var_renamer}{function to rename variables for plot titles}

\item{title_size}{font size of title}

\item{layout}{any ggraph layout, default is "tree"}

\item{layout_seed}{any numeric, required if a stochastic igraph layout is named}

\item{circular}{should the layout be circular?}

\item{node_sort}{sort nodes by "increasing" or "decreasing" size? NULL for no sorting.
Use \code{tax_sort()} before \code{taxatree_plots()} for finer control.}

\item{add_circles}{add faint concentric circles to plot, behind each rank?}

\item{drop_ranks}{drop ranks from tree if not included in stats dataframe}

\item{l1}{Luminance value at the scale endpoints, NULL for palette's default}

\item{l2}{Luminance value at the scale midpoint, NULL for palette's default}

\item{colour_na}{colour for NA values in tree.
(if unused ranks are not dropped, they will have NA values for colour_stat)}
}
\value{
list of ggraph ggplots
}
\description{
\itemize{
\item Uses a ps_extra object to make a tree graph structure from the taxonomic table.
\item Then adds statistical results stored in "taxatree_stats" of ps_extra data
\item You must use \code{taxatree_models()} first to generate statistical model results.
\item You can adjust p-values with \code{taxatree_stats_p_adjust()}
}
}
\details{
\code{taxatree_plotkey} plots same layout as \code{taxatree_plots}, but in a fixed colour

See website article for more examples of use:
https://david-barnett.github.io/microViz/articles/web-only/modelling-taxa.html

Uses ggraph, see help for main underlying graphing function with \code{?ggraph::ggraph}

It is possible to provide multiple significance markers for multiple thresholds,
by passing vectors to the sig_shape, sig_threshold, etc. arguments.
It is critically important that the thresholds are provided in decreasing
order of severity, e.g. sig_threshold = c(0.001, 0.01, 0.1) and you must provide
a shape value for each of them.
}
\examples{
# Limited examples, see website article for more

library(dplyr)
library(ggplot2)

data(dietswap, package = "microbiome")
ps <- dietswap

# create some binary variables for easy visualisation
ps <- ps \%>\% ps_mutate(
  female = if_else(sex == "female", 1, 0, NaN),
  african = if_else(nationality == "AFR", 1, 0, NaN)
)

# This example dataset has some taxa with the same name for phylum and family...
# We can fix problems like this with the tax_prepend_ranks function
# (This will always happen with Actinobacteria!)
ps <- tax_prepend_ranks(ps)

# filter out rare taxa
ps <- ps \%>\% tax_filter(
  min_prevalence = 0.5, prev_detection_threshold = 100
)

# delete the Family rank as we will not use it for this small example
# this is necessary as taxatree_plots can only plot consecutive ranks
ps <- ps \%>\% tax_mutate(Family = NULL)

# specify variables used for modelling
models <- taxatree_models(
  ps = ps, type = corncob::bbdml, ranks = c("Phylum", "Genus"),
  formula = ~ female + african, verbose = TRUE
)
# models list stored as attachment in ps_extra
models

# get stats from models
stats <- taxatree_models2stats(models, param = "mu")
stats

plots <- taxatree_plots(
  data = stats, colour_trans = "identity",
  size_stat = list(mean = mean),
  size_guide = "legend", node_size_range = c(1, 6)
)

# if you change the size_stat for the plots, do the same for the key!!
key <- taxatree_plotkey(
  data = stats,
  rank == "Phylum" | p.value < 0.05, # labelling criteria
  .combine_label = all, # label only taxa where criteria met for both plots
  size_stat = list(mean = mean),
  node_size_range = c(2, 7), size_guide = "none",
  taxon_renamer = function(x) {
    stringr::str_remove_all(x, "[PG]: | [ae]t rel.")
  }
)

# cowplot is powerful for arranging trees and key and colourbar legend
legend <- cowplot::get_legend(plots[[1]])
plot_col <- cowplot::plot_grid(
  plots[[1]] + theme(legend.position = "none"),
  plots[[2]] + theme(legend.position = "none"),
  ncol = 1
)
cowplot::plot_grid(key, plot_col, legend, nrow = 1, rel_widths = c(4, 2, 1))
}
\seealso{
\code{\link[=taxatree_models]{taxatree_models()}} to calculate statistical models for each taxon

\code{\link[=taxatree_plotkey]{taxatree_plotkey()}} to plot the corresponding labelled key

\code{\link[=taxatree_plot_labels]{taxatree_plot_labels()}} and \code{\link[=taxatree_label]{taxatree_label()}} to add labels

\code{\link[=taxatree_stats_p_adjust]{taxatree_stats_p_adjust()}} to adjust p-values
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_paths.R
\name{add_paths}
\alias{add_paths}
\title{Add paths connecting points on a ggplot scatterplot}
\usage{
add_paths(
  ggplot,
  id_var,
  id_values,
  mapping = NULL,
  arrow = grid::arrow(length = grid::unit(2, units = "mm")),
  ...
)
}
\arguments{
\item{ggplot}{ggplot scatterplot such as the output of ord_plot}

\item{id_var}{name of variable used to identify grouping of points}

\item{id_values}{values of id_var variable used to identify groups of points to draw}

\item{mapping}{ggplot aesthetics created by aes(), e.g. aes(colour = ?) - group is already set to id_var internally!}

\item{arrow}{arrowhead to add to path, NULL for none}

\item{...}{additional arguments passed to geom_path}
}
\value{
ggplot with added layer with geom_path
}
\description{
Useful for tracing a few select individuals over time on an ordination plot.
Samples in phyloseq must be arranged in order of timepoint for the path connections to be drawn in the correct order!
You can arrange the samples in timepoint order with ps_arrange.
}
\examples{
library(ggplot2)
data("dietswap", package = "microbiome")

# arrange by timepoint first (or whatever your own time variable is)
dietswap \%>\%
  ps_arrange(timepoint) \%>\%
  tax_fix() \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc(method = "PCA") \%>\%
  ord_plot(colour = "timepoint", alpha = 0.5, size = 2) \%>\%
  add_paths(
    id_var = "subject", id_values = c("azl", "byn"),
    mapping = aes(colour = timepoint), size = 1.5
  )

# paths do NOT connect points in the correct order without arranging first
dietswap \%>\%
  tax_fix() \%>\%
  tax_transform("clr", rank = "Genus") \%>\%
  ord_calc(method = "PCA") \%>\%
  ord_plot(colour = "timepoint", alpha = 0.5) \%>\%
  add_paths(
    id_var = "subject", id_values = c("azl", "byn"),
    mapping = aes(colour = timepoint), size = 1.5
  ) +
  ggtitle("WRONG PATH ORDER", "use ps_arrange first!")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_models2stats.R
\name{taxatree_models2stats}
\alias{taxatree_models2stats}
\title{Extract dataframe of statistics from taxatree_models list}
\usage{
taxatree_models2stats(data, fun = "auto", ..., .keep_models = FALSE)
}
\arguments{
\item{data}{ps_extra with taxatree_models, or just the list of models}

\item{fun}{function to assist extraction of stats dataframe from models}

\item{...}{extra arguments passed to fun}

\item{.keep_models}{should the models list be kept in the ps_extra output?}
}
\value{
data.frame, attached to ps_extra
}
\description{
Extract dataframe of statistics from taxatree_models list
}
\examples{
# This example is an abbreviated excerpt from article on taxon modelling on
# the microViz documentation website

library(dplyr)
data("ibd_phylo", package = "corncob")

# We'll keep only the Ulcerative Colitis and Healthy Control samples, to
# simplify the analyses for this example. We'll also remove the Species
# rank information, as most OTUs in this dataset are not assigned to a
# species. We'll also use `tax_fix` to fill any gaps where the Genus is
# unknown, with the family name or whatever higher rank classification is
# known.

phylo <- ibd_phylo \%>\%
  ps_filter(DiseaseState \%in\% c("UC", "nonIBD")) \%>\%
  tax_mutate(Species = NULL) \%>\%
  tax_fix()

# Let's make some sample data variables that are easier to use and compare
# in the statistical modelling ahead. We will convert dichotomous
# categorical variables into similar binary variables (values: 1 for true,
# or 0 for false). We will also scale and center the numeric variable for
# age.

phylo <- phylo \%>\%
  ps_mutate(
    UC = ifelse(DiseaseState == "UC", yes = 1, no = 0),
    female = ifelse(gender == "female", yes = 1, no = 0),
    antibiotics = ifelse(abx == "abx", yes = 1, no = 0),
    steroids = ifelse(steroids == "steroids", yes = 1, no = 0),
    age_scaled = scale(age, center = TRUE, scale = TRUE)
  )

lm_models <- phylo \%>\%
  tax_fix() \%>\%
  tax_prepend_ranks() \%>\%
  tax_transform("compositional", rank = "Genus") \%>\%
  tax_filter(min_prevalence = 0.1, use_counts = TRUE) \%>\%
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) \%>\%
  taxatree_models(
    type = lm,
    ranks = c("Phylum", "Class", "Genus"),
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  )


lm_stats <- lm_models \%>\% taxatree_models2stats()

# inspect the ps_extra returned, now with taxatree_stats dataframe added
lm_stats

# inspect the dataframe itself
lm_stats$taxatree_stats

# keep the models on the ps_extra object
lm_models \%>\% taxatree_models2stats(.keep_models = TRUE)

# you can adjust the p values with taxatree_stats_p_adjust

# you can plot the results with taxatree_plots
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-taxAnnotation.R
\name{anno_tax_prev}
\alias{anno_tax_prev}
\title{Helper to specify heatmap annotation for showing taxa prevalence as barplot}
\usage{
anno_tax_prev(
  undetected = 0,
  use_counts = TRUE,
  size = grid::unit(20, "mm"),
  baseline = 0,
  border = TRUE,
  bar_width = 0.6,
  gp = grid::gpar(fill = "#CCCCCC"),
  ylim = NULL,
  extend = 0.05,
  axis = TRUE,
  ...,
  data = NULL,
  taxa = NULL,
  which = NULL
)
}
\arguments{
\item{undetected}{the value above which taxa are classed as detected/present in a sample}

\item{use_counts}{try to retrieve counts from data object?}

\item{size}{width or height as a grid unit object}

\item{baseline}{baseline of bars. The value should be "min" or "max", or a numeric value. It is enforced to be zero for stacked barplots.}

\item{border}{Wether draw borders of the annotation region?}

\item{bar_width}{Relative width of the bars. The value should be smaller than one.}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{ylim}{Data ranges. By default it is \code{range(x)} if \code{x} is a vector, or \code{range(rowSums(x))} if \code{x} is a matrix.}

\item{extend}{The extension to both side of \code{ylim}. The value is a percent value corresponding to \code{ylim[2] - ylim[1]}.}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_barplot]{ComplexHeatmap::anno_barplot}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{taxa}{OPTIONAL selection vector of taxa (names, numbers or logical),
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to taxAnnotation(),
which itself is used by cor_heatmap and comp_heatmap as tax_anno argument.
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- tax_filter(ibd_phylo, min_prevalence = 5)
psq <- tax_mutate(psq, Species = NULL)
psq <- tax_fix(psq)
psq <- tax_agg(psq, rank = "Family")
taxa <- tax_top(psq, n = 15, rank = "Family")

# makes a function that takes data, taxa and which (at minimum)
fun <- anno_tax_prev()

# manually specify the prevalence barplot function by giving it data etc.
heatmapAnnoFunction <- fun(data = psq, which = "row", taxa = taxa)

# draw the barplot without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)

grid::grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)

# let's change some style options and specify the data up front
grid::grid.newpage()
pushViewport(vp)
anno_tax_prev(
  data = psq, taxa = taxa, which = "column",
  gp = grid::gpar(fill = "red", lwd = 3, alpha = 0.5),
  border = FALSE, bar_width = 1
) \%>\%
  draw()

# clear drawings
grid::grid.newpage()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_extra-accessors.R
\name{ps_extra-accessors}
\alias{ps_extra-accessors}
\alias{ps_get}
\alias{dist_get}
\alias{ord_get}
\alias{info_get}
\alias{perm_get}
\alias{bdisp_get}
\alias{otu_get}
\alias{tt_get}
\alias{samdat_tbl}
\title{Extract elements from ps_extra class}
\usage{
ps_get(ps_extra)

dist_get(ps_extra)

ord_get(ps_extra)

info_get(ps_extra)

perm_get(ps_extra)

bdisp_get(ps_extra)

otu_get(data, taxa = NA, samples = NA, counts = FALSE)

tt_get(data)

samdat_tbl(data, sample_names_col = ".sample_name")
}
\arguments{
\item{ps_extra}{ps_extra class object}

\item{data}{phyloseq or ps_extra}

\item{taxa}{subset of taxa to return, NA for all (default)}

\item{samples}{subset of samples to return, NA for all (default)}

\item{counts}{should otu_get ensure it returns counts? if present in object}

\item{sample_names_col}{name of column where sample_names are put.
if NA, return data.frame with rownames (sample_names)}
}
\value{
element of ps_extra class object (or NULL)
}
\description{
\itemize{
\item \code{ps_get}     returns phyloseq
\item \code{info_get}   returns ps_extra_info object
\item \code{dist_get}   returns distance matrix (or NULL)
\item \code{ord_get}    returns ordination object (or NULL)
\item \code{perm_get}   returns adonis2() permanova model (or NULL)
\item \code{bdisp_get}  returns results of betadisper() (or NULL)
\item \code{otu_get}    returns phyloseq otu_table matrix with taxa as columns
\item \code{tt_get}     returns phyloseq tax_table
\item \code{samdat_tbl} returns phyloseq sample_data as a tibble,
with sample_names as new first column called .sample_name
}
}
\examples{
data("dietswap", package = "microbiome")
psx <- tax_transform(dietswap, "identity", rank = "Genus")
psx

ps_get(psx)
info_get(psx)

dist_get(psx) # this ps_extra has no dist_calc result
ord_get(psx) # this ps_extra has no ord_calc result
perm_get(psx) # this ps_extra has no dist_permanova result
bdisp_get(psx) # this ps_extra has no dist_bdisp result

# these can be returned from phyloseq objects too
otu_get(psx)[1:6, 1:4]
tt_get(psx) \%>\% head()
samdat_tbl(psx) \%>\% head()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_select.R
\name{ps_select}
\alias{ps_select}
\title{Select phyloseq sample_data using dplyr::select syntax}
\usage{
ps_select(ps, ...)
}
\arguments{
\item{ps}{phyloseq with sample_data}

\item{...}{passed straight to dplyr::select}
}
\value{
phyloseq object
}
\description{
Simple selection of phyloseq sample_data variables, might be useful for printing reduced sample_data, or inside other functions
}
\examples{
library(phyloseq)
library(dplyr)
data("enterotype", package = "phyloseq")

head(sample_data(enterotype))

enterotype \%>\%
  ps_select(!contains("Sample")) \%>\%
  sample_data() \%>\%
  head()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-varAnnotation.R
\name{anno_var_density}
\alias{anno_var_density}
\title{Helper to specify heatmap annotation for variable distribution density plot}
\usage{
anno_var_density(
  fun = identity,
  size = grid::unit(30, "mm"),
  type = c("lines", "violin", "heatmap"),
  xlim = NULL,
  heatmap_colors = c("white", "forestgreen"),
  joyplot_scale = 1.5,
  border = TRUE,
  gp = grid::gpar(fill = "lightgrey"),
  axis = TRUE,
  ...,
  data = NULL,
  vars = NULL,
  which = NULL
)
}
\arguments{
\item{fun}{function applied to all variables, with apply()}

\item{size}{width or height as a grid unit object}

\item{type}{Type of graphics to represent density distribution. "lines" for normal density plot; "violine" for violin plot and "heatmap" for heatmap visualization of density distribution.}

\item{xlim}{Range on x-axis.}

\item{heatmap_colors}{A vector of colors for interpolating density values.}

\item{joyplot_scale}{Relative height of density distribution. A value higher than 1 increases the height of the density distribution and the plot will represented as so-called "joyplot".}

\item{border}{Wether draw borders of the annotation region?}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_density]{ComplexHeatmap::anno_density}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{vars}{OPTIONAL selection vector of variable names,
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to varAnnotation(),
which itself is used by cor_heatmap var_anno argument.
}
\examples{
library(ComplexHeatmap)
set.seed(123)
fakeData <- as.data.frame.matrix(matrix(rnorm(500, 10, 3), ncol = 10))
names(fakeData) <- paste0("var_", 1:10)

# draw the plots without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(
  anno_var_density(data = fakeData, vars = names(fakeData), which = "row")
)

grid.newpage()
pushViewport(vp)
draw(
  anno_var_density(
    data = fakeData, fun = function(x) log(x + 1),
    vars = rev(names(fakeData)), type = "heatmap",
    which = "column"
  )
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq_validate.R
\name{phyloseq_validate}
\alias{phyloseq_validate}
\title{Check for (and fix) common problems with phyloseq objects}
\usage{
phyloseq_validate(
  ps,
  remove_undetected = FALSE,
  min_tax_length = 4,
  verbose = TRUE
)
}
\arguments{
\item{ps}{phyloseq object}

\item{remove_undetected}{if TRUE, removes taxa that sum to zero across all samples}

\item{min_tax_length}{minimum number of characters to not consider a tax_table entry suspiciously short}

\item{verbose}{print informative messages if true}
}
\value{
possibly modified phyloseq object
}
\description{
\itemize{
\item It checks for, and messages about, common uninformative entries in the tax_table, which often cause unwanted results
\item If there is no sample_data, it creates a sample_data dataframe with the sample_names (as "SAMPLE" variable)
\item If there is no tax_table, it creates a 1-column tax_table matrix with the taxa_names, and calls the rank "unique"
\item If remove_undetected = TRUE, it removes taxa where \code{phyloseq::taxa_sums()} is equal to zero, with a warning
}
}
\examples{
data(dietswap, package = "microbiome")

# expect warning about taxa summing to zero
phyloseq_validate(dietswap, remove_undetected = TRUE, verbose = TRUE)

# verbose = FALSE will suppress messages and warnings but still:
# replace NULL sample_data and remove taxa that sum to 0 across all samples
# (if remove_undetected = TRUE)
phyloseq_validate(dietswap, verbose = FALSE)

# Sometimes you might have a phyloseq with no sample_data
# This isn't compatible with some microViz functions, like comp_barplot
# So some functions internally use phyloseq_validate to fix this
dietswap@sam_data <- NULL
phyloseq_validate(dietswap)

# Sometimes you might have a phyloseq with no tax_table
# This isn't compatible with some microViz functions, like tax_top,
# so this is another reason to start your analyses with phyloseq_validate!
data("soilrep", package = "phyloseq")
soilrep # has NULL tax_table
phyloseq_validate(soilrep)

# If no messages or warnings are emitted,
# this means no problems were detected, and nothing was changed
# (but only if verbose = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_plot_iris.R
\name{ord_plot_iris}
\alias{ord_plot_iris}
\title{Circular compositional barplot sorted by ordination angle}
\usage{
ord_plot_iris(
  data,
  tax_level,
  axes = 1:2,
  n_taxa = 10,
  ord_plot = "none",
  taxon_renamer = function(x) identity(x),
  palette = distinct_palette(n_taxa),
  anno_colour = NULL,
  anno_colour_style = list(),
  anno_binary = NULL,
  anno_binary_style = list(),
  keep_all_vars = FALSE,
  scaling = 2,
  count_warn = TRUE,
  ...
)
}
\arguments{
\item{data}{ps_extra list output of ord_calc}

\item{tax_level}{taxonomic aggregation level (from rank_names(ps))}

\item{axes}{which 2 axes of ordination to use for ordering bars}

\item{n_taxa}{how many taxa to colour show distinct colours for (all other taxa grouped into "other").}

\item{ord_plot}{add a matching ordination plot to your iris plot
('list' returns separate plots in a list, 'above'/'below' uses patchwork to pair plots together into one)}

\item{taxon_renamer}{function to rename taxa in the legend}

\item{palette}{colour palette}

\item{anno_colour}{name of sample_data variable to use for colouring geom_segment annotation ring}

\item{anno_colour_style}{list of further arguments passed to geom_segment e.g. size}

\item{anno_binary}{name(s) of binary sample_data variable(s) (levels T/F or 1/0) to use for filtered geom_point annotation ring(s) (annotates at TRUE values)}

\item{anno_binary_style}{list of further arguments passed to geom_point e.g. colour, size, y, etc.}

\item{keep_all_vars}{slows down processing but is required for any post-hoc plot customisation options}

\item{scaling}{Type 2, or type 1 scaling. For more info, see \url{https://sites.google.com/site/mb3gustame/constrained-analyses/rda}.
Either "species" or "site" scores are scaled by (proportional) eigenvalues, and the other set of scores is left unscaled (from ?vegan::scores.cca)}

\item{count_warn}{warn if count data are not available? i.e. phyloseq otu_table is not positive integers and ps_extra counts slot is NULL}

\item{...}{extra args passed to comp_barplot e.g. bar_width}
}
\value{
ggplot
}
\description{
Use with \code{ord_calc} output as data argument.
Order of samples extracted from ordination axes in data.
Best paired with ordination plot made from same \code{ord_calc} output.
}
\details{
data must also contain counts table if taxa were transformed (e.g. for clr PCA ordination)
(i.e. you must have used \code{tax_transform} with keep_counts = TRUE, if transformation was not "identity")

You cannot set a variable fill aesthetic (only fixed) for the annotation points,
as the fill is used for the taxonomic composition bars
}
\examples{
library(dplyr)
library(ggplot2)
data("dietswap", package = "microbiome")

# although these iris plots are great for 100s of samples
# we'll take a subset of the data (for speed in this example)
ps <- dietswap \%>\%
  ps_filter(timepoint \%in\% c(1, 2)) \%>\%
  # copy an otu to the sample data
  ps_otu2samdat("Prevotella melaninogenica et rel.") \%>\%
  # create a couple of useful variables
  ps_mutate(
    female = sex == "female",
    african = nationality == "AFR",
    log_P.melaninogenica = log10(`Prevotella melaninogenica et rel.` + 1)
  )

# define a function for taking the end off the long genus names in this dataset
tax_renamer <- function(tax) {
  stringr::str_remove(tax, " [ae]t rel.")
}

ord <- ps \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("aitchison") \%>\%
  ord_calc(method = "PCoA")

# ordination plot for comparison
ord \%>\% ord_plot(color = "log_P.melaninogenica", size = 3)

ord_plot_iris(
  data = ord,
  tax_level = "Genus",
  n_taxa = 10,
  anno_colour = "nationality",
  anno_colour_style = list(size = 3),
  anno_binary = "female",
  anno_binary_style = list(shape = "F", size = 2.5),
  taxon_renamer = tax_renamer
) +
  scale_colour_brewer(palette = "Dark2")

# It is also possible to use comp_barplot customisation arguments
# like bar_width and bar_outline_colour, and to make interactive iris plots
# using ggiraph:

if (interactive()) {
  hover_over_me <- ord_plot_iris(
    data = ord,
    tax_level = "Genus",
    n_taxa = 10,
    anno_colour = "nationality",
    anno_colour_style = list(size = 3),
    anno_binary = "female",
    anno_binary_style = list(shape = "F", size = 2.5),
    taxon_renamer = tax_renamer,
    interactive = TRUE,
    bar_width = 0.8, bar_outline_colour = "black"
  ) +
    scale_colour_brewer(palette = "Dark2")

  ggiraph::girafe(ggobj = hover_over_me)
}

# Using PCA for ordination after transformations (e.g. clr) means the untransformed taxonomic
# data are only available for plotting as compositions if you transformed with
# tax_transform(keep_counts = TRUE) and your original data were in fact counts.
# Compositional data will also work, and you can set count_warn to FALSE to avoid the warning

clr_pca <- ps \%>\%
  tax_agg("Genus") \%>\%
  tax_transform("clr") \%>\%
  ord_calc(method = "PCA")

# you can generate a simple paired layout of ord_plot and iris plot
# or separately create and pair the plots yourself, for more control

# simple pairing
ord_plot_iris(
  data = clr_pca, n_taxa = 12,
  tax_level = "Genus",
  taxon_renamer = tax_renamer,
  ord_plot = "below",
  bar_width = 0.8, bar_outline_colour = "black",
  anno_binary = "african",
  anno_binary_style = list(
    y = 1.08, colour = "gray50", shape = "circle open", size = 1, stroke = 1.5
  )
)

# manual pairing
plot1 <- clr_pca \%>\% ord_plot(
  plot_taxa = 6:1, tax_vec_length = 0.6,
  colour = "gray50", shape = "nationality",
  taxon_renamer = tax_renamer,
  auto_caption = NA, center = TRUE,
) +
  scale_shape_manual(values = c(AFR = "circle", AAM = "circle open"))

iris <- ord_plot_iris(
  data = clr_pca, n_taxa = 15,
  tax_level = "Genus",
  taxon_renamer = tax_renamer,
  anno_binary = "african",
  anno_binary_style = list(y = 1.05, colour = "gray50", shape = "circle", size = 1)
) +
  # shrink legend text size
  theme(legend.text = element_text(size = 7))

cowplot::plot_grid(plot1, iris, nrow = 1, align = "h", axis = "b", rel_widths = 3:4)

# you can add multiple rings of binary annotations
ord_plot_iris(
  data = clr_pca, n_taxa = 15,
  tax_level = "Genus",
  taxon_renamer = tax_renamer,
  anno_binary = c("african", "female"),
  anno_binary_style = list(
    colour = c("gray50", "coral"),
    shape = c("circle", "F"), size = c(0.5, 2)
  )
) +
  theme(legend.text = element_text(size = 7))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-taxAnnotation.R
\name{taxAnnotation}
\alias{taxAnnotation}
\title{Helper to specify a HeatmapAnnotation for taxa}
\usage{
taxAnnotation(
  ...,
  name,
  annotation_legend_param = list(),
  show_legend = TRUE,
  gp = grid::gpar(col = NA),
  border = FALSE,
  gap = grid::unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_label = NULL,
  annotation_name_gp = grid::gpar(),
  annotation_name_offset = NULL,
  annotation_name_rot = NULL,
  annotation_name_align = TRUE,
  annotation_name_side = "auto",
  .data = NULL,
  .taxa = NULL,
  .side = NULL
)
}
\arguments{
\item{...}{Name-value pairs where the names correspond to annotation names and values
are the output of taxon annotation functions such as anno_tax_prev() or
manually specified AnnotationFunction objects}

\item{name}{Name of the heatmap annotation, optional.}

\item{annotation_legend_param}{A list which contains parameters for annotation legends. See \code{\link[ComplexHeatmap]{color_mapping_legend,ColorMapping-method}} for all possible options.}

\item{show_legend}{Whether show annotation legends. The value can be one single value or a vector.}

\item{gp}{Graphic parameters for simple annotations (with \code{fill} parameter ignored).}

\item{border}{border of single annotations.}

\item{gap}{Gap between annotations. It can be a single value or a vector of \code{\link[grid]{unit}} objects.}

\item{show_annotation_name}{Whether show annotation names? For column annotation, annotation names are drawn either on the left or the right, and for row annotations, names are draw either on top or at the bottom. The value can be a vector.}

\item{annotation_label}{Labels for the annotations. By default it is the same as individual annotation names.}

\item{annotation_name_gp}{Graphic parameters for anntation names. Graphic paramters can be vectors.}

\item{annotation_name_offset}{Offset to the annotation names, a \code{\link[grid]{unit}} object. The value can be a vector.}

\item{annotation_name_rot}{Rotation of the annotation names. The value can be a vector.}

\item{annotation_name_align}{Whether to align the annotation names.}

\item{annotation_name_side}{Side of the annotation names.}

\item{.data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{.taxa}{OPTIONAL selection vector of taxa (names, numbers or logical),
only set this if providing .data argument to override default}

\item{.side}{OPTIONAL string, indicating the side the taxa annotation should be placed:
only set this to override default}
}
\value{
HeatmapAnnotation object
}
\description{
Helper to specify a HeatmapAnnotation for taxa
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- tax_filter(ibd_phylo, min_prevalence = 5)
psq <- tax_mutate(psq, Species = NULL)
psq <- tax_fix(psq)
psq <- tax_agg(psq, rank = "Family")
taxa <- tax_top(psq, n = 15, rank = "Family")

customAxis <- list(labels_rot = 0, at = c(0, 0.5, 1))

# makes a function that takes data, taxa and which (at minimum)
fun <- taxAnnotation(
  gap = grid::unit(2.5, "mm"),
  Prev. = anno_tax_prev(axis_param = customAxis, ylim = c(0, 1), extend = 0),
  `Prop. Abd.` = anno_tax_box(size = unit(40, "mm"), axis_param = customAxis),
  `Log10p Abd.` = anno_tax_density(type = "heatmap")
)

# manually specify the prevalence barplot function by giving it data etc.
heatmapAnnoFunction <- fun(.data = psq, .side = "top", .taxa = taxa)

# draw the annotation without a heatmap, you will never normally do this!
grid.newpage()
vp <- viewport(width = 0.65, height = 0.75)
pushViewport(vp)
draw(heatmapAnnoFunction)

# try again as a row annotation
grid.newpage()
pushViewport(vp)
draw(fun(.data = psq, .side = "right", .taxa = rev(taxa)))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_calc_dominant.R
\name{ps_calc_dominant}
\alias{ps_calc_dominant}
\title{Calculate dominant taxon in each phyloseq sample}
\usage{
ps_calc_dominant(
  ps,
  rank,
  threshold = 0.3,
  n_max = 6,
  var = paste("dominant", rank, sep = "_"),
  none = "none",
  other = "other"
)
}
\arguments{
\item{ps}{phyloseq object}

\item{rank}{taxonomic rank to calculate dominance at}

\item{threshold}{minimum proportion at which to consider a sample dominated by a taxon}

\item{n_max}{maximum number of taxa that can be listed as dominant taxa}

\item{var}{name of variable to add to phyloseq object sample data}

\item{none}{character value to use when no taxon reaches threshold}

\item{other}{character value to use when another taxon (>n_max) dominates}
}
\value{
phyloseq object
}
\description{
Which taxon is most abundant in each sample of your phyloseq object?
This function adds this information as a new variable in your phyloseq
sample_data.
\itemize{
\item If the most abundant taxon is below the proportional abundance threshold,
the dominant taxon will be "none" for that sample
\item If there are more than n_max dominant taxa across all samples (not
including "none") the dominant taxon will be "other" for those samples
}
}
\details{
Thanks to Vitor Heidrich for the idea and a draft implementation
}
\examples{
library(ggplot2)
ps <- corncob::ibd_phylo \%>\%
  tax_filter(min_prevalence = 3) \%>\%
  tax_fix() \%>\%
  phyloseq_validate()

ps \%>\%
  ps_filter(DiseaseState == "CD") \%>\%
  ps_calc_dominant(rank = "Genus") \%>\%
  comp_barplot(tax_level = "Genus", label = "dominant_Genus", n_taxa = 12) +
  coord_flip()

ps \%>\%
  ps_calc_dominant(rank = "Genus") \%>\%
  tax_transform(rank = "Genus", trans = "clr") \%>\%
  ord_calc("PCA") \%>\%
  ord_plot(colour = "dominant_Genus", size = 3, alpha = 0.6) +
  scale_colour_brewer(palette = "Dark2")

# customise function options
ps \%>\%
  ps_calc_dominant(
    rank = "Family", other = "Other", none = "Not dominated",
    threshold = 0.4, n_max = 3
  ) \%>\%
  tax_transform(rank = "Genus", trans = "clr") \%>\%
  ord_calc("PCA") \%>\%
  ord_plot(colour = "dominant_Family", size = 3, alpha = 0.6) +
  scale_colour_manual(values = c(
    Bacteroidaceae = "forestgreen", Lachnospiraceae = "darkblue",
    Ruminococcaceae = "darkorange", Other = "red", "Not dominated" = "grey"
  ))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{heat_grid}
\alias{heat_grid}
\title{set options for drawing gridlines on heatmaps}
\usage{
heat_grid(
  col = "white",
  alpha = 1,
  lty = 1,
  lwd = 0.5,
  lex = 1,
  lineend = "round",
  linejoin = "round"
)
}
\arguments{
\item{col}{Colour for lines and borders.}

\item{alpha}{Alpha channel for transparency}

\item{lty}{Line type}

\item{lwd}{Line width}

\item{lex}{Multiplier applied to line width}

\item{lineend}{Line end style (round, butt, square)}

\item{linejoin}{Line join style (round, mitre, bevel)}
}
\description{
set options for drawing gridlines on heatmaps
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-sampleAnnotation.R
\name{anno_sample_cat}
\alias{anno_sample_cat}
\title{Helper to specify comp_heatmap annotation for categorical sample data}
\usage{
anno_sample_cat(
  var,
  col = distinct_palette(),
  renamer = identity,
  size = grid::unit(5, "mm"),
  legend = TRUE,
  legend_title = "",
  box_col = "white",
  box_lwd = 0.5,
  border_col = NA,
  border_lwd = 1,
  data = NULL,
  samples = NULL,
  which = NULL,
  ...
)
}
\arguments{
\item{var}{name of variable to use for annotation data}

\item{col}{colors vector, at least as long as unique(x), optionally named by x levels}

\item{renamer}{function to rename levels of variable \code{var}}

\item{size}{width or height as a grid unit object}

\item{legend}{generate legend for this annotation
(attached as attribute of heatmap, and not automatically included in plot)}

\item{legend_title}{title for legend, if drawn}

\item{box_col}{colour of boxes around individual cells}

\item{box_lwd}{line width of boxes around individual cells}

\item{border_col}{colour of border around all cells}

\item{border_lwd}{line width of border around all cells}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{samples}{OPTIONAL selection vector of sample names,
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}

\item{...}{
  Arguments passed on to \code{\link[=anno_cat]{anno_cat}}
  \describe{
    \item{\code{x}}{data vector, treated as categorical}
    \item{\code{width}}{grid unit object or NULL}
    \item{\code{height}}{grid unit object or NULL}
  }}
}
\value{
vector of values
}
\description{
Use this as an argument to sampleAnnotation(),
which itself is used by comp_heatmap() as sample_anno argument.
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- ibd_phylo
samples <- phyloseq::sample_names(psq)

# makes a function that takes data, taxa and which (at minimum)
fun <- anno_sample_cat(var = "ibd")

# manually specify the prevalence barplot function by giving it data etc.
heatmapAnnoFunction <- fun(data = psq, which = "row", samples = samples)

# draw the barplot without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)

grid::grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)
# A legend is attached by default to anno_cat() output, let's plot that.
pushViewport(viewport(x = 0.75))
draw(attr(heatmapAnnoFunction, "Legend"))

# change some options and specify the data up front
grid::grid.newpage()
pushViewport(vp)
anno_sample_cat(
  data = psq, var = "DiseaseState", samples = samples, which = "column",
  size = grid::unit(5, "cm"), col = distinct_palette(pal = "kelly")
) \%>\%
  draw()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_otu2samdat.R
\name{ps_otu2samdat}
\alias{ps_otu2samdat}
\title{Copy phyloseq otu_table data to sample_data}
\usage{
ps_otu2samdat(ps, taxa = NULL)
}
\arguments{
\item{ps}{phyloseq with sample_data}

\item{taxa}{list of taxa_names to copy to sample_data, or NULL (which selects all with \code{phyloseq::taxa_names()})}
}
\value{
phyloseq with augmented sample_data
}
\description{
Copy phyloseq otu_table data to sample_data
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")

ps <- dietswap \%>\% ps_otu2samdat("Akkermansia")
sample_variables(ps)

# or if you do not specify any taxa, all are copied
ps_all <- dietswap \%>\% ps_otu2samdat()
sample_variables(ps_all)[1:15]

# this could be useful for colouring ordination plots, for example
ps \%>\%
  ps_mutate(log_akkermansia = log(Akkermansia)) \%>\%
  dist_calc("bray") \%>\%
  ord_calc(method = "PCoA") \%>\%
  ord_plot(
    colour = "log_akkermansia",
    size = 3, shape = "nationality"
  )
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_fix.R
\name{tax_fix}
\alias{tax_fix}
\title{Replace unknown, NA, or short tax_table values}
\usage{
tax_fix(
  ps,
  min_length = 4,
  unknowns = NA,
  suffix_rank = "classified",
  sep = " ",
  anon_unique = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{ps}{phyloseq or tax_table (taxonomyTable)}

\item{min_length}{replace strings shorter than this}

\item{unknowns}{also replace strings matching any in this vector, NA default vector shown in details!}

\item{suffix_rank}{"classified" (default) or "current", when replacing an entry, should the suffix be taken from the lowest classified rank
for that taxon, "classified", or the "current" unclassified rank?}

\item{sep}{character(s) separating new name and taxonomic rank level suffix (see suffix_rank)}

\item{anon_unique}{make anonymous taxa unique by replacing unknowns with taxa_name?
otherwise they are replaced with paste("unknown", first_rank_name),
which is therefore the same for every anonymous taxon, meaning they will be merged if tax_agg is used.
(anonymous taxa are taxa with all unknown values in their tax_table row, i.e. cannot be classified even at highest rank available)}

\item{verbose}{emit warnings when cannot replace with informative name?}
}
\value{
object same class as ps
}
\description{
Identifies phyloseq tax_table values as unknown or uninformative and
replaces them with the first informative value from a higher taxonomic rank.
\itemize{
\item Short values in phyloseq tax_table are typically empty strings or " ", or "g__" etc.
so it is helpful to replace them. (If this is unwanted: set \code{min_length} = 0 to avoid filtering on length.)
\item Values in \code{unknowns} are also removed, even if longer than \code{min_length}.
It is up to the user to specify sensible values in \code{unknowns} if their dataset has other unwanted values.
\item NA values are also replaced.
}

See this article for an extended discussion of tax_table fixing.
\url{https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html}
}
\details{
By default (unknowns = NA), unknowns is set to a vector containing:

's__' 'g__' 'f__' 'o__' 'c__' 'p__' 'k__' 'S__' 'G__' 'F__' 'O__' 'C__' 'P__' 'K__' 'NA' 'NaN' ' ' ''
'unknown' 'Unknown' 's__unknown' 's__Unknown' 's__NA' 'g__unknown' 'g__Unknown' 'g__NA'
'f__unknown' 'f__Unknown' 'f__NA' 'o__unknown' 'o__Unknown' 'o__NA' 'c__unknown' 'c__Unknown' 'c__NA'
'p__unknown' 'p__Unknown' 'p__NA' 'k__unknown' 'k__Unknown' 'k__NA' 'S__unknown' 'S__Unknown' 'S__NA'
'G__unknown' 'G__Unknown' 'G__NA' 'F__unknown' 'F__Unknown' 'F__NA' 'O__unknown' 'O__Unknown' 'O__NA'
'C__unknown' 'C__Unknown' 'C__NA' 'P__unknown' 'P__Unknown' 'P__NA' 'K__unknown' 'K__Unknown' 'K__NA'
}
\examples{
library(dplyr)
library(phyloseq)

data(dietswap, package = "microbiome")
ps <- dietswap

# create unknowns to test filling
tt <- tax_table(ps)
ntax <- ntaxa(ps)
set.seed(123)
g <- sample(1:ntax, 30)
f <- sample(g, 10)
p <- sample(f, 3)
tt[g, 3] <- "g__"
tt[f, 2] <- "f__"
tt[p, 1] <- "p__"
tt[sample(1:ntax, 10), 3] <- "unknown"
# create a row with only NAs
tt[1, ] <- NA

tax_table(ps) <- tax_table(tt)

ps
# tax_fix with defaults should solve most problems
tax_table(ps) \%>\% head(50)

# this will replace "unknown"s as well as short values including "g__" and "f__"
tax_fix(ps) \%>\%
  tax_table() \%>\%
  head(50)

# This will only replace short entries, and so won't replace literal "unknown" values
ps \%>\%
  tax_fix(unknowns = NULL) \%>\%
  tax_table() \%>\%
  head(50)

# Change rank suffix and separator settings
tax_fix(ps, suffix_rank = "current", sep = " - ") \%>\%
  tax_table() \%>\%
  head(50)

# by default, completely unclassified (anonymous) taxa are named by their
# taxa_names / rownames at all ranks.
# This makes anonymous taxa distinct from each other,
# and so they won't be merged on aggregation with tax_agg.
# If you think your anonymous taxa should merge on tax_agg,
# or you just want them to be named the all same for another reason,
# set anon_unique = FALSE (compare the warning messages)
tax_fix(ps, anon_unique = FALSE)
tax_fix(ps, anon_unique = TRUE)

# here's a larger example tax_table shows its still fast with 1000s rows,
# from microbiomeutilities package
# library(microbiomeutilities)
# data("hmp2")
# system.time(tax_fix(hmp2, min_length = 1))
}
\seealso{
\code{\link{tax_fix_interactive}} for interactive tax_fix help
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_transform.R
\name{tax_transform}
\alias{tax_transform}
\title{Transform taxa in phyloseq object and record transformation}
\usage{
tax_transform(
  data,
  trans,
  rank = NA,
  keep_counts = TRUE,
  chain = FALSE,
  zero_replace = 0,
  add = 0,
  transformation = NULL,
  ...
)
}
\arguments{
\item{data}{a phyloseq object or \code{ps_extra} list output from \code{tax_agg}}

\item{trans}{any valid taxa transformation (e.g. from \code{microbiome::transform})}

\item{rank}{If data is phyloseq: data are aggregated at this rank before transforming.
If NA, runs tax_agg(data, rank = NA).
If rank is NA and data is already ps_extra, any preceding aggregation is left as is.}

\item{keep_counts}{if TRUE, store the pre-transformation count data in ps_extra counts slot}

\item{chain}{if TRUE, transforming again is possible when data are already transformed
i.e. multiple transformations can be chained with multiple tax_transform calls}

\item{zero_replace}{Replace any zeros with this value before transforming. Either a numeric, or
"halfmin" which replaces zeros with half of the smallest value across the
entire dataset.
Beware: the choice of zero replacement is not tracked in the ps_extra output.}

\item{add}{Add this value to the otu_table before transforming. If \code{add} != 0,
\code{zero_replace} does nothing. Either a numeric, or "halfmin".
Beware: this choice is not tracked in the ps_extra output.}

\item{transformation}{deprecated, use \code{trans} instead!}

\item{...}{any extra arguments passed to \code{microbiome::transform} or pass
undetected = \verb{a number} when using trans = "binary"}
}
\value{
\code{ps_extra} list including phyloseq object and info
}
\description{
Transform taxa features, and optionally aggregate at specified taxonomic rank beforehand.
You can pipe the results of \code{tax_agg} into \code{tax_transform},
or equivalently set the rank argument in \code{tax_transform}.
}
\details{
This function often uses \code{microbiome::transform} internally and can perform the
same transformations, including many from \code{vegan::decostand} (where the default MARGIN = 2).
See below for notes about some of the available transformations.

\code{tax_transform} returns a \code{ps_extra} list containing the transformed phyloseq object and
extra info (used for annotating \code{ord_plot} ordinations):
\itemize{
\item tax_transform (a string recording the transformation),
\item tax_agg (a string recording the taxonomic aggregation rank if specified here or earlier in \code{tax_agg}).
}

A few commonly used transformations:
\itemize{
\item "clr" performs the centered log ratio transformation using \code{microbiome::transform}
\item "compositional" converts the data into proportions, from 0 to 1.
\item "identity" does not transform the data, and records this choice for \code{ord_plot}
\item "binary" can be used to transform tax abundances into presence/abundance data.
\item "log2" which performs a log base 2 transformation
(don't forget to set zero_replace if there are any zeros in your data)
}
}
\section{clr transformation note}{


If any values are zero, the clr transform routine first adds a small
pseudocount of min(relative abundance)/2 to all values. To avoid this, you
can replace any zeros in advance by setting zero_replace to a number > 0.
}

\section{Binary transformation notes}{


By default, otu_table values of 0 are kept as 0, and all positive values
are converted to 1 (like \code{decostand(method = "pa")}).
You can set a different threshold, by passing e.g. undetected = 10, for
example, in which case all abundances of 10 or below would be converted to 0.
All abundances above 10 would be converted to 1s.

Beware that the choice of detection threshold is not tracked in the ps_extra.
}

\examples{
data("dietswap", package = "microbiome")

# aggregate taxa at Phylum level and center log ratio transform the phyla counts
tax_transform(dietswap, trans = "clr", rank = "Phylum")

# this is equivalent to the two-step method (agg then transform)
tax_agg(dietswap, rank = "Phylum") \%>\% tax_transform("clr")

# does nothing except record tax_agg as "unique" and tax_transform as "identity" in ps_extra info
dietswap \%>\% tax_transform("identity", rank = NA)

# binary transformation (convert abundances to presence/absence or detected/undetected)
tax_transform(dietswap, trans = "binary", rank = "Genus")
# change detection threshold by setting undetected argument (default is 0)
tax_transform(dietswap, trans = "binary", rank = "Genus", undetected = 50) \%>\%
  otu_get() \%>\%
  .[1:6, 1:4]

# log2 transformation after replacing all zeros with a pseudocount of 1
tax_transform(dietswap, trans = "log2", rank = "Family", zero_replace = 1)

# log2 transformation after replacing all zeros with a pseudocount of half
# the minimum non-zero count value in the aggregated dataset
tax_transform(dietswap, trans = "log2", rank = "Family", zero_replace = "halfmin")
}
\seealso{
\code{microbiome::\link[microbiome]{transform}} for some more info on available transformations

\code{vegan::\link[vegan]{decostand}} for even more transformation options

\code{\link{tax_agg}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_scale.R
\name{tax_scale}
\alias{tax_scale}
\title{Mean-center and SD-scale taxa in phyloseq}
\usage{
tax_scale(data, center = TRUE, scale = TRUE, do = NA, keep_counts = TRUE)
}
\arguments{
\item{data}{phyloseq or ps_extra or otu_table}

\item{center}{if TRUE: center each taxon by subtracting its mean}

\item{scale}{if TRUE, divide each centred taxon by its standard deviation (or divide by RMS if not centred!)}

\item{do}{alternative argument that overrides center and scale options! takes "both", "scale", "center" or "neither"}

\item{keep_counts}{if TRUE, retain the original count data in ps_extra counts slot}
}
\description{
Wrapper for applying base scale function to phyloseq otu_table
}
\examples{
data("dietswap", package = "microbiome")
ps <- dietswap
ps \%>\%
  otu_get() \%>\%
  .[1:6, 1:6]

# standard use (mean center and SD scale)
tax_scale(ps) \%>\%
  otu_get() \%>\%
  .[1:6, 1:6] # Aerococcus is NaN as standard deviation = 0 (0 prevalence)

# RMS scale only (directly on otu_table)
otu_get(ps) \%>\%
  tax_scale(center = FALSE) \%>\%
  .[1:6, 1:6] # Aerococcus is NaN as standard deviation = 0 (0 prevalence)

# example using alternative `do` argument (to center only, no scaling)
tax_scale(ps, do = "center") \%>\%
  otu_get() \%>\%
  .[1:6, 1:6]

# preserves existing info
tax_transform(ps, "compositional", rank = "Genus") \%>\% tax_scale()

# drop other ps_extra objects previously calculated with unscaled data
psxDist <- tax_agg(ps, "Genus") \%>\% dist_calc()
psxDist
psxDist \%>\% tax_scale()
tax_scale(psxDist) \%>\% info_get()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_sort_ord.R
\name{tax_sort_ord}
\alias{tax_sort_ord}
\alias{ord_order_taxa}
\title{Order taxa in phyloseq by their loading vectors}
\usage{
tax_sort_ord(ps, ord, axes = 1:2, scaling = 2)

ord_order_taxa(ord, axes = 1:2, scaling = 2)
}
\arguments{
\item{ps}{phyloseq object to be sorted}

\item{ord}{ps_extra with ordination object}

\item{axes}{which axes to use for sorting? numerical vector of length 1 or 2}

\item{scaling}{Type 2, or type 1 scaling. For more info, see \url{https://sites.google.com/site/mb3gustame/constrained-analyses/rda}.
Either "species" or "site" scores are scaled by (proportional) eigenvalues, and the other set of scores is left unscaled (from ?vegan::scores.cca)}
}
\description{
\code{tax_sort_ord} reorders taxa in a phyloseq object based on the relative
length of their taxa scores / "loading" vector lengths on 1 or 2 ordination axes.

\code{ord_order_taxa} gets the taxa names in order from the ordination
contained in a ps_extra list. This is used internally by \code{tax_sort_ord}.
}
\seealso{
\itemize{
\item These functions were created to support ordering of taxa bars on \code{ord_plot_iris}
\item \code{ps_sort_ord} for ordering samples in phyloseq by ordination
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp_barplot.R
\name{comp_barplot}
\alias{comp_barplot}
\title{Plot (grouped and ordered) compositional barplots}
\usage{
comp_barplot(
  ps,
  tax_level,
  n_taxa = 8,
  tax_order = sum,
  merge_other = TRUE,
  taxon_renamer = function(x) identity(x),
  sample_order = "bray",
  order_with_all_taxa = FALSE,
  label = "SAMPLE",
  group_by = NA,
  facet_by = NA,
  bar_width = 1,
  bar_outline_colour = "grey5",
  bar_outline_width = 0.1,
  palette = distinct_palette(n_taxa),
  tax_transform_for_ordering = "identity",
  tax_transform_for_plot = "compositional",
  seriate_method = "OLO_ward",
  keep_all_vars = TRUE,
  interactive = FALSE,
  max_taxa = 10000,
  other_name = "other",
  ...
)
}
\arguments{
\item{ps}{phyloseq object}

\item{tax_level}{taxonomic aggregation level (from rank_names(ps))}

\item{n_taxa}{how many taxa to show distinct colours for (all others grouped into "Other")}

\item{tax_order}{order of taxa within the bars, either a function for tax_sort (e.g. sum),
or a vector of (all) taxa names at tax_level to set order manually}

\item{merge_other}{if FALSE, taxa coloured/filled as "other" remain distinct,
and so can have bar outlines drawn around them}

\item{taxon_renamer}{function that takes taxon names and returns modified names for legend}

\item{sample_order}{vector of sample names;
or any distance measure in dist_calc that doesn't require phylogenetic tree;
or "default" for the order returned by phyloseq::sample_names(ps)}

\item{order_with_all_taxa}{if TRUE, this will always use all taxa (not just the top n_taxa)
to calculate distances for sample ordering}

\item{label}{sample label variable name}

\item{group_by}{splits dataset by this variable (must be categorical)
\itemize{
\item resulting in a list of plots, one for each level of the group_by variable.
}}

\item{facet_by}{facets plots by this variable (must be categorical). If group_by is also set
the faceting will occur separately in the plot for each group.}

\item{bar_width}{default 1 avoids random gapping otherwise seen with many samples
(set to something less than 1 to introduce gaps between fewer samples)}

\item{bar_outline_colour}{line colour separating taxa and samples
(use NA for none)}

\item{bar_outline_width}{width of line separating taxa and samples
(for no outlines set bar_outline_colour = NA)}

\item{palette}{palette for taxa fill colours}

\item{tax_transform_for_ordering}{transformation of taxa values used before ordering samples by similarity}

\item{tax_transform_for_plot}{default "compositional" draws proportions of total counts per sample,
but you could reasonably use another transformation,
e.g. "identity", if you have truly quantitative microbiome profiling data}

\item{seriate_method}{name of any ordering method suitable for distance matrices
(see ?seriation::seriate)}

\item{keep_all_vars}{FALSE may speed up internal melting with ps_melt for large phyloseq objects
but TRUE is required for some post-hoc plot customisation}

\item{interactive}{creates plot suitable for use with ggiraph}

\item{max_taxa}{maximum distinct taxa groups to show
(only really useful for limiting complexity of interactive plots
e.g. within ord_explore)}

\item{other_name}{name for other taxa after N}

\item{...}{extra arguments passed to facet_wrap() (if facet_by is not NA)}
}
\value{
ggplot or list of harmonised ggplots
}
\description{
Stacked barplots showing composition of phyloseq samples for a specified number of coloured taxa.
Normally your phyloseq object should contain counts data,
as by default \code{comp_barplot()} performs the "compositional" taxa transformation for you,
and requires count input for some sample_order methods!
}
\details{
\itemize{
\item sample_order: Either specify a list of sample names to order manually, or the bars/samples can/will be sorted by similarity, according to a specified distance measure (default 'bray'-curtis),
\item seriate_method specifies a seriation/ordering algorithm (default Ward hierarchical clustering with optimal leaf ordering, see seriation::list_seriation_methods())
\item group_by: You can group the samples on distinct plots by levels of a variable in the phyloseq object. The list of ggplots produced can be arranged flexibly with the patchwork package functions. If you want to group by several variables you can create an interaction variable with interaction(var1, var2) in the phyloseq sample_data BEFORE using comp_barplot.
\item facet_by can allow faceting of your plot(s) by a grouping variable. Using this approach is less flexible than using group_by but means you don't have to arrange a list of plots yourself like with the group_by argument. Using facet_by is equivalent to adding a call to facet_wrap(facets = facet_by, scales = "free") to your plot(s). Calling facet_wrap() yourself is itself a more flexible option as you can add other arguments like the number of rows etc. However you must use keep_all_vars = TRUE if you will add faceting manually.
\item bar_width: No gaps between bars, unless you want them (decrease width argument to add gaps between bars).
\item bar_outline_colour: Bar outlines default to "grey5" for almost black outlines. Use NA if you don't want outlines.
\item merge_other: controls whether bar outlines can be drawn around individual (lower abundance) taxa that are grouped in "other" category. If you want to see the diversity of taxa in "other" use merge_taxa = FALSE, or use TRUE if you prefer the cleaner merged look
\item palette: Default colouring is consistent across multiple plots if created with the group_by argument, and the defaults scheme retains the colouring of the most abundant taxa irrespective of n_taxa
}


}
\examples{
library(ggplot2)
data(dietswap, package = "microbiome")

# illustrative simple customised example
dietswap \%>\%
  ps_filter(timepoint == 1) \%>\%
  comp_barplot(
    tax_level = "Family", n_taxa = 8,
    bar_outline_colour = NA,
    sample_order = "bray",
    bar_width = 0.7,
    taxon_renamer = toupper
  ) + coord_flip()

# change colour palette with the distinct_palette() function
# remember to set the number of colours to the same as n_taxa argument!
dietswap \%>\%
  ps_filter(timepoint == 1) \%>\%
  comp_barplot(
    tax_level = "Family", n_taxa = 8,
    bar_outline_colour = NA,
    sample_order = "bray",
    bar_width = 0.7,
    palette = distinct_palette(8, pal = "kelly"),
    taxon_renamer = toupper
  ) + coord_flip()

# Order samples by the value of one of more sample_data variables.
# Use ps_arrange and set sample_order = "default" in comp_barplot.
# ps_mutate is also used here to create an informative variable for axis labelling
dietswap \%>\%
  ps_mutate(subject_timepoint = interaction(subject, timepoint)) \%>\%
  ps_filter(nationality == "AAM", group == "DI", sex == "female") \%>\%
  ps_arrange(desc(subject), desc(timepoint)) \%>\%
  comp_barplot(
    tax_level = "Genus", n_taxa = 12,
    sample_order = "default",
    bar_width = 0.7,
    bar_outline_colour = "black",
    order_with_all_taxa = TRUE,
    label = "subject_timepoint"
  ) + coord_flip()

# how many taxa are in those light grey "other" bars?
# set merge_other, to find out (& remember to set a bar_outline_colour)
dietswap \%>\%
  ps_mutate(subject_timepoint = interaction(subject, timepoint)) \%>\%
  ps_filter(nationality == "AAM", group == "DI", sex == "female") \%>\%
  ps_arrange(desc(subject), desc(timepoint)) \%>\%
  comp_barplot(
    tax_level = "Genus", n_taxa = 12,
    sample_order = "default",
    merge_other = FALSE,
    bar_width = 0.7,
    bar_outline_colour = "black",
    order_with_all_taxa = TRUE,
    label = "subject_timepoint"
  ) + coord_flip()


# Often to compare groups, average compositions are presented
p1 <- phyloseq::merge_samples(dietswap, group = "group") \%>\%
  comp_barplot(
    tax_level = "Genus", n_taxa = 12,
    sample_order = c("ED", "HE", "DI"),
    bar_width = 0.8
  ) +
  coord_flip() + labs(x = NULL, y = NULL)
p1

# However that "group-averaging" approach hides a lot of within-group variation
p2 <- comp_barplot(dietswap,
  tax_level = "Genus", n_taxa = 12, group_by = "group",
  sample_order = "euclidean", bar_outline_colour = NA
) \%>\%
  patchwork::wrap_plots(nrow = 3, guides = "collect") &
  coord_flip() & labs(x = NULL, y = NULL) &
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p2

# Only from p2 you can see that the apparently higher average relative abundance
# of Oscillospira in group DI is probably driven largely by a subgroup
# of DI samples with relatively high Oscillospira.

# make a list of 2 harmonised composition plots (grouped by sex)
p <- comp_barplot(dietswap,
  n_taxa = 15, tax_level = "Genus",
  bar_outline_colour = "black", merge_other = TRUE,
  sample_order = "aitchison", group_by = "sex"
)

# plot them side by side with patchwork package
patch <- patchwork::wrap_plots(p, ncol = 2, guides = "collect")
patch & coord_flip() # make bars in all plots horizontal (note: use & instead of +)

# beautifying tweak #
# modify one plot in place (flip the order of the samples in the 2nd plot)
# notice that the scaling is for the x-axis
# (that's because coord_flip is used afterwards when displaying the plots
patch[[2]] <- patch[[2]] + scale_x_discrete(limits = rev)
# Explainer: rev() function takes current limits and reverses them.
# You could also pass a completely arbitrary order, naming all samples

# you can theme all plots with the & operator
patch & coord_flip() &
  theme(axis.text.y = element_text(size = 5), legend.text = element_text(size = 6))
# See https://patchwork.data-imaginist.com/index.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_stats_p_adjust.R
\name{taxatree_stats_p_adjust}
\alias{taxatree_stats_p_adjust}
\title{Adjust p values in taxatree_stats dataframe}
\usage{
taxatree_stats_p_adjust(
  data,
  method,
  grouping = "rank",
  p = "p.value",
  new_var = paste0("p.adj.", method, ".", grouping)
)
}
\arguments{
\item{data}{ps_extra with taxatree_stats dataframe, or just the dataframe}

\item{method}{any method from \code{stats::p.adjust.methods}}

\item{grouping}{defines grouping of p-values into families for adjustment, see details.}

\item{p}{name of variable containing p values for adjustment}

\item{new_var}{name of new variable created for adjusted p values
(automatically inferred by default)}
}
\value{
ps_extra with dataframe of statistics, or just the data.frame
}
\description{
Apply a p value adjustment method from \code{stats::p.adjust.methods},
such as false-discovery rate adjustment with "BH",
or more conservative family-wise error rate controlling methods such as "holm" or "bonferroni".
}
\details{
Define how to group the p values for adjustment with the \code{grouping} argument.
The default is to adjust the p values in groups at each taxonomic rank,
but you could also adjust per "model" / "taxon" or per "term".
Or even group by a combination of rank and term with c("rank", "term")
}
\examples{
# This example is an abbreviated excerpt from article on taxon modelling on
# the microViz documentation website

library(corncob)
library(dplyr)
data("ibd_phylo", package = "corncob")

# We'll keep only the Ulcerative Colitis and Healthy Control samples, to
# simplify the analyses for this example. We'll also remove the Species
# rank information, as most OTUs in this dataset are not assigned to a
# species. We'll also use `tax_fix` to fill any gaps where the Genus is
# unknown, with the family name or whatever higher rank classification is
# known.

phylo <- ibd_phylo \%>\%
  ps_filter(DiseaseState \%in\% c("UC", "nonIBD")) \%>\%
  tax_mutate(Species = NULL) \%>\%
  tax_fix()

# Let's make some sample data variables that are easier to use and compare
# in the statistical modelling ahead. We will convert dichotomous
# categorical variables into similar binary variables (values: 1 for true,
# or 0 for false). We will also scale and center the numeric variable for
# age.

phylo <- phylo \%>\%
  ps_mutate(
    UC = ifelse(DiseaseState == "UC", yes = 1, no = 0),
    female = ifelse(gender == "female", yes = 1, no = 0),
    antibiotics = ifelse(abx == "abx", yes = 1, no = 0),
    steroids = ifelse(steroids == "steroids", yes = 1, no = 0),
    age_scaled = scale(age, center = TRUE, scale = TRUE)
  )

bb_models <- phylo \%>\%
  tax_fix() \%>\%
  tax_prepend_ranks() \%>\%
  tax_filter(min_prevalence = 0.3) \%>\%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order"),
    variables = c("UC", "female", "antibiotics", "steroids", "age_scaled")
  )

bb_stats <- bb_models \%>\%
  taxatree_models2stats(param = "mu") \%>\%
  taxatree_stats_p_adjust(method = "BH", grouping = "rank")

bb_stats

bb_stats$taxatree_stats

# you can also directly modify the dataframe,
# and choose a different variable name
bb_stats$taxatree_stats \%>\%
  taxatree_stats_p_adjust(
    method = "holm", grouping = "taxon", new_var = "p_adj_holm"
  )

# see all available adjustment methods
stats::p.adjust.methods
}
\seealso{
\code{\link{taxatree_models2stats}}

\code{\link{taxatree_models}}

\code{stats::\link[stats]{p.adjust}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap-utils.R
\name{prev}
\alias{prev}
\title{Calculate prevalence from numeric vector}
\usage{
prev(x, undetected = 0)
}
\arguments{
\item{x}{numeric vector (of taxon counts or proportions)}

\item{undetected}{value above which a taxon is considered present or detected}
}
\value{
numeric value
}
\description{
Useful as helper for taxon prevalence calculation
}
\examples{
prev(c(0, 0, 1, 2, 4))
prev(c(0, 0, 1, 2, 4), undetected = 1.5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_plotkey.R
\name{taxatree_plot_labels}
\alias{taxatree_plot_labels}
\title{Add labels to taxatree plots/key}
\usage{
taxatree_plot_labels(
  p,
  circular = TRUE,
  taxon_renamer = identity,
  fun = ggrepel::geom_text_repel,
  label_var = "label",
  x_nudge = 0.1,
  y_nudge = 0.025,
  fontface = "bold",
  size = 2.5,
  colour = "grey15",
  max.overlaps = Inf,
  min.segment.length = 0,
  segment.size = 0.15,
  segment.color = "grey15",
  point.padding = 0.05,
  box.padding = 0.1,
  seed = NA,
  ...
)
}
\arguments{
\item{p}{taxatree_plotkey or taxatree_plots output plot}

\item{circular}{is the plot layout circular? labels are drawn differently for circular trees}

\item{taxon_renamer}{function that takes taxon names and returns modified names for labels}

\item{fun}{ggrepel labelling function: geom_text_repel or geom_label_repel}

\item{label_var}{name of variable in taxatree_stats that indicates which taxa to label}

\item{x_nudge}{absolute amount by which the initial position of taxon labels is nudged
(relevant only for circular layouts, use nudge_x for other layouts)}

\item{y_nudge}{absolute amount by which the initial position of taxon labels is nudged
(relevant only for circular layouts, use nudge_y for other layouts)}

\item{fontface}{fontface of label text}

\item{size}{size of labels}

\item{colour}{colour of label outlines and text}

\item{max.overlaps}{max number of overlapping labels tolerated}

\item{min.segment.length}{min length of label line segment to bother drawing}

\item{segment.size}{thickness of line segment}

\item{segment.color}{colour of line segment}

\item{point.padding}{padding around node points (for label positioning)}

\item{box.padding}{padding around labels/text (for label positioning)}

\item{seed}{set this for reproducible label positions}

\item{...}{
  Arguments passed on to \code{\link[ggrepel:geom_text_repel]{ggrepel::geom_text_repel}}
  \describe{
    \item{\code{arrow}}{specification for arrow heads, as created by \code{\link[grid]{arrow}}}
    \item{\code{force}}{Force of repulsion between overlapping text labels. Defaults
to 1.}
    \item{\code{force_pull}}{Force of attraction between a text label and its
corresponding data point. Defaults to 1.}
    \item{\code{max.time}}{Maximum number of seconds to try to resolve overlaps.
Defaults to 0.5.}
    \item{\code{max.iter}}{Maximum number of iterations to try to resolve overlaps.
Defaults to 10000.}
    \item{\code{xlim}}{Limits for the x and y axes. Text labels will be constrained
to these limits. By default, text labels are constrained to the entire plot
area.}
    \item{\code{ylim}}{Limits for the x and y axes. Text labels will be constrained
to these limits. By default, text labels are constrained to the entire plot
area.}
    \item{\code{direction}}{"both", "x", or "y" -- direction in which to adjust position of labels}
    \item{\code{verbose}}{If \code{TRUE}, some diagnostics of the repel algorithm are printed}
  }}
}
\description{
Finer control over label drawing for \code{taxatree_plotkey}
(with .draw_label = FALSE),
and label drawing for \code{taxatree_plots} output too.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_calc.R
\name{dist_calc}
\alias{dist_calc}
\title{Calculate distances between pairs of samples in phyloseq object}
\usage{
dist_calc(
  data,
  dist = c("bray", "gunifrac", "unifrac", "wunifrac", "va-wunifrac", "aitchison",
    "euclidean")[1],
  gunifrac_alpha = 0.5,
  ...
)
}
\arguments{
\item{data}{ps_extra object, e.g. output from tax_transform()}

\item{dist}{name of distance to calculate between pairs of samples}

\item{gunifrac_alpha}{setting alpha value only relevant if gunifrac distance used}

\item{...}{optional distance-specific named arguments passed to phyloseq::distance()}
}
\value{
list with distance matrix, phyloseq object, and name of distance used
}
\description{
Can compute various sample-sample distances using the microbiota composition of your samples:
\itemize{
\item Bray Curtis ('bray') or any other ecological distance from phyloseq::distance() / vegan::vegdist()
\item UniFrac distances (using the GUniFrac package)
\itemize{
\item generalised: 'gunifrac' (optionally set weighting alpha in gunifrac alpha)
\item unweighted: 'unifrac'
\item weighted: 'wunifrac'
\item variance adjusted weighted: 'va-wunifrac'
}
\item Aitchison distance (Euclidean distance after centered log ratio transform clr, see details)
\item Euclidean distance
}

Use dist_calc with ps_extra output of tax_transform (or tax_agg).
It returns a ps_extra object containing the phyloseq and the name of the distance used
in addition to the distance matrix itself.
The resulting object is intended to be piped into ord_calc or dist_permanova functions.
Alternatively you can directly access the distance matrix with dist_get().
}
\section{Aitchison distance note}{


You should EITHER:
\enumerate{
\item skip the dist_calc function and call ord_calc(method = "PCA") directly on an object with taxa transformed with tax_transform(trans = "clr")
\item pass an object with untransformed (or 'identity' transformed) taxa to the data argument of dist_calc() and specify dist = "aitchison".
}

If ordination plots with taxon loading vectors are desired, users require option 1.
If the distance matrix is required for permanova, users require option 2.
}

\section{Binary Jaccard distance note}{


Jaccard distance can be computed on abundances, but often in microbiome
research it is the Binary Jaccard distance that is desired. So remember to
first perform a "binary" transformation with \code{tax_transform("binary")},
OR pass an additional argument to \code{dist_calc("jaccard", binary = TRUE)}
}

\examples{
# bray curtis distance on genera-level features
data("dietswap", package = "microbiome")
bc <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("bray")
bc
class(bc)

# gunifrac distance using phyloseq input
data("esophagus", package = "phyloseq")
gunifrac <- esophagus \%>\%
  dist_calc("gunifrac") \%>\%
  dist_get()
class(gunifrac)
}
\seealso{
\code{\link{tax_transform}} for the function to use before dist_calc

\code{\link{ord_calc}}

\code{\link{ord_plot}}

\code{\link{dist_permanova}}

\code{phyloseq::\link[phyloseq:distance]{distance}}

\code{vegan::\link[vegan:vegdist]{vegdist}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_plotkey.R
\name{taxatree_plotkey}
\alias{taxatree_plotkey}
\title{Draw labelled key to accompany taxatree_plots}
\usage{
taxatree_plotkey(
  data,
  ...,
  size_stat = list(prevalence = prev),
  node_size_range = c(1.5, 5),
  edge_width_range = node_size_range * 0.8,
  size_guide = "none",
  size_trans = "identity",
  colour = "lightgrey",
  edge_alpha = 0.7,
  title = "Key",
  title_size = 14,
  taxon_renamer = identity,
  .combine_label = any,
  .draw_label = TRUE,
  .calc_label = TRUE,
  layout = "tree",
  layout_seed = NA,
  circular = identical(layout, "tree"),
  node_sort = NULL,
  add_circles = isTRUE(circular),
  drop_ranks = TRUE
)
}
\arguments{
\item{data}{ps_extra (or phyloseq)}

\item{...}{logical conditions for labelling
e.g. rank == "Phylum", p.value < 0.1 | taxon \%in\% listOfTaxa}

\item{size_stat}{named list of length 1, giving function calculated for each taxon,
to determine the size of nodes (and edges). Name used as size legend title.}

\item{node_size_range}{min and max node sizes, decrease to avoid node overlap}

\item{edge_width_range}{min and max edge widths}

\item{size_guide}{guide for node sizes, try "none", "legend" or ggplot2::guide_legend()}

\item{size_trans}{transformation for size scale
you can use (the name of) any transformer from the scales package,
such as "identity", "log1p", or "sqrt"}

\item{colour}{fixed colour and fill of nodes and edges}

\item{edge_alpha}{fixed alpha value for edges}

\item{title}{title of plot (NULL for none)}

\item{title_size}{font size of title}

\item{taxon_renamer}{function that takes taxon names and returns modified names for labels}

\item{.combine_label}{all or any: function to combine multiple logical "label" values for a taxon
(relevant if taxatree_stats already present in data)}

\item{.draw_label}{should labels be drawn, or just the bare tree, set this to FALSE if you want
to draw customised labels with taxatree_plot_labels afterwards}

\item{.calc_label}{if you already set labels with taxatree_label:
set this to FALSE to use / avoid overwriting that data
(ignores \code{...} if FALSE)}

\item{layout}{any ggraph layout, default is "tree"}

\item{layout_seed}{any numeric, required if a stochastic igraph layout is named}

\item{circular}{should the layout be circular?}

\item{node_sort}{sort nodes by "increasing" or "decreasing" size? NULL for no sorting.
Use \code{tax_sort()} before \code{taxatree_plots()} for finer control.}

\item{add_circles}{add faint concentric circles to plot, behind each rank?}

\item{drop_ranks}{drop ranks from tree if not included in stats dataframe}
}
\value{
ggplot
}
\description{
Draw labelled key to accompany taxatree_plots
}
\examples{
# see taxatree_plots() examples
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap_annotations-deprecated.R
\name{deprecated-heatmap-annotations}
\alias{deprecated-heatmap-annotations}
\alias{tax_anno}
\alias{anno_prev}
\alias{anno_abund}
\alias{var_anno}
\alias{old_anno_var_hist}
\alias{old_anno_var_box}
\title{DEPRECATED Heatmap annotations helpers}
\usage{
tax_anno(
  undetected = 0,
  which = NA,
  prev = 1,
  abund = 2,
  size = 30,
  gap = 2,
  rel_sizes = NA,
  args = NULL,
  ...
)

anno_prev(
  data,
  taxa,
  undetected = 0,
  which = "row",
  size = 15,
  bar_width = 0.6,
  gp = grid::gpar(fill = "grey85"),
  ...
)

anno_abund(
  data,
  taxa,
  undetected = 0,
  which = "row",
  size = 15,
  point_size = 0.75,
  box_width = 0.6,
  gp = grid::gpar(fill = "grey85"),
  ...
)

var_anno(
  annos = "var_box",
  funs = "identity",
  names = NA,
  which = "column",
  size = 15 * length(annos),
  gap = 2,
  rel_sizes = NA,
  args = NULL,
  ...
)

old_anno_var_hist(data, vars = NA, which = "column", size = 15, ...)

old_anno_var_box(data, vars = NA, which = "column", size = 15, ...)
}
\arguments{
\item{undetected}{value above which taxa are considered present/detected in a sample}

\item{which}{"row" or "column" annnotation}

\item{prev}{order in which prevalence annotation shown (number, or NA to not show)}

\item{abund}{order in which abundance annotation shown (number, or NA to not show)}

\item{size}{total size (mm) of annotations (width/height depending on which)}

\item{gap}{gap in mm between annotations}

\item{rel_sizes}{relative sizes of annotations (NA for equal sizes, or same length as annos)}

\item{args}{extra args passed to each annotation: give as list of lists
(one inner list per arg, named, e.g. list(prev = list(whatever = whatever))}

\item{...}{further named args to be passed on (to list)}

\item{data}{phyloseq or ps-extra (or a data.frame or matrix for anno_var_* functions)}

\item{taxa}{names of taxa to plot}

\item{bar_width}{relative width of barchart bars}

\item{gp}{a grid::gpar() object for graphics parameter settings like fill or lwd}

\item{point_size}{size of outlier points in mm}

\item{box_width}{relative width of boxplot boxes}

\item{annos}{name(s) of annotation(s) to show, in order (e.g. 'var_box', 'var_hist')}

\item{funs}{function(s) to transform matrix of variable values before plotting
(length must be 1 or same length as annos)}

\item{names}{names to use for each annotation in annos}

\item{vars}{names of variables to plot}
}
\description{
Functions to easily define ComplexHeatmap annotations for taxa and/or variables
\itemize{
\item tax_anno creates list describing taxa annotation (for cor_heatmap or comp_heatmap)
\item var_anno creates list describing variable annotation (for cor_heatmap)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_arrows.R
\name{Ordination-arrows}
\alias{Ordination-arrows}
\alias{vec_constraint}
\alias{vec_tax_sel}
\alias{vec_tax_all}
\title{Create ordination plot vector styling lists

Used by ord_plot, see examples there.}
\usage{
vec_constraint(
  size = 1,
  alpha = 0.8,
  colour = "brown",
  arrow = grid::arrow(length = grid::unit(0.005, units = "npc"), type = "closed", angle
    = 30),
  lineend = "round",
  linejoin = "mitre",
  ...
)

vec_tax_sel(
  size = 0.5,
  alpha = 1,
  colour = "black",
  arrow = grid::arrow(length = grid::unit(0.005, units = "npc"), type = "closed", angle
    = 30),
  lineend = "round",
  linejoin = "mitre",
  ...
)

vec_tax_all(size = 0.5, alpha = 0.25, arrow = NULL, ...)
}
\arguments{
\item{size}{width of vector}

\item{alpha}{opacity of vector}

\item{colour}{colour of vector}

\item{arrow}{arrow style specified with grid::arrow() or NULL for no arrow}

\item{lineend}{Line end style (round, butt, square).}

\item{linejoin}{Line join style (round, mitre, bevel).}

\item{...}{further arguments passed to geom_segment}
}
\value{
list
}
\description{
Create ordination plot vector styling lists

Used by ord_plot, see examples there.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-sampleAnnotation.R
\name{anno_sample}
\alias{anno_sample}
\title{Helper to specify simple comp_heatmap annotation for other sample data}
\usage{
anno_sample(var, fun = identity, data = NULL, samples = NULL)
}
\arguments{
\item{var}{name of variable to use for annotation data}

\item{fun}{function to transform variable \code{var}}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{samples}{OPTIONAL selection vector of sample names,
only set this if providing data argument to override default}
}
\value{
vector of values
}
\description{
Use this as an argument to sampleAnnotation(),
which itself is used by comp_heatmap() as sample_anno argument.

This creates a vector, which sampleAnnotation() interprets as a
simple annotation, so then you set colours and legend parameters
for each simple annotation as further arguments in sampleAnnotation.
}
\examples{
# see `?sampleAnnotation()`
}
\seealso{
\code{\link[=sampleAnnotation]{sampleAnnotation()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_name.R
\name{tax_name}
\alias{tax_name}
\title{Set unique taxa_names for phyloseq object}
\usage{
tax_name(
  ps,
  prefix = c("tax", "asv", "otu")[1],
  rank = NA,
  pad_number = TRUE,
  sep = "_"
)
}
\arguments{
\item{ps}{phyloseq object}

\item{prefix}{e.g. 'tax', 'asv', or 'otu' (or set your own)}

\item{rank}{name of taxonomic rank from which to use classifications in new names}

\item{pad_number}{should unique numbers have zeros added to the front
(e.g. 001, 002) to be made the same number of characters?}

\item{sep}{character to separate the unique number and any taxonomic classification
info (relevant if rank given)}
}
\value{
phyloseq object
}
\description{
If your current taxa_names aren't what you want (e.g. they are long DNA sequences),
this function will help you set sensible unique names.

It combines:
\itemize{
\item a prefix like tax, asv, or otu (pick an appropriate prefix or set your own)
\item a unique (sequential) number
\item classification information from a chosen taxonomic rank (optional)
}
}
\details{
Don't confuse this with the phyloseq function taxa_names().
}
\examples{
library(phyloseq)
# get example data
data("enterotype")
ps <- enterotype
head(taxa_names(ps)) # these are mostly fine (except the -1), but imagine you wanted new names

# consider storing the original names for reference (e.g. if they are DNA sequences)
old_taxa_names <- taxa_names(ps)

ps <- tax_name(ps)
taxa_names(ps) \%>\% head()

# probably better to include the genus info to make these names more informative
ps <- tax_name(ps, rank = "Genus")
taxa_names(ps) \%>\% head()

# store new names with old names in dataframe for reference
names_df <- tibble::tibble(old = old_taxa_names, new = taxa_names(ps))

# alternative settings
tax_name(ps, pad_number = FALSE) \%>\%
  taxa_names() \%>\%
  head()
tax_name(ps, prefix = "whateveryoulike") \%>\%
  taxa_names() \%>\%
  head()
tax_name(ps, rank = "Genus", sep = "-") \%>\%
  taxa_names() \%>\%
  head()
}
\seealso{
\code{phyloseq::\link[phyloseq]{taxa_names}} for accessing and manually setting names
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_filter.R
\name{tax_filter}
\alias{tax_filter}
\title{Filter rare and/or low abundance taxa from a phyloseq object}
\usage{
tax_filter(
  ps,
  min_prevalence = 1,
  prev_detection_threshold = 1,
  min_total_abundance = 0,
  min_sample_abundance = 0,
  tax_level = NA,
  names_only = FALSE,
  use_counts = TRUE,
  undetected = NULL,
  verbose = TRUE
)
}
\arguments{
\item{ps}{phyloseq or ps_extra (ideally with count data available)}

\item{min_prevalence}{number or proportion of samples that a taxon must be present in}

\item{prev_detection_threshold}{min required counts (or value) for a taxon to be considered present
in that sample (or set undetected arg)}

\item{min_total_abundance}{minimum total readcount of a taxon, summed across all samples
(can be proportion of all counts)}

\item{min_sample_abundance}{taxa must have at least this many reads in one or more samples
(or proportion of that sample's reads)}

\item{tax_level}{if given, aggregates data at named taxonomic rank before filtering,
but returns phyloseq at the ORIGINAL level of aggregation!}

\item{names_only}{if names_only is true return only names of taxa, not the phyloseq}

\item{use_counts}{expect count data in phyloseq otu_table? default is TRUE}

\item{undetected}{e.g. 0, value at (or below) which a taxon is considered not present in that sample.
If set, this overrides prev_detection_threshold.}

\item{verbose}{message about proportional prevalence calculations?}
}
\value{
filtered phyloseq object AT ORIGINAL LEVEL OF AGGREGATION
(not at the level in tax_level)
}
\description{
Removes taxa (from all samples) that do not meet a given criterion or combination of criteria.
If a value for min_prevalence, min_total_abundance or min_sample_abundance is 1 or greater, then
it is treated as an absolute minimum number of samples/reads. If <1, it is treated as proportion of all samples/reads.
This function is designed to work with counts. otu_table must contain counts particularly if you want to set a non-zero value for min_total_abundance.
}
\examples{
data("dietswap", package = "microbiome")
# Dropping rare and low abundance taxa #
# Filter at unique taxa level, keeping only those with a prevalence of 10\% or more
# and at least 10 thousand reads when summed across all samples.
# Then aggregate to Family level taxonomy.
dietswap \%>\%
  tax_filter(min_prevalence = 0.1, min_total_abundance = 10000) \%>\%
  tax_agg("Family")

# Keeping ubiquitous families #
# keep only families that have at least 1000 counts present in 90\% of samples
# then aggregate the remaining taxa at 'Genus' level
dietswap \%>\%
  tax_filter(
    tax_level = "Family", min_prevalence = 0.90,
    prev_detection_threshold = 1000
  ) \%>\%
  tax_agg("Genus")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_seriate.R
\name{ps_seriate}
\alias{ps_seriate}
\title{Arrange samples in a phyloseq by microbiome similarity}
\usage{
ps_seriate(
  ps,
  method = "OLO_ward",
  dist = "bray",
  tax_transform = "identity",
  add_variable = FALSE
)
}
\arguments{
\item{ps}{phyloseq object}

\item{method}{seriation method for ordering samples, from seriation::seriate}

\item{dist}{distance method for dist_calc (only used if required for particular seriation method!)}

\item{tax_transform}{transformation to apply before seriation or any distance calculation}

\item{add_variable}{add a variable to the sample data indicating seriation order}
}
\value{
phyloseq
}
\description{
Uses seriation methods from seriation::seriate and often dist_calc (depending on if seriation method requires a distance matrix)
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")

dietswap \%>\%
  sample_data() \%>\%
  head(8)

dietswap \%>\%
  tax_agg("Genus") \%>\%
  .$ps \%>\%
  ps_seriate(method = "OLO_ward", dist = "bray") \%>\%
  sample_data() \%>\%
  head(8)
}
\seealso{
\code{\link{ps_arrange}} \code{\link{ps_reorder}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_names2rank.R
\name{tax_names2rank}
\alias{tax_names2rank}
\alias{tax_names2tt}
\title{Add taxa_names as last column in phyloseq tax_table}
\usage{
tax_names2rank(data, colname = "unique")

tax_names2tt(data, colname = "unique")
}
\arguments{
\item{data}{phyloseq object, or ps_extra or tax_table (taxonomyTable)}

\item{colname}{name of new rank to add at right side of tax_table}
}
\value{
same class object as passed in to data
}
\description{
The taxa names in your phyloseq may specify a further unique classification
of your taxa, e.g. ASVs, that is not otherwise represented in the tax_table itself.
This function fixes that, and allows you to include this level in taxatree_plots for example.
}
\details{
\code{tax_names2tt} is the old name of \code{tax_names2rank}.
The use of \code{tax_names2tt} is deprecated.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-varAnnotation.R
\name{anno_var_hist}
\alias{anno_var_hist}
\title{Helper to specify heatmap annotation for variable distribution density plot}
\usage{
anno_var_hist(
  fun = identity,
  size = grid::unit(30, "mm"),
  n_breaks = 11,
  border = FALSE,
  gp = grid::gpar(fill = "#CCCCCC"),
  axis = TRUE,
  ...,
  data = NULL,
  vars = NULL,
  which = NULL
)
}
\arguments{
\item{fun}{function applied to all variables, with apply()}

\item{size}{width or height as a grid unit object}

\item{n_breaks}{number of breaks}

\item{border}{Wether draw borders of the annotation region?}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_density]{ComplexHeatmap::anno_density}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{vars}{OPTIONAL selection vector of variable names,
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to varAnnotation(),
which itself is used by cor_heatmap var_anno argument.
}
\examples{
library(ComplexHeatmap)
set.seed(123)
fakeData <- as.data.frame.matrix(matrix(rnorm(500, 10, 3), ncol = 10))
names(fakeData) <- paste0("var_", 1:10)

# draw the histograms without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(
  anno_var_hist(data = fakeData, vars = names(fakeData), which = "row")
)

grid.newpage()
pushViewport(vp)
draw(
  anno_var_hist(
    data = fakeData, fun = sqrt,
    vars = rev(names(fakeData)), n_breaks = 5,
    which = "column", gp = grid::gpar(fill = 2:6, lwd = c(0.9, 2.5))
  )
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_bdisp.R
\name{dist_bdisp}
\alias{dist_bdisp}
\title{Wrapper for vegan::betadisper()}
\usage{
dist_bdisp(
  data,
  variables,
  method = c("centroid", "median")[[1]],
  complete_cases = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{data}{ps_extra output from dist_calc}

\item{variables}{list of variables to use as group}

\item{method}{centroid or median}

\item{complete_cases}{drop samples with NAs in any of the variables listed}

\item{verbose}{sends messages about progress if true}
}
\value{
ps_extra list containing betadisper results
}
\description{
Takes the output of dist_calc function. Or use with the result of the permanova function to ensure the results correspond to exactly the same input data.
Runs betadisper for all categorical variables in variables argument.
See help('betadisper', package = 'vegan').
}
\examples{
library(phyloseq)
library(vegan)
data("dietswap", package = "microbiome")

# add some missings to demonstrate automated removal
sample_data(dietswap)$sex[3:6] <- NA
# create a numeric variable to show it will be skipped with a warning
dietswap <- ps_mutate(dietswap, timepoint = as.numeric(timepoint))

# straight to the betadisp
bd1 <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("aitchison") \%>\%
  dist_bdisp(variables = c("sex", "bmi_group", "timepoint")) \%>\%
  bdisp_get()
bd1$sex
# quick vegan plotting methods
plot(bd1$sex$model, label.cex = 0.5)
boxplot(bd1$sex$model)

# compute distance and use for both permanova and dist_bdisp
testDist <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("bray")

PERM <- testDist \%>\%
  dist_permanova(
    variables = c("sex", "bmi_group"),
    n_processes = 1, n_perms = 99
  )
str(PERM, max.level = 1)

bd <- PERM \%>\% dist_bdisp(variables = c("sex", "bmi_group"))
bd
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-varAnnotation.R
\name{anno_var_box}
\alias{anno_var_box}
\title{Helper to specify heatmap annotation for showing variable distributions}
\usage{
anno_var_box(
  fun = identity,
  size = grid::unit(30, "mm"),
  border = TRUE,
  gp = grid::gpar(fill = "#CCCCCC"),
  ylim = NULL,
  extend = 0.05,
  outline = TRUE,
  box_width = 0.6,
  pch = 1,
  pointsize = grid::unit(0.5, "mm"),
  axis = TRUE,
  ...,
  data = NULL,
  vars = NULL,
  which = NULL
)
}
\arguments{
\item{fun}{function applied to all variables, with apply()}

\item{size}{width or height as a grid unit object}

\item{border}{Wether draw borders of the annotation region?}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{ylim}{Data ranges. By default it is \code{range(x)} if \code{x} is a vector, or \code{range(rowSums(x))} if \code{x} is a matrix.}

\item{extend}{The extension to both side of \code{ylim}. The value is a percent value corresponding to \code{ylim[2] - ylim[1]}.}

\item{outline}{Whether draw outline of boxplots?}

\item{box_width}{Relative width of boxes. The value should be smaller than one.}

\item{pch}{Point style.}

\item{pointsize}{size of outlier points, as grid::unit() object}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_boxplot]{ComplexHeatmap::anno_boxplot}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{vars}{OPTIONAL selection vector of variable names,
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to varAnnotation(),
which itself is used by cor_heatmap as var_anno() argument.
}
\examples{
library(ComplexHeatmap)
set.seed(123)
fakeData <- as.data.frame.matrix(matrix(rnorm(500, 10, 3), ncol = 10))
names(fakeData) <- paste0("var_", 1:10)

# draw the boxplot without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(
  anno_var_box(data = fakeData, vars = names(fakeData), which = "column")
)

grid.newpage()
pushViewport(vp)
draw(
  anno_var_box(
    data = fakeData, fun = function(x) log(x + 1),
    vars = rev(names(fakeData)),
    which = "row"
  )
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-taxAnnotation.R
\name{anno_tax_box}
\alias{anno_tax_box}
\title{Helper to specify heatmap annotation for showing taxa abundance on boxplot}
\usage{
anno_tax_box(
  undetected = 0,
  only_detected = TRUE,
  trans = "compositional",
  zero_replace = 0,
  use_counts = TRUE,
  size = grid::unit(30, "mm"),
  border = TRUE,
  gp = grid::gpar(fill = "#CCCCCC"),
  ylim = NULL,
  extend = 0.05,
  outline = TRUE,
  box_width = 0.6,
  pch = 1,
  pointsize = grid::unit(0.5, "mm"),
  axis = TRUE,
  ...,
  data = NULL,
  taxa = NULL,
  which = NULL
)
}
\arguments{
\item{undetected}{the value above which taxa are classed as detected/present in a sample}

\item{only_detected}{only plot values for samples where the taxon abundance is > undetected}

\item{trans}{name of transformation suitable for tax_transform,
or a function calling tax_transform, and/or tax_scale,
(a function must take a phyloseq or ps_extra, and return one)}

\item{zero_replace}{zero_replace value for for tax_transform, ignored if trans is a function}

\item{use_counts}{try to retrieve counts from data object?}

\item{size}{width or height as a grid unit object}

\item{border}{Wether draw borders of the annotation region?}

\item{gp}{Graphic parameters for points. The length of each graphic parameter can be 1, length of \code{x} if \code{x} is a vector, or number of columns of \code{x} is \code{x} is a matrix.}

\item{ylim}{Data ranges. By default it is \code{range(x)} if \code{x} is a vector, or \code{range(rowSums(x))} if \code{x} is a matrix.}

\item{extend}{The extension to both side of \code{ylim}. The value is a percent value corresponding to \code{ylim[2] - ylim[1]}.}

\item{outline}{Whether draw outline of boxplots?}

\item{box_width}{Relative width of boxes. The value should be smaller than one.}

\item{pch}{Point style.}

\item{pointsize}{size of outlier points, as grid::unit() object}

\item{axis}{Whether to add axis?}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:anno_boxplot]{ComplexHeatmap::anno_boxplot}}
  \describe{
    \item{\code{axis_param}}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  }}

\item{data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{taxa}{OPTIONAL selection vector of taxa (names, numbers or logical),
only set this if providing data argument to override default}

\item{which}{OPTIONAL indicating if it is a 'column' or a 'row' annotation,
only set this if providing data argument to override default}
}
\value{
function or ComplexHeatmap AnnotationFunction object
}
\description{
Use this as an argument to taxAnnotation(),
which itself is used by cor_heatmap and comp_heatmap as tax_anno argument.
}
\examples{
library("ComplexHeatmap")
data("ibd_phylo", package = "corncob")
psq <- tax_filter(ibd_phylo, min_prevalence = 5)
psq <- tax_mutate(psq, Species = NULL)
psq <- tax_fix(psq)
psq <- tax_agg(psq, rank = "Family")
taxa <- tax_top(psq, n = 15, rank = "Family")
# makes a function that takes data, taxa and which (at minimum)
fun <- anno_tax_box()
# manually specify the prevalence barplot function by giving it data etc.
heatmapAnnoFunction <- fun(data = psq, which = "column", taxa = taxa)
# draw the barplot without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)
grid.newpage()
pushViewport(vp)
draw(heatmapAnnoFunction)

# let's change some style options and specify the data up front
grid::grid.newpage()
pushViewport(vp)
draw(anno_tax_box(
  data = psq, taxa = taxa, which = "row", pointsize = grid::unit(1, "mm"),
  gp = grid::gpar(fill = "red"), border = FALSE, box_width = 0.2
))

# clear drawings
grid::grid.newpage()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_select.R
\name{tax_select}
\alias{tax_select}
\title{Subset phyloseq object by (partial) taxa names}
\usage{
tax_select(
  ps,
  tax_list,
  ranks_searched = "all",
  strict_matches = FALSE,
  n_typos = 1,
  deselect = FALSE
)
}
\arguments{
\item{ps}{phyloseq object}

\item{tax_list}{e.g. c('g__Bifidobacterium', 'g__Akkermansia', 'g__Bacteroides', 'g__Streptococcus')}

\item{ranks_searched}{'all' or a list of which taxonomic ranks should be searched for the names in tax_list?}

\item{strict_matches}{only perfect full name matches allowed if TRUE}

\item{n_typos}{how many typos to allow in each name? uses agrep approximate matching if > 0}

\item{deselect}{if TRUE, the matching taxa will be REMOVED instead!}
}
\value{
phyloseq object with fewer taxa
}
\description{
Convenient name-based taxa selection/filtering of phyloseq object, including approximate name matching.
Takes a phyloseq with tax table and a (partial) taxonomic name, or a list/vector of taxonomic names (full or partial matches).
}
\details{
tax_select will also search the otu names/rownames, BUT only for perfect matches.
}
\examples{
# Get example phyloseq object data
data("dietswap", package = "microbiome")
pSeq <- dietswap

# SELECTION EXAMPLES #
a <- pSeq \%>\% tax_select(tax_list = "Bif", n_typos = 0, ranks_searched = "Genus")
b <- pSeq \%>\% tax_select(tax_list = "Bifidobacterium", n_typos = 0)
c <- pSeq \%>\% tax_select(tax_list = "Bif", n_typos = 1)
identical(a, b) # TRUE
identical(a, c) # FALSE

pSeq \%>\% tax_select(tax_list = "Bifidobactrium") # default 1 typo allowed
one <- pSeq \%>\% tax_select(tax_list = "Akkarmensia", n_typos = 2)
two <- pSeq \%>\% tax_select(tax_list = "Akkermansia", n_typos = 0)
identical(one, two) # TRUE

# DESELECTION EXAMPLE # #
pSeq \%>\% tax_select(tax_list = "Bif", strict_matches = FALSE, deselect = TRUE)
# Incorrect example
# pSeq \%>\% tax_select(tax_list = "Bif", strict_matches = TRUE) # fails
}
\seealso{
\code{\link{ps_select}} for selecting variables in phyloseq sample_data

\code{\link{agrep}} for the function that powers the approximate matching in tax_select
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_reorder.R
\name{ps_reorder}
\alias{ps_reorder}
\title{Set order of samples in phyloseq object}
\usage{
ps_reorder(ps, sample_order)
}
\arguments{
\item{ps}{phyloseq}

\item{sample_order}{names or current numerical indices of samples in desired order}
}
\value{
phyloseq
}
\description{
Manually set order of samples by specifying samples names in desired order.
}
\details{
Ordering of samples in a phyloseq is controlled from the otu_table slot!
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")

dietswap \%>\%
  sample_data() \%>\%
  head(8)

new_order <- rev(sample_names(dietswap))
dietswap \%>\%
  ps_reorder(new_order) \%>\%
  sample_data() \%>\%
  head(8)

# random ordering with numbers
set.seed(1000)
random_order <- sample(1:nsamples(dietswap))
dietswap \%>\%
  ps_reorder(random_order) \%>\%
  sample_data() \%>\%
  head(8)
}
\seealso{
\code{\link{ps_arrange}} for arranging samples by sample_data variables (or otu_table)

\code{\link{ps_seriate}} for arranging samples by microbiome similarity

\code{\link{ps_filter}} for keeping only some samples, based on sample_data
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_calc.R
\name{ord_calc}
\alias{ord_calc}
\title{Ordinate samples (arrange by similarity in multiple dimensions)}
\usage{
ord_calc(
  data,
  method = c("auto", "PCoA", "PCA", "CCA", "RDA", "CAP", "NMDS")[1],
  constraints = NULL,
  conditions = NULL,
  scale_cc = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{data}{ps_extra list object: output from dist_calc(), or tax_transform() if no distance calculation required for method e.g. for RDA}

\item{method}{which ordination method to use? "auto" means automatically determined from ps_extra and other args.
If you really know what you want: manually set one of 'PCoA', 'PCA', 'CCA', 'CAP' or 'RDA'}

\item{constraints}{(a vector of) valid sample_data name(s) to constrain analyses, or leave as NULL for unconstrained ordination.
Non-NULL values are compatible with method = "auto"/"RDA"/"CAP"}

\item{conditions}{(a vector of) valid sample_data name(s) to partial these out of analyses with Condition(), or leave as NULL}

\item{scale_cc}{If TRUE (default) ensures any constraints and conditions variables are scaled before use, to ensure their effects are comparable.
If set to FALSE you must ensure you have already set the variables on a similar scale yourself!
If there are no constraints or conditions, scale_cc does nothing.}

\item{verbose}{If TRUE or "max", show any warnings and messages about constraint and conditions scaling and missings etc.
FALSE suppresses warnings!}

\item{...}{optional arguments passed on to phyloseq::ordinate()}
}
\value{
ps_extra list object
}
\description{
Used before plotting with ord_plot() or explorating interactively with ord_explore().
Use method = "auto" to automatically pick an appropriate method from:
\itemize{
\item "PCA" (Principle Components Analysis) combines taxa abundances into new dimensions. The first axes display the greatest variation in your microbial data.
\item "RDA" (Redundancy Analysis) is constrained PCA, roughly speaking. It finds variation in your data that can be explained both by the constraints variables, and the microbial data.
\item "PCoA" (Principle Coordinates Analysis) finds a coordinate system that best preserves the original distances between samples.
\item "CAP" (Constrained Analysis of Principle Coordinates) is also known as distance-based Redundancy Analysis.
}

Alternatively to leaving method = "auto", you can explicitly specify any of the above methods, or choose one of the following:
\itemize{
\item "CCA" (Canonical Correspondence Analysis) - NOT canonical correlation analysis!
\item "NMDS" (Non-metric Multidimensional Scaling)
}

You are strongly recommended to check out this useful website for introductory explanations of these methods
the "GUide to STatistical Analysis in Microbial Ecology": \url{https://sites.google.com/site/mb3gustame/}
}
\details{
Extends functionality of phyloseq::ordinate().
Results can be used directly in ord_plot().
You can extract the ordination object for other analyses with ord_get()
}
\examples{
library(phyloseq)
library(vegan)
data("dietswap", package = "microbiome")

# create a couple of numerical variables to use as constraints
dietswap <- ps_mutate(
  dietswap,
  female = dplyr::if_else(sex == "female", true = 1, false = 0),
  weight = dplyr::recode(bmi_group, obese = 3, overweight = 2, lean = 1)
)

# add a couple of missing values to demo automated dropping of observations with missings
sample_data(dietswap)$female[c(3, 4)] <- NA

# compute ordination
test <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("bray") \%>\%
  ord_calc(constraints = c("weight", "female"))
# familiarise yourself with the structure of the returned ps_extra list object
test
str(test, max.level = 1)

# compute RDA with centre-log-ratio transformed taxa
test2 <- dietswap \%>\%
  tax_agg("Genus") \%>\%
  tax_transform("clr") \%>\%
  ord_calc(constraints = c("weight", "female"))
# plot with vegan package graphics to show it returns a standard ordination object
ord_get(test2) \%>\% vegan::ordiplot()
ord_plot(test2, plot_taxa = 8:1)
# This is equivalent to CAP with "aitchison" distance
# but the latter (below) doesn't allow plotting taxa loadings with ord_plot
dietswap \%>\%
  tax_agg("Genus") \%>\%
  dist_calc("aitchison") \%>\%
  ord_calc(constraints = c("weight", "female")) \%>\%
  ord_plot()
}
\seealso{
\code{\link{dist_calc}} for distance matrix calculation

\code{\link{ord_plot}} and \code{\link{ord_explore}}

\code{phyloseq \link[phyloseq]{ordinate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microViz.R
\docType{package}
\name{microViz}
\alias{microViz}
\title{microViz: microbiome data analysis and visualization}
\description{
microViz provides functions for statistics and visualization of microbiome sequencing data.
microViz wraps, extends and complements popular microbial ecology packages like phyloseq, vegan, and microbiome.

Check out the website for tutorials and illustrated help pages.

\url{https://david-barnett.github.io/microViz/}
}
\author{
David Barnett
(\href{https://orcid.org/0000-0003-1961-7206}{ORCID})
(\href{https://github.com/david-barnett}{GitHub})
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_drop_incomplete.R
\name{ps_drop_incomplete}
\alias{ps_drop_incomplete}
\title{Deselect phyloseq samples with sample_data missings}
\usage{
ps_drop_incomplete(ps, vars = NA, verbose = FALSE)
}
\arguments{
\item{ps}{phyloseq with sample_data}

\item{vars}{vector of variable names to check for missings (or NA, which uses all variables in sample data)}

\item{verbose}{message about number of samples dropped if verbose not FALSE, (and only if > 0 samples dropped)
and message about number of missing per variable in vars if verbose = "max" (and message even if 0 samples dropped)}
}
\value{
phyloseq
}
\description{
Check phyloseq object sample_data for missing values (NAs)
\itemize{
\item specify which variables to check with vars argument, or check all
\item drop samples with any missings
}
}
\details{
This is a wrapper for \code{\link[stats:complete.cases]{stats::complete.cases}} function.
}
\examples{
library(phyloseq)
data("enterotype", package = "phyloseq")

enterotype
ps_drop_incomplete(enterotype)
ps_drop_incomplete(enterotype, vars = "Enterotype", verbose = TRUE)
ps_drop_incomplete(enterotype, vars = "Sample_ID", verbose = TRUE)
ps_drop_incomplete(enterotype, vars = c("Enterotype", "Sample_ID"))
ps_drop_incomplete(enterotype, verbose = "max")
}
\seealso{
\code{\link{ps_filter}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-sampleAnnotation.R
\name{anno_cat}
\alias{anno_cat}
\title{Create colored rectangle annotations for categorical data}
\usage{
anno_cat(
  x,
  which,
  renamer = identity,
  col = distinct_palette(),
  width = NULL,
  height = NULL,
  box_col = "white",
  box_lwd = 0.5,
  border_col = NA,
  border_lwd = 1,
  legend = TRUE,
  legend_title = ""
)
}
\arguments{
\item{x}{data vector, treated as categorical}

\item{which}{Whether it is a column annotation or a row annotation?}

\item{renamer}{function renaming variable values for legend}

\item{col}{colors vector, at least as long as unique(x), optionally named by x levels}

\item{width}{grid unit object or NULL}

\item{height}{grid unit object or NULL}

\item{box_col}{colour of boxes around individual cells}

\item{box_lwd}{line width of boxes around individual cells}

\item{border_col}{colour of border around all cells}

\item{border_lwd}{line width of border around all cells}

\item{legend}{generate legend for this annotation
(attached as attribute of heatmap, and not automatically included in plot)}

\item{legend_title}{title for legend, if drawn}
}
\value{
AnnotationFunction
}
\description{
Similar to anno_simple but with individual boxes!
}
\examples{
library(ComplexHeatmap)
# draw the annotation without a heatmap, you will never normally do this!
vp <- viewport(width = 0.75, height = 0.75)

grid::grid.newpage()
pushViewport(vp)
cats <- letters[1:4]
draw(anno_cat(cats, which = "row"))

grid::grid.newpage()
pushViewport(vp)
draw(
  anno_cat(
    x = cats, col = structure(names = cats, 1:4), which = "column",
    box_col = "black", box_lwd = 5
  )
)

# developer note #
# list of annotations can be split and ordered (adding NULL makes a list)
# https://jokergoo.github.io/ComplexHeatmap-reference/book/a-list-of-heatmaps.html
# (section #4.8 concatenate-only-the-annotations)
grid::grid.newpage()
pushViewport(vp)
annoList <- rowAnnotation(
  hi = anno_cat(cats, which = "row", border_col = "black")
) +
  NULL
draw(object = annoList, row_split = c(1, 1:3), row_order = 4:1)
pushViewport(viewport(x = 0.6))
draw(anno_cat(cats, "row", legend_title = "abcd") \%>\% attr("Legend"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps-varAnnotation.R
\name{varAnnotation}
\alias{varAnnotation}
\title{Helper to specify a HeatmapAnnotation for variables in cor_heatmap}
\usage{
varAnnotation(
  ...,
  name,
  annotation_legend_param = list(),
  show_legend = TRUE,
  gp = grid::gpar(col = NA),
  border = FALSE,
  gap = grid::unit(2, "mm"),
  show_annotation_name = TRUE,
  annotation_label = NULL,
  annotation_name_gp = grid::gpar(),
  annotation_name_offset = NULL,
  annotation_name_rot = NULL,
  annotation_name_align = FALSE,
  annotation_name_side = "auto",
  .data = NULL,
  .vars = NULL,
  .side = NULL
)
}
\arguments{
\item{...}{Name-value pairs where the names correspond to annotation names and values
are the output of variable annotation functions
such as anno_var_box(), or manually specified AnnotationFunction objects}

\item{name}{Name of the heatmap annotation, optional.}

\item{annotation_legend_param}{A list which contains parameters for annotation legends. See \code{\link[ComplexHeatmap]{color_mapping_legend,ColorMapping-method}} for all possible options.}

\item{show_legend}{Whether show annotation legends. The value can be one single value or a vector.}

\item{gp}{Graphic parameters for simple annotations (with \code{fill} parameter ignored).}

\item{border}{border of single annotations.}

\item{gap}{Gap between annotations. It can be a single value or a vector of \code{\link[grid]{unit}} objects.}

\item{show_annotation_name}{Whether show annotation names? For column annotation, annotation names are drawn either on the left or the right, and for row annotations, names are draw either on top or at the bottom. The value can be a vector.}

\item{annotation_label}{Labels for the annotations. By default it is the same as individual annotation names.}

\item{annotation_name_gp}{Graphic parameters for anntation names. Graphic paramters can be vectors.}

\item{annotation_name_offset}{Offset to the annotation names, a \code{\link[grid]{unit}} object. The value can be a vector.}

\item{annotation_name_rot}{Rotation of the annotation names. The value can be a vector.}

\item{annotation_name_align}{Whether to align the annotation names.}

\item{annotation_name_side}{Side of the annotation names.}

\item{.data}{OPTIONAL phyloseq or ps_extra,
only set this to override use of same data as in heatmap}

\item{.vars}{OPTIONAL selection vector of variables (names, numbers or logical),
only set this if providing .data argument to override default}

\item{.side}{OPTIONAL string, indicating the side for the variable annotations:
only set this to override default}
}
\value{
HeatmapAnnotation object
}
\description{
Helper to specify a HeatmapAnnotation for variables in cor_heatmap
}
\seealso{
\code{\link[=taxAnnotation]{taxAnnotation()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ord_explore.R
\name{ord_explore}
\alias{ord_explore}
\title{Interactively explore microbial compositions of ordinated samples}
\usage{
ord_explore(
  data,
  sample_id = NULL,
  seriate_method = "OLO_ward",
  app_options = list(launch.browser = TRUE),
  plot_widths = c(7, 9),
  ...
)
}
\arguments{
\item{data}{ps_extra list output of ord_calc, or phyloseq}

\item{sample_id}{name of sample ID variable to use as default for selecting samples}

\item{seriate_method}{seriation method to order phyloseq samples by similarity}

\item{app_options}{passed to shinyApp() options argument}

\item{plot_widths}{widths of plots in inches, including any legends
(first number is ordination, second is composition barplot)}

\item{...}{additional arguments passed to ord_plot}
}
\value{
nothing, opens default browser
}
\description{
A Shiny app used to create and explore an interactive version of \code{ord_plot()}.
You can select samples on an ordination plot to view their composition with stacked barplots.

The \code{ord_explore()} data argument takes either:
\itemize{
\item the output of \code{ord_calc()} (i.e. a ps_extra with an ordination)
\item a plain phyloseq object: \code{ord_explore()} will help you build an ordination
}

Once the app is running (in your browser), you can:
\enumerate{
\item Create/edit the ordination if required
\itemize{
\item look at the R console error messages if your chosen options don't build
}
\item Style the ordination plot (e.g. choose dimensions; set colour and size; ...)
\itemize{
\item Taxa loading arrows can be added only to PCA, RDA and CCA plots
\item Convex hulls or ellipses can only be drawn if Colour is set to a variable
\item To track individuals over time with the path plotter, your data MUST already be sorted by time (e.g. with ps_arrange)!
}
\item Click on or use the lasso tool to select 2 or more samples to view their compositions
\itemize{
\item By default samples can be selected individually
\item Set the "Select" option to another variable to select by level of that variable
}
\item Style the taxonomic compositions barplot
\itemize{
\item The samples are ordered using the seriate_method argument and the same transformation and distance as used in the ordination plot
\item The app may lag if you select 100s of samples and ungroup the "other" category
\item To avoid this lag: either reduce the number of taxa or samples, or deselect "Interactive" barplot
}
\item Stop the app by clicking the red stop button in the R console
\itemize{
\item Closing the web browser window doesn't stop the app,
(you can find the app again at the local http address shown in the R console)
\item Don't forget to copy the ordination plot code before you close the app
}
}

See the Details section for some known limitations of the app.
Please report any other app problems on the microViz GitHub issues page.
}
\details{
Limitations:
\itemize{
\item If a "Select:" grouping variable is NA for some samples,
then that grouping variable cannot be used to select those samples
\item "Shape:" can only be mapped to variables with maximum 5 distinct levels,
not including NAs. NAs in the shape variable are shown as hollow circles.
}

On some web browsers, e.g. older versions of firefox, the numeric inputs'
buttons are sometimes hard to click.
As a workaround, click the box and type a number or use the arrow keys.
This problem occurs in all Shiny apps, not just microViz.
}
\examples{
# example code only runs in interactive R session
if (interactive()) {
  library(phyloseq)
  library(dplyr)

  # example of quickstart approach with interactive ordination calculation #
  corncob::ibd_phylo \%>\%
    # filtering makes subsequent calculations faster
    tax_filter(min_prevalence = 2) \%>\%
    tax_fix() \%>\%
    ord_explore()

  # simple example with precalculated ordination #
  data("enterotype")
  taxa_names(enterotype)[1] <- "unclassified" # replaces the "-1" taxon name
  ps <- tax_fix(enterotype) # remove NA taxa
  ord1 <- ps \%>\%
    tax_transform("identity", rank = "Genus") \%>\%
    dist_calc("bray") \%>\%
    ord_calc(method = "PCoA")

  ord_explore(data = ord1, auto_caption = 6)

  # constrained ordination example #
  data("dietswap", package = "microbiome")

  # create a couple of numerical variables to use as constraints
  dietswap <- dietswap \%>\%
    ps_mutate(
      weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
      female = if_else(sex == "female", true = 1, false = 0)
    ) \%>\%
    tax_agg("Genus")

  constrained_aitchison_rda <- dietswap \%>\%
    tax_transform("clr") \%>\%
    ord_calc(constraints = c("weight", "female"))

  # label style arguments can be passed to ord_explore
  constrained_aitchison_rda \%>\%
    ord_explore(
      tax_lab_style = list(size = 3),
      constraint_lab_style = list(size = 4), auto_caption = 6
    )
  # Try changing the point colour to bmi_group or similar
  # Style points interactively!
  # (setting colour/shape/etc as arguments doesn't work)

  # dietswap is actually a longitudinal dataset, with multiple samples per
  # subject. If we arrange by timepoint first (!!!), we can use the "paths"
  # additional plot layer from the ord_explore "Add:" menu to track
  # individual subjects over time.
  dietswap \%>\%
    ps_arrange(timepoint) \%>\%
    tax_fix() \%>\%
    ord_explore()


  # Another dataset, where "size" variable drives gradient on PC1
  # Try setting size and/or alpha to correspond to "size"!
  # Then edit the ordination to use "size" as a condition, see what happens
  # hmp2 <- microbiomeutilities::hmp2
  hmp2 \%>\%
    tax_fix() \%>\%
    tax_transform(rank = "Genus", "identity") \%>\%
    dist_calc("aitchison") \%>\%
    ord_calc() \%>\%
    ord_explore()

  # another dataset
  data("soilrep", package = "phyloseq")
  # test auto creation of SAMPLE var
  ps <- soilrep \%>\% ps_select(-Sample)
  # The barplot is actually quite useless with the 16000+ anonymous OTUs
  # in this dataset, but the 1000s of unmerged "other" categories do render
  phyloseq_validate(ps) \%>\%
    tax_fix() \%>\%
    dist_calc("aitchison") \%>\%
    ord_calc() \%>\%
    ord_explore()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_mutate.R
\name{tax_mutate}
\alias{tax_mutate}
\title{Modify or compute new taxonomic ranks in phyloseq}
\usage{
tax_mutate(ps, ...)
}
\arguments{
\item{ps}{phyloseq object with a tax_table, or just a tax_table}

\item{...}{passed straight to dplyr::mutate (see examples and dplyr::mutate help)}
}
\value{
phyloseq object with modified tax_table
}
\description{
Add or overwrite tax_table ranks.
Use dplyr::mutate() syntax.
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")

# compute new rank
tax_mutate(dietswap, loud_genus = toupper(Genus)) \%>\%
  tt_get() \%>\%
  head()

# overwrite a current rank
tax_mutate(dietswap, Genus = toupper(Genus)) \%>\%
  tt_get() \%>\%
  head()

# overwrite all ranks
tax_mutate(dietswap, dplyr::across(.fns = toupper)) \%>\%
  tt_get() \%>\%
  head()

# add a new rank at the beginning
tax_mutate(dietswap, Root = "Bacteria", .before = 1) \%>\%
  tt_get() \%>\%
  head()

# this is an error as ranks can't be any other class than character
# tax_mutate(dietswap, Genus = 1:ntaxa(dietswap))
}
\seealso{
\code{\link[dplyr]{mutate}}

\code{\link{ps_mutate}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxatree_label.R
\name{taxatree_label}
\alias{taxatree_label}
\title{Add logical label column to taxatree_stats dataframe}
\usage{
taxatree_label(
  data,
  ...,
  .label_var = "label",
  .node_fun = list(prevalence = prev)
)
}
\arguments{
\item{data}{ps_extra (or phyloseq)}

\item{...}{REQUIRED logical conditions for labelling
e.g. rank == "Phylum", p.value < 0.1 | taxon \%in\% listOfTaxa}

\item{.label_var}{name of label indicator variable to be created.
If you change this, beware that taxatree_plotkey will not work, you will
need to called taxatree_plot_label with}

\item{.node_fun}{named list of length 1 providing \code{taxatree_nodes} \code{fun} arg.
(name of list iterm is available for use in ...)}
}
\value{
ps_extra with (modified) taxatree_stats dataframe
}
\description{
\code{taxatree_label} is used internally by \code{taxatree_plotkey}, but
can also be used prior to \code{taxatree_plots} to label those plots directly

\code{...} arguments are passed to \code{dplyr::filter()},
so all \code{filter} syntax can be used.
}
\details{
If taxatree_stats missing (or if data is a phyloseq)
it will create a plain taxatree_stats dataframe using only taxatree_nodes

\code{node_fun} can also be a precalculated dataframe (output of taxatree_nodes)
but you should probably not use this option.
This is used internally for efficiency inside \code{taxatree_plotkey()}
}
\examples{
# simple example with plain phyloseq input
data("dietswap", package = "microbiome")
labelled <- dietswap \%>\%
  tax_prepend_ranks() \%>\%
  taxatree_label(rank == "Phylum", prevalence > 0.1)

# Note that "prevalence" column was available in data
# because it is created by `taxatree_nodes()` using the named function
# provided to the `node_fun` argument

# ps_Extra is returned
labelled

# notice how both conditions must be met for label column to be TRUE
labelled$taxatree_stats
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_filter.R
\name{ps_filter}
\alias{ps_filter}
\title{Filter phyloseq samples by sample_data variables}
\usage{
ps_filter(ps, ..., .target = "sample_data", .keep_all_taxa = FALSE)
}
\arguments{
\item{ps}{phyloseq object}

\item{...}{passed directly to dplyr::filter (see examples and ?dplyr::filter)}

\item{.target}{which slot of phyloseq to use for filtering by,
currently only "sample_data" supported}

\item{.keep_all_taxa}{if FALSE (the default),
remove taxa which are no longer present in the dataset after filtering}
}
\value{
phyloseq object (with filtered sample_data)
}
\description{
Keep only samples with sample_data matching one or more conditions.
Use this function as you would use use dplyr::filter(), but with a phyloseq object!
}
\examples{
library(phyloseq)
library(dplyr)

data("enterotype", package = "phyloseq")
enterotype
sample_data(enterotype)[1:10, 1:5]

# keep only samples with seqtech not equal to sanger
ps1 <- ps_filter(enterotype, SeqTech != "Sanger")
ps1
sample_data(ps1)[1:10, 1:5]

# keep only samples with no NAs in any samples
ps2 <- enterotype \%>\% ps_filter(across(everything(), ~ !is.na(.)))
ps2
sample_data(ps2)[1:8, 1:8]

# ps2 is equivalent to dropping samples with incomplete sample_variables and tax_filtering 0s
ps3 <- enterotype \%>\%
  ps_drop_incomplete() \%>\%
  tax_filter(prev_detection_threshold = 1e-20, use_counts = FALSE)
# we needed to set a low detection threshold because this example data is proportions
identical(ps2, ps3) # TRUE

# function will give warning if some of the otu_values are negative
# (which may happen when filtering data that has e.g. clr-transformed taxa abundances)
# as it attempts to discard any taxa that become always absent/0 after filtering (by default)
# set .keep_all_taxa = TRUE to avoid this filtering behaviour, which is unwanted in this case
enterotype \%>\%
  tax_transform("clr") \%>\%
  ps_get() \%>\%
  ps_filter(SeqTech == "Sanger", .keep_all_taxa = TRUE)
}
\seealso{
\code{\link[dplyr]{filter}} explains better how to give arguments to this function

\code{\link{tax_filter}} for filtering taxa (not samples)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_dedupe.R
\name{ps_dedupe}
\alias{ps_dedupe}
\title{De-duplicate phyloseq samples}
\usage{
ps_dedupe(ps, vars, method = "readcount", verbose = TRUE)
}
\arguments{
\item{ps}{phyloseq with sample data}

\item{vars}{names of variables, whose (combined) levels identify groups from which only 1 sample is desired}

\item{method}{keep 1 sample with max "readcount" or the "first" or "last" samples encountered in given sample_data order for each dupe group}

\item{verbose}{message names of samples removed if TRUE (or "debug" for more info, when method = "readcount")}
}
\value{
phyloseq object
}
\description{
Use 1 or more variables in the sample_data to identify and remove duplicate samples (leaving 1 per category).

\strong{methods:}
\itemize{
\item method = "readcount" keeps the one sample in each duplicate group with the highest total number of reads according to sample_sums
\item method = "first" keeps the first sample in each duplicate group encountered in the row order of the sample_data
\item method = "last" keeps the last sample in each duplicate group encountered in the row order of the sample_data
}
}
\details{
What happens when duplicated samples have exactly equal readcounts in method = "readcount"?
The first encountered maximum is kept (in sample_data row order, like method = "first")
}
\examples{
data("dietswap", package = "microbiome")

dietswap
# let's pretend the dietswap data contains technical replicates from each subject
# we want to keep only one of them
ps_dedupe(dietswap, vars = "subject", method = "readcount", verbose = TRUE)

# contrived example to show identifying "duplicates" via the interaction of multiple columns
ps1 <- ps_dedupe(
  ps = dietswap, method = "readcount", verbose = TRUE,
  vars = c("timepoint", "group", "bmi_group")
)
phyloseq::sample_data(ps1)

ps2 <- ps_dedupe(
  ps = dietswap, method = "first", verbose = TRUE,
  vars = c("timepoint", "group", "bmi_group")
)
phyloseq::sample_data(ps2)
}
\seealso{
\code{\link{ps_filter}} for filtering samples by sample_data variables
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scale_shape_girafe_filled.R
\name{scale_shape_girafe_filled}
\alias{scale_shape_girafe_filled}
\title{Filled shapes for ggiraph interactive plots}
\usage{
scale_shape_girafe_filled()
}
\value{
ggplot2 Scale object
}
\description{
Generates a custom ggplot2 shape scale, as used in ord_explore's ordination.
Uses filled shapes, therefore fill aesthetic must be set, in addition to
colour, to have filled shapes.
Points with NA values for the shape variable are shown as hollow circles.
}
\details{
Composite shapes e.g. number 7 "square cross" cause ggiraph interactive
plots to fail when a variable shape and tooltip is set.

Shapes used are, in order: "circle filled", "triangle filled",
"square filled", "diamond filled", and "triangle down filled"
}
\examples{
corncob::ibd_phylo \%>\%
  tax_fix() \%>\%
  phyloseq_validate() \%>\%
  tax_transform(rank = "Genus", trans = "clr") \%>\%
  ord_calc(
    method = "PCA"
  ) \%>\%
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:6,
    colour = "DiseaseState", fill = "DiseaseState",
    shape = "circle", alpha = 0.5,
    size = 3
  ) +
  scale_shape_girafe_filled()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dist_calc_seq.R
\name{dist_calc_seq}
\alias{dist_calc_seq}
\title{Calculate distances between sequential samples in ps_extra/phyloseq object}
\usage{
dist_calc_seq(
  data,
  dist,
  group,
  seq,
  unequal = "warn",
  start_value = NaN,
  return = "data",
  var_name = paste0(dist, "_DistFromLast")
)
}
\arguments{
\item{data}{ps_extra object, e.g. output from tax_transform()}

\item{dist}{name of distance to calculate between pairs of sequential samples}

\item{group}{name of variable in phyloseq sample_data used to define groups of samples}

\item{seq}{name of variable in phyloseq sample_data used to define order of samples
within groups}

\item{unequal}{"error" or "warn" or "ignore" if groups of samples, defined
by group argument, are of unequal size}

\item{start_value}{value returned for the first sample in each group, which
has no preceding sample in the group's sequence, and so has no obvious value}

\item{return}{format of return object: "data" returns ps_extra with sorted samples and
additional variable. "vector" returns only named vector of sequential
distances.}

\item{var_name}{name of variable created in ps_extra if return arg = "data"}
}
\value{
ps_extra object sorted and with new sequential distance variable
or a named vector of that variable
}
\description{
Calculate distances between sequential samples in ps_extra/phyloseq object
}
\examples{
library(ggplot2)
library(dplyr)
data("dietswap", package = "microbiome")

pseq <- dietswap \%>\%
  tax_transform("identity", rank = "Genus") \%>\%
  dist_calc_seq(
    dist = "aitchison", group = "subject", seq = "timepoint",
    # group sizes are unequal because some subjects are missing a timepoint
    unequal = "ignore"
  )

pseq \%>\%
  samdat_tbl() \%>\%
  dplyr::select(1, subject, timepoint, dplyr::last_col())

# ggplot heatmap - unsorted
pseq \%>\%
  samdat_tbl() \%>\%
  filter(timepoint != 1) \%>\%
  ggplot(aes(x = timepoint, y = subject)) +
  geom_tile(aes(fill = aitchison_DistFromLast)) +
  scale_fill_viridis_c(na.value = NA, name = "dist") +
  theme_minimal(base_line_size = NA) +
  scale_y_discrete(limits = rev(levels(samdat_tbl(pseq)$subject)))

# ComplexHeatmap plotting with clustering #
library(tidyr)
library(ComplexHeatmap)

# make data matrix
heatmat <- pseq \%>\%
  samdat_tbl() \%>\%
  filter(timepoint != 1) \%>\%
  pivot_wider(
    id_cols = subject,
    names_from = timepoint, names_prefix = "t",
    values_from = aitchison_DistFromLast
  ) \%>\%
  tibble::column_to_rownames("subject")

heatmat <- as.matrix.data.frame(heatmat)

heatmap <- Heatmap(
  name = "dist",
  matrix = heatmat, col = viridisLite::viridis(12), na_col = "white",
  cluster_columns = FALSE,
  cluster_rows = hclust(dist(heatmat), method = "ward.D"),
  width = unit(1.5, "in"), rect_gp = gpar(col = "black"),
  row_names_side = "left", row_names_gp = gpar(fontsize = 8)
)
heatmap

# comparison with subject tracking on PCA
pca <- pseq \%>\%
  # already sorted data
  dist_calc("aitchison") \%>\%
  ord_calc("PCoA") \%>\%
  ord_plot(alpha = 0.1, shape = "nationality", size = 2) \%>\%
  add_paths(
    mapping = aes(colour = subject, alpha = timepoint, size = timepoint),
    id_var = "subject", id_values = c(
      "eve", "hsf", # low variation
      "vem", # medium
      "ufm", # high variation
      "pku" # starts high
    )
  ) +
  scale_alpha_continuous(range = c(0.3, 1), breaks = c(2, 4, 6)) +
  scale_size_continuous(range = c(1, 2), breaks = c(2, 4, 6))

heatmap
pca
}
\seealso{
\code{\link{dist_calc}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmap-utils.R
\name{adjacent_side}
\alias{adjacent_side}
\title{Simple heatmap helper to get a default adjacent side for another annotation}
\usage{
adjacent_side(side = c("top", "right", "bottom", "left"))
}
\arguments{
\item{side}{one of "top", "right", "bottom", or "left"}
}
\value{
character
}
\description{
Simple heatmap helper to get a default adjacent side for another annotation
}
\examples{
adjacent_side("top")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps_arrange.R
\name{ps_arrange}
\alias{ps_arrange}
\title{Arrange samples in phyloseq by sample_data variables or taxon abundance}
\usage{
ps_arrange(ps, ..., .target = "sample_data")
}
\arguments{
\item{ps}{phyloseq object}

\item{...}{dots passed directly to dplyr::arrange()}

\item{.target}{arrange samples by "sample_data" variables or "otu_table" taxa abundances}
}
\value{
phyloseq
}
\description{
Uses information in the sample_data or tax_table of phyloseq object
to set the order of the samples
(sample_data or tax_table specified by .target arg)

Give this function arguments in the same way you would use dplyr::arrange()
}
\examples{
data("dietswap", package = "microbiome")

dietswap \%>\%
  ps_arrange(subject, timepoint) \%>\%
  phyloseq::sample_data() \%>\%
  head(8)

ps <- dietswap \%>\% ps_arrange(subject, desc(timepoint))
phyloseq::sample_data(ps) \%>\% head(8)
phyloseq::otu_table(ps)[1:8, 1:8]

# you can also arrange samples by the abundances of taxa in the otu tables
pst <- dietswap \%>\% ps_arrange(desc(Akkermansia), .target = "otu_table")
phyloseq::otu_table(pst)[1:8, 1:8]
phyloseq::sample_data(pst) \%>\% head(8)
}
\seealso{
\code{\link[dplyr]{arrange}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{cor_heatmap}
\alias{cor_heatmap}
\title{Microbe-to-sample-data correlation heatmap}
\usage{
cor_heatmap(
  data,
  taxa = NA,
  tax_anno = taxAnnotation(Prev. = anno_tax_prev(), Abun. = anno_tax_box()),
  taxon_renamer = identity,
  vars = NA,
  var_anno = NULL,
  cor = c("pearson", "kendall", "spearman"),
  cor_use = "everything",
  colors = heat_palette(palette = "Blue-Red 2", sym = TRUE),
  numbers = heat_numbers(decimals = 1, col = "black", fontface = "plain"),
  taxa_side = "right",
  vars_side = adjacent_side(taxa_side),
  seriation_method = "OLO_ward",
  seriation_dist = "euclidean",
  seriation_method_col = seriation_method,
  seriation_dist_col = seriation_dist,
  var_fun = "identity",
  grid_col = "white",
  grid_lwd = 0.5,
  anno_tax = NULL,
  anno_vars = NULL,
  ...
)
}
\arguments{
\item{data}{phyloseq or phyloseq extra}

\item{taxa}{list of taxa to include, or NA for all}

\item{tax_anno}{NULL or annotation function for taxa: taxAnnotation() output.}

\item{taxon_renamer}{function to rename taxa before plotting}

\item{vars}{selection of variable names from sample_data}

\item{var_anno}{NULL or annotation function for variables: varAnnotation() output.}

\item{cor}{correlation coefficient. pearson/kendall/spearman,
can be abbreviated (used as legend title)}

\item{cor_use}{passed to cor(use = cor_use)}

\item{colors}{output of heat_palette() to set heatmap fill color scheme}

\item{numbers}{output of heat_numbers() to draw numbers on heatmap cells}

\item{taxa_side}{"top"/"right"/"bottom"/"left": controls heatmap orientation and where any
annotations specified in tax_anno are placed}

\item{vars_side}{which side to place any variable annotations specified in var_anno,
must be an adjacent side to taxa_side}

\item{seriation_method}{method to order the rows (in seriation::seriate)}

\item{seriation_dist}{distance to use in seriation_method (if needed)}

\item{seriation_method_col}{method to order the columns (in seriation::seriate)}

\item{seriation_dist_col}{distance to use in seriation_method_col (if needed)}

\item{var_fun}{a function (or name of) to be applied to columns of a matrix of vars
before correlating (but not used in any variable annotations)}

\item{grid_col}{colour of gridlines, or NA for none}

\item{grid_lwd}{width of gridlines}

\item{anno_tax}{DEPRECATED:
optional annotation of taxa distributions: tax_anno() list output,
or a pre-made ComplexHeatmap HeatmapAnnotation}

\item{anno_vars}{DEPRECATED: use var_anno argument instead.
Optional annotation of variable distributions:
var_anno() list output, or a pre-made ComplexHeatmap HeatmapAnnotation}

\item{...}{
  Arguments passed on to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap}}
  \describe{
    \item{\code{row_dend_width}}{Width of the row dendrogram, should be a \code{\link[grid]{unit}} object.}
    \item{\code{show_row_dend}}{Whether show row dendrogram?}
    \item{\code{row_dend_gp}}{Graphic parameters for the dendrogram segments. If users already provide a \code{\link[stats]{dendrogram}} object with edges rendered, this argument will be ignored.}
    \item{\code{show_row_names}}{Whether show row names.}
    \item{\code{row_names_gp}}{Graphic parameters for row names.}
    \item{\code{row_names_rot}}{Rotation of row names.}
    \item{\code{row_names_centered}}{Should row names put centered?}
    \item{\code{show_heatmap_legend}}{Whether show heatmap legend?}
  }}
}
\description{
Plot correlations between (transformed) microbial abundances and
(selected) numeric-like sample_data variables from a phyloseq object.

Lots of customisation options available through the listed arguments,
and you can pass any other argument from \code{ComplexHeatmap::Heatmap()} too.
}
\details{
Using a data.frame for the data argument is also possible,
in which case the (selected) numeric-like variables will be correlated
with each other, and all arguments relating to taxa will be ignored.
}
\examples{
library(dplyr)
data("dietswap", package = "microbiome")

# create a couple of numerical variables to use
psq <- dietswap \%>\%
  ps_mutate(
    weight = recode(bmi_group, obese = 3, overweight = 2, lean = 1),
    female = if_else(sex == "female", true = 1, false = 0),
    african = if_else(nationality == "AFR", true = 1, false = 0)
  )
psq <- tax_filter(psq, min_prevalence = 1 / 10, min_sample_abundance = 1 / 10)
psq <- tax_agg(psq, "Genus")

# randomly select 20 taxa from the 50 most abundant taxa
set.seed(123)
taxa <- sample(tax_top(psq, n = 50), size = 20)

# NOTE: detection threshold set to 50 as HITchip example data seems to have background noise
ud <- 50

# make simple correlation heatmap with all numeric-like variables
cor_heatmap(
  data = psq, taxa = taxa,
  tax_anno = taxAnnotation(
    Prv. = anno_tax_prev(undetected = ud),
    Abd. = anno_tax_box(undetected = ud)
  )
)

# You can create an annotation object separately in advance
taxAnno <- taxAnnotation(
  Prv. = anno_tax_prev(undetected = ud), Abd. = anno_tax_box(undetected = ud)
)
class(taxAnno) # "function"

# You can select which numeric-like variables to correlate taxa with
cor_heatmap(
  psq, taxa,
  vars = c("african", "female", "weight"), tax_anno = taxAnno
)

# Also you can choose alternative correlation measures
cor_heatmap(psq, taxa, cor = "spearman", tax_anno = taxAnno)

# Annotating variables is possible, and easy with varAnnotation()
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno,
  var_anno = varAnnotation(Val. = anno_var_box(size = grid::unit(2, "cm")))
)

# you can transform the variables before correlating by var_fun
# notice this does not affect the data used for annotations
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno, var_fun = "exp",
  var_anno = varAnnotation(Val. = anno_var_box(size = grid::unit(2, "cm")))
)

# other and multiple annotations
cor_heatmap(
  data = psq, taxa = taxa[1:10], vars = c("african", "weight", "female"),
  tax_anno = taxAnno,
  var_anno = varAnnotation(
    value = anno_var_hist(size = grid::unit(15, "mm")),
    log10p = anno_var_box(function(x) log10(x + 1))
  )
)

# make the same heatmap, but rotated
cor_heatmap(
  data = psq, taxa = taxa[1:10], vars = c("african", "weight", "female"),
  tax_anno = taxAnno, taxa_side = "top",
  var_anno = varAnnotation(
    value = anno_var_hist(size = grid::unit(15, "mm")),
    log10p = anno_var_box(function(x) log10(x + 1))
  )
)

# You can change the colour scheme used, using heat_palette()
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno,
  colors = heat_palette("Green-Orange", rev = TRUE, sym = TRUE)
)

# You can hide or change the style of the numbers with heat_numbers()
cor_heatmap(data = psq, taxa = taxa, tax_anno = taxAnno, numbers = NULL)
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno,
  colors = heat_palette("Berlin", rev = TRUE, sym = TRUE),
  numbers = heat_numbers(decimals = 2, col = "white", fontface = "bold")
)

# You can hide or change the style of the grid lines with grid_col & grid_lwd
cor_heatmap(psq, taxa = taxa, tax_anno = taxAnno, grid_col = NA) # hidden
cor_heatmap(psq, taxa = taxa, tax_anno = taxAnno, grid_lwd = 3) # bigger

# You can pass any other argument from `ComplexHeatmap::Heatmap()` to `...`

# e.g. You can set the absolute width and height of the heatmap body
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno,
  width = grid::unit(40, "mm"), height = grid::unit(10, "cm")
)

# e.g. You can suppress the legend
cor_heatmap(
  data = psq, taxa = taxa, tax_anno = taxAnno, show_heatmap_legend = FALSE,
  width = grid::unit(40, "mm"), height = grid::unit(10, "cm")
)
}
\seealso{
\code{\link[=taxAnnotation]{taxAnnotation()}} \code{\link[=varAnnotation]{varAnnotation()}}

\code{\link[=comp_heatmap]{comp_heatmap()}}

\code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_sort.R
\name{tax_sort}
\alias{tax_sort}
\title{Sort taxa in phyloseq otu_table and tax_table}
\usage{
tax_sort(
  data,
  by = "name",
  at = "names",
  ...,
  tree_warn = TRUE,
  verbose = TRUE,
  use_counts = TRUE
)
}
\arguments{
\item{data}{ps_extra or phyloseq}

\item{by}{how to sort, see description}

\item{at}{"names" or a taxonomic rank to apply sorting method to, as specified in \code{by}.}

\item{...}{used if summary function given, or pass \code{undetected} arg for tax_transform("binary") if by = "prev" or "prevalence"}

\item{tree_warn}{If phylogenetic tree is present in phyloseq phy_tree slot, taxa cannot be reordered.
Default behaviour of tax_sort is to remove the phylogenetic tree and warn about this.
tree_warn = FALSE will suppress the warning message, but still remove the tree!}

\item{verbose}{passed to phyloseq_validate verbose
(if TRUE: message about suspicious values in tax_table, and how to fix)}

\item{use_counts}{use count data if available, instead of transformed data}
}
\value{
sorted phyloseq or ps_extra
}
\description{
Multiple ways of sorting taxa are possible and determined by the \code{by} argument.
The \code{by} argument must be one of:
\itemize{
\item 'rev' to reverse the current order
\item 'name' (sort alphabetically by \code{at})
\item a sample name (descending abundance sorting within that sample)
\item summary stat. function e.g. \code{sum} or \code{mean}
}

The \code{at} argument must be "names" for sorting unique taxa,
or a rank name, for sorting at that rank.
\code{at} is ignored when \code{by} is "rev".
}
\details{
Don't forget to pass \code{na.rm = TRUE} to \code{...}
if using a summary stat function in \code{by}
}
\examples{
library(phyloseq)
data("dietswap", package = "microbiome")
dietswap

# reverse current order
dietswap \%>\%
  tax_sort("rev") \%>\%
  tax_table() \%>\%
  head(30)

# sort alphabetically by a taxonomic rank (or "names" for taxa_names)
dietswap \%>\%
  tax_sort(by = "name", at = "Phylum") \%>\%
  tax_table() \%>\%
  head(30)

# sequentially sorting by higher ranks
# sets tax_table in nested alphabetical order
dietswap \%>\%
  tax_sort(at = "names") \%>\%
  tax_sort(at = "Genus") \%>\%
  tax_sort(at = "Family") \%>\%
  tax_sort(at = "Phylum") \%>\%
  tax_table() \%>\%
  head(30)

# sort by function e.g. median abundance
dietswap \%>\%
  tax_sort(by = median) \%>\%
  taxa_names() \%>\%
  head(20)

# order by descending abundance in a single named sample
dietswap \%>\%
  tax_sort(by = "Sample-1") \%>\%
  otu_table() \%>\%
  .[1:8, 1:4]


# sum order should always equal mean order if non-negative abundances
# don't forget to add na.rm = TRUE if you expect NAs in otu_table somehow
dietswap \%>\%
  tax_sort(by = sum, na.rm = TRUE) \%>\%
  taxa_names() \%>\%
  head(20)

# if your phyloseq object has a phylogenetic tree,
# tax_sort will remove the tree, and warn you about this.

# You can sort by abundance at higher taxonomic ranks,
# without losing lower rank info
# e.g. sort (descending) by phyla abundances
dietswap \%>\%
  tax_sort(by = sum, at = "Phylum") \%>\%
  tax_table() \%>\%
  head()

# You can sort by ascending abundance by reversing afterwards
dietswap \%>\%
  tax_sort(by = "prev", at = "Phylum") \%>\%
  tax_sort(by = "rev") \%>\%
  tax_table() \%>\%
  head()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_chull.R
\name{stat_chull}
\alias{stat_chull}
\title{Draw convex hull for a set of points on a ggplot}
\usage{
stat_chull(
  mapping = NULL,
  data = NULL,
  geom = "polygonHollow",
  position = "identity",
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE,
  ...
)
}
\arguments{
\item{mapping}{Set of aesthetic mappings created by \code{\link[ggplot2:aes]{aes()}} or
\code{\link[ggplot2:aes_]{aes_()}}. If specified and \code{inherit.aes = TRUE} (the
default), it is combined with the default mapping at the top level of the
plot. You must supply \code{mapping} if there is no plot mapping.}

\item{data}{The data to be displayed in this layer. There are three
options:

If \code{NULL}, the default, the data is inherited from the plot
data as specified in the call to \code{\link[ggplot2:ggplot]{ggplot()}}.

A \code{data.frame}, or other object, will override the plot
data. All objects will be fortified to produce a data frame. See
\code{\link[ggplot2:fortify]{fortify()}} for which variables will be created.

A \code{function} will be called with a single argument,
the plot data. The return value must be a \code{data.frame}, and
will be used as the layer data. A \code{function} can be created
from a \code{formula} (e.g. \code{~ head(.x, 10)}).}

\item{geom}{The geometric object to use display the data}

\item{position}{Position adjustment, either as a string, or the result of
a call to a position adjustment function.}

\item{na.rm}{If \code{FALSE}, the default, missing values are removed with
a warning. If \code{TRUE}, missing values are silently removed.}

\item{show.legend}{logical. Should this layer be included in the legends?
\code{NA}, the default, includes if any aesthetics are mapped.
\code{FALSE} never includes, and \code{TRUE} always includes.
It can also be a named logical vector to finely select the aesthetics to
display.}

\item{inherit.aes}{If \code{FALSE}, overrides the default aesthetics,
rather than combining with them. This is most useful for helper functions
that define both data and aesthetics and shouldn't inherit behaviour from
the default plot specification, e.g. \code{\link[ggplot2:borders]{borders()}}.}

\item{...}{Other arguments passed on to \code{\link[ggplot2:layer]{layer()}}. These are
often aesthetics, used to set an aesthetic to a fixed value, like
\code{colour = "red"} or \code{size = 3}. They may also be parameters
to the paired geom/stat.}
}
\description{
Draws a (convex) polygon around the outermost points of a set of points.
Useful as a visual aid for identifying groups of points on a scatterplot,
such as an ordination plot.
}
\details{
This is a ggplot2 extension - slightly modified from the original code found here:

\url{https://CRAN.r-project.org/package=ggplot2/vignettes/extending-ggplot2.html}
}
\examples{
library(ggplot2)
corncob::ibd_phylo \%>\%
  tax_fix() \%>\%
  tax_transform(rank = "Genus", trans = "clr") \%>\%
  ord_calc(method = "PCA") \%>\%
  ord_plot(colour = "DiseaseState", shape = "DiseaseState", alpha = 0.5) +
  stat_chull(aes(colour = DiseaseState))

corncob::ibd_phylo \%>\%
  tax_fix() \%>\%
  tax_transform(rank = "Genus", trans = "clr") \%>\%
  ord_calc(method = "PCA") \%>\%
  ord_plot(colour = "DiseaseState", shape = "DiseaseState", alpha = 0.5) +
  stat_chull(aes(colour = DiseaseState, fill = DiseaseState), alpha = 0.1)
}
\seealso{
\code{ggplot2::\link{stat_ellipse}}

\code{\link{ord_plot}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heatmaps.R
\name{heat_palette}
\alias{heat_palette}
\title{Easy palettes for ComplexHeatmap}
\usage{
heat_palette(
  palette = ifelse(sym, "Blue-Red 3", "Rocket"),
  breaks = "auto",
  range = NA,
  sym = FALSE,
  rev = FALSE
)
}
\arguments{
\item{palette}{named palette from colorspace::hcl_palettes() diverging/sequential
or a vector of colour names/hexcodes}

\item{breaks}{number of breaks, "auto" is 11 for a named palette, or uses palette length}

\item{range}{NA to return palette generating function that takes range
or numeric vector indicating the range, to return a palette}

\item{sym}{makes palette range symmetrical around 0 if TRUE}

\item{rev}{reverse the palette?}
}
\value{
circlize::colorRamp2 palette if range = NA,
or function returning a palette when given a range
}
\description{
Pass a named colorspace hcl palette to circlize::colorRamp2.
\itemize{
\item If you do not specify a range this function returns a function and
the heatmap color palette will use the range of the data automatically
\item If you do specify a range, this returns a colour palette with that range
}
}
