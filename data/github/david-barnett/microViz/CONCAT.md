
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
