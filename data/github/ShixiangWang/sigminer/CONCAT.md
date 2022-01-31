
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Sigminer: Mutational Signature Analysis and Visualization in R <img src="man/figures/logo.png" alt="logo" align="right" height="140" width="120"/>

[![CRAN
status](https://www.r-pkg.org/badges/version/sigminer)](https://cran.r-project.org/package=sigminer)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/ShixiangWang/sigminer/workflows/R-CMD-check/badge.svg)](https://github.com/ShixiangWang/sigminer/actions)
[![Coverage
status](https://codecov.io/gh/ShixiangWang/sigminer/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/sigminer?branch=master)
[![](https://cranlogs.r-pkg.org/badges/grand-total/sigminer?color=orange)](https://cran.r-project.org/package=sigminer)
[![Closed
issues](https://img.shields.io/github/issues-closed/ShixiangWang/sigminer.svg)](https://github.com/ShixiangWang/sigminer/issues?q=is%3Aissue+is%3Aclosed)
[![Lines Of
Code](https://tokei.rs/b1/github/ShixiangWang/sigminer?category=code)](https://github.com/ShixiangWang/sigminer)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FShixiangWang%2Fsigminer&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
![install with
bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)
[![check in Biotreasury](https://img.shields.io/badge/Biotreasury-collected-brightgreen)](https://biotreasury.rjmart.cn/#/tool?id=10043) 

## :bar\_chart: Overview

The cancer genome is shaped by various mutational processes over its
lifetime, stemming from exogenous and cell-intrinsic DNA damage, and
error-prone DNA replication, leaving behind characteristic mutational
spectra, termed **mutational signatures**. This package, **sigminer**,
helps users to extract, analyze and visualize signatures from genome
alteration records, thus providing new insight into cancer study.

For pipeline tool, please see its co-evolutionary CLI
[sigflow](https://github.com/ShixiangWang/sigflow).

**SBS signatures**:

<img src="man/figures/README-unnamed-chunk-1-1.png" width="100%" />

**Copy number signatures**:

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

**DBS signatures**:

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

**INDEL (i.e. ID) signatures**:

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

**Genome rearrangement signatures**:

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

### :airplane: Features

-   supports a standard *de novo* pipeline for identification of **5**
    types of signatures: copy number, SBS, DBS, INDEL and RS (genome
    rearrangement signature).
-   supports quantify exposure for one sample based on *known
    signatures*.
-   supports association and group analysis and visualization for
    signatures.
-   supports two types of signature exposures: relative exposure
    (relative contribution of signatures in each sample) and absolute
    exposure (estimated variation records of signatures in each sample).
-   supports basic summary and visualization for profile of mutation
    (powered by **maftools**) and copy number.
-   supports parallel computation by R packages **foreach**, **future**
    and **NMF**.
-   efficient code powered by R packages **data.table** and
    **tidyverse**.
-   elegant plots powered by R packages **ggplot2**, **ggpubr**,
    **cowplot** and **patchwork**.
-   well tested by R package **testthat** and documented by R package
    **roxygen2**, **roxytest**, **pkgdown**, and etc. for both reliable
    and reproducible research.

## :arrow\_double\_down: Installation

You can install the stable release of **sigminer** from CRAN with:

``` r
install.packages("BiocManager")
BiocManager::install("sigminer", dependencies = TRUE)
```

You can install the development version of **sigminer** from Github
with:

``` r
remotes::install_github("ShixiangWang/sigminer", dependencies = TRUE)
# For Chinese users, run 
remotes::install_git("https://gitee.com/ShixiangWang/sigminer", dependencies = TRUE)
```

You can also install **sigminer** from conda `bioconda` channel with

``` sh
# Please note version number of the bioconda release

# You can install an individual environment firstly with
# conda create -n sigminer
# conda activate sigminer
conda install -c bioconda -c conda-forge r-sigminer
```

## :beginner: Usage

A complete documentation of **sigminer** can be read online at
<https://shixiangwang.github.io/sigminer-book/>. All functions are well
organized and documented at
<https://shixiangwang.github.io/sigminer/reference/index.html>. For
usage of a specific function `fun`, run `?fun` in your R console to see
its documentation.

## :paperclip: Citation

If you use **sigminer** in academic field, please cite one of the
following papers.

------------------------------------------------------------------------

-   ***Wang S, Li H, Song M, Tao Z, Wu T, He Z, et al. (2021) Copy
    number signature analysis tool and its application in prostate
    cancer reveals distinct mutational processes and clinical outcomes.
    PLoS Genet 17(5): e1009557.***
    <https://doi.org/10.1371/journal.pgen.1009557>

-   ***Shixiang Wang, Ziyu Tao, Tao Wu, Xue-Song Liu, Sigflow: An
    Automated And Comprehensive Pipeline For Cancer Genome Mutational
    Signature Analysis, Bioinformatics, btaa895***.
    <https://doi.org/10.1093/bioinformatics/btaa895>

------------------------------------------------------------------------

## :arrow\_down: Download Stats

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

## :page\_with\_curl: References

Please properly cite the following references when you are using any
corresponding features. The references are also listed in the function
documentation. Very thanks to the works, **sigminer** cannot be created
without the giants.

1.  Mayakonda, Anand, et al. “Maftools: efficient and comprehensive
    analysis of somatic variants in cancer.” Genome research 28.11
    (2018): 1747-1756.
2.  Gaujoux, Renaud, and Cathal Seoighe. “A Flexible R Package for
    Nonnegative Matrix Factorization.”" BMC Bioinformatics 11, no. 1
    (December 2010).
3.  H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
    Springer-Verlag New York, 2016.
4.  Kim, Jaegil, et al. “Somatic ERCC2 mutations are associated with a
    distinct genomic signature in urothelial tumors.” Nature genetics
    48.6 (2016): 600.
5.  Alexandrov, Ludmil B., et al. “Deciphering signatures of mutational
    processes operative in human cancer.” Cell reports 3.1 (2013):
    246-259.
6.  Degasperi, Andrea, et al. “A practical framework and online tool for
    mutational signature analyses show intertissue variation and driver
    dependencies.” Nature cancer 1.2 (2020): 249-263.
7.  Alexandrov, Ludmil B., et al. “The repertoire of mutational
    signatures in human cancer.” Nature 578.7793 (2020): 94-101.
8.  Macintyre, Geoff, et al. “Copy number signatures and mutational
    processes in ovarian carcinoma.” Nature genetics 50.9 (2018): 1262.
9.  Tan, Vincent YF, and Cédric Févotte. “Automatic relevance
    determination in nonnegative matrix factorization with the/spl
    beta/-divergence.” IEEE Transactions on Pattern Analysis and Machine
    Intelligence 35.7 (2012): 1592-1605.
10. Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG,
    Alexandrov LB: SigProfilerMatrixGenerator: a tool for visualizing
    and exploring patterns of small mutational events. BMC Genomics
    2019, 20:685
    <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2>

## :page\_facing\_up: LICENSE

The software is made available for non commercial research purposes only
under the
[MIT](https://github.com/ShixiangWang/sigminer/blob/master/LICENSE.md).
However, notwithstanding any provision of the MIT License, the software
currently may not be used for commercial purposes without explicit
written permission after contacting Shixiang Wang
<wangshx@shanghaitech.edu.cn> or Xue-Song Liu
<liuxs@shanghaitech.edu.cn>.

MIT © 2019-Present Shixiang Wang, Xue-Song Liu

MIT © 2018 Anand Mayakonda

------------------------------------------------------------------------

[**Cancer Biology Group**](https://github.com/XSLiuLab)
**@ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**

![Alt](https://repobeats.axiom.co/api/embed/7cd2cf8a196dde9d8d1e13c9b23bc2f157d8254e.svg "Repobeats analytics image")
# sigminer 2.1.2

- Enhanced the `read_copynumber_seqz()` to include minor copy number. (Thanks to yancey)
- Added input `range` check in `sig_estimate()`. (#391)

# sigminer 2.1.1

- Expanded `output_*` function by adding option `sig_db`.
- Fixed the error using `sigminer::get_genome_annotation()` before loading it.
- Fixed the bug the `get_pLOH_score()` return nothing for sample without LOH.

# sigminer 2.1.0

- Added `sig_unify_extract()` as an unified signature extractor.
- Fixed error showing reference signature profile for `CNS_TCGA` database.

# sigminer 2.0.5

- Impl `y_limits` option in `show_sig_profile()` (#381).
- Added function `get_pLOH_score()` for representing the genome that displayed LOH.
- Added function `read_copynumber_ascat()` for reading ASCAT result ASCAT object in `.rds` format.
- Added function `get_intersect_size()` for getting overlap size between intervals.
- Added option to `get_Aneuploidy_score()` to remove short arms of chr13/14/15/21/22 from calculation.

# sigminer 2.0.4

- Implemented Cohen-Sharir method-like Aneuploidy Score.
- Enhanced error handling in `show_sig_feature_corrplot()` (#376).
- Fixed INDEL classification.
- Fixed end position determination in `read_vcf()`.
- Updated INDEL adjustment.
- Included TCGA copy number signatures from SigProfiler.
- Updated docs.

# sigminer 2.0.3

- Preprocessed INDELs before labeling them in `sig_tally()` (#370).
- Fixed `sigprofiler_extract()` extracting copy number signatures and
rolled up sigprofiler version (#369).

# sigminer 2.0.2

- Fixed `output_sig()` error in handling exposure plot with >9 signatures (#366).
- Added `limitsize = FALSE` for `ggsave()` or `ggsave2()` for handling big figure.

# sigminer 2.0.1

- Supported `mm9` genome build.
- Removed FTP link as CRAN suggested (#359).
- Updated README.

# sigminer 2.0.0

## BUG REPORTS

- Fixed the SigProfiler installation error due to Python version in conda environment.
- Fixed classification bug due to repeated function name `call_component`.
- Fixed the bug when `read_vcf()` with `##` commented VCF files.

## ENHANCEMENTS

- Added support for latest COSMIC v3.2 as reference signatures. You can obtain them by

```r
for (i in c("latest_SBS_GRCh37", "latest_DBS_GRCh37", "latest_ID_GRCh37",
            "latest_SBS_GRCh38", "latest_DBS_GRCh38",
            "latest_SBS_mm9", "latest_DBS_mm9",
            "latest_SBS_mm10", "latest_DBS_mm10",
            "latest_SBS_rn6", "latest_DBS_rn6")) {
  message(i)
  get_sig_db(i)
}
```

- Updated `keep_only_pass` to `FALSE` at default.
- Added RSS and unexplained variance calculation in `get_sig_rec_similarity()`.
- Added data check and filter in `output_tally()` and `show_catalogue()`.
- Enhanced `show_group_enrichment()` (#353) & added a new option to cluster rows.
- Removed unnecessary CN classifications code in recent development.

## NEW FUNCTIONS

## DEPRECATED

- Dropped copy number "M"" method to avoid misguiding user to use/read wrong signature profile and keep code simple.

# sigminer 1.2.5

## BUG REPORTS

## ENHANCEMENTS

- Modified the default visualization of `bp_show_survey()`.
- Enhanced `torch` check.

## NEW FUNCTIONS

- `read_sv_as_rs()` and `sig_tally.RS()` for simplified genome rearrangement classification matrix generation (experimental).

## DEPRECATED

# sigminer 1.2.4

## BUG REPORTS

- Fixed the assign problem about match pair in `bp_extract_signatures()` 
with `lpSolve` package instead of using my problematic code.

## ENHANCEMENTS

- Supported `mm10` in `read_vcf()`.
- Removed large data files and store them in Zenodo to reduce package size.
- Added cores check.
- Upgraded SP to v1.1.0 (need test).
- Tried installing Torch before SP (need test).

## NEW FUNCTIONS

## DEPRECATED

# sigminer 1.2.2

## BUG REPORTS

- Fixed bug in silhouette calculation in `bp_extract_signatures()` (#332).
PAY ATTENTION: this may affect results.
- Fixed bug using custom signature name in `show_sig_profile_loop()`.

## ENHANCEMENTS

- Subset signatures to plot is available by `sig_names` option.
- sigminer is available in bioconda channel: `https://anaconda.org/bioconda/r-sigminer/`
- Updated `ms` strategy in `sig_auto_extract()` by assigning each signature to its
best matched reference signatures.
- Added `get_shannon_diversity_index()` to get diversity index for signatures (#333).
- Added new method "S" (from Steele et al. 2019) for tallying copy number data (#329).
- Included new (RS) reference signatures (related to #331).
- Updated the internal code for getting relative activity in `get_sig_exposure()`.

## NEW FUNCTIONS

- `bp_get_clustered_sigs()` to get clustered mean signatures.

## DEPRECATED

# sigminer 1.2.1

- Updated author list.

## BUG REPORTS

## ENHANCEMENTS

- Added a quick start vignette.
- A new option `highlight` is added to `show_sig_number_survey()` and `bp_show_survey2()` to highlight a selected number.

## NEW FUNCTIONS

## DEPRECATED

# sigminer 1.2.0

## BUG REPORTS

## ENHANCEMENTS

- A new option `cut_p_value` is added to `show_group_enrichment()` to cut continous p values as binned regions.
- A Python backend for `sig_extract()` is provided.
- User now can directly use `sig_extract()` and `sig_auto_extract()` instead of loading NMF package firstly.
- Added benchmark results for different extraction approaches in README.
- The threshold for `auto_reduce` in `sig_fit()` is modified from 0.99 to 0.95 and similarity update threshold updated from `>0` to `>=0.01`.
- Removed `pConstant` option from `sig_extract()` and `sig_estimate()`. Now a
auto-check function is created for avoiding the error from NMF package due to
no contribution of a component in all samples.

## NEW FUNCTIONS

- `bp_show_survey2()` to plot a simplified version for signature number survey (#330).
- `read_xena_variants()` to read variant data from UCSC Xena as a `MAF` object for signature analysis.
- `get_sig_rec_similarity()` for getting reconstructed profile similarity for `Signature` object (#293).
- Added functions start with `bp_` which are combined to provide a best practice for extracting
signatures in cancer researches. See more details, run `?bp` in your R console.

## DEPRECATED

# sigminer 1.1.0

- Added data simulation.
- Suppressed `future` warnings.
- Fixed p value calculation in bootstrap analysis.
- Fixed typo in `show_cor()`, thanks to @Miachol.
- Added `y_tr` option in `show_sig_profile()` to transform y axis values.
- Optimized default behavior of `read_copynumber()`.
  - Support LOH records when user input minor allele copy number.
  - Set `complement = FALSE` as default.
  - Free dependencies between option `use_all` and `complement`.
- Added visualization support for genome rearrangement signatures (#300).
- Added four database for reference signatures from <https://doi.org/10.1038/s43018-020-0027-5> (#299).
- Added new measure 'CV' for `show_sig_bootstrap()` (#298).
- Added `group_enrichment()` and `show_group_enrichment()` (#277).
- Optimized signature profile visualization (#295).
- Updated `?sigminer` documentation.
- Added `ms` strategy to select optimal solution by maximizing cosine similarity
to reference signatures.
- Added `same_size_clustering()` for same size clustering.
- Added `show_cosmic()` to support reading COSMIC signatures in web browser (#288).
- Changed argument `rel_threshold` behavior in `sig_fit()` and `get_sig_exposure()`.
Made them more consistent and allowed un-assigned signature contribution (#285).
- Updated all COSMIC signatures to v3.1 and their aetiologies (#287).

# sigminer 1.0.19

- Added more specific reference signatures from SigProfiler, e.g. `SBS_mm9`.
- Supported `data.frame` as input object for `sig` in `get_sig_similarity()` and `sig_fit()`.
- Modified `g_label` option in `show_group_distribution()` to better control group names. 
- Added `test` option and variable checking in `show_cor()`.
- Updated `output_sig()` to output signature exposure distribution ([#280](https://github.com/ShixiangWang/sigminer/issues/280)).
- Added `show_cor()` for general association analysis.
- Added options in `show_group_distribution()` to control segments.

# sigminer 1.0.18

- Fixed bugs when outputing only 1 signatures.
- Fixed label inverse bug in `add_labels()`, thanks to TaoTao for reporting.

# sigminer 1.0.17

- Handled `,` seperated indices in show_cosmic_signatures.
- Added option `set_order` in `get_sig_similarity()` (#274).
- Outputed more stats information in `output_sig()`.
- Fixed default y axis title in `show_sig_bootstrap_error()`, now it is "Reconstruction error (L2 norm)"

# sigminer 1.0.16

- Added `auto_reduce` option in `sig_fit*` functions to improve signature fitting.
- Return cosine similarity for sample profile in `sig_fit()`.
- Set default strategy in `sig_auto_extract()` to 'optimal'.
- Supported search reference signature index in `get_sig_cancer_type_index()`.
- Outputed legacy COSMIC similarity for SBS signatures.
- Added new option in `sigprofiler_extract()` to reduce failure in when `refit` is enabled.
- Outputed both relative and absolute signature exposure in `output_sig()`.
- Updated background color in `show_group_distribution()`.
- Modified the default theme for signature profile in COSMIC style.
- Updated the copy number classification method.

# sigminer 1.0.15

- Handled null catalogue.
- Supported ordering the signatures for results from SigProfiler.
- Supported importing refit results from SigProfiler.
- Set `optimize` option in `sig_extract()` and `sig_auto_extract()`.

# sigminer 1.0.14

- Supported signature index separated by `,` in `sig_fit()` and `sig_fit_bootstrap*` functions.
- Added `output_*` functions from [sigflow](https://github.com/ShixiangWang/sigflow).
- Enhanced DBS search and error handling in `sig_tally()`.
- Added option `highlight_genes` in `show_cn_group_profile()` to show gene labels.
- Added `get_sig_cancer_type_index()` to get reference signature index.
- Added `show_group_distribution()` to show group distribution.
- Added options in `show_cn_profile()` to show specified ranges and add copy number value labels.
- Used package `nnls` instead of `pracma` for NNLS implementation in `sig_fit()`.

# sigminer 1.0.13

- Supported `BSgenome.Hsapiens.1000genomes.hs37d5` in `sig_tally()`.
- Remove changing `MT` to `M` in mutation data.
- Fixed bug in extract numeric signature names and signature orderings in `show_sig_exposure()`.
- Added `letter_colors` as an unexported discrete palette.

# sigminer 1.0.12

- Added `transform_seg_table()`.
- Added `show_cn_group_profile()`.
- Added `show_cn_freq_circos()`.
- `sig_orders` option in `show_sig_profile()` function now can select and order signatures to plot.
- Added `show_sig_profile_loop()` for better signature profile visualization.

# sigminer 1.0.11

- Added option to control the SigProfilerExtractor to avoid issue in docker image build.

# sigminer 1.0.10

- Some updates.
- Compatible with SigProfiler 1.0.15

# sigminer 1.0.9

- Tried to speed up joining adjacent segments in `read_copynumber()`, got 200% improvement.

# sigminer 1.0.8

- Tried to speed up joining adjacent segments in `read_copynumber()`, got 20% improvement.
- Added `cosine()` function.
- Added and exported `get_sig_db()` to let users directly load signature database.
- Added `sigprofiler_extract()` and `sigprofiler_import()` to call SigProfiler and import results.
- Added `read_vcf()` for simply reading VCF files.
- Implemented DBS-1248.
- Added `show_sig_profile_heatmap()`.
- Supported mouse genome 'mm10' ([#241](https://github.com/ShixiangWang/sigminer/issues/241)).
- Added `read_copynumber_seqz()` to read sequenza result directory.
- Speed up the annotation process in `read_copynumber()`.

# sigminer 1.0.7

- Fixed bug in OsCN feature calculation.
- Removed useless options in `read_maf()`.
- Modify method 'LS' in `sig_fit()` to 'NNLS' and implement it with **pracma** package ([#216](https://github.com/ShixiangWang/sigminer/issues/216)).
- Made `use_all` option in `read_copynumber()` working correctly.
- Fixed potential problem raised by unordered copy number segments ([#217](https://github.com/ShixiangWang/sigminer/issues/217)).
- Fixed a typo, correct `MRSE` to `RMSE`.
- Added feature in `show_sig_bootstrap_*()` for plotting aggregated values.
- Fixed bug when use `get_groups()` for clustering.
- Fixed bug about using reference components from NatGen 2018 paper.
- Added option `highlight_size` for `show_sig_bootstrap_*()`.
- Fixed bug about signature profile plotting for method 'M'.

# sigminer 1.0.6

- Added "scatter" in `sig_fit()` function to better visualize a few samples.
- Added "highlight" option.
- `lsei` package was removed from CRAN, here I reset default method to 'QP' and tried best to keep the LS usage in sigminer ([#189](https://github.com/ShixiangWang/sigminer/issues/189)).
- Made consistent copy number labels in `show_sig_profile()` and added input checking for this function.
- Fixed unconsistent bootstrap when use `furrr`, solution is from <https://github.com/DavisVaughan/furrr/issues/107>.
- Properly handled null-count sample in `sig_fit()` for methods `QP` and `SA`.
- Supported boxplot or violin in `show_sig_fit()` and `show_sig_bootstrap_*` functions.
- Added job mode for `sig_fit_bootstrap_batch` for more useful in practice.
- Added `show_groups()` to show the signature contribution in each group from `get_groups()`.
- Expanded clustering in `get_groups()` to result of `sig_fit()`.
- Properly handled null-count samples in `sig_fit_bootstrap_batch()`.
- Added strand bias labeling for INDEL.
- Added COSMIC TSB signatures.

# sigminer 1.0.5

- Exported APOBEC result when the mode is 'ALL' in `sig_tally()`.
- Added batch bootstrap analysis feature (#158).
- Supported all common signature plotting.
- Added strand feature to signature profile.

# sigminer 1.0.4

- Added profile plot for DBS and INDEL.
- Fixed error for signature extraction in mode 'DBS' or 'ID'.
- Fixed method 'M' for CN tally cannot work when `cores > 1` (#161).

# sigminer 1.0.3

- Added multiple methods for `sig_fit()`.
- Added feature `sig_fit_bootstrap()` for bootstrap results.
- Added multiple classification method for SBS signature.
- Added strand bias enrichment analysis for SBS signature.
- Moved multiple packages from field `Imports` to `Suggests`.
- Added feature `report_bootstrap_p_value()` to report p values.
- Added common DBS and ID signature.
- Updated citation.

# sigminer 1.0.2

- Added merged transcript info for hg19 and hg38 build, this is availabe by `data()`.
- Added gene info for hg19 and hg38 build to extdata directory.

# sigminer 1.0.1

- Removed `fuzzyjoin` package from dependency.
- Moved `ggalluvial` package to field `suggsets`.

# sigminer 1.0.0

All users, this is a break-through version of **sigminer**,
most of functions have been modified, more features are implemented.
Please read the reference list to see the function groups and their
functionalities.

Please read the vignette for usage.

I Hope it helps your research work and makes a new contribution
to the scientific community.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  warning = FALSE
)
```

# Sigminer: Mutational Signature Analysis and Visualization in R <img src="man/figures/logo.png" alt="logo" align="right" height="140" width="120"/>

[![CRAN status](https://www.r-pkg.org/badges/version/sigminer)](https://cran.r-project.org/package=sigminer) [![lifecycle](https://img.shields.io/badge/lifecycle-stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) [![R-CMD-check](https://github.com/ShixiangWang/sigminer/workflows/R-CMD-check/badge.svg)](https://github.com/ShixiangWang/sigminer/actions) [![Coverage status](https://codecov.io/gh/ShixiangWang/sigminer/branch/master/graph/badge.svg)](https://codecov.io/github/ShixiangWang/sigminer?branch=master) [![](https://cranlogs.r-pkg.org/badges/grand-total/sigminer?color=orange)](https://cran.r-project.org/package=sigminer) [![Closed issues](https://img.shields.io/github/issues-closed/ShixiangWang/sigminer.svg)](https://github.com/ShixiangWang/sigminer/issues?q=is%3Aissue+is%3Aclosed) [![Lines Of Code](https://tokei.rs/b1/github/ShixiangWang/sigminer?category=code)](https://github.com/ShixiangWang/sigminer) [![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2FShixiangWang%2Fsigminer&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com) ![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square) [![check in Biotreasury](https://img.shields.io/badge/Biotreasury-collected-brightgreen)](https://biotreasury.rjmart.cn/#/tool?id=10043) 


## :bar_chart: Overview

The cancer genome is shaped by various mutational processes over its lifetime, stemming from exogenous and cell-intrinsic DNA damage, and error-prone DNA replication, leaving behind characteristic mutational spectra, termed **mutational signatures**. This package, **sigminer**, helps users to extract, analyze and visualize signatures from genome alteration records, thus providing new insight into cancer study.

For pipeline tool, please see its co-evolutionary CLI [sigflow](https://github.com/ShixiangWang/sigflow).

**SBS signatures**:

```{r, fig.width=12, fig.height=5, echo=FALSE, message=FALSE}
library(sigminer)
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p1 <- show_sig_profile(sig2, mode = "SBS", style = "cosmic", x_label_angle = 90)
p1

```

**Copy number signatures**:

```{r fig.width=12, fig.height=5, echo=FALSE}
load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p2 <- show_sig_profile(sig,
  style = "cosmic",
  mode = "copynumber",
  method = "W",
  normalize = "feature"
)
p2
```

**DBS signatures**:

```{r fig.width=12, fig.height=5, echo=FALSE}
DBS = system.file("extdata", "DBS_signatures.rds",
                        package = "sigminer", mustWork = TRUE)
DBS = readRDS(DBS)

# Show signature profile
p3 <- show_sig_profile(DBS$db[, 1:3] %>% as.matrix(), mode = "DBS", style = "cosmic", check_sig_names = F)
p3
```

**INDEL (i.e. ID) signatures**:

```{r fig.width=12, fig.height=5, echo=FALSE}
ID = system.file("extdata", "ID_signatures.rds",
                        package = "sigminer", mustWork = TRUE)
ID = readRDS(ID)

# Show signature profile
p4 <- show_sig_profile(ID$db[, 4:6] %>% as.matrix(), mode = "ID", style = "cosmic", check_sig_names = F)
p4
```

**Genome rearrangement signatures**:

```{r fig.width=10, fig.height=5, echo=FALSE, message=FALSE}
p5 <- show_cosmic_sig_profile(sig_index = c("R1", "R2", "R3"), sig_db = "RS_Nik_lab", style = "cosmic", show_index = FALSE)
p5
```


### :airplane: Features

-   supports a standard *de novo* pipeline for identification of **5** types of signatures: copy number, SBS, DBS, INDEL and RS (genome rearrangement signature).
-   supports quantify exposure for one sample based on *known signatures*.
-   supports association and group analysis and visualization for signatures.
-   supports two types of signature exposures: relative exposure (relative contribution of signatures in each sample) and absolute exposure (estimated variation records of signatures in each sample).
-   supports basic summary and visualization for profile of mutation (powered by **maftools**) and copy number.
-   supports parallel computation by R packages **foreach**, **future** and **NMF**.
-   efficient code powered by R packages **data.table** and **tidyverse**.
-   elegant plots powered by R packages **ggplot2**, **ggpubr**, **cowplot** and **patchwork**.
-   well tested by R package **testthat** and documented by R package **roxygen2**, **roxytest**, **pkgdown**, and etc. for both reliable and reproducible research.


## :arrow_double_down: Installation

You can install the stable release of **sigminer** from CRAN with:

```{r, eval=FALSE}
install.packages("BiocManager")
BiocManager::install("sigminer", dependencies = TRUE)
```

You can install the development version of **sigminer** from Github with:

```{r, eval=FALSE}
remotes::install_github("ShixiangWang/sigminer", dependencies = TRUE)
# For Chinese users, run 
remotes::install_git("https://gitee.com/ShixiangWang/sigminer", dependencies = TRUE)
```

You can also install **sigminer** from conda `bioconda` channel with

```sh
# Please note version number of the bioconda release

# You can install an individual environment firstly with
# conda create -n sigminer
# conda activate sigminer
conda install -c bioconda -c conda-forge r-sigminer
```

## :beginner: Usage

A complete documentation of **sigminer** can be read online at <https://shixiangwang.github.io/sigminer-book/>. All functions are well organized and documented at <https://shixiangwang.github.io/sigminer/reference/index.html>. For usage of a specific function `fun`, run `?fun` in your R console to see its documentation.

## :paperclip: Citation

If you use **sigminer** in academic field, please cite one of the following papers.

------------------------------------------------------------------------

-   ***Wang S, Li H, Song M, Tao Z, Wu T, He Z, et al. (2021) Copy number signature analysis tool and its application in prostate cancer reveals distinct mutational processes and clinical outcomes. PLoS Genet 17(5): e1009557.*** <https://doi.org/10.1371/journal.pgen.1009557>

-   ***Shixiang Wang, Ziyu Tao, Tao Wu, Xue-Song Liu, Sigflow: An Automated And Comprehensive Pipeline For Cancer Genome Mutational Signature Analysis, Bioinformatics, btaa895***. <https://doi.org/10.1093/bioinformatics/btaa895>

------------------------------------------------------------------------

## :arrow_down: Download Stats

```{r, echo=FALSE, warning=FALSE, message=FALSE}
library("ggplot2")
library("cranlogs")
library("dplyr")

x <- cran_downloads("sigminer", from = "2019-07-01", to = "2021-09-01")

if (!is.null(x)) {
   #head(x)
   p <- ggplot(x, aes(date, cumsum(count), group=package, color=package)) +
       geom_line() + geom_point(aes(shape=package)) + 
    theme_minimal() + theme(legend.position = "none") +
    labs(x = "Time", y = "Cumulative downloads")
   p
}
```

## :page_with_curl: References

Please properly cite the following references when you are using any corresponding features. The references are also listed in the function documentation. Very thanks to the works, **sigminer** cannot be created without the giants.

1.  Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
2.  Gaujoux, Renaud, and Cathal Seoighe. "A Flexible R Package for Nonnegative Matrix Factorization."" BMC Bioinformatics 11, no. 1 (December 2010).
3.  H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
4.  Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors." Nature genetics 48.6 (2016): 600.
5.  Alexandrov, Ludmil B., et al. "Deciphering signatures of mutational processes operative in human cancer." Cell reports 3.1 (2013): 246-259.
6.  Degasperi, Andrea, et al. "A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies." Nature cancer 1.2 (2020): 249-263.
7.  Alexandrov, Ludmil B., et al. "The repertoire of mutational signatures in human cancer." Nature 578.7793 (2020): 94-101.
8.  Macintyre, Geoff, et al. "Copy number signatures and mutational processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.
9.  Tan, Vincent YF, and Cédric Févotte. "Automatic relevance determination in nonnegative matrix factorization with the/spl beta/-divergence." IEEE Transactions on Pattern Analysis and Machine Intelligence 35.7 (2012): 1592-1605.
10. Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB: SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics 2019, 20:685 <https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2>

## :page_facing_up: LICENSE

The software is made available for non commercial research purposes only under the [MIT](https://github.com/ShixiangWang/sigminer/blob/master/LICENSE.md). However, notwithstanding any provision of the MIT License, the software currently may not be used for commercial purposes without explicit written permission after contacting Shixiang Wang [wangshx\@shanghaitech.edu.cn](mailto:wangshx@shanghaitech.edu.cn) or Xue-Song Liu [liuxs\@shanghaitech.edu.cn](mailto:liuxs@shanghaitech.edu.cn).

MIT © 2019-Present Shixiang Wang, Xue-Song Liu

MIT © 2018 Anand Mayakonda

------------------------------------------------------------------------

[**Cancer Biology Group**](https://github.com/XSLiuLab) **\@ShanghaiTech**

**Research group led by Xue-Song Liu in ShanghaiTech University**

![Alt](https://repobeats.axiom.co/api/embed/7cd2cf8a196dde9d8d1e13c9b23bc2f157d8254e.svg "Repobeats analytics image")
---
title: "A Quick Start of sigminer Package"
author: "Shixiang Wang ( wangshx@shanghaitech.edu.cn )"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{A Quick Start of sigminer Package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, echo = FALSE, message = FALSE}
library(markdown)
library(knitr)
knitr::opts_chunk$set(
    error = FALSE,
    tidy  = TRUE,
    message = FALSE,
    fig.align = "center",
    collapse = TRUE,
    comment = "#>")
options(width = 100)
options(rmarkdown.html_vignette.check_title = FALSE)
```

Assume you have already gotten a catalog matrix (sample-by-component) like below:

```{r message=TRUE}
library(sigminer)
data("simulated_catalogs")
mat <- t(simulated_catalogs$set1)

mat[1:5, 1:5]
```

Extract signatures with:

```{r eval=FALSE}
# Here I reduce the values for n_bootstrap and n_nmf_run
# for reducing the run time.
# In practice, you should keep default or increase the values
# for better estimation.
#
# The input data here is simulated from 10 mutational signatures
e1 <- bp_extract_signatures(
  mat,
  range = 8:12,
  n_bootstrap = 5,
  n_nmf_run = 10
)
```

```{r include=FALSE}
e1 <- readRDS("e1.rds")
```

Check which signature number is proper:

```{r message=TRUE, fig.width=4, fig.height=3}
bp_show_survey2(e1, highlight = 10)
```

Get the `10` signatures:

```{r}
obj <- bp_get_sig_obj(e1, 10)
```

Show signature profile:

```{r fig.width=10, fig.height=8}
show_sig_profile(obj, mode = "SBS", style = "cosmic")
```
Show signature activity (a.k.a. exposure) profile:

```{r fig.width=8, fig.height=5}
show_sig_exposure(obj, rm_space = TRUE)
```

Calculate the similarity to COSMIC reference signatures:

```{r message=TRUE}
sim <- get_sig_similarity(obj, sig_db = "SBS")
```

```{r fig.width=10, fig.height=6}
if (require(pheatmap)) {
  pheatmap::pheatmap(sim$similarity)
}
```


## More

Please go to [*reference* list](https://shixiangwang.github.io/sigminer/reference/index.html) for well organized functions and documentation.

For more about mutational signature and **sigminer** usage, you can read [*sigminer-book*](https://shixiangwang.github.io/sigminer-book/).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{transcript.mm10}
\alias{transcript.mm10}
\title{Merged Transcript Location at Genome Build mm10}
\format{
A \code{data.table}
}
\source{
from GENCODE release M25.
}
\description{
Merged Transcript Location at Genome Build mm10
}
\examples{
data(transcript.mm10)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_pLOH_score.R
\name{get_pLOH_score}
\alias{get_pLOH_score}
\title{Get proportions of pLOH score from Allele Specific Copy Number Profile}
\usage{
get_pLOH_score(data, rm_chrs = c("chrX", "chrY"), genome_build = "hg19")
}
\arguments{
\item{data}{a CopyNumber object or a \code{data.frame} containing at least
'chromosome', 'start', 'end', 'segVal', "minor_cn", 'sample' these columns.}

\item{rm_chrs}{chromosomes to be removed in calculation. Default is sex
chromosomes (recommended).}

\item{genome_build}{genome build version, should be 'hg19', 'hg38', 'mm9' or 'mm10'.}
}
\value{
A \code{data.frame}
}
\description{
pLOH score represents the genome that displayed LOH.
}
\examples{
# Load toy dataset of absolute copynumber profile
load(system.file("extdata", "toy_segTab.RData",
  package = "sigminer", mustWork = TRUE
))

set.seed(1234)
segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
cn <- read_copynumber(segTabs,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  genome_measure = "wg", complement = TRUE, add_loh = TRUE
)

df <- get_pLOH_score(cn)
df

df2 <- get_pLOH_score(cn@data)
df2

}
\references{
Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." bioRxiv (2021).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cosmic_signature.R
\name{show_cosmic}
\alias{show_cosmic}
\title{Show Signature Information in Web Browser}
\usage{
show_cosmic(x = "home")
}
\arguments{
\item{x}{a string indicating location
("home" for COSMIC signature home, "legacy" for COSMIC v2 signatures,
"SBS" for COSMIC v3 SBS signatures, "DBS" for COSMIC v3 DBS signatures,
"ID" for COSMIC v3 INDEL signatures) or signature index (e.g.
"SBS1", "DBS2", "ID3").}
}
\value{
Nothing.
}
\description{
Show Signature Information in Web Browser
}
\examples{
\dontrun{
show_cosmic()
show_cosmic("legacy")
show_cosmic("SBS")
show_cosmic("DBS")
show_cosmic("ID")
show_cosmic("SBS1")
show_cosmic("DBS2")
show_cosmic("ID3")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_hyper_mutation.R
\name{handle_hyper_mutation}
\alias{handle_hyper_mutation}
\title{Handle Hypermutant Samples}
\usage{
handle_hyper_mutation(nmf_matrix)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}
}
\value{
a \code{matrix}.
}
\description{
This can be used for SNV/INDEL count matrix. For copy number analysis,
please skip it.
}
\references{
Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
Nature genetics 48.6 (2016): 600.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simulated_catalogs}
\alias{simulated_catalogs}
\title{A List of Simulated SBS-96 Catalog Matrix}
\format{
A list of matrix
}
\source{
Generate from code under data_raw/
}
\description{
Data from \doi{10.1038/s43018-020-0027-5}.
5 simulated mutation catalogs are used by the paper but only 4 are available.
The data are simulated from COSMIC mutational signatures 1, 2, 3, 5, 6, 8,
12, 13, 17 and 18. Each sample is a linear combination of 5 randomly selected
signatures with the addiction of Poisson noise. The number of mutation in
each sample is randomly selected between 1,000 and 50,000 mutations, in log
scale so that a lower number of mutations is more likely to be selected.
The proportion of each signature in each sample is also random.
}
\examples{
data(simulated_catalogs)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_fit_bootstrap_batch.R
\name{sig_fit_bootstrap_batch}
\alias{sig_fit_bootstrap_batch}
\title{Exposure Instability Analysis of Signature Exposures with Bootstrapping}
\usage{
sig_fit_bootstrap_batch(
  catalogue_matrix,
  methods = c("QP"),
  n = 100L,
  min_count = 1L,
  p_val_thresholds = c(0.05),
  use_parallel = FALSE,
  seed = 123456L,
  job_id = NULL,
  result_dir = tempdir(),
  ...
)
}
\arguments{
\item{catalogue_matrix}{a numeric matrix \code{V} with row representing components and
columns representing samples, typically you can get \code{nmf_matrix} from \code{sig_tally()} and
transpose it by \code{t()}.}

\item{methods}{a subset of \code{c("NNLS", "QP", "SA")}.}

\item{n}{the number of bootstrap replicates.}

\item{min_count}{minimal exposure in a sample, default is 1. Any patient has total exposure less
than this value will be filtered out.}

\item{p_val_thresholds}{a vector of relative exposure threshold for calculating p values.}

\item{use_parallel}{if \code{TRUE}, use parallel computation based on \strong{furrr} package.
It can also be an integer for specifying cores.}

\item{seed}{random seed to reproduce the result.}

\item{job_id}{a job ID, default is \code{NULL}, can be a string. When not \code{NULL}, all bootstrapped results
will be saved to local machine location defined by \code{result_dir}. This is very useful for running
more than 10 times for more than 100 samples.}

\item{result_dir}{see above, default is temp directory defined by R.}

\item{...}{other common parameters passing to \link{sig_fit_bootstrap}, including
\code{sig}, \code{sig_index}, \code{sig_db}, \code{db_type}, \code{mode}, \code{auto_reduce} etc.}
}
\value{
a \code{list} of \code{data.table}.
}
\description{
Read \link{sig_fit_bootstrap} for more option setting.
}
\examples{
W <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
colnames(W) <- c("sig1", "sig2")
W <- apply(W, 2, function(x) x / sum(x))

H <- matrix(c(2, 5, 3, 6, 1, 9, 1, 2), ncol = 4)
colnames(H) <- paste0("samp", 1:4)

V <- W \%*\% H
V

if (requireNamespace("quadprog")) {
  z10 <- sig_fit_bootstrap_batch(V, sig = W, n = 10)
  z10
}
}
\seealso{
\link{sig_fit}, \link{sig_fit_bootstrap}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cosmic_signature_profile.R
\name{show_cosmic_sig_profile}
\alias{show_cosmic_sig_profile}
\title{Plot Reference (Mainly COSMIC) Signature Profile}
\usage{
show_cosmic_sig_profile(
  sig_index = NULL,
  show_index = TRUE,
  sig_db = "legacy",
  ...
)
}
\arguments{
\item{sig_index}{a vector for signature index. "ALL" for all signatures.}

\item{show_index}{if \code{TRUE}, show valid indices.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}

\item{...}{other arguments passing to \link{show_sig_profile}.}
}
\value{
a \code{ggplot} object
}
\description{
Plot Reference (Mainly COSMIC) Signature Profile
}
\examples{
show_cosmic_sig_profile()
show_cosmic_sig_profile(sig_db = "SBS")
show_cosmic_sig_profile(sig_index = 1:5)
show_cosmic_sig_profile(sig_db = "SBS", sig_index = c("10a", "17a"))

gg <- show_cosmic_sig_profile(sig_index = 1:5)
gg$aetiology
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_profile_loop.R
\name{show_sig_profile_loop}
\alias{show_sig_profile_loop}
\title{Show Signature Profile with Loop Way}
\usage{
show_sig_profile_loop(
  Signature,
  sig_names = NULL,
  ncol = 1,
  nrow = NULL,
  x_lab = "Components",
  ...
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix with row representing components (motifs) and column
representing signatures (column names must start with 'Sig').}

\item{sig_names}{subset signatures or set name of signatures, can be a character vector.
Default is \code{NULL}, prefix 'Sig' plus number is used.}

\item{ncol}{(optional) Number of columns in the plot grid.}

\item{nrow}{(optional) Number of rows in the plot grid.}

\item{x_lab}{x axis lab.}

\item{...}{other parameters but \code{sig_order} passing to \link{show_sig_profile}.}
}
\value{
a \code{ggplot} result from \code{cowplot::plot_grid()}.
}
\description{
Show Signature Profile with Loop Way
}
\examples{

load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p1 <- show_sig_profile_loop(sig2, mode = "SBS")
p1
p2 <- show_sig_profile_loop(sig2, mode = "SBS", style = "cosmic", sig_names = c("A", "B", "C"))
p2
}
\seealso{
\link{show_sig_profile}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_intersect_size.R
\name{get_intersect_size}
\alias{get_intersect_size}
\title{Get Overlap Size between Interval x and y}
\usage{
get_intersect_size(x.start, x.end, y.start, y.end)
}
\arguments{
\item{x.start}{start position of interval x.}

\item{x.end}{start position of interval x.}

\item{y.start}{start position of interval x.}

\item{y.end}{start position of interval x.}
}
\value{
a numeric vector.
}
\description{
Get Overlap Size between Interval x and y
}
\examples{
o1 <- get_intersect_size(1, 5, 3, 20)
o1
o2 <- get_intersect_size(3, 20, 1, 10)
o2
o3 <- get_intersect_size(c(1, 2, 1), c(10, 4, 6), c(4, 2, 5), c(10, 3, 22))
o3

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_tally.R
\name{sig_tally}
\alias{sig_tally}
\alias{sig_tally.CopyNumber}
\alias{sig_tally.RS}
\alias{sig_tally.MAF}
\title{Tally a Genomic Alteration Object}
\usage{
sig_tally(object, ...)

\method{sig_tally}{CopyNumber}(
  object,
  method = "Wang",
  ignore_chrs = NULL,
  indices = NULL,
  add_loh = FALSE,
  feature_setting = sigminer::CN.features,
  cores = 1,
  keep_only_matrix = FALSE,
  ...
)

\method{sig_tally}{RS}(object, keep_only_matrix = FALSE, ...)

\method{sig_tally}{MAF}(
  object,
  mode = c("SBS", "DBS", "ID", "ALL"),
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  genome_build = NULL,
  add_trans_bias = FALSE,
  ignore_chrs = NULL,
  use_syn = TRUE,
  keep_only_matrix = FALSE,
  ...
)
}
\arguments{
\item{object}{a \link{CopyNumber} object or \link{MAF} object or SV object (from \link{read_sv_as_rs}).}

\item{...}{custom setting for operating object. Detail see S3 method for
corresponding class (e.g. \code{CopyNumber}).}

\item{method}{method for feature classification, can be one of
"Wang" ("W"), "S" (for method described in Steele et al. 2019).}

\item{ignore_chrs}{Chromsomes to ignore from analysis. e.g. chrX and chrY.}

\item{indices}{integer vector indicating segments to keep.}

\item{add_loh}{flag to add LOH classifications.}

\item{feature_setting}{a \code{data.frame} used for classification.
\strong{Only used when method is "Wang" ("W")}.
Default is \link{CN.features}. Users can also set custom input with "feature",
"min" and "max" columns available. Valid features can be printed by
\code{unique(CN.features$feature)}.}

\item{cores}{number of computer cores to run this task.
You can use \code{\link[future:re-exports]{future::availableCores()}} function to check how
many cores you can use.}

\item{keep_only_matrix}{if \code{TRUE}, keep only matrix for signature extraction.
For a \code{MAF} object, this will just return the most useful matrix.}

\item{mode}{type of mutation matrix to extract, can be one of 'SBS', 'DBS' and 'ID'.}

\item{ref_genome}{'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
'BSgenome.Mmusculus.UCSC.mm10',  'BSgenome.Mmusculus.UCSC.mm9', etc.}

\item{genome_build}{genome build 'hg19', 'hg38', 'mm9' or "mm10", if not set, guess it by \code{ref_genome}.}

\item{add_trans_bias}{if \code{TRUE}, consider transcriptional bias categories.
'T:' for Transcribed (the variant is on the transcribed strand);
'U:' for Un-transcribed (the variant is on the untranscribed strand);
'B:' for Bi-directional (the variant is on both strand and is transcribed either way);
'N:' for Non-transcribed (the variant is in a non-coding region and is untranslated);
'Q:' for Questionable.
\strong{NOTE}: the result counts of 'B' and 'N' labels are a little different from
SigProfilerMatrixGenerator, the reason is unknown (may be caused by annotation file).}

\item{use_syn}{Logical. If \code{TRUE}, include synonymous variants in analysis.}
}
\value{
a \code{list} contains a \code{matrix} used for NMF de-composition.
}
\description{
Tally a variation object like \link{MAF}, \link{CopyNumber} and return a matrix for NMF de-composition and more.
This is a generic function,
so it can be further extended to other mutation cases.
\strong{Please read details about how to set sex for identifying copy number signatures}.
Please read \url{https://osf.io/s93d5/} for the generation of SBS, DBS and ID (INDEL)
components.
}
\details{
For identifying copy number signatures, we have to derive copy number
features firstly. Due to the difference of copy number values in sex chromosomes
between male and female, we have to do an extra step \strong{if we don't want to
ignore them}.

I create two options to control this, the default values are shown as
the following, you can use the same way to set (per R session).

\code{options(sigminer.sex = "female", sigminer.copynumber.max = NA_integer_)}
\itemize{
\item If your cohort are all females, you can totally ignore this.
\item If your cohort are all males, set \code{sigminer.sex} to 'male' and
\code{sigminer.copynumber.max} to a proper value (the best is consistent
with \link{read_copynumber}).
\item If your cohort contains both males and females, set \code{sigminer.sex}
as a \code{data.frame} with two columns "sample" and "sex". And
set \code{sigminer.copynumber.max} to a proper value (the best is consistent
with \link{read_copynumber}).
}
}
\section{Methods (by class)}{
\itemize{
\item \code{CopyNumber}: Returns copy number features, components and component-by-sample matrix

\item \code{RS}: Returns genome rearrangement sample-by-component matrix

\item \code{MAF}: Returns SBS mutation sample-by-component matrix and APOBEC enrichment
}}

\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))
\donttest{
# Use method designed by Wang, Shixiang et al.
cn_tally_W <- sig_tally(cn, method = "W")
}
# Use method designed by Steele et al.
# See example in read_copynumber
\donttest{
# Prepare SBS signature analysis
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    use_syn = TRUE
  )
  mt_tally$nmf_matrix[1:5, 1:5]

  ## Use strand bias categories
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    use_syn = TRUE, add_trans_bias = TRUE
  )
  ## Test it by enrichment analysis
  enrich_component_strand_bias(mt_tally$nmf_matrix)
  enrich_component_strand_bias(mt_tally$all_matrices$SBS_24)
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}
}
}
\references{
Wang, Shixiang, et al. "Copy number signature analyses in prostate cancer reveal
distinct etiologies and clinical outcomes." medRxiv (2020).

Steele, Christopher D., et al. "Undifferentiated sarcomas develop through
distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.

Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.

Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702.

Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB: SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics 2019, 20:685 https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2
}
\seealso{
\link{sig_estimate} for estimating signature number for \link{sig_extract},
\link{sig_auto_extract} for extracting signatures using automatic relevance determination technique.
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_group_distribution.R
\name{show_group_distribution}
\alias{show_group_distribution}
\title{Show Groupped Variable Distribution}
\usage{
show_group_distribution(
  data,
  gvar,
  dvar,
  fun = stats::median,
  order_by_fun = FALSE,
  alpha = 0.8,
  g_label = "label",
  g_angle = 0,
  g_position = "top",
  point_size = 1L,
  segment_size = 1L,
  segment_color = "red",
  xlab = NULL,
  ylab = NULL,
  nrow = 1L,
  background_color = c("#DCDCDC", "#F5F5F5")
)
}
\arguments{
\item{data}{a \code{data.frame}.}

\item{gvar}{a group variable name/index.}

\item{dvar}{a distribution variable name/index.}

\item{fun}{a function to summarize, default is \link[stats:median]{stats::median}, can also be \link{mean}.}

\item{order_by_fun}{if \code{TRUE}, reorder the groups by summary measure computed
by argument \code{fun}.}

\item{alpha}{alpha for points, range from 0 to 1.}

\item{g_label}{a string 'label' (default) for labeling with sample size,
or 'norm' to show just group name, or a named vector to set facet labels.}

\item{g_angle}{angle for facet labels, default is \code{0}.}

\item{g_position}{position for facet labels, default is 'top', can also
be 'bottom'.}

\item{point_size}{size of point.}

\item{segment_size}{size of segment.}

\item{segment_color}{color of segment.}

\item{xlab}{title for x axis.}

\item{ylab}{title for y axis.}

\item{nrow}{number of row.}

\item{background_color}{background color for plot panel.}
}
\value{
a \code{ggplot} object.
}
\description{
This is a general function, it can be used in any proper analysis.
}
\examples{
set.seed(1234)
data <- data.frame(
  yval = rnorm(120),
  gr = c(rep("A", 50), rep("B", 40), rep("C", 30))
)
p <- show_group_distribution(data,
  gvar = 2, dvar = 1,
  g_label = "norm",
  background_color = "grey"
)
p
p2 <- show_group_distribution(data,
  gvar = "gr", dvar = "yval",
  g_position = "bottom",
  order_by_fun = TRUE,
  alpha = 0.3
)
p2

# Set custom group names
p3 <- show_group_distribution(data,
  gvar = 2, dvar = 1,
  g_label = c("A" = "X", "B" = "Y", "C" = "Z")
)
p3
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_group_comparison.R
\name{show_group_comparison}
\alias{show_group_comparison}
\title{Plot Group Comparison Result}
\usage{
show_group_comparison(
  group_comparison,
  xlab = "group",
  ylab_co = NA,
  legend_title_ca = NA,
  legend_position_ca = "bottom",
  set_ca_sig_yaxis = FALSE,
  set_ca_custom_xlab = FALSE,
  show_pvalue = TRUE,
  ca_p_threshold = 0.01,
  method = "wilcox.test",
  p.adjust.method = "fdr",
  base_size = 12,
  font_size_x = 12,
  text_angle_x = 30,
  text_hjust_x = 0.2,
  ...
)
}
\arguments{
\item{group_comparison}{a \code{list} from result of \link{get_group_comparison} function.}

\item{xlab}{lab name of x axis for all plots. if it is \code{NA}, remove title for x axis.}

\item{ylab_co}{lab name of y axis for plots of continuous type data. Of note,
this argument should be a character vector has same length as \code{group_comparison},
the location for categorical type data should mark with \code{NA}.}

\item{legend_title_ca}{legend title for plots of categorical type data.}

\item{legend_position_ca}{legend position for plots of categorical type data.
Of note,
this argument should be a character vector has same length as \code{group_comparison},
the location for continuous type data should mark with \code{NA}.}

\item{set_ca_sig_yaxis}{if \code{TRUE}, use y axis to show signature proportion instead of
variable proportion.}

\item{set_ca_custom_xlab}{only works when \code{set_ca_sig_yaxis} is \code{TRUE}. If
\code{TRUE}, set x labels using input \code{xlab}, otherwise variable names will be used.}

\item{show_pvalue}{if \code{TRUE}, show p values.}

\item{ca_p_threshold}{a p threshold for categorical variables, default is 0.01.
A p value less than 0.01 will be shown as \code{P < 0.01}.}

\item{method}{a character string indicating which method to be used for comparing means.
It can be 't.test', 'wilcox.test' etc..}

\item{p.adjust.method}{correction method, default is 'fdr'. Run \code{p.adjust.methods} to
see all available options.}

\item{base_size}{overall font size.}

\item{font_size_x}{font size for x.}

\item{text_angle_x}{text angle for x.}

\item{text_hjust_x}{adjust x axis text}

\item{...}{other paramters pass to \code{\link[ggpubr:compare_means]{ggpubr::compare_means()}} or \code{\link[ggpubr:stat_compare_means]{ggpubr::stat_compare_means()}}
according to the specified \code{method}.}
}
\value{
a \code{list} of \code{ggplot} objects.
}
\description{
Using result data from \link{get_group_comparison}, this function plots
genotypes/phenotypes comparison between signature groups using \strong{ggplot2} package and return
a list of \code{ggplot} object contains individual and combined plots. The combined
plot is easily saved to local using \code{\link[cowplot:save_plot]{cowplot::save_plot()}}. Of note, default fisher
test p values are shown for categorical data and fdr values are shown for
continuous data.
}
\examples{
load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
  package = "sigminer", mustWork = TRUE
))

# Assign samples to clusters
groups <- get_groups(sig, method = "k-means")

set.seed(1234)

groups$prob <- rnorm(10)
groups$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(groups), replace = TRUE)

# Compare groups (filter NAs for categorical coloumns)
groups.cmp <- get_group_comparison(groups[, -1],
  col_group = "group",
  cols_to_compare = c("prob", "new_group"),
  type = c("co", "ca"), verbose = TRUE
)

# Compare groups (Set NAs of categorical columns to 'Rest')
groups.cmp2 <- get_group_comparison(groups[, -1],
  col_group = "group",
  cols_to_compare = c("prob", "new_group"),
  type = c("co", "ca"), NAs = "Rest", verbose = TRUE
)

show_group_comparison(groups.cmp)

ggcomp <- show_group_comparison(groups.cmp2)
ggcomp$co_comb
ggcomp$ca_comb
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{CN.features}
\alias{CN.features}
\title{Classification Table of Copy Number Features Devised by Wang et al. for Method 'W'}
\format{
A \code{data.table} with "sigminer.features" class name
}
\source{
Generate from code under data_raw/
}
\description{
Classification Table of Copy Number Features Devised by Wang et al. for Method 'W'
}
\examples{
data(CN.features)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigminer.R
\docType{package}
\name{sigminer}
\alias{sigminer}
\title{sigminer: Extract, Analyze and Visualize Signatures for Genomic Variations}
\description{
\itemize{
\item Author: \href{https://shixiangwang.github.io/home/}{Shixiang Wang} (\href{mailto:w_shixiang@163.com}{w_shixiang@163.com})
\item Please go to \url{https://shixiangwang.github.io/sigminer-doc/} for full vignette.
\item Please go to \url{https://shixiangwang.github.io/sigminer/reference/index.html}
for organized documentation of functions and datasets.
\item Result visualization for \link{MAF} is provide by \strong{maftools} package,
please read its \href{https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html}{vignette}.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_seg_table.R
\name{transform_seg_table}
\alias{transform_seg_table}
\title{Transform Copy Number Table}
\usage{
transform_seg_table(
  data,
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  ref_type = c("cytoband", "gene"),
  values_fill = NA,
  values_fn = function(x, ...) {     round(mean(x, ...)) },
  resolution_factor = 1L
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{ref_type}{annotation data type used for constructing matrix.}

\item{values_fill}{Optionally, a (scalar) value that specifies what each
\code{value} should be filled in with when missing.

This can be a named list if you want to apply different aggregations
to different value columns.}

\item{values_fn}{Optionally, a function applied to the \code{value} in each cell
in the output. You will typically use this when the combination of
\code{id_cols} and \code{value} column does not uniquely identify an observation.

This can be a named list if you want to apply different aggregations
to different value columns.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}
}
\value{
a \code{data.table}.
}
\description{
Transform Copy Number Table
}
\examples{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))
# Compute the mean segVal in each cytoband
x <- transform_seg_table(cn, resolution_factor = 1)
x
# Compute the mean segVal in each half-cytoband
x2 <- transform_seg_table(cn, resolution_factor = 2)
x2
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chromsize.hg38}
\alias{chromsize.hg38}
\title{Chromosome Size of Genome Build hg38}
\format{
A data.frame
}
\source{
Generate from UCSC gold path
}
\description{
Chromosome Size of Genome Build hg38
}
\examples{
data(chromsize.hg38)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_auto_extract.R
\name{sig_auto_extract}
\alias{sig_auto_extract}
\title{Extract Signatures through the Automatic Relevance Determination Technique}
\usage{
sig_auto_extract(
  nmf_matrix = NULL,
  result_prefix = "BayesNMF",
  destdir = tempdir(),
  method = c("L1W.L2H", "L1KL", "L2KL"),
  strategy = c("stable", "optimal", "ms"),
  ref_sigs = NULL,
  K0 = 25,
  nrun = 10,
  niter = 2e+05,
  tol = 1e-07,
  cores = 1,
  optimize = FALSE,
  skip = FALSE,
  recover = FALSE
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{result_prefix}{prefix for result data files.}

\item{destdir}{path to save data runs, default is \code{tempdir()}.}

\item{method}{default is "L1W.L2H", which uses an exponential prior for W and
a half-normal prior for H (This method is used by PCAWG project, see reference #3).
You can also use "L1KL" to set expoential priors for both W and H, and "L2KL" to
set half-normal priors for both W and H. The latter two methods are originally
implemented by \href{https://software.broadinstitute.org/cancer/cga/msp}{SignatureAnalyzer software}.}

\item{strategy}{the selection strategy for returned data. Set 'stable' for getting optimal
result from the most frequent K. Set 'optimal' for getting optimal result from all Ks.
Set 'ms' for getting result with maximum mean cosine similarity with provided reference
signatures. See \code{ref_sigs} option for details.
If you want select other solution, please check \link{get_bayesian_result}.}

\item{ref_sigs}{A Signature object or matrix or string for specifying
reference signatures, only used when \code{strategy = 'ms'}.
See \code{Signature} and \code{sig_db} options in \link{get_sig_similarity} for details.}

\item{K0}{number of initial signatures.}

\item{nrun}{number of independent simulations.}

\item{niter}{the maximum number of iterations.}

\item{tol}{tolerance for convergence.}

\item{cores}{number of cpu cores to run NMF.}

\item{optimize}{if \code{TRUE}, then refit the denovo signatures with QP method, see \link{sig_fit}.}

\item{skip}{if \code{TRUE}, it will skip running a previous stored result. This can be used to
extend run times, e.g. you try running 10 times firstly and then you want to extend it to
20 times.}

\item{recover}{if \code{TRUE}, try to recover result from previous runs based on input \code{result_prefix},
\code{destdir} and \code{nrun}. This is pretty useful for reproducing result. Please use \code{skip} if you want
to recover an unfinished job.}
}
\value{
a \code{list} with \code{Signature} class.
}
\description{
A bayesian variant of NMF algorithm to enable optimal inferences for the
number of signatures through the automatic relevance determination technique.
This functions delevers highly interpretable and sparse representations for
both signature profiles and attributions at a balance between data fitting and
model complexity (this method may introduce more signatures than expected,
especially for copy number signatures (thus \strong{I don't recommend you to use this feature
to extract copy number signatures})). See detail part and references for more.
}
\details{
There are three methods available in this function: "L1W.L2H", "L1KL" and "L2KL".
They use different priors for the bayesian variant of NMF algorithm
(see \code{method} parameter) written by reference #1 and implemented in
\href{https://software.broadinstitute.org/cancer/cga/msp}{SignatureAnalyzer software}
(reference #2).

I copied source code for the three methods from Broad Institute and supplementary
files of reference #3, and wrote this higher function. It is more friendly for users
to extract, visualize and analyze signatures by combining with other powerful functions
in \strong{sigminer} package. Besides, I implemented parallel computation to speed up
the calculation process and a similar input and output structure like \code{\link[=sig_extract]{sig_extract()}}.
}
\examples{
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))
res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
# At default, all run files are stored in tempdir()
dir(tempdir(), pattern = "Test_copynumber")
\donttest{
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
mt_tally <- sig_tally(
  laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE
)

x <- sig_auto_extract(mt_tally$nmf_matrix,
  strategy = "ms", nrun = 3, ref_sigs = "legacy"
)
x
}
}
\references{
Tan, Vincent YF, and Cédric Févotte. "Automatic relevance determination in nonnegative matrix factorization with the/spl beta/-divergence."
IEEE Transactions on Pattern Analysis and Machine Intelligence 35.7 (2012): 1592-1605.

Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
Nature genetics 48.6 (2016): 600.

Alexandrov, Ludmil, et al. "The repertoire of mutational signatures in human cancer." BioRxiv (2018): 322859.
}
\seealso{
\link{sig_tally} for getting variation matrix,
\link{sig_extract} for extracting signatures using \strong{NMF} package, \link{sig_estimate} for
estimating signature number for \link{sig_extract}.
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_profile.R
\name{show_cn_profile}
\alias{show_cn_profile}
\title{Show Sample Copy Number Profile}
\usage{
show_cn_profile(
  data,
  samples = NULL,
  show_n = NULL,
  show_title = FALSE,
  show_labels = NULL,
  chrs = paste0("chr", 1:22),
  position = NULL,
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  ylim = NULL,
  nrow = NULL,
  ncol = NULL,
  return_plotlist = FALSE
)
}
\arguments{
\item{data}{a \link{CopyNumber} object or a \code{data.frame} containing at least 'chromosome', 'start',
'end', 'segVal' these columns.}

\item{samples}{default is NULL, can be a chracter vector representing multiple samples. If \code{data} argument
is a \code{data.frame}, a column called \code{sample} must exist.}

\item{show_n}{number of samples to show, this is used for checking.}

\item{show_title}{if \code{TRUE}, show title for multiple samples.}

\item{show_labels}{one of \code{NULL}, "s" (for labelling short segments < 1e7)
or "a" (all segments).}

\item{chrs}{chromosomes start with 'chr'.}

\item{position}{a position range, e.g. \code{"chr1:3218923-116319008"}. Only data
overlaps with this range will be shown.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{ylim}{limites for y axis.}

\item{nrow}{number of rows in the plot grid when multiple samples are selected.}

\item{ncol}{number of columns in the plot grid when multiple samples are selected.}

\item{return_plotlist}{default is \code{FALSE}, if \code{TRUE}, return a plot list instead of a combined plot.}
}
\value{
a \code{ggplot} object or a \code{list}
}
\description{
Sometimes it is very useful to check details about copy number profile for one or multiple
samples. This function is designed to do this job and can be further modified by \strong{ggplot2}
related packages.
}
\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

p <- show_cn_profile(cn, nrow = 2, ncol = 1)
p
\donttest{
p2 <- show_cn_profile(cn,
  nrow = 2, ncol = 1,
  position = "chr1:3218923-116319008"
)
p2
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\name{subset.CopyNumber}
\alias{subset.CopyNumber}
\title{Subsetting CopyNumber object}
\usage{
\method{subset}{CopyNumber}(x, subset = TRUE, ...)
}
\arguments{
\item{x}{a \link{CopyNumber} object to be subsetted.}

\item{subset}{logical expression indicating rows to keep.}

\item{...}{further arguments to be passed to or from other methods.
Useless here.}
}
\value{
a \link{CopyNumber} object
}
\description{
Subset \code{data} slot of \link{CopyNumber} object, un-selected rows will move to
dropoff.segs slot, annotation slot will update in the same way.
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_copynumber.R
\name{read_copynumber}
\alias{read_copynumber}
\title{Read Absolute Copy Number Profile}
\usage{
read_copynumber(
  input,
  pattern = NULL,
  ignore_case = FALSE,
  seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
  samp_col = "sample",
  add_loh = FALSE,
  loh_min_len = 10000,
  loh_min_frac = 0.05,
  join_adj_seg = TRUE,
  skip_annotation = FALSE,
  use_all = add_loh,
  min_segnum = 0L,
  max_copynumber = 20L,
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  genome_measure = c("called", "wg"),
  complement = FALSE,
  ...
)
}
\arguments{
\item{input}{a \code{data.frame} or a file or a directory contains copy number profile.}

\item{pattern}{an optional regular expression used to select part of files if
\code{input} is a directory, more detail please see \code{\link[=list.files]{list.files()}} function.}

\item{ignore_case}{logical. Should pattern-matching be case-insensitive?}

\item{seg_cols}{four strings used to specify chromosome, start position,
end position and copy number value in \code{input}, respectively.
Default use names from ABSOLUTE calling result.}

\item{samp_col}{a character used to specify the sample column name. If \code{input}
is a directory and cannot find \code{samp_col}, sample names will use file names
(set this parameter to \code{NULL} is recommended in this case).}

\item{add_loh}{if \code{TRUE}, add LOH labels to segments. \strong{NOTE} a column
'minor_cn' must exist to indicate minor allele copy number value.
Sex chromosome will not be labeled.}

\item{loh_min_len}{The length cut-off for labeling a segment as 'LOH'.
Default is \verb{10Kb}.}

\item{loh_min_frac}{When \code{join_adj_seg} set to \code{TRUE}, only the length fraction
of LOH region is larger than this value will be labeled as 'LOH'.
Default is 30\%.}

\item{join_adj_seg}{if \code{TRUE} (default), join adjacent segments with
same copy number value. This is helpful for precisely count the number of breakpoint.
When set \code{use_all=TRUE}, the mean function will be applied to extra numeric columns
and unique string columns will be pasted by comma for joined records.}

\item{skip_annotation}{if \code{TRUE}, skip annotation step, it may affect some analysis
and visualization functionality, but speed up reading data.}

\item{use_all}{default is \code{FALSE}. If \code{True}, use all columns from raw input.}

\item{min_segnum}{minimal number of copy number segments within a sample.}

\item{max_copynumber}{bigger copy number within a sample will be reset to this value.}

\item{genome_build}{genome build version, should be 'hg19', 'hg38', 'mm9' or 'mm10'.}

\item{genome_measure}{default is 'called', can be 'wg' or 'called'.
Set 'called' will use called segments size to compute total size for CNA burden calculation,
this option is useful for WES and target sequencing.
Set 'wg' will use autosome size from genome build, this option is useful for WGS, SNP etc..}

\item{complement}{if \code{TRUE}, complement chromosome (except 'Y') does not show in input data
with normal copy 2.}

\item{...}{other parameters pass to \code{\link[data.table:fread]{data.table::fread()}}}
}
\value{
a \link{CopyNumber} object.
}
\description{
Read \strong{absolute} copy number profile for preparing CNV signature
analysis. See detail part of \code{\link[=sig_tally]{sig_tally()}} to see how to handle sex to get correct
summary.
}
\examples{
# Load toy dataset of absolute copynumber profile
load(system.file("extdata", "toy_segTab.RData",
  package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  genome_build = "hg19", complement = FALSE
)
cn
cn_subset <- subset(cn, sample == "TCGA-DF-A2KN-01A-11D-A17U-01")

# Add LOH
set.seed(1234)
segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
cn <- read_copynumber(segTabs,
  seg_cols = c("chromosome", "start", "end", "segVal"),
  genome_measure = "wg", complement = TRUE, add_loh = TRUE
)
# Use tally method "S" (Steele et al.)
tally_s <- sig_tally(cn, method = "S")

tab_file <- system.file("extdata", "metastatic_tumor.segtab.txt",
  package = "sigminer", mustWork = TRUE
)
cn2 <- read_copynumber(tab_file)
cn2
}
\seealso{
\link{read_maf} for reading mutation data to \link{MAF} object.
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_distribution.R
\name{show_cn_distribution}
\alias{show_cn_distribution}
\title{Show Copy Number Distribution either by Length or Chromosome}
\usage{
show_cn_distribution(
  data,
  rm_normal = TRUE,
  mode = c("ld", "cd"),
  fill = FALSE,
  scale_chr = TRUE,
  base_size = 14
)
}
\arguments{
\item{data}{a \link{CopyNumber} object.}

\item{rm_normal}{logical. Whether remove normal copy (i.e. "segVal" equals 2), default is \code{TRUE}.}

\item{mode}{either "ld" for distribution by CN length or "cd" for distribution by chromosome.}

\item{fill}{when \code{mode} is "cd" and \code{fill} is \code{TRUE}, plot percentage instead of count.}

\item{scale_chr}{logical. If \code{TRUE}, normalize count to per Megabase unit.}

\item{base_size}{overall font size.}
}
\value{
a \code{ggplot} object
}
\description{
Visually summarize copy number distribution either by copy number segment length
or chromosome. Input is a \link{CopyNumber} object, \code{genome_build} option will
read from \code{genome_build} slot of object.
}
\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))
# Plot distribution
p1 <- show_cn_distribution(cn)
p1
p2 <- show_cn_distribution(cn, mode = "cd")
p2
p3 <- show_cn_distribution(cn, mode = "cd", fill = TRUE)
p3
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_adj_p.R
\name{get_adj_p}
\alias{get_adj_p}
\title{Get Adjust P Values from Group Comparison}
\source{
https://github.com/kassambara/ggpubr/issues/143
}
\usage{
get_adj_p(
  data,
  .col,
  .grp = "Sample",
  comparisons = NULL,
  method = "wilcox.test",
  p.adjust.method = "fdr",
  p.digits = 3L,
  ...
)
}
\arguments{
\item{data}{a \code{data.frame} containing column for groups and column for comparison.}

\item{.col}{column name for comparison.}

\item{.grp}{column name for groups.}

\item{comparisons}{Default is \code{NULL}, use all combination in group column.
It can be a list of length-2 vectors. The entries in the vector are either
the names of 2 values on the x-axis or the 2 integers that correspond to the
index of the groups of interest, to be compared.}

\item{method}{a character string indicating which method to be used for comparing means.
It can be 't.test', 'wilcox.test' etc..}

\item{p.adjust.method}{correction method, default is 'fdr'. Run \code{p.adjust.methods} to
see all available options.}

\item{p.digits}{how many significant digits are to be used.}

\item{...}{other arguments passed to \code{\link[ggpubr:compare_means]{ggpubr::compare_means()}}}
}
\value{
a \code{data.frame} containing comparison result
}
\description{
Setting \code{aes(label=..p.adj..)} in \code{\link[ggpubr:compare_means]{ggpubr::compare_means()}} does not
show adjust p values. The returned result of this function can be combined with \code{\link[ggpubr:stat_pvalue_manual]{ggpubr::stat_pvalue_manual()}} to fix
this problem.
}
\details{
More info see \code{\link[ggpubr:compare_means]{ggpubr::compare_means()}}, \code{\link[ggpubr:stat_compare_means]{ggpubr::stat_compare_means()}} and \code{\link[stats:p.adjust]{stats::p.adjust()}}.
}
\examples{
library(ggpubr)
# T-test
stat.test <- compare_means(
  len ~ dose,
  data = ToothGrowth,
  method = "t.test",
  p.adjust.method = "fdr"
)
stat.test
# Create a simple box plot
p <- ggboxplot(ToothGrowth, x = "dose", y = "len")
p

# Add p values
my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
p + stat_compare_means(method = "t.test", comparisons = my_comparisons)

# Try adding adjust p values
# proposed by author of ggpubr
# however it does not work
p + stat_compare_means(aes(label = ..p.adj..), method = "t.test", comparisons = my_comparisons)

# Solution:
# calculate adjust p values and their location
# then use stat_pvalue_manual() function
p_adj <- get_adj_p(ToothGrowth, .col = "len", .grp = "dose")
p_adj
p + stat_pvalue_manual(p_adj, label = "p.adj")

# Show selected comparisons
# Of note, p value is ajusted
# for three comparisons, but only
# two are showed in figure
p_adj <- get_adj_p(ToothGrowth,
  .col = "len", .grp = "dose",
  comparisons = list(c("0.5", "1"), c("1", "2"))
)
p + stat_pvalue_manual(p_adj, label = "p.adj")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_h_arrow.R
\name{add_h_arrow}
\alias{add_h_arrow}
\title{Add Horizontal Arrow with Text Label to a ggplot}
\usage{
add_h_arrow(
  p,
  x,
  y,
  label = "optimal number",
  space = 0.01,
  vjust = 0.3,
  seg_len = 0.1,
  arrow_len = unit(2, "mm"),
  arrow_type = c("closed", "open"),
  font_size = 5,
  font_family = c("serif", "sans", "mono"),
  font_face = c("plain", "bold", "italic")
)
}
\arguments{
\item{p}{a \code{ggplot}.}

\item{x}{position at x axis.}

\item{y}{position at y axis.}

\item{label}{text label.}

\item{space}{a small space between arrow and text.}

\item{vjust}{vertical adjustment, set to 0 to align with the bottom,
0.5 for the middle, and 1 (the default) for the top.}

\item{seg_len}{length of the arrow segment.}

\item{arrow_len}{length of the arrow.}

\item{arrow_type}{type of the arrow.}

\item{font_size}{font size.}

\item{font_family}{font family.}

\item{font_face}{font face.}
}
\value{
a \code{ggplot} object.
}
\description{
Add Horizontal Arrow with Text Label to a ggplot
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_freq_circos.R
\name{show_cn_freq_circos}
\alias{show_cn_freq_circos}
\title{Show Copy Number Variation Frequency Profile with Circos}
\usage{
show_cn_freq_circos(
  data,
  groups = NULL,
  cutoff = 2L,
  resolution_factor = 1L,
  title = c("AMP", "DEL"),
  chrs = paste0("chr", 1:22),
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  cols = NULL,
  plot_ideogram = TRUE,
  track_height = 0.5,
  ideogram_height = 1,
  ...
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{groups}{a named list or a column name for specifying groups.}

\item{cutoff}{copy number value cutoff for splitting data into AMP and DEL.
The values equal to cutoff are discarded. Default is \code{2}, you can also set
a length-2 vector, e.g. \code{c(2, 2)}.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}

\item{title}{length-2 titles for AMP and DEL.}

\item{chrs}{chromosomes start with 'chr'.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{cols}{length-2 colors for AMP and DEL.}

\item{plot_ideogram}{default is \code{TRUE}, show ideogram.}

\item{track_height}{track height in \code{mm} unit.}

\item{ideogram_height}{ideogram height in \code{mm} unit.}

\item{...}{other parameters passing to \link[circlize:circos.genomicLines]{circlize::circos.genomicLines}.}
}
\value{
Nothing.
}
\description{
Show Copy Number Variation Frequency Profile with Circos
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

show_cn_freq_circos(cn)
ss <- unique(cn@data$sample)
show_cn_freq_circos(cn, groups = list(a = ss[1:5], b = ss[6:10]), cols = c("red", "green"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cor.R
\name{show_cor}
\alias{show_cor}
\title{A Simple and General Way for Association Analysis}
\usage{
show_cor(
  data,
  x_vars = colnames(data),
  y_vars = x_vars,
  cor_method = "spearman",
  vis_method = "square",
  lab = TRUE,
  test = TRUE,
  hc_order = FALSE,
  p_adj = NULL,
  ...
)
}
\arguments{
\item{data}{a \code{data.frame}.}

\item{x_vars}{variables/column names shown in x axis.}

\item{y_vars}{variables/column names shown in y axis.}

\item{cor_method}{method for correlation, default is 'spearman'.}

\item{vis_method}{visualization method, default is 'square',
can also be 'circle'.}

\item{lab}{logical value. If TRUE, add correlation coefficient on the plot.}

\item{test}{if \code{TRUE}, run test for correlation and mark significance.}

\item{hc_order}{logical value. If \code{TRUE},
correlation matrix will be hc.ordered using \code{hclust} function.}

\item{p_adj}{p adjust method, see \link[stats:p.adjust]{stats::p.adjust} for details.}

\item{...}{other parameters passing to \code{ggcorrplot::ggcorrplot()}.}
}
\value{
a \code{ggplot} object
}
\description{
All variables must be continuous.
The matrix will be returned as an element of \code{ggplot} object.
This is basically a wrapper of R package
\href{https://github.com/kassambara/ggcorrplot}{ggcorrplot}.
}
\examples{
data("mtcars")
p1 <- show_cor(mtcars)
p2 <- show_cor(mtcars,
  x_vars = colnames(mtcars)[1:4],
  y_vars = colnames(mtcars)[5:8]
)
p3 <- show_cor(mtcars, vis_method = "circle", p_adj = "fdr")
p1
p1$cor
p2
p3

## Auto detect problem variables
mtcars$xx <- 0L
p4 <- show_cor(mtcars)
p4
}
\seealso{
\link{show_sig_feature_corrplot} for specific and more powerful
association analysis and visualization.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_group_profile.R
\name{show_cn_group_profile}
\alias{show_cn_group_profile}
\title{Show Summary Copy Number Profile for Sample Groups}
\usage{
show_cn_group_profile(
  data,
  groups = NULL,
  fill_area = TRUE,
  cols = NULL,
  chrs = paste0("chr", c(1:22, "X")),
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  cutoff = 2L,
  resolution_factor = 1L,
  force_y_limit = TRUE,
  highlight_genes = NULL,
  repel = FALSE,
  nrow = NULL,
  ncol = NULL,
  return_plotlist = FALSE
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{groups}{a named list or a column name for specifying groups.}

\item{fill_area}{default is \code{TRUE}, fill area with colors.}

\item{cols}{length-2 colors for AMP and DEL.}

\item{chrs}{chromosomes start with 'chr'.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{cutoff}{copy number value cutoff for splitting data into AMP and DEL.
The values equal to cutoff are discarded. Default is \code{2}, you can also set
a length-2 vector, e.g. \code{c(2, 2)}.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}

\item{force_y_limit}{default is \code{TRUE}, force multiple plots}

\item{highlight_genes}{gene list to highlight.
have same y ranges. You can also set a length-2 numeric value.}

\item{repel}{if \code{TRUE} (default is \code{FALSE}), repel highlight genes to
avoid overlap.}

\item{nrow}{number of rows in the plot grid when multiple samples are selected.}

\item{ncol}{number of columns in the plot grid when multiple samples are selected.}

\item{return_plotlist}{default is \code{FALSE}, if \code{TRUE}, return a plot list instead of a combined plot.}
}
\value{
a (list of) \code{ggplot} object.
}
\description{
Show Summary Copy Number Profile for Sample Groups
}
\examples{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

p1 <- show_cn_group_profile(cn)
p1
\donttest{
ss <- unique(cn@data$sample)
p2 <- show_cn_group_profile(cn, groups = list(a = ss[1:5], b = ss[6:10]))
p2
p3 <- show_cn_group_profile(cn,
  groups = list(g1 = ss[1:5], g2 = ss[6:10]),
  force_y_limit = c(-1, 1), nrow = 2
)
p3

## Set custom cutoff for custom data
data <- cn@data
data$segVal <- data$segVal - 2L
p4 <- show_cn_group_profile(data,
  groups = list(g1 = ss[1:5], g2 = ss[6:10]),
  force_y_limit = c(-1, 1), nrow = 2,
  cutoff = c(0, 0)
)
p4

## Add highlight gene
p5 <- show_cn_group_profile(cn, highlight_genes = c("TP53", "EGFR"))
p5
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{centromeres.mm10}
\alias{centromeres.mm10}
\title{Location of Centromeres at Genome Build mm10}
\format{
A data.frame
}
\source{
Generate from \url{https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz}
}
\description{
Location of Centromeres at Genome Build mm10
}
\examples{
data(centromeres.mm10)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_groups.R
\name{show_groups}
\alias{show_groups}
\title{Show Signature Contribution in Clusters}
\usage{
show_groups(grp_dt, ...)
}
\arguments{
\item{grp_dt}{a result \code{data.table} from \link{get_groups}.}

\item{...}{parameters passing to \code{\link[=legend]{legend()}}, e.g. \code{x = "topleft"}.}
}
\value{
nothing.
}
\description{
See example section in \code{\link[=sig_fit]{sig_fit()}} for an examples.
}
\seealso{
\link{get_groups}, \link{sig_fit}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\docType{class}
\name{MAF-class}
\alias{MAF-class}
\alias{MAF}
\title{Class MAF}
\description{
S4 class for storing summarized MAF. It is from \code{maftools} package.
}
\details{
More about MAF object please see \href{https://github.com/PoisonAlien/maftools}{maftools}.
}
\section{Slots}{

\describe{
\item{\code{data}}{data.table of MAF file containing all non-synonymous variants.}

\item{\code{variants.per.sample}}{table containing variants per sample}

\item{\code{variant.type.summary}}{table containing variant types per sample}

\item{\code{variant.classification.summary}}{table containing variant classification per sample}

\item{\code{gene.summary}}{table containing variant classification per gene}

\item{\code{summary}}{table with basic MAF summary stats}

\item{\code{maf.silent}}{subset of main MAF containing only silent variants}

\item{\code{clinical.data}}{clinical data associated with each sample/Tumor_Sample_Barcode in MAF.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{transcript.hg38}
\alias{transcript.hg38}
\title{Merged Transcript Location at Genome Build hg38}
\format{
A \code{data.table}
}
\source{
from GENCODE release v33.
}
\description{
Merged Transcript Location at Genome Build hg38
}
\examples{
data(transcript.hg38)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_components.R
\name{show_cn_components}
\alias{show_cn_components}
\title{Show Copy Number Components}
\usage{
show_cn_components(
  parameters,
  method = "Wang",
  show_weights = TRUE,
  log_y = FALSE,
  return_plotlist = FALSE,
  base_size = 12,
  nrow = 2,
  align = "hv",
  ...
)
}
\arguments{
\item{parameters}{a \code{data.frame} contain parameter components, obtain this
from \link{sig_tally} function.}

\item{method}{method for feature classification, can be one of
"Wang" ("W"), "S" (for method described in Steele et al. 2019).}

\item{show_weights}{default is \code{TRUE}, show weights for each component.
Only used when method is "Macintyre".}

\item{log_y}{logical, if \code{TRUE}, show \code{log10} based y axis, only
works for input from "Wang" ("W") method.}

\item{return_plotlist}{if \code{TRUE}, return a list of ggplot objects but a combined plot.}

\item{base_size}{overall font size.}

\item{nrow}{(optional) Number of rows in the plot grid.}

\item{align}{(optional) Specifies whether graphs in the grid should be horizontally ("h") or
vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".}

\item{...}{other options pass to \code{\link[cowplot]{plot_grid}} function of \strong{cowplot} package.}
}
\value{
a \code{ggplot} object
}
\description{
Show classified components ("Wang" ("W") method) for copy number data.
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_similarity.R
\name{get_sig_similarity}
\alias{get_sig_similarity}
\title{Calculate Similarity between Identified Signatures and Reference Signatures}
\usage{
get_sig_similarity(
  Signature,
  Ref = NULL,
  sig_db = c("legacy", "SBS", "DBS", "ID", "TSB", "SBS_Nik_lab", "RS_Nik_lab",
    "RS_BRCA560", "RS_USARC", "CNS_USARC", "CNS_TCGA", "SBS_hg19", "SBS_hg38", "SBS_mm9",
    "SBS_mm10", "DBS_hg19", "DBS_hg38", "DBS_mm9", "DBS_mm10", "SBS_Nik_lab_Organ",
    "RS_Nik_lab_Organ", "latest_SBS_GRCh37", "latest_DBS_GRCh37", "latest_ID_GRCh37",
    "latest_SBS_GRCh38", "latest_DBS_GRCh38", "latest_SBS_mm9", "latest_DBS_mm9",
    "latest_SBS_mm10", "latest_DBS_mm10", "latest_SBS_rn6", "latest_DBS_rn6"),
  db_type = c("", "human-exome", "human-genome"),
  method = "cosine",
  normalize = c("row", "feature"),
  feature_setting = sigminer::CN.features,
  set_order = TRUE,
  pattern_to_rm = NULL,
  verbose = TRUE
)
}
\arguments{
\item{Signature}{a \code{Signature} object or a component-by-signature matrix/\code{data.frame}
(sum of each column is 1) or a normalized component-by-sample matrix/\code{data.frame}
(sum of each column is 1).
More please see examples.}

\item{Ref}{default is \code{NULL}, can be a same object as \code{Signature}.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}

\item{db_type}{only used when \code{sig_db} is enabled.
"" for keeping default, "human-exome" for transforming to exome frequency of component,
and "human-genome" for transforming to whole genome frequency of component.
Currently only works for 'SBS'.}

\item{method}{default is 'cosine' for cosine similarity.}

\item{normalize}{one of "row" and "feature". "row" is typically used
for common mutational signatures. "feature" is designed by me to use when input
are copy number signatures.}

\item{feature_setting}{a \code{data.frame} used for classification.
\strong{Only used when method is "Wang" ("W")}.
Default is \link{CN.features}. Users can also set custom input with "feature",
"min" and "max" columns available. Valid features can be printed by
\code{unique(CN.features$feature)}.}

\item{set_order}{if \code{TRUE}, order the return similarity matrix.}

\item{pattern_to_rm}{patterns for removing some features/components in similarity
calculation. A vector of component name is also accepted.
The remove operation will be done after normalization. Default is \code{NULL}.}

\item{verbose}{if \code{TRUE}, print extra info.}
}
\value{
a \code{list} containing smilarities, aetiologies if available, best match and RSS.
}
\description{
The reference signatures can be either a \code{Signature} object specified by \code{Ref} argument
or known COSMIC signatures specified by \code{sig_db} argument.
Two COSMIC databases are used for comparisons - "legacy" which includes 30 signaures,
and "SBS" - which includes updated/refined 65 signatures. This function is modified
from \code{compareSignatures()} in \strong{maftools} package.
\strong{NOTE}: all reference signatures are generated from gold standard tool:
SigProfiler.
}
\examples{
# Load mutational signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))

s1 <- get_sig_similarity(sig2, Ref = sig2)
s1

s2 <- get_sig_similarity(sig2)
s2
s3 <- get_sig_similarity(sig2, sig_db = "SBS")
s3

# Set order for result similarity matrix
s4 <- get_sig_similarity(sig2, sig_db = "SBS", set_order = TRUE)
s4

## Remove some components
## in similarity calculation
s5 <- get_sig_similarity(sig2,
  Ref = sig2,
  pattern_to_rm = c("T[T>G]C", "T[T>G]G", "T[T>G]T")
)
s5

## Same to DBS and ID signatures
x1 <- get_sig_db("DBS_hg19")
x2 <- get_sig_db("DBS_hg38")
s6 <- get_sig_similarity(x1$db, x2$db)
s6
}
\references{
Alexandrov, Ludmil B., et al. "The repertoire of mutational signatures in human cancer." Nature 578.7793 (2020): 94-101.

Degasperi, Andrea, et al. "A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies." Nature cancer 1.2 (2020): 249-263.

Steele, Christopher D., et al. "Undifferentiated sarcomas develop through distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.

Nik-Zainal, Serena, et al. "Landscape of somatic mutations in 560 breast cancer whole-genome sequences." Nature 534.7605 (2016): 47-54.

Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." bioRxiv (2021).
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cn_ploidy.R
\name{get_cn_ploidy}
\alias{get_cn_ploidy}
\title{Get Ploidy from Absolute Copy Number Profile}
\usage{
get_cn_ploidy(data)
}
\arguments{
\item{data}{a \link{CopyNumber} object or a \code{data.frame} containing at least 'chromosome', 'start',
'end', 'segVal' these columns.}
}
\value{
a value or a \code{data.table}
}
\description{
Get Ploidy from Absolute Copy Number Profile
}
\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

df <- get_cn_ploidy(cn)
df
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_group_enrichment.R
\name{show_group_enrichment}
\alias{show_group_enrichment}
\title{Show Group Enrichment Result}
\usage{
show_group_enrichment(
  df_enrich,
  return_list = FALSE,
  scales = "free",
  add_text_annotation = TRUE,
  fill_by_p_value = TRUE,
  use_fdr = TRUE,
  cut_p_value = FALSE,
  cut_breaks = c(-Inf, -5, log10(0.05), -log10(0.05), 5, Inf),
  cut_labels = c("↓ 1e-5", "↓ 0.05", "non-significant", "↑ 0.05", "↑ 1e-5"),
  fill_scale = scale_fill_gradient2(low = "#08A76B", mid = "white", high = "red",
    midpoint = ifelse(fill_by_p_value, 0, 1)),
  cluster_row = FALSE,
  ...
)
}
\arguments{
\item{df_enrich}{result \code{data.frame} from \link{group_enrichment}.}

\item{return_list}{if \code{TRUE}, return a list of \code{ggplot} object so user
can combine multiple plots by other R packages like \code{patchwork}.}

\item{scales}{Should scales be fixed (\code{"fixed"}, the default),
free (\code{"free"}), or free in one dimension (\code{"free_x"},
\code{"free_y"})?}

\item{add_text_annotation}{if \code{TRUE}, add text annotation in box.
When show p value with filled color, the text indicates relative change;
when show relative change with filled color, the text indicates p value.}

\item{fill_by_p_value}{if \code{TRUE}, show log10 based p values with filled color.
The +/- of p values indicates change direction.}

\item{use_fdr}{if \code{TRUE}, show FDR values instead of raw p-values.}

\item{cut_p_value}{if \code{TRUE}, cut p values into 5 regions for better visualization.
Only works when \code{fill_by_p_value = TRUE}.}

\item{cut_breaks}{when \code{cut_p_value} is \code{TRUE}, this option set the (log10 based) breaks.}

\item{cut_labels}{when \code{cut_p_value} is \code{TRUE}, this option set the labels.}

\item{fill_scale}{a \code{Scale} object generated by \code{ggplot2} package to
set color for continuous values.}

\item{cluster_row}{if \code{TRUE}, cluster rows with Hierarchical Clustering ('complete' method).}

\item{...}{other parameters passing to \link[ggplot2:facet_wrap]{ggplot2::facet_wrap}, only used
when \code{return_list} is \code{FALSE}.}
}
\value{
a (list of) \code{ggplot} object.
}
\description{
See \link{group_enrichment} for examples.
NOTE the box fill and the box text have different meanings.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulation}
\alias{simulation}
\alias{simulate_signature}
\alias{simulate_catalogue}
\alias{simulate_catalogue_matrix}
\title{Simulation Analysis}
\usage{
simulate_signature(x, weights = NULL)

simulate_catalogue(x, n, weights = NULL)

simulate_catalogue_matrix(x)
}
\arguments{
\item{x}{a numeric vector representing a signature/catalog or matrix with rows representing
signatures/samples and columns representing components.}

\item{weights}{a numeric vector for weights.}

\item{n}{an integer indicating mutation number to be generated in a catalog.}
}
\value{
a \code{matrix}.
}
\description{
\itemize{
\item \code{simulate_signature()} - Simulate signatures from signature pool.
\item \code{simulate_catalogue()} - Simulate catalogs from signature/catalog pool.
\item \code{simulate_catalogue_matrix()} - Simulate a bootstrapped catalog matrix.
}
}
\examples{
# Generate a catalog
set.seed(1234)
catalog <- as.integer(table(sample(1:96, 1000, replace = TRUE)))
names(catalog) <- paste0("comp", 1:96)
# Generate a signature
sig <- catalog / sum(catalog)

# Simulate catalogs
x1 <- simulate_catalogue(catalog, 10) # 10 mutations
x1
x2 <- simulate_catalogue(catalog, 100) # 100 mutations
x2
x3 <- simulate_catalogue(catalog, 1000) # 1000 mutations
x3
# Similar with a signature
x4 <- simulate_catalogue(sig, 10) # 10 mutations
x4

# Load SBS signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
s <- t(sig2$Signature.norm)
# Generate a signature from multiple signatures/catalogs
s1 <- simulate_signature(s)
s1
s2 <- simulate_signature(s, weights = 1:3)
s2
# Generate a catalog from multiple signatures/catalogs
c1 <- simulate_catalogue(s, 100, weights = 1:3)
c1
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_fit_bootstrap.R
\name{sig_fit_bootstrap}
\alias{sig_fit_bootstrap}
\title{Obtain Bootstrap Distribution of Signature Exposures of a Certain Tumor Sample}
\usage{
sig_fit_bootstrap(
  catalog,
  sig,
  n = 100L,
  sig_index = NULL,
  sig_db = "legacy",
  db_type = c("", "human-exome", "human-genome"),
  show_index = TRUE,
  method = c("QP", "NNLS", "SA"),
  auto_reduce = FALSE,
  SA_not_bootstrap = FALSE,
  type = c("absolute", "relative"),
  rel_threshold = 0,
  mode = c("SBS", "DBS", "ID", "copynumber"),
  find_suboptimal = FALSE,
  suboptimal_ref_error = NULL,
  suboptimal_factor = 1.05,
  ...
)
}
\arguments{
\item{catalog}{a named numeric vector or a numeric matrix with dimension Nx1.
N is the number of component, 1 is the sample.}

\item{sig}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix/\code{data.frame} with row representing components (motifs) and
column representing signatures.}

\item{n}{the number of bootstrap replicates.}

\item{sig_index}{a vector for signature index. "ALL" for all signatures.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}

\item{db_type}{only used when \code{sig_db} is enabled.
"" for keeping default, "human-exome" for transforming to exome frequency of component,
and "human-genome" for transforming to whole genome frequency of component.
Currently only works for 'SBS'.}

\item{show_index}{if \code{TRUE}, show valid indices.}

\item{method}{method to solve the minimazation problem.
'NNLS' for non-negative least square; 'QP' for quadratic programming; 'SA' for simulated annealing.}

\item{auto_reduce}{if \code{TRUE}, try reducing the input reference signatures to increase
the cosine similarity of reconstructed profile to observed profile.}

\item{SA_not_bootstrap}{if \code{TRUE}, directly run 'SA' multiple times with original input instead of
bootstrap samples.}

\item{type}{'absolute' for signature exposure and 'relative' for signature relative exposure.}

\item{rel_threshold}{numeric vector, a signature with relative exposure
lower than (equal is included, i.e. \code{<=}) this value will be set to 0
(both absolute exposure and relative exposure).
In this case, sum of signature contribution may not equal to 1.}

\item{mode}{signature type for plotting, now supports 'copynumber', 'SBS',
'DBS', 'ID' and 'RS' (genome rearrangement signature).}

\item{find_suboptimal}{logical, if \code{TRUE}, find suboptimal decomposition with
slightly higher error than the optimal solution by method 'SA'. This is useful
to explore hidden dependencies between signatures. More see reference.}

\item{suboptimal_ref_error}{baseline error used for finding suboptimal solution.
if it is \code{NULL}, then use 'SA' method to obtain the optimal error.}

\item{suboptimal_factor}{suboptimal factor to get suboptimal error, default is \code{1.05},
i.e., suboptimal error is \code{1.05} times baseline error.}

\item{...}{control parameters passing to argument \code{control} in \code{GenSA} function when use method 'SA'.}
}
\value{
a \code{list}
}
\description{
This can be used to obtain the confidence of signature exposures or search
the suboptimal decomposition solution.
}
\examples{
W <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
colnames(W) <- c("sig1", "sig2")
W <- apply(W, 2, function(x) x / sum(x))

H <- matrix(c(2, 5, 3, 6, 1, 9, 1, 2), ncol = 4)
colnames(H) <- paste0("samp", 1:4)

V <- W \%*\% H
V

if (requireNamespace("quadprog", quietly = TRUE)) {
  H_bootstrap <- sig_fit_bootstrap(V[, 1], W, n = 10, type = "absolute")
  ## Typically, you have to run many times to get close to the answer
  boxplot(t(H_bootstrap$expo))
  H[, 1]

  ## Return P values
  ## In practice, run times >= 100
  ## is recommended
  report_bootstrap_p_value(H_bootstrap)
  ## For multiple samples
  ## Input a list
  report_bootstrap_p_value(list(samp1 = H_bootstrap, samp2 = H_bootstrap))

  #   ## Find suboptimal decomposition
  #   H_suboptimal <- sig_fit_bootstrap(V[, 1], W,
  #     n = 10,
  #     type = "absolute",
  #     method = "SA",
  #     find_suboptimal = TRUE
  #   )
}
}
\references{
Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604
}
\seealso{
\link{report_bootstrap_p_value}, \link{sig_fit}, \link{sig_fit_bootstrap_batch}
}
\keyword{bootstrap}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chromsize.mm10}
\alias{chromsize.mm10}
\title{Chromosome Size of Genome Build mm10}
\format{
A data.frame
}
\source{
Generate from UCSC gold path \url{http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes}
}
\description{
Chromosome Size of Genome Build mm10
}
\examples{
data(chromsize.mm10)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_exposure.R
\name{get_sig_exposure}
\alias{get_sig_exposure}
\title{Get Signature Exposure from 'Signature' Object}
\usage{
get_sig_exposure(
  Signature,
  type = c("absolute", "relative"),
  rel_threshold = 0.01
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw exposure matrix with column representing samples (patients) and row
representing signatures.}

\item{type}{'absolute' for signature exposure and 'relative' for signature relative exposure.}

\item{rel_threshold}{only used when type is 'relative', relative exposure less
than (\code{<=}) this value will be set to 0 and thus all signature exposures
may not sum to 1. This is similar to this argument in \link{sig_fit}.}
}
\value{
a \code{data.table}
}
\description{
The expected number of mutations (or copy number segment records) with each signature was
determined after a scaling transformation V ~ WH = W'H' where W' = WU' and H' = UH.
The scaling matrix U is a KxK diagnal matrix (K is signature number, U' is the inverse of U)
with the element corresponding to the L1-norm of column vectors of W
(ie. the sum of the elements of the vector). As a result, the k-th row vector of the final
matrix H' represents the absolute exposure (activity) of the k-th process across samples
(e.g., for SBS, the estimated (or expected) number of mutations generated by the k-th process).
Of note, for copy number signatures, only components of feature CN was used for calculating H'.
}
\examples{
# Load mutational signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Get signature exposure
expo1 <- get_sig_exposure(sig2)
expo1
expo2 <- get_sig_exposure(sig2, type = "relative")
expo2
}
\references{
Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
Nature genetics 48.6 (2016): 600.
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_fit.R
\name{show_sig_fit}
\alias{show_sig_fit}
\title{Show Signature Fit Result}
\usage{
show_sig_fit(
  fit_result,
  samples = NULL,
  signatures = NULL,
  plot_fun = c("boxplot", "violin", "scatter"),
  palette = "aaas",
  title = NULL,
  xlab = FALSE,
  ylab = "Signature exposure",
  legend = "none",
  width = 0.3,
  outlier.shape = NA,
  add = "jitter",
  add.params = list(alpha = 0.3),
  ...
)
}
\arguments{
\item{fit_result}{result object from \link{sig_fit}.}

\item{samples}{samples to show, if \code{NULL}, all samples are used.}

\item{signatures}{signatures to show.}

\item{plot_fun}{set the plot function.}

\item{palette}{the color palette to be used for coloring or filling by groups.
Allowed values include "grey" for grey color palettes; brewer palettes e.g.
"RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and
scientific journal palettes from ggsci R package, e.g.: "npg", "aaas",
"lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".}

\item{title}{plot main title.}

\item{xlab}{character vector specifying x axis labels. Use xlab = FALSE to
hide xlab.}

\item{ylab}{character vector specifying y axis labels. Use ylab = FALSE to
hide ylab.}

\item{legend}{character specifying legend position. Allowed values are one of
c("top", "bottom", "left", "right", "none"). To remove the legend use
legend = "none". Legend position can be also specified using a numeric
vector c(x, y); see details section.}

\item{width}{numeric value between 0 and 1 specifying box width.}

\item{outlier.shape}{point shape of outlier. Default is 19. To hide outlier,
specify \code{outlier.shape = NA}. When jitter is added, then outliers will
be automatically hidden.}

\item{add}{character vector for adding another plot element (e.g.: dot plot or
error bars). Allowed values are one or the combination of: "none",
"dotplot", "jitter", "boxplot", "point", "mean", "mean_se", "mean_sd",
"mean_ci", "mean_range", "median", "median_iqr", "median_hilow",
"median_q1q3", "median_mad", "median_range"; see ?desc_statby for more
details.}

\item{add.params}{parameters (color, shape, size, fill, linetype) for the
argument 'add'; e.g.: add.params = list(color = "red").}

\item{...}{other arguments to be passed to
\code{\link[ggplot2]{geom_boxplot}}, \code{\link[ggpubr]{ggpar}} and
\code{\link[ggpubr]{facet}}.}
}
\value{
a \code{ggplot} object.
}
\description{
See \link{sig_fit} for examples.
}
\seealso{
\link{sig_fit}, \link{show_sig_bootstrap_exposure}, \link{sig_fit_bootstrap}, \link{sig_fit_bootstrap_batch}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chromsize.hg19}
\alias{chromsize.hg19}
\title{Chromosome Size of Genome Build hg19}
\format{
A data.frame
}
\source{
Generate from UCSC gold path
}
\description{
Chromosome Size of Genome Build hg19
}
\examples{
data(chromsize.hg19)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_exposure.R
\name{show_sig_exposure}
\alias{show_sig_exposure}
\title{Plot Signature Exposure}
\usage{
show_sig_exposure(
  Signature,
  sig_names = NULL,
  groups = NULL,
  grp_order = NULL,
  grp_size = NULL,
  cutoff = NULL,
  style = c("default", "cosmic"),
  palette = use_color_style(style),
  base_size = 12,
  font_scale = 1,
  rm_space = FALSE,
  rm_grid_line = TRUE,
  rm_panel_border = FALSE,
  hide_samps = TRUE,
  legend_position = "top"
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw \strong{absolute} exposure matrix with column representing samples (patients) and row
representing signatures (signature names must end with different digital numbers,
e.g. Sig1, Sig10, x12). If you named signatures with letters,
you can specify them by \code{sig_names} parameter.}

\item{sig_names}{set name of signatures, can be a character vector.}

\item{groups}{sample groups, default is \code{NULL}.}

\item{grp_order}{order of groups, default is \code{NULL}.}

\item{grp_size}{font size of groups.}

\item{cutoff}{a cutoff value to remove hyper-mutated samples.}

\item{style}{plot style, one of 'default' and 'cosmic', works when
parameter \code{set_gradient_color} is \code{FALSE}.}

\item{palette}{palette used to plot, default use a built-in palette
according to parameter \code{style}.}

\item{base_size}{overall font size.}

\item{font_scale}{a number used to set font scale.}

\item{rm_space}{default is \code{FALSE}. If \code{TRUE}, it will remove border color
and expand the bar width to 1. This is useful when the sample size is big.}

\item{rm_grid_line}{default is \code{FALSE}, if \code{TRUE}, remove grid lines of plot.}

\item{rm_panel_border}{default is \code{TRUE} for style 'cosmic',
remove panel border to keep plot tight.}

\item{hide_samps}{if \code{TRUE}, hide sample names.}

\item{legend_position}{position of legend, default is 'top'.}
}
\value{
a \code{ggplot} object
}
\description{
Currently support copy number signatures and mutational signatures.
}
\examples{
# Load mutational signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature exposure
p1 <- show_sig_exposure(sig2)
p1

# Load copy number signature
load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature exposure
p2 <- show_sig_exposure(sig)
p2
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_maf.R
\name{read_maf}
\alias{read_maf}
\title{Read MAF Files}
\usage{
read_maf(maf, verbose = TRUE)
}
\arguments{
\item{maf}{tab delimited MAF file. File can also be gz compressed. Required. Alternatively, you can also provide already read MAF file as a dataframe.}

\item{verbose}{TRUE logical. Default to be talkative and prints summary.}
}
\description{
This function is a wrapper of \link[maftools:read.maf]{maftools::read.maf}.
Useless options in \link[maftools:read.maf]{maftools::read.maf} are dropped here.
You can also use \link[maftools:read.maf]{maftools::read.maf} to read the data.
All reference alleles and mutation alleles should be recorded in
positive strand format.
}
\examples{
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
if (!require("R.utils")) {
  message("Please install 'R.utils' package firstly")
} else {
  laml <- read_maf(maf = laml.maf)
  laml
}
}
\seealso{
\link{read_copynumber} for reading copy number data to \link{CopyNumber} object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_labels.R
\name{add_labels}
\alias{add_labels}
\title{Add Text Labels to a ggplot}
\usage{
add_labels(
  p,
  x,
  y,
  y_end = NULL,
  n_label = NULL,
  labels = NULL,
  revert_order = FALSE,
  font_size = 5,
  font_family = "serif",
  font_face = c("plain", "bold", "italic"),
  ...
)
}
\arguments{
\item{p}{a \code{ggplot}.}

\item{x}{position at x axis.}

\item{y}{position at y axis.}

\item{y_end}{end position of y axis when \code{n_label} is set.}

\item{n_label}{the number of label, when this is set,
the position of labels at y axis is auto-generated
according to \code{y} and \code{y_end}.}

\item{labels}{text labels or a \code{similarity} object from \link{get_sig_similarity}.}

\item{revert_order}{if \code{TRUE}, revert label order.}

\item{font_size}{font size.}

\item{font_family}{font family.}

\item{font_face}{font face.}

\item{...}{other parameters passing to \link[ggplot2:annotate]{ggplot2::annotate}.}
}
\value{
a \code{ggplot} object.
}
\description{
Add text labels to a ggplot object, such as the result
from \link{show_sig_profile}.
}
\examples{
# Load mutational signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p <- show_sig_profile(sig2, mode = "SBS")

# Method 1
p1 <- add_labels(p,
  x = 0.75, y = 0.3, y_end = 0.9, n_label = 3,
  labels = paste0("text", 1:3)
)
p1

# Method 2
p2 <- add_labels(p,
  x = c(0.15, 0.6, 0.75), y = c(0.3, 0.6, 0.9),
  labels = paste0("text", 1:3)
)
p2

# Method 3
sim <- get_sig_similarity(sig2)
p3 <- add_labels(p,
  x = c(0.15, 0.6, 0.75), y = c(0.25, 0.55, 0.8),
  labels = sim, font_size = 2
)
p3
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{transcript.hg19}
\alias{transcript.hg19}
\title{Merged Transcript Location at Genome Build hg19}
\format{
A \code{data.table}
}
\source{
from GENCODE release v33.
}
\description{
Merged Transcript Location at Genome Build hg19
}
\examples{
data(transcript.hg19)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cytobands.mm9}
\alias{cytobands.mm9}
\title{Location of Chromosome Cytobands at Genome Build mm9}
\format{
A data.frame
}
\source{
from UCSC \url{http://hgdownload.cse.ucsc.edu/goldenpath/mm9/database/cytoBand.txt.gz}
}
\description{
Location of Chromosome Cytobands at Genome Build mm9
}
\examples{
data(cytobands.mm9)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-tidy-eval.R
\name{tidyeval}
\alias{tidyeval}
\alias{expr}
\alias{enquo}
\alias{enquos}
\alias{sym}
\alias{syms}
\alias{.data}
\alias{.env}
\alias{:=}
\alias{as_name}
\alias{as_label}
\title{Tidy eval helpers}
\description{
\itemize{
\item \code{\link[rlang]{sym}()} creates a symbol from a string and
\code{\link[rlang:sym]{syms}()} creates a list of symbols from a
character vector.
\item \code{\link[rlang:nse-defuse]{enquo}()} and
\code{\link[rlang:nse-defuse]{enquos}()} delay the execution of one or
several function arguments. \code{enquo()} returns a single quoted
expression, which is like a blueprint for the delayed computation.
\code{enquos()} returns a list of such quoted expressions.
\item \code{\link[rlang:nse-defuse]{expr}()} quotes a new expression \emph{locally}. It
is mostly useful to build new expressions around arguments
captured with \code{\link[=enquo]{enquo()}} or \code{\link[=enquos]{enquos()}}:
\code{expr(mean(!!enquo(arg), na.rm = TRUE))}.
\item \code{\link[rlang]{as_name}()} transforms a quoted variable name
into a string. Supplying something else than a quoted variable
name is an error.

That's unlike \code{\link[rlang]{as_label}()} which also returns
a single string but supports any kind of R object as input,
including quoted function calls and vectors. Its purpose is to
summarise that object into a single label. That label is often
suitable as a default name.

If you don't know what a quoted expression contains (for instance
expressions captured with \code{enquo()} could be a variable
name, a call to a function, or an unquoted constant), then use
\code{as_label()}. If you know you have quoted a simple variable
name, or would like to enforce this, use \code{as_name()}.
}

To learn more about tidy eval and how to use these tools, visit
\url{https://tidyeval.tidyverse.org} and the
\href{https://adv-r.hadley.nz/metaprogramming.html}{Metaprogramming
section} of \href{https://adv-r.hadley.nz}{Advanced R}.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_unify_extract.R
\name{sig_unify_extract}
\alias{sig_unify_extract}
\title{An Unified Interface to Extract Signatures}
\usage{
sig_unify_extract(
  nmf_matrix,
  range = 2:5,
  nrun = 10,
  approach = c("bayes_nmf", "repeated_nmf", "bootstrap_nmf", "sigprofiler"),
  cores = 1L,
  ...
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{range}{signature number range, i.e. \code{2:5}.}

\item{nrun}{the number of iteration to be performed to extract each signature number.}

\item{approach}{approach name.
\itemize{
\item "repeated_nmf" - \link{sig_extract}
\item "bayes_nmf" - \link{sig_auto_extract}
\item "bootstrap_nmf" - \link{bp_extract_signatures}
\item "sigprofiler" - \link{sigprofiler}
}}

\item{cores}{number of cores used for computation.}

\item{...}{other parameters passing to signature extractor based
on the \code{approach} setting.}
}
\value{
Result dependent on the \code{approach} setting.
}
\description{
This function provides an unified interface to signature extractor
implemented in \strong{sigminer}. If you determine a specific \code{approach},
please also read the documentation of corresponding extractor.
See "Arguments" part.
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Extract signatures
# It is same as sig_extract(cn_tally_W$nmf_matrix, 2, nrun = 1)
res <- sig_unify_extract(cn_tally_W$nmf_matrix, 2,
  nrun = 1,
  approach = "repeated_nmf"
)
# Auto-extract signatures based on bayesian NMF
res2 <- sig_unify_extract(cn_tally_W$nmf_matrix,
  nrun = 1,
  approach = "bayes_nmf"
)
}
}
\seealso{
\link{sig_extract}, \link{sig_auto_extract}, \link{bp_extract_signatures},
\link{sigprofiler}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_fit.R
\name{sig_fit}
\alias{sig_fit}
\title{Fit Signature Exposures with Linear Combination Decomposition}
\usage{
sig_fit(
  catalogue_matrix,
  sig,
  sig_index = NULL,
  sig_db = c("legacy", "SBS", "DBS", "ID", "TSB", "SBS_Nik_lab", "RS_Nik_lab",
    "RS_BRCA560", "RS_USARC", "CNS_USARC", "CNS_TCGA", "SBS_hg19", "SBS_hg38", "SBS_mm9",
    "SBS_mm10", "DBS_hg19", "DBS_hg38", "DBS_mm9", "DBS_mm10", "SBS_Nik_lab_Organ",
    "RS_Nik_lab_Organ", "latest_SBS_GRCh37", "latest_DBS_GRCh37", "latest_ID_GRCh37",
    "latest_SBS_GRCh38", "latest_DBS_GRCh38", "latest_SBS_mm9", "latest_DBS_mm9",
    "latest_SBS_mm10", "latest_DBS_mm10", "latest_SBS_rn6", "latest_DBS_rn6"),
  db_type = c("", "human-exome", "human-genome"),
  show_index = TRUE,
  method = c("QP", "NNLS", "SA"),
  auto_reduce = FALSE,
  type = c("absolute", "relative"),
  return_class = c("matrix", "data.table"),
  return_error = FALSE,
  rel_threshold = 0,
  mode = c("SBS", "DBS", "ID", "copynumber"),
  true_catalog = NULL,
  ...
)
}
\arguments{
\item{catalogue_matrix}{a numeric matrix \code{V} with row representing components and
columns representing samples, typically you can get \code{nmf_matrix} from \code{sig_tally()} and
transpose it by \code{t()}.}

\item{sig}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix/\code{data.frame} with row representing components (motifs) and
column representing signatures.}

\item{sig_index}{a vector for signature index. "ALL" for all signatures.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}

\item{db_type}{only used when \code{sig_db} is enabled.
"" for keeping default, "human-exome" for transforming to exome frequency of component,
and "human-genome" for transforming to whole genome frequency of component.
Currently only works for 'SBS'.}

\item{show_index}{if \code{TRUE}, show valid indices.}

\item{method}{method to solve the minimazation problem.
'NNLS' for non-negative least square; 'QP' for quadratic programming; 'SA' for simulated annealing.}

\item{auto_reduce}{if \code{TRUE}, try reducing the input reference signatures to increase
the cosine similarity of reconstructed profile to observed profile.}

\item{type}{'absolute' for signature exposure and 'relative' for signature relative exposure.}

\item{return_class}{string, 'matrix' or 'data.table'.}

\item{return_error}{if \code{TRUE}, also return sample error (Frobenius norm) and cosine
similarity between observed sample profile (asa. spectrum) and reconstructed profile. NOTE:
it is better to obtain the error when the type is 'absolute', because the error is
affected by relative exposure accuracy.}

\item{rel_threshold}{numeric vector, a signature with relative exposure
lower than (equal is included, i.e. \code{<=}) this value will be set to 0
(both absolute exposure and relative exposure).
In this case, sum of signature contribution may not equal to 1.}

\item{mode}{signature type for plotting, now supports 'copynumber', 'SBS',
'DBS', 'ID' and 'RS' (genome rearrangement signature).}

\item{true_catalog}{used by \link{sig_fit_bootstrap}, user never use it.}

\item{...}{control parameters passing to argument \code{control} in \code{GenSA} function when use method 'SA'.}
}
\value{
The exposure result either in \code{matrix} or \code{data.table} format.
If \code{return_error} set \code{TRUE}, a \code{list} is returned.
}
\description{
The function performs a signatures decomposition of a given mutational
catalogue \code{V} with known signatures \code{W} by solving the minimization problem
\verb{min(||W*H - V||)} where W and V are known.
}
\details{
The method 'NNLS' solves the minimization problem with nonnegative least-squares constraints.
The method 'QP' and 'SA' are modified from SignatureEstimation package.
See references for details.
Of note, when fitting exposures for copy number signatures, only components of
feature CN is used.
}
\examples{
W <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
colnames(W) <- c("sig1", "sig2")
W <- apply(W, 2, function(x) x / sum(x))

H <- matrix(c(2, 5, 3, 6, 1, 9, 1, 2), ncol = 4)
colnames(H) <- paste0("samp", 1:4)

V <- W \%*\% H
V

if (requireNamespace("quadprog", quietly = TRUE)) {
  H_infer <- sig_fit(V, W, method = "QP")
  H_infer
  H

  H_dt <- sig_fit(V, W, method = "QP", auto_reduce = TRUE, return_class = "data.table")
  H_dt

  ## Show results
  show_sig_fit(H_infer)
  show_sig_fit(H_dt)

  ## Get clusters/groups
  H_dt_rel <- sig_fit(V, W, return_class = "data.table", type = "relative")
  z <- get_groups(H_dt_rel, method = "k-means")
  show_groups(z)
}

# if (requireNamespace("GenSA", quietly = TRUE)) {
#   H_infer <- sig_fit(V, W, method = "SA")
#   H_infer
#   H
#
#   H_dt <- sig_fit(V, W, method = "SA", return_class = "data.table")
#   H_dt
#
#   ## Modify arguments to method
#   sig_fit(V, W, method = "SA", maxit = 10, temperature = 100)
#
#   ## Show results
#   show_sig_fit(H_infer)
#   show_sig_fit(H_dt)
# }
}
\references{
Daniel Huebschmann, Zuguang Gu and Matthias Schlesner (2019). YAPSA: Yet Another Package for Signature Analysis. R package version 1.12.0.

Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604

Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
Nature genetics 48.6 (2016): 600.
}
\seealso{
\link{sig_extract}, \link{sig_auto_extract}, \link{sig_fit_bootstrap}, \link{sig_fit_bootstrap_batch}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{output_bootstrap}
\alias{output_bootstrap}
\title{Output Signature Bootstrap Fitting Results}
\usage{
output_bootstrap(x, result_dir, mut_type = "SBS", sig_db = mut_type)
}
\arguments{
\item{x}{result from \link{sig_fit_bootstrap_batch}.}

\item{result_dir}{a result directory.}

\item{mut_type}{one of 'SBS', 'DBS', 'ID' or 'CN'.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}
}
\value{
Nothing.
}
\description{
Output Signature Bootstrap Fitting Results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/best_practice.R
\name{bp}
\alias{bp}
\alias{bp_extract_signatures}
\alias{bp_extract_signatures_iter}
\alias{bp_cluster_iter_list}
\alias{bp_get_clustered_sigs}
\alias{bp_get_sig_obj}
\alias{bp_get_stats}
\alias{bp_get_rank_score}
\alias{bp_show_survey2}
\alias{bp_show_survey}
\alias{bp_attribute_activity}
\title{A Best Practice for Signature Extraction and Exposure (Activity) Attribution}
\usage{
bp_extract_signatures(
  nmf_matrix,
  range = 2:5,
  n_bootstrap = 20L,
  n_nmf_run = 50,
  RTOL = 0.001,
  min_contribution = 0,
  cores = min(4L, future::availableCores()),
  cores_solution = min(cores, length(range)),
  seed = 123456L,
  handle_hyper_mutation = TRUE,
  report_integer_exposure = FALSE,
  only_core_stats = nrow(nmf_matrix) > 100,
  cache_dir = file.path(tempdir(), "sigminer_bp"),
  keep_cache = FALSE,
  pynmf = FALSE,
  use_conda = TRUE,
  py_path = "/Users/wsx/anaconda3/bin/python"
)

bp_extract_signatures_iter(
  nmf_matrix,
  range = 2:5,
  sim_threshold = 0.95,
  max_iter = 10L,
  n_bootstrap = 20L,
  n_nmf_run = 50,
  RTOL = 0.001,
  min_contribution = 0,
  cores = min(4L, future::availableCores()),
  cores_solution = min(cores, length(range)),
  seed = 123456L,
  handle_hyper_mutation = TRUE,
  report_integer_exposure = FALSE,
  only_core_stats = nrow(nmf_matrix) > 100,
  cache_dir = file.path(tempdir(), "sigminer_bp"),
  keep_cache = FALSE,
  pynmf = FALSE,
  use_conda = FALSE,
  py_path = "/Users/wsx/anaconda3/bin/python"
)

bp_cluster_iter_list(x, k = NULL, include_final_iteration = TRUE)

bp_get_clustered_sigs(SigClusters, cluster_label)

bp_get_sig_obj(obj, signum = NULL)

bp_get_stats(obj)

bp_get_rank_score(obj)

bp_show_survey2(
  obj,
  x = "signature_number",
  left_y = "silhouette",
  right_y = "L2_error",
  left_name = left_y,
  right_name = right_y,
  left_color = "black",
  right_color = "red",
  left_shape = 16,
  right_shape = 18,
  shape_size = 4,
  highlight = NULL
)

bp_show_survey(
  obj,
  add_score = FALSE,
  scales = c("free_y", "free"),
  fixed_ratio = TRUE
)

bp_attribute_activity(
  input,
  sample_class = NULL,
  nmf_matrix = NULL,
  method = c("bt", "stepwise"),
  bt_use_prop = FALSE,
  return_class = c("matrix", "data.table"),
  use_parallel = FALSE,
  cache_dir = file.path(tempdir(), "sigminer_attribute_activity"),
  keep_cache = FALSE
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{range}{a \code{numeric} vector containing the ranks of factorization to try. Note that duplicates are removed
and values are sorted in increasing order. The results are notably returned in this order.}

\item{n_bootstrap}{number of bootstrapped (resampling) catalogs used.
When it is \code{0}, the original (input) mutation catalog is used for NMF decomposition,
this is not recommended, just for testing, user should not set it to \code{0}.}

\item{n_nmf_run}{number of NMF runs for each bootstrapped or original catalog.
At default, in total n_bootstrap x n_nmf_run (i.e. 1000) NMF runs are used
for the task.}

\item{RTOL}{a threshold proposed by Nature Cancer paper to control how to
filter solutions of NMF. Default is \verb{0.1\%} (from reference #2),
only NMF solutions with KLD (KL deviance) <= \verb{100.1\%} minimal KLD are kept.}

\item{min_contribution}{a component contribution threshold to filer out small
contributed components.}

\item{cores}{number of cpu cores to run NMF.}

\item{cores_solution}{cores for processing solutions, default is equal to argument \code{cores}.}

\item{seed}{a random seed to make reproducible result.}

\item{handle_hyper_mutation}{default is \code{TRUE}, handle hyper-mutant samples.}

\item{report_integer_exposure}{if \code{TRUE}, report integer signature
exposure by bootstrapping technique.}

\item{only_core_stats}{if \code{TRUE}, only calculate the core stats for signatures and samples.}

\item{cache_dir}{a directory for keep temp result files.}

\item{keep_cache}{if \code{TRUE}, keep cache results.}

\item{pynmf}{if \code{TRUE}, use Python NMF driver \href{http://nimfa.biolab.si/index.html}{Nimfa}.
The seed currently is not used by this implementation, so the only way to reproduce
your result is setting \code{keep_cache = TRUE}.}

\item{use_conda}{if \code{TRUE}, create an independent conda environment to run NMF.}

\item{py_path}{path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'. In my
test, it is more stable than \code{use_conda=TRUE}. You can install the Nimfa package by yourself
or set \code{use_conda} to \code{TRUE} to install required Python environment, and then set this option.}

\item{sim_threshold}{a similarity threshold for selecting samples to auto-rerun
the extraction procedure (i.e. \code{bp_extract_signatures()}), default is \code{0.95}.}

\item{max_iter}{the maximum iteration size, default is 10, i.e., at most run
the extraction procedure 10 times.}

\item{x}{result from \code{\link[=bp_extract_signatures_iter]{bp_extract_signatures_iter()}} or a list of
\code{Signature} objects.}

\item{k}{an integer sequence specifying the cluster number to get silhouette.}

\item{include_final_iteration}{if \code{FALSE}, exclude final iteration result
from clustering for input from \code{\link[=bp_extract_signatures_iter]{bp_extract_signatures_iter()}}, not applied
if input is a list of \code{Signature} objects.}

\item{SigClusters}{result from \code{\link[=bp_cluster_iter_list]{bp_cluster_iter_list()}}.}

\item{cluster_label}{cluster labels for a specified cluster number, obtain it
from \code{SigClusters$sil_df}.}

\item{obj}{a \code{ExtractionResult} object from \code{\link[=bp_extract_signatures]{bp_extract_signatures()}}.}

\item{signum}{a integer vector to extract the corresponding \code{Signature} object(s).
If it is \code{NULL} (default), all will be returned.}

\item{left_y}{column name for left y axis.}

\item{right_y}{column name for right y axis.}

\item{left_name}{label name for left y axis.}

\item{right_name}{label name for right y axis.}

\item{left_color}{color for left axis.}

\item{right_color}{color for right axis.}

\item{left_shape}{shape setting.}

\item{right_shape}{shape setting.}

\item{shape_size}{shape setting.}

\item{highlight}{a \code{integer} to highlight a \code{x}.}

\item{add_score}{if \code{FALSE}, don't show score and label optimal points by
rank score.}

\item{scales}{one of "free_y" (default) and "free" to control the scales
of plot facet.}

\item{fixed_ratio}{if \code{TRUE} (default), make the x/y axis ratio fixed.}

\item{input}{result from \code{\link[=bp_extract_signatures]{bp_extract_signatures()}} or a Signature object.}

\item{sample_class}{a named string vector whose names are sample names
and values are class labels (i.e. cancer subtype). If it is \code{NULL} (the default),
treat all samples as one group.}

\item{method}{one of 'bt' (use bootstrap exposure median, from reference #2,
\strong{the most recommended way in my personal view}) or stepwise'
(stepwise reduce and update signatures then do signature fitting
with last signature sets, from reference #2, the result tends to assign
the contribution of removed signatures to the remaining signatures,
\strong{maybe I misunderstand the paper method? PAY ATTENTION}).}

\item{bt_use_prop}{this parameter is only used for \code{bt} method to reset
low contributing signature activity (relative activity \verb{<0.01}). If \code{TRUE},
use empirical P value calculation way (i.e. proportion, used by reference \verb{#2}),
otherwise a \code{t.test} is applied.}

\item{return_class}{string, 'matrix' or 'data.table'.}

\item{use_parallel}{if \code{TRUE}, use parallel computation based on \strong{furrr} package.
It can also be an integer for specifying cores.}
}
\value{
It depends on the called function.
}
\description{
These functions are combined to provide a best practice for optimally
identifying mutational signatures and attributing their activities (exposures)
in tumor samples. They are listed in order to use.
\itemize{
\item \code{bp_extract_signatures()} for extracting signatures.
\item \code{bp_show_survey()} for showing measures change under different
signature numbers to help user select optimal signature number.
At default, an aggregated score (named score) is generated to
suggest the best solution.
\item \code{bp_show_survey2()} for showing simplified signature number survey like
\code{\link[=show_sig_number_survey]{show_sig_number_survey()}}.
\item \code{bp_get_sig_obj()} for get a (list of) \code{Signature} object which is common
used in \strong{sigminer} for analysis and visualization.
\item \code{bp_attribute_activity()} for optimizing signature activities (exposures).
NOTE: the activities from extraction step may be better!
You can also use \link{sig_extract} to get optimal NMF result from multiple NMF runs.
Besides, you can use \link{sig_fit} to quantify exposures based on signatures extracted
from \code{bp_extract_signatures()}.
\item \code{bp_extract_signatures_iter()} for extracting signature in a iteration way.
\item \code{bp_cluster_iter_list()} for clustering (\code{hclust} with average linkage)
iterated signatures to help collapse
multiple signatures into one. The result cluster can be visualized by
\code{plot()} or \code{factoextra::fviz_dend()}.
\item \code{bp_get_clustered_sigs()} for getting clustered (grouped) mean signatures from signature clusters.
\item Extra: \code{bp_get_stats}() for obtaining stats for signatures and samples of a solution.
These stats are aggregated (averaged) as the stats for a solution
(specific signature number).
\item Extra: \code{bp_get_rank_score()} for obtaining rank score for all signature numbers.
}
}
\details{
The signature extraction approach is adopted from reference #1, #2, and
the whole best practice is adopted from the pipeline used by reference #3.
I implement the whole procedure with R code based on the method description
of papers. The code is well organized, tested and documented so user will
find it pretty simple and useful. Besides, the structure of the results is
very clear to see and also visualize like other approaches provided by \strong{sigminer}.
}
\section{Measure Explanation in Survey Plot}{

The survey plot provides a pretty good way to facilitate the signature number
selection. A \code{score} measure is calculated as the weighted mean of selected
measures and visualized as the first sub-plot. The optimal number is highlighted
with red color dot and the best values for each measures are also
highlighted with orange color dots. The detail of 6 measures shown in plot are
explained as below.
\itemize{
\item \code{score} - an aggregated score based on rank scores from selected measures below.
The higher, the better. When two signature numbers have the same score,
the larger signature number is preferred (this is a rare situation, you
have to double check other measures).
\item \code{silhouette} - the average silhouette width for signatures, also named as ASW in reference #2.
The signature number with silhouette decreases sharply is preferred.
\item \code{distance} - the average sample reconstructed cosine distance, the lower value is better.
\item \code{error} - the average sample reconstructed error calculated with L2 formula
(i.e. L2 error). This lower value is better. This measure represents a
similar concept like \code{distance} above, they are all used to quantify how well
sample mutation profiles can be reconstructed from signatures, but \code{distance}
cares the whole mutation profile similarity while \code{error} here cares value difference.
\item \verb{pos cor} - the average positive signature exposure correlation coefficient.
The lower value is better. This measure is constructed based on my understanding
about signatures: mutational signatures are typically treated as independent
recurrent patterns, so their activities are less correlated.
\item \code{similarity} - the average similarity within in a signature cluster.
Like \code{silhouette}, the point decreases sharply is preferred.
In the practice, results from multiple NMF runs are clustered
with "clustering with match" algorithm proposed by reference #2. This value
indicates if the signature profiles extracted from different NMF runs are similar.
}
}

\examples{
data("simulated_catalogs")
\donttest{
# Here I reduce the values for n_bootstrap and n_nmf_run
# for reducing the run time.
# In practice, you should keep default or increase the values
# for better estimation.
#
# The input data here is simulated from 10 mutational signatures

# e1 <- bp_extract_signatures(
#   t(simulated_catalogs$set1),
#   range = 8:12,
#   n_bootstrap = 5,
#   n_nmf_run = 10
# )
#
# To avoid computation in examples,
# Here just load the result
# (e1$signature and e1$exposure set to NA to reduce package size)
load(system.file("extdata", "e1.RData", package = "sigminer"))


# See the survey for different signature numbers
# The suggested solution is marked as red dot
# with highest integrated score.
p1 <- bp_show_survey(e1)
p1
# You can also exclude plotting and highlighting the score
p2 <- bp_show_survey(e1, add_score = FALSE)
p2

# You can also plot a simplified version
p3 <- bp_show_survey2(e1, highlight = 10)
p3

# Obtain the suggested solution from extraction result
obj_suggested <- bp_get_sig_obj(e1, e1$suggested)
obj_suggested
# If you think the suggested signature number is not right
# Just pick up the solution you want
obj_s8 <- bp_get_sig_obj(e1, 8)

# Track the reconstructed profile similarity
rec_sim <- get_sig_rec_similarity(obj_s8, t(simulated_catalogs$set1))
rec_sim

# After extraction, you can assign the signatures
# to reference COSMIC signatures
# More see ?get_sig_similarity
sim <- get_sig_similarity(obj_suggested)
# Visualize the match result
if (require(pheatmap)) {
  pheatmap::pheatmap(sim$similarity)
}

# You already got the activities of signatures
# in obj_suggested, however, you can still
# try to optimize the result.
# NOTE: the optimization step may not truly optimize the result!
expo <- bp_attribute_activity(e1, return_class = "data.table")
expo$abs_activity
}

\dontrun{
# Iterative extraction:
# This procedure will rerun extraction step
# for those samples with reconstructed catalog similarity
# lower than a threshold (default is 0.95)
e2 <- bp_extract_signatures_iter(
  t(simulated_catalogs$set1),
  range = 9:11,
  n_bootstrap = 5,
  n_nmf_run = 5,
  sim_threshold = 0.99
)
e2
# When the procedure run multiple rounds
# you can cluster the signatures from different rounds by
# the following command
# bp_cluster_iter_list(e2)

## Extra utilities
rank_score <- bp_get_rank_score(e1)
rank_score
stats <- bp_get_stats(e2$iter1)
# Get the mean reconstructed similarity
1 - stats$stats_sample$cosine_distance_mean
}
}
\references{
Alexandrov, Ludmil B., et al. "Deciphering signatures of mutational processes operative in human cancer." Cell reports 3.1 (2013): 246-259.

Degasperi, Andrea, et al. "A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies." Nature cancer 1.2 (2020): 249-263.

Alexandrov, Ludmil B., et al. “The repertoire of mutational signatures in human cancer.” Nature 578.7793 (2020): 94-101.
}
\seealso{
See \link{sig_estimate}, \link{sig_extract}, \link{sig_auto_extract},
\link{sigprofiler_extract} for other approaches.
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{chromsize.mm9}
\alias{chromsize.mm9}
\title{Chromosome Size of Genome Build mm9}
\format{
A data.frame
}
\source{
Generate from UCSC gold path \url{http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes}
}
\description{
Chromosome Size of Genome Build mm9
}
\examples{
data(chromsize.mm9)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{use_color_style}
\alias{use_color_style}
\title{Set Color Style for Plotting}
\usage{
use_color_style(
  style,
  mode = c("SBS", "copynumber", "DBS", "ID", "RS"),
  method = "Wang"
)
}
\arguments{
\item{style}{one of 'default' and 'cosmic'.}

\item{mode}{only used when the \code{style} is 'cosmic', can be one of
"SBS", "copynumber", "DBS", "ID".}

\item{method}{used to set a more custom palette for different methods.}
}
\value{
color values.
}
\description{
Set Color Style for Plotting
}
\examples{
use_color_style("default")
use_color_style("cosmic")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_cancer_type_index.R
\name{get_sig_cancer_type_index}
\alias{get_sig_cancer_type_index}
\title{Obtain Signature Index for Cancer Types}
\usage{
get_sig_cancer_type_index(
  sig_type = c("legacy", "SBS", "DBS", "ID"),
  seq_type = c("WGS", "WES"),
  source = c("PCAWG", "TCGA", "nonPCAWG"),
  keyword = NULL
)
}
\arguments{
\item{sig_type}{signature type.}

\item{seq_type}{sequencing type.}

\item{source}{data source.}

\item{keyword}{keyword to search in the signature index database.}
}
\value{
a \code{list}.
}
\description{
Obtain Signature Index for Cancer Types
}
\examples{
l1 <- get_sig_cancer_type_index()
l2 <- get_sig_cancer_type_index(sig_type = "SBS")
l3 <- get_sig_cancer_type_index(sig_type = "DBS", source = "PCAWG", seq_type = "WGS")
l4 <- get_sig_cancer_type_index(sig_type = "ID")
l5 <- get_sig_cancer_type_index(keyword = "breast")
l1
l2
l3
l4
l5
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_copynumber_seqz.R
\name{read_copynumber_seqz}
\alias{read_copynumber_seqz}
\title{Read Absolute Copy Number Profile from Sequenza Result Directory}
\usage{
read_copynumber_seqz(target_dir, return_df = FALSE, ...)
}
\arguments{
\item{target_dir}{a directory path.}

\item{return_df}{if \code{TRUE}, return a \code{data.frame} directly, otherwise return a
\link{CopyNumber} object.}

\item{...}{other parameters passing to \code{\link[=read_copynumber]{read_copynumber()}}.}
}
\value{
a \code{data.frame} or a \code{CopyNumber} object.
}
\description{
Read Absolute Copy Number Profile from Sequenza Result Directory
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_convert.R
\name{sig_convert}
\alias{sig_convert}
\title{Convert Signatures between different Genomic Distribution of Components}
\usage{
sig_convert(sig, from = "human-genome", to = "human-exome")
}
\arguments{
\item{sig}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix/\code{data.frame} with row representing components (motifs) and
column representing signatures.}

\item{from}{either one of "human-genome" and "human-exome" or an opportunity matrix
(repeated \code{n} columns with each row represents the total number of mutations for
a component, \code{n} is the number of signature).}

\item{to}{same as \code{from}.}
}
\value{
a \code{matrix}.
}
\description{
Converts signatures between two representations relative to different sets of mutational opportunities.
Currently, only SBS signature is supported.
}
\details{
The default opportunity matrix for "human-genome" and "human-exome" comes from COSMIC
signature database v2 and v3.
}
\examples{
# Load SBS signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Exome-relative to Genome-relative
sig_converted <- sig_convert(sig2,
  from = "human-exome",
  to = "human-genome"
)
sig_converted

show_sig_profile(sig2, style = "cosmic")
show_sig_profile(sig_converted, style = "cosmic")
}
\references{
\code{convert_signatures} function from sigfit package.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_vcf.R
\name{read_vcf}
\alias{read_vcf}
\title{Read VCF Files as MAF Object}
\usage{
read_vcf(
  vcfs,
  samples = NULL,
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  keep_only_pass = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{vcfs}{VCF file paths.}

\item{samples}{sample names for VCF files.}

\item{genome_build}{genome build version like "hg19".}

\item{keep_only_pass}{if \code{TRUE}, keep only 'PASS' mutation for analysis.}

\item{verbose}{if \code{TRUE}, print extra info.}
}
\value{
a \link{MAF}.
}
\description{
MAF file is more recommended. In this function, we will mimic
the MAF object from the key \code{c(1, 2, 4, 5, 7)} columns of VCF file.
}
\examples{
vcfs <- list.files(system.file("extdata", package = "sigminer"), "*.vcf", full.names = TRUE)
\donttest{
maf <- read_vcf(vcfs)
maf <- read_vcf(vcfs, keep_only_pass = TRUE)
}
}
\seealso{
\link{read_maf}, \link{read_copynumber}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{transcript.mm9}
\alias{transcript.mm9}
\title{Merged Transcript Location at Genome Build mm9}
\format{
A \code{data.table}
}
\source{
from UCSC \url{http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/transcriptome.txt.gz}
}
\description{
Merged Transcript Location at Genome Build mm9
}
\examples{
data(transcript.mm9)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigprofiler.R
\name{sigprofiler}
\alias{sigprofiler}
\alias{sigprofiler_extract}
\alias{sigprofiler_import}
\title{Extract Signatures with SigProfiler}
\usage{
sigprofiler_extract(
  nmf_matrix,
  output,
  range = 2:5,
  nrun = 10L,
  refit = FALSE,
  refit_plot = FALSE,
  is_exome = FALSE,
  init_method = c("nndsvd_min", "random", "alexandrov-lab-custom", "nndsvd", "nndsvda",
    "nndsvdar"),
  cores = -1L,
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  use_conda = FALSE,
  py_path = NULL,
  sigprofiler_version = "1.1.3"
)

sigprofiler_import(
  output,
  order_by_expo = FALSE,
  type = c("suggest", "refit", "all")
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{output}{output directory.}

\item{range}{signature number range, i.e. \code{2:5}.}

\item{nrun}{the number of iteration to be performed to extract each signature number.}

\item{refit}{if \code{TRUE}, then refit the denovo signatures with nnls. Same
meaning as \code{optimize} option in \link{sig_extract} or \link{sig_auto_extract}.}

\item{refit_plot}{if \code{TRUE}, SigProfiler will make
denovo to COSMIC sigantures decompostion plots. However, this may fail due
to some matrix cannot be identified by SigProfiler plot program.}

\item{is_exome}{if \code{TRUE}, the exomes will be extracted.}

\item{init_method}{the initialization algorithm for W and H matrix of NMF.
Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar', 'alexandrov-lab-custom'
and 'nndsvd_min'.}

\item{cores}{number of cores used for computation.}

\item{genome_build}{I think this option is useless when input is \code{matrix}, keep it
in case it is useful.}

\item{use_conda}{if \code{TRUE}, create an independent conda environment to run SigProfiler.}

\item{py_path}{path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'.}

\item{sigprofiler_version}{version of \code{SigProfilerExtractor}. If this
package is not installed, the specified package will be installed.
If this package is installed, this option is useless.}

\item{order_by_expo}{if \code{TRUE}, order the import signatures by their exposures, e.g. the signature
contributed the most exposure in all samples will be named as \code{Sig1}.}

\item{type}{one of 'suggest' (for suggested solution), 'refit' (for refit solution) or 'all' (for all solutions).}
}
\value{
For \code{sigprofiler_extract()}, returns nothing. See \code{output} directory.

For \code{sigprofiler_import()}, a \code{list} containing \code{Signature} object.
}
\description{
This function provides an interface to software SigProfiler.
More please see \url{https://github.com/AlexandrovLab/SigProfilerExtractor}.
Typically, a reference genome is not required because the input is a matrix (my understanding).
}
\examples{
if (FALSE) {
  load(system.file("extdata", "toy_copynumber_tally_W.RData",
    package = "sigminer", mustWork = TRUE
  ))

  reticulate::conda_list()

  sigprofiler_extract(cn_tally_W$nmf_matrix, "~/test/test_sigminer",
    use_conda = TRUE
  )

  sigprofiler_extract(cn_tally_W$nmf_matrix, "~/test/test_sigminer",
    use_conda = FALSE, py_path = "/Users/wsx/anaconda3/bin/python"
  )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_estimate.R
\name{sig_estimate}
\alias{sig_estimate}
\alias{show_sig_number_survey}
\alias{show_sig_number_survey2}
\title{Estimate Signature Number}
\usage{
sig_estimate(
  nmf_matrix,
  range = 2:5,
  nrun = 10,
  use_random = FALSE,
  method = "brunet",
  seed = 123456,
  cores = 1,
  keep_nmfObj = FALSE,
  save_plots = FALSE,
  plot_basename = file.path(tempdir(), "nmf"),
  what = "all",
  verbose = FALSE
)

show_sig_number_survey(
  object,
  x = "rank",
  left_y = "cophenetic",
  right_y = "rss",
  left_name = left_y,
  right_name = toupper(right_y),
  left_color = "black",
  right_color = "red",
  left_shape = 16,
  right_shape = 18,
  shape_size = 4,
  highlight = NULL
)

show_sig_number_survey2(
  x,
  y = NULL,
  what = c("all", "cophenetic", "rss", "residuals", "dispersion", "evar", "sparseness",
    "sparseness.basis", "sparseness.coef", "silhouette", "silhouette.coef",
    "silhouette.basis", "silhouette.consensus"),
  na.rm = FALSE,
  xlab = "Total signatures",
  ylab = "",
  main = "Signature number survey using NMF package"
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{range}{a \code{numeric} vector containing the ranks of factorization to try. Note that duplicates are removed
and values are sorted in increasing order. The results are notably returned in this order.}

\item{nrun}{a \code{numeric} giving the number of run to perform for each value in \code{range}, \code{nrun} set to 30~50 is
enough to achieve robust result.}

\item{use_random}{Should generate random data from input to test measurements. Default is \code{TRUE}.}

\item{method}{specification of the NMF algorithm. Use 'brunet' as default.
Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.}

\item{seed}{specification of the starting point or seeding method, which will compute a starting point,
usually using data from the target matrix in order to provide a good guess.}

\item{cores}{number of cpu cores to run NMF.}

\item{keep_nmfObj}{default is \code{FALSE}, if \code{TRUE}, keep NMF objects from runs, and the result may be huge.}

\item{save_plots}{if \code{TRUE}, save signature number survey plot to local machine.}

\item{plot_basename}{when save plots, set custom basename for file path.}

\item{what}{a character vector whose elements partially match one of the following item,
which correspond to the measures computed by \code{summary()} on each – multi-run – NMF result:
'all', 'cophenetic', 'rss', 'residuals', 'dispersion', 'evar', 'silhouette'
(and more specific \verb{*.coef}, \verb{*.basis}, \verb{*.consensus}), 'sparseness'
(and more specific \verb{*.coef}, \verb{*.basis}).
It specifies which measure must be plotted (what='all' plots all the measures).}

\item{verbose}{if \code{TRUE}, print extra message.}

\item{object}{a \code{Survey} object generated from \link{sig_estimate}, or
a \code{data.frame} contains at least rank columns and columns for
one measure.}

\item{x}{a \code{data.frame} or \code{NMF.rank} object obtained from \code{\link[=sig_estimate]{sig_estimate()}}.}

\item{left_y}{column name for left y axis.}

\item{right_y}{column name for right y axis.}

\item{left_name}{label name for left y axis.}

\item{right_name}{label name for right y axis.}

\item{left_color}{color for left axis.}

\item{right_color}{color for right axis.}

\item{left_shape, right_shape, shape_size}{shape setting.}

\item{highlight}{a \code{integer} to highlight a \code{x}.}

\item{y}{for random simulation,
a \code{data.frame} or \code{NMF.rank} object obtained from \code{\link[=sig_estimate]{sig_estimate()}}.}

\item{na.rm}{single logical that specifies if the rank
  for which the measures are NA values should be removed
  from the graph or not (default to \code{FALSE}).  This is
  useful when plotting results which include NAs due to
  error during the estimation process. See argument
  \code{stop} for \code{nmfEstimateRank}.}

\item{xlab}{x-axis label}

\item{ylab}{y-axis label}

\item{main}{main title}
}
\value{
\itemize{
\item sig_estimate: a \code{list} contains information of NMF run and rank survey.
}

\itemize{
\item show_sig_number_survey: a \code{ggplot} object
}

\itemize{
\item show_sig_number_survey2: a \code{ggplot} object
}
}
\description{
Use \strong{NMF} package to evaluate the optimal number of signatures.
This is used along with \link{sig_extract}.
Users should \code{library(NMF)} firstly. If NMF objects are returned,
the result can be further visualized by NMF plot methods like
\code{NMF::consensusmap()} and \code{NMF::basismap()}.

\code{sig_estimate()} shows comprehensive rank survey generated by
\strong{NMF} package, sometimes
it is hard to consider all measures. \code{show_sig_number_survey()} provides a
one or two y-axis visualization method to help users determine
the optimal signature number (showing both
stability ("cophenetic") and error (RSS) at default).
Users can also set custom measures to show.

\code{show_sig_number_survey2()} is modified from \strong{NMF} package to
better help users to explore survey of signature number.
}
\details{
The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
More custom features please directly use \link[NMF:nmfEstimateRank]{NMF::nmfEstimateRank}.
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))
library(NMF)
cn_estimate <- sig_estimate(cn_tally_W$nmf_matrix,
  cores = 1, nrun = 5,
  verbose = TRUE
)

p <- show_sig_number_survey2(cn_estimate$survey)
p

# Show two measures
show_sig_number_survey(cn_estimate)
# Show one measure
p1 <- show_sig_number_survey(cn_estimate, right_y = NULL)
p1
p2 <- add_h_arrow(p, x = 4.1, y = 0.953, label = "selected number")
p2

# Show data from a data.frame
p3 <- show_sig_number_survey(cn_estimate$survey)
p3
# Show other measures
head(cn_estimate$survey)
p4 <- show_sig_number_survey(cn_estimate$survey,
  right_y = "dispersion",
  right_name = "dispersion"
)
p4
p5 <- show_sig_number_survey(cn_estimate$survey,
  right_y = "evar",
  right_name = "evar"
)
p5
}
}
\references{
Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
}
\seealso{
\link{sig_extract} for extracting signatures using \strong{NMF} package, \link{sig_auto_extract} for
extracting signatures using automatic relevance determination technique.

\link{sig_estimate} for estimating signature number for \link{sig_extract},
\link{show_sig_number_survey2} for more visualization method.
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_extract.R
\name{sig_extract}
\alias{sig_extract}
\title{Extract Signatures through NMF}
\usage{
sig_extract(
  nmf_matrix,
  n_sig,
  nrun = 10,
  cores = 1,
  method = "brunet",
  optimize = FALSE,
  pynmf = FALSE,
  use_conda = TRUE,
  py_path = "/Users/wsx/anaconda3/bin/python",
  seed = 123456,
  ...
)
}
\arguments{
\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}

\item{n_sig}{number of signature. Please run \link{sig_estimate} to select a suitable value.}

\item{nrun}{a \code{numeric} giving the number of run to perform for each value in \code{range}, \code{nrun} set to 30~50 is
enough to achieve robust result.}

\item{cores}{number of cpu cores to run NMF.}

\item{method}{specification of the NMF algorithm. Use 'brunet' as default.
Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.}

\item{optimize}{if \code{TRUE}, then refit the denovo signatures with QP method, see \link{sig_fit}.}

\item{pynmf}{if \code{TRUE}, use Python NMF driver \href{http://nimfa.biolab.si/index.html}{Nimfa}.
The seed currently is not used by this implementation.}

\item{use_conda}{if \code{TRUE}, create an independent conda environment to run NMF.}

\item{py_path}{path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'. In my
test, it is more stable than \code{use_conda=TRUE}. You can install the Nimfa package by yourself
or set \code{use_conda} to \code{TRUE} to install required Python environment, and then set this option.}

\item{seed}{specification of the starting point or seeding method, which will compute a starting point,
usually using data from the target matrix in order to provide a good guess.}

\item{...}{other arguments passed to \code{\link[NMF:nmf]{NMF::nmf()}}.}
}
\value{
a \code{list} with \code{Signature} class.
}
\description{
Do NMF de-composition and then extract signatures.
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
res <- sig_extract(cn_tally_W$nmf_matrix, 2, nrun = 1)
}
}
\references{
Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.

Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
}
\seealso{
\link{sig_tally} for getting variation matrix,
\link{sig_estimate} for estimating signature number for \link{sig_extract}, \link{sig_auto_extract} for
extracting signatures using automatic relevance determination technique.
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{centromeres.hg38}
\alias{centromeres.hg38}
\title{Location of Centromeres at Genome Build hg38}
\format{
A data.frame
}
\source{
Generate from Genome Reference Consortium
}
\description{
Location of Centromeres at Genome Build hg38
}
\examples{
data(centromeres.hg38)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cytobands.mm10}
\alias{cytobands.mm10}
\title{Location of Chromosome Cytobands at Genome Build mm10}
\format{
A data.frame
}
\source{
from UCSC \url{http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz}
}
\description{
Location of Chromosome Cytobands at Genome Build mm10
}
\examples{
data(cytobands.mm10)
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
% Please edit documentation in R/get_groups.R
\name{get_groups}
\alias{get_groups}
\title{Get Sample Groups from Signature Decomposition Information}
\usage{
get_groups(
  Signature,
  method = c("consensus", "k-means", "exposure", "samples"),
  n_cluster = NULL,
  match_consensus = TRUE
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract}.
Now it can be used to relative exposure result in \code{data.table} format from \link{sig_fit}.}

\item{method}{grouping method, more see details, could be one of the following:
\itemize{
\item 'consensus' - returns the cluster membership based on the hierarchical clustering of the consensus matrix,
it can only be used for the result obtained by \code{\link[=sig_extract]{sig_extract()}} with multiple runs using \strong{NMF} package.
\item 'k-means' -  returns the clusters by k-means.
\item 'exposure' - assigns a sample into a group whose signature exposure
is dominant.
\item 'samples' - returns the cluster membership based on the contribution of signature to each sample,
it can only be used for the result obtained by \code{\link[=sig_extract]{sig_extract()}} using \strong{NMF} package.
}}

\item{n_cluster}{only used when the \code{method} is 'k-means'.}

\item{match_consensus}{only used when the \code{method} is 'consensus'.
If \code{TRUE}, the result will match order as shown in consensus map.}
}
\value{
a \code{data.table} object
}
\description{
One of key results from signature analysis is to cluster samples into different
groups. This function takes \code{Signature} object as input
and return the membership in each cluster.
}
\details{
Users may find there are bigger differences between using method 'samples' and 'exposure' but
they use a similar idear to find dominant signature, here goes the reason:

Method 'samples' using data directly from NMF decomposition, this means the two matrix
\code{W} (basis matrix or signature matrix) and \code{H} (coefficient matrix or exposure matrix) are
the results of NMF. For method 'exposure', it uses the signature exposure loading matrix.
In this situation, each signture represents a number of mutations (alterations)
about implementation please see source code of \code{\link[=sig_extract]{sig_extract()}} function.
}
\examples{
\donttest{
# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
sig <- sig_extract(cn_tally_W$nmf_matrix, 2,
  nrun = 10
)

# Methods 'consensus' and 'samples' are from NMF::predict()
g1 <- get_groups(sig, method = "consensus", match_consensus = TRUE)
g1
g2 <- get_groups(sig, method = "samples")
g2

# Use k-means clustering
g3 <- get_groups(sig, method = "k-means")
g3
}
}
\seealso{
\code{\link[NMF:predict]{NMF::predict()}}, \link{show_groups}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bayesian_result.R
\name{get_bayesian_result}
\alias{get_bayesian_result}
\title{Get Specified Bayesian NMF Result from Run}
\usage{
get_bayesian_result(run_info)
}
\arguments{
\item{run_info}{a \code{data.frame} with 1 row and two necessary columns \code{Run} and \code{file}.}
}
\value{
a \code{list}.
}
\description{
Sometimes, we may want to use or inspect specified run result from \link{sig_auto_extract}.
This function is designed for this purpose.
}
\examples{
load(system.file("extdata", "toy_copynumber_tally_W.RData",
  package = "sigminer", mustWork = TRUE
))

res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)

# All run info are stored in res$Raw$summary_run
# Obtain result of run 1
res_run1 <- get_bayesian_result(res$Raw$summary_run[1, ])
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_group_comparison.R
\name{get_group_comparison}
\alias{get_group_comparison}
\title{Get Comparison Result between Signature Groups}
\usage{
get_group_comparison(
  data,
  col_group,
  cols_to_compare,
  type = "ca",
  NAs = NA,
  verbose = FALSE
)
}
\arguments{
\item{data}{a \code{data.frame} containing signature groups and genotypes/phenotypes
(including categorical and continuous type data) want to analyze. User need to
construct this \code{data.frame} by him/herself.}

\item{col_group}{column name of signature groups.}

\item{cols_to_compare}{column names of genotypes/phenotypes want to summarize based on groups.}

\item{type}{a characater vector with length same as \code{cols_to_compare},
'ca' for categorical type and 'co' for continuous type.}

\item{NAs}{default is \code{NA}, filter \code{NA}s for categorical columns.
Otherwise a value (either length 1 or length same as \code{cols_to_compare}) fill \code{NA}s.}

\item{verbose}{if \code{TRUE}, print extra information.}
}
\value{
a \code{list} contains data, summary, p value etc..
}
\description{
Compare genotypes/phenotypes based on signature groups (samples are assigned to
several groups). For categorical
type, calculate fisher p value (using \link[stats:fisher.test]{stats::fisher.test}) and count table.
In larger than 2 by 2 tables, compute p-values by Monte Carlo simulation.
For continuous type, calculate anova p value (using \link[stats:aov]{stats::aov}),
summary table and Tukey Honest significant difference (using \link[stats:TukeyHSD]{stats::TukeyHSD}).
The result of this function can be plotted by \code{\link[=show_group_comparison]{show_group_comparison()}}.
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
  package = "sigminer", mustWork = TRUE
))

# Assign samples to clusters
groups <- get_groups(sig, method = "k-means")

set.seed(1234)

groups$prob <- rnorm(10)
groups$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(groups), replace = TRUE)

# Compare groups (filter NAs for categorical coloumns)
groups.cmp <- get_group_comparison(groups[, -1],
  col_group = "group",
  cols_to_compare = c("prob", "new_group"),
  type = c("co", "ca"), verbose = TRUE
)

# Compare groups (Set NAs of categorical columns to 'Rest')
groups.cmp2 <- get_group_comparison(groups[, -1],
  col_group = "group",
  cols_to_compare = c("prob", "new_group"),
  type = c("co", "ca"), NAs = "Rest", verbose = TRUE
)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_shannon_diversity_index.R
\name{get_shannon_diversity_index}
\alias{get_shannon_diversity_index}
\title{Get Shannon Diversity Index for Signatures}
\usage{
get_shannon_diversity_index(rel_expo, cutoff = 0.001)
}
\arguments{
\item{rel_expo}{a \code{data.frame} with numeric columns indicating
\strong{relative} signature exposures for each sample. Typically
this data can be obtained from \code{\link[=get_sig_exposure]{get_sig_exposure()}}.}

\item{cutoff}{a relative exposure cutoff for filtering signatures,
default is \verb{0.1\%}.}
}
\value{
a \code{data.frame}
}
\description{
\deqn{H = - \sum_{i=1}^n{p_i ln(p_i)}}
where \code{n} is the number
of signatures identified in the signature with exposure > \code{cutoff},
and \code{pi} is the normalized exposure of the ith signature with
exposure > \code{cutoff}. Exposures of signatures were normalized to
sum to \code{1}.
}
\examples{
# Load mutational signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Get signature exposure
rel_expo <- get_sig_exposure(sig2, type = "relative")
rel_expo
diversity_index <- get_shannon_diversity_index(rel_expo)
diversity_index
}
\references{
Steele, Christopher D., et al. "Undifferentiated sarcomas develop through distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_bootstrap.R
\name{show_sig_bootstrap}
\alias{show_sig_bootstrap}
\alias{show_sig_bootstrap_exposure}
\alias{show_sig_bootstrap_error}
\alias{show_sig_bootstrap_stability}
\title{Show Signature Bootstrap Analysis Results}
\usage{
show_sig_bootstrap_exposure(
  bt_result,
  sample = NULL,
  signatures = NULL,
  methods = "QP",
  plot_fun = c("boxplot", "violin"),
  agg_fun = c("mean", "median", "min", "max"),
  highlight = "auto",
  highlight_size = 4,
  palette = "aaas",
  title = NULL,
  xlab = FALSE,
  ylab = "Signature exposure",
  width = 0.3,
  dodge_width = 0.8,
  outlier.shape = NA,
  add = "jitter",
  add.params = list(alpha = 0.3),
  ...
)

show_sig_bootstrap_error(
  bt_result,
  sample = NULL,
  methods = "QP",
  plot_fun = c("boxplot", "violin"),
  agg_fun = c("mean", "median"),
  highlight = "auto",
  highlight_size = 4,
  palette = "aaas",
  title = NULL,
  xlab = FALSE,
  ylab = "Reconstruction error (L2 norm)",
  width = 0.3,
  dodge_width = 0.8,
  outlier.shape = NA,
  add = "jitter",
  add.params = list(alpha = 0.3),
  legend = "none",
  ...
)

show_sig_bootstrap_stability(
  bt_result,
  signatures = NULL,
  measure = c("RMSE", "CV", "MAE", "AbsDiff"),
  methods = "QP",
  plot_fun = c("boxplot", "violin"),
  palette = "aaas",
  title = NULL,
  xlab = FALSE,
  ylab = "Signature instability",
  width = 0.3,
  outlier.shape = NA,
  add = "jitter",
  add.params = list(alpha = 0.3),
  ...
)
}
\arguments{
\item{bt_result}{result object from \link{sig_fit_bootstrap_batch}.}

\item{sample}{a sample id.}

\item{signatures}{signatures to show.}

\item{methods}{a subset of \code{c("NNLS", "QP", "SA")}.}

\item{plot_fun}{set the plot function.}

\item{agg_fun}{set the aggregation function when \code{sample} is \code{NULL}.}

\item{highlight}{set the color for optimal solution. Default is "auto", which use the same color as
bootstrap results, you can set it to color like "red", "gold", etc.}

\item{highlight_size}{size for highlighting triangle, default is \code{4}.}

\item{palette}{the color palette to be used for coloring or filling by groups.
Allowed values include "grey" for grey color palettes; brewer palettes e.g.
"RdBu", "Blues", ...; or custom color palette e.g. c("blue", "red"); and
scientific journal palettes from ggsci R package, e.g.: "npg", "aaas",
"lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".}

\item{title}{plot main title.}

\item{xlab}{character vector specifying x axis labels. Use xlab = FALSE to
hide xlab.}

\item{ylab}{character vector specifying y axis labels. Use ylab = FALSE to
hide ylab.}

\item{width}{numeric value between 0 and 1 specifying box width.}

\item{dodge_width}{dodge width.}

\item{outlier.shape}{point shape of outlier. Default is 19. To hide outlier,
specify \code{outlier.shape = NA}. When jitter is added, then outliers will
be automatically hidden.}

\item{add}{character vector for adding another plot element (e.g.: dot plot or
error bars). Allowed values are one or the combination of: "none",
"dotplot", "jitter", "boxplot", "point", "mean", "mean_se", "mean_sd",
"mean_ci", "mean_range", "median", "median_iqr", "median_hilow",
"median_q1q3", "median_mad", "median_range"; see ?desc_statby for more
details.}

\item{add.params}{parameters (color, shape, size, fill, linetype) for the
argument 'add'; e.g.: add.params = list(color = "red").}

\item{...}{other parameters passing to \link[ggpubr:ggboxplot]{ggpubr::ggboxplot} or \link[ggpubr:ggviolin]{ggpubr::ggviolin}.}

\item{legend}{character specifying legend position. Allowed values are one of
c("top", "bottom", "left", "right", "none"). To remove the legend use
legend = "none". Legend position can be also specified using a numeric
vector c(x, y); see details section.}

\item{measure}{measure to estimate the exposure instability, can be one of 'RMSE', 'CV', 'MAE' and 'AbsDiff'.}
}
\value{
a \code{ggplot} object
}
\description{
See details for description.
}
\details{
Functions:
\itemize{
\item \link{show_sig_bootstrap_exposure} - this function plots exposures from bootstrap samples with both dotted boxplot.
The optimal exposure (the exposure from original input) is shown as triangle point. \strong{Only one sample can be plotted}.
\item \link{show_sig_bootstrap_error} - this function plots decomposition errors from bootstrap samples with both dotted boxplot.
The error from optimal solution (the decomposition error from original input) is shown as triangle point. \strong{Only one sample can be plotted}.
\item \link{show_sig_bootstrap_stability} - this function plots the signature exposure instability for specified signatures. Currently,
the instability measure supports 3 types:
\itemize{
\item 'RMSE' for Mean Root Squared Error (default) of bootstrap exposures and original exposures for each sample.
\item 'CV' for  Coefficient of Variation (CV) based on RMSE (i.e. \code{RMSE / btExposure_mean}).
\item 'MAE' for Mean Absolute Error of bootstrap exposures and original exposures for each sample.
\item 'AbsDiff' for Absolute Difference between mean bootstram exposure and original exposure.
}
}
}
\examples{
\donttest{
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
  laml <- read_maf(maf = laml.maf)
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    use_syn = TRUE
  )

  library(NMF)
  mt_sig <- sig_extract(mt_tally$nmf_matrix,
    n_sig = 3,
    nrun = 2,
    cores = 1
  )

  mat <- t(mt_tally$nmf_matrix)
  mat <- mat[, colSums(mat) > 0]
  bt_result <- sig_fit_bootstrap_batch(mat, sig = mt_sig, n = 10)
  ## Parallel computation
  ## bt_result = sig_fit_bootstrap_batch(mat, sig = mt_sig, n = 10, use_parallel = TRUE)

  ## At default, mean bootstrap exposure for each sample has been calculated
  p <- show_sig_bootstrap_exposure(bt_result, methods = c("QP"))
  ## Show bootstrap exposure (optimal exposure is shown as triangle)
  p1 <- show_sig_bootstrap_exposure(bt_result, methods = c("QP"), sample = "TCGA-AB-2802")
  p1
  p2 <- show_sig_bootstrap_exposure(bt_result,
    methods = c("QP"),
    sample = "TCGA-AB-3012",
    signatures = c("Sig1", "Sig2")
  )
  p2

  ## Show bootstrap error
  ## Similar to exposure above
  p <- show_sig_bootstrap_error(bt_result, methods = c("QP"))
  p
  p3 <- show_sig_bootstrap_error(bt_result, methods = c("QP"), sample = "TCGA-AB-2802")
  p3

  ## Show exposure (in)stability
  p4 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"))
  p4
  p5 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "MAE")
  p5
  p6 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "AbsDiff")
  p6
  p7 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "CV")
  p7
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}
}
}
\references{
Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604
}
\seealso{
\link{sig_fit_bootstrap_batch}, \link{sig_fit}, \link{sig_fit_bootstrap}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_copynumber_ascat.R
\name{read_copynumber_ascat}
\alias{read_copynumber_ascat}
\title{Read Copy Number Data from ASCAT Result Files}
\usage{
read_copynumber_ascat(x)
}
\arguments{
\item{x}{one or more \code{.rds} format files which contains \code{ASCAT} object from result of \code{ascat.runAscat()}
in \href{https://github.com/VanLoo-lab/ascat}{\strong{ASCAT}} package.}
}
\value{
a tidy \code{list}.
}
\description{
Note, the result is not a \code{CopyNumber} object, you need to generate it
by yourself.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{output_tally}
\alias{output_tally}
\title{Output Tally Result in Barplots}
\usage{
output_tally(x, result_dir, mut_type = "SBS")
}
\arguments{
\item{x}{a matrix with row representing components (motifs) and column
representing samples.}

\item{result_dir}{a result directory.}

\item{mut_type}{one of 'SBS', 'DBS', 'ID' or 'CN'.}
}
\value{
Nothing.
}
\description{
Output Tally Result in Barplots
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_profile_heatmap.R
\name{show_sig_profile_heatmap}
\alias{show_sig_profile_heatmap}
\title{Show Signature Profile with Heatmap}
\usage{
show_sig_profile_heatmap(
  Signature,
  mode = c("SBS", "DBS"),
  normalize = c("row", "column", "raw"),
  filters = NULL,
  x_lab = NULL,
  y_lab = NULL,
  legend_name = "auto",
  palette = "red",
  x_label_angle = 90,
  x_label_vjust = 1,
  x_label_hjust = 0.5,
  y_label_angle = 0,
  y_label_vjust = 0.5,
  y_label_hjust = 1,
  flip_xy = FALSE,
  sig_names = NULL,
  sig_orders = NULL,
  check_sig_names = TRUE
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix with row representing components (motifs) and column
representing signatures (column names must start with 'Sig').}

\item{mode}{one of "SBS" and "DBS".}

\item{normalize}{one of 'row', 'column', 'raw' and "feature", for row normalization (signature),
column normalization (component), raw data, row normalization by feature, respectively.
Of note, 'feature' only works when the mode is 'copynumber'.}

\item{filters}{a pattern used to select components to plot.}

\item{x_lab}{x label.}

\item{y_lab}{y label.}

\item{legend_name}{name of figure legend.}

\item{palette}{color for value.}

\item{x_label_angle}{angle for x axis text.}

\item{x_label_vjust}{vjust for x axis text.}

\item{x_label_hjust}{hjust for x axis text.}

\item{y_label_angle}{angle for y axis text.}

\item{y_label_vjust}{vjust for y axis text.}

\item{y_label_hjust}{hjust for y axis text.}

\item{flip_xy}{if \code{TRUE}, flip x axis and y axis.}

\item{sig_names}{subset signatures or set name of signatures, can be a character vector.
Default is \code{NULL}, prefix 'Sig' plus number is used.}

\item{sig_orders}{set order of signatures, can be a character vector.
Default is \code{NULL}, the signatures are ordered by alphabetical order.
If an integer vector set, only specified signatures are plotted.}

\item{check_sig_names}{if \code{TRUE}, check signature names when input is
a matrix, i.e., all signatures (colnames) must start with 'Sig'.}
}
\value{
a \code{ggplot} object.
}
\description{
This is a complementary function to \code{\link[=show_sig_profile]{show_sig_profile()}}, it is used for visualizing
some big signatures, i.e. SBS-1536, not all signatures are supported. See details for
current supported signatures.
}
\details{
Support:
\itemize{
\item SBS-24
\item SBS-96
\item SBS-384
\item SBS-1536
\item SBS-6144
\item DBS-78
\item DBS-186
}
}
\examples{
# Load SBS signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p1 <- show_sig_profile_heatmap(sig2, mode = "SBS")
p1
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_consensusmap.R
\name{show_sig_consensusmap}
\alias{show_sig_consensusmap}
\title{Show Signature Consensus Map}
\usage{
show_sig_consensusmap(
  sig,
  main = "Consensus matrix",
  tracks = c("consensus:", "silhouette:"),
  lab_row = NA,
  lab_col = NA,
  ...
)
}
\arguments{
\item{sig}{a \code{Signature} object obtained from \link{sig_extract}.}

\item{main}{Main title as a character string or a grob.}

\item{tracks}{Special additional annotation tracks to
  highlight associations between basis components and
  sample clusters: \describe{ \item{basis}{matches each row
  (resp. column) to the most contributing basis component
  in \code{basismap} (resp. \code{coefmap}). In
  \code{basismap} (resp. \code{coefmap}), adding a track
  \code{':basis'} to \code{annCol} (resp. \code{annRow})
  makes the column (resp. row) corresponding to the
  component being also highlited using the mathcing
  colours.} }}

\item{lab_row}{labels for the rows.}

\item{lab_col}{labels for the columns.}

\item{...}{other parameters passing to \code{NMF::consensusmap()}.}
}
\value{
nothing
}
\description{
This function is a wrapper of \code{NMF::consensusmap()}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_catalogue.R
\name{show_catalogue}
\alias{show_catalogue}
\title{Show Alteration Catalogue Profile}
\usage{
show_catalogue(
  catalogue,
  mode = c("SBS", "copynumber", "DBS", "ID", "RS"),
  method = "Wang",
  normalize = c("raw", "row", "feature"),
  style = c("default", "cosmic"),
  samples = NULL,
  samples_name = NULL,
  x_lab = "Components",
  y_lab = "Counts",
  ...
)
}
\arguments{
\item{catalogue}{result from \link{sig_tally} or a
matrix with row representing components (motifs) and
column representing samples}

\item{mode}{signature type for plotting, now supports 'copynumber', 'SBS',
'DBS', 'ID' and 'RS' (genome rearrangement signature).}

\item{method}{method for copy number feature classification in \link{sig_tally},
can be one of "Wang" ("W"), "S".}

\item{normalize}{normalize method.}

\item{style}{plot style, one of 'default' and 'cosmic'.}

\item{samples}{default is \code{NULL}, show sum of all samples in one row.
If not \code{NULL}, show specified samples.}

\item{samples_name}{set the sample names shown in plot.}

\item{x_lab}{x axis lab.}

\item{y_lab}{y axis lab.}

\item{...}{other arguments passing to \link{show_sig_profile}.}
}
\value{
a \code{ggplot} object
}
\description{
Show Alteration Catalogue Profile
}
\examples{
data("simulated_catalogs")
p <- show_catalogue(simulated_catalogs$set1, style = "cosmic")
p
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{output_fit}
\alias{output_fit}
\title{Output Signature Fitting Results}
\usage{
output_fit(x, result_dir, mut_type = "SBS", sig_db = mut_type)
}
\arguments{
\item{x}{result from \link{sig_fit}.}

\item{result_dir}{a result directory.}

\item{mut_type}{one of 'SBS', 'DBS', 'ID' or 'CN'.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}
}
\value{
Nothing.
}
\description{
Output Signature Fitting Results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{output_sig}
\alias{output_sig}
\title{Output Signature Results}
\usage{
output_sig(sig, result_dir, mut_type = "SBS", sig_db = mut_type)
}
\arguments{
\item{sig}{a \code{Signature} object.}

\item{result_dir}{a result directory.}

\item{mut_type}{one of 'SBS', 'DBS', 'ID' or 'CN'.}

\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}
}
\value{
Nothing.
}
\description{
Output Signature Results
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sv.R
\name{read_sv_as_rs}
\alias{read_sv_as_rs}
\title{Read Structural Variation Data as RS object}
\usage{
read_sv_as_rs(input)
}
\arguments{
\item{input}{a \code{data.frame} or a file with the following columns:
"sample", "chr1", "start1", "end1", "chr2", "start2", "end2", "strand1", "strand2", "svclass".
NOTE: If column "svclass" already exists in input, "strand1" and "strand2" are optional.
If "svclass" is not provided, \code{read_sv_as_rs()} will compute it by
"strand1","strand2"(strand1/strand2),"chr1" and "chr2":
\itemize{
\item translocation, if mates are on different chromosomes.
\item inversion (+/-) and (-/+), if mates on the same chromosome.
\item deletion (+/+), if mates on the same chromosome.
\item tandem-duplication (-/-), if mates on the same chromosome.
}}
}
\value{
a \code{list}
}
\description{
Read Structural Variation Data as RS object
}
\examples{
sv <- readRDS(system.file("extdata", "toy_sv.rds", package = "sigminer", mustWork = TRUE))
rs <- read_sv_as_rs(sv)
# svclass is optional
rs2 <- read_sv_as_rs(sv[, setdiff(colnames(sv), "svclass")])
identical(rs, rs2)
\donttest{
tally_rs <- sig_tally(rs)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_profile.R
\name{show_sig_profile}
\alias{show_sig_profile}
\title{Show Signature Profile}
\usage{
show_sig_profile(
  Signature,
  mode = c("SBS", "copynumber", "DBS", "ID", "RS"),
  method = "Wang",
  by_context = FALSE,
  normalize = c("row", "column", "raw", "feature"),
  y_tr = NULL,
  filters = NULL,
  feature_setting = sigminer::CN.features,
  style = c("default", "cosmic"),
  palette = use_color_style(style, ifelse(by_context, "SBS", mode), method),
  set_gradient_color = FALSE,
  free_space = "free_x",
  rm_panel_border = style == "cosmic",
  rm_grid_line = style == "cosmic",
  rm_axis_text = FALSE,
  bar_border_color = ifelse(style == "default", "grey50", "white"),
  bar_width = 0.7,
  paint_axis_text = TRUE,
  x_label_angle = ifelse(mode == "copynumber" & !(startsWith(method, "T") | method ==
    "X"), 60, 90),
  x_label_vjust = ifelse(mode == "copynumber" & !(startsWith(method, "T") | method ==
    "X"), 1, 0.5),
  x_label_hjust = 1,
  x_lab = "Components",
  y_lab = "auto",
  y_limits = NULL,
  params = NULL,
  show_cv = FALSE,
  params_label_size = 3,
  params_label_angle = 60,
  y_expand = 1,
  digits = 2,
  base_size = 12,
  font_scale = 1,
  sig_names = NULL,
  sig_orders = NULL,
  check_sig_names = TRUE
)
}
\arguments{
\item{Signature}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract},
or just a raw signature matrix with row representing components (motifs) and column
representing signatures (column names must start with 'Sig').}

\item{mode}{signature type for plotting, now supports 'copynumber', 'SBS',
'DBS', 'ID' and 'RS' (genome rearrangement signature).}

\item{method}{method for copy number feature classification in \link{sig_tally},
can be one of "Wang" ("W"), "S".}

\item{by_context}{for specific use.}

\item{normalize}{one of 'row', 'column', 'raw' and "feature", for row normalization (signature),
column normalization (component), raw data, row normalization by feature, respectively.
Of note, 'feature' only works when the mode is 'copynumber'.}

\item{y_tr}{a function (e.g. \code{log10}) to transform y axis before plotting.}

\item{filters}{a pattern used to select components to plot.}

\item{feature_setting}{a \code{data.frame} used for classification.
\strong{Only used when method is "Wang" ("W")}.
Default is \link{CN.features}. Users can also set custom input with "feature",
"min" and "max" columns available. Valid features can be printed by
\code{unique(CN.features$feature)}.}

\item{style}{plot style, one of 'default' and 'cosmic', works when
parameter \code{set_gradient_color} is \code{FALSE}.}

\item{palette}{palette used to plot when \code{set_gradient_color} is \code{FALSE},
default use a built-in palette according to parameter \code{style}.}

\item{set_gradient_color}{default is \code{FALSE}, if \code{TRUE}, use gradient colors
to fill bars.}

\item{free_space}{default is 'free_x'. If "fixed", all panels have the same size.
If "free_y" their height will be proportional to the length of the y scale;
if "free_x" their width will be proportional to the length of the x scale;
or if "free" both height and width will vary.
This setting has no effect unless the appropriate scales also vary.}

\item{rm_panel_border}{default is \code{TRUE} for style 'cosmic',
remove panel border to keep plot tight.}

\item{rm_grid_line}{default is \code{FALSE}, if \code{TRUE}, remove grid lines of plot.}

\item{rm_axis_text}{default is \code{FALSE}, if \code{TRUE}, remove component texts.
This is useful when multiple signature profiles are plotted together.}

\item{bar_border_color}{the color of bar border.}

\item{bar_width}{bar width. By default, set to 70\% of the resolution of the
data.}

\item{paint_axis_text}{if \code{TRUE}, color on text of x axis.}

\item{x_label_angle}{font angle for x label.}

\item{x_label_vjust}{font vjust for x label.}

\item{x_label_hjust}{font hjust for x label.}

\item{x_lab}{x axis lab.}

\item{y_lab}{y axis lab.}

\item{y_limits}{limits to expand in y axis. e.g., \code{0.2}, \code{c(0, 0.3)}.}

\item{params}{params \code{data.frame} of components, obtained from \link{sig_tally}.}

\item{show_cv}{default is \code{FALSE}, if \code{TRUE}, show coefficient of variation when
\code{params} is not \code{NULL}.}

\item{params_label_size}{font size for params label.}

\item{params_label_angle}{font angle for params label.}

\item{y_expand}{y expand height for plotting params of copy number signatures.}

\item{digits}{digits for plotting params of copy number signatures.}

\item{base_size}{overall font size.}

\item{font_scale}{a number used to set font scale.}

\item{sig_names}{subset signatures or set name of signatures, can be a character vector.
Default is \code{NULL}, prefix 'Sig' plus number is used.}

\item{sig_orders}{set order of signatures, can be a character vector.
Default is \code{NULL}, the signatures are ordered by alphabetical order.
If an integer vector set, only specified signatures are plotted.}

\item{check_sig_names}{if \code{TRUE}, check signature names when input is
a matrix, i.e., all signatures (colnames) must start with 'Sig'.}
}
\value{
a \code{ggplot} object
}
\description{
Who don't like to show a barplot for signature profile? This is for it.
}
\examples{
# Load SBS signature
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p1 <- show_sig_profile(sig2, mode = "SBS")
p1

# Use 'y_tr' option to transform values in y axis
p11 <- show_sig_profile(sig2, mode = "SBS", y_tr = function(x) x * 100)
p11

# Load copy number signature from method "W"
load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
  package = "sigminer", mustWork = TRUE
))
# Show signature profile
p2 <- show_sig_profile(sig,
  style = "cosmic",
  mode = "copynumber",
  method = "W",
  normalize = "feature"
)
p2

# Visualize rearrangement signatures
s <- get_sig_db("RS_Nik_lab")
ss <- s$db[, 1:3]
colnames(ss) <- c("Sig1", "Sig2", "Sig3")
p3 <- show_sig_profile(ss, mode = "RS", style = "cosmic")
p3
}
\seealso{
\link{show_sig_profile_loop}, \link{show_sig_profile_heatmap}
}
\author{
Shixiang Wang
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_feature_association.R
\name{get_sig_feature_association}
\alias{get_sig_feature_association}
\title{Calculate Association between Signature Exposures and Other Features}
\usage{
get_sig_feature_association(
  data,
  cols_to_sigs,
  cols_to_features,
  type = "ca",
  method_co = c("spearman", "pearson", "kendall"),
  method_ca = stats::wilcox.test,
  min_n = 0.01,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{data}{a \code{data.frame} contains signature exposures and other features}

\item{cols_to_sigs}{colnames for signature exposure}

\item{cols_to_features}{colnames for other features}

\item{type}{a character vector containing 'ca' for categorical variable and 'co' for continuous variable,
it must have the same length as \code{cols_to_features}.}

\item{method_co}{method for continuous variable, default is "spearman", could also be "pearson" and "kendall".}

\item{method_ca}{method for categorical variable, default is "wilcox.test"}

\item{min_n}{a minimal fraction (e.g. 0.01) or a integer number (e.g. 10) for filtering some variables with few positive events.
Default is 0.01.}

\item{verbose}{if \code{TRUE}, print extra message.}

\item{...}{other arguments passing to test functions, like \code{cor.test}.}
}
\value{
a \code{list}. For 'co' features, 'measure' means correlation coefficient.
For 'ca' features, 'measure' means difference in means of signature exposure.
}
\description{
Association of signature exposures with other features will be performed using one of two procedures:
for a continuous association variable (including ordinal variable), correaltion is performed;
for a binary association variable, samples will be divided into two groups and Mann-Whitney U-test
is performed to test for differences in signature exposure medians between the two groups.
See \link{get_tidy_association} for cleaning association result.
}
\seealso{
\link{get_tidy_association}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/same_size_clustering.R
\name{same_size_clustering}
\alias{same_size_clustering}
\title{Same Size Clustering}
\usage{
same_size_clustering(
  mat,
  diss = FALSE,
  clsize = NULL,
  algo = c("nnit", "hcbottom", "kmvar"),
  method = c("maxd", "random", "mind", "elki", "ward.D", "average", "complete",
    "single")
)
}
\arguments{
\item{mat}{a data/distance matrix.}

\item{diss}{if \code{TRUE}, treat \code{mat} as a distance matrix.}

\item{clsize}{integer, number of sample within a cluster.}

\item{algo}{algorithm.}

\item{method}{method.}
}
\value{
a vector.
}
\description{
This is a wrapper for several implementation that classify samples into
same size clusters, the details please see \href{http://jmonlong.github.io/Hippocamplus/2018/06/09/cluster-same-size/}{this blog}.
The source code is modified based on code from the blog.
}
\examples{
set.seed(1234L)
x <- rbind(
  matrix(rnorm(100, sd = 0.3), ncol = 2),
  matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2)
)
colnames(x) <- c("x", "y")

y1 <- same_size_clustering(x, clsize = 10)
y11 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, diss = TRUE)

y2 <- same_size_clustering(x, clsize = 10, algo = "hcbottom", method = "ward.D")

y3 <- same_size_clustering(x, clsize = 10, algo = "kmvar")
y33 <- same_size_clustering(as.matrix(dist(x)), clsize = 10, algo = "kmvar", diss = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_circos.R
\name{show_cn_circos}
\alias{show_cn_circos}
\title{Show Copy Number Profile in Circos}
\usage{
show_cn_circos(
  data,
  samples = NULL,
  show_title = TRUE,
  chrs = paste0("chr", 1:22),
  genome_build = c("hg19", "hg38", "mm10", "mm9"),
  col = NULL,
  side = "inside",
  ...
)
}
\arguments{
\item{data}{a \link{CopyNumber} object or a \code{data.frame} containing at least 'chromosome', 'start',
'end', 'segVal' these columns.}

\item{samples}{default is \code{NULL}, can be a chracter vector representing multiple samples or
number of samples to show.
If data argument is a \code{data.frame}, a column called sample must exist.}

\item{show_title}{if \code{TRUE} (default), show title with sample ID.}

\item{chrs}{chromosomes start with 'chr'.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{col}{colors for the heatmaps. If it is \code{NULL}, set to
\code{circlize::colorRamp2(c(1, 2, 4), c("blue", "black", "red"))}.}

\item{side}{side of the heatmaps.}

\item{...}{other parameters passing to \link[circlize:circos.genomicHeatmap]{circlize::circos.genomicHeatmap}.}
}
\value{
a circos plot
}
\description{
Another visualization method for copy number profile like \link{show_cn_profile}.
}
\examples{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))
\donttest{
show_cn_circos(cn, samples = 1)
show_cn_circos(cn, samples = "TCGA-99-7458-01A-11D-2035-01")

## Remove title
show_cn_circos(cn, samples = 1, show_title = FALSE)

## Subset chromosomes
show_cn_circos(cn, samples = 1, chrs = c("chr1", "chr2", "chr3"))

## Arrange plots
layout(matrix(1:4, 2, 2))
show_cn_circos(cn, samples = 4)

layout(1) # reset layout
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cytobands.hg38}
\alias{cytobands.hg38}
\title{Location of Chromosome Cytobands at Genome Build hg38}
\format{
A data.frame
}
\source{
from UCSC
}
\description{
Location of Chromosome Cytobands at Genome Build hg38
}
\examples{
data(cytobands.hg38)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_genome_annotation.R
\name{get_genome_annotation}
\alias{get_genome_annotation}
\title{Get Genome Annotation}
\usage{
get_genome_annotation(
  data_type = c("chr_size", "centro_loc", "cytobands", "transcript"),
  chrs = paste0("chr", c(1:22, "X", "Y")),
  genome_build = c("hg19", "hg38", "mm10", "mm9")
)
}
\arguments{
\item{data_type}{'chr_size' for chromosome size,
'centro_loc' for location of centromeres,
'cytobands' for location of chromosome cytobands
and 'transcript' for location of transcripts.}

\item{chrs}{chromosomes start with 'chr'}

\item{genome_build}{one of 'hg19', 'hg38'}
}
\value{
a \code{data.frame} containing annotation data
}
\description{
Get Genome Annotation
}
\examples{
df1 <- get_genome_annotation()
df1

df2 <- get_genome_annotation(genome_build = "hg38")
df2

df3 <- get_genome_annotation(data_type = "centro_loc")
df3

df4 <- get_genome_annotation(data_type = "centro_loc", genome_build = "hg38")
df4

df5 <- get_genome_annotation(data_type = "cytobands")
df5

df6 <- get_genome_annotation(data_type = "cytobands", genome_build = "hg38")
df6
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_enrichment.R
\name{group_enrichment}
\alias{group_enrichment}
\title{General Group Enrichment Analysis}
\usage{
group_enrichment(
  df,
  grp_vars = NULL,
  enrich_vars = NULL,
  cross = TRUE,
  co_method = c("t.test", "wilcox.test")
)
}
\arguments{
\item{df}{a \code{data.frame}.}

\item{grp_vars}{character vector specifying group variables to split samples
into subgroups (at least 2 subgroups, otherwise this variable will be skipped).}

\item{enrich_vars}{character vector specifying measure variables to be compared.
If variable is not numeric, only binary cases are accepted in the form of
\code{TRUE/FALSE} or \code{P/N} (P for positive cases and N for negative cases).
Of note, \code{NA} values set to negative cases.}

\item{cross}{logical, default is \code{TRUE}, combine all situations provided by
\code{grp_vars} and \code{enrich_vars}. For examples, \code{c('A', 'B')} and \code{c('C', 'D')}
will construct 4 combinations(i.e. "AC", "AD", "BC" and "BD"). A variable can
not be in both \code{grp_vars} and \code{enrich_vars}, such cases will be automatically
drop. If \code{FALSE}, use pairwise combinations, see section "examples" for use cases.}

\item{co_method}{test method for continuous variable, default is 't.test'.}
}
\value{
a \code{data.table} with following columns:
\itemize{
\item \code{grp_var}: group variable name.
\item \code{enrich_var}: enrich variable (variable to be compared) name.
\item \code{grp1}: the first group name, should be a member in \code{grp_var} column.
\item \code{grp2}: the remaining samples, marked as 'Rest'.
\item \code{grp1_size}: sample size for \code{grp1}.
\item \code{grp1_pos_measure}: for binary variable, it stores the proportion of
positive cases in \code{grp1}; for continuous variable, it stores mean value.
\item \code{grp2_size}: sample size for \code{grp2}.
\item \code{grp2_pos_measure}: same as \code{grp1_pos_measure} but for \code{grp2}.
\item \code{measure_observed}: for binary variable, it stores odds ratio;
for continuous variable, it stores scaled mean ratio.
\item \code{measure_tested}: only for binary variable, it stores
estimated odds ratio and its 95\% CI from \code{fisher.test()}.
\item \code{p_value}: for binary variable, it stores p value from \code{fisher.test()};
for continuous variable, it stores value from \code{wilcox.test()} or \code{t.test()}.
\item \code{type}: one of "binary" and "continuous".
\item \code{method}: one of "fish.test", "wilcox.test" and "t.test".
}
}
\description{
This function takes a \code{data.frame} as input, compares proportion of positive
cases or mean measure in one subgroup and the remaining samples.
}
\examples{
set.seed(1234)
df <- dplyr::tibble(
  g1 = factor(abs(round(rnorm(99, 0, 1)))),
  g2 = rep(LETTERS[1:4], c(50, 40, 8, 1)),
  e1 = sample(c("P", "N"), 99, replace = TRUE),
  e2 = rnorm(99)
)

print(str(df))
print(head(df))

# Compare g1:e1, g1:e2, g2:e1 and g2:e2
x1 <- group_enrichment(df, grp_vars = c("g1", "g2"), enrich_vars = c("e1", "e2"))
x1

# Only compare g1:e1, g2:e2
x2 <- group_enrichment(df,
  grp_vars = c("g1", "g2"),
  enrich_vars = c("e1", "e2"),
  co_method = "wilcox.test",
  cross = FALSE
)
x2

# Visualization
p1 <- show_group_enrichment(x1, fill_by_p_value = TRUE)
p1
p2 <- show_group_enrichment(x1, fill_by_p_value = FALSE)
p2
p3 <- show_group_enrichment(x1, return_list = TRUE)
p3
}
\seealso{
\link{show_group_enrichment}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cosine.R
\name{cosine}
\alias{cosine}
\title{Calculate Cosine Measures}
\usage{
cosine(x, y)
}
\arguments{
\item{x}{a numeric vector or matrix with column representing vector to calculate similarity.}

\item{y}{must be same format as \code{x}.}
}
\value{
a numeric value or \code{matrix}.
}
\description{
Calculate Cosine Measures
}
\examples{
x <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
y <- c(0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0)
z1 <- cosine(x, y)
z1
z2 <- cosine(matrix(x), matrix(y))
z2
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_vcf.R
\name{read_xena_variants}
\alias{read_xena_variants}
\title{Read UCSC Xena Variant Format Data as MAF Object}
\usage{
read_xena_variants(path)
}
\arguments{
\item{path}{a path to variant file.}
}
\value{
a \code{MAF} object.
}
\description{
Read UCSC Xena Variant Format Data as MAF Object
}
\examples{
\donttest{
if (requireNamespace("UCSCXenaTools")) {
  library(UCSCXenaTools)
  options(use_hiplot = TRUE)
  example_file <- XenaGenerate(subset = XenaDatasets == "mc3/ACC_mc3.txt") \%>\%
    XenaQuery() \%>\%
    XenaDownload()
  x <- read_xena_variants(example_file$destfiles)
  x@data
  y <- sig_tally(x)
  y
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature_obj_operation.R
\name{sig_operation}
\alias{sig_operation}
\alias{sig_names}
\alias{sig_modify_names}
\alias{sig_number}
\alias{sig_attrs}
\alias{sig_signature}
\alias{sig_exposure}
\title{Obtain or Modify Signature Information}
\usage{
sig_names(sig)

sig_modify_names(sig, new_names)

sig_number(sig)

sig_attrs(sig)

sig_signature(sig, normalize = c("row", "column", "raw", "feature"))

sig_exposure(sig, type = c("absolute", "relative"))
}
\arguments{
\item{sig}{a \code{Signature} object obtained either from \link{sig_extract} or \link{sig_auto_extract}.}

\item{new_names}{new signature names.}

\item{normalize}{one of 'row', 'column', 'raw' and "feature", for row normalization (signature),
column normalization (component), raw data, row normalization by feature, respectively.}

\item{type}{one of 'absolute' and 'relative'.}
}
\value{
a \code{Signature} object or data.
}
\description{
Obtain or Modify Signature Information
}
\examples{
## Operate signature names
load(system.file("extdata", "toy_mutational_signature.RData",
  package = "sigminer", mustWork = TRUE
))
sig_names(sig2)
cc <- sig_modify_names(sig2, new_names = c("Sig2", "Sig1", "Sig3"))
sig_names(cc)

# The older names are stored in tags.
print(attr(cc, "tag"))
## Get signature number
sig_number(sig2)
## Get signature attributes
sig_number(sig2)
## Get signature matrix
z <- sig_signature(sig2)
z <- sig_signature(sig2, normalize = "raw")
## Get exposure matrix
## Of note, this is different from get_sig_exposure()
## it returns a matrix instead of data table.
z <- sig_exposure(sig2) # it is same as sig$Exposure
z <- sig_exposure(sig2, type = "relative") # it is same as sig2$Exposure.norm
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{centromeres.mm9}
\alias{centromeres.mm9}
\title{Location of Centromeres at Genome Build mm9}
\format{
A data.frame
}
\source{
Generate from \url{https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/}
with code:\if{html}{\out{<div class="sh">}}\preformatted{for i in $(seq 1 19) X Y;
do
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/chr$\{i\}_gap.txt.gz
done
}\if{html}{\out{</div>}}
}
\description{
Location of Centromeres at Genome Build mm9
}
\examples{
data(centromeres.mm9)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tidy_association.R
\name{get_tidy_association}
\alias{get_tidy_association}
\title{Get Tidy Signature Association Results}
\usage{
get_tidy_association(cor_res, p_adjust = FALSE, method = "fdr")
}
\arguments{
\item{cor_res}{data returned by \code{\link[=get_sig_feature_association]{get_sig_feature_association()}}}

\item{p_adjust}{logical, if \code{TRUE}, adjust p values by data type.}

\item{method}{p value correction method, see \link[stats:p.adjust]{stats::p.adjust} for
more detail.}
}
\value{
a \code{data.frame}
}
\description{
Get Tidy Signature Association Results
}
\seealso{
\link{get_sig_feature_association}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cytobands.hg19}
\alias{cytobands.hg19}
\title{Location of Chromosome Cytobands at Genome Build hg19}
\format{
A data.frame
}
\source{
from UCSC
}
\description{
Location of Chromosome Cytobands at Genome Build hg19
}
\examples{
data(cytobands.hg19)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class.R
\docType{class}
\name{CopyNumber-class}
\alias{CopyNumber-class}
\alias{CopyNumber}
\title{Class CopyNumber}
\description{
S4 class for storing summarized absolute copy number profile.
}
\section{Slots}{

\describe{
\item{\code{data}}{data.table of absolute copy number calling.}

\item{\code{summary.per.sample}}{data.table of copy number variation summary per sample.}

\item{\code{genome_build}}{genome build version, should be one of 'hg19' or 'hg38'.}

\item{\code{genome_measure}}{Set 'called' will use autosomo called segments size to compute total size
for CNA burden calculation, this option is useful for WES and target sequencing.
Set 'wg' will autosome size from genome build, this option is useful for WGS, SNP etc..}

\item{\code{annotation}}{data.table of annotation for copy number segments.}

\item{\code{dropoff.segs}}{data.table of copy number segments dropped from raw input.}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_features.R
\name{show_cn_features}
\alias{show_cn_features}
\title{Show Copy Number Feature Distributions}
\usage{
show_cn_features(
  features,
  method = "Wang",
  rm_outlier = FALSE,
  ylab = NULL,
  log_y = FALSE,
  return_plotlist = FALSE,
  base_size = 12,
  nrow = 2,
  align = "hv",
  ...
)
}
\arguments{
\item{features}{a feature \code{list} generate from \link{sig_tally} function.}

\item{method}{method for feature classification, can be one of
"Wang" ("W"), "S" (for method described in Steele et al. 2019).}

\item{rm_outlier}{default is \code{FALSE}, if \code{TRUE}, remove outliers. Only
works when method is "Wang" ("W").}

\item{ylab}{lab of y axis.}

\item{log_y}{logical, if \code{TRUE}, show \code{log10} based y axis, only
works for input from "Wang" ("W") method.}

\item{return_plotlist}{if \code{TRUE}, return a list of ggplot objects but a combined plot.}

\item{base_size}{overall font size.}

\item{nrow}{(optional) Number of rows in the plot grid.}

\item{align}{(optional) Specifies whether graphs in the grid should be horizontally ("h") or
vertically ("v") aligned. Options are "none" (default), "hv" (align in both directions), "h", and "v".}

\item{...}{other options pass to \code{\link[cowplot]{plot_grid}} function of \code{cowplot} package.}
}
\value{
a \code{ggplot} object
}
\description{
Show Copy Number Feature Distributions
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{centromeres.hg19}
\alias{centromeres.hg19}
\title{Location of Centromeres at Genome Build hg19}
\format{
A data.frame
}
\source{
Generate from UCSC gold path
}
\description{
Location of Centromeres at Genome Build hg19
}
\examples{
data(centromeres.hg19)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_reconstructed_similarity.R
\name{get_sig_rec_similarity}
\alias{get_sig_rec_similarity}
\title{Get Reconstructed Profile Cosine Similarity, RSS, etc.}
\usage{
get_sig_rec_similarity(Signature, nmf_matrix)
}
\arguments{
\item{Signature}{a \code{Signature} object.}

\item{nmf_matrix}{a \code{matrix} used for NMF decomposition with rows indicate samples and columns indicate components.}
}
\value{
a \code{data.table}.
}
\description{
See \link{bp_extract_signatures} for examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoring.R
\name{scoring}
\alias{scoring}
\title{Score Copy Number Profile}
\usage{
scoring(object, TD_size_cutoff = c(1000, 1e+05, 2e+06), TD_cn_cutoff = Inf)
}
\arguments{
\item{object}{a object of \link{CopyNumber}.}

\item{TD_size_cutoff}{a length-3 numeric vector used to specify the start, midpoint, end
segment size for determining tandem duplication size range, midpoint is used to split
TD into short TD and long TD. Default is 1Kb to 100Kb for short TD, 100Kb to 2Mb for long
TD.}

\item{TD_cn_cutoff}{a number defining the maximum copy number of TD,
default is \code{Inf}, i.e. no cutoff.}
}
\value{
a \code{data.table} with following scores:
\itemize{
\item cnaBurden: CNA burden representing the altered genomic fraction as previously reported.
\item cnaLoad: CNA load representing the quantity of copy number alteration.
\item MACN: mean altered copy number (MACN) reflecting the property of altered copy number segments,
calculated as
\deqn{MACN = \frac{\sum_{i} CN_i}{N_{cnv}}}
where \eqn{CN_i} is the copy number of altered segment \eqn{i}, \eqn{N_{cnv}} is
the number of CNV.
\item weightedMACN: same as MACN but weighted with segment length.
\deqn{MACN_{weighted} = \frac{\sum_{i} (CN_i \times L_{i})}{ \sum_{i} L_{i} }}
where \eqn{L_{i}} is the length of altered copy number segment \eqn{i}.
\item Ploidy: ploidy, the formula is same as \code{weightedMACN} but using all copy number segments instead of
altered copy number segments.
\item TDP_pnas: tandem duplication phenotype score from \url{https://www.pnas.org/content/113/17/E2373},
the threshold \code{k} in reference is omitted.
\deqn{TDP = - \frac{\sum_{chr} |TD_{obs}-TD_{exp}|}{TD_{total}}}
where \eqn{TD_{total}} is the number of TD, \eqn{TD_{obs}} and
\eqn{TD_exp} are observed number of TD and expected number of TD for each chromosome.
\item TDP: tandem duplication score used defined by our group work,
TD represents segment with copy number greater than 2.
\deqn{TD = \frac{TD_{total}}{\sum_{chr} |TD_{obs}-TD_{exp}|+1}}
\item sTDP: TDP score for short TD.
\item lTDP: TDP score for long TD.
\item TDP_size : TDP region size (Mb).
\item sTDP_size: sTDP region size (Mb).
\item lTDP_size: lTDP region size(Mb).
\item Chromoth_state: chromothripsis state score,
according to reference \doi{10.1016/j.cell.2013.02.023},
chromothripsis frequently leads to massive loss of segments on
the affected chromosome with segmental losses being interspersed with regions displaying
normal (disomic) copy-number (e.g., copy-number states oscillating between
copy-number = 1 and copy-number = 2), form tens to hundreds of locally clustered DNA rearrangements.
Most of methods use both SV and CNV to infer chromothripsis, here we roughly quantify it with
\deqn{\sum_{chr}{N_{OsCN}^2}}
where \eqn{N_{OsCN}} is the number of oscillating copy number pattern "2-1-2" for each chromosome.
}
}
\description{
Returns quantification of copy number profile and events including
tandem duplication and Chromothripisis etc.
Only copy number data from autosome is used here.
\strong{Some of the quantification methods are rough,
you use at your risk}. You should do some extra work to check the
result scores.
}
\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

d <- scoring(cn)
d

d2 <- scoring(cn, TD_cn_cutoff = 4L)
d2
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{hello}
\alias{hello}
\title{Say Hello to Users}
\usage{
hello()
}
\description{
Say Hello to Users
}
\examples{
hello()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sig_fit_bootstrap.R
\name{report_bootstrap_p_value}
\alias{report_bootstrap_p_value}
\title{Report P Values from bootstrap Results}
\usage{
report_bootstrap_p_value(x, thresholds = c(0.01, 0.05, 0.1))
}
\arguments{
\item{x}{a (list of) result from \link{sig_fit_bootstrap}.}

\item{thresholds}{a vector of relative exposure threshold for calculating p values.}
}
\value{
a (list of) \code{matrix}
}
\description{
See examples in \link{sig_fit_bootstrap}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_Aneuploidy_score.R
\name{get_Aneuploidy_score}
\alias{get_Aneuploidy_score}
\title{Get Aneuploidy Score from Copy Number Profile}
\usage{
get_Aneuploidy_score(
  data,
  ploidy_df = NULL,
  genome_build = "hg19",
  rm_black_arms = FALSE
)
}
\arguments{
\item{data}{a CopyNumber object or a \code{data.frame} containing at least
'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{ploidy_df}{default is \code{NULL}, compute ploidy by segment-size weighted copy number
aross autosome, see \link{get_cn_ploidy}. You can also provide a \code{data.frame} with 'sample'
and 'ploidy' columns.}

\item{genome_build}{genome build version, should be 'hg19', 'hg38', 'mm9' or 'mm10'.}

\item{rm_black_arms}{if \code{TRUE}, remove short arms of chr13/14/15/21/22 from calculation
as documented in reference #3.}
}
\value{
A \code{data.frame}
}
\description{
This implements a Cohen-Sharir method (see reference) like "Aneuploidy Score" computation.
You can read the source code to see how it works. Basically, it follows
the logic of Cohen-Sharir method but with some difference in detail implementation.
Their results should be counterpart, but with no data validation for now.
\strong{Please raise an issue if you find problem/bugs in this function}.
}
\examples{
# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

df <- get_Aneuploidy_score(cn)
df

df2 <- get_Aneuploidy_score(cn@data)
df2

df3 <- get_Aneuploidy_score(cn@data,
  ploidy_df = get_cn_ploidy(cn@data)
)
df3
}
\references{
\itemize{
\item Cohen-Sharir, Y., McFarland, J. M., Abdusamad, M., Marquis, C., Bernhard, S. V., Kazachkova, M., ... & Ben-David, U. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature, 1-6.
\item Logic reference: \url{https://github.com/quevedor2/aneuploidy_score/}.
\item Taylor, Alison M., et al. "Genomic and functional approaches to understanding cancer aneuploidy." Cancer cell 33.4 (2018): 676-689.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrich_component_strand_bias.R
\name{enrich_component_strand_bias}
\alias{enrich_component_strand_bias}
\title{Performs Strand Bias Enrichment Analysis for a Given Sample-by-Component Matrix}
\usage{
enrich_component_strand_bias(mat)
}
\arguments{
\item{mat}{a sample-by-component matrix from \link{sig_tally} with strand bias labels "T:" and "B:".}
}
\value{
a \code{data.table} sorted by \code{p_value}.
}
\description{
See \link{sig_tally} for examples.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_group_mapping.R
\name{show_group_mapping}
\alias{show_group_mapping}
\title{Map Groups using Sankey}
\usage{
show_group_mapping(
  data,
  col_to_flow,
  cols_to_map,
  include_sig = FALSE,
  fill_na = FALSE,
  title = NULL,
  xlab = NULL,
  ylab = NULL,
  custom_theme = cowplot::theme_minimal_hgrid()
)
}
\arguments{
\item{data}{a \code{data.frame} containing signature group and other categorical groups.}

\item{col_to_flow}{length-1 character showing the column to flow, typically a signature group.}

\item{cols_to_map}{character vector showing colnames of other groups.}

\item{include_sig}{default if \code{FALSE}, if \code{TRUE}, showing signature group.}

\item{fill_na}{length-1 string to fill NA, default is \code{FALSE}.}

\item{title}{the title.}

\item{xlab}{label for x axis.}

\item{ylab}{label for y axis.}

\item{custom_theme}{theme for plotting, default is \code{cowplot::theme_minimal_hgrid()}.}
}
\value{
a \code{ggplot} object
}
\description{
This feature is designed for signature analysis. However, users can also use
it in other similar situations.
}
\examples{
data <- dplyr::tibble(
  Group1 = rep(LETTERS[1:5], each = 10),
  Group2 = rep(LETTERS[6:15], each = 5),
  zzzz = c(rep("xx", 20), rep("yy", 20), rep(NA, 10))
)
p1 <- show_group_mapping(data, col_to_flow = "Group1", cols_to_map = colnames(data)[-1])
p1

p2 <- show_group_mapping(data,
  col_to_flow = "Group1", cols_to_map = colnames(data)[-1],
  include_sig = TRUE
)
p2
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cn_freq_table.R
\name{get_cn_freq_table}
\alias{get_cn_freq_table}
\title{Get CNV Frequency Table}
\usage{
get_cn_freq_table(
  data,
  genome_build = "hg19",
  cutoff = 2L,
  resolution_factor = 1L
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{cutoff}{copy number value cutoff for splitting data into AMP and DEL.
The values equal to cutoff are discarded. Default is \code{2}, you can also set
a length-2 vector, e.g. \code{c(2, 2)}.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}
}
\value{
a \code{data.table}.
}
\description{
Get CNV Frequency Table
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_sig_db.R
\name{get_sig_db}
\alias{get_sig_db}
\title{Get Curated Reference Signature Database}
\usage{
get_sig_db(sig_db = "legacy")
}
\arguments{
\item{sig_db}{default 'legacy', it can be 'legacy' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/}{COSMIC v2 'SBS'}),
'SBS', 'DBS', 'ID' and 'TSB' (for \href{https://cancer.sanger.ac.uk/cosmic/signatures/}{COSMIV v3.1 signatures})
for small scale mutations.
For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
"SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
Degasperi, Andrea, et al. (2020) (reference #2);
"RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
"CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA.
\strong{UPDATE}, the latest version of reference version can be automatically
downloaded and loaded from \url{https://cancer.sanger.ac.uk/signatures/downloads/}
when a option with \code{latest_} prefix is specified (e.g. "latest_SBS_GRCh37").
\strong{Note}: the signature profile for different genome builds are basically same.
And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
For all available options, check the parameter setting.}
}
\value{
a \code{list}.
}
\description{
Reference mutational signatures and their aetiologies,
mainly obtained from COSMIC database
(SigProfiler results) and cleaned before saving into
\strong{sigminer} package. You can obtain:
\itemize{
\item COSMIC legacy SBS signatures.
\item COSMIC v3 SBS signatures.
\item COSMIC v3 DBS signatures.
\item COSMIC v3 ID (indel) signatures.
\item SBS and RS (rearrangement) signatures from Nik lab 2020 Nature Cancer paper.
\item RS signatures from BRCA560 and USARC cohorts.
\item Copy number signatures from USARC cohort and TCGA.
}
}
\examples{
s1 <- get_sig_db()
s2 <- get_sig_db("SBS")
s3 <- get_sig_db("DBS")
s4 <- get_sig_db("DBS_mm10")
s5 <- get_sig_db("SBS_Nik_lab")
s6 <- get_sig_db("ID")
s7 <- get_sig_db("RS_BRCA560")
s8 <- get_sig_db("RS_USARC")
s9 <- get_sig_db("RS_Nik_lab")
s10 <- get_sig_db("CNS_USARC")
s11 <- get_sig_db("CNS_TCGA")
s1
s2
s3
s4
s5
s6
s7
s8
s9
s10
s11
}
\references{
\itemize{
\item Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." bioRxiv (2021).
\item Alexandrov, Ludmil B., et al. "The repertoire of mutational signatures in human cancer." Nature 578.7793 (2020): 94-101.
\item Steele, Christopher D., et al. "Undifferentiated sarcomas develop through distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.
}
}
\seealso{
\link{get_sig_similarity}, \link{sig_fit} and \link{show_cosmic_sig_profile}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_sig_feature_corrplot.R
\name{show_sig_feature_corrplot}
\alias{show_sig_feature_corrplot}
\title{Draw Corrplot for Signature Exposures and Other Features}
\usage{
show_sig_feature_corrplot(
  tidy_cor,
  feature_list,
  sort_features = FALSE,
  sig_orders = NULL,
  drop = TRUE,
  return_plotlist = FALSE,
  p_val = 0.05,
  xlab = "Signatures",
  ylab = "Features",
  co_gradient_colors = scale_color_gradient2(low = "blue", mid = "white", high = "red",
    midpoint = 0),
  ca_gradient_colors = co_gradient_colors,
  plot_ratio = "auto",
  breaks_count = NULL
)
}
\arguments{
\item{tidy_cor}{data returned by \link{get_tidy_association}.}

\item{feature_list}{a character vector contains features want to be plotted.
If missing, all features will be used.}

\item{sort_features}{default is \code{FALSE}, use feature order obtained from the previous
step. If \code{TRUE}, sort features as \code{feature_list}.}

\item{sig_orders}{signature levels for ordering.}

\item{drop}{if \code{TRUE}, when a feature has no association with all signatures
(p value larger than threshold set by \code{p_val}), this feature will be removed
from the plot. Otherwise, this feature (a row) will keep with all blank white.}

\item{return_plotlist}{if \code{TRUE}, return as a list of \code{ggplot} objects.}

\item{p_val}{p value threshold. If p value larger than this threshold,
the result becomes blank white.}

\item{xlab}{label for x axis.}

\item{ylab}{label for y axis.}

\item{co_gradient_colors}{a Scale object representing gradient colors used to plot for continuous features.}

\item{ca_gradient_colors}{a Scale object representing gradient colors used to plot for categorical features.}

\item{plot_ratio}{a length-2 numeric vector to set the height/width ratio.}

\item{breaks_count}{breaks for sample count. If set it to \code{NULL},
ggplot \code{bin} scale will be used to automatically determine the
breaks. If set it to \code{NA}, \code{aes} for sample will be not used.}
}
\value{
a \code{ggplot2} object
}
\description{
This function is for association visualization. Of note,
the parameters \code{p_val} and \code{drop} will affect the visualization
of association results under p value threshold.
}
\examples{

# The data is generated from Wang, Shixiang et al.
load(system.file("extdata", "asso_data.RData",
  package = "sigminer", mustWork = TRUE
))

p <- show_sig_feature_corrplot(
            tidy_data.seqz.feature,
            p_val = 0.05,
            breaks_count = c(0L,200L, 400L, 600L, 800L, 1020L))
p
}
\seealso{
\link{get_tidy_association} and \link{get_sig_feature_association}
}
