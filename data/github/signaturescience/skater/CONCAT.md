
<!-- README.md is generated from README.Rmd. Please edit that file -->

# skater

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/skater)](https://CRAN.R-project.org/package=skater)

[![biorXiv](https://img.shields.io/badge/biorXiv-10.1101%2F2021.07.21.453083-red)](https://www.biorxiv.org/content/10.1101/2021.07.21.453083v1)

[![DOI](https://zenodo.org/badge/339462170.svg)](https://zenodo.org/badge/latestdoi/339462170)

[![R-CMD-check-stable](https://github.com/signaturescience/skater/workflows/R-CMD-check-stable/badge.svg)](https://github.com/signaturescience/skater/actions)

[![R-CMD-check-dev](https://github.com/signaturescience/skater/workflows/R-CMD-check-dev/badge.svg)](https://github.com/signaturescience/skater/actions)

<!-- badges: end -->

**S**NP-based **K**inship **A**nalysis **T**esting and **E**valuation:
miscellaneous **R** data analysis utilties.

## Installation

Install stable release from CRAN:

``` r
install.packages("skater")
```

Install development version from GitHub:

``` r
remotes::install_github("signaturescience/skater", build_vignettes=TRUE)
```

## Usage

The [“Basic Usage”
vignette](https://signaturescience.github.io/skater/articles/basic_usage.html)
steps through the primary functionality of the `skater` package:

``` r
vignette("basic_usage", package = "skater")
```

Full documentation: <https://signaturescience.github.io/skater/>.
# skater 0.1.1

- Updated vignette to handle future changes in tidyr (thanks @DavisVaughan, #56)
- Updated DESCRIPTION with `URL` and `BugReports` links.

# skater 0.1.0

- Initial release
## Test environments

- Local MacOS install, R 4.0.4
- R hub
    - Fedora Linux, R-devel
    - Ubuntu Linux 20.04.1 LTS, R-release
    - Windows Server 2008 R2 SP1, R-devel

## R CMD check results

- Local `R CMD check`: Status OK, 0 errors, 0 warnings, 0 notes
- R hub: 
    - NOTE, New submission
    - NOTE possibly mis-spelled words in description: IBD, et, al, benchmarking, polymorphism. IBD is defined in the DESCRIPTION as "identical by descent"; benchmarking and polymorphism are spelled correctly; and "et al." is used in a reference before linking to the doi with `<doi:...>`.

## Revisions after initial CRAN inspection

- Added more detailed description about package functionality in DESCRIPTION.
- Defined acronyms in DESCRIPTION.
- Added Signature Science, LLC as `cph` in DESCRIPTION Authors.
- Added reference to Description field of DESCRIPTION in the form: `authors (year) <doi:...>` with reference to preprint describing methods.
- Better explanation of identical by descent (IBD) segment to kinship coefficient math in function documentation.
- Stopped exporting two internal functions, removed examples, clarified documentation.
- Added a return value for `plot_pedigree()` (called for side effects).
- Updated exported functions to ensure `@return` `\value` notes class of the output value and what it means.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# skater

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/skater)](https://CRAN.R-project.org/package=skater)

[![biorXiv](https://img.shields.io/badge/biorXiv-10.1101%2F2021.07.21.453083-red)](https://www.biorxiv.org/content/10.1101/2021.07.21.453083v1)

[![DOI](https://zenodo.org/badge/339462170.svg)](https://zenodo.org/badge/latestdoi/339462170)

[![R-CMD-check-stable](https://github.com/signaturescience/skater/workflows/R-CMD-check-stable/badge.svg)](https://github.com/signaturescience/skater/actions)

[![R-CMD-check-dev](https://github.com/signaturescience/skater/workflows/R-CMD-check-dev/badge.svg)](https://github.com/signaturescience/skater/actions)

<!-- badges: end -->

**S**NP-based **K**inship **A**nalysis **T**esting and **E**valuation: miscellaneous **R** data analysis utilties.

## Installation

Install stable release from CRAN:

``` r
install.packages("skater")
```

Install development version from GitHub:

``` r
remotes::install_github("signaturescience/skater", build_vignettes=TRUE)
```


## Usage

The ["Basic Usage" vignette](https://signaturescience.github.io/skater/articles/basic_usage.html) steps through the primary functionality of the `skater` package:

``` r
vignette("basic_usage", package = "skater")
```

Full documentation: <https://signaturescience.github.io/skater/>.
---
title: '**skater**: An R package for SNP-based Kinship Analysis, Testing, and Evaluation'
author:
- name: Stephen D. Turner
  affiliation: Signature Science, LLC., Austin, TX 78759, USA.
- name: V. P. Nagraj
  affiliation: Signature Science, LLC., Austin, TX 78759, USA.
- name: Matthew Scholz
  affiliation: Signature Science, LLC., Austin, TX 78759, USA.
- name: Shakeel Jessa
  affiliation: Signature Science, LLC., Austin, TX 78759, USA.
- name: Carlos Acevedo
  affiliation: Signature Science, LLC., Austin, TX 78759, USA.
- name: Jianye Ge
  affiliation: Center for Human Identification, Department of Microbiology, Immunology,
    and Genetics, University of North Texas Health Science Center, Fort Worth, TX
    76107, USA.
- name: August E. Woerner
  affiliation: Center for Human Identification, Department of Microbiology, Immunology,
    and Genetics, University of North Texas Health Science Center, Fort Worth, TX
    76107, USA.
- name: Bruce Budowle
  affiliation: Center for Human Identification, Department of Microbiology, Immunology,
    and Genetics, University of North Texas Health Science Center, Fort Worth, TX
    76107, USA.
abstract: | 
    | 
    | **Motivation:** SNP-based kinship analysis with genome-wide relationship estimation and IBD segment analysis methods produces results that often require further downstream processing and manipulation. A dedicated software package that consistently and intuitively implements this analysis functionality is needed.  
    | **Results:**  Here we present the skater R package for **S**NP-based **k**inship **a**nalysis, **t**esting, and **e**valuation with **R**. The skater package contains a suite of well-documented tools for importing, parsing, and analyzing pedigree data, performing relationship degree inference, benchmarking relationship degree classification, and summarizing IBD segment data.  
    | **Availability:** The skater package is implemented as an R package and is released under the MIT license at https\:\/\/github.com/signaturescience/skater. Documentation is available at https\:\/\/signaturescience.github.io/skater.
output:
  BiocWorkflowTools::f1000_article: default
bibliography: bibliography.bib
keywords: bioinformatics, kinship, R, genealogy, SNPs, single nucleotide polymorphisms, relatedness
---

**R version**: `r R.version.string`

<!-- **Bioconductor version**: `r BiocManager::version()` -->

<!--  Update the name of the package below to be the name of the workflow package you are working on  -->

**skater package version**: `r packageVersion("skater")`

----

# Introduction

Inferring familial relationships between individuals using genetic data is a common practice in population genetics, medical genetics, and forensics. There are multiple approaches to estimating relatedness between samples, including genome-wide measures, such as those implemented in Plink [@purcell2007] or KING [@manichaikul2010], and methods that rely on identity by descent (IBD) segment detection, such as GERMLINE [@gusev2009], hap-IBD [@zhou2020], and IBIS [@seidman2020]. Recent efforts focusing on benchmarking these methods [@ramstetter2017] have been aided by tools for simulating pedigrees and genome-wide SNP data [@caballero2019]. Analyzing results from genome-wide SNP-based kinship analysis or comparing analyses to simulated data for benchmarking have to this point required writing one-off analysis functions or utility scripts that are seldom distributed with robust documentation, test suites, or narrative examples of usage. There is a need in the field for a well-documented software package with a consistent design and API that contains functions to assist with downstream manipulation, benchmarking, and analysis of SNP-based kinship assessment methods. Here we present the skater package for **S**NP-based **k**inship **a**nalysis, **t**esting, and **e**valuation with **R**.

# Methods

## Implementation

<!-- For software tool papers, this section should address how the tool works and any relevant technical details required for implementation of the tool by other developers. -->

The skater package provides an intuitive collection of analysis and utility functions for SNP-based kinship analysis. Functions in the package include tools for importing, parsing, and analyzing pedigree data, performing relationship degree inference, benchmarking relationship degree classification, and summarizing IBD segment data, described in full in the _Use Cases_ section below. The package adheres to "tidy" data analysis principles, and builds upon the tools released under the tidyverse R ecosystem [@Wickham2019].

The skater package is hosted in the Comprehensive R Archive Network (CRAN) which is the main repository for R packages: <http://CRAN.R-project.org/package=skater>. Users can install skater in R by executing the following code:

```{r, eval=FALSE}
install.packages("skater")
```

Alternatively, the development version of skater is available on GitHub at <https://github.com/signaturescience/skater>. The development version may contain new features which are not yet available in the version hosted on CRAN. This version can be installed using the `install_github()` function in the devtools package:

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("signaturescience/skater", build_vignettes=TRUE)
```

When installing skater, other packages which skater depends on are automatically installed, including magritr, tibble, dplyr, tidyr, readr, purrr, kinship2, corrr, rlang, and others.

## Operation

<!-- This part of the methods should include the minimal system requirements needed to run the software and an overview of the workflow for the tool for users of the tool. -->

Minimal system requirements for installing and using skater include R (version 3.0.0 or higher) and several tidyverse packages [@Wickham2019] that many R users will already have installed. Use cases are demonstrated in detail below. In summary, the skater package has functions for:

- Reading in various output files produced by commonly used tools in SNP-based kinship analysis
- Pedigree parsing, manpulation, and analysis
- Relationship degree inference
- Benchmarking and assessing relationship classification accuracy
- IBD segment analysis post-processing

A comprehensive reference for all the functions in the skater package is available at <https://signaturescience.github.io/skater/>. 

<!-- # Results  -->

<!-- Optional - only if novel data or analyses are included -->

<!-- This section is only required if the paper includes novel data or analyses, and should be written as a traditional results section. -->


# Use Cases 

The `skater` package provides a collection of analysis and utility functions for **S**NP-based **k**inship **a**nalysis, **t**esting, and **e**valuation as an **R** package. Functions in the package include tools for working with pedigree data, performing relationship degree inference, assessing classification accuracy, and summarizing IBD segment data.

```{r}
library(skater)
```

<!-- Optional - only if NO new datasets are included -->

<!-- This section is required if the paper does not include novel data or analyses.  -->
<!-- Examples of input and output files should be provided with some explanatory context.  Any novel or complex variable parameters should also be explained in sufficient detail to allow users to understand and use the tool's functionality. -->

## Pedigree parsing, manipulation, and analysis

<!-- The commented line below was used in the shorter JOSS/Bioinformatics submission, expanded here -->
<!-- The skater package has several functions for importing, parsing, and analyzing pedigree data. Pedigrees define familial relationships in a hierarchical structure. Many genomics tools for working with pedigrees start with a .fam file, which is a tabular format with one row per individual and columns for unique IDs of the mother, father, and the family unit. The skater package contains the function `read_fam()` to read in a PLINK-formatted .fam file and another function `fam2ped()` to convert the content into a pedigree object as a nested tibble with one row per family. All pedigree processing from skater internally leverages a data structure from the kinship2 package [@sinnwell2014]. Further functions such as `plot_pedigree()` produce a multi-page PDF drawing a diagram of the pedigree for each family, while `ped2kinpair()` produces a pairwise list of relationships between all individuals in the data with the expected kinship coefficients for each pair (see skater package vignette). -->

Pedigrees define familial relationships in a hierarchical structure. One of the common formats used by PLINK [@purcell2007] and other genetic analysis tools is the `.fam` file. A `.fam` file is a tabular format with one row per individual and columns for unique IDs of the mother, father, and the family unit. The package includes `read_fam()` to read files in this format:

```{r}
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
fam <- read_fam(famfile)
fam
```

Family structures imported from `.fam` formated files can then be translated to the `pedigree` structure used by the `kinship2` package [@sinnwell2014]. The "fam" format may include multiple families, and the `fam2ped()` function will collapse them all into a `tibble` with one row per family:

```{r}
peds <- fam2ped(fam)
peds
```

In the example above, the resulting `tibble` is nested by family ID. The `data` column contains the individual family information, while the `ped` column contains the pedigree object for that family. Using standard tidyverse operations, the resulting tibble can be unnested for any particular family:

```{r}
peds %>% 
  dplyr::filter(fid=="testped1") %>% 
  tidyr::unnest(cols=data)
```

A single pedigree can also be inspected or visualized (standard base R plot arguments such as `mar` or `cex` can be used to adjust aesthetics):

```{r}
peds$ped[[1]]
```

```{r plotped, fig.width=6, fig.height=3, fig.cap="Pedigree diagram for the first family in the pedigree shown in the code above."}
plot(peds$ped[[1]], mar=c(1,4,1,4), cex=.7)
```

The `plot_pedigree()` function from `skater` will iterate over a list of pedigree objects, writing a multi-page PDF, with each page containing a pedigree from family:

```{r, eval=FALSE}
plot_pedigree(peds$ped, file="3gens.ped.pdf")
```

The `ped2kinpair()` function takes a pedigree object and produces a pairwise list of relationships between all individuals in the data with the expected kinship coefficients for each pair.

The function can be run on a single family:

```{r}
ped2kinpair(peds$ped[[1]])
```

This function can also be mapped over all families in the pedigree:

```{r}
kinpairs <- 
  peds %>% 
  dplyr::mutate(pairs=purrr::map(ped, ped2kinpair)) %>% 
  dplyr::select(fid, pairs) %>% 
  tidyr::unnest(cols=pairs)
kinpairs
```

Note that this maps `ped2kinpair()` over all `ped` objects in the input `tibble`, and that relationships are not shown for between-family relationships.

## Relationship degree inference and benchmarking

The skater package includes functions to translate kinship coefficients to relationship degrees. The kinship coefficients could come from `ped2kinpair()` or other kinship estimation software.

<!-- The commented line below was used in the shorter JOSS/Bioinformatics submission, expanded here -->
<!-- The `dibble()` function creates a **d**egree **i**nference t**ibble**, with degrees up to the specified maximum degree resolution, expected kinship coefficient, and lower and upper inference ranges as defined in @manichaikul2010. The `kin2degree()` function infers the relationship degree given a kinship coefficient and a maximum degree resolution (e.g., 7th-degree relatives) up to which anything more distant is classified as unrelated. -->

The `dibble()` function creates a **d**egree **i**nference `tibble`, with degrees up to the specified `max_degree` (default=3), expected kinship coefficient, and lower (`l`) and upper (`u`) inference ranges as defined in @manichaikul2010. Degree 0 corresponds to self / identity / monozygotic twins, with an expected kinship coefficient of 0.5, with inference range >=0.354. Anything beyond the maximum degree resolution is considered unrelated (degree `NA`). Note also that while the theoretical upper boundary for the kinship coefficient is 0.5, the inference range for 0-degree (same person or identical twins) extends to 1 to allow for floating point arithmetic and stochastic effects resulting in kinship coefficients above 0.5.

```{r}
dibble()
```

The degree inference `max_degree` default is 3. Change this argument to allow more granular degree inference ranges:

```{r}
dibble(max_degree = 5)
```

Note that the distance between relationship degrees becomes smaller as the relationship degree becomes more distant. The `dibble()` function will emit a warning with `max_degree` >=10, and will stop with an error at >=12.

The `kin2degree()` function infers the relationship degree given a kinship coefficient and a `max_degree` up to which anything more distant is treated as unrelated. Example first degree relative:

```{r}
kin2degree(.25, max_degree=3)
```

Example 4th degree relative, but using the default max_degree resolution of 3:

```{r}
kin2degree(.0312, max_degree=3)
```

Example 4th degree relative, but increasing the degree resolution:

```{r}
kin2degree(.0312, max_degree=5)
```

The `kin2degree()` function is vectorized over values of `k`, so it can be used inside of a `mutate` on a `tibble` of kinship coefficients:

```{r}
# Get two pairs from each type of relationship we have in kinpairs:
kinpairs_subset <- 
  kinpairs %>% 
  dplyr::group_by(k) %>% 
  dplyr::slice(1:2)
kinpairs_subset

# Infer degree out to third degree relatives:
kinpairs_subset %>% 
  dplyr::mutate(degree=kin2degree(k, max_degree=3))
```



## Benchmarking Degree Classification

<!-- The commented line below was used in the shorter JOSS/Bioinformatics submission, expanded here -->
<!-- Once estimated kinship is converted to degree, it may be of interest to compare the inferred degree to known degrees of relatedness. When aggregated over many relationships and inferences, this can help benchmark performance of a particular kinship analysis method. The skater package adapts a `confusion_matrix()` function from @clark2021 to provide standard contingency table metrics (e.g. sensitivity, specificity, PPV, precision, recall, F1, etc.) with a new reciprocal RMSE (R-RMSE) metric. The R-RMSE metric is defined more thoroughly in the skater package vignette and may be a preferable measure of classification accuracy when benchmarking relationship degree estimation. In many kinship benchmarking analyses, classification error is treated in a categorical manner (exact match plus or minus one degree), neglecting the true amount of sharing as a real number. Taking the reciprocal of the target and predicted degree in a typical RMSE calculation results in larger penalties for more egregious misclassifications (e.g., classifying a first-degree relative pair as second-degree) than misclassifications at more distant relationships (e.g., classifying a fourth-degree relative pair as fifth-degree). -->

Once estimated kinship is converted to degree, it may be of interest to compare the inferred degree to truth. When aggregated over many relationships and inferences, this approach can help benchmark performance of a particular kinship analysis method.

The skater package adapts a `confusion_matrix()` function from @clark2021 to provide standard contingency table metrics (e.g. sensitivity, specificity, PPV, precision, recall, F1, etc.) with a new reciprocal RMSE (R-RMSE) metric. The `confusion_matrix()` function on its own outputs a list with four objects:

1. A `tibble` with calculated accuracy, lower and upper bounds, the guessing rate and p-value of the accuracy vs. the guessing rate. 
2. A `tibble` with contingency table statistics calculated for each class. Details on the statistics calculated for each class can be reviewed on the help page for `?confusion_matrix`.
3. A `matrix` with the contingency table object itself.
4. A `vector` with the reciprocal RMSE (R-RMSE). The R-RMSE represents an alternative to classification accuracy when benchmarking relationship degree estimation and is calculated using the formula in (1). Taking the reciprocal of the target and predicted degree results in larger penalties for more egregious misclassifications (e.g., classifying a first-degree relative pair as second degree) than misclassifications at more distant relationships (e.g., misclassifying a fourth-degree relative pair as fifth-degree). The +0.5 adjustment prevents division-by-zero when a 0th-degree (identical) relative pair is introduced.

\begin{equation}
\sqrt{\frac{\sum_{i=1}^{k}(\frac{1}{\text{Target}+0.5}-\frac{1}{\text{Predicted}+0.5})^2}{k}}
\end{equation}

To illustrate the usage, this example will start with the `kinpairs` data from above and randomly flip ~20% of the true relationship degrees:

```{r}
# Function to randomly flip levels of a factor (at 20%, by default)
randomflip <- function(x, p=.2) ifelse(runif(length(x))<p, sample(unique(x)), x)

# Infer degree (truth/target) using kin2degree, then randomly flip 20% of them
set.seed(42)
kinpairs_inferred <- kinpairs %>% 
  dplyr::mutate(degree_truth=kin2degree(k, max_degree=3)) %>% 
  dplyr::mutate(degree_truth=tidyr::replace_na(degree_truth, "unrelated")) %>% 
  dplyr::mutate(degree_inferred=randomflip(degree_truth))
kinpairs_inferred
```

Next, running the `confusion_matrix()` function will return all four objects noted above:

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth)
```

Standard tidyverse functions such as `purrr::pluck()` can be used to isolate just the contingency table:

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth) %>% 
  purrr::pluck("Table")
```

The `confusion_matrix()` function includes an argument to output in a tidy (`longer=TRUE`) format, and the example below illustrates how to spread contingency table statistics by class:

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth, 
                 longer = TRUE) %>% 
  purrr::pluck("Other") %>% 
  tidyr::spread(Class, Value) %>% 
  dplyr::relocate(Average, .after=dplyr::last_col()) %>% 
  dplyr::mutate_if(rlang::is_double, signif, 2)
```

## IBD segment analysis

Tools such as hap-IBD [@zhou2020], and IBIS [@seidman2020] detect shared IBD segments between individuals. The skater package includes functionality to take those IBD segments, compute shared genomic centimorgan (cM) length, and converts that shared cM to a kinship coefficient. In addition to inferred segments, these functions can estimate "truth" kinship from simulated IBD segments [@caballero2019]. The `read_ibd()` function reads pairwise IBD segments from IBD inference tools and from simulated IBD segments. The `read_map()` function reads in genetic map in a standard format which is required to translate the total centimorgans shared IBD to a kinship coefficient using the `ibd2kin()` function. See `?read_ibd` and `?read_map` for additional details on expected format.

The `read_ibd()` function reads in the pairwise IBD segment format. Input to this function can either be inferred IBD segments from hap-IBD (`source="hapibd"`) or simulated segments (`source="pedsim"`). The first example below uses data in the `hap-ibd` output format:

```{r, message=FALSE}
hapibd_filepath <- system.file("extdata", "GBR.sim.ibd.gz", 
                               package="skater")
hapibd_seg <- read_ibd(hapibd_filepath, source = "hapibd")
hapibd_seg
```

In order to translate the shared genomic cM length to a kinship coefficient, a genetic map must first be read in with `read_map()`. Software for IBD segment inference and simulation requires a genetic map. The map loaded for kinship estimation should be the same one used for creating the shared IBD segment output. The example below uses a minimal genetic map that ships with `skater`:

```{r}
gmap_filepath <- system.file("extdata", "sexspec-avg-min.plink.map", 
                             package="skater")
gmap <- read_map(gmap_filepath)
gmap
```

The `ibd2kin()` function takes the segments and map file and outputs a `tibble` with one row per pair of individuals and columns for individual 1 ID, individual 2 ID, and the kinship coefficient for the pair:

```{r}
ibd_dat <- ibd2kin(.ibd_data=hapibd_seg, .map=gmap)
ibd_dat
```

<!-- # Discussion  -->

<!-- Optional - only if novel data or analyses are included -->

<!-- This section is only required if the paper includes novel data or analyses, and should be written in the same style as a traditional discussion section. -->

<!-- Please include a brief discussion of allowances made (if any) for controlling bias or unwanted sources of variability, and the limitations of any novel datasets. -->


<!-- # Conclusions  -->

<!-- Optional - only if novel data or analyses are included -->

<!-- This section is only required if the paper includes novel data or analyses, and should be written as a traditional conclusion. -->


# Summary 

<!-- Optional - only if NO new datasets are included -->

The skater R package provides a robust software package for data import, manipulation, and analysis tasks typically encountered when working with SNP-based kinship analysis tools. All package functions are internally documented with examples, and the package contains a vignette demonstrating usage, inputs, outputs, and interpretation of all key functions. The package contains internal tests that are automatically run with continuous integration via GitHub Actions whenever the package code is updated. The skater package is permissively licensed (MIT) and is easily extensible to accommodate outputs from new genome-wide relatedness and IBD segment methods as they become available.


<!-- # Data availability  -->

<!-- Optional - only if novel data or analyses are included -->

<!-- Please add details of where any datasets that are mentioned in the paper, and that have not have not previously been formally published, can be found.  If previously published datasets are mentioned, these should be cited in the references, as per usual scholarly conventions. -->


# Software availability

<!-- This section will be generated by the Editorial Office before publication. Authors are asked to provide some initial information to assist the Editorial Office, as detailed below. -->

<!-- 1. URL link to where the software can be downloaded from or used by a non-coder (AUTHOR TO PROVIDE; optional) -->
<!-- 2. URL link to the author's version control system repository containing the source code (AUTHOR TO PROVIDE; required) -->
<!-- 3. Link to source code as at time of publication (*F1000Research* TO GENERATE) -->
<!-- 4. Link to archived source code as at time of publication (*F1000Research* TO GENERATE) -->
<!-- 5. Software license (AUTHOR TO PROVIDE; required) -->

1. Software available from: <http://CRAN.R-project.org/package=skater>.
1. Source code available from: <https://github.com/signaturescience/skater>.
1. Archived source code at time of publication: <https://doi.org/10.5281/zenodo.5761996>.
1. Software license: MIT License.

# Author information

SDT, VPN, and MBS developed the R package.

All authors contributed to method development.

SDT wrote the first draft of the manuscript.

All authors assisted with manuscript revision.

All authors read and approved the final manuscript.


# Competing interests

<!-- All financial, personal, or professional competing interests for any of the authors that could be construed to unduly influence the content of the article must be disclosed and will be displayed alongside the article. If there are no relevant competing interests to declare, please add the following: 'No competing interests were disclosed'. -->

No competing interests were disclosed.

# Grant information

<!-- Please state who funded the work discussed in this article, whether it is your employer, a grant funder etc. Please do not list funding that you have that is not relevant to this specific piece of research. For each funder, please state the funder's name, the grant number where applicable, and the individual to whom the grant was assigned. If your work was not funded by any grants, please include the line: 'The author(s) declared that no grants were involved in supporting this work.' -->

This work was supported in part by award 2019-DU-BX-0046 (Dense DNA Data for Enhanced Missing Persons Identification) to B.B., awarded by the National Institute of Justice, Office of Justice Programs, U.S. Department of Justice and by internal funds from the Center for Human Identification. The opinions, findings, and conclusions or recommendations expressed are those of the authors and do not necessarily reflect those of the U.S. Department of Justice.


<!-- # Acknowledgments -->

<!-- This section should acknowledge anyone who contributed to the research or the article but who does not qualify as an author based on the criteria provided earlier (e.g. someone or an organization that provided writing assistance). Please state how they contributed; authors should obtain permission to acknowledge from all those mentioned in the Acknowledgments section. -->

<!-- Please do not list grant funding in this section. -->






<!--

# USING R MARKDOWN

Some examples of commonly used markdown syntax are listed below, to help you get started.

## Cross-references

For portability between different output formats, use the syntax introduced by *bookdown*, such as `(\#label)` for labels and `\@ref(label)` for cross-references. The following sections provide examples of referencing tables, figures, and equations.

## Citations

You can include references in a standard Bibtex file.  The name of this file is given in the header of the markdown document (in our case it is *sample.bib*).  References to entries in the Bibtex file are made using square brackets and use an @ plus the key for the entry you are referencing [@Smith:2012qr].  You can combine multiple entries by separating them with a semi-colon [@Smith:2012qr; @Smith:2013jd].
The default bibliography style uses numerical citations.  For superscript or author-year citations set the header metadata field `natbiboptions` to either `super` or `round`, respectively.

If you specify a figure caption to a code chunk using the chunk option `fig.cap`, the plot will be automatically labeled and numbered. The figure label is generated from the label of the code chunk by prefixing it with `fig:`, e.g., see Figure \@ref(fig:plot).

## Tables

Markdown syntax tends to lack some of the more sophisticated formatting features available in LaTeX, so you may need to edit the tables later to get the desired format.

| First name  | Last Name | Grade |
| ----------- | --------- | ----- |
| John        | Doe       |   7.5 |
| Richard     | Miles     |     2 |

Table: Caption to table.

Just like figures, tables with captions will also be numbered and can be referenced. Captions are entered as a paragraph beginning with the string "Table:" (or just ":"), which may appear either before or after the table. A label for the table should appear in the beginning of the caption in the form of `(\#tab:label)`, e.g., see Table \@ref(tab:table).

: (\#tab:table) A table with text justification.

| First name  | Last Name | Grade |
| ----------- | :-------: | ----: |
| John        | Doe       |   7.5 |
| Richard     | Miles     |     2 |

## Figures

You can include static figures (i.e. no generated by code) using the `include_graphics()` function from the **knitr** package, in a standard code chunk.

You can again use the `fig.cap` option to provide the figure caption, and reference the image based on the code chunk label.  You can also use options such as `fig.align` and `fig.width` to adjust the position and size of the image within the final document, e.g. Figure \@ref(fig:frog-picture) is a frog.

Alternatively, you can use the standard markdown syntax like so:

![This is a smaller version of the same picture, inserted using the standard markdown syntax](frog.jpg){width=2cm height=2cm}

Please give figures appropriate filenames, e.g.: figure1.pdf, figure2.png.

Figure legends should briefly describe the key messages of the figure such that the figure can stand alone from the main text. However, all figures should also be discussed in the article text. Each legend should have a concise title of no more than 15 words. The legend itself should be succinct, while still explaining all symbols and abbreviations. Avoid lengthy descriptions of methods.  

For any figures reproduced from another publication (as long as appropriate permission has been obtained from the copyright holder -see under the heading 'Submission'), please include a line in the legend to state that: 'This figure has been reproduced with kind permission from [include original publication citation]'.

-->
---
title: "Basic Usage"
output: 
  rmarkdown::html_vignette:
      toc: true
vignette: >
  %\VignetteIndexEntry{Basic Usage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)
options(rmarkdown.html_vignette.check_title = FALSE)
```


## Overview

The `skater` package provides a collection of analysis and utility functions for **S**NP-based **k**inship **a**nalysis, **t**esting, and **e**valuation as an **R** package. Functions in the package include tools for working with pedigree data, performing relationship degree inference, assessing classification accuracy, and summarizing IBD segment data.

```{r}
library(skater)
```

## Pedigree parsing and manipulation

Pedigrees define familial relationships in a hierarchical structure. 

One of the formats used by PLINK and other genetic analysis tools is the `.fam` file.[^plink-fam] A `.fam` file is a tabular format with one row per individual and columns for unique IDs of the mother, father, and the family unit. The package includes `read_fam()` to read files in this format:

[^plink-fam]: https://www.cog-genomics.org/plink/1.9/formats#fam

```{r}
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
fam <- read_fam(famfile)
fam
```

Family structures imported from ".fam" formated files can then be translated to the `pedigree` structure used by the `kinship2` package.[^kinship2-ref] The "fam" format may include multiple families, and the `fam2ped()` function will collapse them all into a `tibble` with one row per family:

[^kinship2-ref]: Sinnwell, Jason P., Terry M. Therneau, and Daniel J. Schaid. "The kinship2 R package for pedigree data." _Human heredity_ 78.2 (2014): 91-93.

```{r}
peds <- fam2ped(fam)
peds
```

In the example above, the resulting `tibble` is nested by family ID. The `data` column contains the individual family information, while the `ped` column contains the pedigree object for that family. You can unnest any particular family:

```{r}
peds %>% 
  dplyr::filter(fid=="testped1") %>% 
  tidyr::unnest(cols=data)
```

You can also look at a single pedigree:

```{r}
peds$ped[[1]]
```

Or plot that pedigree:

```{r plotped, fig.width=8, fig.height=8, fig.align="center"}
plot(peds$ped[[1]], mar=c(1,4,1,4))
```

The `plot_pedigree()` function from `skater` will iterate over a list of pedigree objects, writing a multi-page PDF, with each page containing a pedigree from family:

```{r, eval=FALSE}
plot_pedigree(peds$ped, file="3gens.ped.pdf")
```

The `ped2kinpair()` function takes a pedigree object and produces a pairwise list of relationships between all individuals in the data with the expected kinship coefficients for each pair.

The function can be run on a single family:

```{r}
ped2kinpair(peds$ped[[1]])
```

Or mapped over all families in the pedigree

```{r}
kinpairs <- 
  peds %>% 
  dplyr::mutate(pairs=purrr::map(ped, ped2kinpair)) %>% 
  dplyr::select(fid, pairs) %>% 
  tidyr::unnest(cols=pairs)
kinpairs
```

Note that this maps `ped2kinpair()` over all `ped` objects in the input `tibble`, and that relationships are not shown for between-family relationships (which should all be zero).

## Degree Inference

The `skater` package includes functions to translate kinship coefficients to relationship degrees. The kinship coefficients could come from `ped2kinpair()` or other kinship estimation software.

The `dibble()` function creates a **d**egree **i**nference `tibble`, with degrees up to the specified `max_degree` (default=3), expected kinship coefficient, and lower (`l`) and upper (`u`) inference ranges as defined in the KING paper.[^manichaikul] Degree 0 corresponds to self / identity / monozygotic twins, with an expected kinship coefficient of 0.5, with inference range >=0.354. Anything beyond the maximum degree resolution is considered unrelated (degree `NA`), with expected kinship coefficient of 0.

[^manichaikul]: Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W. M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics (Oxford, England), 26(22), 2867–2873. https://doi.org/10.1093/bioinformatics/btq559

```{r}
dibble()
```

The degree inference `max_degree` default is 3. Change this argument to allow more granular degree inference ranges:

```{r}
dibble(max_degree = 5)
```

Note that the distance between relationship degrees becomes smaller as the relationship degree becomes more distant. The `dibble()` function will throw a warning with `max_degree` >=10, and will stop with an error at >=12.

The `kin2degree()` function infers the relationship degree given a kinship coefficient and a `max_degree` up to which anything more distant is treated as unrelated. Example first degree relative:

```{r}
kin2degree(.25, max_degree=3)
```

Example 4th degree relative, but using the default max_degree resolution of 3:

```{r}
kin2degree(.0312, max_degree=3)
```

Example 4th degree relative, but increasing the degree resolution:

```{r}
kin2degree(.0312, max_degree=5)
```

The `kin2degree()` function is vectorized over values of `k`, so it can be used inside of a `mutate` on a `tibble` of kinship coefficients:

```{r}
# Get two pairs from each type of relationship we have in kinpairs:
kinpairs_subset <- 
  kinpairs %>% 
  dplyr::group_by(k) %>% 
  dplyr::slice(1:2)
kinpairs_subset

# Infer degree out to third degree relatives:
kinpairs_subset %>% 
  dplyr::mutate(degree=kin2degree(k, max_degree=3))
```

## Benchmarking Degree Classification

Once estimated kinship is converted to degree, it may be of interest to compare the inferred degree to truth. When aggregated over many relationships and inferences, this method can help benchmark performance of a particular kinship analysis method.

The `skater` package adapts functionality from the `confusionMatrix` package[^confusionMatrix] in the `confusion_matrix()` function.

[^confusionMatrix]: https://github.com/m-clark/confusionMatrix

The `confusion_matrix()` function on its own outputs a list with three objects:

1. A `tibble` with calculated accuracy, lower and upper bounds, the guessing rate and p-value of the accuracy vs. the guessing rate. 
2. A `tibble` with the following statistics (for each class):
    - Sensitivity = A/(A+C)
    - Specificity = D/(B+D)
    - Prevalence = (A+C)/(A+B+C+D)
    - PPV = (sensitivity * prevalence)/((sensitivity * prevalence) + ((1-specificity) * (1-prevalence)))
    - NPV = (specificity * (1-prevalence))/(((1-sensitivity) * prevalence) + ((specificity) * (1-prevalence)))
    - Detection Rate = A/(A+B+C+D)
    - Detection Prevalence = (A+B)/(A+B+C+D)
    - Balanced Accuracy = (sensitivity+specificity)/2
    - Precision = A/(A+B)
    - Recall = A/(A+C)
    - F1 = harmonic mean of precision and recall
    - False Discovery Rate = 1 - PPV
    - False Omission Rate = 1 - NPV
    - False Positive Rate = 1 - Specificity
    - False Negative Rate = 1 - Sensitivity
3. A `matrix` with the contingency table object itself.
4. A `vector` with the reciprocal RMSE (R-RMSE). The R-RMSE is calculated as `sqrt(mean((1/(Target+.5)-1/(Predicted+.5))^2)))`, and is a superior measure to classification accuracy when benchmarking relationship degree estimation. Taking the reciprocal of the target and predicted degree results in larger penalties for more egregious misclassifications (e.g., classifying a first-degree relative pair as second degree) than misclassifications at more distant relationships (e.g., misclassifying a fourth-degree relative pair as fifth-degree). The +0.5 adjustment prevents division-by-zero when a 0th-degree (identical) relative pair is introduced.

To illustrate the usage, first take the `kinpairs` data from above and randomly flip ~20% of the true relationship degrees.

```{r}
# Function to randomly flip levels of a factor (at 20%, by default)
randomflip <- function(x, p=.2) ifelse(runif(length(x))<p, sample(unique(x)), x)

# Infer degree (truth/target) using kin2degree, then randomly flip 20% of them
set.seed(42)
kinpairs_inferred <- kinpairs %>% 
  dplyr::mutate(degree_truth=kin2degree(k, max_degree=3)) %>% 
  dplyr::mutate(degree_truth=as.character(degree_truth)) %>%
  dplyr::mutate(degree_truth=tidyr::replace_na(degree_truth, "unrelated")) %>% 
  dplyr::mutate(degree_inferred=randomflip(degree_truth))
kinpairs_inferred
```

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth)
```

You can use `purrr::pluck()` to isolate just the contingency table:

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth) %>% 
  purrr::pluck("Table")
```

Or optionally output in a tidy (`longer=TRUE`) format, then spread stats by class:

```{r}
confusion_matrix(prediction = kinpairs_inferred$degree_inferred, 
                 target = kinpairs_inferred$degree_truth, 
                 longer = TRUE) %>% 
  purrr::pluck("Other") %>% 
  tidyr::spread(Class, Value) %>% 
  dplyr::relocate(Average, .after=dplyr::last_col()) %>% 
  dplyr::mutate_if(rlang::is_double, signif, 2) %>% 
  knitr::kable()
```

## IBD Segment Analysis

Tools such as `hap-ibd`[^hap-ibd] are capable of inferring shared IBD segments between individuals. The `skater` package includes functionality to take those IBD segments, compute shared genomic centimorgan (cM) length, and convert that shared cM to a kinship coefficient. In addition to inferred segments, these functions can estimate "truth" kinship from data simulated by `ped-sim`.[^ped-sim]

[^hap-ibd]: https://github.com/browning-lab/hap-ibd#output-files
[^ped-sim]: https://github.com/williamslab/ped-sim#output-ibd-segments-file

The `read_ibd()` function reads in the pairwise IBD segment format. Input to this function can either be inferred IBD segments from hap-IBD (`source="hapibd"`) or simulated segments (`source="pedsim"`). The first example below uses data in the `hap-ibd` output format:

```{r}
hapibd_fp <- system.file("extdata", "GBR.sim.ibd.gz", package="skater", mustWork=TRUE)
hapibd_seg <- read_ibd(hapibd_fp, source = "hapibd")
hapibd_seg
```

In order to translate the shared genomic cM length to a kinship coefficient, you must load a genetic map with `read_map()`. Software for IBD segment inference and simulation requires a genetic map. The map loaded for kinship estimation should be the same one used for creating the shared IBD segment output. The example below uses a minimal genetic map created with `min_map`[^min_map] that ships with `skater`:

[^min_map]: https://github.com/williamslab/min_map

```{r}
gmapfile <- system.file("extdata", "sexspec-avg-min.plink.map", package="skater", mustWork=TRUE)
gmap <- read_map(gmapfile)
gmap
```

The `ibd2kin()` function takes the segments and map file and outputs a `tibble` with one row per pair of individuals and columns for individual 1 ID, individual 2 ID, and the kinship coefficient for the pair:

```{r}
ibd_dat <- ibd2kin(.ibd_data=hapibd_seg, .map=gmap)
ibd_dat
```

As noted above, the IBD segment kinship estimation can be performed on simulated segments. The package includes an example of IBD data in that format:

```{r}
pedsim_fp <- system.file("extdata", "GBR.sim.seg.gz", package="skater", mustWork=TRUE)
pedsim_seg <- read_ibd(pedsim_fp, source = "pedsim")
pedsim_seg
```

Notably, `ped-sim` differentiates IBD1 and IBD2 segments. Given that IBD1 and IBD2 segments are weighted differently in kinship calculation, this should be accounted for in processing. In the example below the shared IBD is calculated separately for IBD1 and IBD2 with `type="IBD1"` and `type="IBD2"` respectively. You can then combine those results and sum the IBD1 and IBD2 kinship coefficients to get the overall kinship coefficient: 

```{r}
ibd1_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD1, .map=gmap, type="IBD1")
ibd2_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD2, .map=gmap, type="IBD2")
dplyr::bind_rows(ibd1_dat,ibd2_dat) %>%
  dplyr::group_by(id1,id2) %>%
  dplyr::summarise(kinship = sum(kinship), .groups = "drop")
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibd.R
\name{ibd2kin}
\alias{ibd2kin}
\title{Compute kinship coefficient from IBD segments}
\usage{
ibd2kin(.ibd_data, .map, type = NULL)
}
\arguments{
\item{.ibd_data}{Tibble with IBD segments created using the \link[skater]{read_ibd} function}

\item{.map}{Tibble with the genetic map data created using the \link[skater]{read_map} function}

\item{type}{Type of IBD to use for kinship coefficient calculation; must be \code{'IBD1'}, \code{'IBD2'}, or \code{NULL} (both IBD1 and IBD2 will be treated the same); default is \code{NULL}}
}
\value{
Tibble with three columns:
\enumerate{
\item id1 (sample identifier 1)
\item id2 (sample identifier 2)
\item kinship (kinship coefficent derived from shared segments)
}
}
\description{
This function is used to retrieve a relatedness measure from IBD segments.
The relatedness value returned is the kinship coefficient.
}
\details{
The input data should be pairwise IBD segments prepared via
\link[skater]{read_ibd}. The function will internally loop over each
chromosome, and use a specified genetic map to convert shared segments to
genetic units. After doing so, the function converts the shared length to a
kinship coefficient by summing \eqn{0.5*IBD2 + 0.25*IBD1}.

Note that the data read in by \link{read_ibd} when \code{source="pedsim"} returns a
list with separate tibbles for IBD1 and IBD2 segments. The current
implementation of this function requires running this function independently
on IBD1 and IBD2 segments, then summarizing (adding) the corresponding
proportions. See examples.
}
\examples{
pedsim_fp <- system.file("extdata", "GBR.sim.seg.gz", package="skater", mustWork=TRUE)
pedsim_seg <- read_ibd(pedsim_fp, source = "pedsim")
gmapfile <- system.file("extdata", "sexspec-avg-min.plink.map", package="skater", mustWork=TRUE)
gmap <- read_map(gmapfile)
ibd1_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD1, .map=gmap, type="IBD1")
ibd2_dat <- ibd2kin(.ibd_data=pedsim_seg$IBD2, .map=gmap, type="IBD2")
dplyr::bind_rows(ibd1_dat,ibd2_dat) \%>\%
  dplyr::group_by(id1,id2) \%>\%
  dplyr::summarise(kinship = sum(kinship), .groups = "drop")

}
\references{
http://faculty.washington.edu/sguy/ibd_relatedness.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{arrange_ids}
\alias{arrange_ids}
\title{Order IDs across two columns}
\usage{
arrange_ids(.data, .id1, .id2)
}
\arguments{
\item{.data}{A tibble with two ID columns to arrange.}

\item{.id1}{Unquoted name of the "id1" column. See examples.}

\item{.id2}{Unquoted name of the "id2" column. See examples.}
}
\value{
A tibble with id1 and id2 rearranged alphanumerically.
}
\description{
Some types of data or results are indexed by two identifiers in
two different columns corresponding to data points for \emph{pairs} of
observations. E.g., you may have columns called \code{id1} and \code{id2} that index
the tibble for all possible pairs of results between samples A, B, and C.
If you attempt to join two tibbles with \code{by=c("id1", "id2")}, the join will
fail if samples are flipped from one dataset to another. E.g., one tibble
may have id1=A and id2=B while the other has id1=B and id2=A. This function
ensures that id1 is alphanumerically first while id2 is alphanumerically
second. See examples.
}
\examples{
d1 <- tibble::tribble(
  ~id1, ~id2, ~results1,
  "a",  "b",       10L,
  "a",  "c",       20L,
  "c",  "b",       30L
)
d2 <- tibble::tribble(
  ~id1, ~id2,  ~results2,
  "b",  "a",       101L,
  "c",  "a",       201L,
  "b",  "c",       301L
)
# Inner join fails because id1!=id2.
dplyr::inner_join(d1, d2, by=c("id1", "id2"))
# Arrange IDs
d1 \%>\% arrange_ids(id1, id2)
d2 \%>\% arrange_ids(id1, id2)
# Inner join
dplyr::inner_join(arrange_ids(d1, id1, id2), arrange_ids(d2, id1, id2), by=c("id1", "id2"))
# Recursively, if you had more than two tibbles
list(d1, d2) \%>\%
  purrr::map(arrange_ids, id1, id2) \%>\%
  purrr::reduce(dplyr::inner_join, by=c("id1", "id2"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ibd.R
\name{interpolate}
\alias{interpolate}
\title{Interpolate over segments}
\usage{
interpolate(ibd_bp, chromgpos)
}
\arguments{
\item{ibd_bp}{Base pair for the IBD segment over which to interpolate}

\item{chromgpos}{Genetic map data for a specific chromosome}
}
\value{
Numeric vector with the genetic distance shared at the segment.
}
\description{
This is an unexported helper used in in \link[skater]{ibd2kin}. The function interpolates over segments to apply genetic length to the segment. It is inspired by Python code distributed by the Browning lab (\href{http://faculty.washington.edu/sguy/ibd_relatedness.html}{documentation}).
}
\references{
http://faculty.washington.edu/sguy/ibd_relatedness.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_fam}
\alias{read_fam}
\title{Read PLINK-formatted .fam file}
\usage{
read_fam(file)
}
\arguments{
\item{file}{Input file path}
}
\value{
A tibble containing the 6 columns from the fam file.
}
\description{
Reads in a \href{https://www.cog-genomics.org/plink/1.9/formats#fam}{PLINK-formatted .fam file}. Input \code{file} must have six columns:
\enumerate{
\item Family ID
\item Individual ID
\item Father ID
\item Mother ID
\item Sex
\item Affected Status
}
}
\examples{
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
fam <- read_fam(famfile)
fam

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confusion_matrix.R
\name{confusion_matrix}
\alias{confusion_matrix}
\title{Calculate various statistics from a confusion matrix}
\usage{
confusion_matrix(
  prediction,
  target,
  positive = NULL,
  prevalence = NULL,
  dnn = c("Predicted", "Target"),
  longer = FALSE,
  ...
)
}
\arguments{
\item{prediction}{A vector of predictions}

\item{target}{A vector of target values}

\item{positive}{The positive class for a 2-class setting. Default is
\code{NULL}, which will result in using the first level of \code{target}.}

\item{prevalence}{Prevalence rate.  Default is \code{NULL}.}

\item{dnn}{The row and column headers for the contingency table returned.
Default is 'Predicted' for rows and 'Target' for columns.}

\item{longer}{Transpose the output to long form.  Default is FALSE (requires
\code{tidyr 1.0}).}

\item{...}{Other parameters, not currently used.}
}
\value{
A list of tibble(s) with the associated statistics and possibly the
frequency table as list column of the first element. If classes contain >1
numeric class and a single non-numeric class (e.g., "1", "2", "3", and
"Unrelated", the RMSE of the reciprocal of the Targets + 0.5 will also be
returned.)
}
\description{
Given a vector of predictions and target values, calculate
numerous statistics of interest. Modified from
\href{https://github.com/m-clark/confusionMatrix}{m-clark/confusion_matrix}.
}
\details{
This returns accuracy, agreement, and other statistics. See the
functions below to find out more. Originally inspired by the
\code{confusionMatrix} function from the \code{caret} package.
}
\examples{
prediction = c(0,1,1,0,0,1,0,1,1,1)
target     = c(0,1,1,1,0,1,0,1,0,1)
confusion_matrix(prediction, target, positive = '1')

set.seed(42)
prediction = sample(letters[1:4], 250, replace = TRUE, prob = 1:4)
target     = sample(letters[1:4], 250, replace = TRUE, prob = 1:4)
confusion_matrix(prediction, target)

prediction = c(rep(1, 50), rep(2, 40), rep(3, 60))
target     = c(rep(1, 50), rep(2, 50), rep(3, 50))
confusion_matrix(prediction, target)
confusion_matrix(prediction, target) \%>\% purrr::pluck("Table")
confusion_matrix(prediction, target, longer=TRUE)
confusion_matrix(prediction, target, longer=TRUE) \%>\%
  purrr::pluck("Other") \%>\%
  tidyr::spread(Class, Value)

# Prediction with an unrelated class
prediction = c(rep(1, 50), rep(2, 40), rep(3, 60), rep("Unrelated", 55))
target     = c(rep(1, 50), rep(2, 50), rep(3, 55), rep("Unrelated", 50))
confusion_matrix(prediction, target)
# Prediction with two unrelated classes
prediction = c(rep(1, 50), rep(2, 40), rep("Third", 60), rep("Unrelated", 55))
target     = c(rep(1, 50), rep(2, 50), rep("Third", 55), rep("Unrelated", 50))
confusion_matrix(prediction, target)

}
\references{
Kuhn, M., & Johnson, K. (2013). Applied predictive modeling.
}
\seealso{
\code{\link{calc_accuracy}} \code{\link{calc_stats}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ped.R
\name{plot_pedigree}
\alias{plot_pedigree}
\title{Plot pedigree}
\usage{
plot_pedigree(ped, file = NULL, width = 10, height = 8)
}
\arguments{
\item{ped}{List of pedigree objects from \link{fam2ped}}

\item{file}{Output file path (must end in ".pdf")}

\item{width}{Width of output PDF}

\item{height}{Height of output PDF}
}
\value{
No return value, called for side effects.
}
\description{
Plot pedigree
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_ibis}
\alias{read_ibis}
\title{Read IBIS coef output file}
\usage{
read_ibis(file)
}
\arguments{
\item{file}{Input file path}
}
\value{
A tibble containing the 6 columns from the ibis file.
}
\description{
Reads in an \code{ibis} \href{https://github.com/williamslab/ibis}{results file}. Input \code{file} must have six columns, whitespace delimited:
\enumerate{
\item id1 (member 1)
\item id2 (member 2)
\item Kinship Coefficient
\item IBD2 (ratio of IBD2/All SNPS)
\item Segment count
\item Kinship Degree
}
}
\examples{
ibisFile <- system.file("extdata", "3gens.ibis.coef", package="skater", mustWork=TRUE)
ibis <- read_ibis(ibisFile)
ibis

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_map}
\alias{read_map}
\title{Read genetic map file}
\usage{
read_map(file)
}
\arguments{
\item{file}{Input file path}
}
\value{
A tibble containing 3 columns:
\enumerate{
\item chr (chromosome)
\item value (genetic length within the physical position boundary)
\item bp (physical position boundary)
}
}
\description{
This function reads in the content from a genetic map file to translate physical distance to genetic units (i.e. cM). Regardless of the source, the input file must be sex-averaged and in a tab-separated "Plink" format (\href{http://zzz.bwh.harvard.edu/plink/data.shtml#map}{documentation}) with the following four columns and no header (i.e. no column names):
\enumerate{
\item Chromosome
\item Identifier (ignored in \code{read_map()})
\item Length (genetic length within the physical position boundary)
\item Position (physical position boundary)
}

The columns must be in the order above. Note that only the first, third, and fourth columns are used in the function.
}
\details{
The genetic map could come from different sources. One source is the HapMap map distributed by the Browning Lab (\href{http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/}{documentation}). If this map file is used, the non-sex chromosomes can be downloaded and concatenated to a single file as follows:\preformatted{wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
unzip plink.GRCh37.map.zip
cat *chr[0-9]*GRCh37.map | sort -k1,1 -k4,4 --numeric-sort > plink.allchr.GRCh37.map
}

Another source is a sex-specific map ("bherer") originally published by Bherer et al and recommended by the developers of \code{ped-sim} for simulating IBD segments (\href{https://github.com/williamslab/ped-sim#map-file}{documentation}). To retrieve and prep this map file for simulation:\preformatted{# Get the refined genetic map and extract
wget --no-check-certificate https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination/raw/master/Refined_genetic_map_b37.tar.gz
tar xvfpz Refined_genetic_map_b37.tar.gz

# Format for ped-sim as per https://github.com/williamslab/ped-sim#map-file-
printf "#chr\\tpos\\tmale_cM\\tfemale_cM\\n" > sexspec.pedsim.map
for chr in \{1..22\}; do
  paste Refined_genetic_map_b37/male_chr$chr.txt Refined_genetic_map_b37/female_chr$chr.txt \\
    | awk -v OFS="\\t" 'NR > 1 && $2 == $6 \{print $1,$2,$4,$8\}' \\
    | sed 's/^chr//' >> sexspec.pedsim.map;
done

# Clean up
rm -rf Refined_genetic_map_b37*
}

After this, the \code{sexspec.pedsim.map} file is ready for use in simulation. However, it must be averaged and reformatted to "Plink format" to use here:\preformatted{cat sexspec.pedsim.map | grep -v "^#" | awk -v OFS="\\t" '\{print $1,".",($3+$4)/2,$2\}' > sexspec-avg.plink.map
}

#' The genetic maps created above are in the tens of megabytes size range. This is trivial to store for most systems but a reduced version would increase portability and ease testing. This "minimum viable genetic map" could be used for testing and as installed package data in an R package for example analysis. Read more about minimum viable genetic maps at:
\itemize{
\item Blog post: \url{https://hapi-dna.org/2020/11/minimal-viable-genetic-maps/}
\item Github repo with python code: \url{https://github.com/williamslab/min_map}
}

The code as written below reduces the averaged sex-specific genetic map from 833776 to 28726 positions (~30X reduction!).\preformatted{# Get minmap script from github
wget https://raw.githubusercontent.com/williamslab/min_map/main/min_map.py

# Create empty minmap
echo -n > sexspec-avg-min.plink.map

# For each autosome...
for chr in \{1..22\}; do
  echo "Working on chromosome $chr..."
  # First pull out just one chromosome
  grep "^$\{chr\}[[:space:]]" sexspec-avg.plink.map > tmp.$\{chr\}
  # Run the python script on that chromosome.
  # The genetic map column is 3rd column (2nd in 0-start). Physical position is last column (3 in 0-based)
  python3 min_map.py -mapfile tmp.$\{chr\} -chr $\{chr\} -genetcol 2 -physcol 3 -noheader -error 0.05
  # Strip out the header and reformat back to plink format, and append to minmap file
  cat min_viable_map$\{chr\}.txt | grep -v "^#" | awk -v OFS="\\t" '\{print $1,".",$4,$2\}' >> sexspec-avg-min.plink.map
  # Clean up
  rm -f min_viable_map$\{chr\}.txt tmp.$\{chr\}
done
}

This averaged version of the Bherer sex-specific map, reduced to a minimum viable genetic map with at most 5\% error, in Plink format, is available as installed package data (see examples). This is useful for testing code, but the full genetic map should be used for most analysis operations.
}
\examples{
gmapfile <- system.file("extdata", "sexspec-avg-min.plink.map", package="skater", mustWork=TRUE)
gmap <- read_map(gmapfile)

}
\references{
\url{http://zzz.bwh.harvard.edu/plink/data.shtml#map}

\url{http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/}

\url{https://github.com/williamslab/ped-sim#map-file}

\url{https://www.nature.com/articles/ncomms14994}

\url{https://www.nature.com/articles/ncomms14994}

\url{https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confusion_matrix.R
\name{calc_stats}
\alias{calc_stats}
\title{Calculate various statistics from a confusion matrix}
\usage{
calc_stats(tabble, prevalence = NULL, positive, ...)
}
\arguments{
\item{tabble}{A frequency table created with \code{\link{table}}}

\item{prevalence}{Prevalence value. Default is \code{NULL}}

\item{positive}{Positive class}

\item{...}{Other, not currently used}
}
\value{
A tibble with (at present) columns for sensitivity, specificity, PPV, NPV, F1 score, detection rate, detection prevalence, balanced accuracy, FDR, FOR, FPR, FNR.  For
more than 2 classes, these statistics are provided for each class.
}
\description{
Given a frequency table of predictions versus target values,
calculate numerous statistics of interest.
}
\details{
Used within confusion_matrix to calculate various confusion matrix
metrics. This is called by \code{confusion_matrix}, but if this is all you
want you can simply supply the table.

Suppose a 2x2 table with notation

\tabular{rcc}{ \tab target \tab \cr Predicted \tab Event \tab No Event
\cr Event \tab A \tab B \cr No Event \tab C \tab D \cr }

The formulas used here are:
\deqn{Sensitivity = A/(A+C)}
\deqn{Specificity = D/(B+D)}
\deqn{Prevalence = (A+C)/(A+B+C+D)}
\deqn{Positive Predictive Value = (sensitivity * prevalence)/((sensitivity*prevalence) + ((1-specificity)*(1-prevalence)))}
\deqn{Negative Predictive Value = (specificity * (1-prevalence))/(((1-sensitivity)*prevalence) + ((specificity)*(1-prevalence)))} \deqn{Detection Rate = A/(A+B+C+D)}
\deqn{Detection Prevalence = (A+B)/(A+B+C+D)}
\deqn{Balanced Accuracy = (sensitivity+specificity)/2}
\deqn{Precision = A/(A+B)}
\deqn{Recall = A/(A+C)}
\deqn{F1 = harmonic mean of precision and recall = (1+beta^2)*precision*recall/((beta^2 * precision)+recall)}
where \code{beta = 1} for this function.
\deqn{False Discovery Rate = 1 - Positive Predictive Value}
\deqn{False Omission Rate = 1 - Negative Predictive Value}
\deqn{False Positive Rate = 1 - Specificity}
\deqn{False Negative Rate = 1 - Sensitivity}
\deqn{D' = qnorm(Sensitivity) - qnorm(1 - Specificity)}
\deqn{AUC ~= pnorm(D'/sqrt(2))}

See the references for discussions of the first five formulas.
Abbreviations:
\describe{
\item{Positive Predictive Value: PPV}{}
\item{Negative Predictive Value: NPV}{}
\item{False Discovery Rate: FDR}{}
\item{False Omission Rate: FOR}{}
\item{False Positive Rate: FPR}{}
\item{False Negative Rate: FNR}{}
}
}
\note{
Different names are used for the same statistics.
\describe{
\item{Sensitivity: True Positive Rate, Recall, Hit Rate, Power}{}
\item{Specificity: True Negative Rate}{}
\item{Positive Predictive Value: Precision}{}
\item{False Negative Rate: Miss Rate, Type II error rate, beta}{}
\item{False Positive Rate: Fallout, Type I error rate, alpha}{}
}

This function is called by \code{confusion_matrix}, but if this is all you
want, you can simply supply the table to this function.
}
\references{
Kuhn, M. (2008), "Building predictive models in R using the
caret package, " \emph{Journal of Statistical Software},
(\url{https://www.jstatsoft.org/article/view/v028i05}).

Altman, D.G., Bland, J.M. (1994) "Diagnostic tests 1: sensitivity and
specificity", \emph{British Medical Journal}, vol 308, 1552.

Altman, D.G., Bland, J.M. (1994) "Diagnostic tests 2: predictive values,"
\emph{British Medical Journal}, vol 309, 102.

Velez, D.R., et. al. (2008) "A balanced accuracy function for epistasis
modeling in imbalanced datasets using multifactor dimensionality
reduction.," \emph{Genetic Epidemiology}, vol 4, 306.
}
\author{
Michael Clark (see \href{https://github.com/m-clark/confusionMatrix}{m-clark/confusion_matrix}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_ibd}
\alias{read_ibd}
\title{Read IBD segment file}
\usage{
read_ibd(file, source)
}
\arguments{
\item{file}{Input file path}

\item{source}{Source of the input file; must be one of \code{"hapibd"} or \code{"pedsim"}}
}
\value{
if \code{source="hapibd"}, a tibble is returned.
If \code{source="pedsim"}, a list with two tibble elements, \code{IBD1} and \code{IBD2} is returned.
Both the \code{hapibd} tibble, and the two \code{pedsim} tibbles contain six columns:
\enumerate{
\item id1 (sample identifier 1)
\item id2 (sample identifier 2)
\item chr (chromosome)
\item start (segment bp start coordinate)
\item end (segment bp end coordinate)
\item length (shared segment length in genetic units, cM)
}
}
\description{
Reads in the inferred IBD segments from \code{hapibd} (\href{https://github.com/browning-lab/hap-ibd#output-files}{documentation}) or IBD segment file generated by \code{ped-sim} (\href{https://github.com/williamslab/ped-sim#output-ibd-segments-file}{documentation}).

If reading a \code{hapibd} segment file, the input data should have the following columns:
\enumerate{
\item First sample identifier
\item First sample haplotype index (1 or 2)
\item Second sample identifier
\item Second sample haplotype index (1 or 2)
\item Chromosome
\item Base coordinate of first marker in segment
\item Base coordinate of last marker in segment
\item cM length of IBD segment
}

If read a \code{pedsim} segment file, the input data should have the following columns:
\enumerate{
\item First sample identifier
\item Second sample identifer
\item Chromosome
\item Physical position start
\item Physical position end
\item IBD type
\item Genetic position start
\item Genetic position end
\item Genetic length (end - start)
}
}
\examples{
hapibd_fp <- system.file("extdata", "GBR.sim.ibd.gz", package="skater", mustWork=TRUE)
hapibd_seg <- read_ibd(hapibd_fp, source = "hapibd")
pedsim_fp <- system.file("extdata", "GBR.sim.seg.gz", package="skater", mustWork=TRUE)
pedsim_seg <- read_ibd(pedsim_fp, source = "pedsim")
}
\references{
\url{https://github.com/browning-lab/hap-ibd#output-files}

\url{https://github.com/williamslab/ped-sim#output-ibd-segments-file}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{kin2degree}
\alias{kin2degree}
\title{Kinship coefficient to degree}
\usage{
kin2degree(k, max_degree = 3L)
}
\arguments{
\item{k}{Kinship coefficient (numeric, typically between 0 and .5, although KING can produce values <0).}

\item{max_degree}{Max degree resolution (default 3). Used to seed
\link[=dibble]{dibble}. Anything below the inference range of
\code{max_degree} will report \code{NA}. See \link[=dibble]{dibble}.}
}
\value{
A vector with inferred degree, up to the maximum degree in \code{dibble} (anything more distant is \code{NA}, i.e., unrelated).
}
\description{
Infers relationship degree given a kinship coefficient.
}
\examples{
kin2degree(0.5)
kin2degree(0.25)
kin2degree(0.125)
kin2degree(0.0625)
kin2degree(0.03125)
kin2degree(0.03125, max_degree=5)
kin2degree(-0.05)
k <- seq(.02, .5, .03)
kin2degree(k)
kin2degree(k, max_degree=5)
tibble::tibble(k=k) \%>\% dplyr::mutate(degree=kin2degree(k))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_akt}
\alias{read_akt}
\title{Read AKT kin output file}
\usage{
read_akt(file)
}
\arguments{
\item{file}{Input file path}
}
\value{
A tibble containing the 7 columns from the akt file.
}
\description{
Reads in an \verb{akt kin} \href{https://illumina.github.io/akt/#kin}{results file}. Input \code{file} must have seven columns, whitespace delimited:
\enumerate{
\item id1 (member 1)
\item id2 (member 2)
\item IBD0 (ratio of IBD0/All SNPS)
\item IBD1 (ratio of IBD1/All SNPS)
\item Kinship Coefficient
\item NSNPS
}
}
\examples{
aktFile <- system.file("extdata", "3gens.akt", package="skater", mustWork=TRUE)
akt <- read_akt(aktFile)
akt

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{kin2cm}
\alias{kin2cm}
\title{Kinship coefficient to cM}
\usage{
kin2cm(k)
}
\arguments{
\item{k}{Kinship coefficient (numeric, typically between 0 and .5, although KING can produce values <0).}
}
\value{
A vector of numeric estimated cM, ranging from 0-3560.
}
\description{
"Converts" a kinship coefficient to put on the same scale as shared cM using the formula
\eqn{cm <- pmin(3560, 4*pmax(0, k)*3560)}.
}
\examples{
kin2cm(.25)
kin2cm(.125)
kin2cm(.0625)
dibble(9) \%>\% dplyr::mutate(cm=kin2cm(k))

}
\references{
\url{https://dnapainter.com/tools/sharedcmv4}.

\url{https://www.ancestry.com/dna/resource/whitePaper/AncestryDNA-Matching-White-Paper.pdf}.

\url{https://verogen.com/wp-content/uploads/2021/03/snp-typing-uas-kinship-estimation-gedmatch-pro-tech-note-vd2020058-a.pdf}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ped.R
\name{fam2ped}
\alias{fam2ped}
\title{Fam to pedigree}
\usage{
fam2ped(fam)
}
\arguments{
\item{fam}{A tibble with six columns of PLINK .fam data as read in by \link{read_fam}.}
}
\value{
A tibble with new listcol \code{ped} containing pedigrees from \code{kinship2::pedigree}.
}
\description{
Converts a \href{https://www.cog-genomics.org/plink/1.9/formats#fam}{PLINK-formatted fam file} to a pedigree object using \link[kinship2:pedigree]{kinship2::pedigree}.
}
\examples{
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
fam <- read_fam(famfile)
fam2ped(fam)

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
See \verb{magrittr::[\\\%>\\\%][magrittr::pipe]} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ped.R
\name{ped2kinpair}
\alias{ped2kinpair}
\title{Pedigree to pairwise kinship}
\usage{
ped2kinpair(ped)
}
\arguments{
\item{ped}{A "pedigree" class object from \link[=fam2ped]{fam2ped}.}
}
\value{
A tibble containing all pairwise kinship coefficients from the input pedigree.
}
\description{
Converts a pedigree class object from \link[=fam2ped]{fam2ped} to a pairwise list of relationships and their expected/theoretical kinship coefficient.
}
\examples{
famfile <- system.file("extdata", "3gens.fam", package="skater", mustWork=TRUE)
famfile \%>\%
  read_fam() \%>\%
  fam2ped() \%>\%
  dplyr::mutate(kinpairs=purrr::map(ped, ped2kinpair)) \%>\%
  dplyr::select(fid, kinpairs) \%>\%
  tidyr::unnest(cols=kinpairs)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confusion_matrix.R
\name{calc_accuracy}
\alias{calc_accuracy}
\title{Calculate Accuracy}
\usage{
calc_accuracy(tabble)
}
\arguments{
\item{tabble}{A frequency table created with \code{\link{table}}}
}
\value{
A tibble with the corresponding statistics
}
\description{
Calculates accuracy and related metrics.
}
\details{
Calculates accuracy, lower and upper bounds, the guessing rate and
p-value of the accuracy vs. the guessing rate. This function is called by
\code{confusion_matrix}, but if this is all you want, you can simply supply
the table to this function.
}
\seealso{
\code{\link{binom.test}}
}
\author{
Michael Clark (see \href{https://github.com/m-clark/confusionMatrix}{m-clark/confusion_matrix}).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{read_plink2_king}
\alias{read_plink2_king}
\title{Read PLINK KING table}
\usage{
read_plink2_king(file)
}
\arguments{
\item{file}{Input file path}
}
\value{
A tibble containing the 6 columns from the \code{plink2 --make-king-table} output.
}
\description{
Reads in the output from \code{plink2 --make-king-table} (\href{https://www.cog-genomics.org/plink/2.0/distance#make_king}{documentation}).
Input \code{file} must have six columns, tab delimited:
\enumerate{
\item id1 (member 1)
\item id2 (member 2)
\item nsnps
\item hethet: proportion of sites where both are heterozygous
\item k: Kinship Coefficient
}
}
\examples{
plink2kingFile <- system.file("extdata", "plink2-king-table.tsv", package="skater", mustWork=TRUE)
plink2king <- read_plink2_king(plink2kingFile)
plink2king
plink2king \%>\% dplyr::filter(k>0.01)

}
\references{
\url{https://www.cog-genomics.org/plink/2.0/distance#make_king}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/skater-package.R
\docType{package}
\name{skater-package}
\alias{skater}
\alias{skater-package}
\title{skater: Utilities for SNP-Based Kinship Analysis}
\description{
Utilities for single nucleotide polymorphism (SNP) based kinship analysis testing and evaluation. The 'skater' package contains functions for importing, parsing, and analyzing pedigree data, performing relationship degree inference, benchmarking relationship degree classification, and summarizing identity by descent (IBD) segment data. Package functions and methods are described in Turner et al. (2021) "skater: An R package for SNP-based Kinship Analysis, Testing, and Evaluation" <doi:10.1101/2021.07.21.453083>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/signaturescience/skater}
  \item Report bugs at \url{https://github.com/signaturescience/skater/issues}
}

}
\author{
\strong{Maintainer}: Stephen Turner \email{vustephen@gmail.com} (\href{https://orcid.org/0000-0001-9140-9028}{ORCID})

Authors:
\itemize{
  \item Matthew Scholz (\href{https://orcid.org/0000-0003-3686-1227}{ORCID})
  \item VP Nagraj (\href{https://orcid.org/0000-0003-0060-566X}{ORCID})
}

Other contributors:
\itemize{
  \item Signature Science, LLC. [copyright holder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dibble}
\alias{dibble}
\title{Degree tibble}
\usage{
dibble(max_degree = 3L)
}
\arguments{
\item{max_degree}{The most distant degree you want to measure (usually between 3-9, default 3).}
}
\value{
A tibble containing the degree, expected kinship coefficient (\code{k}),
lower (\code{l}) and upper (\code{u}) inference bounds.
}
\description{
Creates a tibble with degree, expected kinship coefficient, and inference boundaries.

Rows will be created up to the \code{max_degree}, with an additional row for any
relationship more distant than \code{max_degree}. The \code{degree} value for the final
row will be \code{NA}. This represents inference criteria for "unrelated"
individuals. See examples.
}
\examples{
dibble(3)
dibble(10)

}
