
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<table class="table">
<thead>
<tr class="header">
<th align="left">
Usage
</th>
<th align="left">
Release
</th>
<th align="left">
Development
</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">
<a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://img.shields.io/badge/licence-GPL--3-blue.svg" alt="Licence"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/package=UCSCXenaTools"><img src="https://www.r-pkg.org/badges/version/UCSCXenaTools" alt="CRAN"></a>
</td>
<td align="left">
<a href="https://travis-ci.org/ropensci/UCSCXenaTools"><img src="https://travis-ci.org/ropensci/UCSCXenaTools.svg?branch=master" alt="Travis"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://cran.r-project.org/"><img src="https://img.shields.io/badge/R%3E%3D-3.5.0-blue.svg" alt="minimal R version"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/web/checks/check_results_UCSCXenaTools.html"><img src="https://cranchecks.info/badges/summary/UCSCXenaTools" alt="cran-checks"></a>
</td>
<td align="left">
<a href="https://ci.appveyor.com/project/ShixiangWang/UCSCXenaTools"><img src="https://ci.appveyor.com/api/projects/status/github/ropensci/UCSCXenaTools?branch=master&svg=true" alt="AppVeyor"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://tinyverse.netlify.com/badge/UCSCXenaTools"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/software-review/issues/315"><img src="https://badges.ropensci.org/315_status.svg" alt="rOpenSci"></a>
</td>
<td align="left">
<a href="https://codecov.io/github/ShixiangWang/UCSCXenaTools?branch=master"><img src="https://codecov.io/github/ShixiangWang/UCSCXenaTools/coverage.svg?branch=master" alt="Codecov"></a>
</td>
</tr>
<tr class="even">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://cranlogs.r-pkg.org/badges/grand-total/UCSCXenaTools" alt="downloads"></a>
</td>
<td align="left">
<a href="https://zenodo.org/badge/latestdoi/178662770"><img src="https://zenodo.org/badge/178662770.svg" alt="DOI"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/UCSCXenaTools/issues?q=is%3Aissue+is%3Aclosed"><img src="https://img.shields.io/github/issues-closed/ropensci/UCSCXenaTools.svg" alt="Closed issues"></a>
</td>
</tr>
<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://cranlogs.r-pkg.org/badges/UCSCXenaTools" alt="month-downloads"></a>
</td>
<td align="left">
<a href="https://doi.org/10.21105/joss.01627"><img src="https://joss.theoj.org/papers/10.21105/joss.01627/status.svg" alt="JOSS" >
</a>
</td>
<td align="left">
<a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Project Status: Active – The project has reached a stable, usable state and is being actively developed." /></a>
</td>
</tr>
</tbody>
</table>

<br> <!-- badges: end -->

# UCSCXenaTools <img src='man/figures/logo.png' align="right" height="200" alt="logo"/>

**UCSCXenaTools** is an R package for accessing genomics data from UCSC
Xena platform, from cancer multi-omics to single-cell RNA-seq. Public
omics data from UCSC Xena are supported through [**multiple turn-key
Xena Hubs**](https://xenabrowser.net/datapages/), which are a collection
of UCSC-hosted public databases such as TCGA, ICGC, TARGET, GTEx, CCLE,
and others. Databases are normalized so they can be combined, linked,
filtered, explored and downloaded.

**Who is the target audience and what are scientific applications of
this package?**

-   Target Audience: cancer and clinical researchers, bioinformaticians
-   Applications: genomic and clinical analyses

## Table of Contents

-   [Installation](#installation)
-   [Data Hub List](#data-hub-list)
-   [Basic usage](#basic-usage)
-   [Citation](#citation)
-   [How to contribute](#how-to-contribute)
-   [Acknowledgment](#acknowledgment)

## Installation

Install stable release from CRAN with:

``` r
install.packages("UCSCXenaTools")
```

You can also install devel version of **UCSCXenaTools** from github
with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/UCSCXenaTools")
```

If you want to build vignette in local, please add two options:

``` r
remotes::install_github("ropensci/UCSCXenaTools", build_vignettes = TRUE, dependencies = TRUE)
```

## Data Hub List

All datasets are available at <https://xenabrowser.net/datapages/>.

Currently, **UCSCXenaTools** supports the following data hubs of UCSC
Xena.

-   UCSC Public Hub: <https://ucscpublic.xenahubs.net/>
-   TCGA Hub: <https://tcga.xenahubs.net/>
-   GDC Xena Hub: <https://gdc.xenahubs.net/>
-   ICGC Xena Hub: <https://icgc.xenahubs.net/>
-   Pan-Cancer Atlas Hub: <https://pancanatlas.xenahubs.net/>
-   UCSC Toil RNAseq Recompute Compendium Hub:
    <https://toil.xenahubs.net/>
-   PCAWG Xena Hub: <https://pcawg.xenahubs.net/>
-   ATAC-seq Hub: <https://atacseq.xenahubs.net/>
-   Singel Cell Xena Hub: <https://singlecellnew.xenahubs.net/>
-   Kids First Xena Hub: <https://kidsfirst.xenahubs.net/>
-   Treehouse Xena Hub: <https://xena.treehouse.gi.ucsc.edu:443/>

Users can update dataset list from the newest version of UCSC Xena by
hand with `XenaDataUpdate()` function, followed by restarting R and
`library(UCSCXenaTools)`.

If any url of data hub is changed or a new data hub is online, please
remind me by emailing to <w_shixiang@163.com> or [opening an issue on
GitHub](https://github.com/ropensci/UCSCXenaTools/issues).

## Basic usage

Download UCSC Xena datasets and load them into R by **UCSCXenaTools** is
a workflow with `generate`, `filter`, `query`, `download` and `prepare`
5 steps, which are implemented as `XenaGenerate`, `XenaFilter`,
`XenaQuery`, `XenaDownload` and `XenaPrepare` functions, respectively.
They are very clear and easy to use and combine with other packages like
`dplyr`.

To show the basic usage of **UCSCXenaTools**, we will download clinical
data of LUNG, LUAD, LUSC from TCGA (hg19 version) data hub. Users can
learn more about **UCSCXenaTools** by running
`browseVignettes("UCSCXenaTools")` to read vignette.

### XenaData data.frame

**UCSCXenaTools** uses a `data.frame` object (built in package)
`XenaData` to generate an instance of `XenaHub` class, which records
information of all datasets of UCSC Xena Data Hubs.

You can load `XenaData` after loading `UCSCXenaTools` into R.

``` r
library(UCSCXenaTools)
#> =========================================================================================
#> UCSCXenaTools version 1.4.6
#> Project URL: https://github.com/ropensci/UCSCXenaTools
#> Usages: https://cran.r-project.org/web/packages/UCSCXenaTools/vignettes/USCSXenaTools.html
#> 
#> If you use it in published research, please cite:
#> Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
#>   from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq.
#>   Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627
#> =========================================================================================
#>                               --Enjoy it--
data(XenaData)

head(XenaData)
#> # A tibble: 6 x 17
#>   XenaHosts XenaHostNames XenaCohorts XenaDatasets SampleCount DataSubtype Label
#>   <chr>     <chr>         <chr>       <chr>              <int> <chr>       <chr>
#> 1 https://… publicHub     Breast Can… ucsfNeve_pu…          51 gene expre… Neve…
#> 2 https://… publicHub     Breast Can… ucsfNeve_pu…          57 phenotype   Phen…
#> 3 https://… publicHub     Glioma (Ko… kotliarov20…         194 copy number Kotl…
#> 4 https://… publicHub     Glioma (Ko… kotliarov20…         194 phenotype   Phen…
#> 5 https://… publicHub     Lung Cance… weir2007_pu…         383 copy number CGH  
#> 6 https://… publicHub     Lung Cance… weir2007_pu…         383 phenotype   Phen…
#> # … with 10 more variables: Type <chr>, AnatomicalOrigin <chr>,
#> #   SampleType <chr>, Tags <chr>, ProbeMap <chr>, LongTitle <chr>,
#> #   Citation <chr>, Version <chr>, Unit <chr>, Platform <chr>
```

### Workflow

Select datasets.

``` r
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "clinical") %>% 
  XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo

df_todo
#> class: XenaHub 
#> hosts():
#>   https://tcga.xenahubs.net
#> cohorts() (3 total):
#>   TCGA Lung Cancer (LUNG)
#>   TCGA Lung Adenocarcinoma (LUAD)
#>   TCGA Lung Squamous Cell Carcinoma (LUSC)
#> datasets() (3 total):
#>   TCGA.LUNG.sampleMap/LUNG_clinicalMatrix
#>   TCGA.LUAD.sampleMap/LUAD_clinicalMatrix
#>   TCGA.LUSC.sampleMap/LUSC_clinicalMatrix
```

Query and download.

``` r
XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
```

**For researchers in China, now Hiplot team has deployed several Xena
mirror sites (`https://xena.hiplot.com.cn/`) at Shanghai. You can set an
option `options(use_hiplot = TRUE)` before querying data step to speed
up both data querying and downloading.**

``` r
options(use_hiplot = TRUE)

XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
#> The hiplot server may down, we will not use it for now.
#> This will check url status, please be patient.
#> All downloaded files will under directory /var/folders/bj/nw1w4g1j37ddpgb6zmh3sfh80000gn/T//RtmpoHs0nL.
#> The 'trans_slash' option is FALSE, keep same directory structure as Xena.
#> Creating directories for datasets...
#> Downloading TCGA.LUNG.sampleMap/LUNG_clinicalMatrix
#> Downloading TCGA.LUAD.sampleMap/LUAD_clinicalMatrix
#> Downloading TCGA.LUSC.sampleMap/LUSC_clinicalMatrix
```

Prepare data into R for analysis.

``` r
cli = XenaPrepare(xe_download)
class(cli)
#> [1] "list"
names(cli)
#> [1] "LUNG_clinicalMatrix" "LUAD_clinicalMatrix" "LUSC_clinicalMatrix"
```

## More to read

-   [Introduction and basic usage of
    UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-intro/)
-   [UCSCXenaTools: Retrieve Gene Expression and Clinical Information
    from UCSC Xena for Survival
    Analysis](https://shixiangwang.github.io/home/en/post/ucscxenatools-201908/)
-   [Obtain RNAseq Values for a Specific Gene in Xena
    Database](https://shixiangwang.github.io/home/en/post/2020-07-22-ucscxenatools-single-gene/)
-   [UCSC Xena Access APIs in
    UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-api/)

## Citation

Cite me by the following paper.

    Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
      from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq. 
      Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627

    # For BibTex
      
    @article{Wang2019UCSCXenaTools,
        journal = {Journal of Open Source Software},
        doi = {10.21105/joss.01627},
        issn = {2475-9066},
        number = {40},
        publisher = {The Open Journal},
        title = {The UCSCXenaTools R package: a toolkit for accessing genomics data from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq},
        url = {https://dx.doi.org/10.21105/joss.01627},
        volume = {4},
        author = {Wang, Shixiang and Liu, Xuesong},
        pages = {1627},
        date = {2019-08-05},
        year = {2019},
        month = {8},
        day = {5},
    }

Cite UCSC Xena by the following paper.

    Goldman, Mary, et al. "The UCSC Xena Platform for cancer genomics data 
        visualization and interpretation." BioRxiv (2019): 326470.

## How to contribute

For anyone who wants to contribute, please follow the guideline:

-   Clone project from GitHub
-   Open `UCSCXenaTools.Rproj` with RStudio
-   Modify source code
-   Run `devtools::check()`, and fix all errors, warnings and notes
-   Create a pull request

## Acknowledgment

This package is based on [XenaR](https://github.com/mtmorgan/XenaR),
thanks [Martin Morgan](https://github.com/mtmorgan) for his work.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# UCSCXenaTools 1.4.7

- 3 datasets were added into `XenaData`.

```r
[1] "TCGA_survival_data_2.txt"                                  
[2] "clinical_CellLinePolyA_21.06_2021-06-15.tsv"               
[3] "CellLinePolyA_21.06_hugo_log2tpm_58581genes_2021-06-15.tsv"
```

- Added support for gene symbol checking and data cache in `fetch()`.

# UCSCXenaTools 1.4.6

- Added code to check hiplot server status.
- Fixed check warnings to follow CRAN policy (#36).

# UCSCXenaTools 1.4.5

- Fixed the download bug for pan-cancer data hub due to unvalid URL. `url_encode()`
is added internally to transform reserved characters (`/` is kept).
- Fixed the download bug because UCSCXenaShiny mutate the result of `XenaQuery()`
and thus change the column number. This bug will not affect XenaTools itself.

# UCSCXenaTools 1.4.4

- Fixed the download bug because UCSCXenaShiny mutate the result of `XenaQuery()`
and thus change the column order. This bug will not affect XenaTools itself.

# UCSCXenaTools 1.4.3

- Implemented using hiplot for `fetch()` functions.

# UCSCXenaTools 1.4.1 (1.4.2)

- Removed unpublished Xena hub TDI.

# UCSCXenaTools 1.4.0

- Supported downloading data from Hiplot mirror site (`https://xena.hiplot.com.cn/`).

# UCSCXenaTools 1.3.6

- Fixed a bug about try times for data download. 
- Make sure a message instead of an error will be returned if download process failed.

# UCSCXenaTools 1.3.5

- Added TDI data Hub containing 9 new datasets.

# UCSCXenaTools 1.3.4

* Updated UCSC Xena datasets, now 1670 datasets available.

# UCSCXenaTools 1.3.3

* Added `fetch_sparse_values()` function.
* Updated treehouse URL.
* Added treehouse datasets.

# UCSCXenaTools 1.3.2

* Fixed bug about an error happened when querying mutation data.
* Dropped "Treehouse" data hub.
* Updated code to update Xena hub datasets.

# UCSCXenaTools 1.3.1

* Added `max_try` option in query and download functions, so they can handle internet connection issue better

# UCSCXenaTools 1.3.0

* Added a new data hub: PCAWG Xena Hub (#24). 
* Added a new data hub: Kids First Xena Hub (#24).
* Updated data update function more robustly.
---
title: 'The UCSCXenaTools R package: a toolkit for accessing genomics data from UCSC Xena platform,
  from cancer multi-omics to single-cell RNA-seq'
authors:
- affiliation: '1, 2, 3'
  name: Shixiang Wang
  orcid: 0000-0001-9855-7357
- affiliation: 1
  name: Xuesong Liu
  orcid: 0000-0002-7736-0077
date: "24 July 2019"
bibliography: paper.bib
tags:
- R
- cancer genomics
- data access
affiliations:
- index: 1
  name: School of Life Science and Technology, ShanghaiTech University
- index: 2
  name: Shanghai Institute of Biochemistry and Cell Biology, Chinese Academy of Sciences
- index: 3
  name: University of Chinese Academy of Sciences
---

# Summary

UCSC Xena platform (https://xenabrowser.net/) provides unprecedented resource for public omics data [@goldman2019ucsc]
from big projects like The Cancer Genome Atlas (TCGA) [@weinstein2013cancer], 
International Cancer Genome Consortium Data Portal (ICGC) [@zhang2011international],
The Cancer Cell Line Encyclopedia (CCLE) [@barretina2012cancer], or reserach groups like @mullighan2008genomic, @puram2017single.
All available data types include single-nucleotide variants (SNVs), small insertions and deletions (INDELs), large structural variants, copy number variation (CNV), expression, DNA methylation, ATAC-seq signals, and phenotypic annotations. 

Despite UCSC Xena platform itself allows users to explore and analyze data, it is hard
for users to incorporate multiple datasets or data types, integrate the selected data with 
popular analysis tools or homebrewed code, and reproduce analysis procedures.
R language is well established and extensively used standard in statistical and bioinformatics research.
Here, we introduce an R package UCSCXenaTools for enabling data retrieval, analysis integration and 
reproducible research for omics data from UCSC Xena platform.

Currently, UCSCXenaTools supports downloading over 1600 datasets from 10 data hubs of UCSC Xena platform
as shown in the following table. Typically, downloading UCSC Xena datasets and loading them into R by UCSCXenaTools 
is a workflow with generate, filter, query, download and prepare 5 steps, which are implemented as functions.
They are very clear and easy to use and combine with other packages like dplyr [@wickham2015dplyr].
Besides, UCSCXenaTools can also query and download subset of a target dataset, 
this is particularly useful when
user focus on studying one object like gene or protein. The key features are summarized in Figure 1.


|Data hub       | Dataset count|URL                                |
|:--------------|-------------:|:----------------------------------|
|tcgaHub        |           879|https://tcga.xenahubs.net          |
|gdcHub         |           449|https://gdc.xenahubs.net           |
|publicHub      |           104|https://ucscpublic.xenahubs.net    |
|pcawgHub       |            53|https://pcawg.xenahubs.net         |
|toilHub        |            50|https://toil.xenahubs.net          |
|singlecellHub  |            45|https://singlecell.xenahubs.net    |
|icgcHub        |            23|https://icgc.xenahubs.net          |
|pancanAtlasHub |            19|https://pancanatlas.xenahubs.net   |
|treehouseHub   |            15|https://xena.treehouse.gi.ucsc.edu |
|atacseqHub     |             9|https://atacseq.xenahubs.net       |

![Overview of UCSCXenaTools](overview.png)

# Acknowledgements

We thank Christine Stawitz and Carl Ganz for their constructive comments.
This package is based on R package [XenaR](https://github.com/mtmorgan/XenaR), thanks [Martin Morgan](https://github.com/mtmorgan) for his work.

# References
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

<!-- badges: start -->

<table class="table">
<thead>
<tr class="header">
<th align="left">
Usage
</th>
<th align="left">
Release
</th>
<th align="left">
Development
</th>
</tr>
</thead>
<tbody>

<tr class="odd">
<td align="left">
<a href="https://www.gnu.org/licenses/gpl-3.0.en.html"><img src="https://img.shields.io/badge/licence-GPL--3-blue.svg" alt="Licence"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/package=UCSCXenaTools"><img src="https://www.r-pkg.org/badges/version/UCSCXenaTools" alt="CRAN"></a>
</td>
<td align="left">
<a href="https://travis-ci.org/ropensci/UCSCXenaTools"><img src="https://travis-ci.org/ropensci/UCSCXenaTools.svg?branch=master" alt="Travis"></a>
</td>
</tr>

<tr class="even">
<td align="left">
<a href="https://cran.r-project.org/"><img src="https://img.shields.io/badge/R%3E%3D-3.5.0-blue.svg" alt="minimal R version"></a>
</td>
<td align="left">
<a href="https://cran.r-project.org/web/checks/check_results_UCSCXenaTools.html"><img src="https://cranchecks.info/badges/summary/UCSCXenaTools" alt="cran-checks"></a>
</td>
<td align="left">
<a href="https://ci.appveyor.com/project/ShixiangWang/UCSCXenaTools"><img src="https://ci.appveyor.com/api/projects/status/github/ropensci/UCSCXenaTools?branch=master&svg=true" alt="AppVeyor"></a>
</td>
</tr>

<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://tinyverse.netlify.com/badge/UCSCXenaTools"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/software-review/issues/315"><img src="https://badges.ropensci.org/315_status.svg" alt="rOpenSci"></a>
</td>
<td align="left">
<a href="https://codecov.io/github/ShixiangWang/UCSCXenaTools?branch=master"><img src="https://codecov.io/github/ShixiangWang/UCSCXenaTools/coverage.svg?branch=master" alt="Codecov"></a>
</td>
</tr>

<tr class="even">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://cranlogs.r-pkg.org/badges/grand-total/UCSCXenaTools" alt="downloads"></a>
</td>
<td align="left">
<a href="https://zenodo.org/badge/latestdoi/178662770"><img src="https://zenodo.org/badge/178662770.svg" alt="DOI"></a>
</td>
<td align="left">
<a href="https://github.com/ropensci/UCSCXenaTools/issues?q=is%3Aissue+is%3Aclosed"><img src="https://img.shields.io/github/issues-closed/ropensci/UCSCXenaTools.svg" alt="Closed issues"></a>
</td>

</tr>
<tr class="odd">
<td align="left">
<a href="https://CRAN.R-project.org/package=UCSCXenaTools"><img src="https://cranlogs.r-pkg.org/badges/UCSCXenaTools" alt="month-downloads"></a>
</td>
<td align="left">
<a href="https://doi.org/10.21105/joss.01627"><img src="https://joss.theoj.org/papers/10.21105/joss.01627/status.svg" alt="JOSS" >
</a>
</td>
<td align="left">
<a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Project Status: Active – The project has reached a stable, usable state and is being actively developed." /></a>
</td>
</tr>
</tbody>
</table>
<br>
<!-- badges: end -->

# UCSCXenaTools <img src='man/figures/logo.png' align="right" height="200" alt="logo"/>

**UCSCXenaTools** is an R package for accessing genomics data from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq. 
Public omics data from UCSC Xena are supported through [**multiple turn-key Xena Hubs**](https://xenabrowser.net/datapages/), which are a collection of UCSC-hosted public databases such as TCGA, ICGC, TARGET, GTEx, CCLE, and others. Databases are normalized so they can be combined, linked, filtered, explored and downloaded.

**Who is the target audience and what are scientific applications of this package?**

* Target Audience: cancer and clinical researchers, bioinformaticians
* Applications: genomic and clinical analyses

## Table of Contents

* [Installation](#installation)
* [Data Hub List](#data-hub-list)
* [Basic usage](#basic-usage)
* [Citation](#citation)
* [How to contribute](#how-to-contribute)
* [Acknowledgment](#acknowledgment)

## Installation

Install stable release from CRAN with:

```{r, eval=FALSE}
install.packages("UCSCXenaTools")
```

You can also install devel version of **UCSCXenaTools** from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/UCSCXenaTools")
```

If you want to build vignette in local, please add two options:

```{r, eval=FALSE}
remotes::install_github("ropensci/UCSCXenaTools", build_vignettes = TRUE, dependencies = TRUE)
```

## Data Hub List

All datasets are available at <https://xenabrowser.net/datapages/>.

Currently, **UCSCXenaTools** supports the following data hubs of UCSC Xena.

* UCSC Public Hub: <https://ucscpublic.xenahubs.net/>
* TCGA Hub: <https://tcga.xenahubs.net/>
* GDC Xena Hub: <https://gdc.xenahubs.net/>
* ICGC Xena Hub: <https://icgc.xenahubs.net/>
* Pan-Cancer Atlas Hub: <https://pancanatlas.xenahubs.net/>
* UCSC Toil RNAseq Recompute Compendium Hub: <https://toil.xenahubs.net/>
* PCAWG Xena Hub: <https://pcawg.xenahubs.net/>
* ATAC-seq Hub: <https://atacseq.xenahubs.net/>
* Singel Cell Xena Hub: <https://singlecellnew.xenahubs.net/>
* Kids First Xena Hub: <https://kidsfirst.xenahubs.net/>
* Treehouse Xena Hub: <https://xena.treehouse.gi.ucsc.edu:443/>

Users can update dataset list from the newest version of UCSC Xena by hand with `XenaDataUpdate()` function, followed
by restarting R and `library(UCSCXenaTools)`.

If any url of data hub is changed or a new data hub is online, please remind me by emailing to <w_shixiang@163.com> or [opening an issue on GitHub](https://github.com/ropensci/UCSCXenaTools/issues).


## Basic usage

Download UCSC Xena datasets and load them into R by **UCSCXenaTools** is a workflow with `generate`, `filter`, `query`, `download` and `prepare` 5 steps, which are implemented as `XenaGenerate`, `XenaFilter`, `XenaQuery`, `XenaDownload` and `XenaPrepare` functions, respectively. They are very clear and easy to use and combine with other packages like `dplyr`.

To show the basic usage of **UCSCXenaTools**, we will download clinical data of LUNG, LUAD, LUSC from TCGA (hg19 version) data hub. Users can learn more about **UCSCXenaTools** by running `browseVignettes("UCSCXenaTools")` to read vignette.

### XenaData data.frame

**UCSCXenaTools** uses a `data.frame` object (built in package) `XenaData` to generate an instance of `XenaHub` class, which records information of all datasets of UCSC Xena Data Hubs.

You can load `XenaData` after loading `UCSCXenaTools` into R.

```{r}
library(UCSCXenaTools)
data(XenaData)

head(XenaData)
```

### Workflow

Select datasets.

```{r}
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "clinical") %>% 
  XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo

df_todo
```

Query and download.

```{r, eval=FALSE}
XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
```

**For researchers in China, now Hiplot team has deployed several Xena mirror sites (`https://xena.hiplot.com.cn/`) at Shanghai. You can set an option `options(use_hiplot = TRUE)` before querying data step to speed up both data querying and downloading.**

```{r}
options(use_hiplot = TRUE)

XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
```

Prepare data into R for analysis.

```{r}
cli = XenaPrepare(xe_download)
class(cli)
names(cli)
```

## More to read

- [Introduction and basic usage of UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-intro/)
- [UCSCXenaTools: Retrieve Gene Expression and Clinical Information from UCSC Xena for Survival Analysis](https://shixiangwang.github.io/home/en/post/ucscxenatools-201908/)
- [Obtain RNAseq Values for a Specific Gene in Xena Database](https://shixiangwang.github.io/home/en/post/2020-07-22-ucscxenatools-single-gene/)
- [UCSC Xena Access APIs in UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-api/)

## Citation

Cite me by the following paper.

```
Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
  from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq. 
  Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627

# For BibTex
  
@article{Wang2019UCSCXenaTools,
	journal = {Journal of Open Source Software},
	doi = {10.21105/joss.01627},
	issn = {2475-9066},
	number = {40},
	publisher = {The Open Journal},
	title = {The UCSCXenaTools R package: a toolkit for accessing genomics data from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq},
	url = {https://dx.doi.org/10.21105/joss.01627},
	volume = {4},
	author = {Wang, Shixiang and Liu, Xuesong},
	pages = {1627},
	date = {2019-08-05},
	year = {2019},
	month = {8},
	day = {5},
}
```

Cite UCSC Xena by the following paper. 

```
Goldman, Mary, et al. "The UCSC Xena Platform for cancer genomics data 
    visualization and interpretation." BioRxiv (2019): 326470.
```

## How to contribute

For anyone who wants to contribute, please follow the guideline:

* Clone project from GitHub
* Open `UCSCXenaTools.Rproj` with RStudio
* Modify source code 
* Run `devtools::check()`, and fix all errors, warnings and notes
* Create a pull request

## Acknowledgment

This package is based on [XenaR](https://github.com/mtmorgan/XenaR), thanks [Martin Morgan](https://github.com/mtmorgan) for his work.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
```{r}
library(UCSCXenaTools)
```


```{r}
df =  XenaData %>%
    dplyr::group_by(XenaHostNames) %>%
    dplyr::summarise(count = dplyr::n())
df
```

```{r}
urls = dplyr::tibble(
    name = as.character(.xena_hosts),
    url = names(.xena_hosts)
)
urls
```

```{r}
data = dplyr::left_join(
    df, urls, by = c("XenaHostNames"="name")
) %>% dplyr::arrange(dplyr::desc(count))


colnames(data) = c("Data hub", "Dataset count", "URL")
data
```


```{r}
knitr::kable(data, format = "markdown")
```


---
title: "UCSCXenaTools: an R package for Accessing Genomics Data from UCSC Xena platform, from Cancer Multi-omics to Single-cell RNA-seq"
author: "Shixiang Wang \\

        ShanghaiTech University"
date: "`r Sys.Date()`"

output:
  prettydoc::html_pretty:
    toc: true
    theme: cayman
    highlight: github
  pdf_document:
    toc: true
vignette: >
  %\VignetteIndexEntry{Basic usage}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


**UCSCXenaTools** is an R package for accessing genomics data from UCSC Xena platform, 
from cancer multi-omics to single-cell RNA-seq. 
Public omics data from UCSC Xena are supported through [**multiple turn-key Xena Hubs**](https://xenabrowser.net/datapages/), which are a collection of UCSC-hosted public databases such as TCGA, ICGC, TARGET, GTEx, CCLE, and others. Databases are normalized so they can be combined, linked, filtered, explored and downloaded.

**Who is the target audience and what are scientific applications of this package?**

* Target Audience: cancer and clinical researchers, bioinformaticians
* Applications: genomic and clinical analyses

## Installation

Install stable release from CRAN with:

```{r, eval=FALSE}
install.packages("UCSCXenaTools")
```

You can also install devel version of **UCSCXenaTools** from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/UCSCXenaTools")
```

If you want to build vignette in local, please add two options:

```{r, eval=FALSE}
remotes::install_github("ropensci/UCSCXenaTools", build_vignettes = TRUE, dependencies = TRUE)
```

The minimum versions to run the vignette is `1.2.4`. 
[GitHub Issue](https://github.com/ropensci/UCSCXenaTools/issues) is a place for discussing any problem.

## Data Hub List

All datasets are available at <https://xenabrowser.net/datapages/>.

Currently, **UCSCXenaTools** supports the following data hubs of UCSC Xena.

* UCSC Public Hub: <https://ucscpublic.xenahubs.net/>
* TCGA Hub: <https://tcga.xenahubs.net/>
* GDC Xena Hub: <https://gdc.xenahubs.net/>
* ICGC Xena Hub: <https://icgc.xenahubs.net/>
* Pan-Cancer Atlas Hub: <https://pancanatlas.xenahubs.net/>
* UCSC Toil RNAseq Recompute Compendium Hub: <https://toil.xenahubs.net/>
* PCAWG Xena Hub: <https://pcawg.xenahubs.net/>
* ATAC-seq Hub: <https://atacseq.xenahubs.net/>
* Singel Cell Xena Hub: <https://singlecellnew.xenahubs.net/>
* Kids First Xena Hub: <https://kidsfirst.xenahubs.net/>
* Treehouse Xena Hub: <https://xena.treehouse.gi.ucsc.edu:443/>

Users can update dataset list from the newest version of UCSC Xena by hand with `XenaDataUpdate()` function, followed
by restarting R and `library(UCSCXenaTools)`.

If any url of data hub is changed or a new data hub is online, please remind me by emailing to <w_shixiang@163.com> or [opening an issue on GitHub](https://github.com/ropensci/UCSCXenaTools/issues).


## Usage

Download UCSC Xena datasets and load them into R by **UCSCXenaTools** is a workflow with `generate`, `filter`, `query`, `download` and `prepare` 5 steps, which are implemented as `XenaGenerate`, `XenaFilter`, `XenaQuery`, `XenaDownload` and `XenaPrepare` functions, respectively. They are very clear and easy to use and combine with other packages like `dplyr`.

To show the basic usage of **UCSCXenaTools**, we will download clinical data of LUNG, LUAD, LUSC from TCGA (hg19 version) data hub.

### XenaData data.frame

**UCSCXenaTools** uses a `data.frame` object (built in package) `XenaData` to generate an instance of `XenaHub` class, which records information of all datasets of UCSC Xena Data Hubs.

You can load `XenaData` after loading `UCSCXenaTools` into R.

```{r}
library(UCSCXenaTools)
data(XenaData)

head(XenaData)
```

### Workflow

Select datasets.

```{r}
# The options in XenaFilter function support Regular Expression
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>% 
  XenaFilter(filterDatasets = "clinical") %>% 
  XenaFilter(filterDatasets = "LUAD|LUSC|LUNG") -> df_todo

df_todo
```

Sometimes we only know some keywords, `XenaScan()` can be used to scan all rows to detect if 
the keywords exist in `XenaData`.


```{r}
x1 = XenaScan(pattern = 'Blood')
x2 = XenaScan(pattern = 'LUNG', ignore.case = FALSE)

x1 %>%
    XenaGenerate()
x2 %>%
    XenaGenerate()
```

Query and download.

```{r, eval=FALSE}
XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
```

**For researchers in China, now Hiplot team has deployed several Xena mirror sites (`https://xena.hiplot.com.cn/`) at Shanghai. You can set an option `options(use_hiplot = TRUE)` before querying data step to speed up both data querying and downloading.**

```{r}
options(use_hiplot = TRUE)

XenaQuery(df_todo) %>%
  XenaDownload() -> xe_download
```

Prepare data into R for analysis.

```{r}
cli = XenaPrepare(xe_download)
class(cli)
names(cli)
```

### Browse datasets

Create two XenaHub objects:

* `to_browse` - a XenaHub object containing a cohort and a dataset.
* `to_browse2` - a XenaHub object containing 2 cohorts and 2 datasets.

```{r}
XenaGenerate(subset = XenaHostNames=="tcgaHub") %>%
    XenaFilter(filterDatasets = "clinical") %>%
    XenaFilter(filterDatasets = "LUAD") -> to_browse

to_browse

XenaGenerate(subset = XenaHostNames=="tcgaHub") %>%
    XenaFilter(filterDatasets = "clinical") %>%
    XenaFilter(filterDatasets = "LUAD|LUSC") -> to_browse2

to_browse2
```

`XenaBrowse()` function can be used to browse dataset/cohort links using your default web browser.
At default, this function limits one dataset/cohort for preventing user to open too many links at once. 

```{r,eval=FALSE}
# This will open you web browser
XenaBrowse(to_browse)

XenaBrowse(to_browse, type = "cohort")
```

```{r, error=TRUE}
# This will throw error
XenaBrowse(to_browse2)

XenaBrowse(to_browse2, type = "cohort")
```

When you make sure you want to open multiple links, you can set `multiple` option to `TRUE`.

```{r, eval=FALSE}
XenaBrowse(to_browse2, multiple = TRUE)
XenaBrowse(to_browse2, type = "cohort", multiple = TRUE)
```

## More usages

The core functionality has been described above. 
I write more usages about this package in my website but not here
because sometimes package check will fail due to internet problem.

- [Introduction and basic usage of UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-intro/) - [PDF](https://shixiangwang.github.io/home/en/tools/ucscxenatools-intro.pdf)
- [APIs of UCSCXenaTools](https://shixiangwang.github.io/home/en/tools/ucscxenatools-api/) - [PDF](https://shixiangwang.github.io/home/en/tools/ucscxenatools-api.pdf)

Read [Obtain RNAseq Values for a Specific Gene in Xena Database](https://shixiangwang.github.io/home/en/tools/ucscxenatools-single-gene/) to see how to get values for single gene. A use case for survival analysis based on single gene expression has been published on rOpenSci, please read
[UCSCXenaTools: Retrieve Gene Expression and Clinical Information from UCSC Xena for Survival Analysis](https://ropensci.org/technotes/2019/09/06/ucscxenatools-surv/).

## QA

### How to resume file from breakpoint

Thanks to the UCSC Xena team, the new feature 'resume from breakpoint' is added and 
can be done by **XenaDownload()** with the `method` and `extra` flags specified.

Of note, the corresponding `wget` or `curl` command must be installed by your OS
and can be found by R.

The folliwng code gives a test example, the data can be viewed on [web page](https://xenabrowser.net/datapages/?dataset=TcgaTargetGtex_expected_count&host=https%3A%2F%2Ftoil.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

```r
library(UCSCXenaTools)
xe = XenaGenerate(subset = XenaDatasets == "TcgaTargetGtex_expected_count")
xe
xq = XenaQuery(xe)
# You cannot resume from breakpoint in default mode
XenaDownload(xq, destdir = "~/test/", force = TRUE)
# You can do it with 'curl' command
XenaDownload(xq, destdir = "~/test/", method = "curl", extra = "-C -", force = TRUE)
# You can do it with 'wget' command
XenaDownload(xq, destdir = "~/test/", method = "wget", extra = "-c", force = TRUE)
```

## Citation

Cite me by the following paper.

```
Wang et al., (2019). The UCSCXenaTools R package: a toolkit for accessing genomics data
  from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq. 
  Journal of Open Source Software, 4(40), 1627, https://doi.org/10.21105/joss.01627

# For BibTex
  
@article{Wang2019UCSCXenaTools,
	journal = {Journal of Open Source Software},
	doi = {10.21105/joss.01627},
	issn = {2475-9066},
	number = {40},
	publisher = {The Open Journal},
	title = {The UCSCXenaTools R package: a toolkit for accessing genomics data from UCSC Xena platform, from cancer multi-omics to single-cell RNA-seq},
	url = {http://dx.doi.org/10.21105/joss.01627},
	volume = {4},
	author = {Wang, Shixiang and Liu, Xuesong},
	pages = {1627},
	date = {2019-08-05},
	year = {2019},
	month = {8},
	day = {5},
}
```

Cite UCSC Xena by the following paper. 

```
Goldman, Mary, et al. "The UCSC Xena Platform for cancer genomics data 
    visualization and interpretation." BioRxiv (2019): 326470.
```

## Acknowledgments

This package is based on [XenaR](https://github.com/mtmorgan/XenaR), thanks [Martin Morgan](https://github.com/mtmorgan) for his work.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fetch.R
\name{fetch}
\alias{fetch}
\alias{fetch_dense_values}
\alias{fetch_sparse_values}
\alias{fetch_dataset_samples}
\alias{fetch_dataset_identifiers}
\alias{has_probeMap}
\title{Fetch Data from UCSC Xena Hosts}
\usage{
fetch(host, dataset)

fetch_dense_values(
  host,
  dataset,
  identifiers = NULL,
  samples = NULL,
  check = TRUE,
  use_probeMap = FALSE,
  time_limit = 30
)

fetch_sparse_values(host, dataset, genes, samples = NULL, time_limit = 30)

fetch_dataset_samples(host, dataset, limit = NULL)

fetch_dataset_identifiers(host, dataset)

has_probeMap(host, dataset, return_url = FALSE)
}
\arguments{
\item{host}{a UCSC Xena host, like "https://toil.xenahubs.net".
All available hosts can be printed by \code{\link[=xena_default_hosts]{xena_default_hosts()}}.}

\item{dataset}{a UCSC Xena dataset, like "tcga_RSEM_gene_tpm".
All available datasets can be printed by running \code{XenaData$XenaDatasets} or
obtained from \href{https://xenabrowser.net/datapages/}{UCSC Xena datapages}.}

\item{identifiers}{Identifiers could be probe (like "ENSG00000000419.12"),
gene (like "TP53") etc.. If it is \code{NULL}, all identifiers in the dataset will be used.}

\item{samples}{ID of samples, like "TCGA-02-0047-01".
If it is \code{NULL}, all samples in the dataset will be used. However, it is better to download
the whole datasets if you query many samples and genes.}

\item{check}{if \code{TRUE}, check whether specified \code{identifiers} and \code{samples} exist the dataset
(all failed items will be filtered out). However, if \code{FALSE}, the code is much faster.}

\item{use_probeMap}{if \code{TRUE}, will check if the dataset has ProbeMap firstly.
When the dataset you want to query has a identifier-to-gene mapping, identifiers can be
gene symbols even the identifiers of dataset are probes or others.}

\item{time_limit}{time limit for getting response in seconds.}

\item{genes}{gene names.}

\item{limit}{number of samples, if \code{NULL}, return all samples.}

\item{return_url}{if \code{TRUE}, returns the info of probeMap
instead of a logical value when the result exists.}
}
\value{
a \code{matirx} or character vector or a \code{list}.
}
\description{
When you want to query just data for several genes/samples from UCSC Xena datasets, a better way
is to use these \code{fetch_} functions instead of downloading a whole dataset. Details about functions
please see the following sections.
}
\details{
There are three primary data types: dense matrix (samples by probes (or say identifiers)),
sparse (sample, position, variant), and segmented (sample, position, value).

Dense matrices can be genotypic or phenotypic, it is a sample-by-identifiers matrix.
Phenotypic matrices have associated field metadata (descriptive names, codes, etc.).
Genotypic matricies may have an associated probeMap, which maps probes to genomic locations.
If a matrix has hugo probeMap, the probes themselves are gene names. Otherwise, a probeMap is
used to map a gene location to a set of probes.
}
\section{Functions}{
\itemize{
\item \code{fetch_dense_values}: fetches values from a dense matrix.

\item \code{fetch_sparse_values}: fetches values from a sparse \code{data.frame}.

\item \code{fetch_dataset_samples}: fetches samples from a dataset

\item \code{fetch_dataset_identifiers}: fetches identifies from a dataset.

\item \code{has_probeMap}: checks if a dataset has ProbeMap.
}}

\examples{
library(UCSCXenaTools)

host <- "https://toil.xenahubs.net"
dataset <- "tcga_RSEM_gene_tpm"
samples <- c("TCGA-02-0047-01", "TCGA-02-0055-01", "TCGA-02-2483-01", "TCGA-02-2485-01")
probes <- c("ENSG00000282740.1", "ENSG00000000005.5", "ENSG00000000419.12")
genes <- c("TP53", "RB1", "PIK3CA")

\donttest{
# Fetch samples
fetch_dataset_samples(host, dataset, 2)
# Fetch identifiers
fetch_dataset_identifiers(host, dataset)
# Fetch expression value by probes
fetch_dense_values(host, dataset, probes, samples, check = FALSE)
# Fetch expression value by gene symbol (if the dataset has probeMap)
has_probeMap(host, dataset)
fetch_dense_values(host, dataset, genes, samples, check = FALSE, use_probeMap = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{downloadTCGA}
\alias{downloadTCGA}
\title{Easily Download TCGA Data by Several Options}
\usage{
downloadTCGA(
  project = NULL,
  data_type = NULL,
  file_type = NULL,
  destdir = tempdir(),
  force = FALSE,
  ...
)
}
\arguments{
\item{project}{default is \code{NULL}. Should be one or more of TCGA project id (character vector) provided by Xena.
See all available project id, please use \code{availTCGA("ProjectID")}.}

\item{data_type}{default is \code{NULL}. Should be a character vector specify data type.
See all available data types by \code{availTCGA("DataType")}.}

\item{file_type}{default is \code{NULL}. Should be a character vector specify file type.
See all available file types by \code{availTCGA("FileType")}.}

\item{destdir}{specify a location to store download data. Default is system temp directory.}

\item{force}{logical. if \code{TRUE}, force to download data no matter whether files exist.
Default is \code{FALSE}.}

\item{...}{other argument to \code{download.file} function}
}
\value{
same as \code{XenaDownload()} function result.
}
\description{
TCGA is a very useful database and here we provide this function to
download TCGA (include TCGA Pancan) datasets in human-friendly way. Users who are not
familiar with R operation will benefit from this.
}
\details{
All availble information about datasets of TCGA can access vis \code{availTCGA()} and
check with \code{showTCGA()}.
}
\examples{
\dontrun{
# download RNASeq data (use UVM as example)
downloadTCGA(project = "UVM",
                 data_type = "Gene Expression RNASeq",
                 file_type = "IlluminaHiSeq RNASeqV2")
}
}
\seealso{
\code{\link[=XenaQuery]{XenaQuery()}},
\code{\link[=XenaFilter]{XenaFilter()}},
\code{\link[=XenaDownload]{XenaDownload()}},
\code{\link[=XenaPrepare]{XenaPrepare()}},
\code{\link[=availTCGA]{availTCGA()}},
\code{\link[=showTCGA]{showTCGA()}}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaHub-class.R
\name{xena_default_hosts}
\alias{xena_default_hosts}
\title{UCSC Xena Default Hosts}
\usage{
xena_default_hosts()
}
\value{
A character vector include current defalut hosts
}
\description{
Return Xena default hosts
}
\seealso{
\code{\link[=XenaHub]{XenaHub()}}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{availTCGA}
\alias{availTCGA}
\title{Get or Check TCGA Available ProjectID, DataType and FileType}
\usage{
availTCGA(which = c("all", "ProjectID", "DataType", "FileType"))
}
\arguments{
\item{which}{a character of \code{c("All", "ProjectID", "DataType", "FileType")}}
}
\description{
Get or Check TCGA Available ProjectID, DataType and FileType
}
\examples{
\donttest{
availTCGA("all")
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaFilter.R
\name{XenaFilter}
\alias{XenaFilter}
\title{Filter a XenaHub Object}
\usage{
XenaFilter(
  x,
  filterCohorts = NULL,
  filterDatasets = NULL,
  ignore.case = TRUE,
  ...
)
}
\arguments{
\item{x}{a \link{XenaHub} object}

\item{filterCohorts}{default is \code{NULL}. A character used to filter cohorts,
regular expression is supported.}

\item{filterDatasets}{default is \code{NULL}. A character used to filter datasets,
regular expression is supported.}

\item{ignore.case}{if \code{FALSE}, the pattern matching is case sensitive
and if \code{TRUE}, case is ignored during matching.}

\item{...}{other arguments except \code{value} passed to \code{\link[base:grep]{base::grep()}}.}
}
\value{
a \code{XenaHub} object
}
\description{
One of main functions in \strong{UCSCXenatools}. It is used to filter
\code{XenaHub} object according to cohorts, datasets. All datasets can be found
at \url{https://xenabrowser.net/datapages/}.
}
\examples{
# operate TCGA datasets
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub")
xe
# get all names of clinical data
xe2 = XenaFilter(xe, filterDatasets = "clinical")
datasets(xe2)
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{showTCGA}
\alias{showTCGA}
\title{Show TCGA data structure by Project ID or ALL}
\usage{
showTCGA(project = "all")
}
\arguments{
\item{project}{a character vector. Can be "all" or one or more of TCGA Project IDs.}
}
\value{
a \code{data.frame} including project data structure information.
}
\description{
This can used to check if data type or file type exist in one or more projects by hand.
}
\examples{
\donttest{
showTCGA("all")
}
}
\seealso{
\code{\link[=availTCGA]{availTCGA()}}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simplify.R
\name{getTCGAdata}
\alias{getTCGAdata}
\title{Get TCGA Common Data Sets by Project ID and Property}
\usage{
getTCGAdata(
  project = NULL,
  clinical = TRUE,
  download = FALSE,
  forceDownload = FALSE,
  destdir = tempdir(),
  mRNASeq = FALSE,
  mRNAArray = FALSE,
  mRNASeqType = "normalized",
  miRNASeq = FALSE,
  exonRNASeq = FALSE,
  RPPAArray = FALSE,
  ReplicateBaseNormalization = FALSE,
  Methylation = FALSE,
  MethylationType = c("27K", "450K"),
  GeneMutation = FALSE,
  SomaticMutation = FALSE,
  GisticCopyNumber = FALSE,
  Gistic2Threshold = TRUE,
  CopyNumberSegment = FALSE,
  RemoveGermlineCNV = TRUE,
  ...
)
}
\arguments{
\item{project}{default is \code{NULL}. Should be one or more of TCGA project id (character vector) provided by Xena.
See all available project id, please use \code{availTCGA("ProjectID")}.}

\item{clinical}{logical. if \code{TRUE}, download clinical information. Default is \code{TRUE}.}

\item{download}{logical. if \code{TRUE}, download data, otherwise return a result list include data
information. Default is \code{FALSE}. You can set this to \code{FALSE} if you want to check what you will download or
use other function provided by \code{UCSCXenaTools} to filter result datasets you want to download.}

\item{forceDownload}{logical. if \code{TRUE}, force to download files no matter if exist. Default is \code{FALSE}.}

\item{destdir}{specify a location to store download data. Default is system temp directory.}

\item{mRNASeq}{logical. if \code{TRUE}, download mRNASeq data. Default is \code{FALSE}.}

\item{mRNAArray}{logical. if \code{TRUE}, download mRNA microarray data. Default is \code{FALSE}.}

\item{mRNASeqType}{character vector. Can be one, two or three
in \code{c("normalized", "pancan normalized", "percentile")}.}

\item{miRNASeq}{logical. if \code{TRUE}, download miRNASeq data. Default is \code{FALSE}.}

\item{exonRNASeq}{logical. if \code{TRUE}, download exon RNASeq data. Default is \code{FALSE}.}

\item{RPPAArray}{logical. if \code{TRUE}, download RPPA data. Default is \code{FALSE}.}

\item{ReplicateBaseNormalization}{logical. if \code{TRUE}, download RPPA data by Replicate Base
Normalization (RBN). Default is \code{FALSE}.}

\item{Methylation}{logical. if \code{TRUE}, download DNA Methylation data. Default is \code{FALSE}.}

\item{MethylationType}{character vector. Can be one or two in \code{c("27K", "450K")}.}

\item{GeneMutation}{logical. if \code{TRUE}, download gene mutation data. Default is \code{FALSE}.}

\item{SomaticMutation}{logical. if \code{TRUE}, download somatic mutation data. Default is \code{FALSE}.}

\item{GisticCopyNumber}{logical. if \code{TRUE}, download Gistic2 Copy Number data. Default is \code{FALSE}.}

\item{Gistic2Threshold}{logical. if \code{TRUE}, download Threshold Gistic2 data. Default is \code{TRUE}.}

\item{CopyNumberSegment}{logical. if \code{TRUE}, download Copy Number Segment data. Default is \code{FALSE}.}

\item{RemoveGermlineCNV}{logical. if \code{TRUE}, download Copy Number Segment data which has removed
germline copy number variation. Default is \code{TRUE}.}

\item{...}{other argument to \code{download.file} function}
}
\value{
if \code{download=TRUE}, return \code{data.frame} from \code{XenaDownload},
otherwise return a list including \code{XenaHub} object and datasets information
}
\description{
This is the most useful function for user to download common
TCGA datasets, it is similar to \code{getFirehoseData} function in \code{RTCGAToolbox}
package.
}
\details{
TCGA Common Data Sets are frequently used for biological analysis.
To make easier to achieve these data, this function provide really easy
options to choose datasets and behavior. All availble information about
datasets of TCGA can access vis \code{availTCGA()} and check with \code{showTCGA()}.
}
\examples{
###### get data, but not download

# 1 choose project and data types you wanna download
getTCGAdata(project = "LUAD", mRNASeq = TRUE, mRNAArray = TRUE,
mRNASeqType = "normalized", miRNASeq = TRUE, exonRNASeq = TRUE,
RPPAArray = TRUE, Methylation = TRUE, MethylationType = "450K",
GeneMutation = TRUE, SomaticMutation = TRUE)

# 2 only choose 'LUAD' and its clinical data
getTCGAdata(project = "LUAD")
\dontrun{
###### download datasets

# 3 download clinical datasets of LUAD and LUSC
getTCGAdata(project = c("LUAD", "LUSC"), clinical = TRUE, download = TRUE)

# 4 download clinical, RPPA and gene mutation datasets of LUAD and LUSC
# getTCGAdata(project = c("LUAD", "LUSC"), clinical = TRUE, RPPAArray = TRUE, GeneMutation = TRUE)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{XenaData}
\alias{XenaData}
\title{Xena Hub Information}
\format{
A \code{tibble}.
}
\source{
Generated from UCSC Xena Data Hubs.
}
\description{
This \code{data.frame} is very useful for selecting datasets fastly and
independent on APIs of UCSC Xena Hubs.
}
\examples{
data(XenaData)
str(XenaData)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaPrepare.R
\name{XenaPrepare}
\alias{XenaPrepare}
\title{Prepare (Load) Downloaded Datasets to R}
\usage{
XenaPrepare(
  objects,
  objectsName = NULL,
  use_chunk = FALSE,
  chunk_size = 100,
  subset_rows = TRUE,
  select_cols = TRUE,
  callback = NULL,
  comment = "#",
  na = c("", "NA", "[Discrepancy]"),
  ...
)
}
\arguments{
\item{objects}{a object of character vector or data.frame. If \code{objects} is data.frame,
it should be returned object of \link{XenaDownload} function. More easier way is
that objects can be character vector specify local files/directory and download urls.}

\item{objectsName}{specify names for elements of return object, i.e. names of list}

\item{use_chunk}{default is \code{FALSE}. If you want to select subset of original data, please set it to
\code{TRUE} and specify corresponding arguments: \code{chunk_size}, \code{select_direction}, \code{select_names},
\code{callback}.}

\item{chunk_size}{the number of rows to include in each chunk}

\item{subset_rows}{logical expression indicating elements or rows to keep:
missing values are taken as false. \code{x} can be a representation of data frame
you wanna do subset operation. Of note, the first colname of most of datasets
in Xena will be set to "sample", you can use it to select rows.}

\item{select_cols}{expression, indicating columns to select from a data frame.
'x' can be a representation of data frame you wanna do subset operation,
e.g. \code{select_cols = colnames(x)[1:3]} will keep only first to third column.}

\item{callback}{a function to call on each chunk, default is \code{NULL},
this option will overvide operations of subset_rows and select_cols.}

\item{comment}{a character specify comment rows in files}

\item{na}{a character vectory specify \code{NA} values in files}

\item{...}{other arguments transfer to \code{read_tsv} function or
\code{read_tsv_chunked} function (when \code{use_chunk} is \code{TRUE}) of \code{readr} package.}
}
\value{
a list contains file data, which in way of tibbles
}
\description{
Prepare (Load) Downloaded Datasets to R
}
\examples{
\dontrun{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub")
hosts(xe)
xe_query = XenaQuery(xe)

xe_download = XenaDownload(xe_query)
dat = XenaPrepare(xe_download)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaQuery.R
\name{XenaQuery}
\alias{XenaQuery}
\title{Query URL of Datasets before Downloading}
\usage{
XenaQuery(x)
}
\arguments{
\item{x}{a \link{XenaHub} object}
}
\value{
a \code{data.frame} contains hosts, datasets and url
}
\description{
Query URL of Datasets before Downloading
}
\examples{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub")
hosts(xe)
\dontrun{
xe_query = XenaQuery(xe)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaHub-class.R
\name{XenaDataUpdate}
\alias{XenaDataUpdate}
\title{Get or Update Newest Data Information of UCSC Xena Data Hubs}
\usage{
XenaDataUpdate(saveTolocal = TRUE)
}
\arguments{
\item{saveTolocal}{logical. Whether save to local R package data directory for permanent use
or Not.}
}
\value{
a \code{data.frame} contains all datasets information of Xena.
}
\description{
Get or Update Newest Data Information of UCSC Xena Data Hubs
}
\examples{
\dontrun{
XenaDataUpdate()
XenaDataUpdate(saveTolocal = TRUE)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_higher.R
\name{cohorts}
\alias{cohorts}
\title{Get cohorts of XenaHub object}
\usage{
cohorts(x)
}
\arguments{
\item{x}{a \link{XenaHub} object}
}
\value{
a character vector contains cohorts
}
\description{
Get cohorts of XenaHub object
}
\examples{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub"); cohorts(xe)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_higher.R
\name{datasets}
\alias{datasets}
\title{Get datasets of XenaHub object}
\usage{
datasets(x)
}
\arguments{
\item{x}{a \link{XenaHub} object}
}
\value{
a character vector contains datasets
}
\description{
Get datasets of XenaHub object
}
\examples{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub"); datasets(xe)
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
See \code{magrittr::\link[magrittr]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaDownload.R
\name{XenaDownload}
\alias{XenaDownload}
\title{Download Datasets from UCSC Xena Hubs}
\usage{
XenaDownload(
  xquery,
  destdir = tempdir(),
  download_probeMap = FALSE,
  trans_slash = FALSE,
  force = FALSE,
  max_try = 3L,
  ...
)
}
\arguments{
\item{xquery}{a tibble object generated by \link{XenaQuery} function.}

\item{destdir}{specify a location to store download data. Default is system temp directory.}

\item{download_probeMap}{if \code{TRUE}, also download ProbeMap data, which used for id mapping.}

\item{trans_slash}{logical, default is \code{FALSE}. If \code{TRUE}, transform slash '/' in dataset id
to '__'. This option is for backwards compatibility.}

\item{force}{logical. if \code{TRUE}, force to download data no matter whether files exist.
Default is \code{FALSE}.}

\item{max_try}{time limit to try downloading the data.}

\item{...}{other argument to \code{download.file} function}
}
\value{
a \code{tibble}
}
\description{
Avaliable datasets list: \url{https://xenabrowser.net/datapages/}
}
\examples{
\dontrun{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub")
hosts(xe)
xe_query = XenaQuery(xe)
xe_download = XenaDownload(xe_query)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaHub-class.R
\docType{class}
\name{XenaHub-class}
\alias{XenaHub-class}
\alias{.XenaHub}
\title{Class XenaHub}
\description{
a S4 class to represent UCSC Xena Data Hubs
}
\section{Slots}{

\describe{
\item{\code{hosts}}{hosts of data hubs}

\item{\code{cohorts}}{cohorts of data hubs}

\item{\code{datasets}}{datasets of data hubs}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shiny.R
\name{XenaShiny}
\alias{XenaShiny}
\title{Xena Shiny App}
\usage{
XenaShiny()
}
\description{
Xena Shiny App
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaScan.R
\name{XenaScan}
\alias{XenaScan}
\title{Scan all rows according to user input by a regular expression}
\usage{
XenaScan(
  XenaData = UCSCXenaTools::XenaData,
  pattern = NULL,
  ignore.case = TRUE
)
}
\arguments{
\item{XenaData}{a \code{data.frame}. Default is \code{data(XenaData)}.
The input of this option can only be \code{data(XenaData)} or its subset.}

\item{pattern}{character string containing a \link[base]{regular expression}
    (or character string for \code{fixed = TRUE}) to be matched
    in the given character vector.  Coerced by
    \code{\link[base]{as.character}} to a character string if possible.  If a
    character vector of length 2 or more is supplied, the first element
    is used with a warning.  Missing values are allowed except for
    \code{regexpr}, \code{gregexpr} and \code{regexec}.}

\item{ignore.case}{if \code{FALSE}, the pattern matching is \emph{case
      sensitive} and if \code{TRUE}, case is ignored during matching.}
}
\value{
a \code{data.frame}
}
\description{
\code{XenaScan()} is a function can be used before \code{\link[=XenaGenerate]{XenaGenerate()}}.
}
\examples{

x1 <- XenaScan(pattern = "Blood")
x2 <- XenaScan(pattern = "LUNG", ignore.case = FALSE)

x1 \%>\%
  XenaGenerate()
x2 \%>\%
  XenaGenerate()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_higher.R
\name{hosts}
\alias{hosts}
\title{Get hosts of XenaHub object}
\usage{
hosts(x)
}
\arguments{
\item{x}{a \link{XenaHub} object}
}
\value{
a character vector contains hosts
}
\description{
Get hosts of XenaHub object
}
\examples{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub"); hosts(xe)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_higher.R
\name{samples}
\alias{samples}
\title{Get Samples of a XenaHub object according to 'by' and 'how' action arguments}
\usage{
samples(
  x,
  i = character(),
  by = c("hosts", "cohorts", "datasets"),
  how = c("each", "any", "all")
)
}
\arguments{
\item{x}{a \link{XenaHub} object}

\item{i}{default is a empty character, it is used to specify
the host, cohort or dataset by \code{by} option otherwise
info will be automatically extracted by code}

\item{by}{a character specify \code{by} action}

\item{how}{a character specify \code{how} action}
}
\value{
a list include samples
}
\description{
One is often interested in identifying samples or features present in each data set,
or shared by all data sets, or present in any of several data sets.
Identifying these samples, including samples in arbitrarily chosen data sets.
}
\examples{
\dontrun{
xe = XenaHub(cohorts = "Cancer Cell Line Encyclopedia (CCLE)")
# samples in each dataset, first host
x = samples(xe, by="datasets", how="each")[[1]]
lengths(x)        # data sets in ccle cohort on first (only) host
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaHub-class.R
\name{XenaHub}
\alias{XenaHub}
\title{Generate a XenaHub Object}
\usage{
XenaHub(
  hosts = xena_default_hosts(),
  cohorts = character(),
  datasets = character(),
  hostName = c("publicHub", "tcgaHub", "gdcHub", "icgcHub", "toilHub",
    "pancanAtlasHub", "treehouseHub", "pcawgHub", "atacseqHub", "singlecellHub",
    "kidsfirstHub")
)
}
\arguments{
\item{hosts}{a character vector specify UCSC Xena hosts, all available hosts can be
found by \code{xena_default_hosts()} function. \code{hostName} is a more recommend option.}

\item{cohorts}{default is empty character vector, all cohorts will be returned.}

\item{datasets}{default is empty character vector, all datasets will be returned.}

\item{hostName}{name of host, available options can be accessed by \code{.xena_hosts}
This is an easier option for user than \code{hosts} option. Note, this option
will overlap \code{hosts}.}
}
\value{
a \link{XenaHub} object
}
\description{
It is used to generate original
\code{XenaHub} object according to hosts, cohorts, datasets or hostName.
If these arguments not specified, all hosts and corresponding datasets
will be returned as a \code{XenaHub} object. All datasets can be found
at \url{https://xenabrowser.net/datapages/}.
}
\examples{
\dontrun{
#1 query all hosts, cohorts and datasets
xe = XenaHub()
xe
#2 query only TCGA hosts
xe = XenaHub(hostName = "tcgaHub")
xe
hosts(xe)     # get hosts
cohorts(xe)   # get cohorts
datasets(xe)  # get datasets
samples(xe)   # get samples
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/api_xq.R
\name{to_snake}
\alias{to_snake}
\title{Convert camel case to snake case}
\usage{
to_snake(name)
}
\arguments{
\item{name}{a character vector}
}
\value{
same length as \code{name} but with snake case
}
\description{
Convert camel case to snake case
}
\examples{
to_snake("sparseDataRange")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaGenerate.R
\name{XenaGenerate}
\alias{XenaGenerate}
\title{Generate and Subset a XenaHub Object from 'XenaData'}
\usage{
XenaGenerate(XenaData = UCSCXenaTools::XenaData, subset = TRUE)
}
\arguments{
\item{XenaData}{a \code{data.frame}. Default is \code{data(XenaData)}.
The input of this option can only be \code{data(XenaData)} or its subset.}

\item{subset}{logical expression indicating elements or rows to keep.}
}
\value{
a \link{XenaHub} object.
}
\description{
Generate and Subset a XenaHub Object from 'XenaData'
}
\examples{
# 1 get all datasets
XenaGenerate()
# 2 get TCGA BRCA
XenaGenerate(subset = XenaCohorts == "TCGA Breast Cancer (BRCA)")
# 3 get all datasets containing BRCA
XenaGenerate(subset = grepl("BRCA", XenaCohorts))
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaBrowse.R
\name{XenaBrowse}
\alias{XenaBrowse}
\title{View Info of Dataset or Cohort at UCSC Xena Website Using Web browser}
\usage{
XenaBrowse(x, type = c("dataset", "cohort"), multiple = FALSE)
}
\arguments{
\item{x}{a \link{XenaHub} object.}

\item{type}{one of "dataset" and "cohort".}

\item{multiple}{if \code{TRUE}, browse multiple links instead of throwing error.}
}
\description{
This will open dataset/cohort link of UCSC Xena
in user's default browser.
}
\examples{
\donttest{
XenaGenerate(subset = XenaHostNames == "tcgaHub") \%>\%
  XenaFilter(filterDatasets = "clinical") \%>\%
  XenaFilter(filterDatasets = "LUAD") -> to_browse
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/XenaQueryProbeMap.R
\name{XenaQueryProbeMap}
\alias{XenaQueryProbeMap}
\title{Query ProbeMap URL of Datasets}
\usage{
XenaQueryProbeMap(x)
}
\arguments{
\item{x}{a \link{XenaHub} object}
}
\value{
a \code{data.frame} contains hosts, datasets and url
}
\description{
If dataset has no ProbeMap, it will be ignored.
}
\examples{
xe = XenaGenerate(subset = XenaHostNames == "tcgaHub")
hosts(xe)
\dontrun{
xe_query = XenaQueryProbeMap(xe)
}
}
\author{
Shixiang Wang \href{mailto:w_shixiang@163.com}{w_shixiang@163.com}
}
