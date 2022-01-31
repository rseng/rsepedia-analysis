[![R-CMD-check](https://github.com/timoast/signac/workflows/R-CMD-check/badge.svg)](https://github.com/timoast/signac/actions)
[![CRAN Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

# Signac

Signac is an extension of [Seurat](https://github.com/satijalab/seurat) for the analysis, interpretation, and exploration of single-cell chromatin datasets.

## Features

Signac is designed for the analysis of single-cell chromatin data, including scATAC-seq,
single-cell targeted tagmentation methods such as scCUT&Tag and scACT-seq,
and multimodal datasets that jointly measure chromatin state alongside other
modalities.

Signac currently supports the following features:

* Calling peaks
* Quantifying per-cell counts in different genomic regions
* Calculating single-cell QC metrics
* Dimensional reduction, visualization, and clustering
* Identifying cell-type-specific peaks
* Visualizing 'pseudo-bulk' coverage tracks
* Integration of multiple single-cell datasets
* Integration with single-cell RNA-seq datasets
* Sequence motif enrichment analysis
* Transcription factor footprinting analysis
* Linking peaks to correlated genes
* Parallelization through the [future](https://cran.r-project.org/package=future) package
* Seamless interface with [Seurat](https://satijalab.org/seurat), [SeuratWrappers](https://github.com/satijalab/seurat-wrappers), [SeuratDisk](https://github.com/mojaveazure/seurat-disk), and [SeuratData](https://github.com/satijalab/seurat-data) functionality
* Interoperability with [Bioconductor](https://bioconductor.org/) tools

Please see the Signac [vignettes](articles/overview.html) page for examples.

For installation instructions see the [install](articles/install.html) page.


# Signac

[![R-CMD-check](https://github.com/timoast/signac/workflows/R-CMD-check/badge.svg)](https://github.com/timoast/signac/actions)
[![CRAN
Version](https://www.r-pkg.org/badges/version/Signac)](https://cran.r-project.org/package=Signac)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/Signac)](https://cran.r-project.org/package=Signac)

## Overview

Signac is a comprehensive R package for the analysis of single-cell
chromatin data. Signac includes functions for quality control,
normalization, dimension reduction, clustering, differential activity,
and more.

Documentation and tutorials can be found at
<https://satijalab.org/signac/>

## Installation

Signac requires that [Bioconductor](https://www.bioconductor.org/) is
installed:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
setRepositories(ind=1:2)
```

To install the latest release of Signac from CRAN:

``` r
install.packages("Signac")
```

To release the latest develop version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("timoast/signac", ref = "develop")
```

## Release notes

For a changelog please see the [NEWS
file](https://github.com/timoast/signac/blob/develop/NEWS.md), also
available on the [Signac
website](https://satijalab.org/signac/news/index.html).

## Getting help

If you encounter a bug or have a feature request, please open an
[issue](https://github.com/timoast/signac/issues).

If you would like to discuss questions related to single-cell analysis,
you can open a
[discussion](https://github.com/timoast/signac/discussions).

## Citing Signac

If you use the Signac package in your work please cite [Stuart et
al. 2021](https://doi.org/10.1038/s41592-021-01282-5)

```
@ARTICLE{signac,
  title     = "Single-cell chromatin state analysis with Signac",
  author    = "Stuart, Tim and Srivastava, Avi and Madad, Shaista and Lareau,
               Caleb A and Satija, Rahul",
  journal   = "Nat. Methods",
  publisher = "Nature Publishing Group",
  pages     = "1--9",
  month     =  nov,
  year      =  2021,
  url       = "https://www.nature.com/articles/s41592-021-01282-5",
  language  = "en"
}
```

## Related packages

-   [Seurat](https://github.com/satijalab/seurat)
-   [SeuratObject](https://github.com/mojaveazure/seurat-object)
-   [SeuratDisk](https://github.com/mojaveazure/seurat-disk)
-   [SeuratData](https://github.com/satijalab/seurat-data)
-   [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
-   [Azimuth](https://github.com/satijalab/azimuth)
# Signac 1.5.0

Bug fixes:

* Fixed bug in `FeatureMatrix()` when cells information not present in Fragment object ([#803](https://github.com/timoast/signac/issues/803))
* Fixed bug in object merging ([#804](https://github.com/timoast/signac/issues/804))
* Add ability to run `LinkPeaks()` using Ensembl IDs ([#858](https://github.com/timoast/signac/issues/858))
* Fix issue in `GeneActivity()` when gene names are `NA` ([#865](https://github.com/timoast/signac/issues/865))
* Fix bug in `FeatureMatrix()` when only one region supplied
* Allow negative values in `ExpressionPlot()` when using scaled data ([#893](https://github.com/timoast/signac/issues/893))

Other changes:

* Added `idf` parameter to `RunTFIDF()` to use precomputed IDF vector
* Added `gene.id` parameter to `GeneActivity()` to allow output genes named using gene ID ([#837](https://github.com/timoast/signac/issues/837))
* Added `sep` parameter to `ConnectionsToLinks()` ([#841](https://github.com/timoast/signac/issues/841))

# Signac 1.4.0

Bug fixes:

* Fixed bug in `FindMotifs()` when using only one region as input ([#732](https://github.com/timoast/signac/issues/732))
* Add check for correct number of columns in fragment file ([#748](https://github.com/timoast/signac/issues/748))
* Fixed gene lookup when annotations contain NA values ([#771](https://github.com/timoast/signac/issues/771))
* Fixed error in `ClosestFeature()` when query contained regions on contigs not present in gene annotation ([#758](https://github.com/timoast/signac/issues/758))
* Fixed bug in `TSSEnrichment()` when using multiple fragment files ([#783](https://github.com/timoast/signac/issues/783))
* Fixed bug in `CallPeaks()` when multiple fragment files used as input
* Fixed bug in `CallPeaks()` to account for 0-based starts in called peaks
* Fixed bug in gene name lookup when gene names contain `-` characters ([#759](https://github.com/timoast/signac/issues/759))

Other changes:

* Updated documentation for `genome` parameter in `AddMotifs()` and `RunChromVAR()` ([#712](https://github.com/timoast/signac/issues/712))
* Updated the `FoldChange()` function to use normalized counts rather than raw counts ([#795](https://github.com/timoast/signac/issues/795))
* Improved error checking in `GeneActivity()` ([#797](https://github.com/timoast/signac/issues/797))
* Added `format` parameter to `CallPeaks()` ([#682](https://github.com/timoast/signac/issues/682))

# Signac 1.3.0 

Bug fixes:

* Fixed `LinkPeaks()` function when running on a single gene ([#629](https://github.com/timoast/signac/issues/629))
* Added `fragment.tempdir` parameter to `CallPeaks()` to enable setting directory
that split fragment files are written to during peak calling ([#579](https://github.com/timoast/signac/issues/579))
* Fixed error in `FeatureMatrix()` when setting `sep` parameter ([#626](https://github.com/timoast/signac/discussions/626))
* Fixed peak calling error when group names contain special characters
* Fixed issue with `RenameCells()` when cell information not present in Fragment object ([#704](https://github.com/timoast/signac/issues/704))

Other changes: 

* Improved error checking for `GeneActivity()` ([#625](https://github.com/timoast/signac/issues/625))
* Added `FoldChange()` method for `ChromatinAssay()` object that sets proper parameters for 
chromatin data. This fixes the calculation of fold changes when running `Seurat::FindMarkers()` on
single-cell chromatin data.

# Signac 1.2.1

New functionality:

* Added `head()` method for `Fragment`-class objects.

Bug fixes:

* Fixed bug in `ChromatinAssay` merging ([#596](https://github.com/timoast/signac/pull/596))

Other changes:

* Added support for fragment files containing headers (cellranger-atac v2; [#609](https://github.com/timoast/signac/issues/609))

# Signac 1.2.0

New functionality:

* Added `BigwigTrack()` function to plot data from bigWig files
* Added `bigwig` and `bigwig.type` arguments to `CoveragePlot()` to
include bigWig files in `CoveragePlot()`
* Added `region.highlight` parameter to `CoveragePlot()`
* Added `biotypes` parameter to `GeneActivity()` and `GetTSSPositions()` functions
* Added `max.width` parameter to `GeneActivity()`
* Added `min.distance` parameter to `LinkPeaks()` ([#561](https://github.com/timoast/signac/pull/561))

Bug fixes:

* Fixed fragment file reading when only one fragment found in requested region ([#474](https://github.com/timoast/signac/issues/474))
* Fixed `standard.chromosomes` parameter in `GetGRangesFromEnsDb()` ([#513](https://github.com/timoast/signac/issues/513))
* Fixed `group.by` parameter in `PlotFootprint()` ([#522](https://github.com/timoast/signac/issues/522))
* Fixed bug that would cause some gene coordinates used by `GeneActivity()` to be 
incorrect ([#521](https://github.com/timoast/signac/issues/521))
* Fixed error message in `FindMotifs()` ([#549](https://github.com/timoast/signac/issues/549))
* Fixed bug in `CountsInRegion()` ([#563](https://github.com/timoast/signac/issues/563))

Other changes:

* Improved speed of ChromatinAssay merging
* Improved error message for `TSSEnrichment()` ([#485](https://github.com/timoast/signac/issues/485))
* Improved error messages when trying to run `ChromatinAssay`-specific functions
on non-`ChromatinAssay` assays
* Performance improvements
* Changed default value for `n` in `NucleosomeSignal()`
* Enabled parallization in `TSSEnrichment()` when `fast=TRUE`
* Added early error checking in `LinkPeaks()` ([#550](https://github.com/timoast/signac/pull/550))
* Change to sparse matrix correlation in `LinkPeaks()` ([#550](https://github.com/timoast/signac/pull/550))
* Moved `biovizBase` and `Biostrings` to suggested packages
* Removed `ggbio` dependency
* Re-implemented `AnnotationPlot()`

# Signac 1.1.1 

New functionality:

* Added `group.by` parameter to `PeakPlot()` to allow coloring plotted genomic 
ranges by metadata variables.
* Added `peaks.group.by` and `ranges.group.by` parameters to `CoveragePlot()` to
allow coloring plotted genomic ranges in `CoveragePlot()` to be colored by metadata
variables.

Bug fixes:

* Update meta feature information (overall peak accessibility) when subsetting 
objects to avoid counts becoming inaccurate ([#332](https://github.com/timoast/signac/issues/332))
* Prevent dropping features when creating a merged ChromatinAssay ([#340](https://github.com/timoast/signac/pull/340))
* Fix compilation error when using g++ version <5 ([#326](https://github.com/timoast/signac/issues/326))
* Retain motif positions during subset ([#364](https://github.com/timoast/signac/issues/364))
* Fix `assay` parameter in `CoveragePlot()`
* Fix error when merging ChromatinAssay object ([#355](https://github.com/timoast/signac/issues/355))
* Add more informative error message when all features or cells removed by parameter choices in `CreateChromatinAssay()` ([#387](https://github.com/timoast/signac/issues/387))
* Fix bug in `CreateChromatinAssay()` when setting both `min.cells` and `min.features` arguments ([#390](https://github.com/timoast/signac/issues/390))
* Improved support for remote fragment files
* Fixed bug in `PlotFootprint()` when only one cell in an identity class ([#406](https://github.com/timoast/signac/issues/406))

Other changes:

* Added citation information to the package
* Added `SeuratObject` dependency

# Signac 1.1.0

New functionality:

* Added `CallPeaks()` function to call peaks using MACS2. Peaks can be called
for different groups of cells separately by setting the `group.by` parameter
* Added `LinkPeaks()` function to link peaks to correlated genes.
* Added `AddMotifs()` function to add motif information to a Seurat object or ChromatinAssay.
* Added `AggregateTiles()` function to combine adjacent genome tiles
* Added `ranges` parameter to `CoveragePlot()` to plot addition sets of genomic ranges
* Added `show.bulk` parameter to `CoveragePlot()` to plot accessibility of all cells combined
* Added ability to remove `Fragment` objects and modify the file path for existing
fragment objects ([#206](https://github.com/timoast/signac/issues/206))

Bug fixes: 

* Fixed bugs in `AlleleFreq()` ([#196](https://github.com/timoast/signac/issues/196)
and [#260](https://github.com/timoast/signac/issues/260))
* Fixed bug in `FeatureMatrix()` ([#205](https://github.com/timoast/signac/issues/205), [#291](https://github.com/timoast/signac/issues/291))
* Fixed bug in `CreateChromatinAssay()` when setting `min.features` argument ([#194](https://github.com/timoast/signac/issues/194))
* Fixed bug in `CreateChromatinAssay()` when setting `min.cells` argument ([#292](https://github.com/timoast/signac/issues/292))
* Fixed bug in `TSSEnrichment()` when cell information not set for fragment files ([#203](https://github.com/timoast/signac/issues/203))
* Fixed bug in `TSSEnrichment()` when no fragments present in TSS region ([#244](https://github.com/timoast/signac/issues/244))
* Removed `qvalue` calculation from `FindMotifs()` ([#223](https://github.com/timoast/signac/issues/223))
* Fixed bug in `SetAssayData()` when setting the `scale.data` slot

Other changes:

* Improved feature matching in `MatchRegionStats()` function when matching distribution of multiple features (eg, GC content and overall accessibility)
* Changed parameter names in `MatchRegionStats()`

# Signac 1.0.0

This release includes major updates to the Signac package, including new
functionality, performance improvements, and new data structures.

The entire package has been updated to use the new `ChromatinAssay` class for the
storage of single-cell chromatin data. This is an extension of the standard 
Seurat `Assay` that adds additional slots needed for the analysis of chromatin 
data, including genomic ranges, genome information, fragment file information,
motifs, gene annotations, and genomic links.

In addition, we have defined a new `Fragment` class to store information 
relating to a fragment file. This makes use of the fragment files within Signac
more robust, as checks are now performed to verify that the expected cells are
present in the fragment file, and that the fragment file or index are not
modified on disk.

Key new functionality:

* **Store multiple fragment files**: you can now store as many fragment
files as needed in a single object, and all functions that use the fragment file
will pull data from each of the files. Cell barcodes in the fragment files do
_not_ need to match the cell barcodes in the object.
* **Use remote fragment files**: you can now use all the same functionality with
fragment files hosted on remote servers accessible through `http` or `ftp`.
* **Transcription factor footprinting**: New `Footprint()` and `PlotFootprint()`
functions for TF footprinting analysis.
* **Bioconductor methods**: call `granges()`, `findOverlaps()`, `seqinfo()`, and
other Bioconductor generic functions directly on the `ChromatinAssay` or
`Seurat` object.
* **New multi-modal visualization methods**: Jointly visualize RNA expression
and chromatin accessibility using the `CoveragePlot()` function.
* **New interactive visualizations**: Interactively browse the genome using the
`CoverageBrowser()` function.
* **Mitochondrial lineage tracing**: New functions to identify informative
mitochondrial alleles, find clonotypes, and predict cell lineage relationships
using mitochondrial mutations.

Other changes:

* Updates to `NucleosomeSignal()`: we have greatly improved the scalability of
`NucleosomeSignal()`, and fixed a bug present in previous versions. The score
computed by `NucleosomeSignal()` in 1.0.0 will be different to that computed by
previous versions of Signac.
* New `CountFragments()` function: a fast, memory-efficient function implemented
in C++ that counts the total number of fragments for each cell barcode present
in a fragment file.
* New `fast` option in the `TSSEnrichment()` function. Setting this to `TRUE`
will compute the TSS enrichment score per cell without storing the entire
cell by TSS position matrix. This can significantly reduce memory requirements
for large datasets, but does not allow subsequent plotting of the TSS signal
for different groups of cells.
* New `TilePlot()` function and `tile` parameter for `CoveragePlot()` to plot
Tn5 integration events in a genomic region for individual cells.
* Performance improvements for `FeatureMatrix()`, `CoveragePlot()`, and
`TSSEnrichment()`
* Added the manually curated hg38 genomic blacklist regions curated by Anshul
Kundaje and Anna Shcherbina. These are available as the `blacklist_hg38_unified`
object.
* Updated the `FRiP()` function to use total fragment counts per cell stored
in object metadata.

# Signac 0.2.5

* New `DepthCor` function to compute the correlation between sequencing depth and
reduced dimension components. 
* Performance improvements for `RunTFIDF`. 
* Removed option to use EnsDb object in `ClosestFeatures` and `CoveragePlot`. Use GRanges instead. 
* Removed `ucsc` parameter from `CoveragePlot`. 
* Fixed bug in FeatureMatrix that would cause fragments to be counted multiple
times if `nchunk` was greater than the number of features used. 
* Fixed bug in `CoveragePlot` that would prevent plotting multiple regions when
using `GRanges`. 
* Fixed bug in `CoveragePlot` that would prevent plotting when a different 
assay was active. 
* Removed dependencies: GenomicFeatures
* Moved dependencies to suggests: Biostrings, BSgenome
* Removed from suggests: BSgenome.Hsapiens.UCSC.hg19, EnsDb.Hsapiens.v75, JASPAR2018

# Signac 0.2.4

* First CRAN release.
* New `SubsetMatrix` function to subset a matrix based on number of non-zero elements 
in the rows or columns.
* Removed `seed.use` parameter from `RunSVD`.

# Signac 0.2.3

* New `UnifyPeaks` function to create a merged set of peaks from multiple samples.

# Signac 0.2.2

* Bug fix for `RunSVD`: previously, scaling was applied to each cell rather than each component.
Now, mean centering and SD scaling are applied to the cell embeddings within a component.
* Added `scale.embeddings` option to `RunSVD` to control whether embeddings are scaled
and centered.
* Added `irlba.work` parameter to `RunSVD`.
* Update to allow comment characters in fragment file cell names

# Signac 0.2.1

* Removed `SingleCoveragePlot` from exported functions
* Added executable examples for all functions
* Store raw SVD output in DimReduc misc slot in `RunSVD`
* Fixed strand orientation for gene plot in `CoveragePlot`
* Fix missing x-axis when plotting peaks but not genes in `CoveragePlot`

# Signac 0.2.0

* Removed dependency on TFBSTools, motifmatchr, AnnotationDbi, ggbio, AnnotationFilter
* Renamed `PeriodPlot` to `FragmentHistogram`
* Removed motif dimension reduction functions
* Removed motif clustering functions
* Removed `neighbors` and `reductions` slots from `motif` class
* Added `motif.names` slot to `motif` class
* Added ability to plot peak ranges in `CoveragePlot`
* Added ability to plot gene annotations from `GRanges` object
* Changed gene plot style in `CoveragePlot`
* Allow passing additional arguments to `FilterFragments`
* Add inst/extdata 
* Change DESCRIPTION file so that Bioconductor dependencies are automatically installed

# Signac 0.1.6

* Bug fix for `GetCellsInRegion`
* Improve documentation

# Signac 0.1.5

* New `TSSEnrichment` and `TSSPlot` functions for TSS enrichment scoring
* New `InsertionBias` function
* New options in `CoveragePlot` for scaling tracks 
* Major speed improvements for `CoveragePlot` 
* Improved documentation (added examples)

# Signac 0.1.4

* Updates to `CoveragePlot`: now plots a Tn5 integration score per base, rather than the whole fragment.

# Signac 0.1.3

* New `GetIntersectingFeatures` function to find overlapping peaks between objects  
* New `MergeWithRegions` function to perform region-aware Seurat object merging  

# Signac 0.1.2

* New `RunChromVAR` function to run chromVAR through Signac  
* New `RegionStats` function to add statistics about peak sequences to the feature metadata  
* Improvements to `FindMotifs`: now selects a set of background peaks matching the sequence characteristics of the input

# Signac 0.1.1

* Added `IntersectMatrix`  
* Added unit tests  
* Bug fixes for `ChunkGRanges`

# Signac 0.1.0

* This is the first release of Signac!
---
name: "\U0001F4DADocumentation"
about: An issue related to https://satijalab.org/signac/
title: ''
labels: documentation
assignees: ''

---

<!-- A clear description of what content at https://satijalab.org/signac or in the Signac function man pages is an issue. -->
---
name: "\U0001F41BBug report"
about: Describe a bug you've seen with a reproducible example
title: ''
labels: bug
assignees: ''

---

<!-- Briefly describe your problem and what output you expect. If you have a question, please use the analysis question template instead. -->

<!-- Before posting an issue, ensure that the bug is reproducible by re-running the code that produced the issue in a new R session.-->

<!-- Please include a minimal reproducible code example. You can use the small test data included in Signac (`atac_small`) to demonstrate the issue, or a public dataset (for example, a dataset used in the Signac vignettes: https://satijalab.org/signac/articles/). If you cannot reproduce the issue using a public dataset, please still provide code that reproduces the issue on your data and we will try to address it. -->

<!-- Please include the output of `sessionInfo()` and your operating system in your issue. -->

```r
# insert reproducible example here
```
---
name: "\U0001F680Feature request"
about: Request a new feature
title: ''
labels: enhancement
assignees: ''

---

<!-- clearly explain a new feature you would like and why this would be useful. -->
---
title: "Installation"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

To use Signac first make sure Bioconductor is installed:

```{r}
# Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# To automatically install Bioconductor dependencies
setRepositories(ind=1:2)
```

## Current release

```{r}
install.packages("Signac")
```

## Development version

Unreleased versions of Signac can be installed from the GitHub repository using
the devtools package:

```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("timoast/signac", ref = "develop")
```

## Old versions

Older Signac releases can be installed from the CRAN archive using devtools:

```{r}
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# replace "0.2.5" with the version that you want to install
devtools::install_version(package = 'Signac', version = package_version('0.2.5'))
```

## Docker

We provide docker images for Signac via [dockerhub](https://hub.docker.com/r/timoast/signac).

To pull the latest image from the command line:

```sh
docker pull timoast/signac:latest
```

To use as a base image in a new Dockerfile:

```sh
FROM timoast/signac:latest
```

## Conda

Signac can also be installed using conda. Note that if you use conda, you 
should install all packages through conda rather than R itself. Make sure to 
set up the conda channels first: http://bioconda.github.io/user/install.html#set-up-channels

```{bash}
conda install -c bioconda r-signac
```

## Installing genome assembly and gene annotation packages

It can also be useful (but not essential) to install species-specific packages 
containing genome and gene annotation information from Bioconductor.

[This](https://useast.ensembl.org/info/website/archives/assembly.html) table
from Ensembl provides a mapping of genome assembly to the corresponding gene
annotation version.

### Human hg19

```{r}
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg19', 'EnsDb.Hsapiens.v75'))
```

### Human hg38

```{r}
BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
```

### Mouse mm10

```{r}
BiocManager::install(c('BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79'))
```
---
title: "scATAC-seq data integration"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Here we demonstrate the integration of multiple single-cell chromatin datasets
derived from human PBMCs. One dataset was generated using the 10x Genomics
multiome technology, and includes DNA accessibility and gene expression information
for each cell. The other dataset was profiled using 10x Genomics scATAC-seq,
and includes DNA accessibility data only.

We will integrate the two datasets together using the shared DNA accessibility
assay, using tools available in the Seurat package. Furthermore, we will
demonstrate transferring both continuous (gene expression) and categorical
(cell labels) information from a reference to a query single-cell chromatin
dataset.

<details>
  <summary>**View data download code**</summary>

The PBMC multiome and scATAC-seq data can be downloaded from the 10x website:

```{bash, eval=FALSE}
# multiome
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

# scATAC
wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz.tbi
```

</details>

## Preprocessing

Here we'll load the PBMC multiome data pre-processed in our
[multiome vignette](articles/pbmc_multiomic.html), and create a new object from
the scATAC-seq data:

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(ggplot2)

# load the pre-processed multiome data
pbmc.multi <- readRDS("../vignette_data/pbmc_multiomic.rds")

# process the scATAC data
# first count fragments per cell
fragpath <- "../vignette_data/atac_pbmc_10k_nextgem_fragments.tsv.gz"
fragcounts <- CountFragments(fragments = fragpath)
atac.cells <- fragcounts[fragcounts$frequency_count > 2000, "CB"]

# create the fragment object
atac.frags <- CreateFragmentObject(path = fragpath, cells = atac.cells)
```

An important first step in any integrative analysis of single-cell chromatin data
is to ensure that the same features are measured in each dataset. Here, we
quantify the multiome peaks in the ATAC dataset to ensure that there are common features
across the two datasets.

```{r message=FALSE, warning=FALSE, cache=TRUE}
# quantify multiome peaks in the scATAC-seq dataset
counts <- FeatureMatrix(
  fragments = atac.frags,
  features = granges(pbmc.multi),
  cells = atac.cells
)

# create object
atac.assay <- CreateChromatinAssay(
  counts = counts,
  min.features = 1000,
  fragments = atac.frags
)
pbmc.atac <- CreateSeuratObject(counts = atac.assay, assay = "peaks")
pbmc.atac <- subset(pbmc.atac, nCount_peaks > 2000 & nCount_peaks < 30000)

# compute LSI
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 10)
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- RunSVD(pbmc.atac)
```

Next we can merge the multiome and scATAC datasets together and observe that
there is a difference between them that appears to be due to the batch
(experiment and technology-specific variation).

```{r message=FALSE, warning=FALSE}
# first add dataset-identifying metadata
pbmc.atac$dataset <- "ATAC"
pbmc.multi$dataset <- "Multiome"

# merge
pbmc.combined <- merge(pbmc.atac, pbmc.multi)

# process the combined dataset
pbmc.combined <- FindTopFeatures(pbmc.combined, min.cutoff = 10)
pbmc.combined <- RunTFIDF(pbmc.combined)
pbmc.combined <- RunSVD(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, reduction = "lsi", dims = 2:30)
p1 <- DimPlot(pbmc.combined, group.by = "dataset")
```

## Integration

To find integration anchors between the two datasets, we need to project them into
a shared low-dimensional space. To do this, we'll use reciprocal LSI projection
(projecting each dataset into the others LSI space) by setting `reduction="rlsi"`.
For more information about the data integration methods in Seurat, see our recent
[paper](https://doi.org/10.1016/j.cell.2019.05.031)
and the [Seurat website](https://satijalab.org/seurat/).

Rather than integrating the normalized data matrix, as is typically done for 
scRNA-seq data, we'll integrate the low-dimensional cell embeddings (the LSI
coordinates) across the datasets using the `IntegrateEmbeddings()` function.
This is much better suited to scATAC-seq data,
as we typically have a very sparse matrix with a large number of features. Note 
that this requires that we first compute an uncorrected LSI embedding using the
merged dataset (as we did above).

```{r message=FALSE, warning=FALSE}
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(pbmc.multi, pbmc.atac),
  anchor.features = rownames(pbmc.multi),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = pbmc.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated, group.by = "dataset")
```

Finally, we can compare the results of the merged and integrated datasets, and
find that the integration has successfully removed the technology-specific variation
in the dataset while retaining the cell-type-specific (biological) variation.

```{r message=FALSE, warning=FALSE, fig.width=12}
(p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
```

Here we've demonstrated the integration method using two datasets, but the same
workflow can be applied to integrate any number of datasets.

## Reference mapping

In cases where we have a large, high-quality dataset, or a dataset containing unique
information not present in other datasets (cell type annotations or additional 
data modalities, for example), we often want to use that dataset as 
a reference and map queries onto it so that we can interpret these query datasets
in the context of the existing reference.

To demonstrate how to do this using single-cell chromatin reference and query
datasets, we'll treat the PBMC multiome dataset here as a reference and map the
scATAC-seq dataset to it using the `FindTransferAnchors()` and `MapQuery()`
functions from Seurat.

```{r message=FALSE, warning=FALSE}
# compute UMAP and store the UMAP model
pbmc.multi <- RunUMAP(pbmc.multi, reduction = "lsi", dims = 2:30, return.model = TRUE)

# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.multi,
  query = pbmc.atac,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

# map query onto the reference dataset
pbmc.atac <- MapQuery(
  anchorset = transfer.anchors,
  reference = pbmc.multi,
  query = pbmc.atac,
  refdata = pbmc.multi$predicted.id,
  reference.reduction = "lsi",
  new.reduction.name = "ref.lsi",
  reduction.model = 'umap'
)
```

<details>
  <summary>**What is `MapQuery()` doing?**</summary>

`MapQuery()` is a wrapper function that runs `TransferData()`, `IntegrateEmbeddings()`,
and `ProjectUMAP()` for a query dataset, and sets sensible default parameters based
on how the anchor object was generated. For finer control over the parameters used 
by each of these functions, you can pass parameters through `MapQuery()` to each function
using the `transferdata.args`, `integrateembeddings.args`, and `projectumap.args` arguments
for `MapQuery()`, or you can run each of the functions yourself. For example:

```{r, eval=FALSE}
pbmc.atac <- TransferData(
  anchorset = transfer.anchors, 
  reference = pbmc.multi,
  weight.reduction = "lsiproject",
  query = pbmc.atac,
  refdata = list(
    celltype = "predicted.id",
    predicted_RNA = "RNA")
)
pbmc.atac <- IntegrateEmbeddings(
  anchorset = transfer.anchors,
  reference = pbmc.multi,
  query = pbmc.atac, 
  reductions = "lsiproject",
  new.reduction.name = "ref.lsi"
)
pbmc.atac <- ProjectUMAP(
  query = pbmc.atac, 
  query.reduction = "ref.lsi",
  reference = pbmc.multi, 
  reference.reduction = "lsi",
  reduction.model = "umap"
)
```

</details>

By running `MapQuery()`, we have mapped the scATAC-seq dataset onto the the
multimodal reference, and enabled cell type labels to be transferred from reference
to query. We can visualize these reference mapping results and the cell type 
labels now associated with the scATAC-seq dataset:

```{r message=FALSE, warning=FALSE, fig.width=12, fig.height=6}
p1 <- DimPlot(pbmc.multi, reduction = "umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Reference")
p2 <- DimPlot(pbmc.atac, reduction = "ref.umap", group.by = "predicted.id", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("Query")

p1 | p2
```

For more information about multimodal reference mapping, see the [Seurat vignette](https://satijalab.org/seurat/articles/multimodal_reference_mapping.html).

## RNA imputation

Above we transferred categorical information (the cell labels) and mapped the 
query data onto an existing reference UMAP. We can also transfer continuous data
from the reference to the query in the same way. Here we demonstrate transferring
the gene expression values from the PBMC multiome dataset (that measured DNA
accessibility and gene expression in the same cells) to the PBMC scATAC-seq 
dataset (that measured DNA accessibility only). Note that we could also transfer
these values using the `MapQuery()` function call above by setting the `refdata`
parameter to a list of values.

```{r message=FALSE, warning=FALSE}
# predict gene expression values
rna <- TransferData(
  anchorset = transfer.anchors,
  refdata = GetAssayData(pbmc.multi, assay = "RNA", slot = "data"),
  weight.reduction = pbmc.atac[["lsi"]],
  dims = 2:30
)

# add predicted values as a new assay
pbmc.atac[["predicted"]] <- rna
```

We can look at some immune marker genes and see that the predicted expression
patterns match our expectation based on known expression patterns.

```{r message=FALSE, warning=FALSE, fig.width=12, fig.height=8}
DefaultAssay(pbmc.atac) <- "predicted"

FeaturePlot(
  object = pbmc.atac,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  reduction = "ref.umap",
  ncol = 3
)
```

```{r include=FALSE}
saveRDS(object = pbmc.atac, file = "../vignette_data/pbmc_atac_integration.rds")
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Data structures and object interaction"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

The Signac package is an extension of [Seurat](https://satijalab.org/seurat/)
designed for the analysis of genomic single-cell assays. This includes any assay
that generates signal mapped to genomic coordinates, such as scATAC-seq,
scCUT&Tag, scACT-seq, and other methods.

As the analysis of these single-cell chromatin datasets presents some unique 
challenges in comparison to the analysis of scRNA-seq data, we have created an
extended `Assay` class to store the additional information needed, including:

* Genomic ranges associated with the features (eg, peaks or genomic bins)
* Gene annotations
* Genome information 
* TF motifs
* Genome-wide signal in a disk-based format (fragment files)
* TF footprinting data
* Tn5 insertion bias data
* Linked genomic regions

A major advantage of the Signac design is its interoperability with existing
functions in the Seurat package, and other packages that are able to use the
Seurat object. This enables straightforward analysis of multimodal single-cell
data through the addition of different assays to the Seurat object.

Here we outline the design of each class defined in the Signac package, and
demonstrate methods that can be run on each class.

## The `ChromatinAssay` Class

The `ChromatinAssay` class extends the standard Seurat `Assay` class and adds
several additional slots for data useful for the analysis of single-cell
chromatin datasets. The class includes all the slots present in a standard 
Seurat [Assay](https://github.com/satijalab/seurat/wiki/Assay), 
with the following additional slots:

* `ranges`: A [`GRanges`](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class)
object containing the genomic coordinates of each feature in the `data` matrix.
* `motifs`: A `Motif` object
* `fragments`: A list of `Fragment` objects
* `seqinfo`: A [`Seqinfo`](https://www.rdocumentation.org/packages/GenomeInfoDb/versions/1.8.3/topics/Seqinfo-class)
object containing information about the genome that the data was mapped to
* `annotation`: A [`GRanges`](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class)
object containing gene annotations
* `bias`: A vector containing Tn5 integration bias information (the frequency of
Tn5 integration at different hexamers)
* `positionEnrichment`: A named list of matrices containing positional
enrichment scores for Tn5 integration (for example, enrichment at the TSS or at
different TF motifs)
* `links`: A [`GRanges`](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class)
object describing linked genomic positions, such as co-accessible sites or
enhancer-gene regulatory relationships.

```{r message=FALSE, warning=FALSE}
library(Seurat)
library(Signac)
```

### Constructing the `ChromatinAssay`

A `ChromatinAssay` object can be constructed using the
`CreateChromatinAssay()` function.

```{r}
# get some data to use in the following examples
counts <- GetAssayData(atac_small, slot = "counts")
```

```{r}
# create a standalone ChromatinAssay object
chromatinassay <- CreateChromatinAssay(counts = counts, genome = "hg19")
```

Here the `genome` parameter can be used to set the `seqinfo` slot. We can pass
the name of a genome present in UCSC (e.g., "hg19" or "mm10"), or we can pass
a `Seqinfo`-class object.

To create a Seurat object that contains a `ChromatinAssay` rather than a
standard `Assay`, we can initialize the object using the `ChromatinAssay` rather
than a count matrix. Note that this feature was added in Seurat 3.2.

```{r}
# create a Seurat object containing a ChromatinAssay
object <- CreateSeuratObject(counts = chromatinassay)
```

### Adding a `ChromatinAssay` to a `Seurat` object

To add a new `ChromatinAssay` object to an existing Seurat object, we can use
the standard assignment operation used for adding standard `Assay` objects and 
other data types to the Seurat object.

```{r}
# create a chromatin assay and add it to an existing Seurat object
object[["peaks"]] <- CreateChromatinAssay(counts = counts, genome = "hg19")
```

### Getting and setting `ChromatinAssay` data

We can get/set data for the `ChromatinAssay` in much the same way we do for a
standard `Assay` object: using the `GetAssayData` and `SetAssayData` functions
defined in `Seurat`. For example:

```{r}
## Getting

# access the data slot, found in standard Assays and ChromatinAssays
data <- GetAssayData(atac_small, slot = "data")

# access the bias slot, unique to the ChromatinAssay
bias <- GetAssayData(atac_small, slot = "bias")

## Setting

# set the data slot
atac_small <- SetAssayData(atac_small, slot = "data", new.data = data)

# set the bias slot
bias <- rep(1, 100)  # create a dummy bias vector
atac_small <- SetAssayData(atac_small, slot = "bias", new.data = bias)
```

We also have a variety of convenience functions defined for getting/setting 
data in specific slots. This includes the `Fragments()`, `Motifs()`, `Links()`,
and `Annotation()` functions. For example, to get or set gene annotation data we
can use the `Annotation()` getter and `Annotation<-` setter functions:

```{r message=FALSE, warning=FALSE}
# first get some gene annotations for hg19
library(EnsDb.Hsapiens.v75)

# convert EnsDb to GRanges
gene.ranges <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# convert to UCSC style
seqlevelsStyle(gene.ranges) <- "UCSC"

# set gene annotations
Annotation(atac_small) <- gene.ranges

# get gene annotation information
Annotation(atac_small)
```

The `Fragments()`, `Motifs()`, and `Links()` functions are demonstrated in other sections
below.

### Other `ChromatinAssay` methods

As the `ChromatinAssay` object uses Bioconductor objects like [`GRanges`](https://www.rdocumentation.org/packages/GenomicRanges/versions/1.24.1/topics/GRanges-class)
and [`Seqinfo`](https://www.rdocumentation.org/packages/GenomeInfoDb/versions/1.8.3/topics/Seqinfo-class)
, we can also call standard Bioconductor functions defined in the 
`IRanges`, `GenomicRanges`, and `GenomeInfoDb` packages on the `ChromatinAssay`
object (or a Seurat object with a `ChromatinAssay` as the default assay).

The following methods use the genomic ranges stored in a `ChromatinAssay` object.

```{r message=FALSE, warning=FALSE}
# extract the genomic ranges associated with each feature in the data matrix
granges(atac_small)

# find the nearest range
nearest(atac_small, subject = Annotation(atac_small))

# distance to the nearest range
distanceToNearest(atac_small, subject = Annotation(atac_small))

# find overlaps with another set of genomic ranges
findOverlaps(atac_small, subject = Annotation(atac_small))
```

Many other methods are defined, see the documentation for `nearest-methods`,
`findOverlaps-methods`, `inter-range-methods`, and `coverage` in Signac for a
full list.

The following methods use the `seqinfo` data stored in a `ChromatinAssay` object.

```{r}
# get the full seqinfo information
seqinfo(atac_small)

# get the genome information
genome(atac_small)

# find length of each chromosome
seqlengths(atac_small)

# find name of each chromosome
seqnames(atac_small)

# assign a new genome
genome(atac_small) <- "hg19"
```

Again, several other methods are available that are not listed here. See the 
documentation for `seqinfo-methods` in Signac for a full list.

For a full list of methods for the `ChromatinAssay` class run:

```{r}
methods(class = 'ChromatinAssay')
```

### Subsetting a `ChromatinAssay`

We can use the standard `subset()` function or the `[` operator to subset Seurat
object containing `ChromatinAssay`s. This works the same way as for standard
`Assay` objects.

```{r}
# subset using the subset() function
# this is meant for interactive use
subset.obj <- subset(atac_small, subset = nCount_peaks > 100)

# subset using the [ extract operator
# this can be used programmatically
subset.obj <- atac_small[, atac_small$nCount_peaks > 100]
```

### Converting between `Assay` and `ChromatinAssay`

To convert from a `ChromatinAssay` to a standard `Assay` use the `as()` function

```{r}
# convert a ChromatinAssay to an Assay
assay <- as(object = atac_small[["peaks"]], Class = "Assay")
assay
```

To convert from a standard `Assay` to a `ChromatinAssay` we use the 
`as.ChromatinAssay()` function. This takes a standard assay object, as well as
information to fill the additional slots in the `ChromatinAssay` class.

```{r}
# convert an Assay to a ChromatinAssay
chromatinassay <- as.ChromatinAssay(assay, seqinfo = "hg19")
chromatinassay
```

## The `Fragment` Class

The `Fragment` class is designed for storing and interacting with a
[fragment file](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments)
commonly used for single-cell chromatin data. It contains the path to an indexed 
fragment file on disk, a MD5 hash for the fragment file and the fragment file
index, and a vector of cell names contained in the fragment file. Importantly,
this is a named vector where the elements of the vector are the cell names as 
they appear in the fragment file, and the name of each element is the cell
name as it appears in the `ChromatinAssay` object storing the `Fragment`
object. This allows a mapping of cell names on disk to cell names in R, and 
avoids the need to alter fragment files on disk. This path can also be a remote
file accessible by `http` or `ftp`.

### Constructing the `Fragment` class

A `Fragment` object can be constructed using the `CreateFragmentObject()`
function. 

```{r}
frag.path <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = frag.path,
  cells = colnames(atac_small), 
  validate.fragments = TRUE
)
```

The `validate.fragments` parameter controls whether the file is inspected to
check whether the expected cell names are present. This can help avoid assigning
the wrong fragment file to the object. If you're sure that the file is correct,
you can set this value to `FALSE` to skip this step and save some time. This
check is typically only run once when the `Fragment` object is created, and is 
not normally run on existing `Fragment` files.

### Inspecting the fragment file

To extract the first few lines of a fragment file on-disk, we can use the
`head()` method defined for `Fragment` objects. This is useful for quickly
checking the chromosome naming style in our fragment file, or checking how the 
cell barcodes are named:

```{r}
head(fragments)
```

### Adding a `Fragment` object to the `ChromatinAssay`

A `ChromatinAssay` object can contain a list of `Fragment` objects. This avoids
the need to merge fragment files on disk and simplifies processes of merging
or integrating different Seurat objects containing `ChromatinAssay`s. To add 
a new `Fragment` object to a `ChromatinAssay`, or a Seurat object containing a 
`ChromatinAssay`, we can use the `Fragments<-` assignment function. This will
do a few things:

1. Re-compute the MD5 hash for the fragment file and index and verify that it
matches the hash computed when the `Fragment` object was created.
2. Check that none of the cells contained in the `Fragment` object being added
are already contained in another `Fragment` object stored in the
`ChromatinAssay`. All fragments from a cell must be present in only one fragment
file.
3. Append the `Fragment` object to the list of `Fragment` objects stored in the
`ChromatinAssay`.

```{r}
Fragments(atac_small) <- fragments
```

The `show()` method for `Fragment`-class objects prints the number of cells that
the `Fragment` object contains data for.

```{r}
fragments
```

Alternatively, we can initialize the `ChromatinAssay` with a `Fragment` object
in a couple of ways. We can either pass a vector of `Fragment` objects to the
`fragments` parameter in `CreateChromatinAssay()`, or pass the path to a
single fragment file. If we pass the path to a fragment file we assume that the
file contains fragments for all cells in the `ChromatinAssay` and that the cell
names are the same in the fragment file on disk and in the `ChromatinAssay`.
For example:

```{r}
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  genome = "hg19",
  fragments = frag.path
)
object <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks"
)
```

This will create a Seurat object containing a `ChromatinAssay`, with a single
`Fragment` object.

### Removing a `Fragment` object from the `ChromatinAssay`

All the `Fragment` objects associated with a `ChromatinAssay` can be removed
by assigning `NULL` using the `Fragment<-` assignment function. For example:

```{r}
Fragments(chrom_assay) <- NULL
Fragments(chrom_assay)
```

To remove a subset of `Fragment` object from the list of `Fragment` objects
stored in the `ChromatinAssay`, you will need to extract the list of `Fragment`
objects using the `Fragments()` function, subset the list of objects, then
assign the subsetted list to the assay using the `Seurat::SetAssayData()`
function. For example:

```{r}
chrom_assay <- SetAssayData(chrom_assay, slot = "fragments", new.data = fragments)
Fragments(chrom_assay)
```

### Changing the fragment file path in an existing `Fragment` object

The path to the fragment file can be updated using the `UpdatePath()` function.
This can be useful if you move the fragment file to a new directory, or if you
copy a stored Seurat object containing a `ChromatinAssay` to a different server.

```{r message=FALSE}
fragments <- UpdatePath(fragments, new.path = frag.path)
```

To change the path to fragment files in an object, you will need to remove
the fragment objects, update the paths, and then add the fragment objects
back to the object. For example:

```{r}
frags <- Fragments(object)  # get list of fragment objects
Fragments(object) <- NULL  # remove fragment information from assay

# create a vector with all the new paths, in the correct order for your list of fragment objects
# In this case we only have 1
new.paths <- list(frag.path)
for (i in seq_along(frags)) {
  frags[[i]] <- UpdatePath(frags[[i]], new.path = new.paths[[i]]) # update path
}

Fragments(object) <- frags # assign updated list back to the object
Fragments(object)
```

### Using remote fragment files

Fragment files hosted on remote servers accessible via http or ftp can also be
added to the `ChromatinAssay` in the same way as for locally-hosted fragment 
files. This can enable the exploration of large single-cell datasets without
the need for downloading large files. For example, we can create a Fragment
object using a file hosted on the 10x Genomics website:

```{r}
fragments <- CreateFragmentObject(
  path = "http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz"
)
fragments
```

When files are hosted remotely, the checks described in the section above (MD5
hash and expected cells) are not performed.

### Getting and setting `Fragment` data

To access the cell names stored in a `Fragment` object, we can use the `Cells()`
function. **Importantly**, this returns the cell names as they appear in the 
`ChromatinAssay`, rather than as they appear in the fragment file itself.

```{r}
fragments <- CreateFragmentObject(
  path = frag.path,
  cells = colnames(atac_small), 
  validate.fragments = TRUE
)

cells <- Cells(fragments)
head(cells)
```

Similarly, we can set the cell name information in a `Fragment` object using the
`Cells<-` assignment function. This will set the named vector of cells stored in
the `Fragment` object. Here we must supply a named vector.

```{r}
names(cells) <- cells
Cells(fragments) <- cells
```

To extract any of the data stored in a `Fragment` object we can also use the
`GetFragmentData()` function. For example, we can find the path to the fragment
file on disk:

```{r}
GetFragmentData(object = fragments, slot = "path")
```

For a full list of methods for the `Fragment` class run:

```{r}
methods(class = 'Fragment')
```

## The `Motif` Class

The `Motif` class stores information needed for DNA sequence motif analysis, and
has the following slots:

* `data`: a sparse feature by motif matrix, where entries are 1 if the feature
contains the motif, and 0 otherwise
* `pwm`: A named list of position weight or position frequency matrices
* `motif.names`: a list of motif IDs and their common names
* `positions`: A `GRangesList` object containing the exact positions of each
motif
* `meta.data`: Additional information about the motifs 

Many of these slots are optional and do not need to be filled, but are only
required when running certain functions. For example, the `positions` slot
will be needed if running TF footprinting.

### Constructing the `Motif` class

A `Motif` object can be constructed using the `CreateMotifObject()` function.
Much of the data needed for constructing a `Motif` object can be generated using
functions from the [TFBSTools](https://www.bioconductor.org/packages/release/bioc/html/TFBSTools.html)
and [motifmatchr](https://www.bioconductor.org/packages/release/bioc/html/motifmatchr.html)
packages. Position frequency matrices for motifs can be loaded using the
[JASPAR](http://jaspar.genereg.net/) packages on
[Bioconductor](https://bioconductor.org/packages/release/data/annotation/html/JASPAR2020.html)
or the [chromVARmotifs](https://github.com/GreenleafLab/chromVARmotifs) package.
For example:

```{r message=FALSE, warning=FALSE}
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606) # 9606 is the species code for human
)

# Scan the DNA sequence of each peak for the presence of each motif
motif.matrix <- CreateMotifMatrix(
  features = granges(atac_small),
  pwm = pfm,
  genome = 'hg19'
)

# Create a new Mofif object to store the results
motif <- CreateMotifObject(
  data = motif.matrix,
  pwm = pfm
)
```

The `show()` method for the `Motif` class prints the total number of motifs and
regions included in the object:

```{r}
motif
```

### Adding a `Motif` object to the `ChromatinAssay`

We can add a `Motif` object to the `ChromatinAssay`, or a Seurat object
containing a `ChromatinAssay` using the `Motifs<-` assignment operator.

```{r}
Motifs(atac_small) <- motif
```

### Getting and setting `Motif` data

Data stored in a `Motif` object can be accessed using the `GetMotifData()` and 
`SetMotifData()` functions.

```{r}
# extract data from the Motif object
pfm <- GetMotifData(object = motif, slot = "pwm")

# set data in the Motif object
motif <- SetMotifData(object = motif, slot = "pwm", new.data = pfm)
```

We can access the set of motifs and set of features used in the `Motif` object
using the `colnames()` and `rownames()` functions:

```{r}
# look at the motifs included in the Motif object
head(colnames(motif))
```

```{r}
# look at the features included in the Motif object
head(rownames(motif))
```

To quickly convert between motif IDs (like `MA0497.1`) and motif common names 
(like MEF2C), we can use the `ConvertMotifID()` function. For example:

```{r}
# convert ID to common name
ids <- c("MA0025.1","MA0030.1","MA0031.1","MA0051.1","MA0056.1","MA0057.1")
names <- ConvertMotifID(object = motif, id = ids)
names
```

```{r}
# convert names to IDs
ConvertMotifID(object = motif, name = names)
```

For a full list of methods for the `Motif` class run:

```{r warning=FALSE, message=FALSE}
methods(class = 'Motif')
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>

---
title: "Parallel and distributed processing"
output: html_document
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
---

Parallel computing is supported in Signac through the [future](https://cran.r-project.org/package=future) package, making it easy to specify different parallelization options.
Here we demonstrate parallelization of the `FeatureMatrix` function and show some benchmark results to get a sense for the amount of speedup you might expect.

The [Seurat](https://satijalab.org/seurat/) package also uses `future` for 
parallelization, and you can see the Seurat [vignette](https://satijalab.org/seurat/v3.1/future_vignette.html) for more information.

## How to enable parallelization in Signac

Parallelization can be enabled simply by importing the `future` package and setting the `plan`.

```{r}
library(future)
plan()
```

By default the plan is set to sequential processing (no parallelization). We can change this
to `multicore`, `multiprocessor`, or `multisession` to get asynchronous processing, and set the 
number of workers to change the number of cores used.

```{r, warning=FALSE, message=FALSE}
plan("multicore", workers = 10)
plan()
```

You might also need to increase the maximum memory usage:

```{r}
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 50 Gb RAM
```

Note that as of `future` version [1.14.0](https://cran.r-project.org/web/packages/future/NEWS), forked processing is disabled when running in RStudio. To enable parallel computing in RStudio, you will need to select the "multisession" option.

## Benchmarking

Here we demonstrate the runtime of `FeatureMatrix` run on 90,686 peaks for 10,247 human PBMCs under different parallelization options:

<details>
  <summary>**View benchmarking code**</summary>
  
The following code was run on Ubuntu 20.04 LTS with 56 Intel Xeon Platinum 8280L CPU @ 2.70GHz

```{bash, eval = FALSE}
# download data
wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_peaks.bed
```
  
```{r, eval=FALSE}
library(Signac)

# load data
fragments <- "../vignette_data/atac_pbmc_10k_nextgem_fragments.tsv.gz"
peaks.10k <- read.table(
  file = "../vignette_data/atac_pbmc_10k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks.10k)
cells <- readLines("../vignette_data/pbmc10k/cells.txt")

fragments <- CreateFragmentObject(path = fragments, cells = cells, validate.fragments = FALSE)

# set number of replicates
nrep <- 5
results <- data.frame()
process_n <- 2000

# run sequentially
timing.sequential <- c()
for (i in seq_len(nrep)) {
  start <- Sys.time()
  fmat <- FeatureMatrix(fragments = fragments, features = peaks, cells = cells, process_n = process_n)
  timing.sequential <- c(timing.sequential, as.numeric(Sys.time() - start, units = "secs"))
}
res <- data.frame(
  "setting" = rep("Sequential", nrep),
  "cores" = rep(1, nrep),
  "replicate" = seq_len(nrep),
  "time" = timing.sequential
)
results <- rbind(results, res)

# 4 core
library(future)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)

timing.4core <- c()
for (i in seq_len(nrep)) {
  start <- Sys.time()
  fmat <- FeatureMatrix(fragments = fragments, features = peaks, cells = cells, process_n = process_n)
  timing.4core <- c(timing.4core, as.numeric(Sys.time() - start, units = "secs"))
}
res <- data.frame(
  "setting" = rep("Parallel", nrep),
  "cores" = rep(4, nrep),
  "replicate" = seq_len(nrep),
  "time" = timing.4core
)
results <- rbind(results, res)

# 10 core
plan("multicore", workers = 10)

timing.10core <- c()
for (i in seq_len(nrep)) {
  start <- Sys.time()
  fmat <- FeatureMatrix(fragments = fragments, features = peaks, cells = cells, process_n = process_n)
  timing.10core <- c(timing.10core, as.numeric(Sys.time() - start, units = "secs"))
}
res <- data.frame(
  "setting" = rep("Parallel", nrep),
  "cores" = rep(10, nrep),
  "replicate" = seq_len(nrep),
  "time" = timing.10core
)
results <- rbind(results, res)

# save results
write.table(
  x = results,
  file = paste0("../vignette_data/pbmc10k/timings_", Sys.Date(), ".tsv"),
  quote = FALSE,
  row.names = FALSE
)
```

</details>

```{r, message=FALSE, warning=FALSE, echo=FALSE}
library(ggplot2)
results <- read.table("../vignette_data/pbmc10k/timings_2021-09-20.tsv", header = TRUE)
results$cores <- factor(results$cores, levels = c("1", "4", "10"))

p <- ggplot(results, aes(x = cores, y = time/60, color = cores)) +
  geom_jitter(width = 1/10) +
  ylim(c(0, 13)) +
  ylab("Time (min)") +
  xlab("Number of cores") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("FeatureMatrix runtime")
p
```

<details>
  <summary>**Session Info**</summary>

```{r}
sessionInfo()
```

</details>
---
title: "Archive"
output: html_document
---

To install past versions of Signac see the
[installation instructions](install.html#old-versions)

Past versions of the Signac documentation can be found here:

## [0.2](../0.2/index.html)

## [1.2.0](../1.2.0/index.html)
---
title: "Calling peaks"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

In this tutorial, we demonstrate how to call peaks on a single-cell ATAC-seq
dataset using [MACS2](https://github.com/macs3-project/MACS).

To use the peak calling functionality in Signac you will first need to install
MACS2. This can be done using [pip](https://pypi.org/project/MACS2/) or
[conda](https://anaconda.org/bioconda/macs2), or by building the package
from [source](https://github.com/macs3-project/MACS).

In this demonstration we use scATAC-seq data for human PBMCs. See our
[vignette](pbmc_vignette.Rmd) for the code used to generate this object,
and links to the raw data. First, load the required packages and the
pre-computed Seurat object:

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)

pbmc <- readRDS("../vignette_data/pbmc.rds")
DimPlot(pbmc)
```

Peak calling can be performed using the `CallPeaks()` function, and can either
be done separately for different groups of cells, or performed using data from
all the cells. To call peaks on each annotated cell type, we can use the 
`group.by` argument:

```{r message=FALSE, warning=FALSE, cache=TRUE}
peaks <- CallPeaks(
  object = pbmc,
  group.by = "predicted.id",
  macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2"
)
```

The results are returned as a `GRanges` object, with an additional metadata column
listing the cell types that each peak was identified in:

```{r echo=FALSE}
knitr::kable(head(as.data.frame(peaks)))
```

To quantify counts in each peak, you can use the `FeatureMatrix()` function.

We can visualize the cell-type-specific MACS2 peak calls alongside the 10x Cellranger
peak calls (currently being used in the `pbmc` object) with the `CoveragePlot()`
function. Here the Cellranger peaks are shown in grey and the MACS2 peaks in red:

```{r fig.height=6, fig.width=10}
CoveragePlot(
  object = pbmc,
  region = "CD8A",
  ranges = peaks,
  ranges.title = "MACS2"
)
```
<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Finding co-accessible networks with Cicero"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

In this vignette we will demonstrate how to find cis-co-accessible networks with 
Cicero using single-cell ATAC-seq data. Please see the
Cicero [website](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/)
for information about Cicero.

To facilitate conversion between the Seurat (used by Signac) and CellDataSet
(used by Cicero) formats, we will use a conversion function in the 
[SeuratWrappers](https://github.com/satijalab/seurat-wrappers) package available
on GitHub.

## Data loading

We will use a single-cell ATAC-seq dataset containing human CD34+ hematopoietic
stem and progenitor cells published by Satpathy and Granja et al. (2019, Nature
Biotechnology). The processed datasets are available on NCBI GEO here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785

This is the same dataset we used in the [trajectory](monocle.html)
vignette, and we'll start by loading the dataset that was created in that
vignette. See the [trajectory](monocle.html) vignette for the code used to
create the object from raw data.

First we will load their dataset and perform some standard preprocessing using
Signac.

```{r include=FALSE}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

if (!requireNamespace("SeuratWrappers", quietly = TRUE))
    remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("leidenbase", quietly = TRUE))
    remotes::install_github("cole-trapnell-lab/leidenbase")

if (!requireNamespace("monocle3", quietly = TRUE))
    remotes::install_github("cole-trapnell-lab/monocle3")

if (!requireNamespace("cicero", quietly = TRUE))
    remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
```

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
```

```{r}
# load the object created in the Monocle 3 vignette
bone <- readRDS("../vignette_data/cd34.rds")
```

## Create the Cicero object

We can find cis-co-accessible networks (CCANs) using [Cicero](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/).

The Cicero developers have developed a separate branch of the package that works
with a Monocle 3 `CellDataSet` object. We will first make sure this branch is
installed, then convert our `Seurat` object for the whole bone marrow dataset to
`CellDataSet` format.

```{r, eval=FALSE}
# Install Cicero
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
```

```{r message=FALSE, warning=FALSE}
library(cicero)
```

```{r message=FALSE, warning=FALSE}
# convert to CellDataSet format and make the cicero object
bone.cds <- as.cell_data_set(x = bone)
bone.cicero <- make_cicero_cds(bone.cds, reduced_coordinates = reducedDims(bone.cds)$UMAP)
```

## Find Cicero connections

We'll demonstrate running Cicero here using just one chromosome to save some time,
but the same workflow can be applied to find CCANs for the whole genome.

Here we demonstrate the most basic workflow for running Cicero. This workflow
can be broken down into several steps, each with parameters that can be changed
from their defaults to fine-tune the Cicero algorithm depending on your data.
We highly recommend that you explore the [Cicero website](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/),
[paper](https://doi.org/10.1016/j.molcel.2018.06.044), and documentation for
more information.

```{r, cache=TRUE}
# get the chromosome sizes from the Seurat object
genome <- seqlengths(bone)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(bone.cicero, genomic_coords = genome.df, sample_num = 100)
```

```{r}
head(conns)
```

## Find cis-co-accessible networks (CCANs)

Now that we've found pairwise co-accessibility scores for each peak, we can now
group these pairwise connections into larger co-accessible networks using the
`generate_ccans()` function from Cicero.

```{r}
ccans <- generate_ccans(conns)
```

```{r}
head(ccans)
```

## Add links to a Seurat object

We can add the co-accessible links found by Cicero to the `ChromatinAssay`
object in Seurat. Using the `ConnectionsToLinks()` function in Signac we can
convert the outputs of Cicero to the format needed to store in the `links` slot
in the `ChromatinAssay`, and add this to the object using the `Links<-`
assignment function.

```{r}
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(bone) <- links
```

We can now visualize these links along with DNA accessibility information by
running `CoveragePlot()` for a region:

```{r fig.height=8, fig.width=8, message=FALSE, warning=FALSE}
CoveragePlot(bone, region = "chr1-40189344-40252549")
```

```{r include=FALSE}
saveRDS(object = bone, file = "../vignette_data/cd34.rds")
```

## Acknowledgements

Thanks to the developers of Cicero, especially Cole Trapnell, Hannah Pliner, and
members of the [Trapnell lab](https://cole-trapnell-lab.github.io/). If you use
Cicero please cite the [Cicero paper](https://doi.org/10.1016/j.molcel.2018.06.044).

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Building trajectories with Monocle 3"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

In this vignette we will demonstrate how to construct cell trajectories with 
Monocle 3 using single-cell ATAC-seq data. Please see the
Monocle 3 [website](https://cole-trapnell-lab.github.io/monocle3/) for 
information about installing Monocle 3.

To facilitate conversion between the Seurat (used by Signac) and CellDataSet
(used by Monocle 3) formats, we will use a conversion function in the 
[SeuratWrappers](https://github.com/satijalab/seurat-wrappers) package available
on GitHub.

## Data loading

We will use a single-cell ATAC-seq dataset containing human CD34+ hematopoietic
stem and progenitor cells published by
[Satpathy and Granja et al. (2019, Nature Biotechnology)](https://doi.org/10.1038/s41587-019-0206-z).

The processed dataset is available on NCBI GEO here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129785

Note that the fragment file is present inside the `GSE129785_RAW.tar` archive,
and the index for the fragment file is not supplied. You can index the file
yourself using [tabix](https://www.htslib.org/doc/tabix.html), for example:
`tabix -p bed <fragment_file>`.

First we will load the dataset and perform some standard preprocessing using
Signac.

```{r include=FALSE}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

if (!requireNamespace("SeuratWrappers", quietly = TRUE))
    remotes::install_github('satijalab/seurat-wrappers')

if (!requireNamespace("leidenbase", quietly = TRUE))
    remotes::install_github("cole-trapnell-lab/leidenbase")

if (!requireNamespace("monocle3", quietly = TRUE))
    remotes::install_github("cole-trapnell-lab/monocle3")

if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE))
    BiocManager::install("EnsDb.Hsapiens.v75")
```

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)
```

```{r}
filepath <- "../vignette_data/GSE129785/GSE129785_scATAC-Hematopoiesis-CD34"

peaks <- read.table(paste0(filepath, ".peaks.txt.gz"), header = TRUE)
cells <- read.table(paste0(filepath, ".cell_barcodes.txt.gz"), header = TRUE, stringsAsFactors = FALSE)
rownames(cells) <- make.unique(cells$Barcodes)

mtx <- readMM(file = paste0(filepath, ".mtx"))
mtx <- as(object = mtx, Class = "dgCMatrix")
colnames(mtx) <- rownames(cells)
rownames(mtx) <- peaks$Feature
```

```{r message=FALSE, warning=FALSE}
bone_assay <- CreateChromatinAssay(
  counts = mtx,
  min.cells = 5,
  fragments = "../vignette_data/GSE129785/GSM3722029_CD34_Progenitors_Rep1_fragments.tsv.gz",
  sep = c("_", "_"),
  genome = "hg19"
)
bone <- CreateSeuratObject(
  counts = bone_assay,
  meta.data = cells,
  assay = "ATAC"
)

# The dataset contains multiple cell types
# We can subset to include just one replicate of CD34+ progenitor cells
bone <- bone[, bone$Group_Barcode == "CD34_Progenitors_Rep1"]

# add cell type annotations from the original paper
cluster_names <- c("HSC",	"MEP",	"CMP-BMP",	"LMPP",	"CLP",	"Pro-B",	"Pre-B",	"GMP",
                  "MDP",	"pDC",	"cDC",	"Monocyte-1",	"Monocyte-2",	"Naive-B",	"Memory-B",
                  "Plasma-cell",	"Basophil",	"Immature-NK",	"Mature-NK1",	"Mature-NK2",	"Naive-CD4-T1",
                  "Naive-CD4-T2",	"Naive-Treg",	"Memory-CD4-T",	"Treg",	"Naive-CD8-T1",	"Naive-CD8-T2",
                  "Naive-CD8-T3",	"Central-memory-CD8-T",	"Effector-memory-CD8-T",	"Gamma delta T")
num.labels <- length(cluster_names)
names(cluster_names) <- paste0( rep("Cluster", num.labels), seq(num.labels) )
bone$celltype <- cluster_names[as.character(bone$Clusters)]

bone[["ATAC"]]
```

Next we can add gene annotations for the hg19 genome to the object. This will 
be useful for computing quality control metrics (TSS enrichment score) and 
plotting.

```{r message=FALSE, warning=FALSE}
library(EnsDb.Hsapiens.v75)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(bone) <- annotations
```

## Quality control

We'll compute TSS enrichment, nucleosome signal score, and the percentage of 
counts in genomic blacklist regions for each cell, and use these metrics to 
help remove low quality cells from the datasets.

```{r message=FALSE, warning=FALSE}
bone <- TSSEnrichment(bone)
bone <- NucleosomeSignal(bone)
bone$blacklist_fraction <- FractionCountsInRegion(bone, regions = blacklist_hg19)
```

```{r fig.height=6, fig.width=10, message=FALSE, warning=FALSE}
VlnPlot(
  object = bone,
  features = c("nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 4
)
```

```{r}
bone <- bone[, (bone$nCount_ATAC < 50000) &
               (bone$TSS.enrichment > 2) & 
               (bone$nucleosome_signal < 5)]
```

## Dataset preprocessing

Next we can run a standard scATAC-seq analysis pipeline using Signac to perform
dimension reduction, clustering, and cell type annotation.

```{r message=FALSE, warning=FALSE}
bone <- FindTopFeatures(bone, min.cells = 10)
bone <- RunTFIDF(bone)
bone <- RunSVD(bone, n = 100)
DepthCor(bone)
```

```{r message=FALSE, warning=FALSE}
bone <- RunUMAP(
  bone,
  reduction = "lsi",
  dims = 2:50,
  reduction.name = "UMAP"
)
```

```{r message=FALSE, warning=FALSE}
bone <- FindNeighbors(bone, dims = 2:50, reduction = "lsi")
bone <- FindClusters(bone, resolution = 0.8, algorithm = 3)
```

```{r}
DimPlot(bone, label = TRUE) + NoLegend()
```

Assign each cluster to the most common cell type based on the original
annotations from the paper.

```{r}
for(i in levels(bone)) {
  cells_to_reid <- WhichCells(bone, idents = i)
  newid <- names(sort(table(bone$celltype[cells_to_reid]),decreasing=TRUE))[1]
  Idents(bone, cells = cells_to_reid) <- newid
}
bone$assigned_celltype <- Idents(bone)
```

```{r}
DimPlot(bone, label = TRUE)
```

Next we can subset the different lineages and create a trajectory for each
lineage. Another way to build the trajectories is to use the whole dataset and 
build separate pseudotime trajectories for the different cell partitions found
by Monocle 3.

```{r}
DefaultAssay(bone) <- "ATAC"

erythroid <- bone[,  bone$assigned_celltype %in% c("HSC", "MEP", "CMP-BMP")]
lymphoid <- bone[, bone$assigned_celltype %in% c("HSC", "LMPP", "GMP", "CLP", "Pro-B", "pDC", "MDP", "GMP")]
```

## Building trajectories with Monocle 3

We can convert the Seurat object to a CellDataSet object using the
`as.cell_data_set()` function from [SeuratWrappers](https://github.com/satijalab/seurat-wrappers)
and build the trajectories using Monocle 3. We'll do this separately for 
erythroid and lymphoid lineages, but you could explore other strategies building
a trajectory for all lineages together.

```{r message=FALSE, warning=FALSE, results='hide'}
erythroid.cds <- as.cell_data_set(erythroid)
erythroid.cds <- cluster_cells(cds = erythroid.cds, reduction_method = "UMAP")
erythroid.cds <- learn_graph(erythroid.cds, use_partition = TRUE)

lymphoid.cds <- as.cell_data_set(lymphoid)
lymphoid.cds <- cluster_cells(cds = lymphoid.cds, reduction_method = "UMAP")
lymphoid.cds <- learn_graph(lymphoid.cds, use_partition = TRUE)
```

To compute pseudotime estimates for each trajectory we need to decide what the
start of each trajectory is. In our case, we know that the hematopoietic stem
cells are the progenitors of other cell types in the trajectory, so we can set
these cells as the root of the trajectory. Monocle 3 includes an interactive 
function to select cells as the root nodes in the graph. This function will be 
launched if calling `order_cells()` without specifying the `root_cells` parameter.
Here we've pre-selected some cells as the root, and saved these to a file for 
reproducibility. This file can be downloaded [here](https://www.dropbox.com/s/w5jbokcj9u6iq04/hsc_cells.txt).

```{r}
# load the pre-selected HSCs
hsc <- readLines("../vignette_data/hsc_cells.txt")
```

```{r message=FALSE, warning=FALSE}
# order cells
erythroid.cds <- order_cells(erythroid.cds, reduction_method = "UMAP", root_cells = hsc)
lymphoid.cds <- order_cells(lymphoid.cds, reduction_method = "UMAP", root_cells = hsc)

# plot trajectories colored by pseudotime
plot_cells(
  cds = erythroid.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

plot_cells(
  cds = lymphoid.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
```

Extract the pseudotime values and add to the Seurat object

```{r}
bone <- AddMetaData(
  object = bone,
  metadata = erythroid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Erythroid"
)

bone <- AddMetaData(
  object = bone,
  metadata = lymphoid.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Lymphoid"
)
```

```{r fig.height=4, fig.width=8, message=FALSE, warning=FALSE}
FeaturePlot(bone, c("Erythroid", "Lymphoid"), pt.size = 0.1) & scale_color_viridis_c()
```

```{r, include=FALSE}
saveRDS(object = bone, file = "../vignette_data/cd34.rds")
```

## Acknowledgements

Thanks to the developers of Monocle 3, especially Cole Trapnell, Hannah Pliner,
and members of the [Trapnell lab](https://cole-trapnell-lab.github.io/). If you
use Monocle please cite the
[Monocle papers](https://cole-trapnell-lab.github.io/monocle3/docs/citations/).

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Frequently asked questions"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## What is Signac?

Signac is an extension of Seurat for the analysis of single-cell chromatin data
(DNA-based single-cell assays). We have extended the Seurat object to include
information about the genome sequence and genomic coordinates of sequenced
fragments per cell, and include functions needed for the analysis of single-cell
chromatin data.

## How do I interact with the object?

Signac uses the Seurat object structure, and so all the Seurat commands can be
used when analysing data with Signac. See the Data Structures and Object
Interaction [vignette](data_structures.html) for an explanation of the classes defined in Signac
and how to use them. See the Seurat documentation for more information about the
Seurat object: https://satijalab.org/seurat/

## How do I merge objects with Signac?

See the [merge](merging.html) and [integration](integration.html) vignettes for
information on combining multiple single-cell chromatin datasets.

## How can I create a fragment file for my dataset?

The fragment file is provided in the output of popular single-cell data
data processing tools such as [chromap](https://github.com/haowenz/chromap), [cellranger-atac](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac),
and [cellranger-arc](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc).

If you are using another method that does not provide a fragment file as output,
you can use the [sinto](https://github.com/timoast/sinto) package to generate
a fragment file from the BAM file. See here for more information on using Sinto
to generate a fragment file:
https://timoast.github.io/sinto/basic_usage.html#create-scatac-seq-fragments-file

## How should I decide on the number of dimensions to use?

Choosing the dimensionality is a general problem in single-cell analysis for
which there is no simple solution. There has been discussion about this for
scRNA-seq, and you can read our recommendations for scRNA-seq in the Seurat
vignettes: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html (see
"Determine the ‘dimensionality’ of the dataset").

Here are some general tips/suggestions that might help guide you in the choice
for number of dimensions:

* the number of dimensions needed will generally scale with the size and
complexity of the dataset
* you can try varying the number of dimensions used and observing how the
resulting clusters or UMAP changes
* it is usually better to choose values that are higher rather than too low
* having a good understanding of the biology will help a lot in knowing whether
the clusters make sense, or if the dimensionality might be too high/low

## An annotation or genome sequence for my organism is not available on Bioconductor, what do I do?

If you are studying an organism that does not have a `BSgenome` genome package
or `EnsDB` annotation package available on BioConductor, you can still use your
own GTF file or FASTA files with Signac.

To create your own `BSgenome` data package, see this
[vignette](https://www.bioconductor.org/packages/devel/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf).

To use a GTF file, you can import it using `rtracklayer`, for example:

```{r, eval=FALSE}
gtf <- rtracklayer::import('genes.gtf')
gene.coords <- gtf[gtf$type == 'gene']
seqlevelsStyle(gene.coords) <- 'UCSC'
gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
```

## How should I cite Signac?

If you use Signac, please cite [Stuart et al., 2021](https://doi.org/10.1038/s41592-021-01282-5):

```{r, eval=FALSE}
@ARTICLE{signac,
  title     = "Single-cell chromatin state analysis with Signac",
  author    = "Stuart, Tim and Srivastava, Avi and Madad, Shaista and Lareau,
               Caleb A and Satija, Rahul",
  journal   = "Nat. Methods",
  publisher = "Nature Publishing Group",
  pages     = "1--9",
  month     =  nov,
  year      =  2021,
  url       = "https://www.nature.com/articles/s41592-021-01282-5",
  language  = "en"
}
```

Signac is an extension of [Seurat](https://satijalab.org/seurat/), and uses the
Seurat object structure, so you should consider citing the
[Seurat paper](https://doi.org/10.1016/j.cell.2019.05.031) if you have used Signac.

---
title: "Vignettes overview"
output: html_document
---

# Guided analyses

The following guided analyses demonstrate a standard end-to-end analysis 
pipeline for different types of single-cell chromatin data.

## [Human peripheral blood mononuclear cells](pbmc_vignette.html)

In this tutorial we analyze a human peripheral blood mononuclear cell (PBMC)
dataset of ~7,000 cells.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/pbmc.png")
```

## [Mouse cortical brain cells](mouse_brain_vignette.html)

In this tutorial we analyze a dataset of ~3,500 cortical neurons from the adult mouse
brain.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/mouse.png")
```

## [Joint scRNA-seq and scATAC-seq analysis: 10x Multiomic](pbmc_multiomic.html)

In this tutorial we demonstrate a joint analysis of combined gene expression and
DNA accessibility data, measured in the same human PBMCs using the 10x Genomics
multiomic kit.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/multiomic.png")
```


## [Joint scRNA-seq and scATAC-seq analysis: SNARE-seq](snareseq.html)

In this tutorial we demonstrate strategies to analyze a SNARE-seq dataset where
we have paired measurements of gene expression and DNA accessibility from the 
same mouse brain nuclei.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/snareseq.png")
```

## [Joint single-cell mitochondrial DNA genotyping and DNA accessibility analysis](mito.html)

In this tutorial we identify clonotypes using mitochondrial DNA mutations identified
from scATAC-seq data, and jointly analyze clonal cellular relationships and DNA
accessibility patterns in a human colorectal cancer sample.

```{r echo=FALSE, out.width="50%"}
knitr::include_graphics("./assets/mito.png")
```

# How-to

The following short vignettes demonstrate how to perform more specialized 
analysis tasks.

## [Peak calling](peak_calling.html)

In this vignette we demonstrate how to perform cell-type-specific peak
calling for scATAC-seq data.

```{r echo=FALSE, out.height="50%"}
knitr::include_graphics("./assets/peaks.png")
```

## [Merging datasets](merging.html)

This vignette outlines strategies for merging different single-cell chromatin
datasets together.

```{r echo=FALSE, out.width="50%"}
knitr::include_graphics("./assets/merge.png")
```

## [Integration and label transfer](integrate_atac.html)

Here we demonstrate the integration of multiple single-cell chromatin datasets,
as well as label transfer from a reference dataset to an unlabeled query dataset.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/integration.png")
```

## [DNA sequence motif enrichment analysis](motif_vignette.html)

In this vignette we demonstrate how to perform DNA sequence motif enrichment
analysis using Signac.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/motifs.png")
```

## [Transcription factor footprinting analysis](footprint.html)

In this vignette we demonstrate how to perform motif footprinting analysis,
using a human hematopoietic stem cell dataset as an example.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/footprint.png")
```

## [Building trajectories with Monocle 3](monocle.html)

Here we demonstrate how to build trajectories using scATAC-seq data with the 
Monocle 3 package and conversion functions present in SeuratWrappers.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/monocle.png")
```

## [Finding co-accessible sites with Cicero](cicero.html)

Here we demonstrate how to find co-accessible peaks in scATAC-seq data using the
Cicero package and conversion functions present in SeuratWrappers.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/cicero.png")
```

## [Data visualization](visualization.html)

Here we demonstrate how to create genome browser-style plots using single-cell
chromatin data.

```{r echo=FALSE, out.height="30%"}
knitr::include_graphics("./assets/viz.png")
```

# Object interaction

The following vignettes demonstrate how to interact with the Seurat object and 
object classes defined in the Signac package.

## [Data structures and object interaction](data_structures.html)

This vignette details each class defined in Signac, the methods that operate on
each class, and provides some examples of how to interact with these objects to
perform common analysis tasks.

## [Parallel and distributed computing](future.html)

This vignette demonstrates how to enable parallel computing in Signac and Seurat,
and gives an example of the amount of speedup that might be expected from 
enabling parallelization. 
---
title: "Analyzing PBMC scATAC-seq"
output: html_document
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
---
  
```{r init, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

For this tutorial, we will be analyzing a single-cell ATAC-seq dataset of human
peripheral blood mononuclear cells (PBMCs) provided by 10x Genomics. The
following files are used in this vignette, all available through the 10x
Genomics website:

* The [Raw data](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5)  
* The [Metadata](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv)  
* The [fragments file](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz)
* The fragments file [index](https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi)

<details>
  <summary>**View data download code**</summary>

To download all the required files, you can run the following lines in a shell:

```{sh, eval=FALSE}
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi
```

</details>

First load in Signac, Seurat, and some other packages we will be using for
analyzing human data.

```{r setup, message=FALSE}
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)
```

## Pre-processing workflow

When pre-processing chromatin data, Signac uses information from two related
input files, both of which can be created using CellRanger:
  
  * **Peak/Cell matrix**. This is analogous to the gene expression count matrix used to analyze single-cell RNA-seq. However, instead of genes, each row of the matrix represents a region of the genome (a peak), that is predicted to represent a region of open chromatin. Each value in the matrix represents the number of Tn5 integration sites for each single barcode (i.e. a cell) that map within each peak. You can find more detail on the [10X Website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/matrices).

  * **Fragment file**. This represents a full list of all unique fragments across all single cells. It is a substantially larger file, is slower to work with, and is stored on-disk (instead of in memory). However, the advantage of retaining this file is that it contains all fragments associated with each single cell, as opposed to only fragments that map to peaks. More information about the fragment file can be found on the 10x Genomics [website](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments) or on the [sinto website](https://timoast.github.io/sinto/basic_usage.html#create-scatac-seq-fragments-file).

We start by creating a Seurat object using the peak/cell matrix and cell
metadata generated by `cellranger-atac`, and store the path to the fragment
file on disk in the Seurat object:
  
```{r warning=FALSE, message=FALSE}
counts <- Read10X_h5(filename = "../vignette_data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../vignette_data/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)
```

```{r}
pbmc
```

The ATAC-seq data is stored using a custom assay, the `ChromatinAssay`. This
enables some specialized functions for analysing genomic single-cell assays such
as scATAC-seq. By printing the assay we can see some of the additional
information that can be contained in the `ChromatinAssay`, including motif
information, gene annotations, and genome information.

```{r}
pbmc[['peaks']]
```

For example, we can call `granges` on a Seurat object with a `ChromatinAssay`
set as the active assay (or on a `ChromatinAssay`) to see the genomic ranges
associated with each feature in the object. See the
[object interaction vignette](data_structures.html) for more information
about the `ChromatinAssay` class.

```{r}
granges(pbmc)
```

We can also add gene annotations to the `pbmc` object for the human genome.
This will allow downstream functions to pull the gene annotation information
directly from the object.

```{r annotations, message=FALSE, warning=FALSE, cache=FALSE}
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc) <- annotations
```

## Computing QC Metrics

We can now compute some QC metrics for the scATAC-seq experiment. We currently
suggest the following metrics below to assess data quality. As with scRNA-seq,
the expected range of values for these parameters will vary depending on your
biological system, cell viability, and other factors.

* Nucleosome banding pattern: The histogram of DNA fragment sizes (determined from the paired-end sequencing reads) should exhibit a strong nucleosome banding pattern corresponding to the length of DNA wrapped around a single nucleosome. We calculate this per single cell, and quantify the approximate ratio of mononucleosomal to nucleosome-free fragments (stored as `nucleosome_signal`)

* Transcriptional start site (TSS) enrichment score. The [ENCODE project](https://www.encodeproject.org/) has defined an ATAC-seq targeting score based on the ratio of fragments centered at the TSS to fragments in TSS-flanking regions (see https://www.encodeproject.org/data-standards/terms/). Poor ATAC-seq experiments typically will have a low TSS enrichment score. We can compute this metric for each cell with the `TSSEnrichment()` function, and the results are stored in metadata under the column name `TSS.enrichment`.

* Total number of fragments in peaks: A measure of cellular sequencing depth / complexity. Cells with very few reads may need to be excluded due to low sequencing depth. Cells with extremely high levels may represent doublets, nuclei clumps, or other artefacts.

* Fraction of fragments in peaks: Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with low values (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed. Note that this value can be sensitive to the set of peaks used.

* Ratio reads in genomic blacklist regions The [ENCODE project](https://www.encodeproject.org/) has provided a list of [blacklist regions](https://github.com/Boyle-Lab/Blacklist), representing reads which are often associated with artefactual signal. Cells with a high proportion of reads mapping to these areas (compared to reads mapping to peaks) often represent technical artifacts and should be removed. ENCODE blacklist regions for human (hg19 and GRCh38), mouse (mm10), Drosophila (dm3), and C. elegans (ce10) are included in the Signac package.

Note that the last three metrics can be obtained from the output of CellRanger
(which is stored in the object metadata), but can also be calculated for non-10x
datasets using Signac (more information at the end of this document).

```{r message=FALSE, warning=FALSE}
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
```

We can inspect the TSS enrichment scores by grouping the cells based on the
score and plotting the accessibility signal over all TSS sites. Setting the 
`fast=TRUE` option in `TSSEnrichment()` will only compute the TSS enrichment
score without storing the entire cell by position matrix of Tn5 insertion
frequency for each cell, and can save memory. However, setting `fast=TRUE` will
not allow downstream plotting of the TSS enrichment signal for different groups
of cells using the `TSSPlot()` function, shown here:

```{r message=FALSE, warning=FALSE}
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
```

We can also look at the fragment length periodicity for all the cells, and group
by cells with high or low nucleosomal signal strength. You can see that cells
that are outliers for the mononucleosomal / nucleosome-free ratio (based on the
plots above) have different nucleosomal banding patterns. The remaining cells
exhibit a pattern that is typical for a successful ATAC-seq experiment.

```{r message=FALSE, warning=FALSE}
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
```

```{r message=FALSE, warning=FALSE, fig.width=18, fig.height=6}
VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
```

Finally we remove cells that are outliers for these QC metrics.

```{r}
pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc
```

## Normalization and linear dimensional reduction

* Normalization: Signac performs term frequency-inverse document frequency (TF-IDF) normalization. This is a two-step normalization procedure, that both normalizes across cells to correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.

* Feature selection: The low dynamic range of scATAC-seq data makes it challenging to perform variable feature selection, as we do for scRNA-seq. Instead, we can choose to use only the top _n_% of features (peaks) for dimensional reduction, or remove features present in less than _n_ cells with the `FindTopFeatures()` function. Here, we will all features, though we note that we see very similar results when using only a subset of features (try setting min.cutoff to 'q75' to use the top 25% all peaks), with faster runtimes. Features used for dimensional reduction are automatically set as `VariableFeatures()` for the Seurat object by this function.

* Dimension reduction: We next run singular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above. This returns a reduced dimension representation of the object (for users who are more familiar with scRNA-seq, you can think of this as analogous to the output of PCA).

The combined steps of TF-IDF followed by SVD are known as latent semantic
indexing (LSI), and were first introduced for the analysis of scATAC-seq data by
[Cusanovich et al. 2015](https://science.sciencemag.org/content/367/6473/45.full).

```{r message=FALSE, warning=FALSE}
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
```

The first LSI component often captures sequencing depth (technical variation)
rather than biological variation. If this is the case, the component should be
removed from downstream analysis. We can assess the correlation between each LSI
component and sequencing depth using the `DepthCor()` function:

```{r}
DepthCor(pbmc)
```

Here we see there is a very strong correlation between the first LSI component
and the total number of counts for the cell, so we will perform downstream steps
without this component.

## Non-linear dimension reduction and clustering

Now that the cells are embedded in a low-dimensional space, we can use methods
commonly applied for the analysis of scRNA-seq data to perform graph-based
clustering and non-linear dimension reduction for visualization. The functions
`RunUMAP()`, `FindNeighbors()`, and `FindClusters()` all come from the Seurat
package.

```{r message=FALSE, warning=FALSE}
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()
```

## Create a gene activity matrix

The UMAP visualization reveals the presence of multiple cell groups in human 
blood. If you are familiar with our
[scRNA-seq analyses of PBMC](https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html),
you may even recognize the presence of certain myeloid and lymphoid populations
in the scATAC-seq data. However, annotating and interpreting clusters is more
challenging in scATAC-seq data as much less is known about the functional roles
of noncoding genomic regions than is known about protein coding regions (genes). 

However, we can try to quantify the activity of each gene in the genome by
assessing the chromatin accessibility associated with each gene, and create a
new gene activity assay derived from the scATAC-seq data. Here we will use a
simple approach of summing the fragments intersecting the gene body and promoter
region (we also recommend exploring the
[Cicero](https://cole-trapnell-lab.github.io/cicero-release/) tool, which can
accomplish a similar goal, and we provide a vignette showing how to run Cicero
within a Signac workflow [here](cicero.html)). 

To create a gene activity matrix, we extract gene coordinates and extend them
to include the 2 kb upstream region (as promoter accessibility is often
correlated with gene expression). We then count the number of fragments for each
cell that map to each of these regions, using the using the `FeatureMatrix()`
function. These steps are automatically performed by the `GeneActivity()`
function:

```{r, message=FALSE, warning=FALSE}
gene.activities <- GeneActivity(pbmc)
```

```{r, message=FALSE, warning=FALSE}
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)
```

Now we can visualize the activities of canonical marker genes to help interpret
our ATAC-seq clusters. Note that the activities will be much noisier than
scRNA-seq measurements. This is because they represent measurements from sparse
chromatin data, and because they assume a general correspondence between gene
body/promoter accessibility and gene expression which may not always be the
case. Nonetheless, we can begin to discern populations of monocytes, B, T, and
NK cells based on these gene activity profiles. However, further subdivision of
these cell types is challenging based on supervised analysis alone.

```{r, fig.width=12, fig.height=10}
DefaultAssay(pbmc) <- 'RNA'

FeaturePlot(
  object = pbmc,
  features = c('MS4A1', 'CD3D', 'LEF1', 'NKG7', 'TREM1', 'LYZ'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
```

## Integrating with scRNA-seq data

To help interpret the scATAC-seq data, we can classify cells based on an
scRNA-seq experiment from the same biological system (human PBMC). We utilize
methods for cross-modality integration and label transfer, described 
[here](https://doi.org/10.1016/j.cell.2019.05.031), with a more in-depth
tutorial [here](https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html).
We aim to identify shared correlation patterns in the gene activity matrix and
scRNA-seq dataset to identify matched biological states across the two
modalities. This procedure returns a classification score for each cell for each
scRNA-seq-defined cluster label.

Here we load a pre-processed scRNA-seq dataset for human PBMCs, also provided by
10x Genomics. You can download the raw data for this experiment from the 10x 
[website](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_v3),
and view the code used to construct this object on
[GitHub](https://github.com/satijalab/Integration2019/blob/master/preprocessing_scripts/pbmc_10k_v3.R).
Alternatively, you can download the pre-processed Seurat object 
[here](https://www.dropbox.com/s/zn6khirjafoyyxl/pbmc_10k_v3.rds).

```{r warning=FALSE, message=FALSE}
# Load the pre-processed scRNA-seq data for PBMCs
pbmc_rna <- readRDS("../vignette_data/pbmc_10k_v3.rds")

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
```

```{r fig.width=12}
plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2
```

You can see that the scRNA-based classifications are entirely consistent with
the UMAP visualization, computed only on the scATAC-seq data. We can now easily
annotate our scATAC-seq-derived clusters (alternatively, we could use the RNA
classifications themselves). We note that cluster 14 maps to CD4 Memory T cells,
but is a very small cluster with lower QC metrics. As this group is likely
representing low-quality cells, we remove it from downstream analysis. 

```{r, cache=FALSE}
pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)
```

## Find differentially accessible peaks between clusters

To find differentially accessible regions between clusters of cells, we can
perform a differential accessibility (DA) test. We utilize logistic regression
for DA, as suggested by
[Ntranos et al. 2018](https://www.biorxiv.org/content/10.1101/258566v2)
for scRNA-seq data, and add the total number of fragments as a latent variable
to mitigate the effect of differential sequencing depth on the result.
For sparse data (such as scATAC-seq), we find it is often necessary
to lower the `min.pct` threshold in `FindMarkers()` from the default (0.1, which
was designed for scRNA-seq data). Here we
will focus on comparing Naive CD4 cells and CD14 monocytes, but any groups of
cells can be compared using these methods. We can also visualize these marker
peaks on a violin plot, feature plot, dot plot, heat map, or any
[visualization tool in Seurat](https://satijalab.org/seurat/v3.0/visualization_vignette.html).

```{r message=TRUE, warning=FALSE}
# change back to working with peaks instead of gene activities
DefaultAssay(pbmc) <- 'peaks'

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "CD14 Mono",
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)
```

```{r, cache=FALSE, fig.width=12}
plot1 <- VlnPlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("CD4 Naive","CD14 Mono")
)
plot2 <- FeaturePlot(
  object = pbmc,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2
```

Another way to find DA regions between two groups of cells is to look at the
fold change accessibility between two groups of cells. This can be much faster
than running more sophisticated DA tests, but is not able to account for
latent variables such as differences in the total sequencing depth between
cells, and does not perform any statistical test. However, this can still be a
useful way to quickly explore data, and can be performed using the `FoldChange()`
function in Seurat.

```{r message=FALSE, warning=FALSE}
fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14 Mono")
head(fc)
```

Peak coordinates can be difficult to interpret alone. We can find the closest
gene to each of these peaks using the `ClosestFeature()` function. If you
explore the gene lists, you will see that peaks open in Naive T cells are close
to genes such as *BCL11B* and *GATA3* (key regulators of  T cell differentiation
), while peaks open in monocytes are close to genes such as *CEBPB* (a key 
regulator of monocyte differentiation). We could follow up this result further
by doing gene ontology enrichment analysis on the gene sets returned by
`ClosestFeature()`, and there are many R packages that can do this (see the
[`GOstats`](https://bioconductor.org/packages/release/bioc/html/GOstats.html)
package for example).

```{r, warning=FALSE, message=FALSE}
open_cd4naive <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)
```

```{r, warning=FALSE, message=FALSE}
head(closest_genes_cd4naive)
```

```{r, warning=FALSE, message=FALSE}
head(closest_genes_cd14mono)
```

## Plotting genomic regions

We can plot the frequency of Tn5 integration across regions of the genome for
cells grouped by cluster, cell type, or any other metadata stored in the object
for any genomic region using the `CoveragePlot()` function.
These represent pseudo-bulk accessibility tracks, where signal from
all cells within a group have been averaged together to visualize the DNA 
accessibility in a region (thanks to Andrew Hill for giving the inspiration for
this function in his excellent [blog post](http://andrewjohnhill.com/blog/2019/04/12/streamlining-scatac-seq-visualization-and-analysis/)). Alongside these accessibility tracks we can visualize other
important information including gene annotation, peak coordinates, and genomic
links (if they're stored in the object). See the
[visualization vignette](visualization.html) for more information.

```{r message=FALSE, warning=FALSE, out.width='90%', fig.height=6}
# set plotting order
levels(pbmc) <- c("CD4 Naive","CD4 Memory","CD8 Naive","CD8 Effector","DN T","NK CD56bright","NK CD56Dim","pre-B",'pro-B',"pDC","DC","CD14 Mono",'CD16 Mono')

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)
```

We can also create an interactive version of these plots using the
`CoverageBrowser()` function. Here is a recorded demonstration showing how we
can use `CoverageBrowser()` to browser the genome and adjust plotting parameters
interactively. The "Save plot" button in `CoverageBrowser()` will add the
current plot to a list of `ggplot` objects that is returned when the browser
session is ended by pressing the "Done" button, allowing interesting views
to be saved during an interactive session.

<iframe width="560" height="315" src="https://www.youtube.com/embed/S9b5rN32IC8" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<details>
  <summary>**Working with datasets that were not quantified using CellRanger**</summary>

The CellRanger software from 10x Genomics generates several useful QC metrics
per-cell, as well as a peak/cell matrix and an indexed fragments file. In the
above vignette, we utilize the CellRanger outputs, but provide alternative
functions in Signac for many of the same purposes here.

### Generating a peak/cell or bin/cell matrix

The `FeatureMatrix` function can be used to generate a count matrix containing
any set of genomic ranges in its rows. These regions could be a set of peaks, or
bins that span the entire genome.

```{r, eval=FALSE}
# not run
# peak_ranges should be a set of genomic ranges spanning the set of peaks to be quantified per cell
peak_matrix <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peak_ranges
)
```

For convenience, we also include a `GenomeBinMatrix()` function that will
generate a set of genomic ranges spanning the entire genome for you, and run
`FeatureMatrix()` internally to produce a genome bin/cell matrix.

```{r, eval=FALSE}
# not run
bin_matrix <- GenomeBinMatrix(
  fragments = Fragments(pbmc),
  genome = seqlengths(pbmc),
  binsize = 5000
)
```

### Counting fraction of reads in peaks

The function `FRiP()` will count the fraction of reads in peaks for each cell,
given a peak/cell assay and a bin/cell assay. Note that this can be run on a
subset of the genome, so that a bin/cell assay does not need to be computed for
the whole genome. This will return a Seurat object will metadata added
corresponding to the fraction of reads in peaks for each cell.

```{r, eval=FALSE}
# not run
total_fragments <- CountFragments("'../vignette_data/atac_v1_pbmc_10k_fragments.tsv.gz'")
pbmc$fragments <- total_fragments[colnames(pbmc), "frequency_count"]

pbmc <- FRiP(
  object = pbmc,
  assay = 'peaks',
  total.fragments = 'fragments'
)
```

### Counting fragments in genome blacklist regions

The ratio of reads in genomic blacklist regions, that are known to artifactually
accumulate reads in genome sequencing assays, can be diagnostic of low-quality
cells. We provide blacklist region coordinates for several genomes (hg19, hg38,
mm9, mm10, ce10, ce11, dm3, dm6) in the Signac package for convenience. These
regions were provided by the ENCODE consortium, and we encourage users to cite
their [paper](https://doi.org/10.1038/s41598-019-45839-z) if you use the regions
in your analysis. The `FractionCountsInRegion()` function can be used to
calculate the fraction of all counts within a given set of regions per cell.
We can use this function and the blacklist regions to find the fraction of
blacklist counts per cell.

```{r, eval=FALSE}
# not run
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg19
)
```

```{r message=FALSE, warning=FALSE, echo=FALSE}
saveRDS(object = pbmc, file = "../vignette_data/pbmc.rds")
```

</details>

<details>
  <summary>**Session Info**</summary>

```{r}
sessionInfo()
```

</details>
---
title: "Visualization of genomic regions"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

In this vignette we will demonstrate how to visualize single-cell data in 
genome-browser-track style plots with Signac.

To demonstrate we'll use the human PBMC dataset processed in
[this vignette](pbmc_vignette.html).

```{r message=FALSE, warning=FALSE}
library(Signac)

# load PBMC dataset
pbmc <- readRDS("../vignette_data/pbmc.rds")
```

There are several different genome browser style plot types available in Signac,
including accessibility tracks, gene annotations, peak coordinate, genomic links,
and fragment positions.

## Plotting aggregated signal

The main plotting function in Signac is `CoveragePlot()`, and this computes the
averaged frequency of sequenced DNA fragments for different groups of cells
within a given genomic region. 

```{r message=FALSE, warning=FALSE}
cov_plot <- CoveragePlot(
  object = pbmc,
  region = "chr2-87011729-87035519",
  annotation = FALSE,
  peaks = FALSE
)
cov_plot
```

We can also request regions of the genome by gene name. This will use the gene
coordinates stored in the Seurat object to determine which genomic region to
plot

```{r message=FALSE, warning=FALSE}
CoveragePlot(
  object = pbmc,
  region = "CD8A",
  annotation = FALSE,
  peaks = FALSE
)
```

## Plotting gene annotations

Gene annotations within a given genomic region can be plotted using the 
`AnnotationPlot()` function.

```{r message=FALSE, warning=FALSE}
gene_plot <- AnnotationPlot(
  object = pbmc,
  region = "chr2-87011729-87035519"
)
gene_plot
```

## Plotting peak coordinates

Peak coordinates within a genomic region can be plotted using the `PeakPlot()`
function.

```{r message=FALSE, warning=FALSE}
peak_plot <- PeakPlot(
  object = pbmc,
  region = "chr2-87011729-87035519"
)
peak_plot
```

## Plotting genomic links

Relationships between genomic positions can be plotted using the `LinkPlot()`
function. This will display an arc connecting two linked positions, with the
transparency of the arc line proportional to a score associated with the link.
These links could be used to encode different things, including regulatory
relationships (for example, linking enhancers to the genes that they regulate),
or experimental data such as Hi-C.

Just to demonstrate how the function works, we've created a fake link here
and added it to the PBMC dataset.

```{r include=FALSE}
library(GenomicRanges)

# generate fake links and add to object
cd8.coords <- LookupGeneCoords(object = pbmc, gene = "CD8A")

links <- GRanges(
  seqnames = "chr2",
  ranges = IRanges(start = start(cd8.coords)[1] + 100, end = end(cd8.coords)[1] - 5000),
  score = 1, 
  group = 1
)

Links(pbmc) <- links
```

```{r message=FALSE, warning=FALSE}
link_plot <- LinkPlot(
  object = pbmc,
  region = "chr2-87011729-87035519"
)
link_plot
```

## Plotting per-cell fragment abundance

While the `CoveragePlot()` function computes an aggregated signal within a 
genomic region for different groups of cells, sometimes it's also useful to
inspect the frequency of sequenced fragments within a genomic region for
_individual_ cells, without aggregation. This can be done using the `TilePlot()`
function.

```{r message=FALSE, warning=FALSE}
tile_plot <- TilePlot(
  object = pbmc,
  region = "chr2-87011729-87035519",
  idents = c("CD4 Memory", "CD8 Effector")
)
tile_plot
```

By default, this selects the top 100 cells for each group based on the total
number of fragments in the genomic region. The genomic region is then tiled 
and the total fragments in each tile counted for each cell, and the resulting
counts for each position displayed as a heatmap.

## Plotting additional data alongside genomic tracks

Multimodal single-cell datasets generate multiple experimental measurements for
each cell. Several methods now exist that are capable of measuring single-cell
chromatin data (such as chromatin accessibility) alongside other measurements
from the same cell, such as gene expression or mitochondrial genotype. In these
cases it's often informative to visualize the multimodal data together in a
single plot. This can be achieved using the `ExpressionPlot()` function. This
is similar to the `VlnPlot()` function in Seurat, but is designed to be
incorportated with genomic track plots generated by `CoveragePlot()`.

```{r message=FALSE, warning=FALSE}
expr_plot <- ExpressionPlot(
  object = pbmc,
  features = "CD8A",
  assay = "RNA"
)
expr_plot
```

We can create similar plots for multiple genes at once simply by passing a list
of gene names

```{r message=FALSE, warning=FALSE}
ExpressionPlot(
  object = pbmc,
  features = c("CD8A", "CD4"),
  assay = "RNA"
)
```

## Combining genomic tracks

Above we've demonstrated how to generate individual tracks and panels that can
be combined into a single plot for a single genomic region. These panels can 
be easily combined using the `CombineTracks()` function.

```{r fig.height=10, message=FALSE, warning=FALSE}
CombineTracks(
  plotlist = list(cov_plot, tile_plot, peak_plot, gene_plot, link_plot),
  expression.plot = expr_plot,
  heights = c(10, 6, 1, 2, 3),
  widths = c(10, 1)
)
```

The `heights` and `widths` parameters control the relative heights and widths of
the individual tracks, according to the order that the tracks appear in the
`plotlist`. The `CombineTracks()` function ensures that the tracks are aligned
vertically and horizontally, and moves the x-axis labels (describing the 
genomic position) to the bottom of the combined tracks.

## Generating multiple tracks

Above we've shown how to create genomic plot panels individually and how to
combine them. This allows more control over how each panel is constructed and 
how they're combined, but involves multiple steps. For convenience, we've
included the ability to generate and combine different panels automatically in
the `CoveragePlot()` function, through the `annotation`, `peaks`, `tile`, and
`features` arguments. We can generate a similar plot to that shown above in a
single function call:

```{r fig.height=10, message=FALSE, warning=FALSE}
CoveragePlot(
  object = pbmc,
  region = "chr2-87011729-87035519",
  features = "CD8A",
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE,
  links = TRUE
)
```

Notice that in this example we create the tile plot for every group of cells
that is shown in the coverage track, whereas above we were able to create a plot
that showed the aggregated coverage for all groups of cells and the tile plot
for only the CD4 memory cells and the CD8 effector cells. A higher degree of 
customization is possible when creating each track separately.

## Interactive visualization

Above we demonstrated the different types of plots that can be constructed using
Signac. Often when exploring genomic data it's useful to be able to
interactively browse through different regions of the genome and adjust tracks
on the fly. This can be done in Signac using the `CoverageBrowser()` function.
This provides all the same functionality of the `CoveragePlot()` function,
except that we can scroll upstream/downstream, zoom in/out of regions, navigate
to new region, and adjust which tracks are shown or how the cells are grouped. 
In exploring the data interactively, often you will find interesting plots that
you'd like to save for viewing later. We've included a "Save plot" button that 
will add the current plot to a list of plots that is returned when the
interactive session is ended. Here's a recorded demonstration of the 
`CoverageBrowser()` function:

<iframe width="560" height="315" src="https://www.youtube.com/embed/S9b5rN32IC8" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>


---
title: "Joint RNA and ATAC analysis: SNARE-seq"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

In this vignette we will analyse a single-cell co-assay dataset measuring gene
expression and DNA accessibility in the same cells. This vignette is similar to
the [PBMC multiomic vignette](pbmc_multiomic.html), but demonstrates a similar
joint analysis in a different species and with data gathered using a different technology.

This dataset was published
by [Chen, Lake, and Zhang (2019)](https://www.nature.com/articles/s41587-019-0290-0)
and uses a technology called SNARE-seq. We will look at a dataset from the adult
mouse brain, and the raw data can be downloaded from NCBI GEO here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126074 

As the fragment files for this dataset are not publicly available we have
re-mapped the raw data to the mm10 genome and created a fragment file using
[Sinto](https://timoast.github.io/sinto/).

The fragment file can be downloaded
here: https://www.dropbox.com/s/se3uag0cvld4xqg/fragments.sort.bed.gz

The fragment file index can be downloaded here: https://www.dropbox.com/s/kvb8fiewkjz3kzf/fragments.sort.bed.gz.tbi

Code used to create the fragment file from raw data is available here:
https://github.com/timoast/SNARE-seq 

## Data loading

First we create a Seurat object containing two different assays, one containing
the gene expression data and one containing the DNA accessibility data.

To load the count data, we can use the `Read10X()` function from Seurat by first
placing the `barcodes.tsv.gz`, `matrix.mtx.gz`, and `features.tsv.gz` files into
a separate folder.

```{r include=FALSE}
if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
    BiocManager::install("EnsDb.Mmusculus.v79")
```

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

# load processed data matrices for each assay
rna <- Read10X("../vignette_data/snare-seq/GSE126074_AdBrainCortex_rna/", gene.column = 1)
atac <- Read10X("../vignette_data/snare-seq/GSE126074_AdBrainCortex_atac/", gene.column = 1)
fragments <- "../vignette_data/snare-seq/fragments.sort.bed.gz"

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations
```

## Quality control

```{r message=FALSE, warning=FALSE}
DefaultAssay(snare) <- "ATAC"
snare <- TSSEnrichment(snare)
snare <- NucleosomeSignal(snare)
snare$blacklist_fraction <- FractionCountsInRegion(
  object = snare,
  assay = 'ATAC',
  regions = blacklist_mm10
)
```

```{r fig.width=18, message=FALSE, warning=FALSE}
Idents(snare) <- "all"  # group all cells together, rather than by replicate
VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)
```

```{r}
snare <- subset(
  x = snare,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)
snare
```

## Gene expression data processing

Process gene expression data using Seurat

```{r message=FALSE, warning=FALSE}
DefaultAssay(snare) <- "RNA"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 30)
snare <- RunUMAP(snare, dims = 1:30, reduction.name = "umap.rna")
snare <- FindNeighbors(snare, dims = 1:30)
snare <- FindClusters(snare, resolution = 0.5, algorithm = 3)
p1 <- DimPlot(snare, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")
```

## DNA accessibility data processing

Process the DNA accessibility data using Signac

```{r message=FALSE, warning=FALSE}
DefaultAssay(snare) <- 'ATAC'

snare <- FindTopFeatures(snare, min.cutoff = 10)
snare <- RunTFIDF(snare)
snare <- RunSVD(snare)
snare <- RunUMAP(snare, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')
p2 <- DimPlot(snare, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")
```

```{r fig.width=12, message=FALSE, warning=FALSE}
p1 + p2
```

## Integration with scRNA-seq data

Next we can annotate cell types in the dataset by transferring labels from an
existing scRNA-seq dataset for the adult mouse brain, produced by the Allen
Institute.

```{r message=FALSE, warning=FALSE}
# label transfer from Allen brain
allen <- readRDS("../vignette_data/allen_brain.rds")

# use the RNA assay in the SNARE-seq data for integration with scRNA-seq
DefaultAssay(snare) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = snare,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = snare[['pca']],
  dims = 1:30
)

snare <- AddMetaData(object = snare, metadata = predicted.labels)
```

```{r}
# label clusters based on predicted ID
new.cluster.ids <- c(
  "L2/3 IT",
  "L4",
  "L6 IT",
  "L5 CT",
  "L4",
  "L5 PT",
  "Pvalb",
  "Sst",
  "Astro",
  "Oligo",
  "Vip/Lamp5",
  "L6 IT.2",
  "L6b",
  "NP"
)
names(x = new.cluster.ids) <- levels(x = snare)
snare <- RenameIdents(object = snare, new.cluster.ids)
snare$celltype <- Idents(snare)
DimPlot(snare, group.by = 'celltype', label = TRUE, reduction = 'umap.rna')
```

## Jointly visualizing gene expression and DNA accessibility

We can visualize both the gene expression and DNA accessibility information at
the same time using the `CoveragePlot()` function. This makes it easy to compare
DNA accessibility in a given region for different cell types and overlay gene 
expression information for different genes.

```{r message=FALSE, warning=FALSE, fig.width=12, fig.height=8}
DefaultAssay(snare) <- "ATAC"
CoveragePlot(snare, region = "chr2-22620000-22660000", features = "Gad2")
```

```{r include=FALSE}
saveRDS(object = snare, file = "../vignette_data/snare.rds")
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Joint RNA and ATAC analysis: 10x multiomic"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

In this vignette, we'll demonstrate how to jointly analyze a single-cell dataset
measuring both DNA accessibility and gene expression in the same cells using
Signac and Seurat. In this vignette we'll be using a publicly available 10x
Genomic Multiome dataset for human PBMCs.

<details>
  <summary>**View data download code**</summary>

You can download the required data by running the following lines in a shell:

```{sh, eval=FALSE}
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
```

</details>

```{r include=FALSE}
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
    BiocManager::install("EnsDb.Hsapiens.v86")

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

if (!requireNamespace("SeuratDisk", quietly = TRUE))
    remotes::install_github("mojaveazure/seurat-disk")
```

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

set.seed(1234)
```

```{r message=FALSE, warning=FALSE}
# load the RNA and ATAC data
counts <- Read10X_h5("../vignette_data/multiomic/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "../vignette_data/multiomic/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
```

```{r message=FALSE, warning=FALSE}
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
```

```{r}
pbmc
```

## Quality control

We can compute per-cell quality control metrics using the DNA accessibility data
and remove cells that are outliers for these metrics, as well as cells with
low or unusually high counts for either the RNA or ATAC assay.

```{r fig.width=18, message=FALSE, warning=FALSE}
DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(
  object = pbmc,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
```

```{r}
# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
pbmc
```
## Peak calling

The set of peaks identified using Cellranger often merges distinct peaks that are
close together. This can create a problem for certain analyses, particularly motif 
enrichment analysis and peak-to-gene linkage. To identify a more accurate set of peaks,
we can call peaks using MACS2 with the `CallPeaks()` function. Here we call peaks
on all cells together, but we could identify peaks for each group of cells separately
by setting the `group.by` parameter, and this can help identify peaks specific to
rare cell populations. 

```{r message=FALSE, warning=FALSE}
# call peaks using MACS2
peaks <- CallPeaks(pbmc, macs2.path = "/home/stuartt/miniconda3/envs/signac/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)
```

## Gene expression data processing

We can normalize the gene expression data using SCTransform, and reduce the dimensionality
using PCA.

```{r message=FALSE, warning=FALSE, results='hide'}
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)
```

## DNA accessibility data processing

Here we process the DNA accessibility assay the same way we would process a 
scATAC-seq dataset, by performing latent semantic indexing (LSI).

```{r message=FALSE, warning=FALSE}
DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
```

## Annotating cell types

To annotate cell types in the dataset we can transfer cell labels from an existing
PBMC reference dataset using tools in the Seurat package. See the
Seurat reference mapping [vignette](https://satijalab.org/seurat/v4.0/reference_mapping.html)
for more information.

We'll use an annotated PBMC reference dataset from [Hao et al. (2020)](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1),
available for download here: https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat

Note that the SeuratDisk package is required to load the reference dataset.
Installation instructions for SeuratDisk can be found [here](https://github.com/mojaveazure/seurat-disk).

```{r message=FALSE, warning=FALSE}
library(SeuratDisk)

# load PBMC reference
reference <- LoadH5Seurat("../vignette_data/multiomic/pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                 "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                 "NK Proliferating", "gdT",
                 "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                 "CD14 Mono", "CD16 Mono",
                 "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")
```

## Joint UMAP visualization

Using the weighted nearest neighbor methods in [Seurat v4](https://www.biorxiv.org/content/10.1101/2020.10.12.335331v1),
we can compute a joint neighbor graph that represent both the gene expression and 
DNA accessibility measurements.

```{r message=FALSE, warning=FALSE}
# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
```


## Linking peaks to genes

For each gene, we can find the set of peaks that may regulate the gene by by
computing the correlation between gene expression and accessibility at nearby
peaks, and correcting for bias due to GC content, overall accessibility, and 
peak size. See the [Signac paper](https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1)
for a full description of the method we use to link peaks to genes.

Running this step on the whole genome can be time consuming, so here we demonstrate
peak-gene links for a subset of genes as an example. The same function can be used
to find links for all genes by omitting the `genes.use` parameter:

```{r message=FALSE, warning=FALSE}
DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(
  object = pbmc,
  peak.assay = "peaks",
  expression.assay = "SCT",
  genes.use = c("LYZ", "MS4A1")
)
```

We can visualize these links using the `CoveragePlot()` function, or alternatively
we could use the `CoverageBrowser()` function in an interactive analysis:

```{r fig.height=10, message=FALSE, warning=FALSE}
idents.plot <- c("B naive", "B intermediate", "B memory",
                 "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")

p1 <- CoveragePlot(
  object = pbmc,
  region = "MS4A1",
  features = "MS4A1",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000
)

p2 <- CoveragePlot(
  object = pbmc,
  region = "LYZ",
  features = "LYZ",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 8000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)
```

```{r include=FALSE}
saveRDS(object = pbmc, file = "../vignette_data/pbmc_multiomic.rds")
```


<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Joint single-cell mitochondrial DNA genotyping and DNA accessibility analysis"
author: Caleb Lareau and Tim Stuart
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output:
  html_document:
    self_contained: True
---

```{r packages, cache=FALSE, message=FALSE, warning = FALSE, echo = TRUE}
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(EnsDb.Hsapiens.v75)
```

Here, we take a look at two different datasets containing both DNA accessibility
measurements and mitochondrial mutation data in the same cells. One was sampled
from a patient with a colorectal cancer (CRC) tumor, and the other is from a 
polyclonal TF1 cell line. This data was produced by Lareau and Ludwig
et al. (2020), and you can read the original paper
here: https://doi.org/10.1038/s41587-020-0645-6.

Processed data files, including mitochondrial variant data for the CRC and TF1
dataset is available on Zenodo here: https://zenodo.org/record/3977808

Raw sequencing data and DNA accessibility processed files for the these datasets
are available on NCBI GEO here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142745

<details>
  <summary>**View data download code**</summary>

The required files can be downloaded by running the following lines in a shell:

```{bash, eval=FALSE}
# ATAC data
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.filtered_peak_bc_matrix.h5
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.singlecell.csv
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.fragments.tsv.gz
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.fragments.tsv.gz.tbi

# mitochondrial allele data
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.A.txt.gz
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.C.txt.gz
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.G.txt.gz
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.T.txt.gz
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.depthTable.txt
wget https://zenodo.org/record/3977808/files/CRC_v12-mtMask_mgatk.chrM_refAllele.txt
```

</details>


# Colorectal cancer dataset

To demonstrate combined analyses of mitochondrial DNA variants and accessible
chromatin, we'll walk through a vignette analyzing cells from a primary
colorectal adenocarcinoma. The sample contains a mixture of malignant
epithelial cells and tumor infiltrating immune cells. 

## Loading the DNA accessibility data

First we load the scATAC-seq data and create a Seurat object following the 
standard workflow for scATAC-seq data.
  
```{r importData, message=FALSE, warning = FALSE}
# load counts and metadata from cellranger-atac
counts <- Read10X_h5(filename = "../vignette_data/mito/CRC_v12-mtMask_mgatk.filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../vignette_data/mito/CRC_v12-mtMask_mgatk.singlecell.csv",
  header = TRUE,
  row.names = 1
)

# load gene annotations from Ensembl
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# create object
crc_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  annotation = annotations,
  min.cells = 10,
  genome = "hg19",
  fragments = '../vignette_data/mito/CRC_v12-mtMask_mgatk.fragments.tsv.gz'
)
crc <- CreateSeuratObject(
  counts = crc_assay,
  assay = 'peaks',
  meta.data = metadata
)

crc[["peaks"]]
```

## Quality control

We can compute the standard quality control metrics for scATAC-seq and filter
out low-quality cells based on these metrics.

```{r message=FALSE, warning=FALSE}
# Augment QC metrics that were computed by cellranger-atac
crc$pct_reads_in_peaks <- crc$peak_region_fragments / crc$passed_filters * 100
crc$pct_reads_in_DNase <- crc$DNase_sensitive_region_fragments / crc$passed_filters * 100
crc$blacklist_ratio <- crc$blacklist_region_fragments / crc$peak_region_fragments

# compute TSS enrichment score and nucleosome banding pattern
crc <- TSSEnrichment(crc)
crc <- NucleosomeSignal(crc)
```

```{r fig.width=8, fig.height = 8, message=FALSE, warning=FALSE}
# visualize QC metrics for each cell
VlnPlot(crc, c("TSS.enrichment", "nCount_peaks", "nucleosome_signal", "pct_reads_in_peaks", "pct_reads_in_DNase", "blacklist_ratio"), pt.size = 0, ncol = 3)
```

```{r}
# remove low-quality cells
crc <- subset(
  x = crc,
  subset = nCount_peaks > 1000 &
    nCount_peaks < 50000 &
    pct_reads_in_DNase > 40 &
    blacklist_ratio < 0.05 &
    TSS.enrichment > 3 & 
    nucleosome_signal < 4
)
crc
```

## Loading the mitochondrial variant data

Next we can load the mitochondrial DNA variant data for these cells that was 
quantified using [mgatk](https://github.com/caleblareau/mgatk). The `ReadMGATK()`
function in Signac allows the output from `mgatk` to be read directly into R in
a convenient format for downstream analysis with Signac. Here, we load the data
and add it to the Seurat object as a new assay.

```{r process_mito, cache=TRUE, message=FALSE, warning = FALSE, echo = TRUE}
# load mgatk output
mito.data <- ReadMGATK(dir = "../vignette_data/mito/crc/")

# create an assay
mito <- CreateAssayObject(counts = mito.data$counts)

# Subset to cell present in the scATAC-seq assat
mito <- subset(mito, cells = colnames(crc))

# add assay and metadata to the seurat object
crc[["mito"]] <- mito
crc <- AddMetaData(crc, metadata = mito.data$depth, col.name = "mtDNA_depth")
```

We can look at the mitochondrial sequencing depth for each cell, and further
subset the cells based on mitochondrial sequencing depth.

```{r message=FALSE, warning=FALSE}
VlnPlot(crc, "mtDNA_depth", pt.size = 0.1) + scale_y_log10()
```

```{r}
# filter cells based on mitochondrial depth
crc <- subset(crc, mtDNA_depth >= 10)
crc
```

## Dimension reduction and clustering

Next we can run a standard dimension reduction and clustering workflow using the
scATAC-seq data to identify cell clusters.

```{r message=FALSE, warning=FALSE}
crc <- RunTFIDF(crc)
crc <- FindTopFeatures(crc, min.cutoff = 10)
crc <- RunSVD(crc)
crc <- RunUMAP(crc, reduction = "lsi", dims = 2:50)
crc <- FindNeighbors(crc, reduction = "lsi", dims = 2:50)
crc <- FindClusters(crc, resolution = 0.5, algorithm = 3)
```

```{r message=FALSE, warning=FALSE}
DimPlot(crc, label = TRUE) + NoLegend()
```

## Generate gene scores

To help interpret these clusters of cells, and assign a cell type label, we'll
estimate gene activities by summing the DNA accessibility in the gene body and
promoter region.

```{r setup_gene_activity, cache=FALSE, message=FALSE, warning = FALSE, echo = TRUE}
# compute gene accessibility
gene.activities <- GeneActivity(crc)

# add to the Seurat object as a new assay
crc[['RNA']] <- CreateAssayObject(counts = gene.activities)

crc <- NormalizeData(
  object = crc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(crc$nCount_RNA)
)
```

## Visualize interesting gene activity scores

We note the following markers for different cell types in the CRC dataset:

- EPCAM is a marker for epithelial cells
- TREM1 is a meyloid marker
- PTPRC = CD45 is a pan-immune cell marker
- IL1RL1 is a basophil marker
- GATA3 is a Tcell maker

```{r viz_gene_activitites, message=FALSE, warning = FALSE, fig.width=8, fig.height=8}
DefaultAssay(crc) <- 'RNA'

FeaturePlot(
  object = crc,
  features = c('TREM1', 'EPCAM', "PTPRC", "IL1RL1","GATA3", "KIT"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 2
)
```

Using these gene score values, we can assign cluster identities: 

```{r}
crc <- RenameIdents(
  object = crc,
  '0' = 'Epithelial',
  '1' = 'Epithelial',
  '2' = 'Basophil',
  '3' = 'Myeloid_1',
  '4' = 'Myeloid_2',
  '5' = 'Tcell'
)
```

One of the myeloid clusters has a lower percentage of fragments in peaks, as
well as a lower overall mitochondrial sequencing depth and a different
nucleosome banding pattern.

```{r cell_filtering_recap, message=FALSE, warning = FALSE}
p1 <- FeatureScatter(crc, "mtDNA_depth", "pct_reads_in_peaks") + ggtitle("") + scale_x_log10()
p2 <- FeatureScatter(crc, "mtDNA_depth", "nucleosome_signal") + ggtitle("") + scale_x_log10()

p1 + p2 + plot_layout(guides = 'collect')
```

We can see that most of the low FRIP cells were the `myeloid 1` cluster. This is
most likely an intra-tumor granulocyte that has relatively poor accessible
chromatin enrichment. Similarly, the unusual nuclear chromatin packaging of this
cell type yields slightly reduced mtDNA coverage compared to the `myeloid 2`
cluster.

## Find informative mtDNA variants

Next, we can identify sites in the mitochondrial genome that vary across cells, 
and cluster the cells into clonotypes based on the frequency of these variants
in the cells. Signac utilizes the principles established in the original
mtscATAC-seq work of identifying high-quality variants.  

```{r call_variants, message=FALSE, warning = FALSE}
variable.sites <- IdentifyVariants(crc, assay = "mito", refallele = mito.data$refallele)
VariantPlot(variants = variable.sites)
```

The plot above clearly shows a group of variants with a higher VMR and strand
concordance. In principle, a high strand concordance reduces the likelihood of
the allele frequency being driven by sequencing error (which predominately
occurs on one but not the other strand. This is due to the preceding nucleotide 
content and a common error in mtDNA genotyping). On the other hand, variants
that have a high VMR are more likely to be clonal variants as the alternate
alleles tend to aggregate in certain cells rather than be equivalently
dispersed about all cells, which would be indicative of some other artifact. 

We note that variants that have a very low VMR and and very high strand
concordance are homoplasmic variants for this sample. While these may be
interesting in some settings (e.g. donor demultiplexing), for inferring
subclones, these are not particularly useful. 

Based on these thresholds, we can filter out a set of informative
mitochondrial variants that differ across the cells.

```{r look_at_variants, message=FALSE, warning = FALSE}
# Establish a filtered data frame of variants based on this processing
high.conf <- subset(
  variable.sites, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)

high.conf[,c(1,2,5)]
```

A few things stand out. First, 10 out of the 12 variants occur at less than 1% 
allele frequency in the population. However, 16147C>T is present at about 62%.
We'll see that this is a clonal variant marking the epithelial cells.
Additionally, all of the called variants are transitions (A - G or C - T) rather
than transversion mutations (A - T or C - G). This fits what we know about how
these mutations arise in the mitochondrial genome. 

Depending on your analytical question, these thresholds can be adjusted to
identify variants that are more prevalent in other cells. 

## Compute the variant allele frequency for each cell

We currently have information for each strand stored in the mito assay to allow
strand concordance to be assessed. Now that we have our set of high-confidence
informative variants, we can create a new assay containing strand-collapsed
allele frequency counts for each cell for these variants using the `AlleleFreq()`
function.

```{r}
crc <- AlleleFreq(
  object = crc,
  variants = high.conf$variant,
  assay = "mito"
)
crc[["alleles"]]
```

## Visualize the variants

Now that the allele frequencies are stored as an additional assay, we can use
the standard functions in Seurat to visualize how these allele frequencies are
distributed across the cells. Here we visualize a subset of the variants using 
`FeaturePlot()` and `DoHeatmap()`.

```{r visualize variants, message=FALSE, warning = FALSE}
DefaultAssay(crc) <- "alleles"
alleles.view <- c("12889G>A", "16147C>T", "9728C>T", "9804G>A")
FeaturePlot(
  object = crc,
  features = alleles.view,
  order = TRUE,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()
```

```{r message=FALSE, warning=FALSE}
DoHeatmap(crc, features = rownames(crc), slot = "data", disp.max = 1) +
  scale_fill_viridis_c()
```

Here, we can see a few interesting patterns for the selected variants. 16147C>T
is present in essentially all epithelial cells and almost exclusively in
epithelial cells (the edge cases where this isn't true are also cases where the
UMAP and clustering don't full agree). It is at 100% allele frequency-- strongly
suggestive of whatever cell of origin of this tumor had the mutation at 100% and
then expanded. We then see at least 3 variants 1227G>A, 12889G>A, and 9728C>T that
are mostly present specifically in the epithelial cells that define subclones. 
Other variants including 3244G>A, 9804G>A, and 824T>C are found specifically
in immune cell populations, suggesting that these arose from a common hematopoetic
progenitor cell (probably in the bone marrow).

# TF1 cell line dataset

Next we'll demonstrate a similar workflow to identify cell clones in a different
dataset, this time generated from a TF1 cell line. This dataset contains more
clones present at a higher proportion, based on the experimental design.

We'll demonstrate how to identify groups of related cells (clones) by clustering
the allele frequency data and how to relate these clonal groups to 
accessibility differences utilizing the multimodal capabilities of Signac. 

## Data loading

<details>
  <summary>**View data download code**</summary>

To download the data from Zenodo run the following in a shell:

```{bash, eval=FALSE}
# ATAC data
wget https://zenodo.org/record/3977808/files/TF1.filtered.fragments.tsv.gz
wget https://zenodo.org/record/3977808/files/TF1.filtered.fragments.tsv.gz.tbi
wget https://zenodo.org/record/3977808/files/TF1.filtered.narrowPeak.gz

# mitochondrial genome data
wget https://zenodo.org/record/3977808/files/TF1_filtered.A.txt.gz
wget https://zenodo.org/record/3977808/files/TF1_filtered.T.txt.gz
wget https://zenodo.org/record/3977808/files/TF1_filtered.C.txt.gz
wget https://zenodo.org/record/3977808/files/TF1_filtered.G.txt.gz
wget https://zenodo.org/record/3977808/files/TF1_filtered.chrM_refAllele.txt.gz
wget https://zenodo.org/record/3977808/files/TF1_filtered.depthTable.txt.gz
```

</details>

```{r}
# read the mitochondrial data
tf1.data <- ReadMGATK(dir = "../vignette_data/mito/tf1/")

# create a Seurat object
tf1 <- CreateSeuratObject(
  counts = tf1.data$counts,
  meta.data = tf1.data$depth,
  assay = "mito"
)

# load the peak set
peaks <- read.table(
  file = "../vignette_data/mito/TF1.filtered.narrowPeak.gz",
  sep = "\t",
  col.names = c("chrom", "start", "end", "peak", "width", "strand", "x", "y", "z", "w")
)
peaks <- makeGRangesFromDataFrame(peaks)

# create fragment object
frags <- CreateFragmentObject(
  path = "../vignette_data/mito/TF1.filtered.fragments.tsv.gz",
  cells = colnames(tf1)
)

# quantify the DNA accessibility data
counts <- FeatureMatrix(
  fragments = frags,
  features = peaks,
  cells = colnames(tf1)
)

# create assay with accessibility data and add it to the Seurat object
tf1[["peaks"]] <- CreateChromatinAssay(
  counts = counts,
  fragments = frags
)
```

## Quality control

```{r message=FALSE, warning=FALSE}
# add annotations
Annotation(tf1[["peaks"]]) <- annotations
```

```{r message=FALSE, warning=FALSE}
DefaultAssay(tf1) <- "peaks"

tf1 <- NucleosomeSignal(tf1)
tf1 <- TSSEnrichment(tf1)
```

```{r}
VlnPlot(tf1, c("nCount_peaks", "nucleosome_signal", "TSS.enrichment"), pt.size = 0.1)
```
```{r}
tf1 <- subset(
  x = tf1,
  subset = nCount_peaks > 500 &
    nucleosome_signal < 2 &
    TSS.enrichment > 2.5
)
tf1
```

## Identifying variants

```{r}
DefaultAssay(tf1) <- "mito"
variants <- IdentifyVariants(tf1, refallele = tf1.data$refallele)
VariantPlot(variants)
```

```{r}
high.conf <- subset(
  variants, subset = n_cells_conf_detected >= 5 &
    strand_correlation >= 0.65 &
    vmr > 0.01
)
```

```{r}
tf1 <- AlleleFreq(tf1, variants = high.conf$variant, assay = "mito")
tf1[["alleles"]]
```

## Identifying clones

Now that we've identified a set of variable alleles, we can cluster the cells
based on the frequency of each of these alleles using the `FindClonotypes()`
function. This uses the Louvain community detection algorithm implemented in 
Seurat.

```{r message=FALSE, warning=FALSE}
DefaultAssay(tf1) <- "alleles"
tf1 <- FindClonotypes(tf1)
```

```{r}
table(Idents(tf1))
```

Here we see that the clonal clustering has identified 12 different clones in the
TF1 dataset. We can further visualize the frequency of alleles in these clones
using `DoHeatmap()`. The `FindClonotypes()` function also performs hierarchical
clustering on both the clonotypes and the alleles, and sets the factor levels
for the clonotypes based on the hierarchical clustering order, and the order of
variable features based on the hierarchical feature clustering. This allows us
to get a decent ordering of both features and clones automatically:

```{r message=FALSE, warning=FALSE}
DoHeatmap(tf1, features = VariableFeatures(tf1), slot = "data", disp.max = 0.1) +
  scale_fill_viridis_c()
```

## Find differentially accessible peaks between clones

Next we can use the clonal information derived from the mitochondrial assay
to find peaks that are differentially accessible between clones.

```{r message=FALSE, warning=FALSE}
DefaultAssay(tf1) <- "peaks"

# find peaks specific to one clone
markers.fast <- FoldChange(tf1, ident.1 = 2)
head(markers.fast)
```

We can the DNA accessibility in these regions for each clone using
the `CoveragePlot()` function. As you can see, the peaks identified are highly
specific to one clone.

```{r}
CoveragePlot(
  object = tf1,
  region = rownames(markers.fast)[1],
  extend.upstream = 2000,
  extend.downstream = 2000
)
```
<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Merging objects"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

In this vignette we demonstrate how to merge multiple Seurat objects containing 
single-cell chromatin data. To demonstrate, we will use four scATAC-seq PBMC
datasets provided by 10x Genomics:

* [500-cell PBMC](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_pbmc_500_nextgem)
* [1k-cell PBMC](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_pbmc_1k_nextgem)
* [5k-cell PBMC](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_pbmc_5k_nextgem)
* [10k-cell PBMC](https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_pbmc_10k_nextgem)

<details>
  <summary>**View data download code**</summary>

To download the required data, run the following lines in a shell:

```{bash, eval=FALSE}
# 500 cell
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_500_nextgem/atac_pbmc_500_nextgem_filtered_peak_bc_matrix.h5

# 1k cell
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_filtered_peak_bc_matrix.h5

# 5k cell
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5

# 10k cell
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_fragments.tsv.gz.tbi
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_peaks.bed
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_singlecell.csv
wget https://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_pbmc_10k_nextgem/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5
```

</details>

When merging multiple single-cell chromatin datasets, it's important to be aware
that if peak calling was performed on each dataset independently, the peaks are 
unlikely to be exactly the same. We therefore need to create a common set of
peaks across all the datasets to be merged.

To create a unified set of peaks we can use functions from the
[GenomicRanges](https://bioconductor.org/packages/GenomicRanges/) package.
The `reduce` function from GenomicRanges will merge all intersecting peaks. 
Another option is to use the `disjoin` function, that will
create distinct non-overlapping sets of peaks. Here is a visual example to 
illustrate the difference between `reduce` and `disjoin`:

```{r, echo=FALSE, include=FALSE}
library(GenomicRanges)
library(ggplot2)
library(ggbio)
library(patchwork)
```

```{r}
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(20, 70, 300), end = c(120, 200, 400)))
gr
```

```{r fig.height=1, fig.width=6, echo=FALSE}
ggplot(gr) + geom_segment(size = 5) + theme_classic() + ggtitle("Ranges") + ylab("")
ggplot(reduce(gr)) + geom_segment(size = 5) + theme_classic() + ggtitle("Reduce") + ylab("")
ggplot(disjoin(gr)) + geom_segment(size = 5) + theme_classic() + ggtitle("Disjoin") + ylab("")
```

## Creating a common peak set

If the peaks were identified independently in each experiment then they will 
likely not overlap perfectly. We can merge peaks from all the datasets to create
a common peak set, and quantify this peak set in each experiment prior to merging
the objects.

First we'll load the peak coordinates for each experiment and convert them to
genomic ranges, the use the `GenomicRanges::reduce` function to create a common
set of peaks to quantify in each dataset.

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
```

```{r}
# read in peak sets
peaks.500 <- read.table(
  file = "../vignette_data/pbmc500/atac_pbmc_500_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.1k <- read.table(
  file = "../vignette_data/pbmc1k/atac_pbmc_1k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.5k <- read.table(
  file = "../vignette_data/pbmc5k/atac_pbmc_5k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.10k <- read.table(
  file = "../vignette_data/pbmc10k/atac_pbmc_10k_nextgem_peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.500 <- makeGRangesFromDataFrame(peaks.500)
gr.1k <- makeGRangesFromDataFrame(peaks.1k)
gr.5k <- makeGRangesFromDataFrame(peaks.5k)
gr.10k <- makeGRangesFromDataFrame(peaks.10k)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.500, gr.1k, gr.5k, gr.10k))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks
```

## Create Fragment objects

To quantify our combined set of peaks we'll need to create a Fragment object for
each experiment. The Fragment class is a specialized class defined in Signac to 
hold all the information related to a single fragment file.

First we'll load the cell metadata for each experiment so that we know what cell
barcodes are contained in each file, then we can create Fragment objects using
the `CreateFragmentObject` function. The `CreateFragmentObject` function 
performs some checks to ensure that the file is present on disk and that it is
compressed and indexed, computes the MD5 sum for the file and the tabix index so
that we can tell if the file is modified at any point, and checks that the 
expected cells are present in the file.

```{r}
# load metadata
md.500 <- read.table(
  file = "../vignette_data/pbmc500/atac_pbmc_500_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

md.1k <- read.table(
  file = "../vignette_data/pbmc1k/atac_pbmc_1k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.5k <- read.table(
  file = "../vignette_data/pbmc5k/atac_pbmc_5k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.10k <- read.table(
  file = "../vignette_data/pbmc10k/atac_pbmc_10k_nextgem_singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.500 <- md.500[md.500$passed_filters > 500, ]
md.1k <- md.1k[md.1k$passed_filters > 500, ]
md.5k <- md.5k[md.5k$passed_filters > 500, ]
md.10k <- md.10k[md.10k$passed_filters > 1000, ] # sequenced deeper so set higher cutoff

# create fragment objects
frags.500 <- CreateFragmentObject(
  path = "../vignette_data/pbmc500/atac_pbmc_500_nextgem_fragments.tsv.gz",
  cells = rownames(md.500)
)

frags.1k <- CreateFragmentObject(
  path = "../vignette_data/pbmc1k/atac_pbmc_1k_nextgem_fragments.tsv.gz",
  cells = rownames(md.1k)
)

frags.5k <- CreateFragmentObject(
  path = "../vignette_data/pbmc5k/atac_pbmc_5k_nextgem_fragments.tsv.gz",
  cells = rownames(md.5k)
)

frags.10k <- CreateFragmentObject(
  path = "../vignette_data/pbmc10k/atac_pbmc_10k_nextgem_fragments.tsv.gz",
  cells = rownames(md.10k)
)
```

## Quantify peaks in each dataset

We can now create a matrix of peaks x cell for each sample using the
`FeatureMatrix` function. This function is parallelized using the
[`future`](https://cran.r-project.org/package=future) package. See the 
[parallelization](https://satijalab.org/signac/articles/future.html) vignette
for more information about using `future`.

```{r message=FALSE, warning=FALSE, cache=TRUE}
pbmc500.counts <- FeatureMatrix(
  fragments = frags.500,
  features = combined.peaks,
  cells = rownames(md.500)
)

pbmc1k.counts <- FeatureMatrix(
  fragments = frags.1k,
  features = combined.peaks,
  cells = rownames(md.1k)
)

pbmc5k.counts <- FeatureMatrix(
  fragments = frags.5k,
  features = combined.peaks,
  cells = rownames(md.5k)
)

pbmc10k.counts <- FeatureMatrix(
  fragments = frags.10k,
  features = combined.peaks,
  cells = rownames(md.10k)
)
```

## Create the objects

We will now use the quantified matrices to create a Seurat object for each 
dataset, storing the Fragment object for each dataset in the assay.

```{r message=FALSE, warning=FALSE}
pbmc500_assay <- CreateChromatinAssay(pbmc500.counts, fragments = frags.500)
pbmc500 <- CreateSeuratObject(pbmc500_assay, assay = "ATAC", meta.data=md.500)

pbmc1k_assay <- CreateChromatinAssay(pbmc1k.counts, fragments = frags.1k)
pbmc1k <- CreateSeuratObject(pbmc1k_assay, assay = "ATAC", meta.data=md.1k)

pbmc5k_assay <- CreateChromatinAssay(pbmc5k.counts, fragments = frags.5k)
pbmc5k <- CreateSeuratObject(pbmc5k_assay, assay = "ATAC", meta.data=md.5k)

pbmc10k_assay <- CreateChromatinAssay(pbmc10k.counts, fragments = frags.10k)
pbmc10k <- CreateSeuratObject(pbmc10k_assay, assay = "ATAC", meta.data=md.10k)
```

## Merge objects

Now that the objects each contain an assay with the same set of features, we can
use the standard `merge` function to merge the objects. This will also merge all
the fragment objects so that we retain the fragment information for each cell in
the final merged object.

```{r}
# add information to identify dataset of origin
pbmc500$dataset <- 'pbmc500'
pbmc1k$dataset <- 'pbmc1k'
pbmc5k$dataset <- 'pbmc5k'
pbmc10k$dataset <- 'pbmc10k'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = pbmc500,
  y = list(pbmc1k, pbmc5k, pbmc10k),
  add.cell.ids = c("500", "1k", "5k", "10k")
)
combined[["ATAC"]]
```

```{r message=FALSE, warning=FALSE}
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
```

```{r}
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
```

The merged object contains all four fragment objects, and contains an internal 
mapping of cell names in the object to the cell names in each fragment file so
that we can retrieve information from the files without having to change the 
cell names in each fragment file. We can check that functions that pull data 
from the fragment files work as expected on the merged object by plotting a 
region of the genome:

```{r coverageplot, message=FALSE, warning=FALSE, cache=FALSE, out.width='90%', fig.height=4}
CoveragePlot(
  object = combined,
  group.by = 'dataset',
  region = "chr14-99700000-99760000"
)
```

```{r echo=FALSE, message=FALSE, warning=FALSE}
saveRDS(object = combined, file = "../vignette_data/pbmc_combined.rds")
```

## Merging without a common feature set

The above approach requires that we have access to a fragment file for each
dataset. In some cases we may not have this data (although we can create a 
fragment file from the BAM file using 
[sinto](https://timoast.github.io/sinto/basic_usage.html#create-scatac-seq-fragments-file)).
In these cases, we can still create a merged object, with the caveat that the
resulting merged count matrix may not be as accurate.

The `merge` function defined in Signac for `ChromatinAssay` objects will 
consider overlapping peaks as equivalent, and adjust the genomic ranges spanned
by the peak so that the features in each object being merged become equivalent.
**Note that this can result in inaccuracies in the count matrix, as some peaks 
will be extended to cover regions that were not originally quantified**. This is
the best that can be done without re-quantification, and we recommend always 
following the procedure outlined above for object merging whenever possible.

Here we demonstrate merging the same four PBMC datasets without creating a 
common feature set:

```{r message=FALSE, warning=FALSE}
# load the count matrix for each object that was generated by cellranger
counts.500 <- Read10X_h5("../vignette_data/pbmc500/atac_pbmc_500_nextgem_filtered_peak_bc_matrix.h5")
counts.1k <- Read10X_h5("../vignette_data/pbmc1k/atac_pbmc_1k_nextgem_filtered_peak_bc_matrix.h5")
counts.5k <- Read10X_h5("../vignette_data/pbmc5k/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5")
counts.10k <- Read10X_h5("../vignette_data/pbmc10k/atac_pbmc_10k_nextgem_filtered_peak_bc_matrix.h5")

# create objects
pbmc500_assay <- CreateChromatinAssay(counts = counts.500, sep = c(":", "-"), min.features = 500)
pbmc500 <- CreateSeuratObject(pbmc500_assay, assay = "peaks")
pbmc1k_assay <- CreateChromatinAssay(counts = counts.1k, sep = c(":", "-"), min.features = 500)
pbmc1k <- CreateSeuratObject(pbmc1k_assay, assay = "peaks")
pbmc5k_assay <- CreateChromatinAssay(counts = counts.5k, sep = c(":", "-"), min.features = 500)
pbmc5k <- CreateSeuratObject(pbmc5k_assay, assay = "peaks")
pbmc10k_assay <- CreateChromatinAssay(counts = counts.10k, sep = c(":", "-"), min.features = 1000)
pbmc10k <- CreateSeuratObject(pbmc10k_assay, assay = "peaks")

# add information to identify dataset of origin
pbmc500$dataset <- 'pbmc500'
pbmc1k$dataset <- 'pbmc1k'
pbmc5k$dataset <- 'pbmc5k'
pbmc10k$dataset <- 'pbmc10k'

# merge
combined <- merge(
  x = pbmc500,
  y = list(pbmc1k, pbmc5k, pbmc10k),
  add.cell.ids = c("500", "1k", "5k", "10k")
)

# process 
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')
```

```{r}
DimPlot(combined, group.by = 'dataset', pt.size = 0.1)
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Analyzing adult mouse brain scATAC-seq"
output: html_document
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
---
  
```{r init, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

For this tutorial, we will be analyzing a single-cell ATAC-seq dataset of adult
mouse brain cells provided by 10x Genomics. The following files are used in this
vignette, all available through the 10x Genomics website.

<details>
  <summary>**View data download code**</summary>

To download the required data, run the following lines in a shell:

```{bash, eval=FALSE}
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_singlecell.csv
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz
wget http://cf.10xgenomics.com/samples/cell-atac/1.1.0/atac_v1_adult_brain_fresh_5k/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz.tbi
```

</details>

This vignette echoes the commands run in the introductory Signac vignette on
[human PBMC](https://satijalab.org/signac/articles/pbmc_vignette.html). We
provide the same analysis in a different system to demonstrate performance and
applicability to other tissue types, and to provide an example from another
species.

First load in Signac, Seurat, and some other packages we will be using for
analyzing mouse data.

```{r include=FALSE}
if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
    BiocManager::install("EnsDb.Mmusculus.v79")

if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
    BiocManager::install("GenomeInfoDb")
```

```{r setup, message=FALSE}
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
set.seed(1234)
```

## Pre-processing workflow
  
```{r}
counts <- Read10X_h5("../vignette_data/atac_v1_adult_brain_fresh_5k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../vignette_data/atac_v1_adult_brain_fresh_5k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

brain_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = '../vignette_data/atac_v1_adult_brain_fresh_5k_fragments.tsv.gz',
  min.cells = 1
)
brain <- CreateSeuratObject(
  counts = brain_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)
```

We can also add gene annotations to the `brain` object for the mouse genome.
This will allow downstream functions to pull the gene annotation information
directly from the object.

```{r message=FALSE, warning=FALSE}
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(brain) <- annotations
```

## Computing QC Metrics

Next we compute some useful per-cell QC metrics.

```{r message=FALSE, warning=FALSE}
brain <- NucleosomeSignal(object = brain)
```

We can look at the fragment length periodicity for all the cells, and group by
cells with high or low nucleosomal signal strength. You can see that cells which
are outliers for the  mononucleosomal/ nucleosome-free ratio have different
banding patterns. The remaining cells exhibit a pattern that is typical for a
successful ATAC-seq experiment.

```{r message=FALSE, warning=FALSE}
brain$nucleosome_group <- ifelse(brain$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = brain, group.by = 'nucleosome_group', region = 'chr1-1-10000000')
```

The enrichment of Tn5 integration events at transcriptional start sites (TSSs)
can also be an important quality control metric to assess the targeting of Tn5
in ATAC-seq experiments. The ENCODE consortium defined a TSS enrichment score as
the number of Tn5 integration site around the TSS normalized to the number of
Tn5 integration sites in flanking regions. See the ENCODE documentation for more
information about the TSS enrichment score
(https://www.encodeproject.org/data-standards/terms/). We can calculate the TSS
enrichment score for each cell using the `TSSEnrichment()` function in Signac.

```{r message=FALSE, warning=FALSE}
brain <- TSSEnrichment(brain, fast = FALSE)
```

```{r message=FALSE, warning=FALSE}
brain$high.tss <- ifelse(brain$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(brain, group.by = 'high.tss') + NoLegend()
```

```{r message=FALSE, warning=FALSE, fig.width=18, fig.height=6}
brain$pct_reads_in_peaks <- brain$peak_region_fragments / brain$passed_filters * 100
brain$blacklist_ratio <- brain$blacklist_region_fragments / brain$peak_region_fragments

VlnPlot(
  object = brain,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
```

We remove cells that are outliers for these QC metrics.

```{r}
brain <- subset(
  x = brain,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 100000 &
    pct_reads_in_peaks > 40 &
    blacklist_ratio < 0.025 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
brain
```

## Normalization and linear dimensional reduction

```{r message=FALSE, warning=FALSE}
brain <- RunTFIDF(brain)
brain <- FindTopFeatures(brain, min.cutoff = 'q0')
brain <- RunSVD(object = brain)
```

The first LSI component often captures sequencing depth (technical variation)
rather than biological variation. If this is the case, the component should be
removed from downstream analysis. We can assess the correlation between each LSI
component and sequencing depth using the `DepthCor()` function:

```{r}
DepthCor(brain)
```

Here we see there is a very strong correlation between the first LSI component
and the total number of counts for the cell, so we will perform downstream steps
without this component.

## Non-linear dimension reduction and clustering

Now that the cells are embedded in a low-dimensional space, we can use methods
commonly applied for the analysis of scRNA-seq data to perform graph-based
clustering, and non-linear dimension reduction for visualization. The functions
`RunUMAP()`, `FindNeighbors()`, and `FindClusters()` all come from the Seurat
package.

```{r message=FALSE, warning=FALSE}
brain <- RunUMAP(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindNeighbors(
  object = brain,
  reduction = 'lsi',
  dims = 2:30
)
brain <- FindClusters(
  object = brain,
  algorithm = 3,
  resolution = 1.2,
  verbose = FALSE
)

DimPlot(object = brain, label = TRUE) + NoLegend()
```

## Create a gene activity matrix

```{r, message=FALSE, warning=FALSE}
# compute gene activities
gene.activities <- GeneActivity(brain)

# add the gene activity matrix to the Seurat object as a new assay
brain[['RNA']] <- CreateAssayObject(counts = gene.activities)
brain <- NormalizeData(
  object = brain,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(brain$nCount_RNA)
)
```

```{r, fig.width=12, fig.height=10}
DefaultAssay(brain) <- 'RNA'
FeaturePlot(
  object = brain,
  features = c('Sst','Pvalb',"Gad2","Neurod6","Rorb","Syt6"),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
```

## Integrating with scRNA-seq data

To help interpret the scATAC-seq data, we can classify cells based on an
scRNA-seq experiment from the same biological system (the adult mouse brain).
We utilize methods for cross-modality integration and label transfer, described
[here](https://doi.org/10.1016/j.cell.2019.05.031), with a more in-depth
tutorial [here](https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html).

You can download the raw data for this experiment from the Allen Institute
[website](http://celltypes.brain-map.org/api/v2/well_known_file_download/694413985),
and view the code used to construct this object on
[GitHub](https://github.com/satijalab/Integration2019/blob/master/preprocessing_scripts/allen_brain.R). 
Alternatively, you can download the pre-processed Seurat object
[here](https://www.dropbox.com/s/kqsy9tvsklbu7c4/allen_brain.rds).

```{r warning=FALSE, message=FALSE}
# Load the pre-processed scRNA-seq data
allen_rna <- readRDS("../vignette_data/allen_brain.rds")
allen_rna <- FindVariableFeatures(
  object = allen_rna,
  nfeatures = 5000
)

transfer.anchors <- FindTransferAnchors(
  reference = allen_rna,
  query = brain,
  reduction = 'cca',
  dims = 1:40
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_rna$subclass,
  weight.reduction = brain[['lsi']],
  dims = 2:30
)

brain <- AddMetaData(object = brain, metadata = predicted.labels)
```


```{r fig.width=12}
plot1 <- DimPlot(allen_rna, group.by = 'subclass', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot2 <- DimPlot(brain, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
plot1 + plot2
```

<details>
  <summary>**Why did we change default parameters?**</summary>

We changed default parameters for `FindIntegrationAnchors()` and
`FindVariableFeatures()` (including more features and dimensions). You can run
the analysis both ways, and observe very similar results. However, when using
default parameters we mislabel cluster 11 cells as Vip-interneurons, when they
are in fact a Meis2 expressing CGE-derived interneuron population recently
described by [us](https://www.nature.com/articles/nature25999) and
[others](https://www.nature.com/articles/ncomms14219). The reason is that this
subset is exceptionally rare in the scRNA-seq data (0.3%), and so the genes
define this subset (for example, *Meis2*) were too lowly expressed to be
selected in the initial set of variable features. We therefore need more genes
and dimensions to facilitate cross-modality mapping. Interestingly, this subset
is 10-fold more abundant in the scATAC-seq data compared to the scRNA-seq data
(see [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0209648)
for possible explanations.)

</details>

You can see that the RNA-based classifications are entirely consistent with the
UMAP visualization, computed only on the ATAC-seq data. We can now easily
annotate our scATAC-seq derived clusters (alternately, we could use the RNA
classifications themselves). We note three small clusters (13, 20, 21) which
represent subdivisions of the scRNA-seq labels. Try transferring the cluster
label (which shows finer distinctions) from the allen scRNA-seq dataset,
to annotate them!

```{r}
# replace each label with its most likely prediction
for(i in levels(brain)) {
  cells_to_reid <- WhichCells(brain, idents = i)
  newid <- names(sort(table(brain$predicted.id[cells_to_reid]),decreasing=TRUE))[1]
  Idents(brain, cells = cells_to_reid) <- newid
}
```

## Find differentially accessible peaks between clusters

Here, we find differentially accessible regions between excitatory neurons in
different layers of the cortex.

```{r message=TRUE, warning=FALSE}
#switch back to working with peaks instead of gene activities
DefaultAssay(brain) <- 'peaks'

da_peaks <- FindMarkers(
  object = brain,
  ident.1 = c("L2/3 IT"), 
  ident.2 = c("L4", "L5 IT", "L6 IT"),
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)
```

```{r fig.width=12}
plot1 <- VlnPlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("L4","L5 IT","L2/3 IT")
)
plot2 <- FeaturePlot(
  object = brain,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  max.cutoff = 'q95'
)
plot1 | plot2
```

```{r, warning=FALSE, message=FALSE}
open_l23 <- rownames(da_peaks[da_peaks$avg_log2FC > 0.25, ])
open_l456 <- rownames(da_peaks[da_peaks$avg_log2FC < -0.25, ])
closest_l23 <- ClosestFeature(brain, open_l23)
closest_l456 <- ClosestFeature(brain, open_l456)
head(closest_l23)
```

```{r, warning=FALSE, message=FALSE}
head(closest_l456)
```

## Plotting genomic regions

We can also create coverage plots grouped by cluster, cell type, or any other
metadata stored in the object for any genomic region using the `CoveragePlot()`
function. These represent pseudo-bulk accessibility tracks, where signal from
all cells within a group have been averaged together to visualize the DNA 
accessibility in a region.

```{r message=FALSE, warning=FALSE, out.width='90%', fig.height=10}
# set plotting order
levels(brain) <- c("L2/3 IT","L4","L5 IT","L5 PT","L6 CT", "L6 IT","NP","Sst","Pvalb","Vip","Lamp5","Meis2","Oligo","Astro","Endo","VLMC","Macrophage")

CoveragePlot(
  object = brain,
  region = c("Neurod6", "Gad2"),
  extend.upstream = 1000,
  extend.downstream = 1000,
  ncol = 1
)
```

```{r message=FALSE, warning=FALSE, echo=FALSE}
saveRDS(object = brain, file = "../vignette_data/adult_mouse_brain.rds")
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Transcription factor footprinting"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
```

## Data loading

For this vignette we'll use the dataset introduced and pre-processed in the
[trajectory building vignette](monocle.html).

```{r}
bone <- readRDS("../vignette_data/cd34.rds")
DimPlot(bone, label = TRUE)
```

To perform a footprinting analysis we first need to add motif information to the
object, including the exact positions of each motif. This can be done using 
functions from the \code{motifmatchr} and \code{TFBSTools} packages.

```{r include=FALSE}
if (!requireNamespace("JASPAR2020", quietly = TRUE))
    BiocManager::install("JASPAR2020")

if (!requireNamespace("TFBSTools", quietly = TRUE))
    BiocManager::install("TFBSTools")

if (!requireNamespace("motifmatchr", quietly = TRUE))
    BiocManager::install("motifmatchr")

if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE))
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
```

```{r message=FALSE, warning=FALSE}
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)

# extract position frequency matrices for the motifs
pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
bone <- AddMotifs(bone, genome = BSgenome.Hsapiens.UCSC.hg19, pfm = pwm)
```

## Motif footprinting

Now we can footprint any motif that we have positional information for. By
default, this includes every instance of the motif in the genome. We can instead
use the `in.peaks = TRUE` parameter to include only those motifs that fall
inside a peak in the assay. The `Footprint()` function gathers all the required
data and stores it in the assay. We can then plot the footprinted motifs using
the `PlotFootprint()` function.

```{r message=FALSE, warning=FALSE}
# gather the footprinting information for sets of motifs
bone <- Footprint(
  object = bone,
  motif.name = c("GATA2", "CEBPA", "EBF1"),
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(bone, features = c("GATA2", "CEBPA", "EBF1"))
```

```{r fig.height=12, message=FALSE, warning=FALSE}
p2 + patchwork::plot_layout(ncol = 1)
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>
---
title: "Motif analysis with Signac"
date: 'Compiled: `r format(Sys.Date(), "%B %d, %Y")`'
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

In this tutorial, we will perform DNA sequence motif analysis in Signac. We will
explore two complementary options for performing motif analysis: one by finding
overrepresented motifs in a set of differentially accessible peaks, one method
performing differential motif activity analysis between groups of cells.

In this demonstration we use data from the adult mouse brain. See our
[vignette](mouse_brain_vignette.html) for the code used to generate this object,
and links to the raw data. First, load the required packages and the
pre-computed Seurat object:

```{r include=FALSE}
if (!requireNamespace("TFBSTools", quietly = TRUE))
    BiocManager::install("TFBSTools")

if (!requireNamespace("JASPAR2020", quietly = TRUE))
    BiocManager::install("JASPAR2020")

if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10", quietly = TRUE))
    BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")

if (!requireNamespace("chromVAR", quietly = TRUE))
    BiocManager::install("chromVAR")

if (!requireNamespace("ggbio", quietly = TRUE))
    BiocManager::install("ggbio")
```

```{r message=FALSE, warning=FALSE}
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
set.seed(1234)
```

```{r message=FALSE, warning=FALSE}
mouse_brain <- readRDS("../vignette_data/adult_mouse_brain.rds")
mouse_brain
```

```{r message=FALSE, warning=FALSE}
p1 <- DimPlot(mouse_brain, label = TRUE, pt.size = 0.1) + NoLegend()
p1
```

## Adding motif information to the Seurat object

To add the DNA sequence motif information required for motif analyses, we can
run the `AddMotifs()` function:

```{r message=FALSE, warning=FALSE}
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
```

To facilitate motif analysis in Signac, we have create the `Motif` class to
store all the required information, including a list of position weight matrices
(PWMs) or position frequency matrices (PFMs) and a motif occurrence matrix.
Here, the `AddMotifs()` function construct a `Motif` object and adds it to our
mouse brain dataset, along with other information such as the base composition
of each peak. A motif object can be added to any Seurat assay using the `SetAssayData()`
function. See the [object interaction vignette](data_structures.html) for
more information.

## Finding overrepresented motifs

To identify potentially important cell-type-specific regulatory sequences, we
can search for DNA motifs that are overrepresented in a set of peaks that are
differentially accessible between cell types.

Here, we find differentially accessible peaks between Pvalb and Sst inhibitory
interneurons. For sparse data (such as scATAC-seq), we find it is often necessary
to lower the `min.pct` threshold in `FindMarkers()` from the default (0.1, which
was designed for scRNA-seq data). 

We then perform a hypergeometric test to test the probability of
observing the motif at the given frequency by chance, comparing with a
background set of peaks matched for GC content. 

```{r message=FALSE, warning=FALSE}
da_peaks <- FindMarkers(
  object = mouse_brain,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
```

<details>
  <summary>**Optional: choosing a set of background peaks**</summary>
  
Matching the set of background peaks is essential when finding enriched DNA
sequence motifs. By default, we choose a set of peaks matched for GC content,
but it can be sometimes be beneficial to further restrict the background peaks
to those that are accessible in the groups of cells compared when finding
differentially accessible peaks.

The `AccessiblePeaks()` function can be used to find a set of peaks that are 
open in a subset of cells. We can use this function to first restrict the set
of possible background peaks to those peaks that were open in the set of cells
compared in `FindMarkers()`, and then create a GC-content-matched set of peaks
from this larger set using `MatchRegionStats()`.

```{r}
# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(mouse_brain, idents = c("Pvalb", "Sst"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(mouse_brain, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)
```

`peaks.matched` can then be used as the background peak set by setting
`background=peaks.matched` in `FindMotifs()`.

</details>

```{r}
# test enrichment
enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak
)
```

```{r echo=FALSE}
knitr::kable(head(enriched.motifs))
```

We can also plot the position weight matrices for the motifs, so we can
visualize the different motif sequences.

```{r}
MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)
```

We and others have previously shown that Mef-family motifs, particularly
*Mef2c*, are enriched in Pvalb-specific peaks in scATAC-seq data
(https://doi.org/10.1016/j.cell.2019.05.031; https://doi.org/10.1101/615179),
and further shown that *Mef2c* is required for the development of Pvalb
interneurons (https://www.nature.com/articles/nature25999). Here our results are 
consistent with these findings, and we observe a strong enrichment of Mef-family
motifs in the top results from `FindMotifs()`.

## Computing motif activities

We can also compute a per-cell motif activity score by running 
[chromVAR](https://greenleaflab.github.io/chromVAR/index.html). This allows us
to visualize motif activities per cell, and also provides an alternative method
of identifying differentially-active motifs between cell types.

ChromVAR identifies motifs associated with variability in chromatin
accessibility between cells. See the chromVAR
[paper](https://www.nature.com/articles/nmeth.4401) for a complete description
of the method.

```{r message=FALSE, warning=FALSE, fig.width=12}
mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(mouse_brain) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = mouse_brain,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2
```

We can also directly test for differential activity scores between cell types.
This tends to give similar results as performing an enrichment test on
differentially accessible peaks between the cell types (shown above).

When performing differential testing on the chromVAR z-score, we can set 
`mean.fxn=rowMeans` and `fc.name="avg_diff"` in the `FindMarkers()` function so
that the fold-change calculation computes the average difference in z-score
between the groups.

```{r message=FALSE, warning=FALSE}
differential.activity <- FindMarkers(
  object = mouse_brain,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
```

<details>
  <summary>**Session Info**</summary>
  
```{r}
sessionInfo()
```

</details>

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/footprinting.R
\name{InsertionBias}
\alias{InsertionBias}
\alias{InsertionBias.ChromatinAssay}
\alias{InsertionBias.Seurat}
\title{Compute Tn5 insertion bias}
\usage{
InsertionBias(object, ...)

\method{InsertionBias}{ChromatinAssay}(object, genome, region = "chr1-1-249250621", verbose = TRUE, ...)

\method{InsertionBias}{Seurat}(
  object,
  genome,
  assay = NULL,
  region = "chr1-1-249250621",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat or ChromatinAssay object}

\item{...}{Additional arguments passed to \code{\link{StringToGRanges}}}

\item{genome}{A BSgenome object}

\item{region}{Region to use when assessing bias. Default is human chromosome 1.}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use}
}
\value{
Returns a Seurat object
}
\description{
Counts the Tn5 insertion frequency for each DNA hexamer.
}
\examples{
\dontrun{
library(BSgenome.Mmusculus.UCSC.mm10)

region.use <- GRanges(
  seqnames = c('chr1', 'chr2'),
  IRanges(start = c(1,1), end = c(195471971, 182113224))
)

InsertionBias(
 object = atac_small,
 genome = BSgenome.Mmusculus.UCSC.mm10,
 region = region.use
)
}
}
\concept{footprinting}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{PlotFootprint}
\alias{PlotFootprint}
\title{Plot motif footprinting results}
\usage{
PlotFootprint(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  label = TRUE,
  repel = TRUE,
  show.expected = TRUE,
  normalization = "subtract",
  label.top = 3,
  label.idents = NULL
)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A vector of features to plot}

\item{assay}{Name of assay to use}

\item{group.by}{A grouping variable}

\item{idents}{Set of identities to include in the plot}

\item{label}{TRUE/FALSE value to control whether groups are labeled.}

\item{repel}{Repel labels from each other}

\item{show.expected}{Plot the expected Tn5 integration frequency below the
main footprint plot}

\item{normalization}{Method to normalize for Tn5 DNA sequence bias. Options
are "subtract", "divide", or NULL to perform no bias correction.}

\item{label.top}{Number of groups to label based on highest accessibility
in motif flanking region.}

\item{label.idents}{Vector of identities to label. If supplied,
\code{label.top} will be ignored.}
}
\description{
Plot motif footprinting results
}
\concept{footprinting}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/dimension_reduction.R
\name{RunSVD}
\alias{RunSVD}
\alias{RunSVD.default}
\alias{RunSVD.Assay}
\alias{RunSVD.Seurat}
\title{Run singular value decomposition}
\usage{
RunSVD(object, ...)

\method{RunSVD}{default}(
  object,
  assay = NULL,
  n = 50,
  scale.embeddings = TRUE,
  reduction.key = "LSI_",
  scale.max = NULL,
  verbose = TRUE,
  irlba.work = n * 3,
  ...
)

\method{RunSVD}{Assay}(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "LSI_",
  scale.max = NULL,
  verbose = TRUE,
  ...
)

\method{RunSVD}{Seurat}(
  object,
  assay = NULL,
  features = NULL,
  n = 50,
  reduction.key = "LSI_",
  reduction.name = "lsi",
  scale.max = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{assay}{Which assay to use. If NULL, use the default assay}

\item{n}{Number of singular values to compute}

\item{scale.embeddings}{Scale cell embeddings within each component to
mean 0 and SD 1 (default TRUE).}

\item{reduction.key}{Key for dimension reduction object}

\item{scale.max}{Clipping value for cell embeddings.
Default (NULL) is no clipping.}

\item{verbose}{Print messages}

\item{irlba.work}{work parameter for \code{\link[irlba]{irlba}}.
Working subspace dimension, larger values can speed convergence at the
cost of more memory use.}

\item{features}{Which features to use. If NULL, use variable features}

\item{reduction.name}{Name for stored dimension reduction object.
Default 'svd'}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Run partial singular value decomposition using \code{\link[irlba]{irlba}}
}
\examples{
x <- matrix(data = rnorm(100), ncol = 10)
RunSVD(x)
RunSVD(atac_small[['peaks']])
RunSVD(atac_small)
}
\concept{dimension_reduction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{FindTopFeatures}
\alias{FindTopFeatures}
\alias{FindTopFeatures.default}
\alias{FindTopFeatures.Assay}
\alias{FindTopFeatures.Seurat}
\title{Find most frequently observed features}
\usage{
FindTopFeatures(object, ...)

\method{FindTopFeatures}{default}(object, assay = NULL, min.cutoff = "q5", verbose = TRUE, ...)

\method{FindTopFeatures}{Assay}(object, assay = NULL, min.cutoff = "q5", verbose = TRUE, ...)

\method{FindTopFeatures}{Seurat}(object, assay = NULL, min.cutoff = "q5", verbose = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{assay}{Name of assay to use}

\item{min.cutoff}{Cutoff for feature to be included in the VariableFeatures
for the object. This can be a percentile specified as 'q' followed by the
minimum percentile, for example 'q5' to set the top 95\% most common features
as the VariableFeatures for the object. Alternatively, this can be an integer
specifying the minimum number of cells containing the feature for the feature
to be included in the set of VariableFeatures. For example, setting to 10
will include features in >10 cells in the set of VariableFeatures. If NULL,
include all features in VariableFeatures. If NA, VariableFeatures will not be
altered, and only the feature metadata will be updated with the total counts
and percentile rank for each feature.}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Find top features for a given assay based on total number of counts for the
feature. Can specify a minimum cell count, or a lower percentile
bound to determine the set of variable features. Running this function will
store the total counts and percentile rank for each feature in the feature
metadata for the assay. To only compute the feature metadata, without
changing the variable features for the assay, set \code{min.cutoff=NA}.
}
\examples{
FindTopFeatures(object = atac_small[['peaks']][])
FindTopFeatures(object = atac_small[['peaks']])
FindTopFeatures(atac_small)
}
\concept{preprocessing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{CellsPerGroup}
\alias{CellsPerGroup}
\title{Cells per group}
\usage{
CellsPerGroup(object, group.by = NULL)
}
\arguments{
\item{object}{A Seurat object}

\item{group.by}{A grouping variable. Default is the active identities}
}
\value{
Returns a vector
}
\description{
Count the number of cells in each group
}
\examples{
CellsPerGroup(atac_small)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_hg19}
\alias{blacklist_hg19}
\title{Genomic blacklist regions for Human hg19 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_hg19
}
\description{
Genomic blacklist regions for Human hg19 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GRangesToString}
\alias{GRangesToString}
\title{GRanges to String}
\usage{
GRangesToString(grange, sep = c("-", "-"))
}
\arguments{
\item{grange}{A GRanges object}

\item{sep}{Vector of separators to use for genomic string. First element is
used to separate chromosome and coordinates, second separator is used to
separate start and end coordinates.}
}
\value{
Returns a character vector
}
\description{
Convert GRanges object to a vector of strings
}
\examples{
GRangesToString(grange = blacklist_hg19)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{FilterCells}
\alias{FilterCells}
\title{Filter cells from fragment file}
\usage{
FilterCells(
  fragments,
  cells,
  outfile = NULL,
  buffer_length = 256L,
  verbose = TRUE
)
}
\arguments{
\item{fragments}{Path to a fragment file}

\item{cells}{A vector of cells to keep}

\item{outfile}{Name for output file}

\item{buffer_length}{Size of buffer to be read from the fragment file. This
must be longer than the longest line in the file.}

\item{verbose}{Display messages}
}
\description{
Remove all fragments that are not from an allowed set of cell barcodes from
the fragment file. This will create a new file on disk that only contains
fragments from cells specified in the \code{cells} argument. The output file
is block gzip-compressed and indexed, ready for use with Signac functions.
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
tmpf <- tempfile(fileext = ".gz")
FilterCells(
  fragments = fpath,
  cells = head(colnames(atac_small)),
  outfile = tmpf
)
file.remove(tmpf)
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{as.ChromatinAssay}
\alias{as.ChromatinAssay}
\alias{as.ChromatinAssay.Assay}
\title{Convert objects to a ChromatinAssay}
\usage{
as.ChromatinAssay(x, ...)

\method{as.ChromatinAssay}{Assay}(
  x,
  ranges = NULL,
  seqinfo = NULL,
  annotation = NULL,
  motifs = NULL,
  fragments = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  ...
)
}
\arguments{
\item{x}{An object to convert to class \code{\link{ChromatinAssay}}}

\item{...}{Arguments passed to other methods}

\item{ranges}{A GRanges object}

\item{seqinfo}{A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
information about the genome used. Alternatively, the name of a UCSC genome
can be provided and the sequence information will be downloaded from UCSC.}

\item{annotation}{Genomic annotation}

\item{motifs}{A \code{\link{Motif}} object}

\item{fragments}{A list of \code{\link{Fragment}} objects}

\item{bias}{Tn5 integration bias matrix}

\item{positionEnrichment}{A named list of position enrichment matrices.}

\item{sep}{Characters used to separate the chromosome, start, and end
coordinates in the row names of the data matrix}
}
\description{
Convert objects to a ChromatinAssay
}
\concept{assay}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_dm6}
\alias{blacklist_dm6}
\title{Genomic blacklist regions for Drosophila dm6 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_dm6
}
\description{
Genomic blacklist regions for Drosophila dm6 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{TilePlot}
\alias{TilePlot}
\title{Plot integration sites per cell}
\usage{
TilePlot(
  object,
  region,
  sep = c("-", "-"),
  tile.size = 100,
  tile.cells = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  assay = NULL,
  cells = NULL,
  group.by = NULL,
  order.by = "total",
  idents = NULL
)
}
\arguments{
\item{object}{A Seurat object}

\item{region}{A set of genomic coordinates to show. Can be a GRanges object,
a string encoding a genomic position, a gene name, or a vector of strings
describing the genomic coordinates or gene names to plot. If a gene name is
supplied, annotations must be present in the assay.}

\item{sep}{Separators to use for strings encoding genomic coordinates. First
element is used to separate the chromosome from the coordinates, second
element is used to separate the start from end coordinate.}

\item{tile.size}{Size of the sliding window for per-cell fragment tile plot}

\item{tile.cells}{Number of cells to display fragment information for in tile
plot.}

\item{extend.upstream}{Number of bases to extend the region upstream.}

\item{extend.downstream}{Number of bases to extend the region downstream.}

\item{assay}{Name of assay to use}

\item{cells}{Which cells to plot. Default all cells}

\item{group.by}{Name of grouping variable to group cells by. If NULL, use the
current cell identities}

\item{order.by}{Option for determining how cells are chosen from each group.
Options are "total" or "random". "total" will select the top cells based on
total number of fragments in the region, "random" will select randomly.}

\item{idents}{List of cell identities to include in the plot. If NULL, use
all identities.}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Plots the presence/absence of Tn5 integration sites for each cell
within a genomic region.
}
\examples{
\donttest{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
Fragments(atac_small) <- fragments
TilePlot(object = atac_small, region = c("chr1-713500-714500"))
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_dm3}
\alias{blacklist_dm3}
\title{Genomic blacklist regions for Drosophila dm3 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_dm3
}
\description{
Genomic blacklist regions for Drosophila dm3 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_ce10}
\alias{blacklist_ce10}
\title{Genomic blacklist regions for C. elegans ce10 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_ce10
}
\description{
Genomic blacklist regions for C. elegans ce10 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{GetMotifData}
\alias{GetMotifData}
\alias{GetMotifData.Motif}
\alias{GetMotifData.ChromatinAssay}
\alias{GetMotifData.Seurat}
\title{Retrieve a motif matrix}
\usage{
GetMotifData(object, ...)

\method{GetMotifData}{Motif}(object, slot = "data", ...)

\method{GetMotifData}{ChromatinAssay}(object, slot = "data", ...)

\method{GetMotifData}{Seurat}(object, assay = NULL, slot = "data", ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{slot}{Information to pull from object (data, pwm, meta.data)}

\item{assay}{Which assay to use. Default is the current active assay}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Get motif matrix for given assay
}
\examples{
motif.obj <- Seurat::GetAssayData(
  object = atac_small[['peaks']], slot = "motifs"
)
GetMotifData(object = motif.obj)
GetMotifData(object = atac_small)
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/granges-methods.R
\name{granges-methods}
\alias{granges-methods}
\alias{granges}
\alias{granges,ChromatinAssay-method}
\alias{granges,Seurat-method}
\title{Access genomic ranges for ChromatinAssay objects}
\usage{
\S4method{granges}{ChromatinAssay}(x, use.names = TRUE, use.mcols = FALSE, ...)

\S4method{granges}{Seurat}(x, use.names = TRUE, use.mcols = FALSE, ...)
}
\arguments{
\item{x}{A \code{\link{ChromatinAssay}} object}

\item{use.names}{Whether the names on the genomic ranges should be
propagated to the returned object.}

\item{use.mcols}{Not supported for \code{\link{ChromatinAssay}} objects}

\item{...}{Additional arguments}
}
\value{
Returns a \code{\link[GenomicRanges]{GRanges}} object
}
\description{
Methods for accessing \code{\link[GenomicRanges]{GRanges}} object
information stored in a \code{\link{ChromatinAssay}} object.
}
\section{Functions}{
\itemize{
\item \code{granges,Seurat-method}: method for Seurat objects
}}

\examples{
granges(atac_small)
}
\seealso{
\itemize{
  \item{\link[GenomicRanges]{granges} in the \pkg{GenomicRanges} package.}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{granges}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{AddChromatinModule}
\alias{AddChromatinModule}
\title{Add chromatin module}
\usage{
AddChromatinModule(object, features, genome, assay = NULL, verbose = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A named list of features to include in each module. The name
of each element in the list will be used to name the modules computed, which
will be stored in the object metadata.}

\item{genome}{A BSgenome object}

\item{assay}{Name of assay to use. If NULL, use the default assay.}

\item{verbose}{Display messages}

\item{...}{Additional arguments passed to \code{RunChromVAR}}
}
\value{
Returns a Seurat object
}
\description{
Compute chromVAR deviations for groups of peaks. The goal of this function is
similar to that of \code{\link[Seurat]{AddModuleScore}} except that it is
designed for single-cell chromatin data. The chromVAR deviations for each
group of peaks will be added to the object metadata.
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{Annotation}
\alias{Annotation}
\alias{Annotation<-}
\alias{Annotation.ChromatinAssay}
\alias{Annotation.Seurat}
\alias{Annotation<-.ChromatinAssay}
\alias{Annotation<-.Seurat}
\title{Annotation}
\usage{
Annotation(object, ...)

Annotation(object, ...) <- value

\method{Annotation}{ChromatinAssay}(object, ...)

\method{Annotation}{Seurat}(object, ...)

\method{Annotation}{ChromatinAssay}(object, ...) <- value

\method{Annotation}{Seurat}(object, ...) <- value
}
\arguments{
\item{object}{A Seurat object or ChromatinAssay object}

\item{...}{Arguments passed to other methods}

\item{value}{A value to set. Can be NULL, to remove the current annotation
information, or a \code{\link[GenomicRanges]{GRanges}} object. If a
\code{GRanges} object is supplied and the genome information is stored in the
assay, the genome of the new annotations must match the genome of the assay.}
}
\value{
Returns a \code{\link[GenomicRanges]{GRanges}} object
if the annotation data is present, otherwise returns NULL
}
\description{
Get the annotation from a ChromatinAssay
}
\examples{
\donttest{
Annotation(atac_small[["peaks"]])
}
\donttest{
Annotation(atac_small)
}
genes <- Annotation(atac_small)
Annotation(atac_small[["peaks"]]) <- genes
genes <- Annotation(atac_small)
Annotation(atac_small) <- genes
}
\concept{assay}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_hg38_unified}
\alias{blacklist_hg38_unified}
\title{Unified genomic blacklist regions for Human GRCh38}
\format{
A GRanges object
}
\source{
\url{https://www.encodeproject.org/files/ENCFF356LFX/}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_hg38_unified
}
\description{
Manually curated genomic blacklist regions for the hg38 genome by Anshul
Kundaje and Anna Shcherbina. See
\url{https://www.encodeproject.org/files/ENCFF356LFX/} for a description of
how this blacklist was curated.
}
\author{
Anshul Kundaje

Anna Shcherbina
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{PeakPlot}
\alias{PeakPlot}
\title{Plot peaks in a genomic region}
\usage{
PeakPlot(
  object,
  region,
  assay = NULL,
  peaks = NULL,
  group.by = NULL,
  color = "dimgrey"
)
}
\arguments{
\item{object}{A \code{\link[SeuratObject]{Seurat}} object}

\item{region}{A genomic region to plot}

\item{assay}{Name of assay to use. If NULL, use the default assay.}

\item{peaks}{A GRanges object containing peak coordinates. If NULL, use
coordinates stored in the Seurat object.}

\item{group.by}{Name of variable in feature metadata (if using ranges in the
Seurat object) or genomic ranges metadata (if using supplied ranges) to color
ranges by. If NULL, do not color by any metadata variable.}

\item{color}{Color to use. If \code{group.by} is not NULL, this can be a
custom color scale (see examples).}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Display the genomic ranges in a \code{\link{ChromatinAssay}} object that fall
in a given genomic region
}
\examples{
\donttest{
# plot peaks in assay
PeakPlot(atac_small, region = "chr1-710000-715000")

# manually set color
PeakPlot(atac_small, region = "chr1-710000-715000", color = "red")

# color by a variable in the feature metadata
PeakPlot(atac_small, region = "chr1-710000-715000", group.by = "count")
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{findOverlaps-methods}
\alias{findOverlaps-methods}
\alias{findOverlaps}
\alias{findOverlaps,Vector,ChromatinAssay-method}
\alias{findOverlaps,ChromatinAssay,Vector-method}
\alias{findOverlaps,ChromatinAssay,ChromatinAssay-method}
\alias{findOverlaps,Vector,Seurat-method}
\alias{findOverlaps,Seurat,Vector-method}
\alias{findOverlaps,Seurat,Seurat-method}
\alias{countOverlaps,Vector,ChromatinAssay-method}
\alias{countOverlaps}
\alias{countOverlaps,ChromatinAssay,Vector-method}
\alias{countOverlaps,ChromatinAssay,ChromatinAssay-method}
\alias{countOverlaps,Seurat,Vector-method}
\alias{countOverlaps,Vector,Seurat-method}
\alias{countOverlaps,Seurat,Seurat-method}
\title{Find overlapping ranges for ChromatinAssay objects}
\usage{
\S4method{findOverlaps}{Vector,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{ChromatinAssay,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{ChromatinAssay,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Vector,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Seurat,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{findOverlaps}{Seurat,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Vector,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{ChromatinAssay,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{ChromatinAssay,ChromatinAssay}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Seurat,Vector}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Vector,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)

\S4method{countOverlaps}{Seurat,Seurat}(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ignore.strand = FALSE
)
}
\arguments{
\item{query, subject}{A \code{\link{ChromatinAssay}} object}

\item{maxgap, minoverlap, type, select, ignore.strand}{See
\code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} and
\pkg{IRanges} packages.}
}
\value{
See \code{\link[GenomicRanges]{findOverlaps}}
}
\description{
The \code{findOverlaps, countOverlaps} methods are available for
\code{\link{ChromatinAssay}} objects. This allows finding overlaps between
genomic ranges and the ranges stored in the ChromatinAssay.
}
\details{
If a ChromatinAssay is set as the default assay in a
\code{\link[SeuratObject]{Seurat}} object, you can also call \code{findOverlaps}
directly on the Seurat object.
}
\section{Functions}{
\itemize{
\item \code{findOverlaps,ChromatinAssay,Vector-method}: method for ChromatinAssay, Vector

\item \code{findOverlaps,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{findOverlaps,Vector,Seurat-method}: method for Vector, Seurat

\item \code{findOverlaps,Seurat,Vector-method}: method for Seurat, Vector

\item \code{findOverlaps,Seurat,Seurat-method}: method for Seurat, Seurat

\item \code{countOverlaps,Vector,ChromatinAssay-method}: method for Vector, ChromatinAssay

\item \code{countOverlaps,ChromatinAssay,Vector-method}: method for ChromatinAssay, Vector

\item \code{countOverlaps,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{countOverlaps,Seurat,Vector-method}: method for Seurat, Vector

\item \code{countOverlaps,Vector,Seurat-method}: method for Vector, Seurat

\item \code{countOverlaps,Seurat,Seurat-method}: method for Seurat, Seurat
}}

\seealso{
\itemize{
  \item{\link[IRanges]{findOverlaps-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{findOverlaps-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{overlaps}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links.R
\name{GetLinkedPeaks}
\alias{GetLinkedPeaks}
\title{Get peaks linked to genes}
\usage{
GetLinkedPeaks(object, features, assay = NULL, min.abs.score = 0)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A list of genes to find linked peaks for}

\item{assay}{Name of assay to use. If NULL, use the default assay}

\item{min.abs.score}{Minimum absolute value of the link score for a link to
be returned}
}
\description{
Find peaks linked to a given set of genes
}
\seealso{
GetLinkedGenes
}
\concept{links}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links.R
\name{GetLinkedGenes}
\alias{GetLinkedGenes}
\title{Get genes linked to peaks}
\usage{
GetLinkedGenes(object, features, assay = NULL, min.abs.score = 0)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A list of peaks to find linked genes for}

\item{assay}{Name of assay to use. If NULL, use the default assay}

\item{min.abs.score}{Minimum absolute value of the link score for a link to
be returned}
}
\description{
Find genes linked to a given set of peaks
}
\seealso{
GetLinkedPeaks
}
\concept{links}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{atac_small}
\alias{atac_small}
\title{A small example scATAC-seq dataset}
\format{
A Seurat object with the following assays
\describe{
  \item{peaks}{A peak x cell dataset}
  \item{bins}{A 5 kb genome bin x cell dataset}
  \item{RNA}{A gene x cell dataset}
}
}
\source{
\url{https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_pbmc_10k}
}
\usage{
atac_small
}
\description{
A subsetted version of 10x Genomics 10k human (hg19) PBMC scATAC-seq dataset
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/motifs.R
\name{AddMotifs}
\alias{AddMotifs}
\alias{AddMotifs.default}
\alias{AddMotifs.ChromatinAssay}
\alias{AddMotifs.Seurat}
\title{Add DNA sequence motif information}
\usage{
AddMotifs(object, ...)

\method{AddMotifs}{default}(object, genome, pfm, verbose = TRUE, ...)

\method{AddMotifs}{ChromatinAssay}(object, genome, pfm, verbose = TRUE, ...)

\method{AddMotifs}{Seurat}(object, genome, pfm, assay = NULL, verbose = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object or ChromatinAssay object}

\item{...}{Additional arguments passed to other methods}

\item{genome}{A \code{BSgenome}, \code{DNAStringSet}, \code{FaFile}, or
string stating the genome build recognized by \code{getBSgenome}.}

\item{pfm}{A \code{PFMatrixList} or \code{PWMatrixList} object containing
position weight/frequency matrices to use}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use. If NULL, use the default assay}
}
\value{
When running on a \code{ChromatinAssay} or \code{Seurat} object,
returns a modified version of the input object. When running on a matrix,
returns a \code{Motif} object.
}
\description{
Construct a \code{\link{Motif}} object containing DNA sequence motif
information and add it to an existing Seurat object or ChromatinAssay.
If running on a Seurat object, \code{AddMotifs} will also run
\code{\link{RegionStats}} to compute the GC content of each peak and store
the results in the feature metadata.
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/motifs.R
\name{RunChromVAR}
\alias{RunChromVAR}
\alias{RunChromVAR.ChromatinAssay}
\alias{RunChromVAR.Seurat}
\title{Run chromVAR}
\usage{
RunChromVAR(object, ...)

\method{RunChromVAR}{ChromatinAssay}(object, genome, motif.matrix = NULL, verbose = TRUE, ...)

\method{RunChromVAR}{Seurat}(
  object,
  genome,
  motif.matrix = NULL,
  assay = NULL,
  new.assay.name = "chromvar",
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Additional arguments passed to
\code{\link[chromVAR]{getBackgroundPeaks}}}

\item{genome}{A \code{BSgenome}, \code{DNAStringSet}, \code{FaFile}, or
string stating the genome build recognized by \code{getBSgenome}.}

\item{motif.matrix}{A peak x motif matrix. If NULL, pull the peak x motif
matrix from a Motif object stored in the assay.}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use}

\item{new.assay.name}{Name of new assay used to store the chromVAR results.
Default is "chromvar".}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object with a new assay
}
\description{
Wrapper to run \code{\link[chromVAR]{chromVAR}} on an assay with a motif
object present. Will return a new Seurat assay with the motif activities
(the deviations in chromatin accessibility across the set of regions) as
a new assay.
}
\details{
See the chromVAR documentation for more information:
\url{https://greenleaflab.github.io/chromVAR/index.html}

See the chromVAR paper: \url{https://www.nature.com/articles/nmeth.4401}
}
\examples{
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
RunChromVAR(object = atac_small[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg19)
}
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
RunChromVAR(object = atac_small, genome = BSgenome.Hsapiens.UCSC.hg19)
}
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{BinarizeCounts}
\alias{BinarizeCounts}
\alias{BinarizeCounts.default}
\alias{BinarizeCounts.Assay}
\alias{BinarizeCounts.Seurat}
\title{Binarize counts}
\usage{
BinarizeCounts(object, ...)

\method{BinarizeCounts}{default}(object, assay = NULL, verbose = TRUE, ...)

\method{BinarizeCounts}{Assay}(object, assay = NULL, verbose = TRUE, ...)

\method{BinarizeCounts}{Seurat}(object, assay = NULL, verbose = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{assay}{Name of assay to use. Can be a list of assays,
and binarization will be applied to each.}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Set counts >1 to 1 in a count matrix
}
\examples{
x <- matrix(data = sample(0:3, size = 25, replace = TRUE), ncol = 5)
BinarizeCounts(x)
BinarizeCounts(atac_small[['peaks']])
BinarizeCounts(atac_small)
}
\concept{preprocessing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{CountFragments}
\alias{CountFragments}
\title{Count fragments}
\usage{
CountFragments(fragments, cells = NULL, max_lines = NULL, verbose = TRUE)
}
\arguments{
\item{fragments}{Path to a fragment file}

\item{cells}{Cells to include. If NULL, include all cells}

\item{max_lines}{Maximum number of lines to read from the fragment file. If
NULL, read all lines in the file.}

\item{verbose}{Display messages}
}
\value{
Returns a data.frame
}
\description{
Count total fragments per cell barcode present in a fragment file.
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
counts <- CountFragments(fragments = fpath)
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{AverageCounts}
\alias{AverageCounts}
\title{Average Counts}
\usage{
AverageCounts(object, assay = NULL, group.by = NULL, verbose = TRUE)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use. Default is the active assay}

\item{group.by}{Grouping variable to use. Default is the active identities}

\item{verbose}{Display messages}
}
\value{
Returns a dataframe
}
\description{
Compute the mean counts per group of cells for a given assay
}
\examples{
AverageCounts(atac_small)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{SubsetMatrix}
\alias{SubsetMatrix}
\title{Subset matrix rows and columns}
\usage{
SubsetMatrix(
  mat,
  min.rows = 1,
  min.cols = 1,
  max.row.val = 10,
  max.col.val = NULL
)
}
\arguments{
\item{mat}{A matrix}

\item{min.rows}{Minimum number of non-zero elements for
the row to be retained}

\item{min.cols}{Minimum number of non-zero elements for
the column to be retained}

\item{max.row.val}{Maximum allowed value in a row for the
row to be retained. If NULL, don't set any limit.}

\item{max.col.val}{Maximum allowed value in a column for
the column to be retained. If NULL, don't set any limit.}
}
\value{
Returns a matrix
}
\description{
Subset the rows and columns of a matrix by removing
rows and columns with less than the specified number of
non-zero elements.
}
\examples{
SubsetMatrix(mat = volcano)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{StringToGRanges}
\alias{StringToGRanges}
\title{String to GRanges}
\usage{
StringToGRanges(regions, sep = c("-", "-"), ...)
}
\arguments{
\item{regions}{Vector of genomic region strings}

\item{sep}{Vector of separators to use for genomic string. First element is
used to separate chromosome and coordinates, second separator is used to
separate start and end coordinates.}

\item{...}{Additional arguments passed to
\code{\link[GenomicRanges]{makeGRangesFromDataFrame}}}
}
\value{
Returns a GRanges object
}
\description{
Convert a genomic coordinate string to a GRanges object
}
\examples{
regions <- c('chr1-1-10', 'chr2-12-3121')
StringToGRanges(regions = regions)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantification.R
\name{GenomeBinMatrix}
\alias{GenomeBinMatrix}
\title{Genome bin matrix}
\usage{
GenomeBinMatrix(
  fragments,
  genome,
  cells = NULL,
  binsize = 5000,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
}
\arguments{
\item{fragments}{Path to tabix-indexed fragments file or a list of
\code{\link{Fragment}} objects}

\item{genome}{A vector of chromosome sizes for the genome. This is used to
construct the genome bin coordinates. The can be obtained by calling
\code{\link[GenomeInfoDb]{seqlengths}} on a
\code{\link[BSgenome]{BSgenome-class}} object.}

\item{cells}{Vector of cells to include. If NULL, include all cells found
in the fragments file}

\item{binsize}{Size of the genome bins to use}

\item{process_n}{Number of regions to load into memory at a time, per thread.
Processing more regions at once can be faster but uses more memory.}

\item{sep}{Vector of separators to use for genomic string. First element is
used to separate chromosome and coordinates, second separator is used to
separate start and end coordinates.}

\item{verbose}{Display messages}
}
\value{
Returns a sparse matrix
}
\description{
Construct a bin x cell matrix from a fragments file.
}
\details{
This function bins the genome and calls \code{\link{FeatureMatrix}} to
construct a bin x cell matrix.
}
\examples{
\donttest{
genome <- 780007
names(genome) <- 'chr1'
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(fpath)
GenomeBinMatrix(
  fragments = fragments,
  genome = genome,
  binsize = 1000
)
}
}
\concept{quantification}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GetCellsInRegion}
\alias{GetCellsInRegion}
\title{Get cells in a region}
\usage{
GetCellsInRegion(tabix, region, cells = NULL)
}
\arguments{
\item{tabix}{Tabix object}

\item{region}{A string giving the region to extract from the fragments file}

\item{cells}{Vector of cells to include in output. If NULL, include all cells}
}
\value{
Returns a list
}
\description{
Extract cell names containing reads mapped within a given genomic region
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
GetCellsInRegion(tabix = fpath, region = "chr1-10245-762629")
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{NucleosomeSignal}
\alias{NucleosomeSignal}
\title{NucleosomeSignal}
\usage{
NucleosomeSignal(
  object,
  assay = NULL,
  n = ncol(object) * 5000,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use. Only required if a fragment path is not
provided. If NULL, use the active assay.}

\item{n}{Number of lines to read from the fragment file. If NULL, read all
lines. Default scales with the number of cells in the object.}

\item{verbose}{Display messages}

\item{...}{Arguments passed to other functions}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object with
added metadata for the ratio of mononucleosomal to nucleosome-free fragments
per cell, and the percentile rank of each ratio.
}
\description{
Calculate the strength of the nucleosome signal per cell.
Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to
fragments < 147 bp (nucleosome-free)
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
Fragments(atac_small) <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  tolerance = 0.5
)
NucleosomeSignal(object = atac_small)
}
\concept{qc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footprinting.R
\name{GetFootprintData}
\alias{GetFootprintData}
\title{Get footprinting data}
\usage{
GetFootprintData(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL
)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A vector of features to extract data for}

\item{assay}{Name of assay to use}

\item{group.by}{A grouping variable}

\item{idents}{Set of identities to group cells by}
}
\value{
Returns a matrix
}
\description{
Extract footprint data for a set of transcription factors or metafeatures.
This function will pull accessibility data for a given feature (eg, a TF),
and perform background normalization for each identity class. This is the
data that's used to create TF footprinting plots with the
\code{PlotFootprint} function.
}
\concept{footprinting}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{Motifs}
\alias{Motifs}
\alias{Motifs<-}
\alias{Motifs.ChromatinAssay}
\alias{Motifs.Seurat}
\alias{Motifs<-.ChromatinAssay}
\alias{Motifs<-.Seurat}
\title{Get or set a motif information}
\usage{
Motifs(object, ...)

Motifs(object, ...) <- value

\method{Motifs}{ChromatinAssay}(object, ...)

\method{Motifs}{Seurat}(object, ...)

\method{Motifs}{ChromatinAssay}(object, ...) <- value

\method{Motifs}{Seurat}(object, ...) <- value
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{value}{A \code{\link{Motif}} object}
}
\description{
Get or set the Motif object for a Seurat object or ChromatinAssay.
}
\examples{
Motifs(atac_small[["peaks"]])
Motifs(atac_small)
motifs <- Motifs(atac_small)
Motifs(atac_small[["peaks"]]) <- motifs
motifs <- Motifs(atac_small)
Motifs(atac_small) <- motifs
}
\concept{assay}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{VariantPlot}
\alias{VariantPlot}
\title{Plot strand concordance vs. VMR}
\usage{
VariantPlot(
  variants,
  min.cells = 2,
  concordance.threshold = 0.65,
  vmr.threshold = 0.01
)
}
\arguments{
\item{variants}{A dataframe containing variant information. This should be
computed using \code{\link{IdentifyVariants}}}

\item{min.cells}{Minimum number of high-confidence cells detected with the
variant for the variant to be displayed.}

\item{concordance.threshold}{Strand concordance threshold}

\item{vmr.threshold}{Mean-variance ratio threshold}
}
\description{
Plot the Pearson correlation between allele frequencies on each strand
versus the log10 mean-variance ratio for the allele.
}
\concept{mito}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/peaks.R
\name{CallPeaks}
\alias{CallPeaks}
\alias{CallPeaks.Seurat}
\alias{CallPeaks.ChromatinAssay}
\alias{CallPeaks.Fragment}
\alias{CallPeaks.default}
\title{Call peaks}
\usage{
CallPeaks(object, ...)

\method{CallPeaks}{Seurat}(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  macs2.path = NULL,
  broad = FALSE,
  format = "BED",
  outdir = tempdir(),
  fragment.tempdir = tempdir(),
  combine.peaks = TRUE,
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = Project(object),
  cleanup = TRUE,
  verbose = TRUE,
  ...
)

\method{CallPeaks}{ChromatinAssay}(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
)

\method{CallPeaks}{Fragment}(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
)

\method{CallPeaks}{default}(
  object,
  macs2.path = NULL,
  outdir = tempdir(),
  broad = FALSE,
  format = "BED",
  effective.genome.size = 2.7e+09,
  extsize = 200,
  shift = -extsize/2,
  additional.args = NULL,
  name = "macs2",
  cleanup = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object, ChromatinAssay object, Fragment object, or the
path to fragment file/s.}

\item{...}{Arguments passed to other methods}

\item{assay}{Name of assay to use}

\item{group.by}{Grouping variable to use. If set, peaks will be called
independently on each group of cells and then combined. Note that to call
peaks using subsets of cells we first split the fragment file/s used, so
using a grouping variable will require extra time to split the files and
perform multiple MACS peak calls, and will store additional files on-disk
that may be large. Note that we store split fragment files in the temp
directory (\code{\link[base]{tempdir}}) by default, and if the program is
interrupted before completing these temporary files will not be removed. If
NULL, peaks are called using all cells together (pseudobulk).}

\item{idents}{List of identities to include if grouping cells (only valid if
also setting the \code{group.by} parameter). If NULL, peaks will be called
for all cell identities.}

\item{macs2.path}{Path to MACS program. If NULL, try to find MACS
automatically.}

\item{broad}{Call broad peaks (\code{--broad} parameter for MACS)}

\item{format}{File format to use. Should be either "BED" or "BEDPE" (see 
MACS documentation).}

\item{outdir}{Path for output files}

\item{fragment.tempdir}{Path to write temporary fragment files. Only used if
\code{group.by} is not NULL.}

\item{combine.peaks}{Controls whether peak calls from different groups of
cells are combined using \code{GenomicRanges::reduce} when calling peaks for
different groups of cells (\code{group.by} parameter). If FALSE, a list of
\code{GRanges} object will be returned. Note that metadata fields such as the
p-value, q-value, and fold-change information for each peak will be lost if
combining peaks.}

\item{effective.genome.size}{Effective genome size parameter for MACS
(\code{-g}). Default is the human effective genome size (2.7e9).}

\item{extsize}{\code{extsize} parameter for MACS. Only relevant if 
format="BED"}

\item{shift}{\code{shift} parameter for MACS. Only relevant if format="BED"}

\item{additional.args}{Additional arguments passed to MACS. This should be a
single character string}

\item{name}{Name for output MACS files. This will also be placed in the
\code{name} field in the GRanges output.}

\item{cleanup}{Remove MACS output files}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[GenomicRanges]{GRanges}} object
}
\description{
Call peaks using MACS. Fragment files linked to the specified assay will be
used to call peaks. If multiple fragment files are present, all will be used
in a single MACS invocation. Returns the \code{.narrowPeak} MACS output as a
\code{GRanges} object.
}
\details{
See \url{https://macs3-project.github.io/MACS/} for MACS documentation.

If you call peaks using MACS2 please cite:
\doi{10.1186/gb-2008-9-9-r137}
}
\concept{quantification}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{Cells.Fragment}
\alias{Cells.Fragment}
\alias{Cells<-.Fragment}
\title{Set and get cell barcode information for a \code{\link{Fragment}} object}
\usage{
\method{Cells}{Fragment}(x, ...)

\method{Cells}{Fragment}(x, ...) <- value
}
\arguments{
\item{x}{A Fragment object}

\item{...}{Arguments passed to other methods}

\item{value}{A vector of cell names to store in the \code{\link{Fragment}}
object}
}
\description{
This returns the names of cells in the object that are contained in the
fragment file. These cell barcodes may not match the barcodes present in the
fragment file. The \code{\link{Fragment}} object contains an internal mapping
of the cell names in the \code{\link{ChromatinAssay}} object to the cell
names in the fragment file, so that cell names can be changed in the
assay without needing to change the cell names on disk.
}
\details{
To access the cell names that are stored in the fragment file itself, use
\code{GetFragmentData(object = x, name = "cells")}.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{BigwigTrack}
\alias{BigwigTrack}
\title{Plot data from BigWig}
\usage{
BigwigTrack(
  region,
  bigwig,
  smooth = 200,
  type = "coverage",
  y_label = "Score",
  max.downsample = 3000,
  downsample.rate = 0.1
)
}
\arguments{
\item{region}{GRanges object specifying region to plot}

\item{bigwig}{Path to a bigwig file}

\item{smooth}{Number of bases to smooth data over (rolling mean). If NULL,
do not apply smoothing.}

\item{type}{Plot type. Can be one of "line", "heatmap", or "coverage"}

\item{y_label}{Y-axis label}

\item{max.downsample}{Minimum number of positions kept when downsampling.
Downsampling rate is adaptive to the window size, but this parameter will set
the minimum possible number of positions to include so that plots do not
become too sparse when the window size is small.}

\item{downsample.rate}{Fraction of positions to retain when downsampling.
Retaining more positions can give a higher-resolution plot but can make the
number of points large, resulting in larger file sizes when saving the plot
and a longer period of time needed to draw the plot.}
}
\description{
Create a BigWig track. Note that this function does not work on windows.
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{Jaccard}
\alias{Jaccard}
\title{Calculate the Jaccard index between two matrices}
\usage{
Jaccard(x, y)
}
\arguments{
\item{x}{The first matrix}

\item{y}{The second matrix}
}
\value{
Returns a matrix
}
\description{
Finds the Jaccard similarity between rows of the two matrices. Note that
the matrices must be binary, and any rows with zero total counts will result
in an NaN entry that could cause problems in downstream analyses.
}
\details{
This will calculate the raw Jaccard index, without normalizing for the
expected similarity between cells due to differences in sequencing depth.
}
\examples{
x <- matrix(data = sample(c(0, 1), size = 25, replace = TRUE), ncol = 5)
Jaccard(x = x, y = x)
}
\concept{dimension_reduction}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{FractionCountsInRegion}
\alias{FractionCountsInRegion}
\title{Fraction of counts in a genomic region}
\usage{
FractionCountsInRegion(object, regions, assay = NULL, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{regions}{A GRanges object containing a set of genomic regions}

\item{assay}{Name of assay to use}

\item{...}{Additional arguments passed to \code{\link{CountsInRegion}}}
}
\value{
Returns a numeric vector
}
\description{
Find the fraction of counts per cell that overlap a given set of genomic
ranges
}
\examples{
\dontrun{
FractionCountsInRegion(
  object = atac_small,
  assay = 'bins',
  regions = blacklist_hg19
)
}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{SetMotifData}
\alias{SetMotifData}
\alias{SetMotifData.Motif}
\alias{SetMotifData.ChromatinAssay}
\alias{SetMotifData.Seurat}
\title{Set motif data}
\usage{
SetMotifData(object, ...)

\method{SetMotifData}{Motif}(object, slot, new.data, ...)

\method{SetMotifData}{ChromatinAssay}(object, slot, new.data, ...)

\method{SetMotifData}{Seurat}(object, assay = NULL, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{slot}{Name of slot to use}

\item{new.data}{motif matrix to add. Should be matrix or sparse matrix class}

\item{assay}{Name of assay whose data should be set}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Set motif matrix for given assay
}
\examples{
motif.obj <- Seurat::GetAssayData(
  object = atac_small[['peaks']], slot = "motifs"
)
SetMotifData(object = motif.obj, slot = 'data', new.data = matrix())
SetMotifData(
  object = atac_small[['peaks']], slot = 'data', new.data = matrix()
)
motif.matrix <- GetMotifData(object = atac_small)
SetMotifData(
object = atac_small, assay = 'peaks', slot = 'data', new.data = motif.matrix
)
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{Fragments}
\alias{Fragments}
\alias{Fragments<-}
\alias{Fragments.ChromatinAssay}
\alias{Fragments.Seurat}
\alias{Fragments<-.ChromatinAssay}
\alias{Fragments<-.Seurat}
\title{Get the Fragment objects}
\usage{
Fragments(object, ...)

Fragments(object, ...) <- value

\method{Fragments}{ChromatinAssay}(object, ...)

\method{Fragments}{Seurat}(object, ...)

\method{Fragments}{ChromatinAssay}(object, ...) <- value

\method{Fragments}{Seurat}(object, ...) <- value
}
\arguments{
\item{object}{A Seurat object or ChromatinAssay object}

\item{...}{Arguments passed to other methods}

\item{value}{A \code{\link{Fragment}} object or list of Fragment objects}
}
\value{
Returns a list of \code{\link{Fragment}} objects. If there are
no Fragment objects present, returns an empty list.
}
\description{
Get the Fragment objects
}
\examples{
Fragments(atac_small[["peaks"]])
Fragments(atac_small)
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
Fragments(atac_small[["bins"]]) <- fragments
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
Fragments(atac_small) <- fragments
}
\concept{assay}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{CombineTracks}
\alias{CombineTracks}
\title{Combine genome region plots}
\usage{
CombineTracks(plotlist, expression.plot = NULL, heights = NULL, widths = NULL)
}
\arguments{
\item{plotlist}{A list of plots to combine. Must be from the same genomic
region.}

\item{expression.plot}{Plot containing gene expression information. If
supplied, this will be placed to the left of the coverage tracks and aligned
with each track}

\item{heights}{Relative heights for each plot. If NULL, the first plot will
be 8x the height of the other tracks.}

\item{widths}{Relative widths for each plot. Only required if adding a gene
expression panel. If NULL, main plots will be 8x the width of the gene
expression panel}
}
\value{
Returns a patchworked ggplot2 object
}
\description{
This can be used to combine coverage plots, peak region plots, gene
annotation plots, and linked element plots. The different tracks are stacked
on top of each other and the x-axis combined.
}
\examples{
\donttest{
p1 <- PeakPlot(atac_small, region = "chr1-29554-39554")
p2 <- AnnotationPlot(atac_small, region = "chr1-29554-39554")
CombineTracks(plotlist = list(p1, p2), heights = c(1, 1))
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantification.R
\name{FeatureMatrix}
\alias{FeatureMatrix}
\title{Feature Matrix}
\usage{
FeatureMatrix(
  fragments,
  features,
  cells = NULL,
  process_n = 2000,
  sep = c("-", "-"),
  verbose = TRUE
)
}
\arguments{
\item{fragments}{A list of \code{\link{Fragment}} objects. Note that if
setting the \code{cells} parameter, the requested cells should be present in
the supplied \code{Fragment} objects. However, if the cells information in
the fragment object is not set (\code{Cells(fragments)} is \code{NULL}), then
the fragment object will still be searched.}

\item{features}{A GRanges object containing a set of genomic intervals.
These will form the rows of the matrix, with each entry recording the number
of unique reads falling in the genomic region for each cell.}

\item{cells}{Vector of cells to include. If NULL, include all cells found
in the fragments file}

\item{process_n}{Number of regions to load into memory at a time, per thread.
Processing more regions at once can be faster but uses more memory.}

\item{sep}{Vector of separators to use for genomic string. First element is
used to separate chromosome and coordinates, second separator is used to
separate start and end coordinates.}

\item{verbose}{Display messages}
}
\value{
Returns a sparse matrix
}
\description{
Construct a feature x cell matrix from a genomic fragments file
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(fpath)
FeatureMatrix(
  fragments = fragments,
  features = granges(atac_small)
)
}
\concept{quantification}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{ValidateCells}
\alias{ValidateCells}
\title{Validate cells present in fragment file}
\usage{
ValidateCells(
  object,
  cells = NULL,
  tolerance = 0.5,
  max.lines = 5e+07,
  verbose = TRUE
)
}
\arguments{
\item{object}{A \code{\link{Fragment}} object}

\item{cells}{A character vector containing cell barcodes to search for.
If NULL, use the cells stored in the Fragment object.}

\item{tolerance}{Fraction of input cells that can be unseen before returning
TRUE. For example, \code{tolerance = 0.01} will return TRUE when 99% of cells
have observed fragments in the file. This can be useful if there are cells
present that have much fewer total counts, and would require extensive
searching before a fragment from those cells are found.}

\item{max.lines}{Maximum number of lines to read in without finding the
required number of cells before returning FALSE. Setting this value avoids
having to search the whole file if it becomes clear that the expected cells
are not present. Setting this value to NULL will enable an exhaustive search
of the entire file.}

\item{verbose}{Display messages}
}
\description{
Search for a fragment from each cell that should exist in the fragment file.
Will iterate through chunks of the fragment file until at least one fragment
from each cell barcode requested is found.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mito.R
\name{FindClonotypes}
\alias{FindClonotypes}
\title{Find clonotypes}
\usage{
FindClonotypes(
  object,
  assay = NULL,
  features = NULL,
  metric = "cosine",
  resolution = 1,
  k = 10,
  algorithm = 3
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use}

\item{features}{Features to include when constructing neighbor graph}

\item{metric}{Distance metric to use}

\item{resolution}{Clustering resolution to use. See
\code{\link[Seurat]{FindClusters}}}

\item{k}{Passed to \code{k.param} argument in
\code{\link[Seurat]{FindNeighbors}}}

\item{algorithm}{Community detection algorithm to use. See
\code{\link[Seurat]{FindClusters}}}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Identify groups of related cells from allele frequency data. This will
cluster the cells based on their allele frequencies, reorder the factor
levels for the cluster identities by hierarchical clustering the collapsed
(pseudobulk) cluster allele frequencies, and set the variable features for
the allele frequency assay to the order of features defined by hierarchical
clustering.
}
\concept{mito}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_hg38}
\alias{blacklist_hg38}
\title{Genomic blacklist regions for Human GRCh38}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_hg38
}
\description{
Genomic blacklist regions for Human GRCh38
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mito.R
\name{ClusterClonotypes}
\alias{ClusterClonotypes}
\title{Find relationships between clonotypes}
\usage{
ClusterClonotypes(object, assay = NULL, group.by = NULL)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use}

\item{group.by}{Grouping variable for cells}
}
\value{
Returns a list containing two objects of class
\code{\link[stats]{hclust}}, one for the cell clustering and one for the
feature (allele) clustering
}
\description{
Perform hierarchical clustering on clonotype data
}
\concept{mito}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{ClosestFeature}
\alias{ClosestFeature}
\title{Closest Feature}
\usage{
ClosestFeature(object, regions, annotation = NULL, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{regions}{A set of genomic regions to query}

\item{annotation}{A GRanges object containing annotation information. If
NULL, use the annotations stored in the object.}

\item{...}{Additional arguments passed to \code{\link{StringToGRanges}}}
}
\value{
Returns a dataframe with the name of each region, the closest feature
in the annotation, and the distance to the feature.
}
\description{
Find the closest feature to a given set of genomic regions
}
\examples{
\donttest{
ClosestFeature(
  object = atac_small,
  regions = head(granges(atac_small))
)
}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/motifs.R
\name{ConvertMotifID}
\alias{ConvertMotifID}
\alias{ConvertMotifID.default}
\alias{ConvertMotifID.Motif}
\alias{ConvertMotifID.ChromatinAssay}
\alias{ConvertMotifID.Seurat}
\title{Convert between motif name and motif ID}
\usage{
ConvertMotifID(object, ...)

\method{ConvertMotifID}{default}(object, name, id, ...)

\method{ConvertMotifID}{Motif}(object, ...)

\method{ConvertMotifID}{ChromatinAssay}(object, ...)

\method{ConvertMotifID}{Seurat}(object, assay = NULL, ...)
}
\arguments{
\item{object}{A Seurat, ChromatinAssay, or Motif object}

\item{...}{Arguments passed to other methods}

\item{name}{A vector of motif names}

\item{id}{A vector of motif IDs. Only one of \code{name} and \code{id} should
be supplied}

\item{assay}{For \code{Seurat} object. Name of assay to use.
If NULL, use the default assay}
}
\value{
Returns a character vector with the same length and order as the
input. Any names or IDs that were not found will be stored as \code{NA}.
}
\description{
Converts from motif name to motif ID or vice versa. To convert common names
to IDs, use the \code{name} parameter. To convert IDs to common names, use
the \code{id} parameter.
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{ChromatinAssay-class}
\alias{ChromatinAssay-class}
\alias{ChromatinAssay}
\title{The ChromatinAssay class}
\description{
The ChromatinAssay object is an extended \code{\link[SeuratObject]{Assay}}
for the storage and analysis of single-cell chromatin data.
}
\section{Slots}{

\describe{
\item{\code{ranges}}{A \code{\link[GenomicRanges]{GRanges}} object describing the
genomic location of features in the object}

\item{\code{motifs}}{A \code{\link{Motif}} object}

\item{\code{fragments}}{A list of \code{\link{Fragment}} objects.}

\item{\code{seqinfo}}{A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
information about the genome sequence used.}

\item{\code{annotation}}{A  \code{\link[GenomicRanges]{GRanges}} object containing
genomic annotations}

\item{\code{bias}}{A vector containing Tn5 integration bias information
(frequency of Tn5 integration at different kmers)}

\item{\code{positionEnrichment}}{A named list of matrices containing positional
enrichment scores for Tn5 integration (for example, enrichment at the TSS)}

\item{\code{links}}{A \code{\link[GenomicRanges]{GRanges}} object describing linked
genomic positions, such as co-accessible sites or enhancer-gene regulatory
relationships. This should be a \code{GRanges} object, where the start and
end coordinates are the two linked genomic positions, and must contain a
"score" metadata column.}
}}

\concept{assay}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GetIntersectingFeatures}
\alias{GetIntersectingFeatures}
\title{Find intersecting regions between two objects}
\usage{
GetIntersectingFeatures(
  object.1,
  object.2,
  assay.1 = NULL,
  assay.2 = NULL,
  distance = 0,
  verbose = TRUE
)
}
\arguments{
\item{object.1}{The first Seurat object}

\item{object.2}{The second Seurat object}

\item{assay.1}{Name of the assay to use in the first object. If NULL, use
the default assay}

\item{assay.2}{Name of the assay to use in the second object. If NULL, use
the default assay}

\item{distance}{Maximum distance between regions allowed for an intersection
to be recorded. Default is 0.}

\item{verbose}{Display messages}
}
\value{
Returns a list of two character vectors containing the row names
in each object that overlap each other.
}
\description{
Intersects the regions stored in the rownames of two objects and
returns a vector containing the names of rows that intersect
for each object. The order of the row names return corresponds
to the intersecting regions, i.e. the nth feature of the first vector
will intersect the nth feature in the second vector. A distance
parameter can be given, in which case features within the given
distance will be called as intersecting.
}
\examples{
GetIntersectingFeatures(
  object.1 = atac_small,
  object.2 = atac_small,
  assay.1 = 'peaks',
  assay.2 = 'bins'
)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mito.R
\name{ReadMGATK}
\alias{ReadMGATK}
\title{Read MGATK output}
\usage{
ReadMGATK(dir, verbose = TRUE)
}
\arguments{
\item{dir}{Path to directory containing MGATK output files}

\item{verbose}{Display messages}
}
\value{
Returns a list containing a sparse matrix (counts) and two dataframes
(depth and refallele).

The sparse matrix contains read counts for each base at each position
and strand.

The depth dataframe contains the total depth for each cell.
The refallele dataframe contains the reference genome allele at each
position.
}
\description{
Read output files from MGATK (\url{https://github.com/caleblareau/mgatk}).
}
\examples{
\dontrun{
data.dir <- system.file("extdata", "test_mgatk", package="Signac")
mgatk <- ReadMGATK(dir = data.dir)
}
}
\concept{mito}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{DepthCor}
\alias{DepthCor}
\title{Plot sequencing depth correlation}
\usage{
DepthCor(object, assay = NULL, reduction = "lsi", n = 10, ...)
}
\arguments{
\item{object}{A \code{\link[SeuratObject]{Seurat}} object}

\item{assay}{Name of assay to use for sequencing depth. If NULL, use the
default assay.}

\item{reduction}{Name of a dimension reduction stored in the
input object}

\item{n}{Number of components to use. If \code{NULL}, use all components.}

\item{...}{Additional arguments passed to \code{\link[stats]{cor}}}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Compute the correlation between total counts and each reduced
dimension component.
}
\examples{
\donttest{
DepthCor(object = atac_small)
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/mito.R
\name{AlleleFreq}
\alias{AlleleFreq}
\alias{AlleleFreq.default}
\alias{AlleleFreq.Assay}
\alias{AlleleFreq.Seurat}
\title{Compute allele frequencies per cell}
\usage{
AlleleFreq(object, ...)

\method{AlleleFreq}{default}(object, variants, ...)

\method{AlleleFreq}{Assay}(object, variants, ...)

\method{AlleleFreq}{Seurat}(object, variants, assay = NULL, new.assay.name = "alleles", ...)
}
\arguments{
\item{object}{A Seurat object, Assay, or matrix}

\item{...}{Arguments passed to other methods}

\item{variants}{A character vector of informative variants to keep. For
example, \code{c("627G>A","709G>A","1045G>A","1793G>A")}.}

\item{assay}{Name of assay to use}

\item{new.assay.name}{Name of new assay to store variant data in}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object with a new assay
containing the allele frequencies for the informative variants.
}
\description{
Collapses allele counts for each strand and normalize by the total number of
counts at each nucleotide position.
}
\concept{mito}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{AccessiblePeaks}
\alias{AccessiblePeaks}
\title{Accessible peaks}
\usage{
AccessiblePeaks(
  object,
  assay = NULL,
  idents = NULL,
  cells = NULL,
  min.cells = 10
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use}

\item{idents}{A set of identity classes to find accessible peaks for}

\item{cells}{A vector of cells to find accessible peaks for}

\item{min.cells}{Minimum number of cells with the peak accessible (>0 counts)
for the peak to be called accessible}
}
\value{
Returns a vector of peak names
}
\description{
Find accessible peaks in a set of cells
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{IntersectMatrix}
\alias{IntersectMatrix}
\title{Intersect genomic coordinates with matrix rows}
\usage{
IntersectMatrix(
  matrix,
  regions,
  invert = FALSE,
  sep = c("-", "-"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{matrix}{A matrix with genomic regions in the rows}

\item{regions}{A set of genomic regions to intersect with regions in the
matrix. Either a vector of strings encoding the genomic coordinates, or a
GRanges object.}

\item{invert}{Discard rows intersecting the genomic regions supplied, rather
than retain.}

\item{sep}{A length-2 character vector containing the separators to be used
for extracting genomic coordinates from a string. The first element will be
used to separate the chromosome name from coordinates, and the second element
used to separate start and end coordinates.}

\item{verbose}{Display messages}

\item{...}{Additional arguments passed to \code{\link[IRanges]{findOverlaps}}}
}
\value{
Returns a sparse matrix
}
\description{
Remove or retain matrix rows that intersect given genomic regions
}
\examples{
counts <- matrix(data = rep(0, 12), ncol = 2)
rownames(counts) <- c("chr1-565107-565550","chr1-569174-569639",
"chr1-713460-714823","chr1-752422-753038",
"chr1-762106-763359","chr1-779589-780271")
IntersectMatrix(matrix = counts, regions = blacklist_hg19)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GeneActivity}
\alias{GeneActivity}
\title{Create gene activity matrix}
\usage{
GeneActivity(
  object,
  assay = NULL,
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 0,
  biotypes = "protein_coding",
  max.width = 5e+05,
  gene.id = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use. If NULL, use the default assay}

\item{features}{Genes to include. If NULL, use all protein-coding genes in
the annotations stored in the object}

\item{extend.upstream}{Number of bases to extend upstream of the TSS}

\item{extend.downstream}{Number of bases to extend downstream of the TTS}

\item{biotypes}{Gene biotypes to include. If NULL, use all biotypes in the
gene annotation.}

\item{max.width}{Maximum allowed gene width for a gene to be quantified.
Setting this parameter can avoid quantifying extremely long transcripts that
can add a relatively long amount of time. If NULL, do not filter genes based
on width.}

\item{gene.id}{Record gene IDs in output matrix rather than gene name.}

\item{verbose}{Display messages}

\item{...}{Additional options passed to \code{\link{FeatureMatrix}}}
}
\value{
Returns a sparse matrix
}
\description{
Compute counts per cell in gene body and promoter region.
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
Fragments(atac_small) <- fragments
GeneActivity(atac_small)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{UnifyPeaks}
\alias{UnifyPeaks}
\title{Unify genomic ranges}
\usage{
UnifyPeaks(object.list, mode = "reduce")
}
\arguments{
\item{object.list}{A list of Seurat objects or ChromatinAssay objects}

\item{mode}{Function to use when combining genomic ranges. Can be "reduce"
(default) or "disjoin".
See \code{\link[GenomicRanges]{reduce}}
and \code{\link[GenomicRanges]{disjoin}}
for more information on these functions.}
}
\value{
Returns a GRanges object
}
\description{
Create a unified set of non-overlapping genomic ranges
from multiple Seurat objects containing single-cell
chromatin data.
}
\examples{
UnifyPeaks(object.list = list(atac_small, atac_small))
}
\concept{preprocessing}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{theme_browser}
\alias{theme_browser}
\title{Genome browser theme}
\usage{
theme_browser(..., legend = TRUE)
}
\arguments{
\item{...}{Additional arguments}

\item{legend}{Display plot legend}
}
\description{
Theme applied to panels in the \code{\link{CoveragePlot}} function.
}
\examples{
\donttest{
PeakPlot(atac_small, region = "chr1-710000-715000") + theme_browser()
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{SplitFragments}
\alias{SplitFragments}
\title{Split fragment file by cell identities}
\usage{
SplitFragments(
  object,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  outdir = getwd(),
  file.suffix = "",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use}

\item{group.by}{Name of grouping variable to group cells by}

\item{idents}{List of identities to include}

\item{outdir}{Directory to write output files}

\item{file.suffix}{Suffix to add to all file names (before file extension).
If splitting multiple fragment files without the \code{append} option set to
TRUE, an additional numeric suffix will be added to each file (eg, .1, .2).}

\item{append}{If splitting multiple fragment files, append cells from the
same group (eg cluster) to the same file. Note that this can cause the output
file to be unsorted.}

\item{buffer_length}{Size of buffer to be read from the fragment file. This
must be longer than the longest line in the file.}

\item{verbose}{Display messages}
}
\description{
Splits a fragment file into separate files for each group of cells. If
splitting multiple fragment files containing common cell types, fragments
originating from different files will be appended to the same file for one
group of cell identities.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{DownsampleFeatures}
\alias{DownsampleFeatures}
\title{Downsample Features}
\usage{
DownsampleFeatures(object, assay = NULL, n = 20000, verbose = TRUE)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of assay to use. Default is the active assay.}

\item{n}{Number of features to retain (default 20000).}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object with
\code{\link[SeuratObject]{VariableFeatures}} set to the randomly sampled features.
}
\description{
Randomly downsample features and assign to VariableFeatures for the object.
This will select n features at random.
}
\examples{
DownsampleFeatures(atac_small, n = 10)
}
\concept{preprocessing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{subset.Motif}
\alias{subset.Motif}
\alias{subset}
\alias{[.Motif}
\title{Subset a Motif object}
\usage{
\method{subset}{Motif}(x, features = NULL, motifs = NULL, ...)

\method{[}{Motif}(x, i, j, ...)
}
\arguments{
\item{x}{A Motif object}

\item{features}{Which features to retain}

\item{motifs}{Which motifs to retain}

\item{...}{Arguments passed to other methods}

\item{i}{Which columns to retain}

\item{j}{Which rows to retain}
}
\value{
Returns a subsetted \code{\link{Motif}} object
}
\description{
Returns a subset of a \code{\link{Motif-class}} object.
}
\examples{
motif.obj <- Seurat::GetAssayData(
  object = atac_small[['peaks']], slot = "motifs"
)
subset(x = motif.obj, features = head(rownames(motif.obj), 10))
motif.obj <- Seurat::GetAssayData(
  object = atac_small, assay = 'peaks', slot = 'motifs'
)
motif.obj[1:10,1:10]
}
\seealso{
\code{\link[base]{subset}}
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_mm10}
\alias{blacklist_mm10}
\title{Genomic blacklist regions for Mouse mm10 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_mm10
}
\description{
Genomic blacklist regions for Mouse mm10 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{package}
\name{Signac-package}
\alias{Signac}
\alias{Signac-package}
\title{Signac: Analysis of Single-Cell Chromatin Data}
\description{
A framework for the analysis and exploration of single-cell chromatin data. The 'Signac' package contains functions for quantifying single-cell chromatin data, computing per-cell quality control metrics, dimension reduction and normalization, visualization, and DNA sequence motif analysis. Reference: Stuart et al. (2021) <doi:10.1038/s41592-021-01282-5>.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/timoast/signac}
  \item \url{https://satijalab.org/signac}
  \item Report bugs at \url{https://github.com/timoast/signac/issues}
}

}
\author{
\strong{Maintainer}: Tim Stuart \email{tstuart@nygenome.org} (\href{https://orcid.org/0000-0002-3044-0897}{ORCID})

Authors:
\itemize{
  \item Avi Srivastava \email{asrivastava@nygenome.org} (\href{https://orcid.org/0000-0001-9798-2079}{ORCID})
}

Other contributors:
\itemize{
  \item Paul Hoffman \email{phoffman@nygenome.org} (\href{https://orcid.org/0000-0002-7693-8957}{ORCID}) [contributor]
  \item Rahul Satija \email{rsatija@nygenome.org} (\href{https://orcid.org/0000-0001-9448-8833}{ORCID}) [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/footprinting.R
\name{Footprint}
\alias{Footprint}
\alias{Footprint.ChromatinAssay}
\alias{Footprint.Seurat}
\title{Transcription factor footprinting analysis}
\usage{
Footprint(object, ...)

\method{Footprint}{ChromatinAssay}(
  object,
  genome,
  motif.name = NULL,
  key = motif.name,
  regions = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  compute.expected = TRUE,
  in.peaks = FALSE,
  verbose = TRUE,
  ...
)

\method{Footprint}{Seurat}(
  object,
  genome,
  regions = NULL,
  motif.name = NULL,
  assay = NULL,
  upstream = 250,
  downstream = 250,
  in.peaks = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat or ChromatinAssay object}

\item{...}{Arguments passed to other methods}

\item{genome}{A \code{\link[BSgenome]{BSgenome}} object}

\item{motif.name}{Name of a motif stored in the assay to footprint. If not
supplied, must supply a set of regions.}

\item{key}{Key to store positional enrichment information under.}

\item{regions}{A set of genomic ranges containing the motif instances}

\item{assay}{Name of assay to use}

\item{upstream}{Number of bases to extend upstream}

\item{downstream}{Number of bases to extend downstream}

\item{compute.expected}{Find the expected number of insertions at each
position given the local DNA sequence context and the insertion bias of Tn5}

\item{in.peaks}{Restrict motifs to those that fall in peaks}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Compute the normalized observed/expected Tn5 insertion frequency
for each position surrounding a set of motif instances.
}
\concept{footprinting}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{Extend}
\alias{Extend}
\title{Extend}
\usage{
Extend(x, upstream = 0, downstream = 0, from.midpoint = FALSE)
}
\arguments{
\item{x}{A range}

\item{upstream}{Length to extend upstream}

\item{downstream}{Length to extend downstream}

\item{from.midpoint}{Count bases from region midpoint,
rather than the 5' or 3' end for upstream and downstream
respectively.}
}
\value{
Returns a \code{\link[GenomicRanges]{GRanges}} object
}
\description{
Resize GenomicRanges upstream and or downstream.
From \url{https://support.bioconductor.org/p/78652/}
}
\examples{
Extend(x = blacklist_hg19, upstream = 100, downstream = 100)
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{GetFragmentData}
\alias{GetFragmentData}
\title{Get Fragment object data}
\usage{
GetFragmentData(object, slot = "path")
}
\arguments{
\item{object}{A \code{\link{Fragment}} object}

\item{slot}{Information to pull from object (path, hash, cells, prefix, suffix)}
}
\description{
Extract data from a \code{\link{Fragment-class}} object
}
\concept{assay}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differential_accessibility.R, R/fragments.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{FoldChange}
\alias{Cells}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{Seurat}{\code{\link[Seurat]{FoldChange}}}

  \item{SeuratObject}{\code{\link[SeuratObject]{Cells}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{CoverageBrowser}
\alias{CoverageBrowser}
\title{Genome browser}
\usage{
CoverageBrowser(object, region, assay = NULL, sep = c("-", "-"), ...)
}
\arguments{
\item{object}{A Seurat object}

\item{region}{A set of genomic coordinates}

\item{assay}{Name of assay to use}

\item{sep}{Separators for genomic coordinates if region supplied as a string
rather than GRanges object}

\item{...}{Parameters passed to \code{\link{CoveragePlot}}}
}
\value{
Returns a list of ggplot objects
}
\description{
Interactive version of the \code{\link{CoveragePlot}} function. Allows
altering the genome position interactively. The current view at any time can
be saved to a list of \code{\link[ggplot2]{ggplot}} objects using the "Save
plot" button, and this list of plots will be returned after ending the
browser by pressing the "Done" button.
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{RunTFIDF}
\alias{RunTFIDF}
\alias{RunTFIDF.default}
\alias{RunTFIDF.Assay}
\alias{RunTFIDF.Seurat}
\title{Compute the term-frequency inverse-document-frequency}
\usage{
RunTFIDF(object, ...)

\method{RunTFIDF}{default}(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 10000,
  idf = NULL,
  verbose = TRUE,
  ...
)

\method{RunTFIDF}{Assay}(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 10000,
  idf = NULL,
  verbose = TRUE,
  ...
)

\method{RunTFIDF}{Seurat}(
  object,
  assay = NULL,
  method = 1,
  scale.factor = 10000,
  idf = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{assay}{Name of assay to use}

\item{method}{Which TF-IDF implementation to use. Choice of:
\itemize{
 \item{1}: The TF-IDF implementation used by Stuart & Butler et al. 2019
 (\doi{10.1101/460147}). This computes
 \eqn{\log(TF \times IDF)}.
 \item{2}: The TF-IDF implementation used by Cusanovich & Hill
 et al. 2018 (\doi{10.1016/j.cell.2018.06.052}). This
 computes \eqn{TF \times (\log(IDF))}.
 \item{3}: The log-TF method used by Andrew Hill.
 This computes \eqn{\log(TF) \times \log(IDF)}.
 \item{4}: The 10x Genomics method (no TF normalization). This computes
 \eqn{IDF}.
}}

\item{scale.factor}{Which scale factor to use. Default is 10000.}

\item{idf}{A precomputed IDF vector to use. If NULL, compute based on the
input data matrix.}

\item{verbose}{Print progress}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Run term frequency inverse document frequency (TF-IDF) normalization on a
matrix.
}
\details{
Four different TF-IDF methods are implemented. We recommend using method 1
(the default).
}
\examples{
mat <- matrix(data = rbinom(n = 25, size = 5, prob = 0.2), nrow = 5)
RunTFIDF(object = mat)
RunTFIDF(atac_small[['peaks']])
RunTFIDF(object = atac_small)
}
\references{
\url{https://en.wikipedia.org/wiki/Latent_semantic_analysis#Latent_semantic_indexing}
}
\concept{preprocessing}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{UpdatePath}
\alias{UpdatePath}
\title{Update the file path for a Fragment object}
\usage{
UpdatePath(object, new.path, verbose = TRUE)
}
\arguments{
\item{object}{A \code{\link{Fragment}} object}

\item{new.path}{Path to the fragment file}

\item{verbose}{Display messages}
}
\description{
Change the path to a fragment file store in a \code{\link{Fragment}}
object. Path must be to the same file that was used to create the fragment
object. An MD5 hash will be computed using the new path and compared to the
hash stored in the Fragment object to verify that the files are the same.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{CreateMotifObject}
\alias{CreateMotifObject}
\title{Create motif object}
\usage{
CreateMotifObject(
  data = NULL,
  pwm = NULL,
  motif.names = NULL,
  positions = NULL,
  meta.data = NULL
)
}
\arguments{
\item{data}{A motif x region matrix}

\item{pwm}{A named list of position weight matrices or position frequency
matrices matching the motif names in \code{data}.
Can be of class PFMatrixList.}

\item{motif.names}{A named list of motif names. List element names
must match the names given in \code{pwm}. If NULL, use the names from the
list of position weight or position frequency matrices. This can be used to
set a alternative common name for the motif. If a PFMatrixList is passed to
\code{pwm}, it will pull the motif name from the PFMatrixList.}

\item{positions}{A \code{\link[GenomicRanges]{GRangesList}} object containing
exact positions of each motif.}

\item{meta.data}{A data.frame containing metadata}
}
\value{
Returns a \code{\link{Motif}} object
}
\description{
Create a \code{\link{Motif-class}} object.
}
\examples{
motif.matrix <- matrix(
  data = sample(c(0,1),
    size = 100,
    replace = TRUE),
  ncol = 5
)
motif <- CreateMotifObject(data = motif.matrix)
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{Fragment-class}
\alias{Fragment-class}
\alias{Fragment}
\title{The Fragment class}
\description{
The Fragment class is designed to hold information needed for working with
fragment files.
}
\section{Slots}{

\describe{
\item{\code{path}}{Path to the fragment file on disk.
See \url{https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments}}

\item{\code{hash}}{A vector of two md5sums: first element is the md5sum of the
fragment file, the second element is the md5sum of the index.}

\item{\code{cells}}{A named vector of cells where each element is the cell barcode
as it appears in the fragment file, and the name of each element is the
corresponding cell barcode as stored in the ChromatinAssay object.}
}}

\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links.R
\name{LinkPeaks}
\alias{LinkPeaks}
\title{Link peaks to genes}
\usage{
LinkPeaks(
  object,
  peak.assay,
  expression.assay,
  expression.slot = "data",
  gene.coords = NULL,
  distance = 5e+05,
  min.distance = NULL,
  min.cells = 10,
  method = "pearson",
  genes.use = NULL,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{object}{A Seurat object}

\item{peak.assay}{Name of assay containing peak information}

\item{expression.assay}{Name of assay containing gene expression information}

\item{expression.slot}{Name of slot to pull expression data from}

\item{gene.coords}{GRanges object containing coordinates of genes in the
expression assay. If NULL, extract from gene annotations stored in the assay.}

\item{distance}{Distance threshold for peaks to include in regression model}

\item{min.distance}{Minimum distance between peak and TSS to include in
regression model. If NULL (default), no minimum distance is used.}

\item{min.cells}{Minimum number of cells positive for the peak and gene
needed to include in the results.}

\item{method}{Which correlation coefficient to compute. Can be "pearson"
(default), "spearman", or "kendall".}

\item{genes.use}{Genes to test. If NULL, determine from expression assay.}

\item{n_sample}{Number of peaks to sample at random when computing the null
distribution.}

\item{pvalue_cutoff}{Minimum p-value required to retain a link. Links with a
p-value equal or greater than this value will be removed from the output.}

\item{score_cutoff}{Minimum absolute value correlation coefficient for a link
to be retained}

\item{gene.id}{Set to TRUE if genes in the expression assay are named
using gene IDs rather than gene names.}

\item{verbose}{Display messages}
}
\value{
Returns a Seurat object with the \code{Links} information set. This is
a \code{\link[GenomicRanges]{granges}} object accessible via the \code{\link{Links}}
function, with the following information:
\itemize{
  \item{score: the correlation coefficient between the accessibility of the
  peak and expression of the gene}
  \item{zscore: the z-score of the correlation coefficient, computed based on
  the distribution of correlation coefficients from a set of background peaks}
  \item{pvalue: the p-value associated with the z-score for the link}
  \item{gene: name of the linked gene}
  \item{peak: name of the linked peak}
}
}
\description{
Find peaks that are correlated with the expression of nearby genes.
For each gene, this function computes the correlation coefficient between
the gene expression and accessibility of each peak within a given distance
from the gene TSS, and computes an expected correlation coefficient for each
peak given the GC content, accessibility, and length of the peak. The expected
coefficient values for the peak are then used to compute a z-score and p-value.
}
\details{
This function was inspired by the method originally described by SHARE-seq
(Sai Ma et al. 2020, Cell). Please consider citing the original SHARE-seq
work if using this function: \doi{10.1016/j.cell.2020.09.056}
}
\concept{links}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{inter-range-methods}
\alias{inter-range-methods}
\alias{range}
\alias{range,ChromatinAssay-method}
\alias{range,Seurat-method}
\alias{reduce,ChromatinAssay-method}
\alias{reduce}
\alias{reduce,Seurat-method}
\alias{gaps,ChromatinAssay-method}
\alias{gaps}
\alias{gaps,Seurat-method}
\alias{disjoin,ChromatinAssay-method}
\alias{disjoin}
\alias{disjoin,Seurat-method}
\alias{isDisjoint,ChromatinAssay-method}
\alias{isDisjoint}
\alias{isDisjoint,Seurat-method}
\alias{disjointBins,ChromatinAssay-method}
\alias{disjointBins}
\alias{disjointBins,Seurat-method}
\title{Inter-range transformations for ChromatinAssay objects}
\usage{
\S4method{range}{ChromatinAssay}(x, ..., with.revmap = FALSE, na.rm = FALSE)

\S4method{range}{Seurat}(x, ..., with.revmap = FALSE, na.rm = FALSE)

\S4method{reduce}{ChromatinAssay}(x, drop.empty.ranges = FALSE, ...)

\S4method{reduce}{Seurat}(x, drop.empty.ranges = FALSE, ...)

\S4method{gaps}{ChromatinAssay}(x, start = NA, end = NA)

\S4method{gaps}{Seurat}(x, start = NA, end = NA)

\S4method{disjoin}{ChromatinAssay}(x, ...)

\S4method{disjoin}{Seurat}(x, ...)

\S4method{isDisjoint}{ChromatinAssay}(x, ...)

\S4method{isDisjoint}{Seurat}(x, ...)

\S4method{disjointBins}{ChromatinAssay}(x, ...)

\S4method{disjointBins}{Seurat}(x, ...)
}
\arguments{
\item{x}{A \code{\link{ChromatinAssay}} object}

\item{...}{Additional arguments}

\item{with.revmap}{See \code{\link[IRanges]{inter-range-methods}} in the
\pkg{IRanges} packages}

\item{na.rm}{Ignored}

\item{drop.empty.ranges}{See \code{?\link{IRanges}{inter-range-methods}}}

\item{start, end}{See \code{?\link{IRanges}{inter-range-methods}}}
}
\description{
The \code{range, reduce, gaps, disjoin, isDisjoint, disjointBins} methods
are available for \code{\link{ChromatinAssay}} objects.
}
\section{Functions}{
\itemize{
\item \code{range,Seurat-method}: method for Seurat objects

\item \code{reduce,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{reduce,Seurat-method}: method for Seurat objects

\item \code{gaps,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{gaps,Seurat-method}: method for Seurat objects

\item \code{disjoin,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{disjoin,Seurat-method}: method for Seurat objects

\item \code{isDisjoint,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{isDisjoint,Seurat-method}: method for Seurat objects

\item \code{disjointBins,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{disjointBins,Seurat-method}: method for Seurat objects
}}

\seealso{
\itemize{
  \item{\link[IRanges]{inter-range-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{inter-range-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{inter_range}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{MotifPlot}
\alias{MotifPlot}
\title{Plot DNA sequence motif}
\usage{
MotifPlot(object, motifs, assay = NULL, use.names = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{motifs}{A list of motifs to plot}

\item{assay}{Name of the assay to use}

\item{use.names}{Use motif names stored in the motif object}

\item{...}{Additional parameters passed to \code{\link[ggseqlogo]{ggseqlogo}}}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot position weight matrix or position frequency matrix for different DNA
sequence motifs.
}
\examples{
\donttest{
motif.obj <- Seurat::GetAssayData(atac_small, slot = "motifs")
MotifPlot(atac_small, motifs = head(colnames(motif.obj)))
}
}
\concept{motifs}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{FRiP}
\alias{FRiP}
\title{Calculate fraction of reads in peaks per cell}
\usage{
FRiP(object, assay, total.fragments, col.name = "FRiP", verbose = TRUE)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of the assay containing a peak x cell matrix}

\item{total.fragments}{Name of a metadata column containing the total number
of sequenced fragments for each cell. This can be computed using the
\code{\link{CountFragments}} function.}

\item{col.name}{Name of column in metadata to store the FRiP information.}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Calculate fraction of reads in peaks per cell
}
\examples{
FRiP(object = atac_small, assay = 'peaks', total.fragments = "fragments")
}
\concept{qc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{ExpressionPlot}
\alias{ExpressionPlot}
\title{Plot gene expression}
\usage{
ExpressionPlot(
  object,
  features,
  assay = NULL,
  group.by = NULL,
  idents = NULL,
  slot = "data"
)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A list of features to plot}

\item{assay}{Name of the assay storing expression information}

\item{group.by}{A grouping variable to group cells by. If NULL, use the
current cell identities}

\item{idents}{A list of identities to include in the plot. If NULL, include
all identities}

\item{slot}{Which slot to pull expression data from}
}
\description{
Display gene expression values for different groups of cells and different
genes. Genes will be arranged on the x-axis and different groups stacked on
the y-axis, with expression value distribution for each group shown as a
violin plot. This is designed to work alongside a genomic coverage track,
and the plot will be able to be aligned with coverage tracks for the same
groups of cells.
}
\examples{
\donttest{
ExpressionPlot(atac_small, features = "TSPAN6", assay = "RNA")
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/links.R
\name{ConnectionsToLinks}
\alias{ConnectionsToLinks}
\title{Cicero connections to links}
\usage{
ConnectionsToLinks(conns, ccans = NULL, threshold = 0, sep = c("-", "-"))
}
\arguments{
\item{conns}{A dataframe containing co-accessible elements. This would
usually be the output of \code{run_cicero} or
\code{assemble_connections}. Specifically, this should be a
dataframe where the first column contains the genomic coordinates of the
first element in the linked pair of elements, with chromosome, start, end
coordinates separated by "-" characters. The second column should be the
second element in the linked pair, formatted in the same way as the first
column. A third column should contain the co-accessibility scores.}

\item{ccans}{This is optional, but if supplied should be a dataframe
containing the cis-co-accessibility network (CCAN) information generated
by \code{generate_ccans}. Specifically, this should be a
dataframe containing the name of the peak in the first column, and the
CCAN that it belongs to in the second column.}

\item{threshold}{Threshold for retaining a coaccessible site. Links with
a value less than or equal to this threshold will be discarded.}

\item{sep}{Separators to use for strings encoding genomic coordinates.
First element is used to separate the chromosome from the coordinates, second
element is used to separate the start from end coordinate.}
}
\value{
Returns a \code{\link[GenomicRanges]{GRanges}} object
}
\description{
Convert the output of Cicero connections to a set of genomic ranges where
the start and end coordinates of the range are the midpoints of the linked
elements. Only elements on the same chromosome are included in the output.
}
\details{
See the Cicero package for more information:
\url{https://bioconductor.org/packages/cicero/}
}
\concept{links}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{CreateFragmentObject}
\alias{CreateFragmentObject}
\title{Create a Fragment object}
\usage{
CreateFragmentObject(
  path,
  cells = NULL,
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{path}{A path to the fragment file. The file should contain a tabix
index in the same directory.}

\item{cells}{A named character vector containing cell barcodes contained in
the fragment file. This does not need to be all cells in the fragment file,
but there should be no cells in the vector that are not present in the
fragment file. A search of the file will be performed until at least one
fragment from each cell is found. If NULL, don't check for expected cells.

Each element of the vector should be a cell barcode that appears in the
fragment file, and the name of each element should be the corresponding cell
name in the object.}

\item{validate.fragments}{Check that expected cells are present in the
fragment file.}

\item{verbose}{Display messages}

\item{...}{Additional arguments passed to \code{ValidateCells}}
}
\description{
Create a \code{Fragment} object to store fragment file information.
This object stores a 32-bit MD5 hash of the fragment file and the fragment
file index so that any changes to the files on-disk can be detected. A check
is also performed to ensure that the expected cells are present in the
fragment file.
}
\examples{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
cells <- colnames(x = atac_small)
names(x = cells) <- paste0("test_", cells)
frags <- CreateFragmentObject(path = fpath, cells = cells, verbose = FALSE, tolerance = 0.5)
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{head.Fragment}
\alias{head.Fragment}
\title{Return the first rows of a fragment file}
\usage{
\method{head}{Fragment}(x, n = 6L, ...)
}
\arguments{
\item{x}{a \code{Fragment} object}

\item{n}{an integer specifying the number of rows to return from the fragment
file}

\item{...}{additional arguments passed to \code{\link[utils]{read.table}}}
}
\value{
The first \code{n} rows of a fragment file as a \code{data.frame}
with the following columns: chrom, start, end, barcode, readCount.
}
\description{
Returns the first \code{n} rows of a fragment file. This allows the content
of a fragment file to be inspected.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{CreateChromatinAssay}
\alias{CreateChromatinAssay}
\title{Create ChromatinAssay object}
\usage{
CreateChromatinAssay(
  counts,
  data,
  min.cells = 0,
  min.features = 0,
  max.cells = NULL,
  ranges = NULL,
  motifs = NULL,
  fragments = NULL,
  genome = NULL,
  annotation = NULL,
  bias = NULL,
  positionEnrichment = NULL,
  sep = c("-", "-"),
  validate.fragments = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{counts}{Unnormalized data (raw counts)}

\item{data}{Normalized data; if provided, do not pass counts}

\item{min.cells}{Include features detected in at least this many cells.
Will subset the counts matrix as well.
To reintroduce excluded features, create a new object with a lower cutoff.}

\item{min.features}{Include cells where at least this many features are
detected.}

\item{max.cells}{Include features detected in less than this many cells.
Will subset the counts matrix as well.
To reintroduce excluded features, create a new object with a higher cutoff.
This can be useful for chromatin assays where certain artefactual loci
accumulate reads in all cells. A percentage cutoff can also be set using
'q' followed by the percentage of cells, for example 'q90' will discard
features detected in 90 percent of cells.
If NULL (default), do not apply any maximum value.}

\item{ranges}{A set of \code{\link[GenomicRanges]{GRanges}} corresponding to
the rows of the input matrix}

\item{motifs}{A Motif object (not required)}

\item{fragments}{Path to a tabix-indexed fragments file for the data
contained in the input matrix. If multiple fragment files are required,
you can add additional \code{\link{Fragment}} object to the assay after it is
created using the \code{\link{CreateFragmentObject}} and
\code{\link{Fragments}} functions. Alternatively, a list of
\code{\link{Fragment}} objects can be provided.}

\item{genome}{A \code{\link[GenomeInfoDb]{Seqinfo}} object containing basic
information about the genome used. Alternatively, the name of a UCSC genome
can be provided and the sequence information will be downloaded from UCSC.}

\item{annotation}{A set of \code{\link[GenomicRanges]{GRanges}} containing
annotations for the genome used}

\item{bias}{A Tn5 integration bias matrix}

\item{positionEnrichment}{A named list of matrices containing positional
signal enrichment information for each cell. Should be a cell x position
matrix, centered on an element of interest (for example, TSS sites).}

\item{sep}{Separators to use for strings encoding genomic coordinates.
First element is used to separate the chromosome from the coordinates,
second element is used to separate the start from end coordinate. Only
used if \code{ranges} is NULL.}

\item{validate.fragments}{Check that cells in the assay are present in the
fragment file.}

\item{verbose}{Display messages}

\item{...}{Additional arguments passed to \code{\link{CreateFragmentObject}}}
}
\description{
Create a \code{\link{ChromatinAssay}} object from a count matrix or
normalized data matrix. The expected format of the input matrix is features x
cells. A set of genomic ranges must be supplied along with the matrix, with
the length of the ranges equal to the number of rows in the matrix. If a set
of genomic ranges are not supplied, they will be extracted from the
row names of the matrix.
}
\concept{assay}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/motifs.R
\name{FindMotifs}
\alias{FindMotifs}
\title{FindMotifs}
\usage{
FindMotifs(
  object,
  features,
  background = 40000,
  assay = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{features}{A vector of features to test for enrichments over background}

\item{background}{Either a vector of features to use as the background set,
or a number specify the number of features to randomly select as a background
set. If a number is provided, regions will be selected to match the sequence
characteristics of the query features. To match the sequence characteristics,
these characteristics must be stored in the feature metadata for the assay.
This can be added using the
 \code{\link{RegionStats}} function. If NULL, use all features in the assay.}

\item{assay}{Which assay to use. Default is the active assay}

\item{verbose}{Display messages}

\item{...}{Arguments passed to \code{\link{MatchRegionStats}}.}
}
\value{
Returns a data frame
}
\description{
Find motifs over-represented in a given set of genomic features.
Computes the number of features containing the motif (observed) and
compares this to the total number of features containing the
motif (background) using the hypergeometric test.
}
\examples{
de.motif <- head(rownames(atac_small))
bg.peaks <- tail(rownames(atac_small))
FindMotifs(
  object = atac_small,
  features = de.motif,
  background = bg.peaks
)
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{Cells<-}
\alias{Cells<-}
\title{Set and get cell barcode information for a Fragment object}
\usage{
Cells(x, ...) <- value
}
\arguments{
\item{x}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{value}{A character vector of cell barcodes}
}
\description{
Set and get cell barcode information for a Fragment object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/quantification.R
\name{AggregateTiles}
\alias{AggregateTiles}
\alias{AggregateTiles.Seurat}
\alias{AggregateTiles.ChromatinAssay}
\alias{AggregateTiles.default}
\title{Quantify aggregated genome tiles}
\usage{
AggregateTiles(object, ...)

\method{AggregateTiles}{Seurat}(
  object,
  genome,
  assay = NULL,
  new.assay.name = "tiles",
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
)

\method{AggregateTiles}{ChromatinAssay}(
  object,
  genome,
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
)

\method{AggregateTiles}{default}(
  object,
  genome,
  cells = NULL,
  min_counts = 5,
  binsize = 5000,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{object}{A Seurat object or ChromatinAssay object}

\item{...}{Additional arguments passed to other methods}

\item{genome}{genome A vector of chromosome sizes for the genome. This is
used to construct the genome bin coordinates. The can be obtained by calling
seqlengths on a BSgenome-class object.}

\item{assay}{Name of assay to use}

\item{new.assay.name}{Name of new assay to create containing aggregated
genome tiles}

\item{min_counts}{Minimum number of counts for a tile to be retained prior to
aggregation}

\item{binsize}{Size of the genome bins (tiles) in base pairs}

\item{verbose}{Display messages}

\item{cells}{Cells to include}
}
\value{
When running on a Seurat object, returns the Seurat object with a new
\code{\link{ChromatinAssay}} added.

When running on a \code{\link{ChromatinAssay}}, returns a new
\code{ChromatinAssay} containing the aggregated genome tiles.

When running on a fragment file, returns a sparse region x cell
matrix.
}
\description{
Quantifies fragment counts per cell in fixed-size genome bins across the
whole genome, then removes bins with less than a desired minimum number of
counts in the bin, then merges adjacent tiles into a single region.
}
\concept{quantification}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{LinkPlot}
\alias{LinkPlot}
\title{Plot linked genomic elements}
\usage{
LinkPlot(object, region, min.cutoff = 0)
}
\arguments{
\item{object}{A \code{\link[SeuratObject]{Seurat}} object}

\item{region}{A genomic region to plot}

\item{min.cutoff}{Minimum absolute score for link to be plotted.}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Display links between pairs of genomic elements within a given region of the
genome.
}
\concept{links}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GetGRangesFromEnsDb}
\alias{GetGRangesFromEnsDb}
\title{Extract genomic ranges from EnsDb object}
\usage{
GetGRangesFromEnsDb(
  ensdb,
  standard.chromosomes = TRUE,
  biotypes = c("protein_coding", "lincRNA", "rRNA", "processed_transcript"),
  verbose = TRUE
)
}
\arguments{
\item{ensdb}{An EnsDb object}

\item{standard.chromosomes}{Keep only standard chromosomes}

\item{biotypes}{Biotypes to keep}

\item{verbose}{Display messages}
}
\description{
Pulls the transcript information for all chromosomes from an EnsDb object.
This wraps \code{\link[biovizBase]{crunch}} and applies the extractor
function to all chromosomes present in the EnsDb object.
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomeinfodb-methods.R
\name{seqinfo-methods}
\alias{seqinfo-methods}
\alias{seqinfo}
\alias{seqinfo,ChromatinAssay-method}
\alias{seqinfo<-,ChromatinAssay-method}
\alias{seqlevels,ChromatinAssay-method}
\alias{seqlevels}
\alias{seqlevels<-,ChromatinAssay-method}
\alias{seqnames,ChromatinAssay-method}
\alias{seqnames}
\alias{seqnames<-,ChromatinAssay-method}
\alias{seqlengths,ChromatinAssay-method}
\alias{seqlengths}
\alias{seqlengths<-,ChromatinAssay-method}
\alias{genome,ChromatinAssay-method}
\alias{genome}
\alias{genome<-,ChromatinAssay-method}
\alias{isCircular,ChromatinAssay-method}
\alias{isCircular}
\alias{isCircular<-,ChromatinAssay-method}
\alias{seqinfo,Seurat-method}
\alias{seqinfo<-,Seurat-method}
\alias{seqlevels,Seurat-method}
\alias{seqlevels<-,Seurat-method}
\alias{seqnames,Seurat-method}
\alias{seqnames<-,Seurat-method}
\alias{seqlengths,Seurat-method}
\alias{seqlengths<-,Seurat-method}
\alias{genome,Seurat-method}
\alias{genome<-,Seurat-method}
\alias{isCircular,Seurat-method}
\alias{isCircular<-,Seurat-method}
\title{Access and modify sequence information for ChromatinAssay objects}
\usage{
\S4method{seqinfo}{ChromatinAssay}(x)

\S4method{seqinfo}{ChromatinAssay}(x) <- value

\S4method{seqlevels}{ChromatinAssay}(x)

\S4method{seqlevels}{ChromatinAssay}(x) <- value

\S4method{seqnames}{ChromatinAssay}(x)

\S4method{seqnames}{ChromatinAssay}(x) <- value

\S4method{seqlengths}{ChromatinAssay}(x)

\S4method{seqlengths}{ChromatinAssay}(x) <- value

\S4method{genome}{ChromatinAssay}(x)

\S4method{genome}{ChromatinAssay}(x) <- value

\S4method{isCircular}{ChromatinAssay}(x)

\S4method{isCircular}{ChromatinAssay}(x) <- value

\S4method{seqinfo}{Seurat}(x)

\S4method{seqinfo}{Seurat}(x) <- value

\S4method{seqlevels}{Seurat}(x)

\S4method{seqlevels}{Seurat}(x) <- value

\S4method{seqnames}{Seurat}(x)

\S4method{seqnames}{Seurat}(x) <- value

\S4method{seqlengths}{Seurat}(x)

\S4method{seqlengths}{Seurat}(x) <- value

\S4method{genome}{Seurat}(x)

\S4method{genome}{Seurat}(x) <- value

\S4method{isCircular}{Seurat}(x)

\S4method{isCircular}{Seurat}(x) <- value
}
\arguments{
\item{x}{A \code{\link{ChromatinAssay}} object}

\item{value}{A \code{\link[GenomeInfoDb]{Seqinfo}} object or name of a UCSC
genome to store in the \code{\link{ChromatinAssay}}}
}
\description{
Methods for accessing and modifying
\code{\link[GenomeInfoDb]{Seqinfo}} object information stored in a
\code{\link{ChromatinAssay}} object.
}
\section{Functions}{
\itemize{
\item \code{seqinfo<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{seqlevels,ChromatinAssay-method}: get method for ChromatinAssay objects

\item \code{seqlevels<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{seqnames,ChromatinAssay-method}: get method for ChromatinAssay objects

\item \code{seqnames<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{seqlengths,ChromatinAssay-method}: get method for ChromatinAssay objects

\item \code{seqlengths<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{genome,ChromatinAssay-method}: get method for ChromatinAssay objects

\item \code{genome<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{isCircular,ChromatinAssay-method}: get method for ChromatinAssay objects

\item \code{isCircular<-,ChromatinAssay-method}: set method for ChromatinAssay objects

\item \code{seqinfo,Seurat-method}: get method for Seurat objects

\item \code{seqinfo<-,Seurat-method}: set method for Seurat objects

\item \code{seqlevels,Seurat-method}: get method for Seurat objects

\item \code{seqlevels<-,Seurat-method}: set method for Seurat objects

\item \code{seqnames,Seurat-method}: get method for Seurat objects

\item \code{seqnames<-,Seurat-method}: set method for Seurat objects

\item \code{seqlengths,Seurat-method}: get method for Seurat objects

\item \code{seqlengths<-,Seurat-method}: set method for Seurat objects

\item \code{genome,Seurat-method}: get method for Seurat objects

\item \code{genome<-,Seurat-method}: set method for Seurat objects

\item \code{isCircular,Seurat-method}: get method for Seurat objects

\item \code{isCircular<-,Seurat-method}: set method for Seurat objects
}}

\seealso{
\itemize{
  \item{\link[GenomeInfoDb]{seqinfo} in the \pkg{GenomeInfoDb} package.}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{seqinfo}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{TSSEnrichment}
\alias{TSSEnrichment}
\title{Compute TSS enrichment score per cell}
\usage{
TSSEnrichment(
  object,
  tss.positions = NULL,
  n = NULL,
  fast = TRUE,
  assay = NULL,
  cells = NULL,
  process_n = 2000,
  verbose = TRUE
)
}
\arguments{
\item{object}{A Seurat object}

\item{tss.positions}{A GRanges object containing the TSS positions. If NULL,
use the genomic annotations stored in the assay.}

\item{n}{Number of TSS positions to use. This will select the first _n_
TSSs from the set. If NULL, use all TSSs (slower).}

\item{fast}{Just compute the TSS enrichment score, without storing the
base-resolution matrix of integration counts at each site. This reduces the
memory required to store the object but does not allow plotting the
accessibility profile at the TSS.}

\item{assay}{Name of assay to use}

\item{cells}{A vector of cells to include. If NULL (default), use all cells
in the object}

\item{process_n}{Number of regions to process at a time if using \code{fast}
option.}

\item{verbose}{Display messages}
}
\value{
Returns a \code{\link[SeuratObject]{Seurat}} object
}
\description{
Compute the transcription start site (TSS) enrichment score for each cell,
as defined by ENCODE:
\url{https://www.encodeproject.org/data-standards/terms/}.
}
\details{
The computed score will be added to the object metadata as "TSS.enrichment".
}
\examples{
\dontrun{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
Fragments(atac_small) <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  tolerance = 0.5
)
TSSEnrichment(object = atac_small)
}
}
\concept{qc}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/preprocessing.R
\name{RegionStats}
\alias{RegionStats}
\alias{RegionStats.default}
\alias{RegionStats.ChromatinAssay}
\alias{RegionStats.Seurat}
\title{Compute base composition information for genomic ranges}
\usage{
RegionStats(object, ...)

\method{RegionStats}{default}(object, genome, verbose = TRUE, ...)

\method{RegionStats}{ChromatinAssay}(object, genome, verbose = TRUE, ...)

\method{RegionStats}{Seurat}(object, genome, assay = NULL, verbose = TRUE, ...)
}
\arguments{
\item{object}{A Seurat object, Assay object, or set of genomic ranges}

\item{...}{Arguments passed to other methods}

\item{genome}{A BSgenome object}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use}
}
\value{
Returns a dataframe
}
\description{
Compute the GC content, region lengths, and dinucleotide base frequencies
for regions in the assay and add to the feature metadata.
}
\examples{
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
RegionStats(
  object = rownames(atac_small),
  genome = BSgenome.Hsapiens.UCSC.hg19
)
}
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
RegionStats(
  object = atac_small[['peaks']],
  genome = BSgenome.Hsapiens.UCSC.hg19
)
}
\dontrun{
library(BSgenome.Hsapiens.UCSC.hg19)
RegionStats(
  object = atac_small,
  assay = 'bins',
  genome = BSgenome.Hsapiens.UCSC.hg19
)
}
}
\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{AnnotationPlot}
\alias{AnnotationPlot}
\title{Plot gene annotations}
\usage{
AnnotationPlot(object, region)
}
\arguments{
\item{object}{A \code{\link[SeuratObject]{Seurat}} object}

\item{region}{A genomic region to plot}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Display gene annotations in a given region of the genome.
}
\examples{
\donttest{
AnnotationPlot(object = atac_small, region = c("chr1-29554-39554"))
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{coverage,ChromatinAssay-method}
\alias{coverage,ChromatinAssay-method}
\alias{coverage}
\alias{coverage,Seurat-method}
\title{Coverage of a ChromatinAssay object}
\usage{
\S4method{coverage}{ChromatinAssay}(
  x,
  shift = 0L,
  width = NULL,
  weight = 1L,
  method = c("auto", "sort", "hash")
)

\S4method{coverage}{Seurat}(
  x,
  shift = 0L,
  width = NULL,
  weight = 1L,
  method = c("auto", "sort", "hash")
)
}
\arguments{
\item{x}{A \code{\link{ChromatinAssay}} object}

\item{shift}{How much each range should be shifted before coverage is
computed. See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.}

\item{width}{Specifies the length of the returned coverage vectors.
See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.}

\item{weight}{Assigns weight to each range in \code{x}.
See \code{\link[IRanges]{coverage}} in the \pkg{IRanges} package.}

\item{method}{See \code{\link[IRanges]{coverage}} in the \pkg{IRanges}
package}
}
\description{
This is the \code{coverage} method for \code{\link{ChromatinAssay}} objects.
}
\section{Functions}{
\itemize{
\item \code{coverage,ChromatinAssay-method}: method for ChromatinAssay objects

\item \code{coverage,Seurat-method}: method for Seurat objects
}}

\seealso{
\itemize{
  \item{\link[IRanges]{coverage-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{coverage-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{coverage}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{FragmentHistogram}
\alias{FragmentHistogram}
\title{Plot fragment length histogram}
\usage{
FragmentHistogram(
  object,
  assay = NULL,
  region = "chr1-1-2000000",
  group.by = NULL,
  cells = NULL,
  log.scale = FALSE,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Which assay to use. Default is the active assay.}

\item{region}{Genomic range to use. Default is fist two megabases of
chromosome 1. Can be a GRanges object, a string, or a vector
of strings.}

\item{group.by}{Name of one or more metadata columns to group (color) the
cells by. Default is the current cell identities}

\item{cells}{Which cells to plot. Default all cells}

\item{log.scale}{Display Y-axis on log scale. Default is FALSE.}

\item{...}{Arguments passed to other functions}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot the frequency that fragments of different lengths are present for
different groups of cells.
}
\examples{
\donttest{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
Fragments(atac_small) <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
FragmentHistogram(object = atac_small, region = "chr1-10245-780007")
}
}
\concept{qc}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{GetTSSPositions}
\alias{GetTSSPositions}
\title{Find transcriptional start sites}
\usage{
GetTSSPositions(ranges, biotypes = "protein_coding")
}
\arguments{
\item{ranges}{A GRanges object containing gene annotations.}

\item{biotypes}{Gene biotypes to include. If NULL, use all biotypes in the
supplied gene annotation.}
}
\description{
Get the TSS positions from a set of genomic ranges containing gene positions.
Ranges can contain exons, introns, UTRs, etc, rather than the whole
transcript. Only protein coding gene biotypes are included in output.
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{TSSPlot}
\alias{TSSPlot}
\title{Plot signal enrichment around TSSs}
\usage{
TSSPlot(object, assay = NULL, group.by = NULL, idents = NULL)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of the assay to use. Should have the TSS enrichment
information for each cell
already computed by running \code{\link{TSSEnrichment}}}

\item{group.by}{Set of identities to group cells by}

\item{idents}{Set of identities to include in the plot}
}
\value{
Returns a \code{\link[ggplot2]{ggplot2}} object
}
\description{
Plot the normalized TSS enrichment score at each position relative to the
TSS. Requires that \code{\link{TSSEnrichment}} has already been run on the
assay.
}
\concept{qc}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\docType{class}
\name{Motif-class}
\alias{Motif-class}
\alias{Motif}
\title{The Motif class}
\description{
The Motif class is designed to store DNA sequence motif information,
including motif PWMs or PFMs, motif positions, and metadata.
}
\section{Slots}{

\describe{
\item{\code{data}}{A sparse, binary, feature x motif matrix. Columns
correspond to motif IDs, rows correspond to genomic features
(peaks or bins). Entries in the matrix should be 1 if the
genomic feature contains the motif, and 0 otherwise.}

\item{\code{pwm}}{A named list of position weight matrices}

\item{\code{motif.names}}{A list containing the name of each motif}

\item{\code{positions}}{A \code{\link[GenomicRanges]{GRangesList}} object containing
exact positions of each motif.}

\item{\code{meta.data}}{A dataframe for storage of additional
information related to each motif. This could include the
names of proteins that bind the motif.}
}}

\concept{motifs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{CountsInRegion}
\alias{CountsInRegion}
\title{Counts in region}
\usage{
CountsInRegion(object, assay, regions, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{assay}{Name of a chromatin assay in the object to use}

\item{regions}{A GRanges object}

\item{...}{Additional arguments passed to \code{\link[IRanges]{findOverlaps}}}
}
\value{
Returns a numeric vector
}
\description{
Count reads per cell overlapping a given set of regions
}
\examples{
\donttest{
CountsInRegion(
  object = atac_small,
  assay = 'bins',
  regions = blacklist_hg19
)
}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{CoveragePlot}
\alias{CoveragePlot}
\title{Plot Tn5 insertion frequency over a region}
\usage{
CoveragePlot(
  object,
  region,
  features = NULL,
  assay = NULL,
  show.bulk = FALSE,
  expression.assay = "RNA",
  expression.slot = "data",
  annotation = TRUE,
  peaks = TRUE,
  peaks.group.by = NULL,
  ranges = NULL,
  ranges.group.by = NULL,
  ranges.title = "Ranges",
  region.highlight = NULL,
  links = TRUE,
  tile = FALSE,
  tile.size = 100,
  tile.cells = 100,
  bigwig = NULL,
  bigwig.type = "coverage",
  heights = NULL,
  group.by = NULL,
  window = 100,
  extend.upstream = 0,
  extend.downstream = 0,
  scale.factor = NULL,
  ymax = NULL,
  cells = NULL,
  idents = NULL,
  sep = c("-", "-"),
  max.downsample = 3000,
  downsample.rate = 0.1,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{region}{A set of genomic coordinates to show. Can be a GRanges object,
a string encoding a genomic position, a gene name, or a vector of strings
describing the genomic coordinates or gene names to plot. If a gene name is
supplied, annotations must be present in the assay.}

\item{features}{A vector of features present in another assay to plot
alongside accessibility tracks (for example, gene names).}

\item{assay}{Name of the assay to plot}

\item{show.bulk}{Include coverage track for all cells combined (pseudo-bulk).
Note that this will plot the combined accessibility for all cells included in
the plot (rather than all cells in the object).}

\item{expression.assay}{Name of the assay containing expression data to plot
alongside accessibility tracks. Only needed if supplying \code{features}
argument.}

\item{expression.slot}{Name of slot to pull expression data from. Only needed
if supplying the \code{features} argument.}

\item{annotation}{Display gene annotations}

\item{peaks}{Display peaks}

\item{peaks.group.by}{Grouping variable to color peaks by. Must be a variable
present in the feature metadata. If NULL, do not color peaks by any variable.}

\item{ranges}{Additional genomic ranges to plot}

\item{ranges.group.by}{Grouping variable to color ranges by. Must be a
variable present in the metadata stored in the \code{ranges} genomic ranges.
If NULL, do not color by any variable.}

\item{ranges.title}{Y-axis title for ranges track. Only relevant if
\code{ranges} parameter is set.}

\item{region.highlight}{Region to highlight on the plot. Should be a GRanges
object containing the coordinates to highlight. By default, regions will be
highlighted in grey. To change the color of the highlighting, include a
metadata column in the GRanges object named "color" containing the color to
use for each region.}

\item{links}{Display links}

\item{tile}{Display per-cell fragment information in sliding windows.}

\item{tile.size}{Size of the sliding window for per-cell fragment tile plot}

\item{tile.cells}{Number of cells to display fragment information for in tile
plot.}

\item{bigwig}{List of bigWig file paths to plot data from. Files can be
remotely hosted. The name of each element in the list will determine the
y-axis label given to the track.}

\item{bigwig.type}{Type of track to use for bigWig files ("line", "heatmap",
or "coverage"). Should either be a single value, or a list of values giving
the type for each individual track in the provided list of bigwig files.}

\item{heights}{Relative heights for each track (accessibility, gene
annotations, peaks, links).}

\item{group.by}{Name of one or more metadata columns to group (color) the
cells by. Default is the current cell identities}

\item{window}{Smoothing window size}

\item{extend.upstream}{Number of bases to extend the region upstream.}

\item{extend.downstream}{Number of bases to extend the region downstream.}

\item{scale.factor}{Scaling factor for track height. If NULL (default),
use the median group scaling factor determined by total number of fragments
sequences in each group.}

\item{ymax}{Maximum value for Y axis. If NULL (default) set to the highest
value among all the tracks.}

\item{cells}{Which cells to plot. Default all cells}

\item{idents}{Which identities to include in the plot. Default is all
identities.}

\item{sep}{Separators to use for strings encoding genomic coordinates. First
element is used to separate the chromosome from the coordinates, second
element is used to separate the start from end coordinate.}

\item{max.downsample}{Minimum number of positions kept when downsampling.
Downsampling rate is adaptive to the window size, but this parameter will set
the minimum possible number of positions to include so that plots do not
become too sparse when the window size is small.}

\item{downsample.rate}{Fraction of positions to retain when downsampling.
Retaining more positions can give a higher-resolution plot but can make the
number of points large, resulting in larger file sizes when saving the plot
and a longer period of time needed to draw the plot.}

\item{...}{Additional arguments passed to \code{\link[patchwork]{wrap_plots}}}
}
\value{
Returns a \code{\link[ggplot2]{ggplot}} object
}
\description{
Plot frequency of Tn5 insertion events for different groups of cells within
given regions of the genome.
}
\details{
Thanks to Andrew Hill for providing an early version of this function.
}
\examples{
\donttest{
fpath <- system.file("extdata", "fragments.tsv.gz", package="Signac")
fragments <- CreateFragmentObject(
  path = fpath,
  cells = colnames(atac_small),
  validate.fragments = FALSE
)
Fragments(atac_small) <- fragments

# Basic coverage plot
CoveragePlot(object = atac_small, region = c("chr1-713500-714500"))

# Show additional ranges
ranges.show <- StringToGRanges("chr1-713750-714000")
CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), ranges = ranges.show)

# Highlight region
CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), region.highlight = ranges.show)

# Change highlight color
ranges.show$color <- "orange"
CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), region.highlight = ranges.show)

# Show expression data
CoveragePlot(object = atac_small, region = c("chr1-713500-714500"), features = "ELK1")
}
}
\concept{visualization}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{ValidateHash}
\alias{ValidateHash}
\title{Validate hashes for Fragment object}
\usage{
ValidateHash(object, verbose = TRUE)
}
\arguments{
\item{object}{A \code{\link{Fragment}} object}

\item{verbose}{Display messages}
}
\description{
Validate hashes for Fragment object
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iranges-methods.R
\name{nearest-methods}
\alias{nearest-methods}
\alias{precede}
\alias{precede,ANY,ChromatinAssay-method}
\alias{precede,ChromatinAssay,ANY-method}
\alias{precede,ChromatinAssay,ChromatinAssay-method}
\alias{precede,ANY,Seurat-method}
\alias{precede,Seurat,ANY-method}
\alias{precede,Seurat,Seurat-method}
\alias{follow,ANY,ChromatinAssay-method}
\alias{follow}
\alias{follow,ChromatinAssay,ANY-method}
\alias{follow,ChromatinAssay,ChromatinAssay-method}
\alias{follow,ANY,Seurat-method}
\alias{follow,Seurat,ANY-method}
\alias{follow,Seurat,Seurat-method}
\alias{nearest,ANY,ChromatinAssay-method}
\alias{nearest}
\alias{nearest,ChromatinAssay,ANY-method}
\alias{nearest,ChromatinAssay,ChromatinAssay-method}
\alias{nearest,ANY,Seurat-method}
\alias{nearest,Seurat,ANY-method}
\alias{nearest,Seurat,Seurat-method}
\alias{distance,ANY,ChromatinAssay-method}
\alias{distance}
\alias{distance,ChromatinAssay,ANY-method}
\alias{distance,ChromatinAssay,ChromatinAssay-method}
\alias{distance,ANY,Seurat-method}
\alias{distance,Seurat,ANY-method}
\alias{distance,Seurat,Seurat-method}
\alias{distanceToNearest,ANY,ChromatinAssay-method}
\alias{distanceToNearest}
\alias{distanceToNearest,ChromatinAssay,ANY-method}
\alias{distanceToNearest,ChromatinAssay,ChromatinAssay-method}
\alias{distanceToNearest,ANY,Seurat-method}
\alias{distanceToNearest,Seurat,ANY-method}
\alias{distanceToNearest,Seurat,Seurat-method}
\title{Find the nearest range neighbors for ChromatinAssay objects}
\usage{
\S4method{precede}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{precede}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{follow}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ANY,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ChromatinAssay,ChromatinAssay}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{ANY,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{Seurat,ANY}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{nearest}{Seurat,Seurat}(x, subject, select = c("arbitrary", "all"), ignore.strand = FALSE)

\S4method{distance}{ANY,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ANY}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ChromatinAssay,ChromatinAssay}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{ANY,Seurat}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{Seurat,ANY}(x, y, ignore.strand = FALSE, ...)

\S4method{distance}{Seurat,Seurat}(x, y, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ANY,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ANY}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ChromatinAssay,ChromatinAssay}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{ANY,Seurat}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{Seurat,ANY}(x, subject, ignore.strand = FALSE, ...)

\S4method{distanceToNearest}{Seurat,Seurat}(x, subject, ignore.strand = FALSE, ...)
}
\arguments{
\item{x}{A query \code{\link{ChromatinAssay}} object}

\item{subject}{The subject \code{\link[GenomicRanges]{GRanges}} or
\code{\link{ChromatinAssay}} object. If missing, \code{x} is used as the
subject.}

\item{select}{Logic for handling ties.
See \code{\link[GenomicRanges]{nearest-methods}} in the \pkg{GenomicRanges}
package.}

\item{ignore.strand}{Logical argument controlling whether strand information
should be ignored.}

\item{y}{For the \code{distance} method, a
\code{\link[GenomicRanges]{GRanges}} object or a \code{\link{ChromatinAssay}}
object}

\item{...}{Additional arguments for methods}
}
\description{
The \code{precede, follow, nearest, distance, distanceToNearest} methods
are available for \code{\link{ChromatinAssay}} objects.
}
\section{Functions}{
\itemize{
\item \code{precede,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{precede,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{precede,ANY,Seurat-method}: method for ANY, Seurat

\item \code{precede,Seurat,ANY-method}: method for Seurat, ANY

\item \code{precede,Seurat,Seurat-method}: method for Seurat, Seurat

\item \code{follow,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{follow,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{follow,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{follow,ANY,Seurat-method}: method for ANY, Seurat

\item \code{follow,Seurat,ANY-method}: method for Seurat, ANY

\item \code{follow,Seurat,Seurat-method}: method for Seurat, Seurat

\item \code{nearest,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{nearest,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{nearest,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{nearest,ANY,Seurat-method}: method for ANY, Seurat

\item \code{nearest,Seurat,ANY-method}: method for Seurat, ANY

\item \code{nearest,Seurat,Seurat-method}: method for Seurat, Seurat

\item \code{distance,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{distance,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{distance,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{distance,ANY,Seurat-method}: method for ANY, Seurat

\item \code{distance,Seurat,ANY-method}: method for Seurat, ANY

\item \code{distance,Seurat,Seurat-method}: method for Seurat, Seurat

\item \code{distanceToNearest,ANY,ChromatinAssay-method}: method for ANY, ChromatinAssay

\item \code{distanceToNearest,ChromatinAssay,ANY-method}: method for ChromatinAssay, ANY

\item \code{distanceToNearest,ChromatinAssay,ChromatinAssay-method}: method for ChromatinAssay, ChromatinAssay

\item \code{distanceToNearest,ANY,Seurat-method}: method for ANY, Seurat

\item \code{distanceToNearest,Seurat,ANY-method}: method for Seurat, ANY

\item \code{distanceToNearest,Seurat,Seurat-method}: method for Seurat, Seurat
}}

\seealso{
\itemize{
  \item{\link[IRanges]{nearest-methods} in the \pkg{IRanges} package.}
  \item{\link[GenomicRanges]{nearest-methods} in the \pkg{GenomicRanges}
  package}
  \item{\link{ChromatinAssay-class}}
 }
}
\concept{nearest}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/mito.R
\name{IdentifyVariants}
\alias{IdentifyVariants}
\alias{IdentifyVariants.default}
\alias{IdentifyVariants.Assay}
\alias{IdentifyVariants.Seurat}
\title{Identify mitochondrial variants}
\usage{
IdentifyVariants(object, ...)

\method{IdentifyVariants}{default}(
  object,
  refallele,
  stabilize_variance = TRUE,
  low_coverage_threshold = 10,
  verbose = TRUE,
  ...
)

\method{IdentifyVariants}{Assay}(object, refallele, ...)

\method{IdentifyVariants}{Seurat}(object, refallele, assay = NULL, ...)
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{refallele}{A dataframe containing reference alleles for the
mitochondrial genome.}

\item{stabilize_variance}{Stabilize variance}

\item{low_coverage_threshold}{Low coverage threshold}

\item{verbose}{Display messages}

\item{assay}{Name of assay to use. If NULL, use the default assay.}
}
\value{
Returns a dataframe
}
\description{
Identify mitochondrial variants present in single cells.
}
\examples{
\dontrun{
data.dir <- "path/to/data/directory"
mgatk <- ReadMGATK(dir = data.dir)
variant.df <- IdentifyVariants(
  object = mgatk$counts,
  refallele = mgatk$refallele
)
}
}
\concept{mito}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fragments.R
\name{ValidateFragments}
\alias{ValidateFragments}
\title{Validate Fragment object}
\usage{
ValidateFragments(object, verbose = TRUE, ...)
}
\arguments{
\item{object}{A \code{\link{Fragment}} object}

\item{verbose}{Display messages}

\item{...}{Additional parameters passed to \code{\link{ValidateCells}}}
}
\description{
Verify that the cells listed in the object exist in the fragment file
and that the fragment file or index have not changed since creating the
fragment object.
}
\concept{fragments}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R, R/objects.R
\name{Links}
\alias{Links}
\alias{Links<-}
\alias{Links.ChromatinAssay}
\alias{Links.Seurat}
\alias{Links<-.ChromatinAssay}
\alias{Links<-.Seurat}
\title{Get or set links information}
\usage{
Links(object, ...)

Links(object, ...) <- value

\method{Links}{ChromatinAssay}(object, ...)

\method{Links}{Seurat}(object, ...)

\method{Links}{ChromatinAssay}(object, ...) <- value

\method{Links}{Seurat}(object, ...) <- value
}
\arguments{
\item{object}{A Seurat object}

\item{...}{Arguments passed to other methods}

\item{value}{A \code{\link[GenomicRanges]{GRanges}} object}
}
\description{
Get or set the genomic link information for a Seurat object or ChromatinAssay
}
\examples{
Links(atac_small[["peaks"]])
Links(atac_small)
links <- Links(atac_small)
Links(atac_small[["peaks"]]) <- links
links <- Links(atac_small)
Links(atac_small) <- links
}
\concept{assay}
\concept{links}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{MatchRegionStats}
\alias{MatchRegionStats}
\title{Match DNA sequence characteristics}
\usage{
MatchRegionStats(
  meta.feature,
  query.feature,
  features.match = c("GC.percent"),
  n = 10000,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{meta.feature}{A dataframe containing DNA sequence information for
features to choose from}

\item{query.feature}{A dataframe containing DNA sequence information for
features to match.}

\item{features.match}{Which features of the query to match when selecting a
set of regions. A vector of column names present in the feature metadata can
be supplied to match multiple characteristics at once. Default is GC content.}

\item{n}{Number of regions to select, with characteristics matching the query}

\item{verbose}{Display messages}

\item{...}{Arguments passed to other functions}
}
\value{
Returns a character vector
}
\description{
Return a vector if genomic regions that match the distribution of a set of
query regions for any given set of characteristics, specified in the input
\code{meta.feature} dataframe.
}
\details{
For each requested feature to match, a density
distribution is estimated using the \code{\link[stats]{density}} function,
and a set of weights for each feature in the dataset estimated based on the
density distribution. If multiple features are to be matched (for example,
GC content and overall accessibility), a joint density distribution is then
computed by multiplying the individual feature weights. A set of features
with characteristics matching the query regions is then selected using the
\code{\link[base]{sample}} function, with the probability of randomly
selecting each feature equal to the joint density distribution weight.
}
\examples{
metafeatures <- Seurat::GetAssayData(
  object = atac_small[['peaks']], slot = 'meta.features'
)
query.feature <- metafeatures[1:10, ]
features.choose <- metafeatures[11:nrow(metafeatures), ]
MatchRegionStats(
  meta.feature = features.choose,
  query.feature = query.feature,
  features.match = "percentile",
  n = 10
)
}
\concept{motifs}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{blacklist_ce11}
\alias{blacklist_ce11}
\title{Genomic blacklist regions for C. elegans ce11 (0-based)}
\format{
A GRanges object
}
\source{
\url{https://github.com/Boyle-Lab/Blacklist}

\doi{10.1038/s41598-019-45839-z}
}
\usage{
blacklist_ce11
}
\description{
Genomic blacklist regions for C. elegans ce11 (0-based)
}
\concept{data}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{LookupGeneCoords}
\alias{LookupGeneCoords}
\title{Get gene coordinates}
\usage{
LookupGeneCoords(object, gene, assay = NULL)
}
\arguments{
\item{object}{A Seurat object}

\item{gene}{Name of a gene to extract}

\item{assay}{Name of assay to use}
}
\description{
Extract the coordinates of the longest transcript for a gene stored in the
annotations within an object.
}
\examples{
LookupGeneCoords(atac_small, gene = "MIR1302-10")
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{CreateMotifMatrix}
\alias{CreateMotifMatrix}
\title{Create motif matrix}
\usage{
CreateMotifMatrix(
  features,
  pwm,
  genome,
  score = FALSE,
  use.counts = FALSE,
  sep = c("-", "-"),
  ...
)
}
\arguments{
\item{features}{A GRanges object containing a set of genomic features}

\item{pwm}{A \code{\link[TFBSTools]{PFMatrixList}} or
\code{\link[TFBSTools]{PWMatrixList}}
object containing position weight/frequency matrices to use}

\item{genome}{Any object compatible with the \code{genome} argument
in \code{\link[motifmatchr]{matchMotifs}}}

\item{score}{Record the motif match score, rather than presence/absence
(default FALSE)}

\item{use.counts}{Record motif counts per region. If FALSE (default),
record presence/absence of motif. Only applicable if \code{score=FALSE}.}

\item{sep}{A length-2 character vector containing the separators to be used
when constructing matrix rownames from the GRanges}

\item{...}{Additional arguments passed to
\code{\link[motifmatchr]{matchMotifs}}}
}
\value{
Returns a sparse matrix
}
\description{
Create a motif x feature matrix from a set of genomic ranges,
the genome, and a set of position weight matrices.
}
\details{
Requires that motifmatchr is installed
\url{https://www.bioconductor.org/packages/motifmatchr/}.
}
\examples{
\dontrun{
library(JASPAR2018)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg19)

pwm <- getMatrixSet(
  x = JASPAR2018,
  opts = list(species = 9606, all_versions = FALSE)
)
motif.matrix <- CreateMotifMatrix(
  features = granges(atac_small),
  pwm = pwm,
  genome = BSgenome.Hsapiens.UCSC.hg19
)
}
}
\concept{motifs}
\concept{preprocessing}
