# SoupX

An R package for the estimation and removal of cell free mRNA contamination in droplet based single cell RNA-seq data.

The problem this package attempts to solve is that all droplet based single cell RNA-seq experiments also capture ambient mRNAs present in the input solution along with cell specific mRNAs of interest.  This contamination is ubiquitous and can vary hugely between experiments (2% - 50%), although around 10% seems reasonably common.

There's no way to know in advance what the contamination is in an experiment, although solid tumours and low-viability cells tend to produce higher contamination fractions.  As the source of the contaminating mRNAs is lysed cells in the input solution, the profile of the contamination is experiment specific and produces a batch effect. 

Even if you decide you don't want to use the SoupX correction methods for whatever reason, you should at least want to know how contaminated your data are.

**NOTE:** From v1.3.0 onward SoupX now includes an option to automatically estimate the contamination fraction.  It is anticipated that this will be the preferred way of using the method for the vast majority of users.  This function (`autoEstCont`) depends on clustering information being provided.  If you are using 10X data mapped with cellranger, this will be loaded automatically, but otherwise it must be provided explicitly by the user using `setClusters`.

## Installation

The latest stable release can be installed from CRAN in the usual way by running,

```R
install.packages('SoupX')
```

If you want to use the latest development version, install it by running,

```R
devtools::install_github("constantAmateur/SoupX",ref='devel')
```

Finally, if you want to use the per-cell contamination estimation (which you almost certainly won't need to), install the branch STAN

```R
devtools::install_github("constantAmateur/SoupX",ref='STAN')
```

If you encounter errors saying `multtest` is unavailable, please install this manually from bioconductor with:

```R
BiocManager::install('multtest')
```

## Quickstart

Decontaminate one channel of 10X data mapped with cellranger by running:

```R
sc = load10X('path/to/your/cellranger/outs/folder')
sc = autoEstCont(sc)
out = adjustCounts(sc)
```

or to manually load decontaminate any other data

```R
sc = SoupChannel(table_of_droplets,table_of_counts)
sc = setClusters(sc,cluster_labels)
sc = autoEstCont(sc)
out = adjustCounts(sc)
```

`out` will then contain a corrected matrix to be used in place of the original table of counts in downstream analyses.


## Documentation

The methodology implemented in this package is explained in detail in [this paper](https://doi.org/10.1093/gigascience/giaa151).  

A detailed vignette is provided with the package and can be viewed [here](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html).  

# Citing SoupX

If you use SoupX in your work, please cite: "Young, M.D., Behjati, S. (2020). SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data, GigaScience, Volume 9, Issue 12, December 2020, giaa151bioRxiv, 303727, https://doi.org/10.1093/gigascience/giaa151"

## Frequently Asked Questions

### I'm getting errors from `autoEstCont` or unrealistic estimates

The automatic estimation of the contamination implemented in `autoEstCont` makes the assumption that there is sufficient diversity in the raw data to identify marker genes (as such genes are commonly useful for estimating the contamination).  If your data is either extremely homogenous (i.e., all one cell type, for example a cell line) or your number of cells is very low (a few hundred or less), then this assumption is unlikely to hold.  In such situations you should think hard about if you really want to include data with such severe limitations.  But if you're sure you do, the best approach is probably to manually specify a contamination fraction in line with what you would expect from similar experiments.

### My data still looks contaminated.  Why didn't SoupX work?

The first thing to do is check that you are providing clustering information, either by doing clustering yourself and running `setClusters` before `adjustCounts` or by loading it automatically from `load10X`.  Cluster information allows far more contamination to be identified and safely removed.

The second thing to consider is if the contamination rate estimate looks plausible.  As estimating the contamination rate is the part of the method that requires the most user input, it can be prone to errors. Generally a contamination rate of 2% or less is low, 5% is usual, 10% moderate and 20% or above very high.  Of course your experience may vary and these expectations are based on fresh tissue experiments on the 10X 3' platform.

Finally, note that SoupX has been designed to try and err on the side of not throwing out real counts.  In some cases it is more important to remove contamination than be sure you've retained all the true counts.  This is particularly true as "over-removal" will not remove all the expression from a truly expressed gene unless you set the over-removal to something extreme.  If this describes your situation you may want to try manually increasing the contamination rate by setting `setContaminationFraction` and seeing if this improves your results.

### I can't find a good set of genes to estimate the contamination fraction.

Generally the gene sets that work best are sets of genes highly specific to a cell type that is present in your data at low frequency.  Think HB genes and erythrocytes, IG genes and B-cells, TPSB2/TPSAB1 and Mast cells, etc.  Before trying anything more esoteric, it is usually a good idea to at least try out the most commonly successful gene sets, particularly HB genes.  If this fails, the `plotMarkerDistribution` function can be used to get further inspiration as described in the vignette.  If all of this yields nothing, we suggest trying a range of corrections to see what effect this has on your downstream analysis.  In our experience most experiments have somewhere between 2-10% contamination.

### `estimateNonExpressingCells` can't find any cells to use to estimate contamination.

At this point we assume that you have chosen a set (or sets) of genes to use to estimate the contamination.  The default behaviour (with 10X data) is to look for cells with strong evidence of endogenous expression of these gene sets in all cells, then exclude any cluster with a cell that has strong evidence of endogenous expression.  This conservative behaviour is designed to stop the over-estimation of the contamination fraction, but can sometimes make estimation difficult.  If all clusters have at least one cell that "looks bad" you have 3 options.
1. Recluster the data to produce more clusters with fewer cells per cluster.  This is the preferred option, but requires more work on the users part.
2. Make the criteria for declaring a cell to be genuinely expressing a gene set less strict.  This seldom works, as usually when a cell is over the threshold, it's over by a lot.  But in some cases tweaking the values `maximumContamination` and/or `pCut` can yield usable results.
3. Set `clusters=FALSE` to force `estimateNonExpressingCells` to consider each cell independently.  If you are going to do this, it is worth making the criteria for excluding a cell more permissive by decreasing `maximumContamination` as much as is reasonable.


## Changelog

### v1.5.0

`load10X` now requires the version of `Seurat::Read10X` that does **not** strip out the numeric suffix.

### v1.4.5 

First CRAN version of the code.  The one significant change other than tweaks to reach CRAN compatibility is that the correction algorithm has been made about 20 times faster.  As such, the parallel option was no longer needed and has been removed. Also includes some other minor tweaks.

### v1.3.6

Addition of `autoEstCont` function to automatically estimate the contamination fraction without the need to specify a set of genes to use for estimation.  A number of other tweaks and bug fixes.

### v1.2.1

Some bug fixes from v1.0.0.  Added some helper functions for integrating metadata into SoupChannel object.  Further integration of cluster information in estimation of contamination and calculation of adjusted counts.  Make the `adjustCounts` routine parallel.

### v1.0.0

Review of method, with focus on simplification of code.  Functions that were being used to "automate" selection of genes for contamination estimation have been removed as they were being misused.  Clustering is now used to guide selection of cells where a set of genes is not expressed.  Default now set to use global estimation of rho.  A hierarchical bayes routine has been added to share information between cells when the user does use cell specific estimation.  See NOTE for further details.

### v0.3.0

Now passes R CMD check without warnings or errors.  Added extra vignette on estimating contamination correctly.  Changed the arguments for the interpolateCellContamination function and made monotonically decreasing lowess the default interpolation method.  A number of other plotting improvements.

### v0.2.3

Added lowess smoothing to interpolation and made it the default.  Modified various functions to allow single channel processing in a more natural way.  Some minor bug fixes.

### v0.2.2

Integrated estimateSoup into class construction to save memory when loading many channels.
Added function to use tf-idf to quickly estimate markers.
Some minor bug fixes and documentation updates.

### v0.2.1

Update documentation and modify plot functions to return source data.frame.

### v0.2.0

A fairly major overhaul of the data structures used by the package.  Not compatible with previous versions.

### v0.1.1

Some bug fixes to plotting routines.

# License

```
Copyright (c) 2018 Genome Research Ltd. 
Author: Matthew Young <my4@sanger.ac.uk> 
 
This program is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License version 3 
as published by the Free Software Foundation. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details <http://www.gnu.org/licenses/>. 
```
## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Resubmission

This is a resubmission.  In this version I have:

* Corrected the issue that caused SoupX to be archived.  Specifically, I had not set the LazyDataCompression flag in DESCRIPTION to 'xz', causing the package to fail size checks.  This flag has been set.

* Removed examples from non-exported functions.

---
title: "SoupX PBMC Demonstration"
author: "Matthew Daniel Young"
date: "`r Sys.Date()`"
fig_width: 8
fig_height: 6
output: 
  pdf_document: default
  html_document: default
vignette: >
  %\VignetteIndexEntry{PBMC Demonstration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r global_options, include=FALSE}
library(knitr)
opts_chunk$set(tidy=TRUE)
```

# Introduction

Before we get started with the specifics of example data sets and using the R package, it is worth understanding at a broad level what the problem this package aims to solve is and how it goes about doing it.  Of course, the best way of doing this is by [reading the pre-print](https://www.biorxiv.org/content/10.1101/303727v1), it's not long I promise.  But if you can't be bothered doing that or just want a refresher, I'll try and recap the main points.

In droplet based, single cell RNA-seq experiments, there is always a certain amount of background mRNAs present in the dilution that gets distributed into the droplets with cells and sequenced along with them.  The net effect of this is to produce a background contamination that represents expression not from the cell contained within a droplet, but the solution that contained the cells.

This collection of cell free mRNAs floating in the input solution (henceforth referred to as "the soup") is created from cells in the input solution being lysed.  Because of this, the soup looks different for each input solution and strongly resembles the expression pattern obtained by summing all the individual cells.

The aim of this package is to provide a way to estimate the composition of this soup, what fraction of UMIs are derived from the soup in each droplet and produce a corrected count table with the soup based expression removed.

The method to do this consists of three parts:

1. Calculate the profile of the soup.
2. Estimate the cell specific contamination fraction.
3. Infer a corrected expression matrix. 

In previous versions of SoupX, the estimation of the contamination fraction (step 2) was the part that caused the most difficulty for the user. The contamination fraction is parametrised as `rho` in the code, with `rho=0` meaning no contamination and `rho=1` meaning 100% of UMIs in a droplet are soup. 

From version 1.3.0 onwards, an automated routine for estimating the contamination fraction is provided, which should be suitable is most circumstances.  However, this vignette will still spend a lot of effort explaining how to calculate the contamination fraction "manually".  This is because there are still circumstances where manually estimating `rho` is preferable or the only option and it is important to understanding how the method works and how it can fail.

While it is possible to run SoupX without clustering information, you will get far better results if some basic clustering is provided.  Therefore, it is **strongly** recommended that you provide some clustering information to SoupX.  If you are using 10X data mapped with cellranger, the default clustering produced by cellranger is automatically loaded and used.  The results are not strongly sensitive to the clustering used.  Seurat with default parameters will also yield similar results.

# Quickstart

If you have some 10X data which has been mapped with cellranger, the typical SoupX work flow would be.

```{r quick_start, eval=FALSE}
install.packages('SoupX')
library(SoupX)
#Load data and estimate soup profile
sc = load10X('Path/to/cellranger/outs/folder/')
#Estimate rho
sc = autoEstCont(sc)
#Clean the data
out = adjustCounts(sc)
```

which would produce a new matrix that has had the contaminating reads removed.  This can then be used in any downstream analysis in place of the original matrix.  Note that by default `adjustCounts` will return non-integer counts.  If you require integers for downstream processing either pass out through `round` or set `roundToInt=TRUE` when running `adjustCounts`.

# Getting started

You install this package like any other R package.  The simplest way is to use the CRAN version by running,

```{r install_CRAN,eval=FALSE}
install.packages('SoupX')
```

If you want to use the latest experimental features, you can install the development version from github using the [devtools](https://devtools.r-lib.org/) `install_github` function as follows:

```{r install, eval=FALSE}
devtools::install_github("constantAmateur/SoupX",ref='devel')
```

Once installed, you can load the package in the usual way,

```{r load}
library(SoupX)
```

# PBMC dataset

Like every other single cell tool out there, we are going to use one of the 10X PBMC data sets to demonstrate how to use this package.  Specifically, we will use this [PBMC dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k).  The starting point is to download the [raw](https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz) and [filtered](https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz) cellranger output and extract them to a temporary folder as follows.

```{r download,eval=FALSE}
tmpDir = tempdir(check=TRUE)
download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))
download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))
untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)
untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)
```

## Loading the data

SoupX comes with a convenience function for loading 10X data processed using cellranger.  If you downloaded the data as above you can use it to get started by running,

```{r load_data,eval=FALSE}
sc = load10X(tmpDir)
```

This will load the 10X data into a `SoupChannel` object.  This is just a list with some special properties, storing all the information associated with a single 10X channel.  A `SoupChannel` object can also be created manually by supplying a table of droplets and a table of counts.  Assuming you have followed the above code to download the PBMC data, you could manually construct a `SoupChannel` by running,

```{r load_data_manual,eval=FALSE}
toc = Seurat::Read10X(file.path(tmpDir,'filtered_gene_bc_matrices','GRCh38'))
tod = Seurat::Read10X(file.path(tmpDir,'raw_gene_bc_matrices','GRCh38'))
sc = SoupChannel(tod,toc)
```

To avoid downloading or including large data files, this vignette will use a pre-loaded and processed object `PBMC_sc`.

```{r load_saved}
data(PBMC_sc)
sc = PBMC_sc
sc
```

## Profiling the soup

Having loaded our data, the first thing to do is to estimate what the expression profile of the soup looks like.  This is actually done for us automatically by the object construction function `SoupChannel` called by `load10X`.  Usually, the default estimation is fine, but it can be done explicitly by setting `calcSoupProfile=FALSE` as follows

```{r estimateSoup, eval=FALSE}
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
sc = estimateSoup(sc)
```

Note that we cannot perform this operation using our pre-saved `PBMC_sc` data as the table of droplets is dropped once the soup profile has been generated to save memory.  Generally, we don't need the full table of droplets once we have determined what the soup looks like.

Usually the only reason to not have `estimateSoup` run automatically is if you want to change the default parameters or have some other way of calculating the soup profile.  One case where you may want to do the latter is if you only have the table of counts available and not the empty droplets.  In this case you can proceed by running

```{r estimateNoDrops, eval=TRUE}
library(Matrix)
toc = sc$toc
scNoDrops = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#Calculate soup profile
soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)
```

In this case the `setSoupProfile` command is used instead of `estimateSoup` and directly adds the custom estimation of the soup profile to the `SoupChannel` object.  Note that we have loaded the `Matrix` library to help us manipulate the sparse matrix `toc`.

## Adding extra meta data to the SoupChannel object

We have some extra meta data that it is essential we include in our `SoupChannel` object.  In general you can add any meta data by providing a `data.frame` with row names equal to the column names of the `toc` when building the `SoupChannel` object.  

However, there are some bits of meta data that are so essential that they have their own special loading functions.  The most essential is clustering information.  Without it, SoupX will still work, but you won't be able to automatically estimate the contamination fraction and the correction step will be far less effective.  Metadata associated with our PBMC dataset is also bundled with SoupX.  We can use it to add clustering data by running,

```{r set_clustering}
data(PBMC_metaData)
sc = setClusters(sc,setNames(PBMC_metaData$Cluster,rownames(PBMC_metaData)))
```

It can also be very useful to be able to visualise our data by providing some kind of dimension reduction for the data.  We can do this by running,

```{r add_DR}
sc = setDR(sc,PBMC_metaData[colnames(sc$toc),c('RD1','RD2')])
```

This is usually not needed when using the `load10X` function as the cellranger produced values are automatically loaded.

## Visual sanity checks

It is often the case that really what you want is to get a rough sense of whether the expression of a gene (or group of genes) in a set of cells is derived from the soup or not.  At this stage we already have enough information to do just this.  Before proceeding, we will briefly discuss how to do this.

Let's start by getting a general overview of our PBMC data by plotting it with the provided annotation.

```{r plot_annot}
library(ggplot2)
dd = PBMC_metaData[colnames(sc$toc),]
mids = aggregate(cbind(RD1,RD2) ~ Annotation,data=dd,FUN=mean)
gg = ggplot(dd,aes(RD1,RD2)) + 
  geom_point(aes(colour=Annotation),size=0.2) +
  geom_label(data=mids,aes(label=Annotation)) +
  ggtitle('PBMC 4k Annotation') +
  guides(colour = guide_legend(override.aes = list(size=1)))
plot(gg)
```

SoupX does not have any of its own functions for generating tSNE (or any other reduced dimension) co-ordinates, so it is up to us to generate them using something else.  In this case I have run [Seurat](https://satijalab.org/seurat/) in a standard way and produced a tSNE map of the data (see `?PBMC`).

Suppose that we are interested in the expression of the gene _IGKC_, a key component immunoglobulins (i.e., antibodies) highly expressed by B-cells.  We can quickly visualise which cells express _IGKC_ by extracting the counts for it from the `SoupChannel` object.

```{r plot_IGKC}
dd$IGKC = sc$toc['IGKC',]
gg = ggplot(dd,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)
```

Wow!  We know from prior annotation that the cells in the cluster at the bottom are B-cells so should express _IGKC_.  But the cluster on the right is a T-cell population.  Taken at face value, we appear to have identified a scattered population of T-cells that are producing antibodies!  Start preparing the nature paper!

Before we get too carried away though, perhaps it's worth checking if the expression of _IGKC_ in these scattered cells is more than we would expect by chance from the soup.  To really answer this properly, we need to know how much contamination is present in each cell, which will be the focus of the next sections.  But we can get a rough idea just by calculating how many counts we would expect for _IGKC_ in each cell, by assuming that cell contained nothing but soup.  The function `soupMarkerMap` allows you to visualise the ratio of observed counts for a gene (or set of genes) to this expectation value.  Let's try it out,

```{r sanity_check}
gg = plotMarkerMap(sc,'IGKC')
plot(gg)
```

There is no need to pass the tSNE coordinates to this function as we stored them in the `sc` object when we ran `setDR` above.  Looking at the resulting plot, we see that the cells in the B-cell cluster have a reddish colour, indicating that they are expressed far more than we would expect by chance, even if the cell was nothing but background.  Our paradigm changing, antibody producing T-cells do not fare so well.  They all have a decidedly bluish hue, indicating that is completely plausible that the expression of _IGKC_ in these cells is due to contamination from the soup.  Those cells that are shown as dots have zero expression for _IGKC_.

We have made these plots assuming each droplet contains nothing but background contamination, which is obviously not true.  Nevertheless, this can still be a useful quick and easy sanity check to perform.

## Estimating the contamination fraction

Probably the most difficult part of using SoupX is accurately estimating the level of background contamination (represented as `rho`) in each channel.  There are two ways to do this: using the automatic `autoEstCont` method, or manually providing a list of "non expressed genes".  This vignette will demonstrate both methods, but we anticipate that the automatic method will be used in most circumstances.  Before that we will describe the idea that underpins both approaches; identifying genes that are not expressed by some cells in our data and the expression that we observe for these genes in these cells must be due to contamination.

This is the most challenging part of the method to understand and we have included a lot of detail here. But successfully applying SoupX does not depend on understanding all these details.  The key thing to understand is that the contamination fraction estimate is the fraction of your data that will be discarded.  If this value is set too low, your "corrected" data will potentially still be highly contaminated.  If you set it too high, you will discard real data, although there are good reasons to want to do this at times (see section below).  If the contamination fraction is in the right ball park, SoupX will remove most of the contamination.  It will generally not matter if this number if off by a few percent.

Note that all modes of determining the contamination fraction add an entry titled `fit` to the `SoupChannel` object which contains details of how the final estimate was reached. 

### Manually specifying the contamination fraction

It is worth considering simply manually fixing the contamination fraction at a certain value.  This seems like a bad thing to do intuitively, but there are actually good reasons you might want to.  When the contamination fraction is set too high, true expression will be removed from your data.  However, this is done in such a way that the counts that are most specific to a subset of cells (i.e., good marker genes) will be the absolute last thing to be removed.  Because of this, it can be a sensible thing to set a high contamination fraction for a set of experiments and be confident that the vast majority of the contamination has been removed.

Even when you have a good estimate of the contamination fraction, you may want to set the value used artificially higher.  SoupX has been designed to be conservative in the sense that it errs on the side of retaining true expression at the cost of letting some contamination to creep through.  Our tests show that a well estimated contamination fraction will remove 80-90% of the contamination (i.e. the soup is reduced by an order of magnitude).  For most applications this is sufficient.  However, in cases where complete removal of contamination is essential, it can be worthwhile to increase the contamination fraction used by SoupX to remove a greater portion of the contamination.

Our experiments indicate that adding 5% extra removes 90-95% of the soup, 10% gets rid of 95-98% and 20% removes 99% or more.

Explicitly setting the contamination fraction can be done by running,

```{r set_rho}
sc = setContaminationFraction(sc,0.2)
```

to set the contamination fraction to 20% for all cells.

### Genes to estimate the contamination fraction

To estimate the contamination fraction, we need a set of genes that we know are not expressed in a set of cells, so by measuring how much expression we observe we can infer the contamination fraction. That is, we need a set of genes that we know are not expressed by cells of a certain type, so that in these cells the only source of expression is the soup.  The difficulty is in identifying these sets of genes and the cells in which they can be assumed to be not expressed.

Note that the purpose of this set of genes is to estimate the contamination fraction and nothing else.  These genes play no special role in the actual removal of background associated counts.  They are categorically **not** a list of genes to be removed or anything of that sort.

Furthermore, if no good set of genes can be provided, and/or the automatic method fails, it is reasonable to consider setting the contamination fraction manually to something and seeing how your results are effected.  A contamination rate of around 0.1 is appropriate for many datasets, but of course every experiment is different.

To make this concrete, let us consider an example.  The genes _HBB_,_HBA2_ are both haemoglobin genes and so should only be expressed in red blood cells and nowhere else.  _IGKC_ is an antibody gene produced only by B cells.  Suppose we're estimating the contamination then using two sets of genes: HB genes (_HBB_ and _HBA2_) and IG genes (_IGKC_).  Let's now look at what happens in a few hypothetical cells:

Cell 1 - Is a red blood cell so expresses _HBB_ and _HBA2_, but should not express _IGKC_. For this cell we want to use _IGKC_ to estimate the contamination fraction but not _HBB_,_HBA2_.

Cell 2 - Is a B-Cell so should express _IGKC_, but not _HBB_ or _HBA2_. For this cell we want to use _HBB_ and _HBA2_ to estimate the contamination fraction, but not _IGKC_.

Cell 3 - Is an endothelial cell, so should not express any of _HBB_,_HBA2_ or _IGKC_. So we want to use all three to estimate the contamination fraction.

Basically we are trying to identify in each cell, a set of genes we know the cell does not express so we can estimate the contamination fraction using the expression we do see.

Now obviously the method doesn't know anything about the biology and we haven't told it what's a B cell, a RBC or anything else. There is nothing stopping you supplying that information if you do have it and that will of course give the best results.

But absent this information, the trick is to use the expression level of these genes in each cell to identify when not to use a gene to estimate the contamination fraction.  This is why the best genes for estimating the contamination fraction are those that are highly expressed in the cells that do use them (like HB or IG genes).  Then we can be confident that observing a low level of expression of a set of genes in a cell is due to background contamination, not a low level of mRNA production by the cell.

Given a set of genes that we suspect may be useful, the function `plotMarkerDistribution` can be used to visualise how this gene's expression is distributed across cells. To continue our example:

Cell 1 - The measured expression of _HBB_ and _HBA2_ is 10 times what we'd expect if the droplet was filled with soup, so the method will not use either of these genes to calculate `rho`. On the other hand _IGKC_ is about .05 times the value we'd get for pure soup, so that is used.

Cell 2 - _HBB_/_HBA2_ have values around .05 times the soup. _IGKC_ is off the charts at 100 times what we'd expect in the soup. So the method concludes that this cell is expressing _IGKC_ and so uses only _HBB_/_HBA2_ to estimate `rho`.

Cell 3 - All three are at around .05, so all are used to estimate `rho`.

To prevent accidentally including cells that genuinely express one of the estimation genes, SoupX will by default exclude any cluster where even one gene has evidence that it expresses a gene.  So in the example above, SoupX would not use HB genes to estimate the contamination rate in Cell 1, or any of the cells belonging to the same cluster as Cell 1.  This very conservative behaviour is to prevent over-estimation of the contamination fraction.

Clustering is beyond the scope of SoupX, so must be supplied by the user.  For 10X data mapped using cellranger, SoupX will automatically pull the graph based clustering produced by cellranger and use that by default.

As indicated above, to get a more accurate estimate, groups with a similar biological function are grouped together so they're either used or excluded as a group. This is why the parameter `nonExpressedGeneList` is given as a list. Each entry in the list is a group of genes that are grouped biologically. So in our example we would set it like:

```{r genes1}
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC'))
```

in this example we'd probably want to include other IG genes and Haemoglobin genes even through they're not particularly highly expressed in general, as they should correlate biologically. That is,

```{r genes2}
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))
```

or something similar.

### The automated method

Estimating the contamination fraction using the automated is as simple as running:

```{r auto_est}
sc = autoEstCont(sc)
```

This will produce a mysterious looking plot with two distributions and a red line.  To understand this plot, we need to understand a little bit about what `autoEstCont` is doing.  The basic idea is that it tries to aggregate evidence from many plausible estimators of `rho` and assigns the true contamination fraction to the one that occurs the most often.  The reason this works is that incorrect estimates should have no preferred value, while true estimates should cluster around the same value.  The solid curve shows something like the frequency of different estimates of `rho`, with a red line indicating its peak, which gives the estimate of `rho`.  If you are using the default values for `priorRho` and `priorRhoStdDev` (which you probably should be) you can ignore the dashed line.

For those wanting a more detailed explanation, here is what happens.  First the function tries to identify genes that are very specific to one cluster of cells (using `quickMarkers`).  The determination of how specific is "very specific" is based on the gene's tf-idf value for the cluster it is specific to.  See the `quickMarkers` help or [this](https://constantamateur.github.io/2020-04-10-scDE/) for an explanation of what this means.  The default of `tfidfMin=1` demands that genes by reasonably specific, so if you are getting a low number of genes for estimation you can consider decreasing this value.  This list is further reduced by keeping only genes that are "highly expressed" in the soup (as these give more accurate estimates of `rho`), where highly expressed is controlled by `soupQuantile`.  The default value sounds strict, but in practice many genes with tf-idf over 1 tend to pass it.

Each of these genes is used to independently estimate `rho` in the usual way.  That is, the clusters for which the gene can be confidently said to not be expressed (as determined by `estimateNonExpressingCells`) are used to estimate `rho`.  We could then just create a histogram of these estimates and pick the most common value.  This would work reasonably well, but we can do better by considering that each estimate has uncertainty associated with it.  So instead we calculate the posterior distribution of the contamination fraction for each gene/cluster pair and determine the best estimate of `rho` by finding the peak of the average across all these distributions.  This is what is shown as the solid curve in the above plot.

The posterior distribution is calculated using a Poisson likelihood with a gamma distribution prior, parametrised by its mean `priorRho` and standard deviation `priorRhoStdDev`.  The dotted line in the above plot shows the prior distribution.  The default parameters have been calibrated to be fairly non-specific with a slight preference towards values of rho in the 0% to 10% range which is most commonly seen for fresh (i.e. not nuclear) single cell experiments.

The default values place only a very weak constraint, as can be seen by setting a uniform prior

```{r auto_est_unif_prior}
sc = autoEstCont(sc,priorRhoStdDev=0.3)
```

which gives the same answer up to two significant figures.  Of course you can break things if you start setting strong, badly motivated priors, so please don't do this.

### The manual way

The alternative to the automatic method is to manually specify which sets of genes to use to estimate the contamination fraction.  These genes need to be such that we are as certain they will not be expressed in each cell.  See the section above on "Genes to estimate the contamination fraction" for an example which may make this clearer.

For some experiments, such as solid tissue studies where red cell lysis buffer has been used, it is obvious what genes to use for this purpose.  In the case of bloody solid tissue, haemoglobin genes will be a ubiquitous contaminant and are not actually produced by any cell other than red blood cells in most contexts.  If this is the case, you can skip the next section and proceed straight to estimating contamination.

#### Picking soup specific genes

However, some times it is not obvious in advance which genes are highly specific to just one population of cells.  This is the case with our PBMC data, which is not a solid tissue biopsy and so it is not clear which gene sets to use to estimate the contamination.  In general it is up to the user to pick sensible genes, but there are a few things that can be done to aid in this selection process.  Firstly, the genes that are the most useful are those expressed most highly in the background.  We can check which genes these are by running:

```{r topSoupGenes}
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)
```

Unfortunately most of the most highly expressed genes in this case are ubiquitously expressed (_RPL_/_RPS_ genes or mitochondrial genes).  So we need some further criteria to aid our selection process.

The function `plotMarkerDistribution` is used to visualise the distribution of expression (relative to what would be expected were each cell pure background) across all cells in the data set.  When no geneset is provided, the function will try and guess which genes might be useful.  

```{r inferNonExpressed}
plotMarkerDistribution(sc)
```

The plot shows the distribution of log10 ratios of observed counts to expected if the cell contained nothing but soup.  A guess at which cells definitely express each gene is made and the background contamination is calculated.  The red line shows the global estimate (i.e., assuming the same contamination fraction for all cells) of the contamination fraction using just that gene.  This "guessing" is done using the `quickMarkers` function to find marker genes of each cluster (see "Automatic method" section).  As such, it will fail if no clusters have been provided.

Note that this is a heuristic set of genes that is intended to help develop your biological intuition.  It absolutely **must not** be used to automatically select a set of genes to estimate the background contamination fraction.  For this reason, the function will not return a list of genes.  **If you select the top N genes from this list and use those to estimate the contamination, you will over-estimate the contamination fraction!**

Note too that the decision of what genes to use to estimate the contamination must be made on a channel by channel basis.  We will find that B-cell specific genes are useful for estimating the contamination in this channel.  If we had another channel with only T-cells, these markers would be of no use.

Looking at this plot, we observe that there are two immunoglobulin genes from the constant region (_IGKC_ and _IGHM_) present and they give a consistent estimate of the contamination fraction of around 10% (-1 on the log10 scale).  As we know that it is reasonable to assume that immunoglobulin genes are expressed only in B-cells, we will decide to use their expression in non B-cells to estimate the contamination fraction.

But there's no reason to just use the genes `quickMarkers` flagged for us.  So let's define a list of all the constant immunoglobulin genes, 

```{r igGenes}
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')
```

it doesn't matter if some of these are not expressed in our data, they will then just not contribute to the estimate.

#### Estimating non-expressing cells

Having decided on a set of genes with which to estimate the contamination, we next need to decide which cells genuinely express these genes and should not be used for estimating the contamination, and which do not and should.  This is done as follows,

```{r calculateNullMatrix}
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes),clusters=FALSE)
```

Which produces a matrix indicating which cells (rows) should use which sets of genes (columns) to estimate the contamination.  You will notice that the function returned a warning about cluster information not being provided.  As discussed above, SoupX tries to be conservative and prevents estimation both from cells with high expression of a gene set (`igGenes` in this case) and any cell that falls in the same cluster.  When no clustering information is given, it cannot do this so defaults to just excluding those cells that are obviously not suitable.  We can visualise which cells have been marked to use for estimation,

```{r visNullMatrix}
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)
```

You'll notice that above we set `clusters=FALSE` which stops SoupX from using clustering information.  We provided this information earlier when we ran `setClusters`.  Let's see how things change if we let clustering information be used.

```{r calcNullMatrixWithClustering}
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes))
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)
```

As you can see the set of cells to be used for estimation with the `igGenes` set has decreased.  In this case it makes not much difference, but in general it is better to provide clustering and be conservative.

It is worth noting one final thing about the specification of `nonExpressedGeneList`.  It seems odd that we have specified `nonExpressedGeneList = list(IG=igGenes)` instead of just `nonExpressedGeneList = igGenes`.  This is because `nonExpressedGeneList` expects sets of genes that are biologically related and expected to be present or not present as a set (e.g. IG genes, HB genes).


#### Calculating the contamination fraction

At this point all the hard work has been done.  To estimate the contamination fraction you need only pass your set of genes and which cells in which to use those sets of genes to `calculateContaminationFraction`.

```{r calcContamination}
sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst)
```

This function will modify the `metaData` table of `sc` object to add a table giving the contamination fraction estimate.  This approach gives a contamination fraction very close to the automatic method.

```{r viewCont}
head(sc$metaData)
```

## Correcting expression profile

We have now calculated or set the contamination fraction for each cell and would like to use this to remove the contamination from the original count matrix.  As with estimating the contamination, this procedure is made much more robust by providing clustering information.  This is because there is much more power to separate true expression from contaminating expression when counts are aggregated into clusters.  Furthermore, the process of redistributing corrected counts from the cluster level to individual cells automatically corrects for variation in the cell specific contamination rate (see the paper for details).

We have already loaded clustering information into our `sc` object with `setClusters`, which will be used by default.  So we can just run.

```{r decontaminate}
out = adjustCounts(sc)
```

The recommended mode of operation will produce a non-integer (although still sparse) matrix where the original counts have been corrected for background expression.  See the help, code, and paper for details of how this is done.  

You should not change the method parameter unless you have a strong reason to do so.  When you need integer counts for downstream analyses, setting `roundToInt=TRUE`, stochastically rounds up with probability equal to the fraction part of the number.  For example, if a cell has 1.2 corrected counts it will be assigned a value of 1 80% of the time and 2 20% of the time.

### Investigating changes in expression

Before proceeding let's have a look at what this has done.  We can get a sense for what has been the most strongly decreased by looking at the fraction of cells that were non-zero now set to zero after correction.

```{r mostZeroed}
cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed
```

Notice that a number of the genes on this list are highly specific markers of one cell type or group of cells (_CD74_/_HLA-DRA_ antigen presenting cells, _IGKC_ B-cells) and others came up on our list of potential cell specific genes.  Notice also the presence of the mitochondrial gene _MT-ND3_.

If on the other hand we focus on genes for which there is a quantitative difference,

```{r mostReduced}
tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=20)
```

we find genes associated with metabolism and translation.  This is often the case as mitochondrial genes are over represented in the background compared to cells, presumably as a result of the soup being generated from distressed cells.

### Visualising expression distribution

Way back at the start, we did a quick visualisation to look at how the ratio of _IGKC_ expression to pure soup was distributed.  Now that we've corrected our data, we can see how that compares to our corrected data.  The function `plotChangeMap` can help us with this.  By default it plots the fraction of expression in each cell that has been deemed to be soup and removed.

```{r IGKC_change}
plotChangeMap(sc,out,'IGKC')
```

which shows us that the expression has been heavily decreased in the areas where it was very surprising to observe it before.  

The interpretation of which cells are expressing which genes can change quite dramatically when we correct for soup contamination. You should explore this yourself by plotting a variety of different genes and seeing which change and which do not.  For example, you could look at,

```{r change_plots,eval=FALSE}
plotChangeMap(sc,out,'LYZ')
plotChangeMap(sc,out,'CD74')
plotChangeMap(sc,out,'IL32')
plotChangeMap(sc,out,'TRAC')
plotChangeMap(sc,out,'S100A9')
plotChangeMap(sc,out,'NKG7')
plotChangeMap(sc,out,'GNLY')
plotChangeMap(sc,out,'CD4')
plotChangeMap(sc,out,'CD8A')
```

In general, the changes tend to be largest for genes that are highly expressed but only in a specific context.


## Integrating with downstream tools

Of course, the next thing you'll want to do is to load this corrected expression matrix into some downstream analysis tool and further analyse the data.

The corrected matrix can then be used for any downstream analysis in place of the uncorrected raw matrix. If you are using 10X data and would like to save these final counts out in the same format, you can use the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) `write10xCounts` function like this,

```{r writeOut,eval=FALSE}
DropletUtils:::write10xCounts('./strainedCounts',out)
```

### Loading into Seurat

For loading into Seurat or other R packages, there is no need to save the output of SoupX to disk.  For a single channel, you can simply run the standard Seurat construction step on the output of `adjustCounts`.  That is,

```{r seurat,eval=FALSE}
library(Seurat)
srat = CreateSeuratObject(out)
```

If you have multiple channels you want to process in SoupX then load into Seurat, you need to create a combined matrix with cells from all channels.  Assuming `scs` is a list named by experiment (e.g., `scs = list(Experiment1 = scExp1,Experiment2 = scExp2)`), this can be done by running something like:

```{r seuratMulti,eval=FALSE}
library(Seurat)
srat = list()
for(nom in names(scs)){
  #Clean channel named 'nom'
  tmp = adjustCounts(scs[[nom]])
  #Add experiment name to cell barcodes to make them unique
  colnames(tmp) = paste0(nom,'_',colnames(tmp))
  #Store the result
  srat[[nom]] = tmp
}
#Combine all count matricies into one matrix
srat = do.call(cbind,srat)
srat = CreateSeuratObject(srat)
```
---
title: "SoupX PBMC Demonstration"
author: "Matthew Daniel Young"
date: "`r Sys.Date()`"
fig_width: 8
fig_height: 6
output: 
  pdf_document: default
  html_document: default
vignette: >
  %\VignetteIndexEntry{PBMC Demonstration}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r global_options, include=FALSE}
library(knitr)
opts_chunk$set(tidy=TRUE)
```

# Introduction

Before we get started with the specifics of example data sets and using the R package, it is worth understanding at a broad level what the problem this package aims to solve is and how it goes about doing it.  Of course, the best way of doing this is by [reading the pre-print](https://www.biorxiv.org/content/10.1101/303727v1), it's not long I promise.  But if you can't be bothered doing that or just want a refresher, I'll try and recap the main points.

In droplet based, single cell RNA-seq experiments, there is always a certain amount of background mRNAs present in the dilution that gets distributed into the droplets with cells and sequenced along with them.  The net effect of this is to produce a background contamination that represents expression not from the cell contained within a droplet, but the solution that contained the cells.

This collection of cell free mRNAs floating in the input solution (henceforth referred to as "the soup") is created from cells in the input solution being lysed.  Because of this, the soup looks different for each input solution and strongly resembles the expression pattern obtained by summing all the individual cells.

The aim of this package is to provide a way to estimate the composition of this soup, what fraction of UMIs are derived from the soup in each droplet and produce a corrected count table with the soup based expression removed.

The method to do this consists of three parts:

1. Calculate the profile of the soup.
2. Estimate the cell specific contamination fraction.
3. Infer a corrected expression matrix. 

In previous versions of SoupX, the estimation of the contamination fraction (step 2) was the part that caused the most difficulty for the user. The contamination fraction is parametrised as `rho` in the code, with `rho=0` meaning no contamination and `rho=1` meaning 100% of UMIs in a droplet are soup. 

From version 1.3.0 onwards, an automated routine for estimating the contamination fraction is provided, which should be suitable is most circumstances.  However, this vignette will still spend a lot of effort explaining how to calculate the contamination fraction "manually".  This is because there are still circumstances where manually estimating `rho` is preferable or the only option and it is important to understanding how the method works and how it can fail.

While it is possible to run SoupX without clustering information, you will get far better results if some basic clustering is provided.  Therefore, it is **strongly** recommended that you provide some clustering information to SoupX.  If you are using 10X data mapped with cellranger, the default clustering produced by cellranger is automatically loaded and used.  The results are not strongly sensitive to the clustering used.  Seurat with default parameters will also yield similar results.

# Quickstart

If you have some 10X data which has been mapped with cellranger, the typical SoupX work flow would be.

```{r quick_start, eval=FALSE}
install.packages('SoupX')
library(SoupX)
#Load data and estimate soup profile
sc = load10X('Path/to/cellranger/outs/folder/')
#Estimate rho
sc = autoEstCont(sc)
#Clean the data
out = adjustCounts(sc)
```

which would produce a new matrix that has had the contaminating reads removed.  This can then be used in any downstream analysis in place of the original matrix.  Note that by default `adjustCounts` will return non-integer counts.  If you require integers for downstream processing either pass out through `round` or set `roundToInt=TRUE` when running `adjustCounts`.

# Getting started

You install this package like any other R package.  The simplest way is to use the CRAN version by running,

```{r install_CRAN,eval=FALSE}
install.packages('SoupX')
```

If you want to use the latest experimental features, you can install the development version from github using the [devtools](https://devtools.r-lib.org/) `install_github` function as follows:

```{r install, eval=FALSE}
devtools::install_github("constantAmateur/SoupX",ref='devel')
```

Once installed, you can load the package in the usual way,

```{r load}
library(SoupX)
```

# PBMC dataset

Like every other single cell tool out there, we are going to use one of the 10X PBMC data sets to demonstrate how to use this package.  Specifically, we will use this [PBMC dataset](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k).  The starting point is to download the [raw](https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz) and [filtered](https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz) cellranger output and extract them to a temporary folder as follows.

```{r download,eval=FALSE}
tmpDir = tempdir(check=TRUE)
download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))
download.file('https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))
untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)
untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)
```

## Loading the data

SoupX comes with a convenience function for loading 10X data processed using cellranger.  If you downloaded the data as above you can use it to get started by running,

```{r load_data,eval=FALSE}
sc = load10X(tmpDir)
```

This will load the 10X data into a `SoupChannel` object.  This is just a list with some special properties, storing all the information associated with a single 10X channel.  A `SoupChannel` object can also be created manually by supplying a table of droplets and a table of counts.  Assuming you have followed the above code to download the PBMC data, you could manually construct a `SoupChannel` by running,

```{r load_data_manual,eval=FALSE}
toc = Seurat::Read10X(file.path(tmpDir,'filtered_gene_bc_matrices','GRCh38'))
tod = Seurat::Read10X(file.path(tmpDir,'raw_gene_bc_matrices','GRCh38'))
sc = SoupChannel(tod,toc)
```

To avoid downloading or including large data files, this vignette will use a pre-loaded and processed object `PBMC_sc`.

```{r load_saved}
data(PBMC_sc)
sc = PBMC_sc
sc
```

## Profiling the soup

Having loaded our data, the first thing to do is to estimate what the expression profile of the soup looks like.  This is actually done for us automatically by the object construction function `SoupChannel` called by `load10X`.  Usually, the default estimation is fine, but it can be done explicitly by setting `calcSoupProfile=FALSE` as follows

```{r estimateSoup, eval=FALSE}
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
sc = estimateSoup(sc)
```

Note that we cannot perform this operation using our pre-saved `PBMC_sc` data as the table of droplets is dropped once the soup profile has been generated to save memory.  Generally, we don't need the full table of droplets once we have determined what the soup looks like.

Usually the only reason to not have `estimateSoup` run automatically is if you want to change the default parameters or have some other way of calculating the soup profile.  One case where you may want to do the latter is if you only have the table of counts available and not the empty droplets.  In this case you can proceed by running

```{r estimateNoDrops, eval=TRUE}
library(Matrix)
toc = sc$toc
scNoDrops = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#Calculate soup profile
soupProf = data.frame(row.names = rownames(toc),
                      est = rowSums(toc)/sum(toc),
                      counts = rowSums(toc))
scNoDrops = setSoupProfile(scNoDrops,soupProf)
```

In this case the `setSoupProfile` command is used instead of `estimateSoup` and directly adds the custom estimation of the soup profile to the `SoupChannel` object.  Note that we have loaded the `Matrix` library to help us manipulate the sparse matrix `toc`.

## Adding extra meta data to the SoupChannel object

We have some extra meta data that it is essential we include in our `SoupChannel` object.  In general you can add any meta data by providing a `data.frame` with row names equal to the column names of the `toc` when building the `SoupChannel` object.  

However, there are some bits of meta data that are so essential that they have their own special loading functions.  The most essential is clustering information.  Without it, SoupX will still work, but you won't be able to automatically estimate the contamination fraction and the correction step will be far less effective.  Metadata associated with our PBMC dataset is also bundled with SoupX.  We can use it to add clustering data by running,

```{r set_clustering}
data(PBMC_metaData)
sc = setClusters(sc,setNames(PBMC_metaData$Cluster,rownames(PBMC_metaData)))
```

It can also be very useful to be able to visualise our data by providing some kind of dimension reduction for the data.  We can do this by running,

```{r add_DR}
sc = setDR(sc,PBMC_metaData[colnames(sc$toc),c('RD1','RD2')])
```

This is usually not needed when using the `load10X` function as the cellranger produced values are automatically loaded.

## Visual sanity checks

It is often the case that really what you want is to get a rough sense of whether the expression of a gene (or group of genes) in a set of cells is derived from the soup or not.  At this stage we already have enough information to do just this.  Before proceeding, we will briefly discuss how to do this.

Let's start by getting a general overview of our PBMC data by plotting it with the provided annotation.

```{r plot_annot}
library(ggplot2)
dd = PBMC_metaData[colnames(sc$toc),]
mids = aggregate(cbind(RD1,RD2) ~ Annotation,data=dd,FUN=mean)
gg = ggplot(dd,aes(RD1,RD2)) + 
  geom_point(aes(colour=Annotation),size=0.2) +
  geom_label(data=mids,aes(label=Annotation)) +
  ggtitle('PBMC 4k Annotation') +
  guides(colour = guide_legend(override.aes = list(size=1)))
plot(gg)
```

SoupX does not have any of its own functions for generating tSNE (or any other reduced dimension) co-ordinates, so it is up to us to generate them using something else.  In this case I have run [Seurat](https://satijalab.org/seurat/) in a standard way and produced a tSNE map of the data (see `?PBMC`).

Suppose that we are interested in the expression of the gene _IGKC_, a key component immunoglobulins (i.e., antibodies) highly expressed by B-cells.  We can quickly visualise which cells express _IGKC_ by extracting the counts for it from the `SoupChannel` object.

```{r plot_IGKC}
dd$IGKC = sc$toc['IGKC',]
gg = ggplot(dd,aes(RD1,RD2)) +
  geom_point(aes(colour=IGKC>0))
plot(gg)
```

Wow!  We know from prior annotation that the cells in the cluster at the bottom are B-cells so should express _IGKC_.  But the cluster on the right is a T-cell population.  Taken at face value, we appear to have identified a scattered population of T-cells that are producing antibodies!  Start preparing the nature paper!

Before we get too carried away though, perhaps it's worth checking if the expression of _IGKC_ in these scattered cells is more than we would expect by chance from the soup.  To really answer this properly, we need to know how much contamination is present in each cell, which will be the focus of the next sections.  But we can get a rough idea just by calculating how many counts we would expect for _IGKC_ in each cell, by assuming that cell contained nothing but soup.  The function `soupMarkerMap` allows you to visualise the ratio of observed counts for a gene (or set of genes) to this expectation value.  Let's try it out,

```{r sanity_check}
gg = plotMarkerMap(sc,'IGKC')
plot(gg)
```

There is no need to pass the tSNE coordinates to this function as we stored them in the `sc` object when we ran `setDR` above.  Looking at the resulting plot, we see that the cells in the B-cell cluster have a reddish colour, indicating that they are expressed far more than we would expect by chance, even if the cell was nothing but background.  Our paradigm changing, antibody producing T-cells do not fare so well.  They all have a decidedly bluish hue, indicating that is completely plausible that the expression of _IGKC_ in these cells is due to contamination from the soup.  Those cells that are shown as dots have zero expression for _IGKC_.

We have made these plots assuming each droplet contains nothing but background contamination, which is obviously not true.  Nevertheless, this can still be a useful quick and easy sanity check to perform.

## Estimating the contamination fraction

Probably the most difficult part of using SoupX is accurately estimating the level of background contamination (represented as `rho`) in each channel.  There are two ways to do this: using the automatic `autoEstCont` method, or manually providing a list of "non expressed genes".  This vignette will demonstrate both methods, but we anticipate that the automatic method will be used in most circumstances.  Before that we will describe the idea that underpins both approaches; identifying genes that are not expressed by some cells in our data and the expression that we observe for these genes in these cells must be due to contamination.

This is the most challenging part of the method to understand and we have included a lot of detail here. But successfully applying SoupX does not depend on understanding all these details.  The key thing to understand is that the contamination fraction estimate is the fraction of your data that will be discarded.  If this value is set too low, your "corrected" data will potentially still be highly contaminated.  If you set it too high, you will discard real data, although there are good reasons to want to do this at times (see section below).  If the contamination fraction is in the right ball park, SoupX will remove most of the contamination.  It will generally not matter if this number if off by a few percent.

Note that all modes of determining the contamination fraction add an entry titled `fit` to the `SoupChannel` object which contains details of how the final estimate was reached. 

### Manually specifying the contamination fraction

It is worth considering simply manually fixing the contamination fraction at a certain value.  This seems like a bad thing to do intuitively, but there are actually good reasons you might want to.  When the contamination fraction is set too high, true expression will be removed from your data.  However, this is done in such a way that the counts that are most specific to a subset of cells (i.e., good marker genes) will be the absolute last thing to be removed.  Because of this, it can be a sensible thing to set a high contamination fraction for a set of experiments and be confident that the vast majority of the contamination has been removed.

Even when you have a good estimate of the contamination fraction, you may want to set the value used artificially higher.  SoupX has been designed to be conservative in the sense that it errs on the side of retaining true expression at the cost of letting some contamination to creep through.  Our tests show that a well estimated contamination fraction will remove 80-90% of the contamination (i.e. the soup is reduced by an order of magnitude).  For most applications this is sufficient.  However, in cases where complete removal of contamination is essential, it can be worthwhile to increase the contamination fraction used by SoupX to remove a greater portion of the contamination.

Our experiments indicate that adding 5% extra removes 90-95% of the soup, 10% gets rid of 95-98% and 20% removes 99% or more.

Explicitly setting the contamination fraction can be done by running,

```{r set_rho}
sc = setContaminationFraction(sc,0.2)
```

to set the contamination fraction to 20% for all cells.

### Genes to estimate the contamination fraction

To estimate the contamination fraction, we need a set of genes that we know are not expressed in a set of cells, so by measuring how much expression we observe we can infer the contamination fraction. That is, we need a set of genes that we know are not expressed by cells of a certain type, so that in these cells the only source of expression is the soup.  The difficulty is in identifying these sets of genes and the cells in which they can be assumed to be not expressed.

Note that the purpose of this set of genes is to estimate the contamination fraction and nothing else.  These genes play no special role in the actual removal of background associated counts.  They are categorically **not** a list of genes to be removed or anything of that sort.

Furthermore, if no good set of genes can be provided, and/or the automatic method fails, it is reasonable to consider setting the contamination fraction manually to something and seeing how your results are effected.  A contamination rate of around 0.1 is appropriate for many datasets, but of course every experiment is different.

To make this concrete, let us consider an example.  The genes _HBB_,_HBA2_ are both haemoglobin genes and so should only be expressed in red blood cells and nowhere else.  _IGKC_ is an antibody gene produced only by B cells.  Suppose we're estimating the contamination then using two sets of genes: HB genes (_HBB_ and _HBA2_) and IG genes (_IGKC_).  Let's now look at what happens in a few hypothetical cells:

Cell 1 - Is a red blood cell so expresses _HBB_ and _HBA2_, but should not express _IGKC_. For this cell we want to use _IGKC_ to estimate the contamination fraction but not _HBB_,_HBA2_.

Cell 2 - Is a B-Cell so should express _IGKC_, but not _HBB_ or _HBA2_. For this cell we want to use _HBB_ and _HBA2_ to estimate the contamination fraction, but not _IGKC_.

Cell 3 - Is an endothelial cell, so should not express any of _HBB_,_HBA2_ or _IGKC_. So we want to use all three to estimate the contamination fraction.

Basically we are trying to identify in each cell, a set of genes we know the cell does not express so we can estimate the contamination fraction using the expression we do see.

Now obviously the method doesn't know anything about the biology and we haven't told it what's a B cell, a RBC or anything else. There is nothing stopping you supplying that information if you do have it and that will of course give the best results.

But absent this information, the trick is to use the expression level of these genes in each cell to identify when not to use a gene to estimate the contamination fraction.  This is why the best genes for estimating the contamination fraction are those that are highly expressed in the cells that do use them (like HB or IG genes).  Then we can be confident that observing a low level of expression of a set of genes in a cell is due to background contamination, not a low level of mRNA production by the cell.

Given a set of genes that we suspect may be useful, the function `plotMarkerDistribution` can be used to visualise how this gene's expression is distributed across cells. To continue our example:

Cell 1 - The measured expression of _HBB_ and _HBA2_ is 10 times what we'd expect if the droplet was filled with soup, so the method will not use either of these genes to calculate `rho`. On the other hand _IGKC_ is about .05 times the value we'd get for pure soup, so that is used.

Cell 2 - _HBB_/_HBA2_ have values around .05 times the soup. _IGKC_ is off the charts at 100 times what we'd expect in the soup. So the method concludes that this cell is expressing _IGKC_ and so uses only _HBB_/_HBA2_ to estimate `rho`.

Cell 3 - All three are at around .05, so all are used to estimate `rho`.

To prevent accidentally including cells that genuinely express one of the estimation genes, SoupX will by default exclude any cluster where even one gene has evidence that it expresses a gene.  So in the example above, SoupX would not use HB genes to estimate the contamination rate in Cell 1, or any of the cells belonging to the same cluster as Cell 1.  This very conservative behaviour is to prevent over-estimation of the contamination fraction.

Clustering is beyond the scope of SoupX, so must be supplied by the user.  For 10X data mapped using cellranger, SoupX will automatically pull the graph based clustering produced by cellranger and use that by default.

As indicated above, to get a more accurate estimate, groups with a similar biological function are grouped together so they're either used or excluded as a group. This is why the parameter `nonExpressedGeneList` is given as a list. Each entry in the list is a group of genes that are grouped biologically. So in our example we would set it like:

```{r genes1}
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC'))
```

in this example we'd probably want to include other IG genes and Haemoglobin genes even through they're not particularly highly expressed in general, as they should correlate biologically. That is,

```{r genes2}
nonExpressedGeneList = list(HB=c('HBB','HBA2'),IG = c('IGKC','IGHG1','IGHG3'))
```

or something similar.

### The automated method

Estimating the contamination fraction using the automated is as simple as running:

```{r auto_est}
sc = autoEstCont(sc)
```

This will produce a mysterious looking plot with two distributions and a red line.  To understand this plot, we need to understand a little bit about what `autoEstCont` is doing.  The basic idea is that it tries to aggregate evidence from many plausible estimators of `rho` and assigns the true contamination fraction to the one that occurs the most often.  The reason this works is that incorrect estimates should have no preferred value, while true estimates should cluster around the same value.  The solid curve shows something like the frequency of different estimates of `rho`, with a red line indicating its peak, which gives the estimate of `rho`.  If you are using the default values for `priorRho` and `priorRhoStdDev` (which you probably should be) you can ignore the dashed line.

For those wanting a more detailed explanation, here is what happens.  First the function tries to identify genes that are very specific to one cluster of cells (using `quickMarkers`).  The determination of how specific is "very specific" is based on the gene's tf-idf value for the cluster it is specific to.  See the `quickMarkers` help or [this](https://constantamateur.github.io/2020-04-10-scDE/) for an explanation of what this means.  The default of `tfidfMin=1` demands that genes by reasonably specific, so if you are getting a low number of genes for estimation you can consider decreasing this value.  This list is further reduced by keeping only genes that are "highly expressed" in the soup (as these give more accurate estimates of `rho`), where highly expressed is controlled by `soupQuantile`.  The default value sounds strict, but in practice many genes with tf-idf over 1 tend to pass it.

Each of these genes is used to independently estimate `rho` in the usual way.  That is, the clusters for which the gene can be confidently said to not be expressed (as determined by `estimateNonExpressingCells`) are used to estimate `rho`.  We could then just create a histogram of these estimates and pick the most common value.  This would work reasonably well, but we can do better by considering that each estimate has uncertainty associated with it.  So instead we calculate the posterior distribution of the contamination fraction for each gene/cluster pair and determine the best estimate of `rho` by finding the peak of the average across all these distributions.  This is what is shown as the solid curve in the above plot.

The posterior distribution is calculated using a Poisson likelihood with a gamma distribution prior, parametrised by its mean `priorRho` and standard deviation `priorRhoStdDev`.  The dotted line in the above plot shows the prior distribution.  The default parameters have been calibrated to be fairly non-specific with a slight preference towards values of rho in the 0% to 10% range which is most commonly seen for fresh (i.e. not nuclear) single cell experiments.

The default values place only a very weak constraint, as can be seen by setting a uniform prior

```{r auto_est_unif_prior}
sc = autoEstCont(sc,priorRhoStdDev=0.3)
```

which gives the same answer up to two significant figures.  Of course you can break things if you start setting strong, badly motivated priors, so please don't do this.

### The manual way

The alternative to the automatic method is to manually specify which sets of genes to use to estimate the contamination fraction.  These genes need to be such that we are as certain they will not be expressed in each cell.  See the section above on "Genes to estimate the contamination fraction" for an example which may make this clearer.

For some experiments, such as solid tissue studies where red cell lysis buffer has been used, it is obvious what genes to use for this purpose.  In the case of bloody solid tissue, haemoglobin genes will be a ubiquitous contaminant and are not actually produced by any cell other than red blood cells in most contexts.  If this is the case, you can skip the next section and proceed straight to estimating contamination.

#### Picking soup specific genes

However, some times it is not obvious in advance which genes are highly specific to just one population of cells.  This is the case with our PBMC data, which is not a solid tissue biopsy and so it is not clear which gene sets to use to estimate the contamination.  In general it is up to the user to pick sensible genes, but there are a few things that can be done to aid in this selection process.  Firstly, the genes that are the most useful are those expressed most highly in the background.  We can check which genes these are by running:

```{r topSoupGenes}
head(sc$soupProfile[order(sc$soupProfile$est,decreasing=TRUE),],n=20)
```

Unfortunately most of the most highly expressed genes in this case are ubiquitously expressed (_RPL_/_RPS_ genes or mitochondrial genes).  So we need some further criteria to aid our selection process.

The function `plotMarkerDistribution` is used to visualise the distribution of expression (relative to what would be expected were each cell pure background) across all cells in the data set.  When no geneset is provided, the function will try and guess which genes might be useful.  

```{r inferNonExpressed}
plotMarkerDistribution(sc)
```

The plot shows the distribution of log10 ratios of observed counts to expected if the cell contained nothing but soup.  A guess at which cells definitely express each gene is made and the background contamination is calculated.  The red line shows the global estimate (i.e., assuming the same contamination fraction for all cells) of the contamination fraction using just that gene.  This "guessing" is done using the `quickMarkers` function to find marker genes of each cluster (see "Automatic method" section).  As such, it will fail if no clusters have been provided.

Note that this is a heuristic set of genes that is intended to help develop your biological intuition.  It absolutely **must not** be used to automatically select a set of genes to estimate the background contamination fraction.  For this reason, the function will not return a list of genes.  **If you select the top N genes from this list and use those to estimate the contamination, you will over-estimate the contamination fraction!**

Note too that the decision of what genes to use to estimate the contamination must be made on a channel by channel basis.  We will find that B-cell specific genes are useful for estimating the contamination in this channel.  If we had another channel with only T-cells, these markers would be of no use.

Looking at this plot, we observe that there are two immunoglobulin genes from the constant region (_IGKC_ and _IGHM_) present and they give a consistent estimate of the contamination fraction of around 10% (-1 on the log10 scale).  As we know that it is reasonable to assume that immunoglobulin genes are expressed only in B-cells, we will decide to use their expression in non B-cells to estimate the contamination fraction.

But there's no reason to just use the genes `quickMarkers` flagged for us.  So let's define a list of all the constant immunoglobulin genes, 

```{r igGenes}
igGenes = c('IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM',
            'IGLC1','IGLC2','IGLC3','IGLC4','IGLC5','IGLC6','IGLC7',
            'IGKC')
```

it doesn't matter if some of these are not expressed in our data, they will then just not contribute to the estimate.

#### Estimating non-expressing cells

Having decided on a set of genes with which to estimate the contamination, we next need to decide which cells genuinely express these genes and should not be used for estimating the contamination, and which do not and should.  This is done as follows,

```{r calculateNullMatrix}
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes),clusters=FALSE)
```

Which produces a matrix indicating which cells (rows) should use which sets of genes (columns) to estimate the contamination.  You will notice that the function returned a warning about cluster information not being provided.  As discussed above, SoupX tries to be conservative and prevents estimation both from cells with high expression of a gene set (`igGenes` in this case) and any cell that falls in the same cluster.  When no clustering information is given, it cannot do this so defaults to just excluding those cells that are obviously not suitable.  We can visualise which cells have been marked to use for estimation,

```{r visNullMatrix}
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)
```

You'll notice that above we set `clusters=FALSE` which stops SoupX from using clustering information.  We provided this information earlier when we ran `setClusters`.  Let's see how things change if we let clustering information be used.

```{r calcNullMatrixWithClustering}
useToEst = estimateNonExpressingCells(sc,nonExpressedGeneList = list(IG=igGenes))
plotMarkerMap(sc,geneSet=igGenes,useToEst=useToEst)
```

As you can see the set of cells to be used for estimation with the `igGenes` set has decreased.  In this case it makes not much difference, but in general it is better to provide clustering and be conservative.

It is worth noting one final thing about the specification of `nonExpressedGeneList`.  It seems odd that we have specified `nonExpressedGeneList = list(IG=igGenes)` instead of just `nonExpressedGeneList = igGenes`.  This is because `nonExpressedGeneList` expects sets of genes that are biologically related and expected to be present or not present as a set (e.g. IG genes, HB genes).


#### Calculating the contamination fraction

At this point all the hard work has been done.  To estimate the contamination fraction you need only pass your set of genes and which cells in which to use those sets of genes to `calculateContaminationFraction`.

```{r calcContamination}
sc = calculateContaminationFraction(sc,list(IG=igGenes),useToEst=useToEst)
```

This function will modify the `metaData` table of `sc` object to add a table giving the contamination fraction estimate.  This approach gives a contamination fraction very close to the automatic method.

```{r viewCont}
head(sc$metaData)
```

## Correcting expression profile

We have now calculated or set the contamination fraction for each cell and would like to use this to remove the contamination from the original count matrix.  As with estimating the contamination, this procedure is made much more robust by providing clustering information.  This is because there is much more power to separate true expression from contaminating expression when counts are aggregated into clusters.  Furthermore, the process of redistributing corrected counts from the cluster level to individual cells automatically corrects for variation in the cell specific contamination rate (see the paper for details).

We have already loaded clustering information into our `sc` object with `setClusters`, which will be used by default.  So we can just run.

```{r decontaminate}
out = adjustCounts(sc)
```

The recommended mode of operation will produce a non-integer (although still sparse) matrix where the original counts have been corrected for background expression.  See the help, code, and paper for details of how this is done.  

You should not change the method parameter unless you have a strong reason to do so.  When you need integer counts for downstream analyses, setting `roundToInt=TRUE`, stochastically rounds up with probability equal to the fraction part of the number.  For example, if a cell has 1.2 corrected counts it will be assigned a value of 1 80% of the time and 2 20% of the time.

### Investigating changes in expression

Before proceeding let's have a look at what this has done.  We can get a sense for what has been the most strongly decreased by looking at the fraction of cells that were non-zero now set to zero after correction.

```{r mostZeroed}
cntSoggy = rowSums(sc$toc>0)
cntStrained = rowSums(out>0)
mostZeroed = tail(sort((cntSoggy-cntStrained)/cntSoggy),n=10)
mostZeroed
```

Notice that a number of the genes on this list are highly specific markers of one cell type or group of cells (_CD74_/_HLA-DRA_ antigen presenting cells, _IGKC_ B-cells) and others came up on our list of potential cell specific genes.  Notice also the presence of the mitochondrial gene _MT-ND3_.

If on the other hand we focus on genes for which there is a quantitative difference,

```{r mostReduced}
tail(sort(rowSums(sc$toc>out)/rowSums(sc$toc>0)),n=20)
```

we find genes associated with metabolism and translation.  This is often the case as mitochondrial genes are over represented in the background compared to cells, presumably as a result of the soup being generated from distressed cells.

### Visualising expression distribution

Way back at the start, we did a quick visualisation to look at how the ratio of _IGKC_ expression to pure soup was distributed.  Now that we've corrected our data, we can see how that compares to our corrected data.  The function `plotChangeMap` can help us with this.  By default it plots the fraction of expression in each cell that has been deemed to be soup and removed.

```{r IGKC_change}
plotChangeMap(sc,out,'IGKC')
```

which shows us that the expression has been heavily decreased in the areas where it was very surprising to observe it before.  

The interpretation of which cells are expressing which genes can change quite dramatically when we correct for soup contamination. You should explore this yourself by plotting a variety of different genes and seeing which change and which do not.  For example, you could look at,

```{r change_plots,eval=FALSE}
plotChangeMap(sc,out,'LYZ')
plotChangeMap(sc,out,'CD74')
plotChangeMap(sc,out,'IL32')
plotChangeMap(sc,out,'TRAC')
plotChangeMap(sc,out,'S100A9')
plotChangeMap(sc,out,'NKG7')
plotChangeMap(sc,out,'GNLY')
plotChangeMap(sc,out,'CD4')
plotChangeMap(sc,out,'CD8A')
```

In general, the changes tend to be largest for genes that are highly expressed but only in a specific context.


## Integrating with downstream tools

Of course, the next thing you'll want to do is to load this corrected expression matrix into some downstream analysis tool and further analyse the data.

The corrected matrix can then be used for any downstream analysis in place of the uncorrected raw matrix. If you are using 10X data and would like to save these final counts out in the same format, you can use the [DropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) `write10xCounts` function like this,

```{r writeOut,eval=FALSE}
DropletUtils:::write10xCounts('./strainedCounts',out)
```

### Loading into Seurat

For loading into Seurat or other R packages, there is no need to save the output of SoupX to disk.  For a single channel, you can simply run the standard Seurat construction step on the output of `adjustCounts`.  That is,

```{r seurat,eval=FALSE}
library(Seurat)
srat = CreateSeuratObject(out)
```

If you have multiple channels you want to process in SoupX then load into Seurat, you need to create a combined matrix with cells from all channels.  Assuming `scs` is a list named by experiment (e.g., `scs = list(Experiment1 = scExp1,Experiment2 = scExp2)`), this can be done by running something like:

```{r seuratMulti,eval=FALSE}
library(Seurat)
srat = list()
for(nom in names(scs)){
  #Clean channel named 'nom'
  tmp = adjustCounts(scs[[nom]])
  #Add experiment name to cell barcodes to make them unique
  colnames(tmp) = paste0(nom,'_',colnames(tmp))
  #Store the result
  srat[[nom]] = tmp
}
#Combine all count matricies into one matrix
srat = do.call(cbind,srat)
srat = CreateSeuratObject(srat)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateNonExpressingCells.R
\name{estimateNonExpressingCells}
\alias{estimateNonExpressingCells}
\title{Calculate which cells genuinely do not express a particular gene or set of genes}
\usage{
estimateNonExpressingCells(
  sc,
  nonExpressedGeneList,
  clusters = NULL,
  maximumContamination = 1,
  FDR = 0.05
)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{nonExpressedGeneList}{A list containing sets of genes which will be used to estimate the contamination fraction.}

\item{clusters}{A named vector indicating how to cluster cells.  Names should be cell IDs, values cluster IDs.  If NULL, we will attempt to load it from sc$metaData$clusters.  If set to FALSE, each cell will be considered individually.}

\item{maximumContamination}{The maximum contamination fraction that you would reasonably expect.  The lower this value is set, the more aggressively cells are excluded from use in estimation.}

\item{FDR}{A Poisson test is used to identify cells to exclude, this is the false discovery rate it uses.  Higher FDR = more aggressive exclusion.}
}
\value{
A matrix indicating which cells to be used to estimate contamination for each set of genes.  Typically passed to the \code{useToEst} parameter of \code{\link{calculateContaminationFraction}} or \code{\link{plotMarkerMap}}.
}
\description{
Given a list of correlated genes (e.g. Haemoglobin genes, Immunoglobulin genes, etc.), make an attempt to estimate which cells genuinely do not express each of these gene sets in turn.  The central idea is that in cells that are not genuinely producing a class of mRNAs (such as haemoglobin genes), any observed expression of these genes must be due to ambient RNA contamination.  As such, if we can identify these cells, we can use the observed level of expression of these genes to estimate the level of contamination.
}
\details{
The ideal way to do this would be to have a prior annotation of your data indicating which cells are (for instance) red blood cells and genuinely expression haemoglobin genes, and which do not and so only express haemoglobin genes due to contamination.  If this is your circumstance, there is no need to run this function, you can instead pass a matrix encoding which cells are haemoglobin expressing and which are not to \code{\link{calculateContaminationFraction}} via the \code{useToEst} parameter.

This function will use a conservative approach to excluding cells that it thinks may express one of your gene sets.  This is because falsely including a cell in the set of non-expressing cells may erroneously inflate your estimated contamination, whereas failing to include a genuine non-expressing cell in this set has no significant effect.

To this end, this function will exclude any cluster of cells in which any cell is deemed to have genuine expression of a gene set.  Clustering of data is beyond the scope of this package, but can be performed by the user.  In the case of 10X data mapped using cellranger and loaded using \code{\link{load10X}}, the cellranger graph based clustering is automatically loaded and used.

To decide if a cell is genuinely expressing a set of genes, a Poisson test is used.  This tests whether the observed expression is greater than \code{maximumContamination} times the expected number of counts for a set of genes, if the cell were assumed to be derived wholly from the background.  This process can be made less conservative (i.e., excluding fewer cells/clusters) by either decreasing the value of the maximum contamination the user believes is plausible (\code{maximumContamination}) or making the significance threshold for the test more strict (by reducing \code{FDR}).
}
\examples{
#Common gene list in real world data
geneList = list(HB=c('HBB','HBA2'))
#Gene list appropriate to toy data
geneList = list(CD7 = 'CD7')
ute = estimateNonExpressingCells(scToy,geneList)
}
\seealso{
calculateContaminationFraction plotMarkerMap
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotSoupCorrelation}
\alias{plotSoupCorrelation}
\title{Plot correlation of expression profiles of soup and aggregated cells}
\usage{
plotSoupCorrelation(sc)
}
\arguments{
\item{sc}{A SoupChannel object.}
}
\value{
A ggplot2 object containing the plot.
}
\description{
Calculates an expression profile by aggregating counts across all cells and plots this (on a log10 scale) against the expression profile of the soup.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{scToy}
\alias{scToy}
\title{Toy SoupChanel object}
\format{
\code{scToy} is a \code{SoupChannel} object.
}
\usage{
data(scToy)
}
\description{
A \code{\link{SoupChannel}} object created from the toy data used in examples.
}
\details{
The toy data is created from a modified version of the extremely reduced \code{Seurat} \code{pbmc_small} dataset.  It includes clusters, tSNE coordinates and a flat estimate of 0.1 contamination.  It includes data for only 226 genes and 62 cells and should not be used for anything other than testing functions as it is not representative of real data in any way.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PBMC_metaData}
\alias{PBMC_metaData}
\title{PBMC 4K meta data}
\format{
\code{PBMC_metaData} is a data.frame with 4 columns: RD1, RD2, Cluster, and Annotation.
}
\source{
\url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
}
\usage{
data(PBMC_metaData)
}
\description{
Collection of bits of meta data relating to the 10X PBMC 4K data.
}
\details{
This data set pertains to the 10X demonstration PBMC 4K data and includes metadata about it in the \code{data.frame} named \code{PBMC_metaData}.

\code{PBMC_metaData} was created using Seurat (v2) to calculate a tSNE representation of the data and cluster cells with these commands.
\itemize{
  \item \code{set.seed(1)}
  \item \code{srat = CreateSeuratObject(sc$toc)}
  \item \code{srat = NormalizeData(srat)}
  \item \code{srat = ScaleData(srat)}
  \item \code{srat = FindVariableGenes(srat)}
  \item \code{srat = RunPCA(srat,pcs.compute=30)}
  \item \code{srat = RunTSNE(srat,dims.use=seq(30))}
  \item \code{srat = FindClusters(srat,dims.use=seq(30),resolution=1)}
  \item \code{PBMC_metaData = as.data.frame(srat@dr$tsne@cell.embeddings)}
  \item \code{colnames(PBMC_metaData) = c('RD1','RD2')}
  \item \code{PBMC_metaData$Cluster = factor(srat@meta.data[rownames(PBMC_metaData),'res.1'])}
  \item \code{PBMC_metaData$Annotation = factor(c('7'='B','4'='B','1'='T_CD4','2'='T_CD4','3'='T_CD8','5'='T_CD8','6'='NK','8'='NK','0'='MNP','9'='MNP','10'='MNP','11'='?')[as.character(PBMC_metaData$Cluster)])}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotMarkerDistribution}
\alias{plotMarkerDistribution}
\title{Plots the distribution of the observed to expected expression for marker genes}
\usage{
plotMarkerDistribution(
  sc,
  nonExpressedGeneList,
  maxCells = 150,
  tfidfMin = 1,
  ...
)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{nonExpressedGeneList}{Which sets of genes to use to estimate soup (see \code{\link{calculateContaminationFraction}}).}

\item{maxCells}{Randomly plot only this many cells to prevent over-crowding.}

\item{tfidfMin}{Minimum specificity cut-off used if finding marker genes (see \code{\link{quickMarkers}}).}

\item{...}{Passed to \code{\link{estimateNonExpressingCells}}}
}
\value{
A ggplot2 object containing the plot.
}
\description{
If each cell were made up purely of background reads, the expression fraction would equal that of the soup.  This plot compares this expectation of pure background to the observed expression fraction in each cell, for each of the groups of genes in \code{nonExpressedGeneList}.  For each group of genes, the distribution of this ratio is plotted across all cells.  A value significantly greater than 1 (0 on log scale) can only be obtained if some of the genes in each group are genuinely expressed by the cell.  That is, the assumption that the cell is pure background does not hold for that gene.
}
\details{
This plot is a useful diagnostic for the assumption that a list of genes is non-expressed in most cell types.  For non-expressed cells, the ratio should cluster around the contamination fraction, while for expressed cells it should be elevated.  The most useful non-expressed gene sets are those for which the genes are either strongly expressed, or not expressed at all.  Such groups of genes will show up in this plot as a bimodal distribution, with one mode containing the cells that do not express these genes around the contamination fraction for this channel and another around a value at some value equal to or greater than 0 (1 on non-log scale) for the expressed cells.

The red line shows the global estimate of the contamination for each group of markers.  This is usually lower than the low mode of the distribution as there will typically be a non-negligible number of cells with 0 observed counts (and hence -infinity log ratio).

If \code{nonExpressedGeneList} is missing, this function will try and find genes that are very specific to different clusters, as these are often the most useful in estimating the contamination fraction.   This is meant only as a heuristic, which can hopefully provide some inspiration as to a class of genes to use to estimation the contamination for your experiment.  Please do **NOT** blindly use the top N genes found in this way to estimate the contamination.  That is, do not feed this list of genes into \code{\link{calculateContaminationFraction}} without any manual consideration or filtering as this *will over-estimate your contamination* (often by a large amount).  For this reason, these gene names are not returned by the function.
}
\examples{
gg = plotMarkerDistribution(scToy,list(CD7='CD7',LTB='LTB'))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjustCounts.R
\name{adjustCounts}
\alias{adjustCounts}
\title{Remove background contamination from count matrix}
\usage{
adjustCounts(
  sc,
  clusters = NULL,
  method = c("subtraction", "soupOnly", "multinomial"),
  roundToInt = FALSE,
  verbose = 1,
  tol = 0.001,
  pCut = 0.01,
  ...
)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{clusters}{A vector of cluster IDs, named by cellIDs.  If NULL clusters auto-loaded from \code{sc}.  If FALSE, no clusters are used.  See details.}

\item{method}{Method to use for correction.  See details.  One of 'multinomial', 'soupOnly', or 'subtraction'}

\item{roundToInt}{Should the resulting matrix be rounded to integers?}

\item{verbose}{Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.}

\item{tol}{Allowed deviation from expected number of soup counts.  Don't change this.}

\item{pCut}{The p-value cut-off used when \code{method='soupOnly'}.}

\item{...}{Passed to expandClusters.}
}
\value{
A modified version of the table of counts, with background contamination removed.
}
\description{
After the level of background contamination has been estimated or specified for a channel, calculate the resulting corrected count matrix with background contamination removed.
}
\details{
This essentially subtracts off the mean expected background counts for each gene, then redistributes any "unused" counts.  A count is unused if its subtraction has no effect.  For example, subtracting a count from a gene that has zero counts to begin with.

As expression data is highly sparse at the single cell level, it is highly recommended that clustering information be provided to allow the subtraction method to share information between cells.  Without grouping cells into clusters, it is difficult (and usually impossible) to tell the difference between a count of 1 due to background contamination and a count of 1 due to endogenous expression.  This ambiguity is removed at the cluster level where counts can be aggregated across cells.  This information can then be propagated back to the individual cell level to provide a more accurate removal of contaminating counts.

To provide clustering information, either set clustering on the SoupChannel object with \code{\link{setClusters}} or explicitly passing the \code{clusters} parameter.  

If \code{roundToInt=TRUE}, this function will round the result to integers.  That is, it will take the floor of the connected value and then round back up with probability equal to the fractional part of the number.

The \code{method} parameter controls how the removal of counts in performed.  This should almost always be left at the default ('subtraction'), which iteratively subtracts counts from all genes as described above.  The 'soupOnly' method will use a p-value based estimation procedure to identify those genes that can be confidently identified as having endogenous expression and removes everything else (described in greater detail below).  Because this method either removes all or none of the expression for a gene in a cell, the correction procedure is much faster.  Finally, the 'multinomial' method explicitly maximises the multinomial likelihood for each cell.  This method gives essentially identical results as 'subtraction' and is considerably slower.

In greater detail, the 'soupOnly' method is done by sorting genes within each cell by their p-value under the null of the expected soup fraction using a Poisson model.  So that genes that definitely do have a endogenous contribution are at the end of the list with p=0.  Those genes for which there is poor evidence of endogenous cell expression are removed, until we have removed approximately nUMIs*rho molecules.  The cut-off to prevent removal of genes above nUMIs*rho in each cell is achieved by calculating a separate p-value for the total number of counts removed to exceed nUMIs*rho, again using a Poisson model.  The two p-values are combined using Fisher's method and the cut-off is applied to the resulting combined p-value calculated using a chi-squared distribution with 4 degrees of freedom.
}
\examples{
out = adjustCounts(scToy)
#Return integer counts only
out = adjustCounts(scToy,roundToInt=TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{expandClusters}
\alias{expandClusters}
\title{Expands soup counts calculated at the cluster level to the cell level}
\usage{
expandClusters(clustSoupCnts, cellObsCnts, clusters, cellWeights, verbose = 1)
}
\arguments{
\item{clustSoupCnts}{Matrix of genes (rows) by clusters (columns) where counts are number of soup counts for that gene/cluster combination.}

\item{cellObsCnts}{Matrix of genes (rows) by cells (columns) giving the observed counts}

\item{clusters}{Mapping from cells to clusters.}

\item{cellWeights}{Weighting to give to each cell when distributing counts.  This would usually be set to the number of expected soup counts for each cell.}

\item{verbose}{Integer giving level of verbosity.  0 = silence, 1 = Basic information, 2 = Very chatty, 3 = Debug.}
}
\value{
A matrix of genes (rows) by cells (columns) giving the number of soup counts estimated for each cell.  Non-integer values possible.
}
\description{
Given a clustering of cells and soup counts calculated for each of those clusters, determines a most likely allocation of soup counts at the cell level.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setProperties.R
\name{setContaminationFraction}
\alias{setContaminationFraction}
\title{Manually set contamination fraction}
\usage{
setContaminationFraction(sc, contFrac, forceAccept = FALSE)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{contFrac}{The contamination fraction.  Either a constant, in which case the same value is used for all cells, or a named vector, in which case the value is set for each cell.}

\item{forceAccept}{A warning or error is usually returned for extremely high contamination fractions.  Setting this to TRUE will turn these into messages and proceed.}
}
\value{
A modified SoupChannel object for which the contamination (rho) has been set.
}
\description{
Manually specify the contamination fraction.
}
\examples{
sc = load10X(system.file('extdata','toyData',package='SoupX'))
sc = setContaminationFraction(sc,0.1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quickMarkers.R
\name{quickMarkers}
\alias{quickMarkers}
\title{Gets top N markers for each cluster}
\usage{
quickMarkers(toc, clusters, N = 10, FDR = 0.01, expressCut = 0.9)
}
\arguments{
\item{toc}{Table of counts.  Must be a sparse matrix.}

\item{clusters}{Vector of length \code{ncol(toc)} giving cluster membership.}

\item{N}{Number of marker genes to return per cluster.}

\item{FDR}{False discover rate to use.}

\item{expressCut}{Value above which a gene is considered expressed.}
}
\value{
data.frame with top N markers (or all that pass the hypergeometric test) and their statistics for each cluster.
}
\description{
Uses tf-idf ordering to get the top N markers of each cluster.  For each cluster, either the top N or all genes passing the hypergeometric test with the FDR specified, whichever list is smallest.
}
\details{
Term Frequency - Inverse Document Frequency is used in natural language processing to identify terms specific to documents.  This function uses the same idea to order genes within a group by how predictive of that group they are.  The main advantage of this is that it is extremely fast and gives reasonable results.

To do this, gene expression is binarised in each cell so each cell is either considered to express or not each gene.  That is, we replace the counts with \code{toc > zeroCut}.  The frequency with which a gene is expressed within the target group is compared to the global frequency to calculate the tf-idf score.  We also calculate a multiple hypothesis corrected p-value based on a hypergeometric test, but this is extremely permissive.
}
\examples{
#Calculate markers of clusters in toy data
mrks = quickMarkers(scToy$toc,scToy$metaData$clusters)
\dontrun{
#Calculate markers from Seurat (v3) object
mrks = quickMarkers(srat@assays$RNA@count,srat@active.ident)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setProperties.R
\name{setSoupProfile}
\alias{setSoupProfile}
\title{Set soup profile}
\usage{
setSoupProfile(sc, soupProfile)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{soupProfile}{A data.frame with columns \code{est} containing the fraction of soup for each gene, \code{counts} containing the total counts for each gene and with row names corresponding to the row names of \code{sc$toc}.}
}
\value{
An updated SoupChannel object with the soup profile set.
}
\description{
Manually sets or updates the soup profile for a SoupChannel object.
}
\examples{
#Suppose only table of counts is available
toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
#Suppress calculating soup profile automatically
sc = SoupChannel(toc,toc,calcSoupProfile=FALSE)
#And add manually
rowSums = Matrix::rowSums
soupProf = data.frame(row.names = rownames(toc),est=rowSums(toc)/sum(toc),counts=rowSums(toc))
sc = setSoupProfile(sc,soupProf)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load10X.R
\name{load10X}
\alias{load10X}
\title{Load a collection of 10X data-sets}
\usage{
load10X(
  dataDir,
  cellIDs = NULL,
  channelName = NULL,
  readArgs = list(),
  includeFeatures = c("Gene Expression"),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{dataDir}{Top level cellranger output directory (the directory that contains the \code{raw_gene_bc_matrices} folder).}

\item{cellIDs}{Barcodes of droplets that contain cells.  If NULL, use the default cellranger set.}

\item{channelName}{The name of the channel to store.  If NULL set to either \code{names(dataDir)} or \code{dataDir} is no name is set.}

\item{readArgs}{A list of extra parameters passed to \code{Seurat::Read10X}.}

\item{includeFeatures}{If multiple feature types are present, keep only the types mentioned here and collapse to a single matrix.}

\item{verbose}{Be verbose?}

\item{...}{Extra parameters passed to \code{SoupChannel} construction function.}
}
\value{
A SoupChannel object containing the count tables for the 10X dataset.
}
\description{
Loads unfiltered 10X data from each data-set and identifies which droplets are cells using the cellranger defaults.
}
\examples{
sc = load10X(system.file('extdata','toyData',package='SoupX'))
}
\seealso{
SoupChannel estimateSoup
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{alloc}
\alias{alloc}
\title{Allocate values to "buckets" subject to weights and constraints}
\usage{
alloc(tgt, bucketLims, ws = rep(1/length(bucketLims), length(bucketLims)))
}
\arguments{
\item{tgt}{Value to distribute between buckets.}

\item{bucketLims}{The maximum value that each bucket can take.  Must be a vector of positive values.}

\item{ws}{Weights to be used for each bucket.  Default value makes all buckets equally likely.}
}
\value{
A vector of the same length as \code{bucketLims} containing values distributed into buckets.
}
\description{
Allocates \code{tgt} of something to \code{length(bucketLims)} different "buckets" subject to the constraint that each bucket has a maximum value of \code{bucketLims} that cannot be exceeded.  By default counts are distributed equally between buckets, but weights can be provided using \code{ws} to have the redistribution prefer certain buckets over others.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/package.R
\docType{package}
\name{SoupX}
\alias{SoupX}
\title{SoupX: Profile, quantify and remove ambient RNA expression from droplet based RNA-seq}
\description{
This package implements the method described in REF.  First a few notes about nomenclature:
soup - Used a shorthand to refer to the ambient RNA which is contained in the input solution to droplet based RNA-seq experiments and ends up being sequenced along with the cell endogenous RNAs that the experiment is aiming to quantify.
channel - This refers to a single run input into a droplet based sequencing platform.  For Chromium 10X 3' sequencing there are currently 8 "channels" per run of the instrument.  Because the profile of the soup depends on the input solution, this is the minimal unit on which the soup should be estimated and subtracted.
}
\details{
The essential step in performing background correction is deciding which genes are not expressed in a reasonable fraction of cells.  This is because SoupX estimates the contamination fraction by comparing the expression of these non-expressed genes in droplets containing cells to the soup defined from empty droplets.  For solid tissue, the set of Haemoglobin genes usually works well.  The key properties a gene should have are:
- it should be easy to identify when it is truly expressed (i.e., when it's expressed, it should be highly expressed) 
- it should be highly specific to a certain cell type or group of cell types so that when the expression level is low, you can be confident that the expression is coming from the soup and not a very low level of expression from the cell

Spike-in RNAs are the best case scenario.  In the case where you do not have spike-ins and haemoglobin genes are not viable estimators, the user should begin by using the \link{plotMarkerDistribution} function to plot those genes with bi-modal distributions that have a pattern of expression across cells that is consistent with high cell-type specificity.  The user should then select a set of genes that can be used for estimation from this list.  One or two high quality genes is usually sufficient to obtain a good estimate for the average contamination level of a channel.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{PBMC_sc}
\alias{PBMC_sc}
\title{SoupChannel from PBMC data}
\format{
\code{PBMC_sc} is a \code{SoupChannel} object with 33,694 genes and 2,170 cells.
}
\source{
\url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k}
}
\usage{
data(PBMC_sc)
}
\description{
\code{\link{SoupChannel}} created from 10X demonstration PBMC 4k data.  The cells have been sub-sampled by a factor of 2 to reduce file size of package.
}
\details{
\code{PBMC_sc} was created by running the following commands.
\itemize{
  \item \code{set.seed(1137)}
  \item \code{tmpDir = tempdir(check=TRUE)}
  \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'tod.tar.gz'))}
  \item \code{download.file('http://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz',destfile=file.path(tmpDir,'toc.tar.gz'))}
  \item \code{untar(file.path(tmpDir,'tod.tar.gz'),exdir=tmpDir)}
  \item \code{untar(file.path(tmpDir,'toc.tar.gz'),exdir=tmpDir)}
  \item \code{library(SoupX)}
  \item \code{PBMC_sc = load10X(tmpDir,calcSoupProfile=FALSE)}
  \item \code{PBMC_sc = SoupChannel(PBMC_sc$tod,PBMC_sc$toc[,sample(ncol(PBMC_sc$toc),round(ncol(PBMC_sc$toc)*0.5))])}
}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classFunctions.R
\name{SoupChannel}
\alias{SoupChannel}
\title{Construct a SoupChannel object}
\usage{
SoupChannel(tod, toc, metaData = NULL, calcSoupProfile = TRUE, ...)
}
\arguments{
\item{tod}{Table of droplets.  A matrix with columns being each droplet and rows each gene.}

\item{toc}{Table of counts.  Just those columns of \code{tod} that contain cells.}

\item{metaData}{Meta data pertaining to the cells.  Optional.  Must be a data-frame with rownames equal to column names of \code{toc}.}

\item{calcSoupProfile}{By default, the soup profile is calculated using \code{\link{estimateSoup}} with default values.  If you want to do something other than the defaults, set this to \code{FALSE} and call \code{\link{estimateSoup}} manually.}

\item{...}{Any other named parameters to store.}
}
\value{
A SoupChannel object.
}
\description{
Creates a SoupChannel object that contains everything related to the soup estimation of a single channel.
}
\examples{
#Load droplet and count tables
tod = Seurat::Read10X(system.file('extdata','toyData','raw_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
#Default calculates soup profile
sc = SoupChannel(tod,toc)
names(sc)
#This can be suppressed
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
names(sc)
}
\seealso{
SoupChannelList estimateSoup setSoupProfile setClusters
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateContaminationFraction.R
\name{calculateContaminationFraction}
\alias{calculateContaminationFraction}
\title{Calculate the contamination fraction}
\usage{
calculateContaminationFraction(
  sc,
  nonExpressedGeneList,
  useToEst,
  verbose = TRUE,
  forceAccept = FALSE
)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{nonExpressedGeneList}{A list containing sets of genes which can be assumed to be non-expressed in a subset of cells (see details).}

\item{useToEst}{A boolean matrix of dimensions ncol(toc) x length(nonExpressedGeneList) indicating which gene-sets should not be assumed to be non-expressed in each cell.  Row names must correspond to the names of \code{nonExpressedGeneList}.  Usually produced by \code{\link{estimateNonExpressingCells}}.}

\item{verbose}{Print best estimate.}

\item{forceAccept}{Passed to \code{\link{setContaminationFraction}}.}
}
\value{
A modified version of \code{sc} with estimates of the contamination (rho) added to the metaData table.
}
\description{
This function computes the contamination fraction using two user-provided bits of information.  Firstly, a list of sets of genes that can be biologically assumed to be absent in at least some cells in your data set.  For example, these might be haemoglobin genes or immunoglobulin genes, which should not be expressed outside of erythroyctes and antibody producing cells respectively.
}
\details{
Secondly, this function needs to know which cells definitely do not express the gene sets described above.  Continuing with the haemoglobin example, which are the erythrocytes that are producing haemoglobin mRNAs and which are non-erythrocytes that we can safely assume produce no such genes.  The assumption made is any expression from a gene set in cell marked as a "non-expressor" for that gene set, must be derived from the soup.  Therefore, the level of contamination present can be estimated from the amount of expression of these genes seen in these cells.

Most often, the genesets are user supplied based on your knowledge of the experiment and the cells in which they are genuinely expressed is estimated using \code{\link{estimateNonExpressingCells}}.  However, they can also be supplied directly if other information is available.

Usually, there is very little variation in the contamination fraction within a channel and very little power to detect the contamination accurately at a single cell level.  As such, the default mode of operation simply estimates one value of the contamination fraction that is applied to all cells in a channel.

The global model fits a simple Poisson glm to the aggregated count data across all cells.

Finally, note that if you are not able to find a reliable set of genes to use for contamination estimation, or you do not trust the values produced, the contamination fraction can be manually set by the user using \code{\link{setContaminationFraction}}.
}
\examples{
#Common gene list in real world data
geneList = list(HB=c('HBB','HBA2'))
#Gene list appropriate to toy data
geneList = list(CD7 = 'CD7')
ute = estimateNonExpressingCells(scToy,geneList)
sc = calculateContaminationFraction(scToy,geneList,ute)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{initProgBar}
\alias{initProgBar}
\title{Create Seurat style progress bar}
\usage{
initProgBar(min, max, ...)
}
\arguments{
\item{min}{Minimum value of parameter.}

\item{max}{Maximum value of parameter.}

\item{...}{Passed to \code{\link{txtProgressBar}}}
}
\value{
A txtProgressBar object to use updating progress.
}
\description{
Creates progress bar that won't ruin log files and shows progress towards 100%.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classFunctions.R
\name{print.SoupChannel}
\alias{print.SoupChannel}
\title{Print method for SoupChannel}
\usage{
\method{print}{SoupChannel}(x, ...)
}
\arguments{
\item{x}{A SoupChannel object.}

\item{...}{Currently unused.}
}
\value{
Nothing.  Prints message to console.
}
\description{
Prints a summary of a SoupChannel object.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setProperties.R
\name{setClusters}
\alias{setClusters}
\title{Sets clustering for SoupChannel}
\usage{
setClusters(sc, clusters)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{clusters}{A named vector, where entries are the cluster IDs and names are cellIDs.  If no names are provided, the order is assumed to match the order in \code{sc$metaData}.}
}
\value{
An updated SoupChannel object with clustering information stored.
}
\description{
Adds or updates clustering information to meta-data table in SoupChannel object.
}
\examples{
sc = load10X(system.file('extdata','toyData',package='SoupX'))
mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
sc = setClusters(sc,mDat$res.1)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/autoEstCont.R
\name{autoEstCont}
\alias{autoEstCont}
\title{Automatically calculate the contamination fraction}
\usage{
autoEstCont(
  sc,
  topMarkers = NULL,
  tfidfMin = 1,
  soupQuantile = 0.9,
  maxMarkers = 100,
  contaminationRange = c(0.01, 0.8),
  rhoMaxFDR = 0.2,
  priorRho = 0.05,
  priorRhoStdDev = 0.1,
  doPlot = TRUE,
  forceAccept = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{sc}{The SoupChannel object.}

\item{topMarkers}{A data.frame giving marker genes.  Must be sorted by decreasing specificity of marker and include a column 'gene' that contains the gene name.  If set to NULL, markers are estimated using \code{\link{quickMarkers}}.}

\item{tfidfMin}{Minimum value of tfidf to accept for a marker gene.}

\item{soupQuantile}{Only use genes that are at or above this expression quantile in the soup.  This prevents inaccurate estimates due to using genes with poorly constrained contribution to the background.}

\item{maxMarkers}{If we have heaps of good markers, keep only the best \code{maxMarkers} of them.}

\item{contaminationRange}{Vector of length 2 that constrains the contamination fraction to lie within this range.  Must be between 0 and 1.  The high end of this range is passed to \code{\link{estimateNonExpressingCells}} as \code{maximumContamination}.}

\item{rhoMaxFDR}{False discovery rate passed to \code{\link{estimateNonExpressingCells}}, to test if rho is less than \code{maximumContamination}.}

\item{priorRho}{Mode of gamma distribution prior on contamination fraction.}

\item{priorRhoStdDev}{Standard deviation of gamma distribution prior on contamination fraction.}

\item{doPlot}{Create a plot showing the density of estimates?}

\item{forceAccept}{Passed to \code{\link{setContaminationFraction}}.  Should we allow very high contamination fractions to be used.}

\item{verbose}{Be verbose?}
}
\value{
A modified SoupChannel object where the global contamination rate has been set.  Information about the estimation is also stored in the slot \code{fit}
}
\description{
The idea of this method is that genes that are highly expressed in the soup and are marker genes for some population can be used to estimate the background contamination.  Marker genes are identified using the tfidf method (see \code{\link{quickMarkers}}).  The contamination fraction is then calculated at the cluster level for each of these genes and clusters are then aggressively pruned to remove those that give implausible estimates.
}
\details{
This set of marker genes is filtered to include only those with tf-idf value greater than \code{tfidfMin}.  A higher tf-idf value implies a more specific marker.  Specifically a cut-off t implies that a marker gene has the property that geneFreqGlobal < exp(-t/geneFreqInClust).  See \code{\link{quickMarkers}}.  It may be necessary to decrease this value for data sets with few good markers.

This set of marker genes is filtered down to include only the genes that are highly expressed in the soup, controlled by the \code{soupQuantile} parameter.  Genes highly expressed in the soup provide a more precise estimate of the contamination fraction.

The pruning of implausible clusters is based on a call to \code{\link{estimateNonExpressingCells}}.  The parameters \code{maximumContamination=max(contaminationRange)} and \code{rhoMaxFDR} are passed to this function.  The defaults set here are calibrated to aggressively prune anything that has even the weakest of evidence that it is genuinely expressed. 

For each cluster/gene pair the posterior distribution of the contamination fraction is calculated (based on gamma prior, controlled by \code{priorRho} and \code{priorRhoStdDev}).  These posterior distributions are aggregated to produce a final estimate of the contamination fraction. The logic behind this is that estimates from clusters that truly estimate the contamination fraction will cluster around the true value, while erroneous estimates will be spread out across the range (0,1) without a 'preferred value'.  The most probable value of the contamination fraction is then taken as the final global contamination fraction.
}
\examples{
#Use less specific markers
scToy = autoEstCont(scToy,tfidfMin=0.8)
#Allow large contamination fractions to be allocated
scToy = autoEstCont(scToy,forceAccept=TRUE)
#Be quiet
scToy = autoEstCont(scToy,verbose=FALSE,doPlot=FALSE)
}
\seealso{
quickMarkers
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateSoup.R
\name{estimateSoup}
\alias{estimateSoup}
\title{Get expression profile of soup}
\usage{
estimateSoup(sc, soupRange = c(0, 100), keepDroplets = FALSE)
}
\arguments{
\item{sc}{A \code{SoupChannel} object.}

\item{soupRange}{Droplets with total UMI count in this range (excluding endpoints) are used to estimate soup.}

\item{keepDroplets}{Storing the full table of counts for all droplets uses a lot of space and is really only used to estimate the soup profile.  Therefore, it is dropped after the soup profile has been estimated unless this is set to \code{TRUE}.}
}
\value{
A modified version of \code{sc} with an extra \code{soupProfile} entry containing a data.frame with the soup profile and confidence limits for all genes.
}
\description{
This is usually called by \code{\link{SoupChannel}}, rather than directly by the user.  Uses the empty droplets in the range provided to calculate the expression profile of the soup under the assumption that these droplets only contain background.
}
\examples{
#Load droplet and count tables
tod = Seurat::Read10X(system.file('extdata','toyData','raw_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
toc = Seurat::Read10X(system.file('extdata','toyData','filtered_gene_bc_matrices','GRCh38',
                                  package='SoupX'))
#Suppress calculation of soup profile automatically on load
sc = SoupChannel(tod,toc,calcSoupProfile=FALSE)
#Retain table of droplets
sc = estimateSoup(sc,keepDroplets=TRUE)
#Or use non-default values
sc = estimateSoup(sc,soupRange=c(60,100))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotMarkerMap}
\alias{plotMarkerMap}
\title{Plot ratio of observed to expected counts on reduced dimension map}
\usage{
plotMarkerMap(
  sc,
  geneSet,
  DR,
  ratLims = c(-2, 2),
  FDR = 0.05,
  useToEst = NULL,
  pointSize = 2,
  pointShape = 21,
  pointStroke = 0.5,
  naPointSize = 0.25
)
}
\arguments{
\item{sc}{SoupChannel object.}

\item{geneSet}{A vector with the names of the genes to aggregate and plot evidence for.}

\item{DR}{A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.  Try and fetch automatically if missing.}

\item{ratLims}{Truncate log ratios at these values.}

\item{FDR}{False Discovery Rate for statistical test of enrichment over background.}

\item{useToEst}{A vector (usually obtained from \code{\link{estimateNonExpressingCells}}), that will be used to mark cells instead of the usual Poisson test.}

\item{pointSize}{Size of points}

\item{pointShape}{Shape of points}

\item{pointStroke}{Stroke size for points}

\item{naPointSize}{Point size for NAs.}
}
\value{
A ggplot2 containing the plot.
}
\description{
Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how likely a set of genes are to be soup derived on that map.  That is, given a set of genes, this function calculates how many counts would be expected if that droplet were nothing but soup and compares that to the observed count.  This is done via a log2 ratio of the two values.  A Poisson test is performed and points that have a statistically significant enrichment over the background (at 5% FDR) are marked.
}
\examples{
gg = plotMarkerMap(scToy,'CD7')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{plotChangeMap}
\alias{plotChangeMap}
\title{Plot maps comparing corrected/raw expression}
\usage{
plotChangeMap(
  sc,
  cleanedMatrix,
  geneSet,
  DR,
  dataType = c("soupFrac", "binary", "counts"),
  logData = FALSE,
  pointSize = 0.5
)
}
\arguments{
\item{sc}{SoupChannel object.}

\item{cleanedMatrix}{A cleaned matrix to compare against the raw one.  Usually the output of \code{\link{adjustCounts}}.}

\item{geneSet}{A vector with the names of the genes to aggregate and plot evidence for.}

\item{DR}{A data.frame, with rows named by unique cell IDs (i.e., <ChannelName>_<Barcode>) the first two columns of which give the coordinates of each cell in some reduced dimension representation of the data.}

\item{dataType}{How should data be represented.  Binary sets each cell to expressed or not, counts converts everything to counts, soupFrac plots the fraction of the observed counts that are identified as contamination (i.e., (old-new)/old) for each cell and is the default.}

\item{logData}{Should we log the thing we plot?}

\item{pointSize}{Size of points}
}
\value{
A ggplot2 containing the plot.
}
\description{
Given some reduced dimensional representation of the data (such as UMAP or tSNE) that has been calculated however you would like, this provides a way to visualise how the expression of a geneSet changes after soup correction.
}
\examples{
out = adjustCounts(scToy)
gg = plotChangeMap(scToy,out,'S100A9')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setProperties.R
\name{setDR}
\alias{setDR}
\title{Manually set dimension reduction for a channel}
\usage{
setDR(sc, DR, reductName = NULL)
}
\arguments{
\item{sc}{A SoupChannel object.}

\item{DR}{The dimension reduction coordinates (e.g., tSNE).  This must be a data.frame, with two columns giving the two dimension reduction coordinates.  The data.frame must either have row names matching the row names of sc$metaData, or be ordered in the same order as sc$metaData.}

\item{reductName}{What to name the reduction (defaults to column names provided).}
}
\value{
A modified SoupChannel object for which the dimension reduction has been set.
}
\description{
Manually specify the dimension reduction
}
\examples{
sc = load10X(system.file('extdata','toyData',package='SoupX'))
mDat = read.table(system.file('extdata','toyData','metaData.tsv',package='SoupX'),sep='\t')
sc = setDR(sc,mDat[,c('tSNE_1','tSNE_2')])
}
