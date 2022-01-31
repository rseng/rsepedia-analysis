# DENDRO & DENDROplan

 **DENDRO**, stands for **D**na based **E**volutio**N**ary tree pre**D**iction by sc**R**na-seq techn**O**logy, is an **R** package, which takes scRNA-seq data for a tumor (or related somatic tissues) and accurately reconstructs its phylogeny, assigning each single cell from the single cell RNA sequencing (scRNA-seq) data to a subclone (Figure 1). Currently there is no phylogenetic reconstruction framework specifically designs for scRNA-seq dataset due to biological dropout (i.e. burstness), sequencing error, and technical dropout. DENDRO perfectly tackles this problem with a Bayesian framework (Beta-Binomial), and achieves high clustering accuracy .

In addition, before conducting a single cell RNA-seq experiment on a tumor sample, it is important to project how subclone detection power depends on the number of cells sequenced and the coverage per cell. To facilitate experiment design, we developed a tool, **DENDROplan** (Figure 2), that  predicts the expected clustering accuracy by DENDRO given sequencing parameters. As a result, researchers can design experiment parameters, such as sequencing depth and number of cells, based on DENDROplan's prediction.


## Manuscript

[DENDRO: genetic heterogeneity profiling and subclone detection by single-cell RNA sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1922-x)


## Questions & Problems

If you have any questions or problems when using DENDRO or DENDROplan, please feel free to open a new issue [here](https://github.com/zhouzilu/DENDRO/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Installation

Install to R/RStudio
Install all packages in the latest version of [R](https://www.r-project.org/).
```r
devtools::install_github("zhouzilu/DENDRO")
```
If you observe error with Biobase try the following and then try reinstall.
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8")
```


## Pipeline overview

This DENDRO package includes two analysis tools: (1) **DENDRO**, a phylogenetic tree construction with real dataset such as tumor and hematopoesis scRNA-seq, and (2) **DENDROplan**, which help design experiment by predicting the accuracy of DENDRO cluster given inferred clonal tree structure, cell number and sequencing depth. Overall pipelines are shown below.

### DENDRO pipeline

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-01.jpg' width='1000' height='600'>
  </p>

  **Figure 1.** A flowchart outlining the procedures of DENDRO. DENDRO starts from scRNA-seq raw data. We recommend STAR 2-pass method for mapping because it is more robust with splicing junction. SNA detection was applied to mapped BAM files. Both counts of total allele reads and counts of alternative allele reads for each cell $c$ at mutation position $g$ are collected. In the next step, a cell-to-cell genetic divergence matrix is calculated using a genetic divergence evaluation function. DENDRO further clusters the cells and pools cells from same cluster together and re-estimate SNA profiles. Based on the re-estimated SNA profiles, DENDRO generates a parsimony tree which shows the evolution relationship between subclones.

### Running DENDRO

  **DENDRO R notebook** with step-by-step demonstration and rich display is available [***here***](http://raw.githack.com/zhouzilu/DENDRO/master/vignette/DENDRO_vignette.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/DENDRO/blob/master/vignette/DENDRO_vignette.Rmd).


### DENDROplan pipeline

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-02.jpg' width='1000' height='600'>
  </p>

  **Figure 2.** The overall DENDROplan pipeline. The analysis starts with a designed tree with an interested clade (purple clade in the example). Based on the tree model, number of cells, sequencing depth and sequencing error rate, we simulate single cell mutation profile. scRNA data was sampled from a reference scRNA-seq dataset given expression level in bulk. A phylogeny computed by DENDRO is further compared with underlining truth, which measured by three statistics - adjust Rand index (global accuracy statistics), capture rate (subclone specific statistic) and purity (subclone specific statistic). 

### Running DENDROplan

  **DENDROplan R notebook** with step-by-step demonstration and rich display is available [***here***](http://raw.githack.com/zhouzilu/DENDRO/master/vignette/DENDROplan_vignette.html). Corresponding **Rmd script** is available [***here***](https://github.com/zhouzilu/DENDRO/blob/master/vignette/DENDROplan_vignette.Rmd).


## Citation

Zhou, Z., Xu, B., Minn, A. et al. DENDRO: genetic heterogeneity profiling and subclone detection by single-cell RNA sequencing. Genome Biol 21, 10 (2020). https://doi.org/10.1186/s13059-019-1922-x

<br>
  Genetic Heterogeneity Profiling by Single Cell RNA Sequencing ([GitHub](https://github.com/zhouzilu/DENDRO))

## Developers & Maintainers

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, University of Pennsylvania

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) 
  <br>
  Department of Statistics, University of Pennsylvania
# DENDRO mutation detection with GATK

 **DENDRO** uses the following approach in the original paper to generate the initial mutation profiles across genes and cells. We write [an example script](https://github.com/zhouzilu/DENDRO/blob/master/script/mutation_detection_mapreduce.sh) modified from one real data analysis. Just want to clarify that I am not an expert with GATK tool. However, given significant amount of tests and failures, the following pipeline works for DENDRO. See [HaplotypeCaller document](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) for more information.


## Questions, Suggestions & Problems

If you have any questions or problems when using DENDRO or DENDROplan, please feel free to open a new issue [here](https://github.com/zhouzilu/DENDRO/issues). You can also email the maintainers of the corresponding packages -- the contact information is shown under Developers & Maintainers.


## Pipeline overview

Initial SNA detection is indeed one of the most important steps in DENDRO pipeline. As stated in the paper, we want to maximize the sensitivity of our detection at this moment. There will be  a lot of false positives calls, but they will be cleaned up in the clustering and pooling steps. Figure 4 shows the detection pipeline. User should check the [example script](https://github.com/zhouzilu/DENDRO/blob/master/script/mutation_detection_mapreduce.sh) together with this figure side-by-side.

<p align="center">
  <img src='https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-04.jpg' width='1000' height='600'>
  </p>

  **Figure 5.** A flowchart outlining the procedures of mutation detection from fastq to final GVCF using GATK [ERC GVCF approach](https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode). The number on the most right corresponding to the steps in the shell script.
  
  
Due to great number of cells, traditional way of calling variants sample by sample is extremely slow. Luckily, GATK has ERC GVCF mode, which utilize a map-reduce like approach. Please check the details [here](https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode). Our script is built upon ERC GVCF mode togehter with the [variant detection best approach on RNAseq data](https://gatkforums.broadinstitute.org/gatk/discussion/3892/the-gatk-best-practices-for-variant-calling-on-rnaseq-in-full-detail).

Example script is attached [here](https://github.com/zhouzilu/DENDRO/blob/master/script/mutation_detection_mapreduce.sh). We also provided an [R script](https://github.com/zhouzilu/DENDRO/blob/master/script/vcf_to_DENDROinput.R) to convert vcf files to DENDRO input.

## Citation

Zhou, Z., Xu, B., Minn, A. et al. DENDRO: genetic heterogeneity profiling and subclone detection by single-cell RNA sequencing. Genome Biol 21, 10 (2020). https://doi.org/10.1186/s13059-019-1922-x

## Developers & Maintainers

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, University of Pennsylvania

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/)
  <br>
  Department of Statistics, University of Pennsylvania
---
title: "DENDROplan Vignette"
author: "Zilu Zhou"
date: "10/25/2019"#"`r Sys.Date()`"
abstract: >
 Before conducting a single cell RNA-seq experiment on a tumor sample, it is important to project how subclone detection power depends on the number of cells sequenced and the coverage per cell. To facilitate experiment design, we developed a tool, **DENDROplan** (Figure 2), that  predicts the expected clustering accuracy by DENDRO given sequencing parameters.  Given a tree structure and a target accuracy, **DENDROplan** computes the necessary read depth and number of cells needed based on the spike-in procedure described above. For more detail, please check our [biorixv preprint](www.rstudio.com)
output:
  rmarkdown::html_document:
    theme: united
    highlight: tango
    toc: true
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: DENDRO.bibtex
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# 1. Installation
Install all packages in the latest version of [R](https://www.r-project.org/).
```{r, eval=FALSE}
devtools::install_github("zhouzilu/DENDRO")
```

The reference dataset from Deng et al. [@Deng2014] is pre-stored in the DENDRO package.

# 2. Questions & issues
If you have any questions or problems when using DENDROplan, please feel free to open a new issue [here](https://github.com/zhouzilu/DENDRO/issues). You can also email the maintainers of the package -- the contact information is below.

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, UPenn

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, UPenn

# 3. DENDROplan pipeline
## 3.1 Overall pipeline

Figure 2 illustrate the overall pipeline. The analysis starts with a designed tree with an interested clade (purple clade in the example). Based on the tree model, number of cells, sequencing depth and sequencing error rate, we simulate single cell mutation profile. scRNA data was sampled from a reference scRNA-seq dataset given expression level in bulk. A phylogeny computed by DENDRO is further compared with underlining truth, which measured by three statistics - adjusted Rand index (global accuracy statistics), capture rate (subclone specific statistic) and purity (subclone specific statistic). 
```{r, out.width = "1000px", fig.align = "center", echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-02.jpg")
```
  **Figure 2**. A flowchart outlining the procedures for DENDROplan. 

## 3.2 Load reference
Our reference dataset is from Deng et al. [@Deng2014] with great sequencing quality and depth. First, let's load `DENDROplan_ref` and check the structure of the reference variable `ref`

```{r, message=FALSE, warning=FALSE}
library(DENDRO)
data("DENDROplan_ref")
str(ref)
```

There are total 22958 potential mutation loci across 133 cells. As a result, when we simulate trees, the maximum number of cells is 133 and the maximum number of mutations is 22958.

## 3.3  Tree simulation
### 3.3.1 Tree type and parameters

The trees that DENDROplan is able to simulate show at Figure 3. The corresponding `k` and `subtype` value are also included. 

```{r, out.width = "1000px", fig.align = "center", echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-03.jpg")
```

  **Figure 3**. Example tree structure DENDRO can simulate with corresponding parameters.

Numbers on branches illustrate index of mutation probability variable `lprob`, and numbers at node illustrate index of cell proportion variable `kprob`. The interested clade is labeled by a *star*. Capture rate and purity are measured on that specific *star* cluster. `k` is the parameter indicating number of subclones. When `k=4`, the subclones can form two types of trees, differetiated by the `subtype` parameter.

### 3.3.2 Accuracy measurements

We measure DENDROplan results by four statistics, explained here:

* Adjusted Rand index: Adjusted Rand index is a measure of the similarity between two data clusterings after adjusted for the chance grouping of elements. For details, see [here](https://en.wikipedia.org/wiki/Rand_index)

* Capture rate: Capture rate is a measure of “false negative rate” of a specific *star* clade. Out of all the cells from the specific clade, how many of them is detected by the algorithm.

* Purity: Purity is a measure of “false positive rate” of a specific *star* clade. Out of all the cells in the “specific cluster” you detected, how many are actually from the true specific clade.

* Observation probability: This is only one single value, which measures the probability of observe all clades in our multipe simulation round.

## 3.4 Practice
### 3.4.1 Exampe I

Assume we only have 100 mutation sites and our tree looks like Figure 3c with mutations and cells uniformly distributed

```{r, message=FALSE, warning=FALSE, fig.width=8}
res=DENDRO.simulation(n=100,ref=ref,k=4,subtype=1)
res
```

Result shows the example tree structure as well as statistics with 95% Confidence Interval (CI).

### 3.4.2 Example II

In the 2nd example, we want to specify the similar tree but with customized cell proportion and mutation loads. We want the four clusters have cell proportion 0.2, 0.2, 0.2 and 0.4, i.e. we want one major subclone. In our code we need to specify this by parameter `kprob`. Also, we want our mutation load unequally distributed with branch 4 (cluster 4 specific) having more mutations `lprob=c(0.15,0.15,0.15,0.25,0.15,0.15)`.

```{r, message=FALSE, warning=FALSE, fig.width=8}
res=DENDRO.simulation(kprob=c(0.2,0.2,0.2,0.4),lprob=c(0.15,0.15,0.15,0.25,0.15,0.15),n=100,ref=ref,k=4,subtype=1)
res
```

Example tree shows consistency with our input.

### 3.4.3 Example III

One important reason for user to use this tool is that it can estimate the cluster accuracy given different sequencing depth. In the original Deng et al. paper [@Deng2014], they collect around 10,000,000 50bp-reads per cell in mice, resulting 46 mapped reads for each mutation site, which is super high depth (Figure 4). 

```{r, out.width = "1000px", fig.align = "center", echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/deng_stat-04.jpg")
```

  **Figure 4**. scRNA-seq library statistics for Deng et al. [@Deng2014]. Modified from [original paper](http://science.sciencemag.org/content/sci/suppl/2014/01/09/343.6167.193.DC1/Deng-SM.pdf).

In real life, usually 1,000,000 reads per cell is pretty good depth. Proportionally, there are around 4.5 mapped reads per mutation site. Let's change the read depth parameter `RD` and see how well DENDRO performs.

```{r}
res=DENDRO.simulation(RD=4.5,n=100,ref=ref,k=4,subtype=1)
res
```

Feel free to play with different parameters. The detailed function document can be found by typing `??DENDRO.simulation` in R.

# 4. Session info

```{r sessionInfo}
sessionInfo()
```

# 5. References

---
title: "DENDRO Vignette"
author: "Zilu Zhou"
date: "10/25/2019"#"`r Sys.Date()`"
abstract: >
 *D*na based *E*volutio*N*ary tree pre*D*iction by sc*R*na-seq techn*O*logy (DENDRO) is a statistical framework that takes scRNA-seq data for a tumor and accurately reconstructs its phylogeny, assigning each single cell from the scRNA-seq data to a subclone. Our approach allows us to (1) cluster cells based on genomic profiles while accounting for transcriptional bursting, technical dropout and sequencing error, as benchmarked on in silico mixture and a spike-in analysis, (2) make robust mutation profile inference in subclonal resolution, and (3) perform DNA-RNA joint analysis with same set of cells and examine the relationship between transcriptomic variation and mutation profile. For more detail, please check our [biorixv preprint, no link yet](www.rstudio.com)
output: 
  rmarkdown::html_document:
    theme: united
    highlight: tango
    toc: true
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: DENDRO.bibtex
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# 1. Installation
Install all packages in the latest version of [R](https://www.r-project.org/).
```{r, eval=FALSE}
devtools::install_github("zhouzilu/DENDRO")
```

# 2. Questions & issues
If you have any questions or problems when using DENDRO, please feel free to open a new issue [here](https://github.com/zhouzilu/DENDRO/issues). You can also email the maintainers of the package -- the contact information is below.

* [Zilu Zhou](https://statistics.wharton.upenn.edu/profile/zhouzilu/) (zhouzilu at pennmedicine dot upenn dot edu)
  <br>
  Genomics and Computational Biology Graduate Group, UPenn

* [Nancy R. Zhang](https://statistics.wharton.upenn.edu/profile/nzh/) (nzh at wharton dot upenn dot edu)
  <br>
  Department of Statistics, UPenn

# 3. DENDRO analysis pipeline
## 3.1 Overall pipeline 

Figure 1 illustrate the overall pipeline. DENDRO starts from scRNA-seq raw data. We recommend STAR 2-pass method for mapping because it is more robust with splicing junction. SNA detection was applied to mapped BAM files. Both counts of total allele reads and counts of alternative allele reads for each cell $c$ at mutation position $g$ are collected. In the next step, a cell-to-cell genetic divergence matrix is calculated using a genetic divergence evaluation function. DENDRO further clusters the cells and polls cells from same cluster together and re-estimate SNA profiles. Based on the re-estimated SNA profiles, DENDRO generates a parsimony tree which shows the evolution relationship between subclones.

```{r, out.width = "1000px", fig.align = "center", echo=FALSE}
knitr::include_graphics("https://raw.githubusercontent.com/zhouzilu/DENDRO/master/figure/Pkg_FIG-01.jpg")
```
  **Figure 1**. A flowchart outlining the procedures for DENDRO. We separate our analysis pipeline into three stages. The subtask is labeled on the right.

## 3.2 Stage I
### 3.2.1 Initial SNA detection with GATK

Starting with scRNA-seq dataset, we first detect primary mutation with GATK tools. Due to large amount of cells, a map-reduce ERC GVCF approach is necessary. An detailed illustration is attached [here](https://github.com/zhouzilu/DENDRO/blob/master/script/vcf_to_DENDROinput.R). After generated the VCF files, DENDRO utilized information of (1) number of alternative allele read counts $X$, (2) number of total allele read counts $N$, and (3) mutation profile matrix $Z$, where $Z=1$ indicates mutations for each cell and loci. $X, N, Z \in R^{M \times C}$. $M$ is total number of mutation loci and $C$ is total number of cells. We also provided an example R script [here]() for $X$,$N$,$Z$ matrix extraction from VCF.

Here we load our demo dataset (200 mutation locus $\times$ 130 cells) generated using spike-in.
```{r, message=FALSE, warning=FALSE}
library(DENDRO)
data("DENDRO_demo")
str(demo)
```

where `Info` indicates mutation information such as chromosome, allele nucleitide and position, and `label` indicates the true label, which we don't have in a real data analysis.

### 3.2.2 Cell and mutation filtering 

Given $X, N$ and $Z$, DENDRO first apply quality control (QC). In default, DENDRO filter out (1) mutations with less than 5% cells detected (too rare) or more than 95% of cells detected (too common), and (2) cells with total read counts deviated than 5 standard deviation. Check `??FilterCellMutation` for more details. 

```{r, message=FALSE, warning=FALSE}
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label,cut.off.VAF = 0.05, cut.off.sd = 5)
```

The above two plots illustrate two distributions: (left) total number of cells for each mutations and (right) total number of mutations for each cell, with the red line indicating filter thresholds.

```{r, message=FALSE, warning=FALSE}
str(demo_qc)
```


## 3.3 Stage II
### 3.3.1 Genetic divergence matrix calculation

Now DENDRO can calculate the genetic divergence matrix. Check `??DENDRO.dist` for more details.

```{r, message=FALSE, warning=FALSE, fig.width=14}
demo_qc$dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
```

### 3.3.2 Clustering with the genetic divergence matrix

Let's apply hierachical clustering and plot out the clustering result colored by known true label: `demo_qc$clade`. To allow user to customize your tree `DENDRO.cluster` inherent arguments from `plot.phylo`, Check `DENDRO.cluster` for more details.

```{r, message=FALSE, warning=FALSE, fig.width=14}
demo_qc$cluster = DENDRO.cluster(demo_qc$dist,label=demo_qc$label,type='fan')
```

Let's decided the optimal number of clusters using an intra-cluster divergence (icd) measurements.

```{r, message=FALSE, warning=FALSE}
demo_qc$icd = DENDRO.icd(demo_qc$dist,demo_qc$cluster)
demo_qc$optK = 3
demo_qc$DENDRO_label = cutree(demo_qc$cluster,demo_qc$optK)
```

We decide the optimal number of cluster by identifying kink or "elbow point" in the icd plot. In this example, `optK = 3`. It is crucial that if there are multiple "elbow point", the *smallest* one is the most robust.

Let's re-plot our data with DENDRO label

```{r, message=FALSE, warning=FALSE, fig.width=14}
demo_qc$cluster = DENDRO.cluster(demo_qc$dist,label=demo_qc$DENDRO_label,type='fan')
```

### 3.3.3 Re-estimate mutation profile within each cluster and QC

DENDRO further re-esimate the subclone-level mutation profiles by pooling all reads within each reads together with a maximum liklihood approach [@li2012likelihood]. Check `??DENDRO.recalculate` for more details.

```{r, message=FALSE, warning=FALSE}
demo_cluster = DENDRO.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name=c('Cluster3','Cluster2','Cluster1'))
```
`cluster.name` specifies the cluster name given the clustering order (1, 2, 3, ...).

## 3.4 Stage III
### 3.4.1 Evolutionary tree construction

Given the filtered cluster-level mutation profiles, we now can construct an neighbor-joining tree using algorithm implemented in package `phangorn`. See `??DENDRO.tree` for more details.

```{r, message=FALSE, warning=FALSE}
DENDRO.tree(demo_cluster$Z)
```

In this phylogenetic tree, Cluster1 has greater genetic divergence compared with Cluster2 and Cluster3, which is consistent with our data generating process. 

### 3.4.2 Other analysis

User could further perform joint differential expression analysis and differential mutation analysis between different subclone groups. Mutation profile across clones is sored at `demo_cluster$Z`.

Differential expression analysis packages are wide-spread. Two methods that I personally preferred are [Seurat MAST implementation](https://satijalab.org/seurat/get_started.html) [@seurat2018] and [scDD](https://bioconductor.org/packages/release/bioc/html/scDD.html) [@scDD2016].

Gene set enrichment analysis is available at [MSigDB, Broad Institute](http://software.broadinstitute.org/gsea/msigdb/) [@gsea2005].


# 4. Session info

```{r sessionInfo}
sessionInfo()
```

# 5. References

\name{DENDRO.simulation}
\alias{DENDRO.simulation}
\title{
Simulation DENDRO performance
}
\description{
DENDRO.simulation assess the DENDRO performance given an imaginary clonal tree with cell numbers, read depth etc. See argument below.
}
\usage{
DENDRO.simulation(kprob = NULL, lprob = NULL, filt = 0, m = 100, n = 1000, epi = 0.001, RD = NULL, ref, k = NULL, subtype = 1, rpt = 100, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kprob}{
  Distribution of cells in each subclone. The order of probability in kprob see tutorial. Default NULL, indicates uniform distribution of cells.
}
  \item{lprob}{
  Distribution of mutations along each branch in the tree. The order of probability in lprob see tutorial. Default NULL, indicates uniform distribution of mutations.
}
  \item{filt}{
  filt threashold indicates mutations with less than filt mutations are filtered out.
}
  \item{m}{
  m indicates number of cells in your simulation. Default 100
}
  \item{n}{
  n indicates number of genes in your simulation. Default 1000
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.001 according to Illunima.
}
  \item{RD}{
  RD indicates read depth of the overall sequencing. Default uses data from Deng et al. 2014, which is 45X with 10,000,000 reads per cell.
}
  \item{ref}{
  The reference dataset, with X1 and X2 matrices, indicating two individual allele read counts.
}
  \item{k}{
  When kprob is NULL, this will decide the total number of subclones. Combined with subtype, they will determine the clonal tree structure. Default NULL.
}
  \item{subtype}{
  When k=4, there are two different clonal tree structure. This will determine which one is it. Default 1
}
  \item{rpt}{
  How many round of simulation do you want to run to generate the estimation and confidence interval. Default 100.
}
  \item{plot}{
  Whether you want to plot out the simulation example tree and results.
}
}
\value{
Summary statistic matrix with mean and confidence interval.
}
\references{
Deng, Qiaolin, et al. "Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells." Science 343.6167 (2014): 193-196.
}
\author{
Zilu Zhou
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("DENDROplan_ref")
res=DENDRO.simulation(RD=4.5,n=100,ref=ref,k=4,subtype=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{DENDRO.dist}
\alias{DENDRO.dist}
\title{
DENDRO specific genetic divergence evaluation function
}
\description{
Calculate the cell-to-cell divergence matrix given variant allele read counts, total allele read counts, estimated mutations and sequencing error rate. This method accounts for transcriptional bursting and sequencing error with a Beta-Binomial framework. This function is linear with number of cells and number of mutations.
}
\usage{
DENDRO.dist(X, N, Z, epi = 0.01, show.progress = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  An matrix contains variants allele read counts across cell (column) and loci (row).
}
  \item{N}{
  An matrix contains total allele read counts across cell (column) and loci (row).
}
  \item{Z}{
  An mutation indicator matrix (1 for mutation, 0 for normal) across cell (column) and loci (row).
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.01 according to Illunima.
}
  \item{show.progress}{
  Whether to show the divergence calculation programs. For large and diverse cell population, this function will take some time and we recommend tracking progress.
}
}
\value{
`dist` returns an object of class `"dist"`.
See https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/dist for more details
}
\references{
https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/dist
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{FilterCellMutation}
\alias{FilterCellMutation}
\title{
Filter out low expressed gene and high dropout cells based on read counts
}
\description{
This function filter out low expressed gene (cut.off.VAF) and high dropout cells (cut.off.sd) based on read counts and plot out the summary distributions.
}
\usage{
FilterCellMutation(X, N, Z, Info = NULL, label = NULL, cut.off.VAF = 0.05, cut.off.sd = 5, plot = TRUE)
}
\arguments{
  \item{X}{
  An matrix contains variants allele read counts across cell (column) and loci (row).
}
  \item{N}{
  An matrix contains total allele read counts across cell (column) and loci (row).
}
  \item{Z}{
  An mutation indicator matrix (1 for mutation, 0 for normal) across cell (column) and loci (row).
}
  \item{Info}{
  Mutation loci information. Matrix with row number same to number of loci and column number > 1. This is optional. Default NULL
}
  \item{label}{
  An integer 1D vector that decide the label color by type. The label can be assigned due to piror information such as individuals or site. Default NULL, where no color label. 
}
  \item{cut.off.VAF}{
  Variant allele frequency (VAF, sum(Z)/length(Z) ) filter threshold. Variants with variants allele frequency less than cut.off.VAF (too rare) or greater than 1-cut.off.VAF (too common) are filtered out. Default 0.05
}
  \item{cut.off.sd}{
  Cell total expression (sum(N)) filter threshold. Cells with total expression less than mean(N)-cut.off.sd*sd(N) (too small) or greater than mean(N)+cut.off.sd*sd(N) (double or triplet) are filtered out. Default 5
}
  \item{plot}{
  TRUE or FALSE, decide whether to plot the distribution and cut off (red line). Default value is TRUE.
}
}
\value{
List of filtered matrix
\item{X}{Input X after filtering}
\item{N}{Input N after filtering}
\item{Z}{Input Z after filtering}
\item{Info}{Input Info after filtering}
\item{label}{Input label after filtering}
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{DENDRO.cluster}
\alias{DENDRO.cluster}
\title{
DENDRO specific cluster method on genetic divergence matrix
}
\description{
Cluster cells based on the genetic divergence matrix using hierachical cluster methods. 
}
\usage{
DENDRO.cluster(dist, method = "ward.D", plot = TRUE, label = NULL, ...)
}
\arguments{
  \item{dist}{
  A dissimilarity structure as produced by `dist`
}
  \item{method}{
  Submethod for `hclust`. Default is "ward.D", because statistical integreity
}
  \item{plot}{
  TRUE or FALSE, decide whether to plot the result. Default value is TRUE.
}
  \item{label}{
  An integer 1D vector that decide the label color. Default NULL, where no color label. 
}  
  \item{...}{
  Other arugments can be inherited from plot.phylo in "ape". For example type='fan' for circle cluster plot.
}
}
\value{
An object of class hclust which describes the tree produced by the clustering process. Check https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust for more detail.

If `plot==TRUE`, a hclust plot will also be displayed.
}
\references{
https://www.rdocumentation.org/packages/stats/versions/3.5.1/topics/hclust
}
\author{
Zilu Zhou
%%  ~~who you are~~
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{DENDRO.tree}
\alias{DENDRO.tree}

\title{
Construct a phylogenetic tree
}
\description{
This method constructs a phylogenetic tree with Neighbor Joining method (implemented in `phangorn`), suggesting the evolutionary relationship between different subclones. 
}
\usage{
DENDRO.tree(Z_cluster, label_cluster = NULL)
}
\arguments{
  \item{Z_cluster}{
  An mutation indicator matrix (1 for mutation, 0 for normal) across subclones (column) and loci (row).
}
  \item{label_cluster}{
  Cluster labels, which will be used to color phylogenetic tree node. Better to track cells and clusters.  This is optional. Default NULL
}
}
\value{
 void. 
 A phylogenetic tree will be plotted.
}
\references{
https://cran.r-project.org/web/packages/phangorn/index.html
}
\author{
Zilu Zhou
}
\examples{
data("DENDRO_demo")
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
icd = DENDRO.icd(dist,cluster)
DENDRO_label = cutree(cluster,3)
demo_cluster = DENDRO.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name=c('Cluster3','Cluster2','Cluster1'))
DENDRO.tree(demo_cluster$Z)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\name{DENDRO.recalculate}
\alias{DENDRO.recalculate}
\title{
Recalculate mutation profile for each cluster
}
\description{
This function calculates the mutation profile for each cluster, after pooling cells within each cluster together. Such that, the cells within each cluster have same mutation profiles and the result is more robust. The estimation is based on a maximum likelihood approach. Loci with mutation observed in all or no subclones are removed.
}
\usage{
DENDRO.recalculate(X, N, Info, DENDRO_label, cluster.name = NULL, top = NULL, epi = 0.001, m = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  An matrix contains variants allele read counts across cell (column) and loci (row).
}
  \item{N}{
  An matrix contains total allele read counts across cell (column) and loci (row).}
\item{Info}{
  Mutation loci information. Matrix with row number same to number of loci and column number > 1. This is optional. Default NULL
}
  \item{DENDRO_label}{
  An integer 1D vector that decide the labeled cluster in each cell. Estimated by DENDRO.dist and DENDRO.cluster.
}
  \item{cluster.name}{
  You can label the name of the cluster after observing the cell composition. This is optional. Default NULL
}
  \item{top}{
  Since many mutations may be observed, `top` ask how many top mutation sites you want to selected for downstream analysis based on its marginal likelihood. This is optional. Default NULL
}
  \item{epi}{
  Sequencing error and rare RNA editing combined rate. Default 0.001 according to Illunima.
}
  \item{m}{
  The polidy in maximum likelihood mutation estimation. Default 2.
}
}
\value{
\item{X}{An matrix contains variants allele read counts across subclones (column) and loci (row).}
\item{N}{An matrix contains total allele read counts across subclones (column) and loci (row).}
\item{Z}{An mutation indicator matrix (1 for mutation, 0 for normal) across subclones (column) and loci (row).}
\item{Info}{Input Info after filtering}
\item{Z_cluster_lg}{An likelihood matrix of mutation inference across subclones (column) and loci (row).}
}
\references{
Li, B., et al., A likelihood-based framework for variant calling and de novo mutation detection in families. PLoS Genet, 2012. 8(10): p. e1002944.
}
\author{
Zilu Zhou
}
\examples{
demo_qc = FilterCellMutation(demo$X,demo$N,demo$Z,demo$Info,demo$label)
dist = DENDRO.dist(demo_qc$X,demo_qc$N,demo_qc$Z,show.progress=FALSE)
cluster=DENDRO.cluster(dist,label=demo_qc$label)
icd = DENDRO.icd(dist,cluster)
DENDRO_label = cutree(cluster,3)
demo_cluster = DENDRO.recalculate(demo_qc$X,demo_qc$N, demo_qc$Info, demo_qc$DENDRO_label, cluster.name=c('Cluster3','Cluster2','Cluster1'))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
