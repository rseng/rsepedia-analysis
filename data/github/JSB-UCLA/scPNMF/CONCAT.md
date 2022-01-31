# scPNMF
A dimensionality reduction method to facilitate gene selection for targeted gene profiling by learning a sparse gene encoding of single cells

## Introduction
`scPNMF` is a method to facilitate gene selection for targeted gene profiling by learning a sparse gene encoding of single cells. Compared with existing gene selection methods, `scPNMF` has two advantages. First, its selected informative genes can better distinguish cell types, with a small number, e.g., < 200 genes. Second, it enables the alignment of new targeted gene profiling data with reference data in a low-dimensional space to help the prediction of cell types in the new data.

## Installation

`scPNMF` can be installed from Github with the following code in `R`:

``` r
install.packages("devtools")
library(devtools)

install_github("JSB-UCLA/scPNMF")
```

## Usage

For detailed info on `scPNMF` method and applications, please check out the package [vignettes](https://htmlpreview.github.io/?https://github.com/JSB-UCLA/scPNMF/blob/main/inst/docs/scPNMF2.html), or with the following code in `R`: 

``` r
install_github("JSB-UCLA/scPNMF", build_vignettes = TRUE)
browseVignettes("scPNMF")
```

## Contact

Any questions or suggestions on `scPNMF` are welcomed! Please report it on [issues](https://github.com/JSB-UCLA/scPNMF/issues), or contact Dongyuan Song (<dongyuansong@ucla.edu>) or Kexin Li (<aileenlikexin@outlook.com>).

## Reference
Dongyuan Song, Kexin Li, Zachary Hemminger, Roy Wollman, Jingyi Jessica Li, scPNMF: sparse gene encoding of single cells to facilitate gene selection for targeted gene profiling, *Bioinformatics*, Volume 37, Issue Supplement_1, July 2021, Pages i358–i366, https://doi.org/10.1093/bioinformatics/btab273
---
title: "Using scPNMF for gene selection and data projection"
author:
  - name: Dongyuan Song
    affiliation: Bioinformatics IDP, University of California, Los Angeles
  - name: Kexin Li
    affiliation: Department of Statistics, University of California, Los Angeles

date: "`r BiocStyle::doc_date()`"
output: 
  BiocStyle::html_document:
    highlight: pygments
    toc: true
    fig_width: 10
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{scPNMF}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r echo=FALSE, include=FALSE}
knitr::opts_chunk$set(tidy = FALSE, cache = TRUE, dev = "png",
                      message = FALSE, error = FALSE, warning = TRUE)
```

# Introduction 

```{r setup}
suppressPackageStartupMessages(library(scPNMF))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(umap))
suppressPackageStartupMessages(library(ggplot2))
```

*scPNMF* is a method to facilitate gene selection for targeted gene profiling by learning a sparse gene encoding of single cells. Compared with existing gene selection methods, *scPNMF* has two advantages. First, its selected informative genes can better distinguish cell types, with a small number, e.g., < 200 gene. Second, it enables the alignment of new targeted gene profiling data with reference data in a low-dimensional space to help the prediction of cell types in the new data.

In this manual, we will demonstrate how to use `r Githubpkg("JSB-UCLA/scPNMF")` package to select informative genes, and project datasets with limited gene numbers (e.g., targeted gene profiling) into a low-dimensional space.

# Setting Up the Data

Here we use `zheng4` dataset as an example. The `zheng4` dataset is stored as an `SingleCellExperiment` object. Unlike other single-cell methods, *scPNMF* **DO NOT** require normalization (e.g., normalized by cell library size). The input data is simply the raw log-count matrix, where rows are genes and columns are cells.
```{r load-data}
data(zheng4, package = "scPNMF")
Input_zheng4 <- logcounts(zheng4)
```

The input data contains `r dim(Input_zheng4)[1]` genes and `r dim(Input_zheng4)[2]` cells.


# Step I: PNMF

The first step is to run the Projective Nonnegative Matrix Factorization (PNMF). Briefly, PNMF takes advantages from both PCA and NMF (as its name!). Here is a table comparing PNMF with PCA and NMF.

Method|Optimization Problem|Non-negativity|Sparsity|Mutually Exclusiveness|New Data Projection
------|--------------------|--------------|--------|----------------------|-------------------
PNMF | $\min\limits_{\mathbf{W}} \|\mathbf{X} - \mathbf{W}\mathbf{W}^T\mathbf{X}\|\ s.t.\ \mathbf{W}\geq 0$ | Yes | Very high | Very high | Yes
PCA | $\min\limits_{\mathbf{W}} \|\mathbf{X} - \mathbf{W}\mathbf{W}^T\mathbf{X}\|\ s.t.\ \mathbf{W}^T\mathbf{W} = \mathbf{I}$ | No | Low | Low | Yes
NMF | $\min\limits_{\mathbf{W}, \mathbf{H}} \|\mathbf{X} - \mathbf{W}\mathbf{H}\|\ s.t.\ \mathbf{W},\mathbf{H}\geq 0$ | Yes | High | High | No

We run the PNMF algorithm. $K$, the number of low dimension, is a key parameter that needs to be pre-specified by users. We proposed a metric $dev.ortho$ (normalized difference between $W^TW$ and $I$) in Section S1.1 of our *scPNMF* paper to select $K$, which is implemented in `K_selection()` function. However, for the purpose of saving time as well as getting a stable result, we recommend the users to directly set $K=20$ on a typical scRNA-seq data. Here we use $K=15$ for demonstration purpose.
```{r PNMF}
res_pnmf <- scPNMF::PNMFfun(X = Input_zheng4,
                            K = 15, method="EucDist", tol=1e-4, maxIter=1000, verboseN = TRUE)
W <- res_pnmf$Weight
S <- res_pnmf$Score
```

The output is a list containing two matrices: the weight matrix (projection matrix) $\mathbf{W}$ and the score matrix (low-dimensional space) $\mathbf{S} = \mathbf{W}^T\mathbf{X}$. The two output matrices are very similar as those from PCA.

# Step II: Basis Selection

The second main step of *scPNMF* is to select informative bases among the $K$ bases found by PNMF in the last step(i.e., columns of $\mathbf{W}$ and rows of $\mathbf{S}$). The main goal is to remove unwanted variations of cells (e.g., variations irrelevant to cell types). There are three main strategies: functional annotations (optional); data-driven basis selection: correlations with cell library sizes and multimodality tests. We have designed functions for each strategy, and finally we show how to perform basis selection with regard to all the strategies by a combined function.

## Functional Annotations

Due to the nice interpretability of PNMF, each basis briefly represents a funtional gene cluster. Users can perform Gene Ontology (GO) analysis to annoate each basis, and decide which bases they believe are useless and remove irrelevant bases. This is an **OPTIONAL** step since it requires prior knowledge from users. Here we demonstrate functional annotations on the first five bases by setting `dim_use = 1:5`.
```{r annotation}
res_annotation <- scPNMF::basisAnnotate(W = W, 
                                        dim_use = 1:5,
                                        id_type = "ENSEMBL",
                                        return_fig = TRUE)
plot(res_annotation$p_comp)
```




## Data-driven Selection (Correlations with Cell Library Sizes and Multimodality Tests)
The default basis selection procedure in *scPNMF* is purely data-driven. We have two basic assumptions:
  1. A basis which is highly correlated with cell library size (here, total log-counts) should be removed since researchers usually treat cell library size as an unwanted variation;
  2. A basis which shows strong multi-modal pattern should be kept since it indicates cell sub-population (e.g., cell type) structure.

The following function will perform the basis check and users can inspect the checking result by setting `return_fig = TRUE`.

```{r tests}
res_test <- scPNMF::basisTest(S = S, X = Input_zheng4,
                              return_fig = TRUE, ncol = 5,
                              mc.cores = 1)
```


## Combined Basis Selection Function

The following wrapper function will select bases automatically.
```{r basisselect}
W_select <- scPNMF::basisSelect(W = W, S = S,
                                X=Input_zheng4, toTest = TRUE, toAnnotate = FALSE, mc.cores = 1)
colnames(W_select)
```



# Application I: Dimensionality Reduction

*scPNMF* can work as a dimensionality reduction method. Empirically, the result is similar to PCA. Here we use UMAP plot to visualize the reduced space. The cells (dots) are colored by the originally given cell type labels.
```{r umap}
set.seed(123)
S_select <- t(Input_zheng4) %*% W_select
cell_type <- colData(zheng4)$phenoid
umap_select <- data.frame(umap(S_select)$layout, cell_type)
colnames(umap_select) <- c("UMAP1", "UMAP2", "cell_type")

ggplot(umap_select, aes(x=UMAP1, y=UMAP2)) + geom_point(size=0.1, alpha=1, aes(color=cell_type)) +
  ggtitle("UMAP Plot") + theme_bw() +
  theme(aspect.ratio=1, plot.title = element_text(face = "bold", size=12, hjust = 0.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=5)))
```



# Application II: Informative Gene Selection 

The most important application of *scPNMF* is to select informative genes (often called highly-variable genes) in an unsupervised way. Compared to other gene selection methods, *scPNMF* works well with a small informative gene number budget. Here we set the pre-specified gene number as $M = 100$.
```{r infogene}
ig <- getInfoGene(W_select, M = 100, by_basis = FALSE, return_trunW = TRUE, dim_use = NULL)
print(ig$InfoGene)
```

The returned object `ig` contains a informative gene vector and a new projection matrix $\mathbf{W}_{S,(M)}$ which only relies on $M$ informative genes. 



# Application III: Data Projection

$\mathbf{W}_{S,(M)}$ can be used to project new data. Since it only uses $M$ genes, this projection works for new datasets with only $M$ genes (e.g., single-cell targeted gene profiling). Here we illustrate how the projection works with a subset of the original `zheng4` data that only keeps $M$ genes to mimic the structure of new targeted gene profiling data.

```{r projection}
S_new <- getProjection(Input_zheng4[ig$InfoGene, ], ig$trunW)
```

We still use a UMAP plot to visualize the low-dimensional space with $M=100$ genes. 
```{r umap2}
set.seed(456)
cell_type <- colData(zheng4)$phenoid
umap_new <- data.frame(umap(S_new)$layout, cell_type)
colnames(umap_new) <- c("UMAP1", "UMAP2", "cell_type")

ggplot(umap_new, aes(x=UMAP1, y=UMAP2)) + geom_point(size=0.1, alpha=1, aes(color=cell_type)) +
  ggtitle("UMAP Plot on Newly Projected Data") + theme_bw() +
  theme(aspect.ratio=1, plot.title = element_text(face = "bold", size=12, hjust = 0.5), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=5)))
```



# Getting help
If you meet problems when using *scPNMF*, please report it on [issues](https://github.com/JSB-UCLA/scPNMF/issues). For more questions, [email us](dongyuansong@ucla.edu)!



# Session information
```{r}
sessionInfo()
```


# References




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getApp.R
\name{getInfoGene}
\alias{getInfoGene}
\title{Get informative genes}
\usage{
getInfoGene(W, M = 100, by_basis = FALSE, return_trunW = FALSE, dim_use = NULL)
}
\arguments{
\item{W}{The weight matrix output by PNMF}

\item{M}{The user-defined informative gene number}

\item{by_basis}{Return informative genes by basis or not}

\item{return_trunW}{Return the truncated weight matrix or not}

\item{dim_use}{The bases (columns) to be used in the selected weight matrix, \code{NULL} value uses all bases}
}
\value{
A vector containing M top informative genes
}
\description{
Get informative genes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisAnnotate.R
\name{basisAnnotate}
\alias{basisAnnotate}
\title{Annotate the bases of scPNMF}
\usage{
basisAnnotate(
  W,
  DE_prop = 0.1,
  dim_use = NULL,
  id_type = "SYMBOL",
  ont = "BP",
  OrgDb = "org.Hs.eg.db",
  pvalueCutoff = 0.1,
  qvalueCutoff = 0.1,
  minGSSize = 20,
  maxGSSize = 100,
  simp = FALSE,
  return_fig = FALSE,
  cat_num = 10,
  word_num = 50
)
}
\arguments{
\item{W}{The weight matrix from scPNMF.}

\item{DE_prop}{The proportion of genes to be treated as "DE" genes for GO analysis.}

\item{dim_use}{Which bases to annotate.}

\item{id_type}{The gene id types. For example, official gene symbol is \code{SYMBOL}. See details in package \code{clusterProfiler}.}

\item{ont}{Which gene ontology to use. One of "BP", "MF" or "CC". The default is "BP".}

\item{OrgDb}{OrgDb to use.}

\item{pvalueCutoff}{P-value cutoff on on enrichment tests to report. Default is 0.1.}

\item{qvalueCutoff}{Q-value cutoff on on enrichment tests to report. Default is 0.1.}

\item{minGSSize}{Minimal size of genes annotated by Ontology term for testing. Default is 20.}

\item{maxGSSize}{Maximal size of genes annotated by Ontology term for testing. Default is 100.}

\item{simp}{If simplifying the GO terms. Default is FALSE.}

\item{return_fig}{If returning the visualization plots. Default is FALSE.}

\item{cat_num}{How many GO terms will be shown for basis comparison plot.}

\item{word_num}{How many GO terms will be used for word cloud plot.}
}
\description{
Perform Gene Ontology (GO) enrichment analysis on each basis to provide biological interpretation
}
\details{
This function annotates the weight matrix \eqn{W}. Each basis of \eqn{W} represents a functional gene cluster. Users can use this annotation to decide which bases should be kept for further analysis.
See details about the GO analysis in \code{clusterProfiler}.
}
\author{
Dongyuan Song
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zheng4}
\alias{zheng4}
\title{A SingelCellExperiment object containing 2,192 genes and 3,994 cells}
\format{
A SingelCellExperiment object with 2,192 rows (genes) and 3,994 cols (cells).
}
\source{
Duo, A., Robinson, M. D., and Soneson, C. (2018). A systematic performance evaluation of clustering methods for single-cell rna-seq data. F1000Research, 7.
}
\usage{
zheng4
}
\description{
A SingelCellExperiment object containing 2,192 genes and 3,994 cells
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pNMF_func.R
\name{PNMFfun}
\alias{PNMFfun}
\title{Fitting PNMF Models}
\usage{
PNMFfun(
  X,
  K = 10,
  tol = 0.001,
  maxIter = 500,
  verboseN = FALSE,
  zerotol = 1e-10,
  method = "EucDist",
  label = NULL,
  mu = 1,
  lambda = 0.01,
  seed = 123
)
}
\arguments{
\item{X}{Input data matrix, where rows represent features (genes), columns represent samples (cells).}

\item{K}{Specification of the factorization rank (number of low dimension).}

\item{tol}{A threshold below which would be considered converged. Default is 1e-3.}

\item{maxIter}{Number of max iteration times. Default is 500.}

\item{verboseN}{A boolean value indicating whether to print number of iterations.}

\item{zerotol}{A threshold on basis loadings below which would be considered zero. Default is 1e-10.}

\item{method}{A character string indicating which method to be used. One of \code{"EucDist"}, \code{"KL"}, or \code{"DPNMF"}.}

\item{label}{A character vector indicating the cluster type for each cell. Required only when \code{method = "DPNMF"}.}

\item{mu}{A numerical value which controls the penalty term. Larger \code{mu} represents haivier penalization of class distances in \code{DPNMF}. Default is 1.}

\item{lambda}{A numerical value which controls the magnituide of within class distances. Larger \code{lambda} represents larger proportion of within class distances in the total penalty term. Default is 0.01.}

\item{seed}{Random seed of the initialization.}
}
\value{
A list with components:
\describe{
  \item{\code{Weight}}{The basis of model fit (\eqn{W}).}
  \item{\code{Score}}{The mapped scores (\eqn{W^TX}), which is the dimension reduced result.}
}
}
\description{
Fast Projective Nonnegative Matrix Factorization Realizatiton based on Euclidean Distance / KL Divergence / Discriminant pNMF.
}
\details{
Given a data matrix (rows as features and columns as samples), this function
computes the Projective Nonnegative Matrix Factorization (PNMF). Based on different objective functions,
the choices are Euclidean distance (\code{"EucDist"}), KL divergence (\code{"KL"}) (Yang, Zhirong, and Erkki Oja. 2010), or Discriminant PNMF (\code{"DPNMF"}) (Guan, Naiyang, et al. 2013). 
\code{"EucDist"} is supposed to be the most common one;
\code{"KL"} is similar to KL-NMF (Poisson-NMF), and may work better for count data; 
\code{"DPNMF"} requires the predefined labels.

The fitting result of PNMF shares characteristics of both PCA and NMF. The model returns
a \code{basis} matrix, which is similar to the loading matrix in PCA. However, notice that
unlike in PCA the first PC always represents the largest variation, each basis vectors are equal in PNMF.
}
\references{
\itemize{
\item Yang, Z., & Oja, E. (2010). Linear and nonlinear projective nonnegative matrix factorization. IEEE Transactions on Neural Networks, 21(5), 734-749.
\item Guan, N., Zhang, X., Luo, Z., Tao, D., & Yang, X. (2013). Discriminant projective non-negative matrix factorization. PloS one, 8(12), e83291.
\item \url{https://github.com/richardbeare/pNMF}
}
}
\author{
Kexin Li, \email{aileenlikexin@outlook.com}

Dongyuan Song, \email{dongyuansong@g.ucla.edu}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisTest.R
\name{basisTest}
\alias{basisTest}
\title{Data-driven basis selection}
\usage{
basisTest(S, X, return_fig = FALSE, adj_method = "BH", mc.cores = 10, ncol = 4)
}
\arguments{
\item{S}{The score matrix of scPNMF.}

\item{X}{The expression matrix of scPNMF.}

\item{return_fig}{If returning the visualization plots. Default is FALSE.}

\item{adj_method}{Multiple test correction method.}

\item{mc.cores}{The number of cores to use.}

\item{ncol}{Columns for facets in plots.}
}
\description{
Perform data-driven basis selection to remove useless bases and keep informative bases. It contains two approaches:
(1) examine bases by correlations with cell library sizes; (2) examine bases by multimodality tests.
}
\author{
Dongyuan Song
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basisSelect.R
\name{basisSelect}
\alias{basisSelect}
\title{Perform basis selection}
\usage{
basisSelect(
  W,
  S = NULL,
  X = NULL,
  toTest = TRUE,
  cor_thres = 0.7,
  pval_thres = 0.01,
  return_fig = FALSE,
  adj_method = "BH",
  mc.cores = 10,
  ncol = 4,
  toAnnotate = FALSE,
  dim_use = NULL
)
}
\arguments{
\item{W}{The weight matrix output by PNMF}

\item{S}{The score matrix output by PNMF}

\item{X}{The input gene by cell logcount matrix}

\item{toTest}{Whether to select bases by Pearson correlation w/ cell library size and test of multimodality}

\item{cor_thres}{Pearson correlation w/ cell library size cutoff. Default is 0.7.}

\item{pval_thres}{Adjusted p-value cutoff on test of multimodality. Default is 0.01.}

\item{return_fig}{Whether to print scatter plot of score vector against cell library size and distribution of score vectors for each basis.}

\item{adj_method}{P-value correction method. Default is "BH".}

\item{mc.cores}{The number of cores to use for function mclapply().}

\item{ncol}{Columns for facets in plots.}

\item{toAnnotate}{Whether to perform Gene Ontology (GO) enrichment analysis on each basis}

\item{dim_use}{The bases (columns) to be used in the selected weight matrix, \code{NULL} value uses all bases}
}
\value{
The selected weight matrix
}
\description{
Wrapper of basis annotation and data-driven basis selection.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getApp.R
\name{getProjection}
\alias{getProjection}
\title{Get data projection}
\usage{
getProjection(X, W_sM)
}
\arguments{
\item{X}{The input gene by cell logcount matrix}

\item{W_sM}{The truncated, selected weight matrix}
}
\value{
A K0 * cells matrix representing low-dimensional data projection
}
\description{
Get data projection
}
