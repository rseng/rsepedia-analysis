# scLINE
## 1	Overview
scLINE is an R package for dimension reduction of single-cell RNA-seq data. scLINE integrates the single-cell RNA-seq data and multiple gene networks collated from public databases, supplementing inter-gene interactions to solve the problem of missing information caused by drop-out events.
![](https://github.com/BMILAB/scLINE/raw/master/figure/fig1.jpg)
            **Figure 1.** Overview of scLINE. (A) Construction of cell-gene networks based on scRNA-seq data. (B) Integrate and build gene networks based on multiple genetic interactions databases. (C) In each iteration, the edges in the cell-gene network are selected, and the corresponding low-dimensional vectors of cell and gene are updated based on the expression information and the integrated gene relationships. Train the model to get the low-dimensional representation of cells and genes.
## 2	Installation
You can install scLINE from github with:
```R
install.packages("devtools")
library(devtools)
install_github("BMILAB/scLINE")
library(scLINE)
```
## 3	Preparations
### 3.1	Single-cell RNA-seq data
The input to scLINE is a matrix of scRNA-seq data in which rows correspond to genes and columns correspond to cells. In this study, we take the dataset Usoskin[1] from human neuronal cells as example.
```R
> load(system.file("data","Usoskin.Rdata",package = "scLINE"))
> exp_mat<-Usoskin$rawdata
> dim(exp_mat)
[1] 25334   622
> exp_mat[6:10,1:5]
      L128_B01 L128_C01 L128_D01 L128_E01 L128_F01
AAMP    55.678    303.9   156.57   80.096   185.63
AANAT    0.000      0.0     0.00    0.000     0.00
AARS    83.517    303.9   156.57    0.000     0.00
ABAT    27.839      0.0     0.00    0.000     0.00
ABCA1    0.000      0.0     0.00    0.000     0.00
```
### 3.2	Gene Network 
Users need to enter a list of gene networks. The network needs to be represented in the form of a triplet, with the following three columns: gene ID, gene ID, score. Gene networks can be derived from existing public databases or inter-gene relationship networks obtained by users through other methods. In this study, we take 3 gene networks integrated from the public database String(v9.1)[2], HumanNet(v1)[3] and IntPath(v2.0)[4] as examples. 
```R
> load(system.file("data","string.Rdata",package = "scLINE"))
> load(system.file("data","ppi.Rdata",package = "scLINE"))
> load(system.file("data","humannet.Rdata",package = "scLINE"))
> gene_network<-list(string = string,ppi = ppi,humannet = humannet)
> str(gene_network)
List of 3
 $ string  :'data.frame':	11353056 obs. of  3 variables:
  ..$ V1: chr [1:11353056] "ARF5" "ARF5" "ARF5" "ARF5" ...
  ..$ V2: chr [1:11353056] "PRKCG" "ZNF148" "PRDX6" "KALRN" ...
  ..$ V3: int [1:11353056] 260 164 159 194 164 189 240 164 164 224 ...
 $ ppi     :'data.frame':	39240 obs. of  3 variables:
  ..$ X1: chr [1:39240] "ALDH1A1" "ITGA7" "PPP1R9A" "SRGN" ...
  ..$ X2: chr [1:39240] "ALDH1A1" "CHRNA1" "ACTG1" "CD44" ...
  ..$ X3: num [1:39240] 1 1 1 1 1 1 1 1 1 1 ...
 $ humannet:'data.frame':	476399 obs. of  3 variables:
  ..$ V1: chr [1:476399] "DLST" "MCM2" "FEN1" "EIF4A2" ...
  ..$ V2: chr [1:476399] "OGDH" "MCM3" "PCNA" "EIF4G1" ...
  ..$ V3: num [1:476399] 4.26 4.25 4.24 4.24 4.24 ...
```
## 4	Standard analysis work-flow
Users can obtain the low-dimensional representation matrix of single cells through the function *scLINE*. The user need to specify the low-dimensional vector dimension *L*, the number of iterations *T*, the number of negative samples *K*, the initial learning rate *rho*, and the choice of whether to set the random weight *Random_weight*.
```R
> lowdim_list<-scLINE(exp_mat, gene_network, L = 100, K = 5, T = 3e7, rho = 0.025, Random_weight = TRUE)
Mathching geneID...
Constructing mapping tables for sampling...
Start iterating...
```
Visualize the obtained low dimensional matrix in 2D through the function *visualize*.
```R
> visualize(lowdim_list$cell_low,Usoskin$label)
```
![](https://github.com/BMILAB/scLINE/raw/master/figure/visualization.jpg)

**Figure 2.** Visualization of Usoskin after scLINE
## 5	Session information
```R
> sessionInfo()
R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)
Matrix products: default
locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     
other attached packages:
[1] scLINE_0.1.0
loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1        rstudioapi_0.10   xml2_1.2.0        roxygen2_6.1.1    magrittr_1.5      usethis_1.5.0     devtools_2.1.0    pkgload_1.0.2    
 [9] R6_2.4.0          rlang_0.3.4       stringr_1.4.0     tools_3.6.0       pkgbuild_1.0.3    sessioninfo_1.1.1 cli_1.1.0         withr_2.1.2      
[17] commonmark_1.7    remotes_2.1.0     assertthat_0.2.1  rprojroot_1.3-2   digest_0.6.19     crayon_1.3.4      processx_3.3.1    purrr_0.3.2      
[25] callr_3.2.0       fs_1.3.1          ps_1.3.0          testthat_2.1.1    memoise_1.1.0     glue_1.3.1        stringi_1.4.3     compiler_3.6.0   
[33] desc_1.2.0        backports_1.1.4   prettyunits_1.0.2
```
## References
1.	Usoskin, D., et al., Unbiased classification of sensory neuron types by large-scale single-cell RNA sequencing. Nature Neuroscience, 2015. 18(1): p. 145-+.
2.	Franceschini, A., et al., STRING v9.1: protein-protein interaction networks, with increased coverage and integration. Nucleic Acids Research, 2013. 41(D1): p. D808-D815.
3.	Lee, I., et al., Prioritizing candidate disease genes by network-based boosting of genome-wide association data. Genome Research, 2011. 21(7): p. 1109-1121.
4.	Zhou, H., et al., IntPath-an integrated pathway gene relationship database for model organisms and important pathogens. Bmc Systems Biology, 2012. 6.

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{visualize}
\alias{visualize}
\title{Visualization}
\usage{
visualize(data, label)
}
\arguments{
\item{data}{a matrix in which each row represents a cell}

\item{label}{true labels of the cell corresponding to  data}
}
\value{
Visualization image of the data
}
\description{
Visualize the low-dimensional vectors obtained by scLINE in 2D
}
\examples{
visualize(low_mat$cell_low,Usoskin$label)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scLINE.R
\name{scLINE}
\alias{scLINE}
\title{Integrated networks embedding for scRNA-seq data}
\usage{
scLINE(mat, network, L, K = 5, T, rho = 0.025, Random_weight)
}
\arguments{
\item{mat}{A matrix of scRNA-seq data in which rows represent gene columns representing cells}

\item{network}{A list of gene networks with 3 columns: geneID, geneID, score}

\item{L}{Dim of the low-dimensional representations}

\item{K}{Number of negative samples, the default is 5}

\item{T}{Number of iterations}

\item{rho}{The initial learning rate, the default is 0.025}

\item{Random_weight}{True of False,Whether to add random weights to the network model}
}
\value{
the low-dimensional matrix of input data
}
\description{
A tool of dimension reduction for scRNA-seq data combined with multiple gene networks embedding model
}
\examples{
load(system.file("data","Usoskin.Rdata",package = "scLINE"))
load(system.file("data","ppi.Rdata",package = "scLINE"))
load(system.file("data","humannet.Rdata",package = "scLINE"))
exp_mat<-Usoskin$rawdata
gene_network<-list(ppi = ppi,humannet = humannet)
lowdim_list<-scLINE(exp_mat, gene_network, L = 20, K = 5, T = 1e7, rho = 0.025, Random_weight = TRUE)
}
