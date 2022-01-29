# SDImpute
A R package.
## SDImpute: A statistical block imputation method based on cell-level and gene-level information for dropouts in single-cell RNA-seq data

## Introduction
SDImpute is an implement block imputation for dropout events in scRNA-seq data. SDImpute combines both cell-level and gene-level information to identify the drop-out events and borrows the information unaffected by dropouts from similar cells to impute the dropouts.

## Installation
You can install SDImpute from GitHub with:

``` r
install.packages("devtools")         
library(devtools)           
install_github("Jinsl-lab/SDImpute")
```

## Quick start
The input data of SDImpute is a gene expression matrix, columns and rows represent cells and genes respectively. You can choose whether or not to preprocess the data.

``` r
library(SDImpute)
# An example scRNA-seq dataset, containing 10,166 genes and 425 cells.
data(data)
imputed_data<-SDImpute(data,do.nor=TRUE,do.log=FALSE,auto.k=FALSE,k=5,M=15,T=0.5)
