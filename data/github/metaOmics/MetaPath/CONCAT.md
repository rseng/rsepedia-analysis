# MetaPath
Meta analysis pathway enrichment

## Required packages
* To install all the required packages, open R console
```{r eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("impute")
BiocManager::install("Biobase")
BiocManager::install("GSEABase")
BiocManager::install("genefilter")
BiocManager::install("ConsensusClusterPlus")
BiocManager::install("AnnotationDbi")
BiocManager::install("Rgraphviz")

devtools::install_github("metaOmics/MetaDE")
devtools::install_github("metaOmics/MetaPath")
`
