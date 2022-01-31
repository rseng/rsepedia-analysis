[![Build Status](https://travis-ci.org/ropensci/cRegulome.svg?branch=master)](https://travis-ci.org/ropensci/cRegulome)
[![codecov](https://codecov.io/gh/ropensci/cRegulome/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/cRegulome)
[![Build status](https://ci.appveyor.com/api/projects/status/gcmojtcsyt7rcwtk?svg=true)](https://ci.appveyor.com/project/ropensci/cregulome)
[![](https://badges.ropensci.org/149_status.svg)](https://github.com/ropensci/onboarding/issues/149)  
[![CRAN version](https://img.shields.io/badge/CRAN-v0.3.0-blue.svg)](https://CRAN.R-project.org/package=cRegulome) 
![downloads](https://cranlogs.r-pkg.org/badges/grand-total/cRegulome)  

# cRegulome  
## Overview  
Transcription factors and microRNAs are important for regulating the gene
expression in normal physiology and pathological conditions. Many
bioinformatics tools were built to predict and identify transcription
factors and microRNA targets and their role in development of diseases
including cancers. The availability of public access high-throughput data
allowed for data-driven predictions and discoveries.
Here, we build on some of these tools and integrative analyses and provide a
tool to access, manage and visualize data from open source databases.
cRegulome provides a programmatic access to the regulome (microRNA and
transcription factor) correlations with target genes in cancer. The package
obtains a local instance of 
[Cistrome Cancer](http://cistrome.org/CistromeCancer/) and 
[miRCancerdb](https://mahshaaban.shinyapps.io/miRCancerdb/) databases and
provides classes and methods to interact with and visualize the correlation
data.  

## What is cRegulome used for?  
cRegulome provides programmatic access to regulome-gene correlation data in 
cancer from different data sources. Researches who are interested in studying 
the role of microRNAs and transcription factors in cancer can use this package 
to construct a small or large scale queries to answer different questions:  

* Which microRNAs and/or transcription factors are associated with a particular
set of genes?  
* What different regulation patterns a microRNA or a transcription factor can 
take in different types of cancer?  
* For a given set of regulatory elements, which genes are likely to be 
regulated by these elements in a certain type of cancer?  

In addition, cRegulome can be used with other R packages like `igraph` to 
study the co-regulation networks in different types of cancer.  
    

## Getting started  
To get starting with cRegulome we show a very quick example. We first start
by downloading a small test database file, make a simple query and convert
the output to a cRegulome object to print and visualize.  

```r
# install the package from CRAN
install.packages('cRegulome')
```

```r
# install the development version from github
devtools::install_github('ropensci/cRegulome')

# install the development version and build vignette from github 
devtools::install_github('ropensci/cRegulome', build_vignettes = TRUE)
```

```{r load_libraries}
# load required libraries
library(cRegulome)
library(RSQLite)
library(ggplot2)
```

```r
if(!file.exists('cRegulome.db')) {
    get_db(test = TRUE)
}

# connect to the db file
conn <- dbConnect(SQLite(), 'cRegulome.db')
```

Or access the same test set file from the package directly  

```r
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)
```

```r
# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
ob <- cmicroRNA(dat)
```

```r
# print object
ob
```

## Documentation

```r
browseVignettes("cRegulome")
```
Alternatively, the vingettes can be found online, [case_study](http://rpubs.com/MahShaaban/cRegulome1) 
and [using_cRegulome](http://rpubs.com/MahShaaban/cRegulome2).

## Citation  

```r
citation("cRegulome")
```

[![](http://www.ropensci.org/public_images/github_footer.png)](http://ropensci.org)

# cRegulome 0.3.2
    
    - Restored on CRAN
    
# cRegulome 0.3.1
    
    - Fixed a CRAN test failure due to changes in a dep pkg

# cRegulome 0.3.0
    
    - Bug fix: since 0.2.0 the argument targets_only did not work properly.
    The bug is fixed and tested in this release.
    - Added the option directed to cor_igraph which allows for constructing
    a directed graph when desired.
    - Added the option to limit the query output of get_tf and get_mir to 
    a predifined set of genes.
  
# cRegulome 0.2.0
    
    - Reduced code dependencies
    - Improved code performance

# cRegulome 0.1.1

    - fix installing in default library tree
    
# cRegulome 0.1.0

    - cRegulome v0.1.0 (2018-02-08) Approved by rOpenSci
    - On CRAN

# cRegulome 0.99.0

    - cRegulome v0.99.0 (2017-09-06) Submit to rOpenSci

# All contributions are most welcomed  

I'd be glad to recieve any comments or ideas to help the package forward.  

## Bugs  

* To report a bug please use the [issue](https://github.com/MahShaaban/cRegulome/issues) page on github  

## Code contributions  

* Fork this repo to your github account  
* Clone the repo to your machine and make changes  
* Push to your account  
* Submit a [pull](https://github.com/MahShaaban/cRegulome/pulls) request at this repo  
    
## Email: [mahmoud.s.fahmy@students.kasralainy.egu.eg](mahmoud.s.fahmy@students.kasralainy.egu.eg)  

## Thank you for contributing!  ## Test environments
* local OS X install, R 3.4.3
* ubuntu 12.04 (on travis-ci), R 3.4.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

## Submissions

* First submission v0.1.0
* Updated the title and description fields
* Updated the description with quoted software names
* Default 'get_db' to write to a temporary directory

---

## Fixed

* Fix installing in default library tree

--- 

## v0.2.0

* Submission v0.2.0
* Reduced code dependencies
* Improved code performance

## v0.3.0

* Submission v0.3.0
* A bug fix
* Minor changes to get_tf, get_mir and cor_igraph
---
title: "Case Study"
author: "Mahmoud Ahmed"
date: "August 22, 2017"
subtitle: Transcription factors and microRNAs in Gastric cancer
vignette: >
    %\VignetteIndexEntry{Case Study} 
    %\VignetteEngine{knitr::rmarkdown} 
    %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE}
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.align = 'center')
```

# Overview   

In this brief study, we attempted to answer questions raised by Shi and colleagues in a paper published in PLOS in 2015. The aim of this case study is to illustrate the use of `cRegulome` to obtain data based on queries of interest. Secondly, to integrate the package in a useful and meaningful workflow and give a context to the kind of information one can get from such databases.  

# Motivation  
Shi et al. studied the transcription factors and microRNA co-regulation of genes involved in gastric cancer to reveal the signaling pathways that drive the development of the disease (Shi et al., 2015). Briefly, they used the previous literature and microarrays to identify the differentially expressed transcription factors (TF) and microRNAs in gastric cancer tissues compared to the normal ones. Then, they identified their target genes using annotation databases and constructed a TF-microRNA gene regulatory network. Finally, they identified the hub-genes and performed a functional analysis namely KEGG pathways enrichment to find which signalling pathways they belong to.  

Here, we tried tackling the same question using `cRegulome` data. We started from the same point like the PLOS paper by using the lists of differentially expressed TF and microRNAs, then used `cRegulome` data to find driving genes and pathways in Stomach and esophageal cancer (STES).  

# PLOS paper data  

## Interesting transcription factors and microRNAs  

We started by obtaining the lists of differentially expressed TFs and microRNAs which were compiled from the previous literature and microarray profiling respectively.  

```{r load_libraries}
library(cRegulome)
library(readxl)
library(ggplot2)
library(RSQLite)
library(R.utils)
library(igraph)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
```

```{r paper_data}
# list of transcription factors
if(!file.exists('tf.xlsx')) 
    download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4393113/bin/pone.0122882.s001.xlsx',
                  destfile = 'tf.xlsx', mode = 'wb')
tf <- read_excel('tf.xlsx', skip = 1)

# list of microRNAs
if(!file.exists('mir.xls')) 
    download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4393113/bin/pone.0122882.s003.xls',
                  destfile = 'mir.xls', mode = 'wb')
mir <- read_excel('mir.xls', skip = 1)
```

Here are the numbers and the first few entries from the lists:  

```{r first_few}
length(unique(tf$SOURCE)); unique(tf$SOURCE) # TFs
length(unique(mir$AccID)); head(unique(mir$AccID), 5) # microRNAs
```

## Use cRegulome to access correlation data  

Here, we show the straight forward way of obtaining correlation/co-expression data using the `cRegulome` package. This is only two simple steps. First, download the database if you are using the package for the first time. And make a query using the TF/microRNAs of interest and limit the output to their known targets. For convenience, the data were included in a test subset of the database file and is used throughout this vignette instead of the full database file.   

```{r cRegulome_data, eval=FALSE}
# download the db file when using it for the first time
destfile = paste(tempdir(), 'cRegulome.db.gz', sep = '/')
if(!file.exists(destfile)) {
    get_db(test = TRUE)
}

# connect to the db file
db_file = paste(tempdir(), 'cRegulome.db', sep = '/')
conn <- dbConnect(SQLite(), db_file)
```

```{r load_testset, include=FALSE}
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)
```

```{r query_database}
# query the database
creg_tf <- get_tf(conn,
                  tf = unique(tf$SOURCE),
                  study = 'STES',
                  targets_only = TRUE)

creg_mir <- get_mir(conn,
                    mir = tolower(unique(mir$AccID)),
                    study = 'STES',
                    targets_only = TRUE)
```

Here is a comparison of the numbers found in the TF/microRNAs previous lists and the query output. 

```{r compare_numbers}
length(unique(creg_mir$mirna_base) %in% unique(tolower(mir$AccID)))
length(unique(creg_tf$tf) %in% unique(tf$SOURCE))
```

# TCGA stomach and esophageal cancer study  

## Transcription factors  

To answer these questions, we first construct a query in cRegulome to get the TF-gene correlations in the STES cancer study. We then look at the numbers of targets, densities and intersections using methods from cRegulome package.  

```{r stes_tf}
# numbers of targets 
table(creg_tf$tf)
```

```{r TF_summary}
# construct a cTF object and plot
ob_tf <- cTF(creg_tf)
cor_joy(ob_tf)
cor_upset(ob_tf)
```

## microRNA   

Similarly, we use the output data.frame of microRNA-gene correlations in the STES study and summarize the numbers, densities and intersections using cRegulome.  

```{r stes_microRNA}
# numbers of targets 
table(creg_mir$mirna_base)
```

```{r microRNA_summary}
# construct a cmicroRNA object and plot
ob_mir <- cmicroRNA(creg_mir)
cor_joy(ob_mir)
cor_upset(ob_mir)
```

## Network construction  

For the purpose of constructing the network, we decided to limit the nodes to the TFs/microRNAs and gene targets with high correlation (absolute Pearson's correlation > 0.3). We first return to `cRegulome` to query the database and tweak the output to be used with the `igraph` package to build the network.  

```{r custom_query}
# query cRegulome to get high correlated targets
creg_tf <- get_tf(conn,
                  tf = unique(tf$SOURCE),
                  study = 'STES',
                  min_abs_cor = .3,
                  targets_only = TRUE)
creg_mir <- get_mir(conn,
                    mir = tolower(unique(mir$AccID)),
                    study = 'STES',
                    min_abs_cor = .3,
                    targets_only = TRUE)
```

First, we construct two separate networks for the TF and the microRNA correlations using the `cor_igraph` function. Then, we combine the two networks and their attributes.  

```{r make_graph}
# make two separate networks
p1 <- cor_igraph(cTF(creg_tf))
p2 <- cor_igraph(cmicroRNA(creg_mir))

# combine networks
p <- graph.union(p1, p2)

# combine attributes
V(p)$type[V(p)$name %in% unique(creg_tf$tf)] <- 'TF'
V(p)$type[V(p)$name %in% unique(creg_mir$mirna_base)] <- 'microRNA'
V(p)$type[is.na(V(p)$type)] <- 'gene'

V(p)$color <- c('lightgreen', 'blue', 'red')[factor(V(p)$type)]

V(p)$label <- ifelse(V(p)$type == 'gene', '', V(p)$name)

E(p)$weight_1[is.na(E(p)$weight_1)] <- E(p)$weight_2[!is.na(E(p)$weight_2)]
```

## Node degrees  

Simple and useful information about the network can be obtained by analyzing the vertices `degree`. A node *degree* is the number of edges it shares with other nodes in the graph. Most of the nodes in the network we constructed have on edge/connection to another node. Most of the gene nodes has one edge and a few genes have 2 to 5 edges. Those are the ones that are regulated by two or more regulatory element (TF/microRNAs).  

```{r node_degrees}
par(mfrow=c(1,2))
deg <- degree(p)

# full network degrees
plot(density(deg), 
     main = 'Full network degrees')

# gene degrees
plot(density(deg[V(p)$type == 'gene']),
     main = 'Gene nodes degrees')
```

Visualizing a dense network may not provide a lot of details, however we notice that the transcription factors (red) and the microRNAs (blue) are in many cases co-regulate one or more gene. So in the following section, we will used a clustering algorithm to capture these connections in sub-communities and used the KEGG enrichment analysis to ask whether they are biologically meaningful.  

```{r plot_network, fig.height=10, fig.width=10}
# plot network
set.seed(123)
par(mfrow=c(1,1))
new_p <- delete.vertices(p, deg < 2)
deg <- degree(new_p)
plot(new_p,
     vertex.size = log2(deg)+1,
     vertex.label.dist = .3,
     vertex.label.cex   = .8,
     vertex.label = V(new_p)$label,
     edge.arrow.size = 0)

legend('bottomleft',
       legend = unique(V(new_p)$type),
       col = unique(V(new_p)$color),
       pch = 19)
```

## Finding clusters  
Here, we tried to find substructures in the network using the fast greedy algorithm. Three clusters were found and are shown in the dendrogram.  

```{r clusters, fig.width=12, fig.height=8}
set.seed(123)
cfg <- cluster_fast_greedy(new_p, weights = E(new_p)$weight_1)
plot_dendrogram(cfg,
                labels = V(new_p)$label,
                mode = 'hclust',
                cex = .5)
```

This is the number of nodes in each cluster.  

```{r numbers_clusters}
clusters <- split(names(membership(cfg)),
                  as.numeric(membership(cfg)))
lengths(clusters)
```

## Pathway enrichment analysis  

```{r kegg_enrichment}
# prepare entrez ids
entrez <- lapply(clusters, function(x) {
    ei <- AnnotationDbi::select(org.Hs.eg.db, x, 'ENTREZID', 'SYMBOL')$ENTREZID
    na.omit(ei)
})

# run kegg enrichment
comp_path <- compareCluster(entrez, fun = 'enrichKEGG', organism = 'hsa')
comp_path@compareClusterResult %>%
    ggplot(aes(x = Description, y = Count)) +
    geom_col() +
    facet_wrap(~Cluster, scales = 'free_y', ncol = 1) +
    coord_flip() +
    labs(x = '')
```

The KEGG pathways enrichment analysis was applied to the `r length(levels(comp_path@compareClusterResult$Cluster))` clusters separately. Clusters `r  as.integer(unique(comp_path@compareClusterResult$Cluster)[1])` and `r  as.integer(unique(comp_path@compareClusterResult$Cluster)[2])` resulted in the enrichment of `r table(as.integer(comp_path@compareClusterResult$Cluster))[1]` and `r table(as.integer(comp_path@compareClusterResult$Cluster))[2]` KEGG pathways respectively. No significant enrichment was found by the other 3 clusters.  

# Remarks  

* In the PLOS paper, the authors used annotation databases to identify the genes co-regulated by TFs and microRNAs and to construct the network. In our analysis, we used a different approach to identify these targets. For TFs, the `targets_only` argument in the call to `get_tf` is based on a ChIP-seq data-driven analysis. For microRNAs, `targets_only` is based on the TargetScan annotations. For both cases, only highly correlated genes (absolute Pearson's correlation > 0.3) were chosen to continue the analysis with.  

* In the published study, the authors used the high `degrees` of the nodes in the networks to define the hubs/or the genes with many outgoing edges to identify the ones which are co-regulated by the TFs and microRNAs of interest. However, we used the `degrees` of the nodes in the network to exclude the genes with few edges. Then, used a clustering algorithm to find highly correlated nodes, then performed the KEGG enrichment analysis on them as separated groups.  

* To sum, we constructed a co-regulation network of TF and microRNA gene targets in stomach and esophageal cancer. We identified a number of genes that are likely to be regulated by these TFs and microRNAs. These genes were clustered in 5 groups. The KEGG pathways enrichment analysis showed a high enrichment of multiple pathways involved in DNA synthesis and repair.  
 

# References  
Shi Y, Wang J, Xin Z, Duan Z, Wang G, Li F. Transcription Factors and microRNA-Co-Regulated Genes in Gastric Cancer Invasion in Ex Vivo. Zhao J-J, editor. PLoS One [Internet]. Public Library of Science; 2015 [cited 2017 Sep 1]; 10: e0122882. doi: 10.1371/journal.pone.0122882.

```{r clean, include=FALSE}
unlink('./*.xls*')
dbDisconnect(conn)
```

---
title: "Using cRegulome"
author: "Mahmoud Ahmed"
date: "August 22, 2017"
vignette: >
    %\VignetteIndexEntry{Using cRegulome}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---

```{r, echo=FALSE}
knitr::opts_chunk$set(message = FALSE, warning = FALSE, fig.align = 'center')
```

# Overview  
Transcription factors and microRNAs are important for regulating the gene
expression in normal physiology and pathological conditions. Many
bioinformatic tools were built to predict and identify transcription
factors and microRNA targets and their role in development of diseases
including cancers. The availability of public access high-throughput data
allowed for data-driven discoveries and  validations of these predictions.
Here, we build on that kind of tools and integrative analyses to provide a
tool to access, manage and visualize data from open source databases.
cRegulome provides a programmatic access to the regulome (microRNA and
transcription factor) correlations with target genes in cancer. The package
obtains a local instance of Cistrome Cancer and miRCancerdb databases and
provides objects and methods to interact with and visualize the correlation
data.  

# Getting started  

To get started with cRegulome, we show a very quick example. We first start
by downloading a small test database file, make a simple query and convert
the output to a cRegulome object to print and visualize.  

```{r load_libraries}
# load required libraries
library(cRegulome)
library(RSQLite)
library(ggplot2)
```

```{r prepare database file, eval=FALSE}
# download the db file when using it for the first time
destfile = paste(tempdir(), 'cRegulome.db.gz', sep = '/')
if(!file.exists(destfile)) {
    get_db(test = TRUE)
}

# connect to the db file
db_file = paste(tempdir(), 'cRegulome.db', sep = '/')
conn <- dbConnect(SQLite(), db_file)
```

```{bash eval=FALSE}
# alternative to downloading the database file
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/9537385/cRegulome.db.gz
gunzip cRegulome.db.gz
```

```{r connect_db, include=FALSE}
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)
```

```{r simple_query}
# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
ob <- cmicroRNA(dat)
```

```{r print_object}
# print object
ob
```

```{r plot_object}
# plot object
cor_plot(ob)
```

# Package Description  

## Data sources  

The two main sources of data used by this package are Cistrome Cancer and
miRCancerdb databases. Cistrome Cancer is based on an integrative analysis of
The Cancer Genome Atlas (TCGA) and public ChIP-seq data. It provides
calculated correlations of (n = 320) transcription factors and their target
genes in (n = 29) cancer study. In addition, Cistrome Cancer provides the
transcription factors regulatory potential to target and non-target genes.
miRCancerdb uses TCGA data and TargetScan annotations to correlate known
microRNAs (n = 750) and target and non-target genes in (n = 25) cancer studies.
  
## Database file  

cRegulome obtains a pre-build SQLite database file of the Cistrome Cancer
and miRCancerdb databases. The details of this build is provided at
[cRegulomedb](https://github.com/MahShaaban/cRegulomedb) in addition to the 
scripts used to pull, format and deposit the data at an on-line repository.
Briefly, the SQLite database consist of 4 tables `cor_mir` and `cor_tf` for 
correlation values and `targets_mir` and `targets_tf` for microRNA miRBase 
ID and transcription factors symbols to genes mappings.  Two indices were 
created to facilitate the database search using the miRBase IDs and 
transcription factors symbols. The database file can be downloaded using the 
function `get_db`.  

To show the details of the database file, the following code connects to
the database and show the names of tables and fields in each of them.  

```{r database_file}
# table names
tabs <- dbListTables(conn)
print(tabs)

# fields/columns in the tables
for(i in seq_along(tabs)) {
  print(dbListFields(conn, tabs[i]))
}
```

## Database query  

To query the database using cRegulome, we provide two main functions;
`get_mir` and `get_tf` for querying microRNA and transcription factors
correlations respectively. Users need to provide the proper IDs for
microRNA, transcription factor symbols and/or TCGA study identifiers.
microRNAs are referred to by the official miRBase IDs, transcription
factors by their corresponding official gene symbols that contains them 
and TCGA studies with their common identifiers. In either cases, the output of 
calling the these functions is a tidy data frame of 4 columns; `mirna_base`/ 
`tf`, `feature`, `cor` and `study` These correspond to the miRBase IDs or
transcription factors symbol, gene symbol, correlation value and the TCGA
study identifier.  

Here we show an example of such a query. Then, we illustrate how this query
is executed on the database using basic `RSQLite` and `dbplyr` which is what
the `get_*` functions are doing.  

```{r database_query}
# query the db for two microRNAs
dat_mir <- get_mir(conn,
                   mir = c('hsa-let-7g', 'hsa-let-7i'),
                   study = 'STES')

# query the db for two transcription factors
dat_tf <- get_tf(conn,
                 tf = c('LEF1', 'MYB'),
                 study = 'STES')

# show first 6 line of each of the data.frames
head(dat_mir); head(dat_tf)
```

## Objects  

Two S3 objects are provided by cRegulome to store and dispatch methods on
the correlation data. cmicroRNA and cTF for microRNA and transcription
factors respectively. The structure of these objects is very similar.
Basically, as all S3 objects, it’s a list of 4 items; microRNA or TF for
the regulome element, features for the gene hits, studies for the TCGA
studies and finally corr is either a `data.frame` when the object has 
data.from a single TCGA study or a named list of data.frames when it has data 
from multiple studies.  Each of these data.frames has the regulome element
(microRNAs or transcription factors) in columns and features/genes in rows.  
To construct these objects, users need to call a constructor function with
the corresponding names on the data.frame output form `get_*`. The reverse
is possible by calling the function `cor_tidy` on the object to get back the 
tidy data.frame.  

```{r cmicroRNA_object}
# explore the cmicroRNA object
ob_mir <- cmicroRNA(dat_mir)
class(ob_mir)
str(ob_mir)
```

```{r cTF_object}
# explore the cTF object
ob_tf <- cTF(dat_tf)
class(ob_tf)
str(ob_tf)
```

## Methods  

cRegulome provides S3 methods to interact a visualize the correlations data
in the cmicroRNA and cTF objects. Table 1 provides an over view of these
functions. These methods dispatch directly on the objects and could be
customized and manipulated in the same way as their generics.  

```{r methods_cmicroRNA}
# cmicroRNA object methods
methods(class = 'cmicroRNA')
```

```{r methods_cTF}
# cTF object methods
methods(class = 'cTF')
```

```{r tidy_method}
# tidy method
head(cor_tidy(ob_mir))
```

```{r cor_hist_method}
# cor_hist method
cor_hist(ob_mir,
     breaks = 100,
     main = '', xlab = 'Correlation')
dev.off()
```


```{r cor_joy_method}
# cor_joy method
cor_joy(ob_mir) +
    labs(x = 'Correlation', y = '')
dev.off()
```

```{r cor_venn_diagram_method}
# cor_venn_diagram method
cor_venn_diagram(ob_mir, cat.default.pos = 'text')
dev.off()
```

```{r cor_upset_method}
# cor_upset method
cor_upset(ob_mir)
dev.off()
```

# Contributions  
Comments, issues and contributions are welcomed at:
[https://github.com/MahShaaban/cRegulome](https://github.com/MahShaaban/cRegulome)  

# Citations  
Please cite:  

```{r citation, eval=FALSE}
citation('cRegulome')
```

```{r clean, echo=FALSE}
dbDisconnect(conn)
unlink('./Venn*')
```

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.tidy.R
\name{cor_tidy}
\alias{cor_tidy}
\title{Tidy \link{cmicroRNA} and \link{cTF} objects}
\usage{
cor_tidy(ob)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}
}
\value{
A tidy \code{data.frame} of four columns. \code{mirna_base} or
\code{tf}is the microRNA miRBase IDs, \code{feature} is the features/genes,
\code{cor} is the corresponding expression correlations and \code{study}
is TCGA study ID.
}
\description{
Tidy \link{cmicroRNA} and \link{cTF} objects
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# convert cmicroRNA object to a tidy data.frame
tidy_cmir <- cor_tidy(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{stat_make_targets}
\alias{stat_make_targets}
\title{Make A SQL statement to extract target features}
\usage{
stat_make_targets(reg, study, type = "mir")
}
\arguments{
\item{reg}{A \code{character} vector of one or more regulator ID.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{type}{A \code{character} string. Either 'mir' of 'tf'. Used to define
columns and tables names.}
}
\value{
A character string
}
\description{
Not meant to be called directly by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.plots.R
\name{cor_plot}
\alias{cor_plot}
\title{Plot method for \link{cmicroRNA} and \link{cTF} objects}
\usage{
cor_plot(ob, study, ...)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{...}{Other options}
}
\value{
A \code{ggplot} object of a dot plot of the correlation values 
between genes and microRNAs or transcription factors in a TCGA study.
}
\description{
A dot plot of microRNA/TF correlation in a single study of TCGA. When the
object \link{cmicroRNA}/\link{cTF} contains more than one TCGA studies, the
argument \code{study} is a requirement.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# print object
cor_plot(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{stat_make}
\alias{stat_make}
\title{Make A SQL statement}
\usage{
stat_make(reg, study, min_abs_cor, max_num, targets, type = "mir")
}
\arguments{
\item{reg}{A \code{character} vector of one or more regulator ID.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{min_abs_cor}{A \code{numeric}, an absolute correlation minimum between 0
and 1 for each \code{mir}.}

\item{max_num}{An \code{integer}, maximum number of \code{features} to show
for each \code{mir} in each \code{study}.}

\item{targets}{A \code{character} vector of gene symbol names.}

\item{type}{A \code{character} string. Either 'mir' of 'tf'. Used to define
columns and tables names.}
}
\value{
A character string
}
\description{
Not meant to be called directly by the user.
}
\examples{
stat_make(reg = 'hsa-let-7g',
          study = 'STES')
          
stat_make(reg = 'hsa-let-7g',
          study = 'STES',
          min_abs_cor = .3)
          
stat_make(reg = 'hsa-let-7g',
          study = 'STES',
          min_abs_cor = .3,
          max_num = 5)
          
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cRegulome.R
\docType{package}
\name{cRegulome}
\alias{cRegulome}
\title{\code{cRegulome} package}
\description{
Download, access and visualize Regulome (microRNA and transcription factors)
data from miRCancer and Cistrome cancer
}
\section{\code{cRegulome} functions to download and query the database file}{

\code{\link{get_db}}
\code{\link{get_tf}}
\code{\link{get_mir}}
}

\section{\code{cRegulome} functions to create S3 objects}{

\code{\link{cTF}}
\code{\link{cmicroRNA}}
}

\section{\code{cRegulome} functions to reshape S3 objects}{

\code{\link{cor_tidy}}
\code{\link{cor_igraph}}
}

\section{\code{cRegulome} functions to visualize data in S3 objects}{

\code{\link{cor_hist}}
\code{\link{cor_joy}}
\code{\link{cor_plot}}
\code{\link{cor_upset}}
\code{\link{cor_venn_diagram}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.plots.R
\name{cor_hist}
\alias{cor_hist}
\title{A histogram of the correlations of microRNA or tf sets}
\usage{
cor_hist(ob, study, ...)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{...}{Other options}
}
\value{
An \code{\link[graphics]{hist}} plot of the correlations values 
between genes a microRNA or a transcription factor in a TCGA study
}
\description{
Plot a \code{\link[graphics]{hist}} of sets of microRNAs or transcription
factors-gene correlations in a TCGA study.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# print object
cor_hist(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.tidy.R
\name{cor_prep}
\alias{cor_prep}
\title{Prepare correlation data for plotting}
\usage{
cor_prep(ob, study, add_dir = TRUE, add_corr = TRUE)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{add_dir}{A \code{logical} default TRUE for whether to add a column
called Direction that has the direction of the correlation; positive or 
negative.}

\item{add_corr}{A \code{logical} default TRUE for whether to add a column
called Correlation that has the absolute value of the correlation}
}
\value{
A \code{data.frame}
}
\description{
Not meant to be called directly by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_tf}
\alias{get_tf}
\title{Get transcription factor correlations from cRegulome.db}
\usage{
get_tf(conn, tf, study, min_abs_cor, max_num, targets_only = FALSE, targets)
}
\arguments{
\item{conn}{A connection such as this returned by 
\code{\link[DBI]{dbConnect}}}

\item{tf}{A required \code{character} vector of the transcription factor of
interest. These are the HUGO official gene symbols of the genes contains the
transcription factor.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{min_abs_cor}{A \code{numeric}, an absolute correlation minimum between 0
and 1 for each \code{mir}.}

\item{max_num}{An \code{integer}, maximum number of \code{features} to show
for each \code{mir} in each \code{study}.}

\item{targets_only}{A \code{logical} whether restrict the output to 
the recognized target features.}

\item{targets}{A \code{character} vector of gene symbol names.}
}
\value{
A tidy \code{data.frame} of four columns. \code{tf} is the official
gene symbols of the genes contains the transcription factor, \code{feature}
is the features/genes, cor is the corresponding expression correlations
and \code{study} is TCGA study ID.
}
\description{
This function access the \code{sqlite} database file which is obtained by
running \link{get_db}. Basically, the function provides ways to query the 
database to the correlation data of the transcription factors of interest. 
The function returns an error if the database file \code{cRegulome.db} is 
not in the working directory.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

\dontrun{
# get transcription factors correlations in all studies
get_tf(conn,
        tf = 'LEF1')
}

# get correlations in a particular study
get_tf(conn,
       tf = 'LEF1',
       study = 'STES')

# enter a custom query with different arguments
get_tf(conn,
       tf = 'LEF1',
       study = 'STES',
       min_abs_cor = .3,
       max_num = 5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.plots.R
\name{cor_venn_diagram}
\alias{cor_venn_diagram}
\title{Venn Diagram of microRNA or transcription factor correlated features}
\usage{
cor_venn_diagram(ob, study, ...)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{...}{Other options}
}
\value{
A venn diagram with a circle or an ellipses for each microRNA and
the number of correlated features.
}
\description{
Count and plot the numbers of microRNA correlated features in
\code{cmicroRNA} object.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# make graph
cor_venn_diagram(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.tidy.R
\name{cor_igraph}
\alias{cor_igraph}
\title{Make an igraph object}
\usage{
cor_igraph(ob, directed = FALSE)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{directed}{A \code{logical} when \code{FALSE} the graph is indirected}
}
\value{
An \code{igraph} object
}
\description{
An \code{igraph} object of from \link{cmicroRNA} or \link{cTF} 
objects.
}
\examples{
# load required libraries
library(RSQLite)
library(cRegulome)

# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- dbConnect(SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# print object
cor_igraph(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.plots.R
\name{cor_joy}
\alias{cor_joy}
\title{A joy plot of correlation of microRNA or tf sets}
\usage{
cor_joy(ob, study, ...)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{...}{Other options}
}
\value{
An \code{\link{ggridges}} plot object
}
\description{
A \code{\link{ggridges}} joy plot of sets of microRNAs or transcription
factors-gene correlations in a TCGA study.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# print object
cor_joy(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_mir}
\alias{get_mir}
\title{Get microRNA correlations from cRegulome.db}
\usage{
get_mir(conn, mir, study, min_abs_cor, max_num, targets_only = FALSE, targets)
}
\arguments{
\item{conn}{A connection such as this returned by 
\code{\link[DBI]{dbConnect}}}

\item{mir}{A required \code{character} vector of the microRNAs of interest.
These are the miRBase ID which are the official identifiers of the
widely used miRBase database, \url{http://www.mirbase.org/}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{min_abs_cor}{A \code{numeric}, an absolute correlation minimum between 0
and 1 for each \code{mir}.}

\item{max_num}{An \code{integer}, maximum number of \code{features} to show
for each \code{mir} in each \code{study}.}

\item{targets_only}{A \code{logical} whether restrict the output to 
the recognized target features.}

\item{targets}{A \code{character} vector of gene symbol names.}
}
\value{
A tidy \code{data.frame} of four columns. \code{mirna_base} is the
microRNA miRBase IDs, \code{feature} is the features/genes, \code{cor}
is the corresponding expression correlations and \code{study} is TCGA
study ID.
}
\description{
This function access the \code{sqlite} database file which is obtained by
running \link{get_db}. Basically, the function provides ways to query the 
database to the correlation data of the microRNAs of interest. The function 
returns an error if the database file \code{cRegulome.db} is not in the 
working directory.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# get microRNA correlations in all studies
get_mir(conn,
        mir = 'hsa-let-7g')

# get correlations in a particular study
get_mir(conn,
        mir = 'hsa-let-7g',
        study = 'STES')

# enter a custom query with different arguments
get_mir(conn,
        mir = 'hsa-let-7g',
        study = 'STES',
        min_abs_cor = .3,
        max_num = 5)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_db}
\alias{get_db}
\title{Get cRegulome.db file}
\usage{
get_db(test = FALSE, destfile, ...)
}
\arguments{
\item{test}{A \code{logical}, default \code{FALSE}. When \code{TRUE}
downloads a database file with the same structure with a subset of
the data for speed.}

\item{destfile}{A character vector for the desired path for the database 
file. By default, when not specified, is constructed by using 
\code{\link{tempdir}} as a directory and the string \code{cRegulome.db.gz}}

\item{...}{Optional arguments passed to \code{\link[utils]{download.file}}}
}
\value{
Downloads a compressed \code{sqlite} file to the current working
directory. The file is named \code{cRegulome.db.gz} by default and it's
not advised to change the name to avoid breaking the other functions
that calls the database.
}
\description{
This function calls \code{\link[utils]{download.file}} to download the
pre-build database file of cRegulome. Additionally, the function checks
the validity of the pre-defined URL and whether the database file exists
in the current working directory to avoid redownloading it. Typically,
users would run this function once at the first time the use the package
or to update the database to the latest version.
}
\examples{
\dontrun{
# download a test set of the database
get_db(test = TRUE)

# download the full database file
get_db(test = FALSE)
}

# load the test db file from shipped with the pacakge
db_file <- system.file("extdata", "cRegulome.db", package = "cRegulome")
file.info(db_file)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{cTF}
\alias{cTF}
\title{Construct cTF object}
\usage{
cTF(dat_tf)
}
\arguments{
\item{dat_tf}{A \code{data.frame} such as this returned by calling
\link{get_tf}.}
}
\value{
An S3 object of class \code{cTF}
}
\description{
Constructs an S3 object called cTF contains data returned by calling
\link{get_tf}. Used to define methods for printing and visualizing
transcription factors-gene expression correlations.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_tf(conn,
              tf = 'LEF1',
              study = 'STES',
              min_abs_cor = .3,
              max_num = 5)

# make a cTF object   
ctf <- cTF(dat)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.plots.R
\name{cor_upset}
\alias{cor_upset}
\title{\code{\link[UpSetR]{upset}} plot of microRNA or tf sets}
\usage{
cor_upset(ob, study, ...)
}
\arguments{
\item{ob}{A \link{cmicroRNA} or \link{cTF} object such as this returned by
calling \link{cmicroRNA} or \link{cTF}.}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{...}{Other options}
}
\value{
An \code{\link[UpSetR]{upset}} plot
}
\description{
\code{\link[UpSetR]{upset}} of sets of microRNAs or transcription
factors and their correlated features in a TCGA study.
}
\examples{
 
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = c('hsa-let-7g', 'hsa-let-7i'),
               study = 'STES')

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

# print object
cor_upset(cmir)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{stat_collect}
\alias{stat_collect}
\title{Collect data from SQLite database}
\usage{
stat_collect(conn, study, stat, type = "mir")
}
\arguments{
\item{conn}{A connection such as this returned by 
\code{\link[DBI]{dbConnect}}}

\item{study}{A \code{character} vector of The Cancer Genome Atlas (TCGA)
study identifiers. To view the available studies in TCGA project,
\url{https://tcga-data.nci.nih.gov/docs/publications/tcga}. When left to
default \code{NULL} all available studies will be included.}

\item{stat}{A string such as this returned by \code{\link{stat_make}}}

\item{type}{A \code{character} string. Either 'mir' of 'tf'. Used to define
columns and tables names.}
}
\value{
A data.frame
}
\description{
Not meant to be called directly by the user.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{cmicroRNA}
\alias{cmicroRNA}
\title{Construct cmicroRNA object}
\usage{
cmicroRNA(dat_mir)
}
\arguments{
\item{dat_mir}{A \code{data.frame} such as this returned by calling
\link{get_mir}.}
}
\value{
An S3 object of class \code{cmicroRNA}
}
\description{
Constructs an S3 object called cmicroRNA contains data returned by calling
\link{get_mir}. Used to define methods for printing and visualizing
microRNA-gene expression correlations.
}
\examples{
# locate the testset file and connect
fl <- system.file('extdata', 'cRegulome.db', package = 'cRegulome')
conn <- RSQLite::dbConnect(RSQLite::SQLite(), fl)

# enter a custom query with different arguments
dat <- get_mir(conn,
               mir = 'hsa-let-7g',
               study = 'STES',
               min_abs_cor = .3,
               max_num = 5)

# make a cmicroRNA object   
cmir <- cmicroRNA(dat)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{stat_collect_targets}
\alias{stat_collect_targets}
\title{Collect target features from SQLite database}
\usage{
stat_collect_targets(conn, stat)
}
\arguments{
\item{conn}{A connection such as this returned by 
\code{\link[DBI]{dbConnect}}}

\item{stat}{A string such as this returned by \code{\link{stat_make}}}
}
\value{
A \code{character} vector
}
\description{
Not meant to be called directly by the user.
}
