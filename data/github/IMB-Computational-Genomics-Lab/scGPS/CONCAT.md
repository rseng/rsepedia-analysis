# _scGPS_ - Single Cell Global fate Potential of Subpopulations 
<p align="center">
	<img src="man/figures/scGPSlogo.png" width="200px">
</p>

The _scGPS_ package website is available at: https://imb-computational-genomics-lab.github.io/scGPS/index.html 

The usage instruction can be found at: https://imb-computational-genomics-lab.github.io/scGPS/articles/vignette.html 
## _scGPS_ general description
_scGPS_ is a complete single cell RNA analysis framework from decomposing a mixed population into clusters (_SCORE_) to analysing the relationship between clusters (_scGPS_). _scGPS_ also performs unsupervised selection of predictive genes defining a subpopulation and/or driving transition between subpopulations. 

The package implements two new algorithms _SCORE_ and _scGPS_.

Key features of the _SCORE_ clustering algorithm

- Unsupervised (no prior number of clusters), stable (with automated selection of stability and resolution parameters through scanning a range of search windows for each run, together with a boostrapping aggregation approach to determine stable clusters), fast (with Rcpp implementation)
- _SCORE_ first builds a reference cluster (the highest resolution) and then runs iterative clustering through 40 windows (or more) in the dendrogram
- Resolution is quantified as the divergence from reference by applying adjusted Rand index
- Stability is the proportional to the number of executive runs without Rand index change while changing the cluster search space
- Optimal resolution is the combination of: stable and high resolution
- Bagging algorithm (bootstrap aggregation) can detect a rare subpopulation, which appears multiple times during different decision tree runs 

Key features of the _scGPS_ algorithm

- Estimates transition scores between any two subpopulations
- _scGPS_ prediction model is based on Elastic Net procedure, which enables to select predictive genes and train interpretable models to predict each subpopulation 
- Genes identified by _scGPS_ perform better than known gene markers in predicting cell subpopulations 
- Transition scores are percents of target cells classified as the same class to the original subpopulation 
- For cell subtype comparision, transition scores are similarity between two subpopulations
- The scores are average values from 100 bootstrap runs
- For comparison, a non-shrinkage procedure with linear discriminant analysis (LDA) is used

## _scGPS_ workflow

_scGPS_ takes scRNA expression dataset(s) from one or more unknown sample(s) to find subpopulations and relationship between these subpopulations. The input dataset(s) contains mixed, heterogeous cells. _scGPS_ first uses _SCORE_ (or _CORE_ V2.0) to identify homogenous subpopulations. _scGPS_ contains a number of functions to verify the subpopulations identified by _SCORE_ (e.g. functions to compare with results from PCA, tSNE and the imputation method CIDR). _scGPS_ also has options to find gene markers that distinguish a subpopulation from the remaining cells and performs pathway enrichment analysis to annotate subpopulation. In the second stage, _scGPS_ applies a machine learning procedure to select optimal gene predictors and to build prediction models that can estimate between-subpopulation transition scores, which are the probability of cells from one subpopulation that can likely transition to the other subpopulation.

<p align="center">
	<img src="man/figures/packagePlan.png" width="450px"> <br>
Figure 1. scGPS workflow. Yellow boxes show inputs, and green boxes show main scGPS analysis.  
</p>




---
title: "scGPS introduction"
description: "scGPS is a single cell RNA analysis framework to decompose a mixed population into clusters (SCORE), followed by analysing the relationship between clusters (scGPS)"
author: "Quan Nguyen and Michael Thompson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{single cell Global fate Potential of Subpopulations}
    %\VignetteEngine{knitr::rmarkdown}
    \usepackage[utf8]{inputenc}
---

```{r setup, eval = TRUE, echo = FALSE}
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)
```

# 1. Installation instruction

```{r installation, eval = FALSE}
# To install scGPS from github (Depending on the configuration of the local
# computer or HPC, possible custom C++ compilation may be required - see
# installation trouble-shootings below)
devtools::install_github("IMB-Computational-Genomics-Lab/scGPS")

# for C++ compilation trouble-shooting, manual download and installation can be
# done from github

git clone https://github.com/IMB-Computational-Genomics-Lab/scGPS

# then check in scGPS/src if any of the precompiled (e.g.  those with *.so and
# *.o) files exist and delete them before recompiling

# then with the scGPS as the R working directory, manually install and load
# using devtools functionality
# Install the package
devtools::install()
#load the package to the workspace 
library(scGPS)

```

# 2. A simple workflow of the scGPS: 
The purpose of this workflow is to solve the following task:

* Given a mixed population with known subpopulations, estimate transition scores between these subpopulation.

## 2.1 Create scGPS objects
```{r scGPS_object, warning = FALSE, message = FALSE}

# load mixed population 1 (loaded from day_2_cardio_cell_sample dataset, 
# named it as day2)
library(scGPS)

day2 <- day_2_cardio_cell_sample
mixedpop1 <- new_scGPS_object(ExpressionMatrix = day2$dat2_counts,
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)

# load mixed population 2 (loaded from day_5_cardio_cell_sample dataset, 
# named it as day5)
day5 <- day_5_cardio_cell_sample
mixedpop2 <- new_scGPS_object(ExpressionMatrix = day5$dat5_counts,
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)

```

## 2.2 Run prediction
```{r prediction, warning = FALSE, message = FALSE}

# select a subpopulation
c_selectID <- 1
# load gene list (this can be any lists of user selected genes)
genes <- training_gene_sample
genes <- genes$Merged_unique
# load cluster information 
cluster_mixedpop1 <- colData(mixedpop1)[,1]
cluster_mixedpop2 <- colData(mixedpop2)[,1]
#run training (running nboots = 3 here, but recommend to use nboots = 50-100)
LSOLDA_dat <- bootstrap_prediction(nboots = 3, mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes = genes, c_selectID  = c_selectID,
    listData = list(), cluster_mixedpop1 = cluster_mixedpop1, 
    cluster_mixedpop2 = cluster_mixedpop2, trainset_ratio = 0.7)
names(LSOLDA_dat)
```

## 2.3 Summarise results 
```{r summarise_results, warning = FALSE, message = FALSE}
# summary results LDA
sum_pred_lda <- summary_prediction_lda(LSOLDA_dat = LSOLDA_dat, nPredSubpop = 4)
# summary results Lasso to show the percent of cells
# classified as cells belonging 
sum_pred_lasso <- summary_prediction_lasso(LSOLDA_dat = LSOLDA_dat,
    nPredSubpop = 4)
# plot summary results 
plot_sum <-function(sum_dat){
    sum_dat_tf <- t(sum_dat)
    sum_dat_tf <- na.omit(sum_dat_tf)
    sum_dat_tf <- apply(sum_dat[, -ncol(sum_dat)],1,
        function(x){as.numeric(as.vector(x))})
    sum_dat$names <- gsub("ElasticNet for subpop","sp",  sum_dat$names )
    sum_dat$names <- gsub("in target mixedpop","in p",  sum_dat$names) 
    sum_dat$names <- gsub("LDA for subpop","sp",  sum_dat$names )
    sum_dat$names <- gsub("in target mixedpop","in p",  sum_dat$names)
    colnames(sum_dat_tf) <- sum_dat$names
    boxplot(sum_dat_tf, las=2)
}
plot_sum(sum_pred_lasso)
plot_sum(sum_pred_lda)
# summary accuracy to check the model accuracy in the leave-out test set 
summary_accuracy(object = LSOLDA_dat)
# summary maximum deviance explained by the model 
summary_deviance(object = LSOLDA_dat)
```


# 3. A complete workflow of the scGPS:
The purpose of this workflow is to solve the following task:

* Given an unknown mixed population, find clusters and estimate relationship between clusters
 
## 3.1 Identify clusters in a dataset using CORE
 (skip this step if clusters are known)
```{r CORE, warning = FALSE, message = FALSE}

# find clustering information in an expresion data using CORE
day5 <- day_5_cardio_cell_sample
cellnames <- colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <-data.frame("Cluster"=cluster, "cellBarcodes" = cellnames)
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts,
                    GeneMetadata = day5$dat5geneInfo, CellMetadata = cellnames)

CORE_cluster <- CORE_clustering(mixedpop2, remove_outlier = c(0), PCA=FALSE)

# to update the clustering information, users can ...
key_height <- CORE_cluster$optimalClust$KeyStats$Height
optimal_res <- CORE_cluster$optimalClust$OptimalRes
optimal_index = which(key_height == optimal_res)

clustering_after_outlier_removal <- unname(unlist(
 CORE_cluster$Cluster[[optimal_index]]))
corresponding_cells_after_outlier_removal <- CORE_cluster$cellsForClustering
original_cells_before_removal <- colData(mixedpop2)[,2]
corresponding_index <- match(corresponding_cells_after_outlier_removal,
                            original_cells_before_removal )
# check the matching
identical(as.character(original_cells_before_removal[corresponding_index]),
         corresponding_cells_after_outlier_removal)
# create new object with the new clustering after removing outliers
mixedpop2_post_clustering <- mixedpop2[,corresponding_index]
colData(mixedpop2_post_clustering)[,1] <- clustering_after_outlier_removal

```

## 3.2 Identify clusters in a dataset using SCORE (Stable Clustering at Optimal REsolution)
(skip this step if clusters are known)

(SCORE aims to get stable subpopulation results by introducing bagging aggregation and bootstrapping to the CORE algorithm)
```{r SCORE with bagging, warning = FALSE, message = FALSE}

# find clustering information in an expresion data using SCORE
day5 <- day_5_cardio_cell_sample
cellnames <- colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <-data.frame("Cluster"=cluster, "cellBarcodes" = cellnames)
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts,
                    GeneMetadata = day5$dat5geneInfo, CellMetadata = cellnames )

SCORE_test <- CORE_bagging(mixedpop2, remove_outlier = c(0), PCA=FALSE,
                                bagging_run = 20, subsample_proportion = .8)

```


## 3.3 Visualise all cluster results in all iterations
```{r visualisation, dev='CairoPDF'}
dev.off()
##3.2.1 plot CORE clustering
p1 <- plot_CORE(CORE_cluster$tree, CORE_cluster$Cluster, 
    color_branch = c("#208eb7", "#6ce9d3", "#1c5e39", "#8fca40", "#154975",
        "#b1c8eb"))
p1
#extract optimal index identified by CORE
key_height <- CORE_cluster$optimalClust$KeyStats$Height
optimal_res <- CORE_cluster$optimalClust$OptimalRes
optimal_index = which(key_height == optimal_res)
#plot one optimal clustering bar
plot_optimal_CORE(original_tree= CORE_cluster$tree,
                 optimal_cluster = unlist(CORE_cluster$Cluster[optimal_index]),
                 shift = -2000)

##3.2.2 plot SCORE clustering
#plot all clustering bars
plot_CORE(SCORE_test$tree, list_clusters = SCORE_test$Cluster)
#plot one stable optimal clustering bar
plot_optimal_CORE(original_tree= SCORE_test$tree,
                 optimal_cluster = unlist(SCORE_test$Cluster[
                    SCORE_test$optimal_index]),
                 shift = -100)

```

## 3.4 Compare clustering results with other dimensional reduction methods (e.g., tSNE)
```{r compare_clustering, fig.width=5, fig.height=5}
t <- tSNE(expression.mat=assay(mixedpop2))
p2 <-plot_reduced(t, color_fac = factor(colData(mixedpop2)[,1]),
                      palletes =1:length(unique(colData(mixedpop2)[,1])))
p2
```

## 3.5 Find gene markers and annotate clusters

```{r find_markers, warning = FALSE, message = FALSE, fig.width=7}
#load gene list (this can be any lists of user-selected genes)
genes <-training_gene_sample
genes <-genes$Merged_unique

#the gene list can also be objectively identified by differential expression
#analysis cluster information is requied for find_markers. Here, we use
#CORE results.

#colData(mixedpop2)[,1] <- unlist(SCORE_test$Cluster[SCORE_test$optimal_index])

suppressMessages(library(locfit))

DEgenes <- find_markers(expression_matrix=assay(mixedpop2),
                            cluster = colData(mixedpop2)[,1],
                            selected_cluster=unique(colData(mixedpop2)[,1]))

#the output contains dataframes for each cluster.
#the data frame contains all genes, sorted by p-values
names(DEgenes)

#you can annotate the identified clusters
DEgeneList_1vsOthers <- DEgenes$DE_Subpop1vsRemaining$id

#users need to check the format of the gene input to make sure they are
#consistent to the gene names in the expression matrix

#the following command saves the file "PathwayEnrichment.xlsx" to the
#working dir
#use 500 top DE genes
suppressMessages(library(DOSE))
suppressMessages(library(ReactomePA))
suppressMessages(library(clusterProfiler))
genes500 <- as.factor(DEgeneList_1vsOthers[seq_len(500)])
enrichment_test <- annotate_clusters(genes, pvalueCutoff=0.05, gene_symbol=TRUE)

#the enrichment outputs can be displayed by running
clusterProfiler::dotplot(enrichment_test, showCategory=10, font.size = 6)

```

# 4. Relationship between clusters within one sample or between two samples
The purpose of this workflow is to solve the following task:

* Given one or two unknown mixed population(s) and clusters in each mixed population, estimate
* Visualise relationship between clusters*
 
## 4.1 Start the scGPS prediction to find relationship between clusters
 
```{r scGPS_prediction, warning = FALSE, message = FALSE}

#select a subpopulation, and input gene list
c_selectID <- 1
#note make sure the format for genes input here is the same to the format
#for genes in the mixedpop1 and mixedpop2
genes = DEgenes$id[1:500]

#run the test bootstrap with nboots = 2 runs

cluster_mixedpop1 <- colData(mixedpop1)[,1]
cluster_mixedpop2 <- colData(mixedpop2)[,1]

LSOLDA_dat <- bootstrap_prediction(nboots = 2, mixedpop1 = mixedpop1,
                             mixedpop2 = mixedpop2, genes = genes, 
                             c_selectID  = c_selectID,
                             listData = list(),
                             cluster_mixedpop1 = cluster_mixedpop1,
                             cluster_mixedpop2 = cluster_mixedpop2)

```

## 4.2 Display summary results for the prediction

```{r summarise_prediction}
#get the number of rows for the summary matrix
row_cluster <-length(unique(colData(mixedpop2)[,1]))

#summary results LDA to to show the percent of cells classified as cells
#belonging by LDA classifier
summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster )

#summary results Lasso to show the percent of cells classified as cells
#belonging by Lasso classifier
summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop = row_cluster)

# summary maximum deviance explained by the model during the model training
summary_deviance(object = LSOLDA_dat)

# summary accuracy to check the model accuracy in the leave-out test set
summary_accuracy(object = LSOLDA_dat)

```

## 4.3 Plot the relationship between clusters in one sample
Here we look at one example use case to find relationship between clusters
within one sample or between two sample

```{r prediction_one_sample, warning = FALSE, message = FALSE, fig.width=5}
#run prediction for 3 clusters
cluster_mixedpop1 <- colData(mixedpop1)[,1]
cluster_mixedpop2 <- colData(mixedpop2)[,1]
#cluster_mixedpop2 <- as.numeric(as.vector(colData(mixedpop2)[,1]))

c_selectID <- 1
#top 200 gene markers distinguishing cluster 1
genes = DEgenes$id[1:200]

LSOLDA_dat1 <- bootstrap_prediction(nboots = 2, mixedpop1 = mixedpop2,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop2,
                        cluster_mixedpop2 = cluster_mixedpop2)

c_selectID <- 2
genes = DEgenes$id[1:200]

LSOLDA_dat2 <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop2,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop2,
                        cluster_mixedpop2 = cluster_mixedpop2)

c_selectID <- 3
genes = DEgenes$id[1:200]
LSOLDA_dat3 <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop2,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop2,
                        cluster_mixedpop2 = cluster_mixedpop2)

c_selectID <- 4
genes = DEgenes$id[1:200]
LSOLDA_dat4 <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop2,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop2,
                        cluster_mixedpop2 = cluster_mixedpop2)


#prepare table input for sankey plot

LASSO_C1S2  <- reformat_LASSO(c_selectID=1, mp_selectID = 2,
                             LSOLDA_dat=LSOLDA_dat1,
                             nPredSubpop = length(unique(colData(mixedpop2)
                                [,1])),
                             Nodes_group ="#7570b3")

LASSO_C2S2  <- reformat_LASSO(c_selectID=2, mp_selectID =2,
                             LSOLDA_dat=LSOLDA_dat2,
                             nPredSubpop = length(unique(colData(mixedpop2)
                                [,1])),
                             Nodes_group ="#1b9e77")

LASSO_C3S2  <- reformat_LASSO(c_selectID=3, mp_selectID =2,
                             LSOLDA_dat=LSOLDA_dat3,
                             nPredSubpop = length(unique(colData(mixedpop2)
                                [,1])),
                             Nodes_group ="#e7298a")

LASSO_C4S2  <- reformat_LASSO(c_selectID=4, mp_selectID =2,
                             LSOLDA_dat=LSOLDA_dat4,
                             nPredSubpop = length(unique(colData(mixedpop2)
                                [,1])),
                             Nodes_group ="#00FFFF")

combined <- rbind(LASSO_C1S2,LASSO_C2S2,LASSO_C3S2, LASSO_C4S2 )
combined <- combined[is.na(combined$Value) != TRUE,]

nboots = 2
#links: source, target, value
#source: node, nodegroup
combined_D3obj <-list(Nodes=combined[,(nboots+3):(nboots+4)],
                     Links=combined[,c((nboots+2):(nboots+1),ncol(combined))])

library(networkD3)

Node_source <- as.vector(sort(unique(combined_D3obj$Links$Source)))
Node_target <- as.vector(sort(unique(combined_D3obj$Links$Target)))
Node_all <-unique(c(Node_source, Node_target))

#assign IDs for Source (start from 0)
Source <-combined_D3obj$Links$Source
Target <- combined_D3obj$Links$Target

for(i in 1:length(Node_all)){
   Source[Source==Node_all[i]] <-i-1
   Target[Target==Node_all[i]] <-i-1
}
# 
combined_D3obj$Links$Source <- as.numeric(Source)
combined_D3obj$Links$Target <- as.numeric(Target)
combined_D3obj$Links$LinkColor <- combined$NodeGroup

#prepare node info
node_df <-data.frame(Node=Node_all)
node_df$id <-as.numeric(c(0, 1:(length(Node_all)-1)))

suppressMessages(library(dplyr))
Color <- combined %>% count(Node, color=NodeGroup) %>% select(2)
node_df$color <- Color$color

suppressMessages(library(networkD3))
p1<-sankeyNetwork(Links =combined_D3obj$Links, Nodes = node_df,
                 Value = "Value", NodeGroup ="color", LinkGroup = "LinkColor", 
                 NodeID="Node", Source="Source", Target="Target", fontSize = 22)
p1

#saveNetwork(p1, file = paste0(path,'Subpopulation_Net.html'))

```

## 4.3 Plot the relationship between clusters in two samples
Here we look at one example use case to find relationship between clusters within one sample or between two sample
 
```{r prediction_two_samples,warning = FALSE, message = FALSE, fig.width=5}
#run prediction for 3 clusters
cluster_mixedpop1 <- colData(mixedpop1)[,1]
cluster_mixedpop2 <- colData(mixedpop2)[,1]
row_cluster <-length(unique(colData(mixedpop2)[,1]))

c_selectID <- 1
#top 200 gene markers distinguishing cluster 1
genes = DEgenes$id[1:200]
LSOLDA_dat1 <- bootstrap_prediction(nboots = 2, mixedpop1 = mixedpop1,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop1,
                        cluster_mixedpop2 = cluster_mixedpop2)


c_selectID <- 2
genes = DEgenes$id[1:200]
LSOLDA_dat2 <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop1,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop1,
                        cluster_mixedpop2 = cluster_mixedpop2)

c_selectID <- 3
genes = DEgenes$id[1:200]
LSOLDA_dat3 <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop1,
                        mixedpop2 = mixedpop2, genes=genes, c_selectID, 
                        listData =list(),
                        cluster_mixedpop1 = cluster_mixedpop1,
                        cluster_mixedpop2 = cluster_mixedpop2)

#prepare table input for sankey plot

LASSO_C1S1  <- reformat_LASSO(c_selectID=1, mp_selectID = 1,
                             LSOLDA_dat=LSOLDA_dat1, nPredSubpop = row_cluster, 
                             Nodes_group = "#7570b3")

LASSO_C2S1  <- reformat_LASSO(c_selectID=2, mp_selectID = 1,
                             LSOLDA_dat=LSOLDA_dat2, nPredSubpop = row_cluster, 
                             Nodes_group = "#1b9e77")

LASSO_C3S1  <- reformat_LASSO(c_selectID=3, mp_selectID = 1,
                             LSOLDA_dat=LSOLDA_dat3, nPredSubpop = row_cluster, 
                             Nodes_group = "#e7298a")


combined <- rbind(LASSO_C1S1,LASSO_C2S1,LASSO_C3S1)

nboots = 2
#links: source, target, value
#source: node, nodegroup
combined_D3obj <-list(Nodes=combined[,(nboots+3):(nboots+4)],
                     Links=combined[,c((nboots+2):(nboots+1),ncol(combined))])
combined <- combined[is.na(combined$Value) != TRUE,]


library(networkD3)

Node_source <- as.vector(sort(unique(combined_D3obj$Links$Source)))
Node_target <- as.vector(sort(unique(combined_D3obj$Links$Target)))
Node_all <-unique(c(Node_source, Node_target))

#assign IDs for Source (start from 0)
Source <-combined_D3obj$Links$Source
Target <- combined_D3obj$Links$Target

for(i in 1:length(Node_all)){
   Source[Source==Node_all[i]] <-i-1
   Target[Target==Node_all[i]] <-i-1
}

combined_D3obj$Links$Source <- as.numeric(Source)
combined_D3obj$Links$Target <- as.numeric(Target)
combined_D3obj$Links$LinkColor <- combined$NodeGroup

#prepare node info
node_df <-data.frame(Node=Node_all)
node_df$id <-as.numeric(c(0, 1:(length(Node_all)-1)))

suppressMessages(library(dplyr))
n <- length(unique(node_df$Node))
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
Color = getPalette(n)
node_df$color <- Color
suppressMessages(library(networkD3))
p1<-sankeyNetwork(Links =combined_D3obj$Links, Nodes = node_df,
                 Value = "Value", NodeGroup ="color", LinkGroup = "LinkColor",
                 NodeID="Node", Source="Source", Target="Target", fontSize = 22)
p1
#saveNetwork(p1, file = paste0(path,'Subpopulation_Net.html'))
```



```{r session_info, include=TRUE, echo=TRUE, results='markup'}
devtools::session_info()
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainLassoLDATraining.R
\name{bootstrap_prediction}
\alias{bootstrap_prediction}
\title{BootStrap runs for both scGPS training and prediction}
\usage{
bootstrap_prediction(
  nboots = 1,
  genes = genes,
  mixedpop1 = mixedpop1,
  mixedpop2 = mixedpop2,
  c_selectID = NULL,
  listData = list(),
  cluster_mixedpop1 = NULL,
  cluster_mixedpop2 = NULL,
  trainset_ratio = 0.5,
  LDA_run = TRUE,
  verbose = FALSE,
  log_transform = FALSE
)
}
\arguments{
\item{nboots}{a number specifying how many bootstraps to be run}

\item{genes}{a gene list to build the model}

\item{mixedpop1}{a \linkS4class{SingleCellExperiment} object from a mixed 
population for training}

\item{mixedpop2}{a \linkS4class{SingleCellExperiment} object from a target 
mixed population for prediction}

\item{c_selectID}{the root cluster in mixedpop1 to becompared to clusters in
mixedpop2}

\item{listData}{a \code{list} object, which contains trained results for the
first mixed population}

\item{cluster_mixedpop1}{a vector of cluster assignment for mixedpop1}

\item{cluster_mixedpop2}{a vector of cluster assignment for mixedpop2}

\item{trainset_ratio}{a number specifying the proportion of cells to be part 
of the training subpopulation}

\item{LDA_run}{logical, if the LDA prediction is added to compare to 
ElasticNet}

\item{verbose}{a logical whether to display additional messages}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a \code{list} with prediction results written in to the index 
\code{out_idx}
}
\description{
ElasticNet and LDA prediction for each of all the 
subpopulations in the new mixed population after training the model for a 
subpopulation in the first mixed population. The number of bootstraps to be 
run can be specified.
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
cluster_mixedpop1 <- colData(mixedpop1)[,1]
cluster_mixedpop2 <- colData(mixedpop2)[,1]
c_selectID <- 2
test <- bootstrap_prediction(nboots = 1, mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes=genes, listData =list(), 
    cluster_mixedpop1 = cluster_mixedpop1, 
    cluster_mixedpop2 = cluster_mixedpop2, c_selectID = c_selectID)
names(test)
test$ElasticNetPredict
test$LDAPredict
}
\seealso{
\code{\link{bootstrap_parallel}} for parallel options
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering_bagging.R
\name{CORE_bagging}
\alias{CORE_bagging}
\title{Main clustering SCORE (CORE V2.0) Stable Clustering
at Optimal REsolution with bagging and bootstrapping}
\usage{
CORE_bagging(
  mixedpop = NULL,
  bagging_run = 20,
  subsample_proportion = 0.8,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  remove_outlier = c(0),
  nRounds = 1,
  PCA = FALSE,
  nPCs = 20,
  ngenes = 1500,
  log_transform = FALSE
)
}
\arguments{
\item{mixedpop}{is a \linkS4class{SingleCellExperiment} object from the train
mixed population.}

\item{bagging_run}{an integer specifying the number of bagging runs to be 
computed.}

\item{subsample_proportion}{a numeric specifying the proportion 
of the tree to be chosen in subsampling.}

\item{windows}{a numeric vector specifying the ranges of each window.}

\item{remove_outlier}{a vector containing IDs for clusters to be removed
the default vector contains 0, as 0 is the cluster with singletons.}

\item{nRounds}{an integer specifying the number rounds to attempt to remove 
outliers.}

\item{PCA}{logical specifying if PCA is used before calculating distance 
matrix.}

\item{nPCs}{an integer specifying the number of principal components to use.}

\item{ngenes}{number of genes used for clustering calculations.}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a \code{list} with clustering results of all iterations, and a 
selected
optimal resolution
}
\description{
CORE is an algorithm to generate reproduciable clustering,
CORE is first implemented in ascend R package. Here, CORE V2.0 uses bagging 
analysis to find a stable clustering result and detect rare clusters mixed
population.
}
\examples{
day5 <- day_5_cardio_cell_sample
cellnames<-colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <- data.frame('cluster' = cluster, 'cellBarcodes' = cellnames)
#day5$dat5_counts needs to be in a matrix format
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test <- CORE_bagging(mixedpop2, remove_outlier = c(0), PCA=FALSE,
    bagging_run = 2, subsample_proportion = .7)
}
\author{
Quan Nguyen, 2018-05-11
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{annotate_clusters}
\alias{annotate_clusters}
\title{annotate_clusters functionally annotates the identified clusters}
\usage{
annotate_clusters(
  DEgeneList,
  pvalueCutoff = 0.05,
  gene_symbol = TRUE,
  species = "human"
)
}
\arguments{
\item{DEgeneList}{is a vector of gene symbols, convertable to ENTREZID}

\item{pvalueCutoff}{is a numeric of the cutoff p value}

\item{gene_symbol}{logical of whether the geneList is a gene symbol}

\item{species}{is the selection of 'human' or 'mouse', default to 'human' 
genes}
}
\value{
write enrichment test output to a file and an enrichment test object 
for plotting
}
\description{
often we need to label clusters with unique biological 
characters. One of the common approach to annotate a cluster is to perform 
functional enrichment analysis. The annotate implements ReactomePA and
clusterProfiler for this analysis type in R. The function require 
installation of several databases as described below.
}
\examples{
genes <-training_gene_sample
genes <-genes$Merged_unique[seq_len(50)]
enrichment_test <- annotate_clusters(genes, pvalueCutoff=0.05, 
    gene_symbol=TRUE, species = 'human')
clusterProfiler::dotplot(enrichment_test, showCategory=15)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{CORE_clustering}
\alias{CORE_clustering}
\title{Main clustering CORE V2.0 updated}
\usage{
CORE_clustering(
  mixedpop = NULL,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  remove_outlier = c(0),
  nRounds = 1,
  PCA = FALSE,
  nPCs = 20,
  ngenes = 1500,
  verbose = FALSE,
  log_transform = FALSE
)
}
\arguments{
\item{mixedpop}{is a \linkS4class{SingleCellExperiment} object from the train
mixed population}

\item{windows}{a numeric specifying the number of windows to test}

\item{remove_outlier}{a vector containing IDs for clusters to be removed
the default vector contains 0, as 0 is the cluster with singletons.}

\item{nRounds}{an integer specifying the number rounds to attempt to remove 
outliers.}

\item{PCA}{logical specifying if PCA is used before calculating distance
matrix}

\item{nPCs}{an integer specifying the number of principal components to use.}

\item{ngenes}{number of genes used for clustering calculations.}

\item{verbose}{a logical whether to display additional messages}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a \code{list} with clustering results of all iterations, and a 
selected optimal resolution
}
\description{
CORE is an algorithm to generate reproduciable clustering,
CORE is first implemented in ascend R package. Here, CORE V2.0 introduces 
several new functionalities, including three key features:
fast (and more memory efficient) implementation with C++ and paralellisation
options allowing clustering of hundreds of thousands of cells 
(ongoing development), outlier revomal important if singletons exist (done),
a number of dimensionality reduction methods including the imputation
implementation (CIDR) for confirming clustering results (done), and an option
to select the number of optimisation tree height windows for increasing
resolution
}
\examples{
day5 <- day_5_cardio_cell_sample
#day5$dat5_counts needs to be in a matrix format
cellnames <- colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <-data.frame('Cluster'=cluster, 'cellBarcodes' = cellnames)
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test <- CORE_clustering(mixedpop2, remove_outlier = c(0), PCA=FALSE, nPCs=20,
    ngenes=1500)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calcDistArma}
\alias{calcDistArma}
\title{Compute Euclidean distance matrix by rows}
\usage{
calcDistArma(x)
}
\arguments{
\item{x}{A numeric matrix}
}
\value{
a distance matrix
}
\description{
Compute Euclidean distance matrix by rows
}
\examples{
mat_test <-matrix(rnbinom(1000,mu=0.01, size=10),nrow=1000)
#library(microbenchmark)
#microbenchmark(calcDistArma(mat_test), dist(mat_test), times=3)


}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{sub_clustering}
\alias{sub_clustering}
\title{sub_clustering for selected cells}
\usage{
sub_clustering(
  object = NULL,
  ngenes = 1500,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  select_cell_index = NULL
)
}
\arguments{
\item{object}{is a \linkS4class{SingleCellExperiment} object from the
train mixed population}

\item{ngenes}{number of genes used for clustering calculations.}

\item{windows}{a numeric vector specifying the ranges of each window.}

\item{select_cell_index}{a vector containing indexes for cells in selected 
clusters to be reclustered}
}
\value{
clustering results
}
\description{
performs 40 clustering runs or more depending on windows
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test_sub_clustering <-sub_clustering(mixedpop2,
    select_cell_index = c(seq_len(100)))
}
\author{
Quan Nguyen, 2018-01-31
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_prediction_summary.R
\name{summary_accuracy}
\alias{summary_accuracy}
\title{get percent accuracy for Lasso model, from \code{n} bootstraps}
\usage{
summary_accuracy(object = NULL)
}
\arguments{
\item{object}{is a list containing the training results from 
the \code{summary_accuracy} summarise \code{n} bootstraps}
}
\value{
a vector of percent accuracy for the selected subpopulation
}
\description{
The training results from \code{training} were written to
}
\examples{
c_selectID<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
LSOLDA_dat <- bootstrap_prediction(nboots = 1,mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
    cluster_mixedpop1 = colData(mixedpop1)[,1],
    cluster_mixedpop2=colData(mixedpop2)[,1])
summary_accuracy(LSOLDA_dat)
summary_deviance(LSOLDA_dat)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_prediction_summary.R
\name{summary_prediction_lda}
\alias{summary_prediction_lda}
\title{get percent deviance explained for LDA model, from \code{n} bootstraps}
\usage{
summary_prediction_lda(LSOLDA_dat = NULL, nPredSubpop = NULL)
}
\arguments{
\item{LSOLDA_dat}{is a list containing the training results from 
\code{training}}

\item{nPredSubpop}{is the number of subpopulations in the target mixed 
population}
}
\value{
a dataframe containg information for the LDA prediction 
results, each column contains prediction results for all subpopulations from
each bootstrap run
}
\description{
the training results from \code{training} were written to
the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises 
prediction explained for \code{n} bootstrap runs and also returns the best
deviance matrix for plotting, as well as the best matrix with Lasso genes
and coefficients
}
\examples{
c_selectID<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
LSOLDA_dat <- bootstrap_prediction(nboots = 1,mixedpop1 = mixedpop1, 
mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
    cluster_mixedpop1 = colData(mixedpop1)[,1],
    cluster_mixedpop2=colData(mixedpop2)[,1])
summary_prediction_lda(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/training_gene_sample.R
\docType{data}
\name{training_gene_sample}
\alias{training_gene_sample}
\title{Input gene list for training \pkg{scGPS}, e.g. differentially 
expressed genes}
\format{
a list instance, containing a count matrix and a vector with 
clustering information.
}
\source{
Dr Joseph Powell's laboratory, IMB, UQ
}
\usage{
training_gene_sample
}
\value{
a vector containing gene symbols
}
\description{
Genes should be ordered from most to least important genes
(1 row per gene)
}
\author{
Quan Nguyen, 2017-11-25
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{CORE_subcluster}
\alias{CORE_subcluster}
\title{sub_clustering (optional) after running CORE 'test'}
\usage{
CORE_subcluster(
  mixedpop = NULL,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  select_cell_index = NULL,
  ngenes = 1500
)
}
\arguments{
\item{mixedpop}{is a \linkS4class{SingleCellExperiment} object from the train
mixed population}

\item{windows}{a numeric specifying the number of windows to test}

\item{select_cell_index}{a vector containing indexes for cells in selected 
clusters to be reclustered}

\item{ngenes}{number of genes used for clustering calculations.}
}
\value{
a \code{list} with clustering results of all iterations, and a 
selected optimal resolution
}
\description{
CORE_subcluster allows re-cluster the CORE clustering 
result
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts,
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test <- CORE_clustering(mixedpop2,remove_outlier= c(0))
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{mean_cpp}
\alias{mean_cpp}
\title{Calculate mean}
\usage{
mean_cpp(x)
}
\arguments{
\item{x}{integer.}
}
\value{
a scalar value
}
\description{
Calculate mean
}
\examples{
mean_cpp(seq_len(10^6))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{subset_cpp}
\alias{subset_cpp}
\title{Subset a matrix}
\usage{
subset_cpp(m1in, rowidx_in, colidx_in)
}
\arguments{
\item{m1in}{an R matrix (expression matrix)}

\item{rowidx_in}{a numeric vector of rows to keep}

\item{colidx_in}{a numeric vector of columns to keep}
}
\value{
a subsetted matrix
}
\description{
Subset a matrix
}
\examples{
mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=100)
subset_mat <- subset_cpp(mat_test, rowidx_in=c(1:10), colidx_in=c(100:500))
dim(subset_mat)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{PrinComp_cpp}
\alias{PrinComp_cpp}
\title{Principal component analysis}
\usage{
PrinComp_cpp(X)
}
\arguments{
\item{X}{an R matrix (expression matrix), rows are genes, columns are cells}
}
\value{
a list with three list pca lists
}
\description{
This function provides significant speed gain if the input
matrix is big
}
\examples{
mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=1000)
#library(microbenchmark)
#microbenchmark(PrinComp_cpp(mat_test), prcomp(mat_test), times=3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/day_5_cardio_cell_sample.R
\docType{data}
\name{day_5_cardio_cell_sample}
\alias{day_5_cardio_cell_sample}
\title{One of the two example single-cell count matrices to be used
for \pkg{scGPS} prediction}
\format{
a list instance, containing a count matrix and a vector with 
clustering information.
}
\source{
Dr Joseph Powell's laboratory, IMB, UQ
}
\usage{
day_5_cardio_cell_sample
}
\value{
a list, with the name day5
}
\description{
The count data set contains counts for 17402 genes for 983 cells
(1 row per gene) randomly subsampled from day-5 cardio-differentiation 
population
The vector of clustering information contains corresponding to cells in the
count matrix
}
\author{
Quan Nguyen, 2017-11-25
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainLassoLDATraining.R
\name{training}
\alias{training}
\title{Main model training function for finding the best model
that characterises a subpopulation}
\usage{
training(
  genes = NULL,
  cluster_mixedpop1 = NULL,
  mixedpop1 = NULL,
  mixedpop2 = NULL,
  c_selectID = NULL,
  listData = list(),
  out_idx = 1,
  standardize = TRUE,
  trainset_ratio = 0.5,
  LDA_run = FALSE,
  log_transform = FALSE
)
}
\arguments{
\item{genes}{a vector of gene names (for ElasticNet shrinkage); gene symbols
must be in the same format with gene names in subpop2. Note that genes are 
listed by the order of importance, e.g. differentially expressed genes that 
are most significan, so that if the gene list contains too many genes, only 
the top 500 genes are used.}

\item{cluster_mixedpop1}{a vector of cluster assignment in mixedpop1}

\item{mixedpop1}{is a \linkS4class{SingleCellExperiment} object from the 
train mixed population}

\item{mixedpop2}{is a \linkS4class{SingleCellExperiment} object from the 
target mixed population}

\item{c_selectID}{a selected number to specify which subpopulation to be used
for training}

\item{listData}{list to store output in}

\item{out_idx}{a number to specify index to write results into the list 
output. This is needed for running bootstrap.}

\item{standardize}{a logical value specifying whether or not to standardize 
the train matrix}

\item{trainset_ratio}{a number specifying the proportion of cells to be part
of the training subpopulation}

\item{LDA_run}{logical, if the LDA run is added to compare to ElasticNet}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a \code{list} with prediction results written in to the indexed 
\code{out_idx}
}
\description{
Training a haft of all cells to find optimal ElasticNet and LDA
models to predict a subpopulation
}
\examples{

c_selectID<-1
out_idx<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts,
GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
listData  <- training(genes, 
    cluster_mixedpop1 = colData(mixedpop1)[, 1],
    mixedpop1 = mixedpop1, mixedpop2 = mixedpop2, c_selectID,
    listData =list(), out_idx=out_idx, trainset_ratio = 0.5)
names(listData)
listData$Accuracy
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainLassoLDATraining.R
\name{bootstrap_parallel}
\alias{bootstrap_parallel}
\title{BootStrap runs for both scGPS training and prediction
with parallel option}
\usage{
bootstrap_parallel(
  ncores = 4,
  nboots = 1,
  genes = genes,
  mixedpop1 = mixedpop1,
  mixedpop2 = mixedpop2,
  c_selectID,
  listData = list(),
  cluster_mixedpop1 = NULL,
  cluster_mixedpop2 = NULL
)
}
\arguments{
\item{ncores}{a number specifying how many cpus to be used for running}

\item{nboots}{a number specifying how many bootstraps to be run}

\item{genes}{a gene list to build the model}

\item{mixedpop1}{a \linkS4class{SingleCellExperiment} object from a mixed 
population for training}

\item{mixedpop2}{a \linkS4class{SingleCellExperiment} object from a target 
mixed population for prediction}

\item{c_selectID}{the root cluster in mixedpop1 to becompared to clusters in
mixedpop2}

\item{listData}{a \code{list} object, which contains trained results for the
first mixed population}

\item{cluster_mixedpop1}{a vector of cluster assignment for mixedpop1}

\item{cluster_mixedpop2}{a vector of cluster assignment for mixedpop2}
}
\value{
a \code{list} with prediction results written in to the index 
\code{out_idx}
}
\description{
same as bootstrap_prediction, but with an multicore option
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
#prl_boots <- bootstrap_parallel(ncores = 4, nboots = 1, genes=genes,
#    mixedpop1 = mixedpop2, mixedpop2 = mixedpop2,  c_selectID=1,
#    listData =list())
#prl_boots[[1]]$ElasticNetPredict
#prl_boots[[1]]$LDAPredict

}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{find_markers}
\alias{find_markers}
\title{find marker genes}
\usage{
find_markers(
  expression_matrix = NULL,
  cluster = NULL,
  selected_cluster = NULL,
  fitType = "local",
  dispersion_method = "per-condition",
  sharing_Mode = "maximum"
)
}
\arguments{
\item{expression_matrix}{is  a normalised expression matrix.}

\item{cluster}{corresponding cluster information in the expression_matrix
by running CORE clustering or using other methods.}

\item{selected_cluster}{a vector of unique cluster ids to calculate}

\item{fitType}{string specifying 'local' or 'parametric' for DEseq dispersion
estimation}

\item{dispersion_method}{one of the options c( 'pooled', 'pooled-CR', 
per-condition', 'blind' )}

\item{sharing_Mode}{one of the options c("maximum", "fit-only", 
"gene-est-only")}
}
\value{
a \code{list} containing sorted DESeq analysis results
}
\description{
Find DE genes from comparing one clust vs remaining
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
# depending on the data, the DESeq::estimateDispersions function requires
# suitable fitType
# and dispersion_method options
DEgenes <- find_markers(expression_matrix=assay(mixedpop1),
                        cluster = colData(mixedpop1)[,1],
                        selected_cluster=c(1), #can also run for more
                        #than one clusters, e.g.selected_cluster = c(1,2)
                        fitType = "parametric", 
                        dispersion_method = "blind",
                        sharing_Mode="fit-only"
                        )
names(DEgenes)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calcDist}
\alias{calcDist}
\title{Compute Euclidean distance matrix by rows}
\usage{
calcDist(x)
}
\arguments{
\item{x}{A numeric matrix}
}
\value{
a distance matrix
}
\description{
Compute Euclidean distance matrix by rows
}
\examples{
mat_test <-matrix(rnbinom(1000,mu=0.01, size=10),nrow=1000)
calcDist(mat_test)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_object.R
\name{new_summarized_scGPS_object}
\alias{new_summarized_scGPS_object}
\title{new_summarized_scGPS_object}
\usage{
new_summarized_scGPS_object(
  ExpressionMatrix = NULL,
  GeneMetadata = NULL,
  CellMetadata = NULL
)
}
\arguments{
\item{ExpressionMatrix}{An expression dataset in matrix format.
Rows should represent a transcript and its normalised counts,
while columns should represent individual cells.}

\item{GeneMetadata}{A data frame or vector containing gene identifiers used 
in the expression matrix. The first column should hold the cell identifiers
you are using in the expression matrix. Other columns contain information 
about the genes, such as their corresponding ENSEMBL transcript identifiers.}

\item{CellMetadata}{A data frame containing cell identifiers 
(usually barcodes) and clustering information (the first column of the data
frame contains clustering information). The column containing clustering
information needs to be named as 'Cluster'. If clustering information is not
available, users can run CORE function and add the information to the scGPS
before running scGPS prediction}
}
\value{
This function generates an scGPS object belonging to the
\linkS4class{SingleCellExperiment}.
}
\description{
\code{\link{new_scGPS_object}} generates a scGPS object in the 
\linkS4class{SingleCellExperiment} class for use with the scGPS package. This
object contains an expression matrix, associated metadata (cells, genes,
clusters). The data are expected to be normalised counts.
}
\examples{
day2 <- day_2_cardio_cell_sample
t <-new_summarized_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
colData(t); show(t); colnames(t)
}
\seealso{
\linkS4class{SingleCellExperiment}
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MainLassoLDATraining.R
\name{predicting}
\alias{predicting}
\title{Main prediction function applying the optimal ElasticNet and LDA models}
\usage{
predicting(
  listData = NULL,
  cluster_mixedpop2 = NULL,
  mixedpop2 = NULL,
  out_idx = NULL,
  standardize = TRUE,
  LDA_run = FALSE,
  c_selectID = NULL,
  log_transform = FALSE
)
}
\arguments{
\item{listData}{a \code{list} object containing trained results for the
selected subpopulation in the first mixed population}

\item{cluster_mixedpop2}{a vector of cluster assignment for mixedpop2}

\item{mixedpop2}{a \linkS4class{SingleCellExperiment} object from the target
mixed population of importance, e.g. differentially expressed genes that are
most significant}

\item{out_idx}{a number to specify index to write results into the list 
output. This is needed for running bootstrap.}

\item{standardize}{a logical of whether to standardize the data}

\item{LDA_run}{logical, if the LDA prediction is added to compare to 
ElasticNet, the LDA model needs to be trained from the training before
inputting to this prediction step}

\item{c_selectID}{a number to specify the trained cluster used for prediction}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a \code{list} with prediction results written in to the index
\code{out_idx}
}
\description{
Predict a new mixed population after training the model for a
subpopulation in the first mixed population.
All subpopulations in the new target mixed population will be predicted, 
where each targeted subpopulation will have a transition score from the 
orginal subpopulation to the new subpopulation.
}
\examples{
c_selectID<-1
out_idx<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
listData  <- training(genes, 
    cluster_mixedpop1 = colData(mixedpop1)[, 1], mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, c_selectID, listData =list(), out_idx=out_idx)
listData  <- predicting(listData =listData,  mixedpop2 = mixedpop2, 
    out_idx=out_idx, cluster_mixedpop2 = colData(mixedpop2)[, 1], 
    c_selectID = c_selectID)

}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_object.R
\name{new_scGPS_object}
\alias{new_scGPS_object}
\title{new_scGPS_object}
\usage{
new_scGPS_object(
  ExpressionMatrix = NULL,
  GeneMetadata = NULL,
  CellMetadata = NULL,
  LogMatrix = NULL
)
}
\arguments{
\item{ExpressionMatrix}{An expression matrix in data.frame or matrix format.
Rows should represent a transcript and its normalised counts,
while columns should represent individual cells.}

\item{GeneMetadata}{A data frame or vector containing gene identifiers used 
in the expression matrix. The first column should hold the gene identifiers
you are using in the expression matrix. Other columns contain information 
about the genes, such as their corresponding ENSEMBL transcript identifiers.}

\item{CellMetadata}{A data frame containing cell identifiers 
(usually barcodes) and an integer representing which batch they belong to.
The column containing clustering information needs to be the first column in 
the CellMetadata dataframe If clustering information is not available, users
can run CORE function and add the information to the scGPS before running
scGPS prediction}

\item{LogMatrix}{optional input for a log matrix of the data. If no log
matrix is supplied one will be created for the object}
}
\value{
This function generates an scGPS object belonging to the 
\linkS4class{SingleCellExperiment}.
}
\description{
\code{\link{new_scGPS_object}} generates a scGPS object in the 
\linkS4class{SingleCellExperiment} class for use with the scGPS package. This
object contains an expression matrix, associated metadata (cells, genes,
clusters). The data are expected to be normalised counts.
}
\examples{
day2 <- day_2_cardio_cell_sample
t <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
colData(t); show(t); colnames(t)
}
\seealso{
\linkS4class{SingleCellExperiment}
}
\author{
Quan Nguyen, 2018-04-06
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_prediction_summary.R
\name{summary_prediction_lasso}
\alias{summary_prediction_lasso}
\title{get percent deviance explained for Lasso model, from \code{n} 
bootstraps}
\usage{
summary_prediction_lasso(LSOLDA_dat = NULL, nPredSubpop = NULL)
}
\arguments{
\item{LSOLDA_dat}{is a list containing the training results from 
\code{training}}

\item{nPredSubpop}{is the number of subpopulations in the target mixed 
population}
}
\value{
a dataframe containg information for the Lasso prediction 
results, each column
contains prediction results for all subpopulations from each bootstrap run
}
\description{
the training results from \code{training} were written to
the object \code{LSOLDA_dat}, the \code{summary_prediction} summarises 
prediction for \code{n} bootstrap runs
}
\examples{
c_selectID<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
LSOLDA_dat <- bootstrap_prediction(nboots = 1,mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
    cluster_mixedpop1 = colData(mixedpop1)[,1],
    cluster_mixedpop2=colData(mixedpop2)[,1])
summary_prediction_lasso(LSOLDA_dat=LSOLDA_dat, nPredSubpop=4)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_prediction_summary.R
\name{summary_deviance}
\alias{summary_deviance}
\title{get percent deviance explained for Lasso model, 
from \code{n} bootstraps}
\usage{
summary_deviance(object = NULL)
}
\arguments{
\item{object}{is a list containing the training results from 
\code{training}}
}
\value{
a \code{list} containing three elements, with a vector of percent
maximum deviance explained, a dataframe containg information for the full 
deviance, and a dataframe containing gene names and coefficients of the best
model
}
\description{
the training results from \code{training} were written to
the object \code{LSOLDA_dat}, the \code{summary_devidance} summarises 
deviance explained for \code{n} bootstrap runs and also returns the best
deviance matrix for plotting, as well as the best matrix with Lasso genes 
and coefficients
}
\examples{
c_selectID<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo,
                    CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
LSOLDA_dat <- bootstrap_prediction(nboots = 2,mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
    cluster_mixedpop1 = colData(mixedpop1)[,1],
    cluster_mixedpop2=colData(mixedpop2)[,1])
summary_deviance(LSOLDA_dat)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scgps_prediction_summary.R
\name{reformat_LASSO}
\alias{reformat_LASSO}
\title{summarise bootstrap runs for Lasso model, from \code{n} bootstraps}
\usage{
reformat_LASSO(
  c_selectID = NULL,
  mp_selectID = NULL,
  LSOLDA_dat = NULL,
  nPredSubpop = NULL,
  Nodes_group = "#7570b3",
  nboots = 2
)
}
\arguments{
\item{c_selectID}{is the original cluster to be projected}

\item{mp_selectID}{is the target mixedpop to project to}

\item{LSOLDA_dat}{is the results from the bootstrap}

\item{nPredSubpop}{is the number of clusters in the target mixedpop 
\code{row_cluster <-length(unique(target_cluster))}}

\item{Nodes_group}{string representation of hexidecimal color code for node}

\item{nboots}{is an integer for how many bootstraps are run}
}
\value{
a dataframe containg information for the Lasso prediction results, 
each column
contains prediction results for all subpopulations from each bootstrap run
}
\description{
the training and prediction results from \code{bootstrap}
were written to the object \code{LSOLDA_dat}, the \code{reformat_LASSO}
summarises prediction for \code{n} bootstrap runs
}
\examples{
c_selectID<-1
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
genes <-training_gene_sample
genes <-genes$Merged_unique
LSOLDA_dat <- bootstrap_prediction(nboots = 2, mixedpop1 = mixedpop1, 
    mixedpop2 = mixedpop2, genes=genes, c_selectID, listData =list(),
    cluster_mixedpop1 = colData(mixedpop1)[,1],
    cluster_mixedpop2=colData(mixedpop2)[,1])
reformat_LASSO(LSOLDA_dat=LSOLDA_dat, 
    nPredSubpop=length(unique(colData(mixedpop2)[,1])), c_selectID = 1, 
    mp_selectID =2, nboots = 2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{var_cpp}
\alias{var_cpp}
\title{Calculate variance}
\usage{
var_cpp(x, bias = TRUE)
}
\arguments{
\item{x}{a vector of gene expression.}

\item{bias}{degree of freedom}
}
\value{
a variance value
}
\description{
Calculate variance
}
\examples{
var_cpp(seq_len(10^6))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{plot_reduced}
\alias{plot_reduced}
\title{plot reduced data}
\usage{
plot_reduced(
  reduced_dat,
  color_fac = NULL,
  dims = c(1, 2),
  dimNames = c("Dim1", "Dim2"),
  palletes = NULL,
  legend_title = "Cluster"
)
}
\arguments{
\item{reduced_dat}{is a matrix with genes in rows and cells in columns}

\item{color_fac}{is a vector of colors corresponding to clusters to determine
colors of scattered plots}

\item{dims}{an integer of the number of dimestions}

\item{dimNames}{a vector of the names of the dimensions}

\item{palletes}{can be a customised color pallete that determine colors for 
density plots, if NULL it will use RColorBrewer 
colorRampPalette(RColorBrewer::brewer.pal(sample_num, 'Set1'))(sample_num)}

\item{legend_title}{title of the plot's legend}
}
\value{
a matrix with the top 20 CIDR dimensions
}
\description{
plot PCA, tSNE, and CIDR reduced datasets
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
#CIDR_dim <-CIDR(expression.matrix=assay(mixedpop1))
#p <- plot_reduced(CIDR_dim, color_fac = factor(colData(mixedpop1)[,1]),
#     palletes = seq_len(length(unique(colData(mixedpop1)[,1]))))
#plot(p)
tSNE_dim <-tSNE(expression.mat=assay(mixedpop1))
p2 <- plot_reduced(tSNE_dim, color_fac = factor(colData(mixedpop1)[,1]),
    palletes = seq_len(length(unique(colData(mixedpop1)[,1]))))
plot(p2)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering_bagging.R
\name{clustering_bagging}
\alias{clustering_bagging}
\title{HC clustering for a number of resolutions}
\usage{
clustering_bagging(
  object = NULL,
  ngenes = 1500,
  bagging_run = 20,
  subsample_proportion = 0.8,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  remove_outlier = c(0),
  nRounds = 1,
  PCA = FALSE,
  nPCs = 20,
  log_transform = FALSE
)
}
\arguments{
\item{object}{is a \linkS4class{SingleCellExperiment} object from the train
mixed population.}

\item{ngenes}{number of genes used for clustering calculations.}

\item{bagging_run}{an integer specifying the number of bagging runs to be 
computed.}

\item{subsample_proportion}{a numeric specifying the proportion of the tree 
to be chosen in subsampling.}

\item{windows}{a numeric vector specifying the rages of each window.}

\item{remove_outlier}{a vector containing IDs for clusters to be removed
the default vector contains 0, as 0 is the cluster with singletons.}

\item{nRounds}{a integer specifying the number rounds to attempt to remove 
outliers.}

\item{PCA}{logical specifying if PCA is used before calculating distance 
matrix.}

\item{nPCs}{an integer specifying the number of principal components to use.}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
a list of clustering results containing each bagging run
as well as the clustering of the original tree and the tree itself.
}
\description{
subsamples cells for each bagging run and performs 40 
clustering runs or more depending on windows.
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test <-clustering_bagging(mixedpop2, remove_outlier = c(0),
    bagging_run = 2, subsample_proportion = .7)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_parallel_distance}
\alias{rcpp_parallel_distance}
\title{distance matrix using C++}
\usage{
rcpp_parallel_distance(mat)
}
\arguments{
\item{mat}{an R matrix (expression matrix), rows are genes, columns are cells}
}
\value{
a distance matrix
}
\description{
This function provides fast and memory efficient distance matrix
calculation
}
\examples{
mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=10000)
#library(microbenchmark)
#microbenchmark(rcpp_parallel_distance(mat_test), dist(mat_test), times=3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{distvec}
\alias{distvec}
\title{Compute Distance between two vectors}
\usage{
distvec(x, y)
}
\arguments{
\item{x}{A numeric vector}

\item{y}{A numeric vector}
}
\value{
a numeric distance
}
\description{
Compute Distance between two vectors
}
\examples{
x <-matrix(rnbinom(1000,mu=0.01, size=10),nrow=1000)
x <- x[1,]
y <-matrix(rnbinom(1000,mu=0.01, size=10),nrow=1000)
y <- y[1,]
distvec(x, y)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{top_var}
\alias{top_var}
\title{select top variable genes}
\usage{
top_var(expression.matrix = NULL, ngenes = 1500)
}
\arguments{
\item{expression.matrix}{is a matrix with genes in rows and cells in columns}

\item{ngenes}{number of genes used for clustering calculations.}
}
\value{
a subsetted expression matrix with the top n most variable genes
}
\description{
subset a matrix by top variable genes
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
SortedExprsMat <-top_var(expression.matrix=assay(mixedpop1))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/day_2_cardio_cell_sample.R
\docType{data}
\name{day_2_cardio_cell_sample}
\alias{day_2_cardio_cell_sample}
\title{One of the two example single-cell count matrices to be used
for training \pkg{scGPS} model}
\format{
a list instance, containing a count matrix and a vector with 
clustering information
}
\source{
Dr Joseph Powell's laboratory, IMB, UQ
}
\usage{
day_2_cardio_cell_sample
}
\value{
a list, with the name day2
}
\description{
The count data set contains counts for 16990 genes for 590 cells
randomly subsampled from day-2 cardio-differentiation population.
The vector of clustering information contains corresponding to cells in the 
count matrix
}
\author{
Quan Nguyen, 2017-11-25
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{find_optimal_stability}
\alias{find_optimal_stability}
\title{Find the optimal cluster}
\usage{
find_optimal_stability(
  list_clusters,
  run_RandIdx,
  bagging = FALSE,
  windows = seq(from = 0.025, to = 1, by = 0.025)
)
}
\arguments{
\item{list_clusters}{is a \code{list} object containing 40 clustering results}

\item{run_RandIdx}{is a \code{data frame} object from iterative clustering 
runs}

\item{bagging}{is a logical that is true if bagging is to be performed, 
changes return}

\item{windows}{a numeric vector specifying the ranges of each window.}
}
\value{
bagging == FALSE => a \code{list} with optimal stability, cluster 
count and summary stats bagging == TRUE => a \code{list} with high res 
cluster count, optimal cluster count and keystats
}
\description{
from calculated stability based on Rand indexes for consecutive
clustering run, find the resolution (window), where the stability is the 
highest
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
cluster_all <-clustering(object=mixedpop2)
stab_df <- find_stability(list_clusters=cluster_all$list_clusters,
                         cluster_ref = cluster_all$cluster_ref)
optimal_stab <- find_optimal_stability(list_clusters = 
    cluster_all$list_clusters, stab_df, bagging = FALSE)

}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimReduction.R
\name{PCA}
\alias{PCA}
\title{PCA}
\usage{
PCA(expression.matrix = NULL, ngenes = 1500, scaling = TRUE, npcs = 50)
}
\arguments{
\item{expression.matrix}{An expression matrix, with genes in rows}

\item{ngenes}{number of genes used for clustering calculations.}

\item{scaling}{a logical of whether we want to scale the matrix}

\item{npcs}{an integer specifying the number of principal components to use.}
}
\value{
a list containing PCA results and variance explained
}
\description{
Select top variable genes and perform prcomp
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
t <-PCA(expression.matrix=assay(mixedpop1))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{find_stability}
\alias{find_stability}
\title{Calculate stability index}
\usage{
find_stability(list_clusters = NULL, cluster_ref = NULL)
}
\arguments{
\item{list_clusters}{is a object from the iterative clustering runs}

\item{cluster_ref}{is a object from the reference cluster}
}
\value{
a \code{data frame} with stability scores and rand_index results
}
\description{
from clustering results, compare similarity between clusters by
adjusted Randindex
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
cluster_all <-clustering(object=mixedpop2)
stab_df <- find_stability(list_clusters=cluster_all$list_clusters,
                         cluster_ref = cluster_all$cluster_ref)
}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DimReduction.R
\name{tSNE}
\alias{tSNE}
\title{tSNE}
\usage{
tSNE(
  expression.mat = NULL,
  topgenes = 1500,
  scale = TRUE,
  thet = 0.5,
  perp = 30
)
}
\arguments{
\item{expression.mat}{An expression matrix, with genes in rows}

\item{topgenes}{number of genes used for clustering calculations.}

\item{scale}{a logical of whether we want to scale the matrix}

\item{thet}{numeric; Speed/accuracy trade-off (increase for less accuracy)}

\item{perp}{numeric; Perplexity parameter
(should not be bigger than 3 * perplexity < nrow(X) - 1, see details for 
interpretation)}
}
\value{
a tSNE reduced matrix containing three tSNE dimensions
}
\description{
calculate tSNE from top variable genes
}
\examples{
day2 <- day_2_cardio_cell_sample
mixedpop1 <-new_scGPS_object(ExpressionMatrix = day2$dat2_counts, 
    GeneMetadata = day2$dat2geneInfo, CellMetadata = day2$dat2_clusters)
t <-tSNE(expression.mat = assay(mixedpop1))

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{plot_optimal_CORE}
\alias{plot_optimal_CORE}
\title{plot one single tree with the optimal clustering result}
\usage{
plot_optimal_CORE(
  original_tree,
  optimal_cluster = NULL,
  shift = -100,
  values = NULL
)
}
\arguments{
\item{original_tree}{a dendrogram object}

\item{optimal_cluster}{a vector of cluster IDs for cells in the dendrogram}

\item{shift}{a numer specifying the gap between the dendrogram and the 
colored}

\item{values}{a vector containing color values of the branches and the
colored bar underneath the tree bar underneath the dendrogram. This
parameter allows better selection of colors for the display.}
}
\value{
a plot with colored braches and bars for the optimal clustering 
result
}
\description{
after an optimal cluster has been identified, users may use this
function to plot the resulting dendrogram with the branch colors represent 
clutering results
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
CORE_cluster <- CORE_clustering(mixedpop2, remove_outlier = c(0))
key_height <- CORE_cluster$optimalClust$KeyStats$Height
optimal_res <- CORE_cluster$optimalClust$OptimalRes
optimal_index = which(key_height == optimal_res)
plot_optimal_CORE(original_tree= CORE_cluster$tree, 
    optimal_cluster = unlist(CORE_cluster$Cluster[optimal_index]),
    shift = -2000)

}
\author{
Quan Nguyen, 2017-11-25
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{tp_cpp}
\alias{tp_cpp}
\title{Transpose a matrix}
\usage{
tp_cpp(X)
}
\arguments{
\item{X}{an R matrix (expression matrix)}
}
\value{
a transposed matrix
}
\description{
Transpose a matrix
}
\examples{
mat_test <-matrix(rnbinom(1000000,mu=0.01, size=10),nrow=100)
tp_mat <- tp_cpp(mat_test)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{rand_index}
\alias{rand_index}
\title{Calculate rand index}
\usage{
rand_index(tab, adjust = TRUE)
}
\arguments{
\item{tab}{a table containing different clustering results in rows}

\item{adjust}{a logical of whether to use the adjusted rand index}
}
\value{
a rand_index value
}
\description{
Comparing clustering results Function for calculating randindex
(adapted from the function by Steve Horvath and Luohua Jiang, UCLA, 2003)
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts,
GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
cluster_all <-clustering(object=mixedpop2)

rand_index(table(unlist(cluster_all$list_clusters[[1]]), 
cluster_all$cluster_ref))

}
\author{
Quan Nguyen and Michael Thompson, 2018-05-11
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{plot_CORE}
\alias{plot_CORE}
\title{Plot dendrogram tree for CORE result}
\usage{
plot_CORE(original.tree, list_clusters = NULL, color_branch = NULL)
}
\arguments{
\item{original.tree}{the original dendrogram before clustering}

\item{list_clusters}{a list containing clustering results for each of the}

\item{color_branch}{is a vector containing user-specified colors (the number
of unique colors should be equal or larger than the number of clusters). This
parameter allows better selection of colors for the display.}
}
\value{
a plot with clustering bars underneath the tree
}
\description{
This function plots CORE and all clustering results underneath
}
\examples{
day5 <- day_5_cardio_cell_sample
cellnames <- colnames(day5$dat5_counts)
cluster <-day5$dat5_clusters
cellnames <-data.frame('Cluster'=cluster, 'cellBarcodes' = cellnames)
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts,
    GeneMetadata = day5$dat5geneInfo, CellMetadata = cellnames)
CORE_cluster <- CORE_clustering(mixedpop2, remove_outlier = c(0))
plot_CORE(CORE_cluster$tree, CORE_cluster$Cluster)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_Eucl_distance_NotPar}
\alias{rcpp_Eucl_distance_NotPar}
\title{Function to calculate Eucledean distance matrix without parallelisation}
\usage{
rcpp_Eucl_distance_NotPar(mat)
}
\arguments{
\item{mat}{an R matrix (expression matrix), with cells in rows and genes
in columns}
}
\value{
a distance matrix
}
\description{
Function to calculate Eucledean distance matrix without parallelisation
}
\examples{
mat_test <-matrix(rnbinom(100000,mu=0.01, size=10),nrow=1000)
#library(microbenchmark)
#microbenchmark(rcpp_Eucl_distance_NotPar(mat_test), dist(mat_test), times=3)

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{add_import}
\alias{add_import}
\title{add_import}
\description{
import packages to namespace
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CORE_clustering.R
\name{clustering}
\alias{clustering}
\title{HC clustering for a number of resolutions}
\usage{
clustering(
  object = NULL,
  ngenes = 1500,
  windows = seq(from = 0.025, to = 1, by = 0.025),
  remove_outlier = c(0),
  nRounds = 1,
  PCA = FALSE,
  nPCs = 20,
  verbose = FALSE,
  log_transform = FALSE
)
}
\arguments{
\item{object}{is a \linkS4class{SingleCellExperiment} object from the
train mixed population}

\item{ngenes}{number of top variable genes to be used}

\item{windows}{a numeric specifying the number of windows to test}

\item{remove_outlier}{a vector containing IDs for clusters to be removed
the default vector contains 0, as 0 is the cluster with singletons}

\item{nRounds}{number of iterations to remove a selected clusters}

\item{PCA}{logical specifying if PCA is used before calculating distance 
matrix}

\item{nPCs}{number of principal components from PCA dimensional reduction to
be used}

\item{verbose}{a logical whether to display additional messages}

\item{log_transform}{boolean whether log transform should be computed}
}
\value{
clustering results
}
\description{
performs 40 clustering runs or more depending on windows
}
\examples{
day5 <- day_5_cardio_cell_sample
mixedpop2 <-new_summarized_scGPS_object(ExpressionMatrix = day5$dat5_counts, 
    GeneMetadata = day5$dat5geneInfo, CellMetadata = day5$dat5_clusters)
test <-clustering(mixedpop2, remove_outlier = c(0))
}
\author{
Quan Nguyen, 2017-11-25
}
