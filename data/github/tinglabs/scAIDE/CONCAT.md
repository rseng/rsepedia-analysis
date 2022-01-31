scAIDE
======================
> scAIDE is an unsupervised clustering framework for single-cell RNA-seq data. We revealed both putative and rare cell types in the 1.3 million neural cell dataset by global clustering. We obtained 64 clusters within 30 minutes, which were then mapped to 19 putative cell types. Three different subpopulations of neural stem/pregenitor cells were identified, each primed for different developmental stages. 

## Overview
![Overview](figures/Overview.png)

## Table of content
- [Introduction](#Introduction)
- [Installation](#Installation)
    - [Required Packages](#--required-installations)
    - [Install](#--Install)
- [Example Usage](#example-usage)
    - [Preprocessing](#--Preprocessing)
    - [AIDE](#--AIDE)
    - [RPH-kmeans](#--RPH-kmeans)
    - [Biological Analysis](#--biological-analysis)
    - [Example Results](#--examples-results)
    - [Scalability](#--scalability)
- [Citation](#citation-&-references)
- [Maintenance](#Maintenance)

## Introduction
There are three main parts to this clustering framework: AIDE embedding, clustering (RPH-kmeans), and biological analysis.

- **AIDE**: autoencoder-imputed distance-preserved embedding (a novel deep learning architecture that learns a good representation of single-cell data accounting for biological noise)
- **RPH-kmeans**: Random projection hashing based k-means algorithm (a novel algorithm which improves the detection of rare cell types)
- **biological analysis**: Biological analytics code are packed into an R package, scAIDE.


## Installation
### - Required Installations
AIDE and RPH-kmeans are implemented in python, biological analytics code in R.

- Python3: tensorflow 1.14, numpy, scipy, scikit_learn, tqdm, seaborn
- R: parallel, gmp, ggplot2

### - Install
- AIDE: please refer to [https://github.com/tinglabs/aide](https://github.com/tinglabs/aide) for details.
- RPH-kmeans: please refer to [https://github.com/tinglabs/rph_kmeans](https://github.com/tinglabs/rph_kmeans) for details.
- scAIDE (R package):

	```r
	require(devtools)
	setwd("where scAIDE folder is located")
	install("scAIDE")
	```

## Example Usage:
A demo is provided in [demo](https://github.com/tinglabs/scAIDE/tree/master/demo) folder, showing details of data preprocessing, embedding with AIDE and clustering with RPH-kmeans. Two datasets (`Mouse bladder` and `Mouse retina`) are also given.

### - Preprocessing
A pre-processed single-cell data is accepted, provided that it is normalized and log-transformed (for optimal performance). 
In Python, the input is configured as n cells (rows) by m genes (columns).


### - AIDE
```python
# Load data:
# For small to medium size datasets (up to few hundred thousands of cells)
import pandas as pd
import numpy as np

# Make sure that the final input is n cells (rows) by m genes (cols)
sc_data = pd.read_csv("single_cell_dataset.csv", index_col=0)
sc_data = sc_data.values.astype('float32') # type = np.ndarray

# Configuring AIDE parameters:
from aide import AIDE, AIDEConfig
config = AIDEConfig()
# We may tune the following 4 parameters, but default values are usually sufficient.
config.pretrain_step_num = 1000 # Pretrain step
config.ae_drop_out_rate = 0.4 # Dropout rate
config.alpha = 12.0 # A parameter that determines the portion of AE loss vs MDS encoder loss
config.early_stop = True # Early stop (maximum step number = 20000, minimum step number = 4000)

# Running AIDE:
encoder = AIDE(name = "sc_test", save_folder = "sc_test_folder")
sc_embedding = encoder.fit_transform(sc_data, config=config)

# save embedding
np.savetxt("~/sc_embedding.txt", sc_embedding)
```

### - RPH-kmeans

```python
from rph_kmeans import RPHKMeans
# In the case that n_clusters is already known:
clt = RPHKMeans(n_init=10, n_clusters=10)
clt_labels = clt.fit_predict(sc_embedding)

# In the case that n_clusters is unknown: In order to automatically detect the number of clusters, 
# we implemented a weighted BIC value that determines the optimal k based on 'kneedle' point.

# Important Note: Please set the parameter max_point to a smaller number for small datasets (i.e. less than 2000 cells). 
# The max_point defaults to 2000 cells, and the RPH algorithm stops when the reduced number of cells is below max_point.
# In other words, RPH is not performed if the dataset is smaller than max_point.

max_point = 50 # Defaults to 2000

from rph_kmeans import select_k_with_bic
kmax = 30 # set the maximum number of k to explore
optimal_k, _, _ = select_k_with_bic(sc_embedding, kmax=kmax, point_reducer_kwargs={'max_point':max_point})

clt = RPHKmeans(n_init=10, n_clusters=optimal_k, max_point = max_point) # run RPH-kmeans with optimal_k to get the clustering results

clt_labels = clt.fit_predict(sc_embedding)

# Output results
np.savetxt("~/clt_labels.txt", clt_labels)
```

### - Biological Analysis

```r
library(scAIDE)

# load original data and clustering results:
sc_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # rows = genes, cols = cells
sc_clusters <- read.table("clt_labels.txt")$V1
# Evaluate wilcox rank sum test and log fold change values
eval_gene_markers <- store_markers(gene_expression_data, sc_clusters, threads = 8)
gene_names <- rownames(gene_expression_data)
# returns the list of markers for each cluster, with your specified threshold
cluster_markers_list <- curate_markers(eval_gene_markers, gene_names, wilcox_threshold=0.001, logfc_threshold=1.5)

# Cell type assignment probability according to the markers in the database
# panglao_marker_list: pre-processed list of markers for neural and immune cell types.
# returns a cluster (rows) by cell types (cols) matrix
celltype_prob <- calculate_celltype_prob(cluster_markers_list, panglao_marker_list, type = "jacc")
celltype_id <- rowMaxs(celltype_prob)

# Enrichment probability (based on hypergeometric distribution), this is to be compared with celltype_id to ensure that the number of marker genes detected is of statistical significance.
n_genes <- nrow(gene_expression_data)
# returns a cluster (rows) by cell types (cols) matrix, with p-value entries
enrichment_prob <- calculate_enrichment_prob(cluster_markers_list, panglao_marker_list, n_genes, type = "jacc")

######################################################################
# Visualizing marker genes:
# example marker list:
selected_marker_genes <- c("SOX2", "ALDOC", "CCND2", "OLIG1", "OLIG2")
gene_expression_subset <- gene_expression_data[match(tolower(selected_marker_genes), tolower(rownmaes(gene_expression_data))), ]
# Process the data for plots
processed_markers <- process_marker_expression(gene_expression_subset, sc_clusters)
# Specify the cell type order, and the gene order that you wish to plot
cell_levels <- unique(sc_clusters)
gene_levels <- selected_marker_genes
marker_plot <- plot_marker_expression(processed_markers, gene_levels=gene_levels, cell_levels=cell_levels)
```

### - Example results
The annotated labels for PBMC and Neural datasets are included in the folder 'Predicted labels'. The .RData files include the predicted annotated labels for these datasets. </br>
</br>
The following figures show the results for the PBMC 68k dataset and the 1.3 million neural dataset. 

<p align="center">
  <img src=figures/pbmc.png alt="pbmc" title="pbmc" align="center" height="300">
  <img src= figures/neural.png alt="neural" title="neural" align="center" height="300">
</p>

### - Scalability
The time taken to cluster 1.3 million cells (with roughly 20,000 genes) is less than 30 minutes, using 7GB of memory.

## File Illustration
- `sc_cluster`: Codes and results of clustering experiments using AIDE and RPH-kmeans.
- `baseline`: Codes and results of clustering experiments using baseline tools (e.g. DCA, MAGIC, scScope, scDeepCluster, ...).
- `demo`: A demo of data preprocessing, embedding with AIDE and clustering with RPH-kmeans.
- `scAIDE`: The R package of biological analysis.
- `figures`: Figures of README

## Citation & References

scAIDE: clustering of large-scale single-cell RNA-seq data reveals putative and rare cell types. NAR Genomics and Bioinformatics 2020.

References:
###### 1. Zheng, G. X. et al. (2017) Massively parallel digital transcriptional profiling of single cells. Nature Communications 8, 14049, doi:10.1038/ncomms14049
###### 2. Genomics, X. J. C. B. 1.3 million brain cells from E18 mice. (2017).
###### 3. Franzen, O., Gan, L. M. & Bjorkegren, J. L. M. PanglaoDB: a web server for exploration of mouse and human single-cell RNA sequencing data. Database (Oxford) 2019.

## Maintenance

If there's any questions / problems regarding scAIDE, please feel free to contact Ken Xie - kk.xie419@gmail.com and Huang Yu - yuhuang-cst@foxmail.com. Thank you!

# sc_cluster
Single cell clustering experiments that use **AIDE** to reduce the dimension of cell-gene matrix and then use **RPH-KMeans** to do clustering.

- [AIDE](https://github.com/tinglabs/aide): Autoencoder-imputed distance-preserved embedding, which combines both **Multidimentional Scaling (MDS)** and **AutoEncoder (AE)** technique, aiming to preserve the distance between imputed data generated by **AE** when reducing dimension.
- [rph-kmeans](https://github.com/tinglabs/rph_kmeans): a variant of kmeans algorithm in which the initial centers is produced by point reduction process using **random projection (RP)**, which is one of the **local sensitive hashing (LSH)** technique.

## File Illustration
```
result
	dim_origin_clt_eval_paper: output of dim_origin_clt_eval.py
	dim_reduce_eval_paper: output of dim_reduce_clt_eval.py
	data_explain_paper: basic infomation of data
	embed_mem_use_paper: memory of running AIDE and PCA
src
	reader.py: read data from raw file
	preprocess.py: data preprocess
	dim_origin_clt_eval.py: data preprocess + clustering
	dim_reduce_clt_eval.py: data preprocess + dimension reduction + clustering
data (not uploaded)
	raw
		1M_neurons
			1M_neurons_filtered_gene_bc_matrices_h5.h5
		sc_brain
			l5_all.loom
		scDeepCluster
			PBMC_68k
				PBMC_68k.h5
			Shekhar_mouse_retina
				Shekhar_mouse_retina.h5
			10X_PBMC
				10X_PBMC.h5
			mouse_bladder_cell
				mouse_bladder_cell.h5
			mouse_ES_cell
				mouse_ES_cell.h5
			worm_neuron_cell
				worm_neuron_cell.h5
``` 

# Demo
- This is a demo for analyzing `Mouse bladder` (small data; 2746 cells *	20670 genes) and `Mouse retina` (large data; 27499 cells * 13166 genes), including the recommended data preprocessing, embedding with AIDE and clustering with RPH-kmeans.
- The two dataset are located in `data` folder.
- The running results are located in `aide_for_bladder` folder and `aide_for_retina` folder respectively.
- For both of the AIDE and RPH-kmeans, default hyper-parameters are used in this demo.

## For small data
The code of analyzing the `Mouse bladder` dataset are provided in `small.py`. Simply run the script (`python3 small.py`) and the running log is as follows:

```
Pretrain begin============================================
Step 50(5.0%): Batch Loss=494.993896484375
Step 100(10.0%): Batch Loss=473.5043029785156
...
Step 950(95.0%): Batch Loss=432.5137634277344
Step 1000(100.0%): Batch Loss=412.89599609375
Pretrain end.============================================
Step 50(0.25%); Global Step 50: Batch Loss=584.4363403320312; [Reconstruct, MDS, L2] Loss = [523.41284, 5.085291, 0.0]
...
Step 5100(25.5%); Global Step 5100: Validation Loss=468.5406188964844; [Reconstruct, MDS, L2] Loss = [451.28278, 1.4381549, 0.0]; Min Val Loss = 467.6856384277344; No Improve = 5; 
Step 5150(25.75%); Global Step 5150: Batch Loss=468.8702392578125; [Reconstruct, MDS, L2] Loss = [450.14258, 1.5606391, 0.0]
Step 5200(26.0%); Global Step 5200: Batch Loss=459.5696716308594; [Reconstruct, MDS, L2] Loss = [443.746, 1.3186402, 0.0]
No improve = 6, early stop!
Training end. Total step = 5200
Type of embedding = <class 'numpy.ndarray'>; Shape of embedding = (2746, 256); Data type of embedding = float32
RPH-KMeans (n_init = 1): ARI = 0.6105 (0.0279), NMI = 0.7604 (0.0086)
RPH-KMeans (n_init = 10): ARI = 0.6679 (0.0630), NMI = 0.7754 (0.0135)
KMeans (init = k-means++, n_init = 1): ARI = 0.5785 (0.0377), NMI = 0.7644 (0.0126)
KMeans (init = k-means++, n_init = 10): ARI = 0.5821 (0.0302), NMI = 0.7648 (0.0070)
```

Here shows the history loss of AIDE:

![history loss](aide_for_bladder/loss.png)

## For large data (e.g. cell_num > 100000)
The code of analyzing the `Mouse retina` dataset are provided in `large.py`. Simply run the script (`python3 large.py `) and the running log is as follows:

```
Pretrain begin============================================
Step 50(5.0%): Batch Loss=458.53656005859375
Step 100(10.0%): Batch Loss=443.490234375
...
Step 950(95.0%): Batch Loss=420.61492919921875
Step 1000(100.0%): Batch Loss=427.2901611328125
Pretrain end.============================================
Step 50(0.25%); Global Step 50: Batch Loss=469.0723876953125; [Reconstruct, MDS, L2] Loss = [445.16647, 1.9921587, 0.0]
...
Step 5700(28.5%); Global Step 5700: Validation Loss=436.14697265625; [Reconstruct, MDS, L2] Loss = [429.66595, 0.54008174, 0.0]; Min Val Loss = 435.47027587890625; No Improve = 5; 
Step 5750(28.75%); Global Step 5750: Batch Loss=434.73870849609375; [Reconstruct, MDS, L2] Loss = [428.11176, 0.5522464, 0.0]
Step 5800(29.0%); Global Step 5800: Batch Loss=439.7626037597656; [Reconstruct, MDS, L2] Loss = [432.93768, 0.56874436, 0.0]
No improve = 6, early stop!
Training end. Total step = 5800
Type of embedding = <class 'numpy.ndarray'>; Shape of embedding = (27499, 256); Data type of embedding = float32
RPH-KMeans (n_init = 1): ARI = 0.8914 (0.0306), NMI = 0.8248 (0.0117)
RPH-KMeans (n_init = 10): ARI = 0.8859 (0.0155), NMI = 0.8246 (0.0054)
KMeans (init = k-means++, n_init = 1): ARI = 0.7556 (0.1552), NMI = 0.7944 (0.0228)
KMeans (init = k-means++, n_init = 10): ARI = 0.6659 (0.1138), NMI = 0.7895 (0.0209)
```

Here shows the history loss of AIDE:

![history loss](aide_for_retina/loss.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{process_panglao_db}
\alias{process_panglao_db}
\title{Pre-processes the Panglao database to curate markers for desired cell types.}
\usage{
process_panglao_db(cell_type_name_list, database)
}
\arguments{
\item{cell_type_name_list}{This is the corresponding cell type names in the database.}

\item{database}{The database table which is loaded.}
}
\description{
This function processes the panglao database to the format of lists, where each list is a cell type containing marker genes.
}
\section{Biological Analysis}{

}

\examples{
panglao_table <-read.csv("panglao_database.csv", header = T, sep = "\t")
type_names <- c("astrocytes", "oligodendrocytes")
preprocessed_pl <- process_panglao_db(type_names, panglao_table)

# Preloaded Panglao database objects:
# panglao_db: the full panglao marker list as of date 7th Feb 2020
# panglao_marker_list: the pre-processed list for all neural and immune cell types.
}
\keyword{calculate_enrichment_prob}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{process_marker_expression}
\alias{process_marker_expression}
\title{Prepares the data for plotting marker genes}
\usage{
process_marker_expression(gene_cell_matrix, cell_labels)
}
\arguments{
\item{gene_cell_matrix}{This is the pre-processed gene expression matrix (normalized and logged). The subset expression from the markers list is used}

\item{cell_labels}{Cell labels which can be annotated or just cluster labels.}
}
\description{
This function calculates the cluster-specific average gene expressions and the percentage of zeros in each cluster. Prepares the data for plotting
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1

# example marker list:
selected_marker_genes <- c("SOX2", "ALDOC", "CCND2", "OLIG1", "OLIG2")
gene_expression_subset <- gene_expression_data[match(tolower(selected_marker_genes), tolower(rownmaes(gene_expression_data))), ]
processed_markers <- process_marker_expression(gene_expression_subset, cluster_vector)
}
\keyword{process_marker_expression}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{curate_markers}
\alias{curate_markers}
\title{Find the expression markers specific to your cluster assignments}
\usage{
curate_markers(
  whole_list,
  gene_names,
  wilcox_threshold = 0.001,
  logfc_threshold = 1.5
)
}
\arguments{
\item{whole_list}{This is the result returned from the function 'store_markers', containing the wilcox test and log fold change values.}

\item{gene_names}{A character vector containing the gene names in the same order as your input gene expression matrix. e.g. c("SOX2", "ALDOC")}

\item{wilcox_threshold}{This is the p-value cut-off for the wilcox test, it is set to 0.001 by default. Recommend a smaller value for a stricter detection}

\item{logfc_threshold}{This is the log fold change threshold, which is 1.5 by default. This can be increased for stricter restriction on sample sizes.}
}
\description{
This function calculates the wilcox rank sum test and log fold change ratio for all genes in each cluster.
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1
eval_gene_markers <- store_markers(gene_expression_data, cluster_vector, threads = 8)
gene_names = rownames(gene_expression_data)
cluster_markers_list <- curate_markers(eval_gene_markers, gene_names, wilcox_threshold=0.001, logfc_threshold=1.5)
}
\keyword{curate_markers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{store_markers}
\alias{store_markers}
\title{Calculate the wilcox rank sum test and log fold change values}
\usage{
store_markers(input_data_matrix, labels_vector, threads = 8)
}
\arguments{
\item{input_data_matrix}{A pre-processed gene expression matrix (ie normalization and log). Matrix input is a genes (row) by cells (col) gene expression matrix with all genes.}

\item{labels_vector}{A numeric vector containing the cell cluster assignment of each cell. If minimum value is 0, it will automatically add 1 to the vector.}

\item{threads}{The number of threads used for parallel computing. Defaults to 8}
}
\description{
This function calculates the wilcox rank sum test and log fold change ratio for all genes in each cluster.
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1
eval_gene_markers <- store_markers(gene_expression_data, cluster_vector, threads = 8)
}
\keyword{store_markers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{calculate_celltype_prob}
\alias{calculate_celltype_prob}
\title{Calculate the cell type probablity assignment according to a markers database}
\usage{
calculate_celltype_prob(clt_marker_list, marker_database_list, type = "jacc")
}
\arguments{
\item{clt_marker_list}{This is the result returned from the function 'curate_markers', a list object containing marker genes for each specific cluster.}

\item{type}{This determines the type of overlap to calculate. Defaults to the Jaccard index, "jacc". Accuracy ("ac") and F1 ("f1") are also available. We recommend using jaccard or accuracy in applications.}

\item{marker_database}{A list object that contains the marker genes from each specific cell type. A pre-processed variable is stored as 'panglao_marker_list', which contains markers for neural and immune cell types.}
}
\description{
This function calculates probability of each cell type in the desired database, according to the number of overlapping genes. We have included the "Panglao" database as of date 7th Feb 2020 in our package.
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1
eval_gene_markers <- store_markers(gene_expression_data, cluster_vector, threads = 8)
gene_names = rownames(gene_expression_data)
cluster_markers_list <- curate_markers(eval_gene_markers, gene_names, wilcox_threshold=0.001, logfc_threshold=1.5)

celltype_prob <- calculate_celltype_prob(cluster_markers_list, panglao_marker_list, type = "jacc")
celltype_id <- rowMaxs(celltype_prob)
}
\keyword{calculate_celltype_prob}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{find_specific_marker}
\alias{find_specific_marker}
\title{Utility function for calculating the overlaps of specific markers vs each cluster}
\usage{
find_specific_marker(gene_name, f_list, type = "jacc")
}
\arguments{
\item{gene_name}{List of gene names that we want to query in the clusters.}

\item{f_list}{A list object that contains the marker genes from each specific cluster.}

\item{type}{This determines the type of overlap to calculate. Defaults to the Jaccard index, "jacc". Accuracy ("ac") and F1 ("f1") are also available. We recommend using jaccard or accuracy in applications.}
}
\description{
This is the utility function behind calculate_celltype_prob, which implements the Jaccard, accuracy and F1 to calculate the value of overlaps between marker genes and cluster specific genes.
}
\section{Biological Analysis}{

}

\keyword{find_specific_marker}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_marker_expressions}
\alias{plot_marker_expressions}
\title{Plot for visualizing cluster-specific markers}
\usage{
plot_marker_expressions(processed_list, gene_levels = NULL, cell_levels = NULL)
}
\arguments{
\item{processed_list}{This is the pre-processed list from the function 'process_marker_expression', which returns an average expression matrix and the percentage matrix (ratio of zeros)}

\item{gene_levels}{This is a character vector containing gene names in the order that you want to display in the plot}

\item{cell_levels}{This is a character vector containing cell type annotations / cluster labels in the order that you want to display in the plot}
}
\description{
This function plots the average expression values (denoted by color of dots) and the percentage of zeros in the cluster (denoted size of dots)
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1

# example marker list:
selected_marker_genes <- c("SOX2", "ALDOC", "CCND2", "OLIG1", "OLIG2")
gene_expression_subset <- gene_expression_data[match(tolower(selected_marker_genes), tolower(rownmaes(gene_expression_data))), ]
processed_markers <- process_marker_expression(gene_expression_subset, cluster_vector)
cell_levels <- unique(cluster_vector)
gene_levels <- selected_marker_genes
marker_plot <- plot_marker_expression(processed_markers, gene_levels=gene_levels, cell_levels=cell_levels)
}
\keyword{plot_marker_expressions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analytical_functions.R
\name{calculate_enrichment_prob}
\alias{calculate_enrichment_prob}
\title{Calculate the enrichment probablity of each cell type based on a hypergeometric distribution.}
\usage{
calculate_enrichment_prob(
  clt_marker_list,
  marker_databse_list,
  total_background_genes = 10000,
  type = "jacc"
)
}
\arguments{
\item{clt_marker_list}{This is the result returned from the function 'curate_markers', a list object containing marker genes for each specific cluster.}

\item{total_background_genes}{This is the number of genes in your dataset}

\item{type}{This determines the type of overlap to calculate. Defaults to the Jaccard index, "jacc". Accuracy ("ac") and F1 ("f1") are also available. We recommend using jaccard or accuracy in applications.}

\item{marker_database}{A list object that contains the marker genes from each specific cell type. A pre-processed variable is stored as 'panglao_marker_list', which contains markers for neural and immune cell types.}
}
\description{
This function calculates enrichment probability for each cell type vs each cluster. In principle, it determines the statistical validity of the cell type assignment, based on the number of overlapping genes.
}
\section{Biological Analysis}{

}

\examples{
gene_expression_data <- read.csv("single_cell_dataset.csv", header = T, row.names = 1) # make sure that genes are in rows
cluster_vector <- read.table("results_from_rpkmeans.txt")$V1
eval_gene_markers <- store_markers(gene_expression_data, cluster_vector, threads = 8)
gene_names = rownames(gene_expression_data)
cluster_markers_list <- curate_markers(eval_gene_markers, gene_names, wilcox_threshold=0.001, logfc_threshold=1.5)

n_genes <- nrow(gene_expression_data)
enrichment_prob <- calculate_enrichment_prob(cluster_markers_list, panglao_marker_list, n_genes, type = "jacc")
}
\keyword{calculate_enrichment_prob}
