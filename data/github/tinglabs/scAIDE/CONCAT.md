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
