# DrivAER for manifold interpretation in scRNA-seq data
**DrivAER** is a method for identification of **Driv**ing transcriptional programs based on **A**uto**E**ncoder derived **R**elevance scores. 
DrivAER infers relevance scores for transcriptional programs with respect to specified outcomes of interest in single-cell RNA sequencing data, such as psuedotemporal ordering or disease status.

See our [manuscript](https://academic.oup.com/gigascience/article/9/12/giaa122/6029835) for more details.

<p align="center"> 
<img src="Figure1.PNG">
</p>

**Workflow** (a) DrivAER iteratively subjects annotated gene sets to unsupervised dimension reduction via DCA. (b) For each gene set the generated two-dimensional data manifold coordinates are used as (c) input features in a random forest model to predict the outcome of interest (i.e. pseudotemporal ordering). (d) The random forest prediction accuracy represents the relevance score. 

## Reproducibility
To reproduce the analysis and figures presented in our manuscript please see the [*Reproducibility*](https://github.com/lkmklsmn/DrivAER/tree/master/Reproducibility) folder.

## Check out our live, interactive tutorial!
The following [Google colab](https://colab.research.google.com/) notebooks allow you to interactively explore DrivAER and can be run within your browser. We have prepared two analysis examples:
1. [Blood development](https://colab.research.google.com/drive/1zrQ7l3Orz7h-eGEX7MHRIBTTXzL_vu9O#scrollTo=VzAzfdHZrOWz)
2. [Interferon stimulation](https://colab.research.google.com/drive/13DA_dYlRjlKma1d9VB65JrfPhBvkEGDC#scrollTo=roa2rIBT1s_R)

## Installation
### via pip
	pip install git+https://github.com/lkmklsmn/DrivAER.git
### via git
	git clone https://github.com/lkmklsmn/DrivAER
	cd DrivAER
	python setup.py install

## Input
1. Raw count expression matrix
2. Outcome of interest (pseudotemporal ordering/cell grouping etc)
3. Gene set annotation

## Output
1. Relevance scores for each annotated transcriptional program
2. Data manifolds derived from each transcriptional program
3. Various visualizations (heatmap, DCA embedding, barplots)

## Usage

### Step 1: Get Gene Set Annotations
DrivAER supports annotations in gmt and csv format, as well as user-defined annotation file.
#### 1.1 From gmt format
The gmt format files for Broad's MSigDB can be downloaded from the [Broad Website](https://www.gsea-msigdb.org/gsea/downloads.jsp).
| Gene set | Source | Gene1 | Gene2 | Gene3|
| ---------- | ---------- |  :----:  |  :----:  |  :----:  | 
| set1 | source | gene1 | gene2 | gene3 |
| set2 | source | gene1 | gene2 | gene3 |
| set3 | source | gene1 | gene2 | gene3 |
	import DrivAER as dv
	C3_mouse = dv.get_anno(filename = "C3.gmt", filetype = "gmt", conv_mouse = True)
#### 1.2 From tsv format
The tsv format files can be downloaded from Trandcription Factor sites, such as [TRRUST](https://www.grnpedia.org/trrust/downloadnetwork.php). DrivAER provides built-in annotations from TRRUST.
| Transcription factor | Target | Type | Source|
| ---------- | ---------- |  :----:  |  :----:  | 
| set1 | gene1 | XX | XX |
| set1 | gene2 | XX | XX |
| set1 | gene3 | XX | XX |
| set2 | gene1 | XX | XX |
	trrust_human = dv.get_anno(filename = "trrust_human.tsv", filetype = "tsv", conv_mouse = False)
#### 1.3 From user-defined gene set annotations
Users can create your own gene set annotations. The format is a pandas series. Index are trandcription factor names or gene set names. Each row contains a list of corresponding genes.

### Step 2: Calculate relevance scores
	res = dv.calc_relevance(count = your_count, pheno = your_pt, tf_targets = C3_mouse, min_targets=5,
                   ae_type="nb-conddisp", epochs=100, early_stop=3, hidden_size=(8, 2, 8), verbose=False)

Additionally, users can replace the DCA with other dimension reduction methods. The commands *calc_relevance_pca*, *calc_relevance_tsne*, *calc_relevance_umap* will perform dimension reduction based on PCA, tSNE and UMAP, respectively.

### Step 3: Generate visualizations
	dv.rank_plot(result, save, path)
	dv.embedding_plot(result, tf_name, pheno, datatype, save, path)
	dv.gene_plot(result, count, tf_name, gene, save, path)
