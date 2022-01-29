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
