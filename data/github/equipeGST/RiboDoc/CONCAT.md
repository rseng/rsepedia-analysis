# Step by step tutorial

Welcome to the RiboDoc tool tutorial !  

>RiboDoc is designed to perform all classical steps of **ribosome profiling** (RiboSeq) data analysis from the FastQ files to the differential expression analysis.


## 1) Install Docker  
First of all, Docker must be present in version 19 or higher.  
If you don’t already have it, now is the time to fix it !  
Docker Engine is available on different OS like macOS and Windows 10 through Docker Desktop and as a static binary installation for a variety of Linux platforms. All are available here : https://docs.docker.com/engine/install/   

>Tips:  
>&emsp;&emsp;&emsp;For Windows, WSL2 and Ubuntu from Microsoft store applications are needed too.  

## 2) Directory preparation  
RiboDoc does not need installation (yipee) but a precise folder architecture is required (boo).  
The first step is the project folder creation. It is named as your project and will be the volume linked to Docker.  
Then, two sub-folders and a file have to be created and completed respectively.   

> **Caution, those steps are majors for the good course of the analysis**  
> **Subfolders don’t have uppercase**    

Folder architecture at this step:  
Project_name  

### a) *fastq* subfolder  
This subfolder, as its name suggests, should contain your FastQs. These must be compressed in .gz.  
**Format of file name must be as following:**  
&emsp;&emsp;&emsp;biological_condition_name.replicat.fastq.gz     
For example, for the first replicate of the wild-type condition the sample can be named *Wild_Type.1.fastq.gz* and *Wild_Type.2.fastq.gz* for the second replicate.   

>Caution, for **Windows**, extensions can be hidden.    

Folder architecture at this step:  
Project_name  
└── fastq   
&emsp;&emsp;&emsp;├── Wild_Type.1.fastq.gz   
&emsp;&emsp;&emsp;├── Wild_Type.2.fastq.gz   
&emsp;&emsp;&emsp;├── Mutant.1.fastq.gz   
&emsp;&emsp;&emsp;└── Mutant.2.fastq.gz   

If you want to try RiboDoc on an example dataset, you can find our data with on GEO : [GEO GSE173856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173856)

### b) *database* subfolder  
In this subfolder, you must put at least the following three files:  
- Your genome fasta file: Whether it's the genome or the transcriptome, it must be your reference fasta file where reads will be aligned. It must, like the other files, be downloaded from the [Ensembl](https://www.ensembl.org/index.html) database.  
For example, for an entire human genome, you can look for [Human genome](https://www.ensembl.org/Homo_sapiens/Info/Index) and download the *Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz* file in the genome assembly list of fasta files, unzip it and place it in the *database* subfolder of your project directory.  

> Note:  
>&emsp;&emsp; If the transcriptome is used, the computer needs less RAM capacity than for a genome.    


- A GFF format annotation file corresponding to the reference genome dropped.  
For example, for the human genome annotations, you can look for [Human genome](https://www.ensembl.org/Homo_sapiens/Info/Index) and download the *Homo_sapiens.GRCh38.103.gff3.gz* file in the genome annotation list of gff files, unzip it and place it in the *database* subfolder of your project directory.  
- out-RNA fasta file: This file must gather together RNA sequences you want to remove from the analysis. As a rule, these are at least rRNA. You can also add mitochondrial RNA, non-chromosomal dna or any other fasta sequence of your choice. You can find the fasta files associated to non-chromosomal dna in the [Ensembl](https://www.ensembl.org/index.html) database in the genome assembly list of fasta files. If you want to remove some specific sequences, you just have to create a file with, for each sequence, one line starting with a ">" where you can add a name for your sequence followed by one line containing the sequence to remove. It gives you this format :  
>&emsp; > Sequence_X
>&emsp; GCTGACACGCTGTCCTCTGGCGACCTGTCGTCGGAGAGGTTGGGCCTCCGGATGCGCGCGGGGCTCTGGC
>&emsp; CTCACGGTGACCGGCTAGCCGGCCGCGCTCCTGCCTTGAGCCGCCTGCCGCGGCCCGCGGGCCTGCTGTT  
>&emsp; > Sequence_Y
>&emsp; CTCTCGCGCGTCCGAGCGTCCCGACTCCCGGTGCCGGCCCGGGTCCGGGTCTCTGACCCACCCGGGGGCG  

If you look for specific sequences, you can find them on the [NCBI website](https://www.ncbi.nlm.nih.gov/). For example, you can find the rRNA fasta file sequences [here](https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta).  
You need to have a file, even if it's empty. So, if you want everything to be used in the analysis, just put an empty file.

If you have them, files containing each annotation length (see next paragraph) are also be dropped into this folder.  

Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── Wild_Type.1.fastq.gz   
│&emsp;&emsp;├── Wild_Type.2.fastq.gz   
│&emsp;&emsp;├── Mutant.1.fastq.gz   
│&emsp;&emsp;└── Mutant.2.fastq.gz  
└── database   
&emsp;&emsp;&emsp;├── reference_genome_sequences.fa  
&emsp;&emsp;&emsp;├── reference_genome_annotations.gff3  
&emsp;&emsp;&emsp;├── RNA_to_remove.fa  
&emsp;&emsp;&emsp;└── annotation_length.txt (if needed/provided)  

### c) [config.yaml](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) file  
Config.yaml file is used to define parameters to tell RiboDoc how to process your data.  
You must download it [here](https://raw.githubusercontent.com/equipeGST/RiboDoc/main/config.yaml) and open it with a text editor.    
It must be carefully completed and be present in the project directory everytime you want to run RiboDoc.  

>Caution  
>&emsp;&emsp;&emsp;Spaces and quotation marks **must not be changed** ! Your information must be entered in quotes and should not have spaces     

#### Project name  
First and easy step, the project name ! You can use the same as your folder.  
*project_name*: principal_directory_name  
#### Name of database subfolder file  
You must enter the full name (**with extensions**) without the path of files added in the database subfolder previously created.   
*fasta*: reference_genome_fasta_file.fa  
*gff*: corresponding_gff_annotation_file.gff3  
*fasta_outRNA*: unwanted_sequences_fasta_file.fa  
#### Pipeline option selection  
During the RiboDoc process, data is trimmed and selected depending on their length.   
*already_trimmed*: If your data contains already trimmed reads, you can set this option on “yes”.   
*adapt_sequence*: If they are not trimmed, you must specify the sequence adapters in quotes on the line here.

You also have to define the range for read length selection. Default values select reads from 25 to 35 bases long.  
*readsLength_min*: minimum read length.   
*readsLength_max*: maximum read length.   

You can also parameter alignments:   
*gff_element_cds*: feature corresponding to CDS in the annotation file. 'CDS' is the default value (can sometimes be ORF).     
*gff_attribut*: attribut to regroup reads during counting. 'ID' is the default value.     
#### Statistical settings  
To be able to perform statistical analyzes, you must define a reference condition as well as your thresholds.   
*reference_condition*: it correspond to the reference biological_condition_name
*transcript_or_gene*: choose whether to perfom the differential analysis on features grouped by "transcripts" or by "gene"  
*p-val*: p-value threshold for the differential analysis. Defaut value is 0.01.  
*logFC*: logFC threshold for the differential analysis. Defaut value is 0 to keep all the genes without logFC filtering.  
#### Window for qualitative test  
During the quality analysis, the periodicity is observed on bases around start and stop.     

>The periodicity must be calculated using a metagene profile. It provides the amount of footprints relative to all annotated start and stop codons in a selected window.   

2 pipelines dedicated to quality controls are available in RiboDoc. The first one uses the [riboWaltz tool](https://github.com/LabTranslationalArchitectomics/riboWaltz) which can need high RAM resources depending on your data. The second pipeline is a series of scripts called TRiP which use a specific gff file format as annotation file. For more details on this format, please check the RiboDoc article.
*qualitative_analysis*: Choose between the 2 qualitative analysis pipeline. Default value is "ribowaltz".  
##############################################################
##### Optional and only for qualitative_analysis: "trip" #####
##############################################################
#### UTR covering option  
*UTR*: This option has to be turned on if you want to compare UTRs coverage against CDS. However, to be realized, this option requires a file with the name of each gene and the length of the associated annotation (one file by annotation: CDS, 3’-UTR and 5’-UTR).  

The full name (**with extensions**) without the path of each file has to be report in the configfile.  
*CDSlength*: complete name of file with CDS length.  
*5primelength*: complete name of file with 5’-UTR length  
*3primelength*: complete name of file with 3’-UTR length  

Elements for UTRs have to be specified only if you put 'yes' for the 'UTR' option  
*gff_element_three_prime_utr*: feature corresponding to 3'UTR in the annotation file. Default value is "three_prime_UTR"  
*gff_element_five_prime_utr*: feature corresponding to 5'UTR in the annotation file. Default value is "five_prime_UTR"  
##############################################################  
The window selected by default is -50/+100 nts and -100/+50 nts around start and stop codons respectively.   
*window_bf*: define your window before start and after stop  
*window_af*: define your window after star and before stop  
#### Number of threads  
Thanks to the use of Snakemake, RiboDoc can analyse multiple samples at the same time. We define that ¼ of available CPUs are necessarily requisitioned for these multiple parallele tasks.  
As some tools used in the analysis pipeline have a multithreaded option (like cutadapt, hisat2, bowtie2 or htseq-count), you can choose the number of threads you allow.    
*threads*: threads factor you allow. Default value is "2" to enable 50% of your CPUs to be used in the analysis.    

>Caution:  
>&emsp;&emsp;&emsp;You can attribute 3 threads maximum.  
>&emsp;&emsp;&emsp;Indeed:  
>&emsp;&emsp;&emsp;Allow 1 thread = Remain on a ¼ of threads used  
>&emsp;&emsp;&emsp;Allow 2 threads = Each sample will have 2 threads at its disposal. 2/4 of the CPUs will therefore be used.   
>&emsp;&emsp;&emsp;Allow 3 threads = ¾ of CPUs are used. You should never go beyond so as not to saturate the computer.     

Folder architecture at this step:  
Project_name  
├── fastq   
│&emsp;&emsp;├── Wild_Type.1.fastq.gz   
│&emsp;&emsp;├── Wild_Type.2.fastq.gz   
│&emsp;&emsp;├── Mutant.1.fastq.gz   
│&emsp;&emsp;└── Mutant.2.fastq.gz  
├── database   
│&emsp;&emsp;├── reference_genome_sequences.fa  
│&emsp;&emsp;├── reference_genome_annotations.gff3  
│&emsp;&emsp;├── RNA_to_remove.fa  
│&emsp;&emsp;└── annotation_length.txt (if needed/provided)  
└── config.yaml  


## 3) Pull RiboDoc  
When the folder architecture is ready, it’s time to start RiboDoc !  
First, open a terminal (this can look scary for some but don't worry, there are only 2 lines to write).  
If you never used RiboDoc on your workstation, you must pull it from Docker hub.  
Copy and past the following command line:    
&emsp;&emsp;&emsp;`docker pull equipegst/ribodoc`  
If you have any error, it might come from a rights problem so you should try to copy and paste this command:    
&emsp;&emsp;&emsp;`sudo docker pull equipegst/ribodoc`     

## 4) Run RiboDoc  
Then, or if you already have used RiboDoc, you can run it thanks to the following command:  
&emsp;&emsp;&emsp;`docker run --rm -v /path/to/working/directory/:/data/ equipegst/ribodoc`  
*/path/to/working/directory/* corresponds to the project_name directory's **full path** (usually starting with a "/").   
*:/data/* **must not be modified in any way** as it corresponds to the path to the project directory inside the container.   

>Caution:   
>&emsp;&emsp;&emsp;The path to your project folder and the "/data/" path must start and finish with a slash "/"   


>For Windows users  
>&emsp;&emsp;&emsp;The path to your project folder has to start at the local disks C or D: C:\your\path\   
>&emsp;&emsp;&emsp;This path to your project folder has to be composed and finished with backslashes "\" (instead of slashes "/")   
>&emsp;&emsp;&emsp;/data/ path does not change in any way !   

## 5) In case of any error   
Managed by snakemake, the pipeline will finish all jobs unrelated to the rule/job that failed before exiting. You can still force the container to stop with Docker Desktop or with the following command lines (might need the "sudo" keyword at the beginning) :   
>&emsp;> docker container ls   

Which provides you the container's ID (*e.g.* 9989909f047d), then :   
>&emsp;> docker stop ID  
Where "ID" is the id obtained with the previous command

If you have issues with the use of Docker, you must refer to their [website](https://docs.docker.com/).     
If the error happens during the use of RiboDoc, the rule (job) which failed is indicated in your terminal. You can then find the error output in the *logs* folder. Each rule have a precise name and a folder related to it with files corresponding to the different steps of this rule. If you can not solve the problem by exploring those files, you can contact us through the ["issues"](https://github.com/equipeGST/RiboDoc/issues) part of our github.
>If you want to send us a request because of an error, the easiest way for us to help you is if you send us an archive named after the rule which failed in your analysis with your config.yaml, logs folder and logsTmp folder to it so we can help you finding what happened and why.

## 6) Understand the results  
Here is the project_name folder architecture after RiboDoc run.  
Initial folders and files are still present and highligth in bold in the tree architecture below.  
**Project_name**  
├── **fastq/**   
│&emsp;&emsp;├── **Wild_Type.1.fastq.gz**   
│&emsp;&emsp;├── **Wild_Type.2.fastq.gz**   
│&emsp;&emsp;├── **Mutant.1.fastq.gz**   
│&emsp;&emsp;└── **Mutant.2.fastq.gz**  
├── database/   
│&emsp;&emsp;├── **reference_genome_sequences.fa**  
│&emsp;&emsp;├── **reference_genome_annotations.gff3**  
│&emsp;&emsp;├── **RNA_to_remove.fa**  
│&emsp;&emsp;└── **annotation_length.txt (if provided)**   
├── **config.yaml**   
├── logs/   
│&emsp;&emsp;└── *one_log_folder_by_job*  
├── logsTmp/   
│&emsp;&emsp;└── *one_file_by_steps_of_interest_for_alignment_stats*  
├── RESULTS/  
│&emsp;&emsp;├── config.yaml     
│&emsp;&emsp;├── BAM/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_bam_by_sample.bam*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_bai_by_bam.bai*  
│&emsp;&emsp;├── DESeq2/  
│&emsp;&emsp;│&emsp;&emsp;├── count_matrix.txt  
│&emsp;&emsp;│&emsp;&emsp;├── complete.txt  
│&emsp;&emsp;│&emsp;&emsp;├── up.txt  
│&emsp;&emsp;│&emsp;&emsp;├── down.txt  
│&emsp;&emsp;│&emsp;&emsp;├── Images/  
│&emsp;&emsp;│&emsp;&emsp;│&emsp;&emsp;├── *all_quantitative_analysis_graphs*  
│&emsp;&emsp;├── fastqc/  
│&emsp;&emsp;│&emsp;&emsp;├── *one_html_by_sample.html*  
│&emsp;&emsp;│&emsp;&emsp;└── *one_zip_by_sample.zip*  
│&emsp;&emsp;├── htseqcount_CDS/  
│&emsp;&emsp;│&emsp;&emsp;└── *one_file_by_sample.txt*  
│&emsp;&emsp;├── riboWaltz/  # If you chose *qualitative_analysis: "ribowaltz"* in the config.yaml file  
│&emsp;&emsp;│&emsp;&emsp;├── *riboWaltz's qualitative analysis results*  
│&emsp;&emsp;├── qualitativeAnalysis/  # If you chose *qualitative_analysis: "trip"* in the config.yaml file  
│&emsp;&emsp;│&emsp;&emsp;├── *TRiP's qualitative analysis results*  
├── dag_all.svg    
└── dag_last_run.svg    
- The *logs/* folder groups together all the error output messages from tools used in RiboDoc analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step to give us feedback.   
- *RESULTS/* folder contains 5 subfolders:   
&emsp;&emsp;i) *BAM/*: it contains a BAM file for each sample (allows visualization on tools such as IGV).  
&emsp;&emsp;ii) *DESeq2/*: it contains the differential analysis html report (*PROJECT_NAME.Final_report.txt*), the count_matrix, the tables and the images related to the DESeq2 analysis.   
&emsp;&emsp;iii) *fastqc/*: it contains raw data quality controls.   
&emsp;&emsp;iv) *htseqcount_CDS/*: it contains htseq output for CDS counts.    
&emsp;&emsp;vii) *qualitativeAnalysis/* or *riboWatlz/*: it contains all files related to qualitative test like periodicity and reads length repartition   
It contains also two files:  
&emsp;&emsp;i) *PROJECT_NAME.Analysis_report.html* gathers standard output of each analysis pipeline tool. It allows to know numbers of reads at each step a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.  
&emsp;&emsp;iii) *config.yaml* to have a parameters backup.     

- The *dag files* which represents the analysis steps with your samples.  

>Last big tip:  
In case that a sample is too variable against other replicats or if new samples sequencing are added to your study, you can delete/move or add them in the *fastq* subfolder, delete/move/rename the subfolder *RESULTS/DESeq2*. Run again RiboDoc on the same *project_name* folder and it only creates missing files to complete the analysis (complete analysis for added samples, new differential analysis with all samples available in the *fastq* subfolder).  

Thank you for using RiboDoc !   
We wish you the best results for your analyzes !  
---
title: "Differential analysis report"
author: "GST team"
date: "`r format (Sys.time(), '%B %d, %Y')`"
output: "html_document"
---

```{r Sys, echo=FALSE}
#Sys.getenv("RSTUDIO_PANDOC")
# Sys.setenv(RSTUDIO_PANDOC="/root/miniconda3/envs/RiboDoc_env/bin/pandoc")
Sys.setenv(RSTUDIO_PANDOC="/RiboDoc/miniconda3/envs/RiboDoc_env/bin/pandoc")
```

```{r rm Global environment, include=FALSE}
rm(list=ls())
```

```{r packages, echo=FALSE, include=FALSE}
library(DESeq2)
library(stringr)
library(FactoMineR)
```

```{r scan_and_paths, echo=FALSE, include=FALSE}
params <- scan(file = "/data/config.yaml",
                   what = "character",
                   sep = ":"
                   )

trans_or_gene <- gsub(" ", "", params[which(params=="transcript_or_gene")+1], fixed = TRUE)
if(trans_or_gene=="gene")
{
  DESeq2_folder <- "/data/RESULTS/DESeq2_by_gene/"
} else {
  DESeq2_folder <- "/data/RESULTS/DESeq2_by_transcript/"
}

pathway_matrix = paste0(DESeq2_folder,"count_matrix.txt")
pathway_names = paste0(DESeq2_folder,"transcript_list.txt")
```

```{r Arguments, echo=FALSE, include=FALSE}
refCond <- gsub(" ", "", params[which(params=="reference_condition")+1], fixed = TRUE)
Var_log2FC <- as.numeric(gsub(" ", "", params[which(params=="logFC")+1], fixed = TRUE))
Var_padj <- as.numeric(gsub(" ", "", params[which(params=="p-val")+1], fixed = TRUE))
```
\
##### *`r paste0("Reference condition : ", refCond)`*
\

```{r mkdir, echo = FALSE, include = FALSE}
dir.create(paste0(DESeq2_folder,"Images/"))
```

```{r setup, echo=FALSE, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(include = TRUE, message = FALSE, warning = FALSE)
knitr::opts_chunk$set(fig.path = paste0(DESeq2_folder,"Images/"), dev=c('png','tiff'))
```


```{r data , echo=FALSE}
expData = round(read.table(pathway_matrix, header = T, row.names = 1, check.names = FALSE))

names_list = read.table(pathway_names, header = T, row.names = 1, check.names = FALSE)
```

```{r Merge by gene, echo=FALSE}
if(trans_or_gene=="gene")
{
  expData_named <- data.frame(merge(expData, names_list, by = "row.names"), row.names = 1)
  expData_named_gene <- expData_named[,-(dim(expData_named)[2]-1)]
  expData_counts_by_gene <- aggregate(expData_named_gene[,-dim(expData_named_gene)[2]], expData_named_gene["Gene_name"], sum)
  expData <- data.frame(expData_counts_by_gene, row.names = 1)
  write.table(expData, paste0(DESeq2_folder,"count_matrix_by_gene.txt"), sep = "\t")
}
```

The present analysis report is part of RiboDoc tool, developed by the GST team, Université Paris-Saclay, CEA, CNRS, Institute for Integrative Biology of the Cell (I2BC), 91198, Gif-sur-Yvette, France.

**OVERVIEW** :  
This report presents the quantitative analysis realized thanks to R software and bioconductor packages including DESeq2, stringr and FactomineR.
First of all the raw data are summarized, then the analysis of the variations between and within biological conditions is carried out. After data normalization, the main results of the differential analysis is presented. 

**RAW DATA** :  
First of all, the count matrix obtained with RiboDoc tool is loaded in R. It contains one column per sample and one row per feature.
\
`r nrow(expData)` features are compared in this study.

\


```{r Determination of the number of samples , include=FALSE }
expData_Sorted <- expData[,sort(colnames(expData))]
Names_Col <- colnames(expData_Sorted)
Cond <- str_replace(Names_Col,regex("[.]{0,1}[[:digit:]]{1,}$"),"")
nb_col_condA <- table(factor(Cond))[[1]]
nb_col_condB <- table(factor(Cond))[[2]]
sample <- Cond[nb_col_condA+1]

if (Cond[1] != refCond) {
  Ord <- c( (1:nb_col_condB)+nb_col_condA , (1:nb_col_condA) )
  expData_Sorted <- expData_Sorted[, Ord]
  
  nb_col_condA <- nb_col_condA + nb_col_condB
  nb_col_condB <- nb_col_condA - nb_col_condB
  nb_col_condA <- nb_col_condA - nb_col_condB
  Names_Col <- colnames(expData_Sorted) 
  sample <- Cond[1]
}
```

*Table 1 : Header of the count matrix*
```{r  results = 'asis', echo = FALSE}
knitr::kable(head(expData_Sorted[,]))
```
##### *`r paste0("Samples : ",nb_col_condB," ",Cond[1]," vs ",nb_col_condA," ",Cond[nb_col_condB+1])`*
 \
 \
Figure 1 shows the library sizes. That is the total number of reads count in each sample. Size must be as similar as possible within conditions.  
\
*Figure 1 : Number of raw reads counts per sample. One color by biological condition.*
```{r Library sizes, fig.height = 3, fig.width = 4, echo=FALSE, dev=c('png','tiff')}
barplot(colSums(expData_Sorted/1000000), ylab = "Total read number (million)",
        main = "Library sizes", col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        cex.names = 0.6,
        las = 2
        )
```
\
\
**VARIATION BETWEEN AND WITHIN BIOLOGICAL CONDITIONS** :  
Quantitative analysis is realized in order to highlight the variability between two (or more) biological conditions. To assess this variability, replicats have to be as close as possible. To check this, a pairwise scatter plot (figure 2), a Spearman correlation (figure 2) and a PCA (figure 3) are produced. The first one shows the number of counts for each gene between two samples. Here these plots are associated with Spearman correlation coefficients, which measure the statistical relationship between two samples. The value must be as close as possible to 1 between replicates and move towards 0 between biological conditions. In the same way, the PCA allows us to visualize variability. The first component (figure 3, x-axe) is expected to clearly separate samples between the different biological conditions and bring together replicats.  
\
*Figure 2 : Pairwise scatter plot for variation analysis*
```{r pairwiseScatter, echo = FALSE ,dev=c('png', 'tiff') }
panel.cor <- function(x, y){
   r <- round(cor(10^x, 10^y, method="spearman"), digits=4 )
    txt <- paste0("R = ", r)
   text(2, 2, txt, cex = 1)
}

upper.panel<-function(x, y){
 points((x),(y), pch = ".", cex =1, col = rgb(0, 0, 0, 0.15))
}

pairs(log10(expData_Sorted), 
      lower.panel= panel.cor,
      upper.panel = upper.panel
      )
```
\
*Figure 3 : PCA for variation analysis*  
```{r PCA,  echo=FALSE  }
par(mfrow = c(2,1))
res.pca <- PCA(t(expData_Sorted) , ncp = 3, graph = FALSE)

plot(res.pca, choix ="ind", autoLab = "yes", axes = c(1,2),width=3, height=3)  

plot(res.pca, choix ="ind", autoLab = "yes", axes = c(1,3))
par(mfrow = c(1,1))
```
\
\
**DATA NORMALIZATION** :  
A DESeqDataSet (DDS) object is created from raw data in terms of conditions. 
Thanks to the “counts” function from DESeq2, we normalize data through the DDS object. It is necessary to erase technical biases and make read counts comparable between samples.  
The DESeq2 normalization uses the median of ratios method. For each feature, a pseudo-reference sample is created (ref=sqrt(featureCount_sampleA*featureCount_sampleB)). Then the ratio of each sample to the reference is calculated (ratio=sample/ref). The normalization factor (or size factor) for each sample corresponds to the median value of all ratios for a given sample (normalization_factor_bySample=median(all_feature_ratio)). One sample raw counts are then divided by its size factor to determine the normalized counts. The median of ratios method is based on the hypothesis that not all features are differentially expressed. So, median is not influenced by outliers which correspond to differentially expressed genes without biological conditions distinction.  
To check that normalization went well, we realized two graphs. The first one (figure 4) shows the library size as figure 1 but after normalization. All samples must have the same size. Boxplots are also generated (figure 5) to show how counts distributions changed between before and after normalization. We expect that normalized counts are nearly the same between all samples unlike raw data.  


*Figure 4 : Number of normalized reads counts per sample.*  


```{r, echo=FALSE}
conds = factor(c(rep("CondA", nb_col_condA), rep("CondB", nb_col_condB)))
colData = data.frame(condition = conds)
ddsObjet = DESeqDataSetFromMatrix(countData = expData_Sorted,
                                  colData   = colData, formula(~ condition))
ddsObjet = estimateSizeFactors(ddsObjet)
Size_Factors <- sizeFactors(ddsObjet)
# Normalised data
normCountData = counts(ddsObjet, normalized = TRUE)
```
\


```{r Library Sizes Norm, fig.height = 3, fig.width =4, echo=FALSE}
barplot(colSums(normCountData/1000000), ylab = "Total read number (million)",
        main = "Library sizes (after normalization)",
        col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        las = 2,
        cex.names = 0.6)
```
\
*Figure 5 : Boxplots of reads distribution. Raw (left) vs normalized (right).*
```{r Boxplot, echo=FALSE}
row_sub_1 = apply(expData_Sorted, 1, function(row) all(row != 0 ))
row_sub_2 = apply(normCountData, 1, function(row) all(row != 0 ))
par(mfrow = c(1,2))

boxplot((expData_Sorted[row_sub_1,]),
        main = "Library sizes",
        log = "y",
        col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        cex.names = 0.6,
        las =2
        ,ylim = c(1,max(expData_Sorted[,]) )
        )
boxplot((normCountData[row_sub_2,]), ylab = "Total read number",
        main = "Library sizes (after normalization)",
        log = "y",
        col = c(rep("yellow", nb_col_condA), rep("blue", nb_col_condB)),
        names = colnames(expData_Sorted),
        cex.names = 0.6,
        las =2,
        yaxt = "n",
        ylim = c(1,max(expData_Sorted[,]) )
        )
```


```{r pagination graph cancellation, include = FALSE}
par(mfrow = c(1,1))
```
\


**DIFFERENTIAL ANALYSIS** :  
The DESeq2 function can be run now on the DDS object. First, the estimateSizeFactors sub-function calculates the relative library depth of each sample. Then, estimateDispersions sub-function estimates the dispersion of counts for each feature. Finally nbinomWaldTest sub-function calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs.
Results are presented in seven graphs. Figure 6 represented the dispersion estimate with the mean of normalized counts by the estimated dispersion. Dots can have three colors : i) black which shows the dispersion estimates by feature. ii) red for representing the fitted dispersion value estimated by model (i.e. the mean-variance relationship function) iii) Blue exhibits the final estimates retained for the statistical test. Outliers are shown with a black dot with a blue circle.  
\
*Figure 6 : Dispersion estimates*
```{r , echo=FALSE, include = FALSE}
ddsEstim = DESeq(ddsObjet)
resDESeq = results(ddsEstim, contrast = c("condition", "CondB", "CondA"))

mcols(resDESeq, use.names = TRUE)
dds = estimateDispersions(ddsObjet)
```


```{r Dispersion, fig.height = 5, fig.width = 5, echo=FALSE}
plotDispEsts(dds)
```

  
*Figure 7 : shows the distribution of logFC frequency. The highest frequency is expected on 0. In fact, the large majority of features must not be differentially expressed.*
```{r ,  echo=FALSE, include = FALSE}
# Control mean
CT_Mean<- data.frame(apply(normCountData[,1:(nb_col_condA)],1,mean))

# Tested condition mean
Mut_Mean<- data.frame(apply(normCountData[,(nb_col_condA+1):dim(normCountData)[2]],1,mean))
```

```{r logFC frequency, fig.height = 3, fig.width =4, echo=FALSE}
hist(resDESeq$log2FoldChange,
     nclass = 200,
     xlab = "logFC values",
     main ="Frequency / LogFC" )
abline(v = 2, col = "red")
abline(v = -2, col = "red")
```
\
Figures 8 and 9 represent the raw and adjusted p-values distribution respectively. The peak around 0 corresponds to the differentially expressed genes. Graphs are expected to be a uniform distribution between 0 to 1.  

*Figure 8 : Distribution of raw p-values*

```{r pValue, fig.height = 4, fig.width =4, echo=FALSE}
hist(resDESeq[,"pvalue"], nclass = 100, xlab = "p-values",
     main = "Histogram of p-values (DESeq2)")
```
  
*Figure 9 : Distribution of adjusted p-values*  

```{r adjusted pValue, fig.height = 4, fig.width =4, echo=FALSE}
hist(resDESeq[,"padj"], nclass = 100, xlab = "padj",
     main = "Histogram adjusted p-values")
```
\
Figure 10 is the MA-plot which shows the mean of normalized counts for each feature by the log ratio of differential expression. Differentially expressed features are represented by red dots. Triangles correspond to features with a too high/low log2FC to be shown on graph.    
   
*Figure 10 : MA-plot*   
```{r MA-plot, fig.height = 4, fig.width =4, echo=FALSE}
plotMA(resDESeq, alpha = Var_padj, ylab = "log2 FC", colSig = rgb(1,0,0,0.5))
```


Figures 11 and 12 are the volcano-plot and its zoom respectively. They represent each feature by its log2FC and its adjusted p-value. Differentially expressed features are still red dots and triangles correspond to outliers features.   
*Figure 11 : Volcano-plot*

```{r volcano plot global, fig.height = 4, fig.width =4, echo=FALSE}
myColors <- ifelse((resDESeq$padj < Var_padj & resDESeq$log2FoldChange > Var_log2FC), rgb(1, 0, 0, 0.25) ,
                   ifelse((resDESeq$padj < Var_padj & resDESeq$log2FoldChange < -Var_log2FC), rgb(1, 0, 0, 0.25) ,
                          rgb(0, 0, 0, 0.25) ) )

plot(resDESeq[, "log2FoldChange"], -log10(resDESeq[, "padj"]),
     pch = 20, cex = 1,
     col=myColors,
     xlab = "log2 FC", ylab = "-log10(padj)",
     main = paste0("Volcano plot - ",Cond[nb_col_condB+1]," vs ",Cond[1])
     )
```
  
*Figure 12 : Zoom of the previous volcano-plot*  

```{r ZOOM volcano plot , fig.height = 4, fig.width =4, echo=FALSE}
limite <- Var_padj/(10^(-log10(Var_padj)*2))


myColors <- ifelse((resDESeq$padj < Var_padj & resDESeq$log2FoldChange > Var_log2FC), rgb(1, 0, 0, 0.25) ,
                  ifelse((resDESeq$padj < Var_padj & resDESeq$log2FoldChange < -Var_log2FC), rgb(1, 0, 0, 0.25) ,
                          rgb(0, 0, 0, 0.25) ) )

Valeurs_Limite_Y <- ifelse((resDESeq$padj < limite), limite ,
              resDESeq$padj) 

Forme_Pixel_Y <- ifelse((resDESeq$padj < limite), 17 ,
                                 20) 


plot(resDESeq[, "log2FoldChange"], -log10(Valeurs_Limite_Y),
    pch = Forme_Pixel_Y,
    col=myColors,
    xlab = "log2 FC", ylab = "-log10(padj)",
    ylim = c(0,-log10(limite)),
    main = paste0("Volcano plot Zoom - ",Cond[nb_col_condB+1]," vs ",Cond[1])
    )
```
\

```{r Final tables creation, include = FALSE}
Bruts_Norm <- cbind(expData_Sorted, normCountData )

Table_Complete <- data.frame(cbind(row.names(Bruts_Norm) , Bruts_Norm, CT_Mean,Mut_Mean , resDESeq) )


colnames(Table_Complete) =  c("ID" ,
                               Names_Col,
                               paste("norm", Names_Col, sep = "_"),
                               paste0(refCond,"_Mean"),
                               paste0(sample, "_Mean"),
                               colnames(resDESeq))

if(trans_or_gene=="gene")
{
  allGenes <- Table_Complete
} else {
  allGenes <- merge(Table_Complete, names_list, by = "row.names")
}

inducedGenes = allGenes[which((allGenes[, "log2FoldChange"] > Var_log2FC) & (allGenes[, "padj"] < Var_padj)),]
dim(inducedGenes)

repressedGenes = allGenes[which((allGenes[, "log2FoldChange"] < -Var_log2FC) & (allGenes[, "padj"] < Var_padj)),]
dim(repressedGenes)

# Sort tables by padj
inducedGenes <- inducedGenes[(order(inducedGenes[,"padj"])),]
repressedGenes <- repressedGenes[(order(repressedGenes[,"padj"])),]
```


```{r, echo = FALSE}
write.table(allGenes,
            file = paste0(DESeq2_folder,"complete.txt"),
            quote = F, sep = "\t", row.names = F)


write.table(inducedGenes,
            file = paste0(DESeq2_folder,"up.txt"),
            quote = F, sep = "\t", row.names = F)

write.table(repressedGenes,
            file = paste0(DESeq2_folder,"down.txt"),
            quote = F, sep = "\t", row.names = F)
```
\
**OUTPUT FILES** :  
All the output files can be found in the directory “RESULTS/DESeq2_by_gene/” or " RESULTS/DESeq2_by_transcript/”  
`r paste0("RESULTS/DESeq2_by_",trans_or_gene,"/complete.txt")`  
`r paste0("RESULTS/DESeq2_by_",trans_or_gene,"/up.txt")  `  
`r paste0("RESULTS/DESeq2_by_",trans_or_gene,"/down.txt")  `   
\
*Explanation of each column* :  

* Genes: unique feature identifier  
* sampleName: raw counts per sample  
* norm.sampleName: rounded normalized counts per sample  
* Mean.sampleName: mean over all samples  
* baseMean: base mean over all samples  
* log2FoldChange: log^2^(FC) reflects the differential expression between Test and Ref
* pvalue: raw p-value from the statistical test  
* padj: adjusted p-value on which the cut-off α is applied  
