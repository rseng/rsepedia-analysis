# Easy RNAseq Analysis with *oneStopRNAseq*
- This is the backend of https://mccb.umassmed.edu/OneStopRNAseq/index.php
- Can be installed on Linux systems (tested on Ubuntu and CentOS, did not test on other flavors of Linux, not working on Mac yet)
- Would need knowledge in basic bash commands to use this workflow

# Installation
- install anaconda by following instructions on https://docs.anaconda.com/anaconda/install/  (installation tested in conda version 4.9.2)
- download OneStopRNASeq by `git clone https://github.com/radio1988/OneStopRNAseq.git`
- `cd OneStopRNAseq/snakemake`
- `conda env create -f workflow/envs/osr-base.yaml`  # create an conda env called 'osr-base'
- `conda activate osr-base`
- `Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'` # install packages for QoRTs
- uncompress example dataset
    - `gunzip example_data/genome/mm10_chr19/mm10.chr19.fa.gz` 
    - `gunzip example_data/genome/mm10_chr19/gencode.vM25.primary_assembly.annotation.ch19.gtf.gz`
- the installation typically takes 15-30 mins

# Running OneStopRNASeq workflow on example datasets
## Start from FASTQ as input
```
mkdir fq_analysis && cd fq_analysis # create workdir
osr_path=$download_path/OneStopRNASeq/snakemake  # this is user specific, e.g. /home/user/git/OneStopRNAseq/snakemake
# put necessary files into workdir
ln -s $osr_path/meta
ln -s $osr_path/example_data
ln -s $osr_path/example_data/fastq_small/ fastq
cp $osr_path/config.fq_example.yaml config.yaml
mkdir -p workflow && cd workflow
ln -s $osr_path/workflow/envs
ln -s $osr_path/workflow/Snakefile
ln -s $osr_path/workflow/osr.py
rsync  -a $osr_path/workflow/script ./
snakemake -j 1 -np   # quick test with a 'dry-run'
snakemake -j 2 -pk  # run the workflow on the example datase with two threads, takes around 30 min for the first run
```

# Quick start on snakemake advanced usage:
- If the workflow did not finish completely, try submitting the jobs again with `snakemake -j 2 -pk`
- If the workdir is locked in a second submission, please kill previously submitted jobs by typing 'snakemake --unlock -j 1', after you make sure the first submission has stopped running
- The workflow can be easily adapted to a LSF system by using `snakemake -pk --cluster 'bsub -q queue_name -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 24:00'` 

# Viewing results
- results are under subfolders in the workdir, e.g. DESeq2, gsea, rMATS.x, fastqc, bam_qc
- interpretation of results can be found here: https://mccb.umassmed.edu/OneStopRNAseq/documents/description_of_output_files.pdf 
- write up for method section in publication can be found here: https://mccb.umassmed.edu/OneStopRNAseq/documents/template_of_method_section.pdf 
# SalmonTE

## Change Logs

* December 19, 2018: `SalmonTE` has been updated to `0.4` with some improvements!
  * Regarding [#14](https://github.com/hyunhwaj/SalmonTE/issues/14), now `SalmonTE` is coupled with the latest version of `snakemake`, so there is no more directory error. Furthermore, `SalmonTE` no longer runs with the older version of `snakemake`, please update the version of the `snakemake` package.
  * Regarding [#22](https://github.com/hyunhwaj/SalmonTE/issues/22), mappability for each file is now reported. (MAPPING_INFO.csv in the quantification output directory)
  * `test` function has been improved, and this will give better representations of the data.
  * Running of `SalmonTE` will not produce massive and enigmatic messages anymore.

* June 4, 2018: Now `SalmonTE` supports `fq` and `fq.gz` extensions.

* May 3, 2018: Update README.txt

* May 2, 2018: A bug fix regarding issue [#10](https://github.com/hyunhwaj/SalmonTE/issues/10)

* January 23, 2018: Added references of *Mus musculus*(mm) and *Danio rerio*(dr), added a function users to allow to build a customized index (`index` mode), and fixed a minor bug.

* December 25, 2017: Fixed a bug for single-end reads dataset.

* December 14, 2017: Fixed issue of the `macOS`. We have figured out there is a problem of old version of `snakemake`. If you already install the package, and the version is not above `4.0.0` (You can check it with `snakemake --version`) then please update version with below command:
`pip3 install snakemake --user --upgrade`

* November 27, 2017: source code of PSB manuscript is out now - our [manuscript](https://github.com/hyunhwaj/SalmonTE-manuscript/) and [response letter](https://github.com/hyunhwaj/SalmonTE-response).

* November 21, 2017: Version 0.2 is out, `SalmonTE` now supports two different type of expressions with `--exprtype` option. Use `--exprtype=TPM` if you want to use TPM values for the statistical analysis. If you want to run differential expression analysis, then I highly recommend to use `--exprtype=count`. [Here is the nice answer why](https://support.bioconductor.org/p/98820/).


## Notice

* June 5, 2018: Our recent paper about TE in Alzhimer's disease has been published in Cell Reports! Please read the paper to see how `SalmonTE` was succesfully applied! [(Cell Reports Link)](https://www.cell.com/action/showAbstract?pii=S2211-1247(18)30722-8)

## What is SalmonTE?
`SalmonTE` is an ultra-Fast and Scalable Quantification Pipeline of Transpose Element (TE) Abundances from Next Generation Sequencing Data. It comes with [Salmon](https://github.com/COMBINE-lab/salmon) which is a fast and accurate transcriptome quantification method. You can read the details of the pipeline and an example of real data study in [my recent published paper in PSB 2018](http://www.worldscientific.com/doi/10.1142/9789813235533_0016).

## What I need to run `SalmonTE`? Why I have to use it?
* You only need to have a set of FASTQ files and phenotype data. Furthermore, **`SalmonTE` automatically decided wether your dataset is paired-ends reads or not.** 
* conditions can be a numeric data or a categorical data. Based on the data type of the conditions of each sample, `SalmonTE` will run differential expression analysis or linear regression.
* Unlikely other TE analysis tools, `SalmonTE` gives you various visualized output. It must be helpful to your research.

## Requirements & Installation
To use `SalmonTE` `python` and `R` must be installed before running it.

~* **Note**: Currently, running `SalmonTE` on MacOS has an issue, and we are try to fix it soon. Thus, we recommend to use `linux` environment to play it.~

* Install `python` and `R` packages

For `python`:

Run following line in your console
```
pip3 install snakemake docopt pandas --user
```

For `R`:
Run following lines in `R` console.

```
install.packages(c("tidyverse", "scales", "WriteXLS", "BiocManager"))
BiocManager::install("DESeq2", version = "3.8")
```


* Clone the repository 

```
git clone https://github.com/hyunhwaj/SalmonTE
```

* Add `PATH` of SalmonTE to your `.bashrc` file:

```
export PATH=$PATH:/PATH_OF_SALMON_TE/
```

* Re log-in to terminal or use `source` command:

```
source ~/.bashrc
```

### Troubleshooting

**Q. I am using `SalmonTE` on `macOS` and `salmonTE` fails to run on `quant` mode with error messages:**

```
CalledProcessError in line xx of SOME_PATH:
Command ' set -euo pipefail;  ROOT_OF_SALMON_TE/SalmonTE/salmon/darwin/bin/salmon quant...' returned non-zero exit status 134.
```

**A.** You may have a problem to run `salmon` which is an essential tool for the pipeline. You may install [Threading Building Blocks library](https://www.threadingbuildingblocks.org/download) to solve the problem. If you are using [homebrew](https://brew.sh) then please use below command:

```
brew install tbb
```



## How to use it?

```
Usage:
    SalmonTE.py index [--ref_name=ref_name] (--input_fasta=fa_file) [--te_only]
    SalmonTE.py quant [--reference=genome] [--outpath=outpath] [--num_threads=numthreads] [--exprtype=exprtype] FILE...
    SalmonTE.py test [--inpath=inpath] [--outpath=outpath] [--tabletype=tabletype] [--figtype=figtype] [--analysis_type=analysis_type] [--conditions=conditions]
    SalmonTE.py (-h | --help)
    SalmonTE.py --version

Options:
    -h --help     Show this screen.
    --version     Show version.
```

## An example of SalmonTE usage with command line 

### Running the `quant` mode to collect TE expressions

#### Parameters

* `--reference`: This will select a reference file, and should the species **identifier** of your data. We are currently supporting references of those species.
    * hs : *Homo Sapiens*
    * mm : *Mus musculus*
    * dm : *Drosophila melanogaster*
    * dr : *Danio rerio*
* `--outpath`: Path to generate quantification results. If you omit this, `SalmonTE_output` in the current path, and output will be stored in the path.
* `--exprtype`: The type of expression measure, and **TPM** or **count** are possible options. If you omit this, then "TPM" is the input of the parameter.
* `--num_threads`: This has to be an integer, and this parameter will configure how many threads will use for the parallization.

After you put your parameters, you can put the directory which includes a list of **FASTQ** files, 

```
SalmonTE.py quant --reference=hs example
```

Or, you can put the list of files like below.

```
SalmonTE.py quant --reference=hs example/CTRL_1_R1.fastq.gz example/CTRL_2_R1.fastq.gz          
```

### Running `test` mode to perform statistical test

Before you run test mode, you should modify <strike>`control.csv`</strike> `condition.csv` file which is stored in the `outpath`. Here are examples of the proper modifications:

For the differential expression analysis, change the file as below. 
**Important**: The control samples has to be labeled as `control`. Other labels will cause errors.

```
SampleID,condition
FASTQ1,control
FASTQ2,control
FASTQ3,treatment
FASTQ4,treatment
```

For the regression analysis, 

```
SampleID,condition
FASTQ1,1.5
FASTQ2,2.1
FASTQ3,3.8
FASTQ4,9.5
```

Once the conditions of every sample has been filled, we can run the test mode like the example commnad-line below:

* `--inpath`: This should be the path which contains output of `quant` mode. 
* `--outpath`: This will be the path to store all outputs for the mode.
* `--tabletype`: The file format of the tables, `csv`, `tsv`, and `xls` are supported. If you omit this, then `xls` formatted file will be generated.
* `--tabletype`: The file format of the figures, `png`, `tiff`, `jpg`, and `pdf` are supported. If you omit this, then `pdf` formated files will be generated.
* `--analysis_type`: The type of the analysis, and **DE** (for a differential analysis) or **LM** (for a linear regression analysis) are possible options. If you omit this, then "DE" is the input of the parameter.
* `--conditions`: The list of conditions will be considered when **DE** has been selected for `--analysis_type.` The input needs to contain two different conditions (written in non-white space characters) and each condition are separated by `,`, and **no white-space characters** are not allowed for the input. i.e. `SalmonTE` does not care input such as `control , treatment`. The first condition of the input will be considered as a normal condition (e.g. healthy condition, wild-type mice) in the study, and the later will be considered as another condition which you are interested (e.g. knock-out mice, treatment).

```
SalmonTE.py test --inpath=SalmonTE_output --outpath=SalmonTE_statistical_test --tabletype=csv --figtype=png --analysis_type=DE --conditions=control,treatment
```

## How to Cite?

```
@inbook{doi:10.1142/9789813235533_0016,
author = {Hyun-Hwan Jeong and Hari Krishna Yalamanchili and Caiwei Guo and Joshua M. Shulman and Zhandong Liu},
title = {An ultra-fast and scalable quantification pipeline for transposable elements from next generation sequencing data},
booktitle = {Biocomputing 2018},
chapter = {},
pages = {168-179},
doi = {10.1142/9789813235533_0016},
URL = {http://www.worldscientific.com/doi/abs/10.1142/9789813235533_0016},
eprint = {http://www.worldscientific.com/doi/pdf/10.1142/9789813235533_0016}
publisher = WORLD SCIENTIFIC
address = 
year = 2017
edition = 
}
```
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at hyunhwaj@bcm.edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
---
title: "DE Analysis and Standard Visualizations"
output:
  html_document:
    toc: yes
    code_folding: hide
    toc_float: yes
params: 
  max_fdr: 0.05
  min_lfc: 0.585
  indfilter: FALSE
  cookscutoff: TRUE
  countFile: '../example_data/counts/featureCounts.mm10.txt' # all three formats need to be tested in development
  annoFile: 'https://raw.githubusercontent.com/radio1988/OneStopRNAseq/master/snakemake/envs/anno_tables/mm10.gencode.vm25.te_included.anno.txt'
  #annoFile: 'https://raw.githubusercontent.com/radio1988/OneStopRNAseq/master/snakemake/envs/anno_tables/hg38.gencode.v34.te_included.anno.txt'
  metaFile: "../meta/meta.xlsx"
  contrastFile:  "../meta/contrast.de.xlsx"
  blackSamples: 'H_B1_1, H_B1_8, C_B1_4, C_B1_1'
---

## Notes:
- Please search for “Warning” in this report for possible problems before proceeding
-	Normalized expression values, differentially expressed genes, and various types of plots such as heatmap and volcano plots are in the DESeq2 folder. For example, TPM.xlsx contains normalized expression values.
- Please use LFC_shrunken for Log2FoldChange, LFC_raw is not recommended
- Please use FDR (padj) instead of raw p-values for identifying significantly differentially expressed genes.
- blackSamples added in params


```{r setup, include=F}
# Bioconductor
library(BiocManager) # 
library(DESeq2)
library(EnhancedVolcano) # hard to install via conda

# CRAN
library(ggrepel) # 
library(ggplot2)
library(ashr) # 
library(MASS)
library(WriteXLS)
library(plyr) # 
library(gdata) #
library(dplyr)
library(RColorBrewer)
library(pheatmap) #
library(PoiClaClu) # 
library(gridExtra) # 
library(grid)
# library(corrplot) 

#dir.create("./rnk", showWarnings = FALSE)
#knitr::opts_knit$set(root.dir = "./")
```


## Input parameter extraction
```{r params}
max_fdr <- params$max_fdr
min_lfc <- params$min_lfc
indfilter <- params$indfilter
cookscutoff <- params$cookscutoff

countFile <- params$countFile
annoFile <- params$annoFile
metaFile <- params$metaFile
contrastFile <-params$contrastFile

blackSamplesStr <- gsub(" " , '', params$blackSamples)
blackSamples <- unlist(strsplit(blackSamplesStr, ','))

paste("Data file name:", countFile)
paste("MetaData file name:", metaFile)
paste("Contrast file name:", contrastFile)
paste("Annotation file name:", annoFile)
paste("FDR cut-off:", max_fdr)
paste("Log2FC cut-off:", min_lfc)
cat('blackSamples to remove:', blackSamples)
```

## Functions
```{r functions}
fix_hyphen <- function(x){
  return(gsub("-", ".", x))
}

readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
  df <- data.frame(df)  #important
  df <- data.frame(lapply(df, fix_hyphen))
  return (df)
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


volcanoplot <- function(res, anno, name='name', 
                        lfcthresh=2, sigthresh=0.05, 
                        labelsig=FALSE, xlim=100, ylim=1000, textcx=1) {
  # If Annotation file is wrong, No annotation match
  res0 <- res
  res<-merge(data.frame(res0), anno, by.x = 0, by.y=1, all.x=T, sort=F)
  if (sum(anno[, 1] %in% row.names(res0)) < 1) {
    warning(
      c("\nThe annotation file and the count filt does have no match in gene-id:\n", 
        "count table gene-id: ", head(row.names(res0), 1), 
        "\nanno table gene-id: ", anno[1:1, 1], 
        "gene-id rather than gene-name used for Volcano-plot"
      ))
    res$Name <- res[, 1]
  }
  
  # remove NA
  res$padj[is.na(res$padj)] <- 1
  # set lim on x, y
  res$padj[res$padj < 10^(-ylim) & !is.na(res$padj)] <- 10^(-ylim) # y-axis top value 50
  res$log2FoldChange[res$log2FoldChange > xlim] <- xlim
  res$log2FoldChange[res$log2FoldChange < -xlim] <- -xlim
  # show num_pos num_neg
  pos <- subset(res, padj<sigthresh & log2FoldChange>lfcthresh)
  neg <- subset(res, padj<sigthresh & log2FoldChange< -lfcthresh)
  pos.n <- dim(pos)[1]
  neg.n <- dim(neg)[1]
  
  #labelcut <- res$padj[order(res$padj)][100] 
  
  p <- EnhancedVolcano(res,
                       lab = res$Name,
                       #selectLab = as.character(res$Name[which(res$padj<labelcut)]), # mark top genes
                       #selectLab = c("FOS", "LDHA"), # mark selected genes
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = name,
                       subtitle = paste("Up:", pos.n, ", Down:", neg.n, sep = ""),
                       xlab = bquote(~Log[2]~ "Fold Change"),
                       ylab = bquote(~-Log[10]~italic(FDR)),
                       #xlim = c(-6, 6),
                       pCutoff = sigthresh,
                       FCcutoff = lfcthresh,
                       #pLabellingCutoff = labelcut,
                       cutoffLineType = 'twodash',
                       cutoffLineWidth = 0.8,
                       # pointSize = 1.0, # todo: fix incomparibility in R3.5 and R4
                       # labSize = 2.0,
                       # DrawConnectors = F,
                       # legend = c("NS","Log2 FC","P","P & Log2 FC"),
                       legendLabels = c('NS', expression(Log[2]~FC),
                                        "FDR", expression(FDR~and~Log[2]~FC)),
                       caption = paste0('Total = ', nrow(res), ' genes'),
                       legendPosition = 'right',
                       legendLabSize = 10,
                       axisLabSize = 10,
                       legendIconSize = 3.0)
  
  # # plot in html
  # grid.arrange(p,
  #              ncol=1,
  #              top = textGrob(' ',
  #                             just = c('center'),
  #                             gp = gpar(fontsize = 32))
  # )
  # grid.rect(gp=gpar(fill=NA))
  
  # save PDF
  pdf(paste(name, "pdf", sep="."), width=8, height=6)
  grid.arrange(p,
               ncol=1,
               top = textGrob(' ',
                              just = c('center'),
                              gp = gpar(fontsize = 32))
  )
  grid.rect(gp=gpar(fill=NA))
  dev.off()
}


maplot <- function (res, thresh=max_fdr, labelsig=FALSE, textcx=1, ...) {
  # todo: improve visualization
  with(res, 
       plot(baseMean, log2FoldChange, col="grey80", pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<max_fdr), 
       points(baseMean, log2FoldChange, col="grey40", pch=20, cex=.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<max_fdr), 
         textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}

zscore <- function(matrix){
  return( t(scale(t(matrix))))
}

rename_num_vector_by_order <- function(l){
  # l have to be a vector of numbers
  # output vector of roman numbers ordered by appearance in the input vector
  # e.g. c(2,3,3,2,1) -> c(I, II, II, I, III)
  # test rename_num_vector_by_order(c(2,3,3,2,1))
  u <- unique(l)
  n=0
  for (i in u){
    n = n+1; 
    l <- replace(l, l==i, as.character(as.roman(n)))
  }
  return(l)
}

Heatmap <- function(df, nclass=2, fname="heatmap", main="title"){
  if (dim(df)[1] > nclass){
    # Heatmap Pre-plot to get hclust
    p <- pheatmap(df, 
                  main = main,
                  cluster_cols = F, 
                  border_color = NA,
                  cutree_rows = nclass, 
                  show_rownames = F)

    ## Print out gene classification (https://www.biostars.org/p/287512/)
    # cutree manually
    gene_classes <- sort(cutree(p$tree_row, k=nclass))
    gene_classes <- data.frame(gene_classes) 
    classification <- merge(gene_classes, df, by = 0)
    row.names(classification) <- classification$Row.names
    # Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
    idx <- rownames(df[p$tree_row[["order"]],])
    classification <- classification[idx,] 
    # rename gene classes
    classification$gene_classes <- rename_num_vector_by_order(classification$gene_classes)
    # save excel
    WriteXLS(classification, row.names = F,
             paste(fname,"gene_class.xlsx", sep = "."))
    # get class label
    annotation_row = data.frame(
      class = classification$gene_classes
    )
    row.names(annotation_row) <- row.names(classification)
    annotation_row$class <- as.character(annotation_row$class)
    # output final plot
    if (dim(df)[1] < 20){ # manage row height
        p2 <- pheatmap(df,
           main = main,
           border_color = NA,
           cluster_cols = F,
           cutree_rows = nclass,
           show_rownames = F,
           annotation_row = annotation_row,
           cellheight = 20,
           filename = paste(fname, "pdf", sep = "."))
        p3 <- pheatmap(df,  # sample clustering
           main = main,
           border_color = NA,
           cluster_cols = T,
           cutree_rows = nclass,
           show_rownames = F,
           annotation_row = annotation_row,
           cellheight = 20,
           filename = paste(fname, "v2.pdf", sep = "."))
    }else{
        p2 <- pheatmap(df,
           main = main,
           border_color = NA,
           cluster_cols = F,
           cutree_rows = nclass,
           show_rownames = F,
           annotation_row = annotation_row,
           filename = paste(fname, "pdf", sep = "."))
        p3 <- pheatmap(df, # sample clustering
           main = main,
           border_color = NA,
           cluster_cols = T,
           cutree_rows = nclass,
           show_rownames = F,
           annotation_row = annotation_row,
           filename = paste(fname, "v2.pdf", sep = "."))
    }

  }else if (dim(df)[1] > 0){
    # output final plot
        if (dim(df)[1] < 20){ # manage row height
            p2 <- pheatmap(df,
               main = main,
               cluster_cols = F,
               cluster_rows = F,
               border_color = NA,
               show_rownames = F,
               cellheight = 20,
               filename = paste(fname, "pdf", sep = "."))
            p3 <- pheatmap(df,
               main = main,
               cluster_cols = T,
               cluster_rows = F,
               border_color = NA,
               show_rownames = F,
               cellheight = 20,
               filename = paste(fname, "v2.pdf", sep = "."))
        }else{
            p2 <- pheatmap(df,
               main = main,
               cluster_cols = F,
               cluster_rows = F,
               border_color = NA,
               show_rownames = F,
               filename = paste(fname, "pdf", sep = "."))
            p3 <- pheatmap(df,
               main = main,
               cluster_cols = T,
               cluster_rows = F,
               border_color = NA,
               show_rownames = F,
               filename = paste(fname, "v2.pdf", sep = "."))
        }
  }else{ print("no sig DEG for heatmap")}
}


process_deseq_res <- function(res="lfcshrink.res", res2="results.res", name='name', anno='anno.df', norm_exp="tpm.df"){
    ## Summary
    print(name)
    print("\n>>> Summary using FDR cut-off only (LFC not used)")
    summary(res, alpha=max_fdr)
    
    print("\n>>> Summary using both FDR and LFC_shrunken cut-off")
    sig_idx <- res$padj<max_fdr & abs(res$log2FoldChange) > min_lfc
    sig_idx[is.na(sig_idx)] <- FALSE
    res_sig <- res[sig_idx,]
    print(table(sig_idx))
    
    up_idx <- res$padj<max_fdr & res$log2FoldChange > min_lfc
    up_idx[is.na(up_idx)] <- FALSE
    res_sig <- res[up_idx,]
    print(table(up_idx))
    
    down_idx <- res$padj<max_fdr & res$log2FoldChange < -min_lfc
    down_idx[is.na(down_idx)] <- FALSE
    res_sig <- res[down_idx,]
    print(table(down_idx))

    # Prep
    res.df <- as.data.frame(res)
    names(res.df)[2] <- "log2FoldChange_shrunken"
    names(res.df)[3] <- "lfcSE_shrunken"

    res2.df <- as.data.frame(res2)
    names(res2.df)[2] <- "log2FoldChange_raw"
    names(res2.df)[3] <- "lfcSE_raw"
    res2.df <- res2.df[, c(2,3)]

    resdata <- merge(res.df, res2.df, by=0, sort=F, all.x=T)
    resdata <- merge(anno, resdata, by.x=1, by.y=1, sort=F, all.y=T)
    resdata <- merge(resdata, norm_exp, by.x=1, by.y=0, all.x=T, sort=F)
    head(resdata)
    sig_idx <- resdata$padj<max_fdr & abs(resdata$log2FoldChange_shrunken) > min_lfc # important to put this line right before output sig.xlsx
    sig_idx[is.na(sig_idx)] <- FALSE
    resdata.sig <- resdata[sig_idx,]
    head(resdata.sig)

    ## Write results
    WriteXLS(x = resdata,
             ExcelFileName = paste(name, 'deseq2.xlsx', sep = '.'),
             row.names = F, SheetNames = 'sheet1', na = 'NA')  # for user

    WriteXLS(x = resdata.sig,
             ExcelFileName = paste(name, 'deseq2.sig.FDR', max_fdr,
                                   'LFC', min_lfc, 'xlsx', sep = '.'),
             row.names = F, SheetNames = 'sheet1', na = 'NA')  # for user
    ## For GSEA
    rnk <- subset(resdata, select = c("Name","log2FoldChange_shrunken"))
    colnames(rnk) <- c("# Name","log2FoldChange_shrunken")
    rnk <- rnk[order(rnk$log2FoldChange_shrunken), ]
    rnk[, 1] <- toupper(rnk[, 1])
    write.table(rnk, 
                paste('rnk/', name, '.rnk', sep = ''), 
                row.names = F, quote = F, sep='\t')


    # # Corrlelation of Length and LFC for Niraj (only after write excel, alters resdata)
    # resdata.sig.cor <- cor.test(resdata.sig$Length, 
    #                        resdata.sig$log2FoldChange_shrunken, 
    #                        method = "spearman")
    # title <- paste("Spearman Cor:", format(resdata.sig.cor$estimate, digits=2, nsmall=2),
    #                "p-value:", format(resdata.sig.cor$p.value, digits=3, nsmall=3),
    #                sep = " ")
    # 
    # resdata.sig$density <- get_density(resdata.sig$Length, resdata.sig$log2FoldChange_shrunken, n = 100)
    # # set lim (careful, only after outputting excel)
    # resdata.sig$log2FoldChange_shrunken[resdata.sig$log2FoldChange_shrunken > 10] <- 10
    # resdata.sig$log2FoldChange_shrunken[resdata.sig$log2FoldChange_shrunken < -10] <- -10
    # resdata.sig$Length[resdata.sig$Length > 20000] <- 20000
    # ggplot(resdata.sig) + 
    #     geom_point(aes(Length, log2FoldChange_shrunken, color = density)) +
    #     scale_color_viridis() +
    #     ggtitle(paste(paste("Sig-DEG for", name), title, sep = "\n") ) + 
    #     ylim(-10, 10) +        
    #     xlim(0,20000)
    
    ##  Plots
    hist(res$pvalue, breaks=50, col="grey80", # todo: improve vis
         main = paste('Histogram of p-values', name, sep = "\n"), 
         xlab = 'pvalues', ylab = 'Frequency')
    
    hist(res$padj, breaks=50, col="grey", 
         main = paste('Histogram of FDR', name, sep = "\n"), 
         xlab = 'FDR', ylab = 'Frequency')
    
    maplot(res, main=paste("MAplot", paste(name, "LFC_shrunken"), sep="\n")) # todo: improve vis
    maplot(res2, main=paste("MAplot", paste(name, "LFC_raw"), sep="\n"))
    
    volcanoplot(res, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
                textcx=.8,  name= paste(name, "LFC_shrunken", sep="."))
    volcanoplot(res2, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
                textcx=.8,name= paste(name, "LFC_raw", sep="."))
    
    n1 <- dim(resdata.sig)[2]
    n2 <- dim(norm_exp)[2]
    zscore.df <- zscore(resdata.sig[, (n1-n2+1):n1])
    rownames(zscore.df) <- resdata.sig[,1]
    colnames(zscore.df) <- gsub (":TPM", "", colnames(zscore.df))
    Heatmap(zscore.df, nclass = 2,
            fname = paste(name, "heatmap", sep="."),
            main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr ))
    # dev.off()


}
```


## Importing data
```{r read_meta_data}
meta.df <- readExcel(metaFile)
meta.df <- meta.df[!(meta.df[, 1] %in% blackSamples), ]

```

```{r read_count_data}
# read
paste('countFile:', countFile)
if (grepl("txt$", countFile)){
  df <- read.table(countFile, 
                   sep="\t", header=TRUE, comment.char = '#', quote="",
                   row.names = 1) # row.name in cts(matrix)
} else if (grepl("xlsx$", countFile)){
  df <- readExcel(countFile)  
} else {
  stop("input count file not xlsx, nor txt format, abort!!!")
}

# auto-detect content format
featureCountsFormat <- F
cleanCountFormat <- F
osrCountFormat <- F

if (colnames(df)[1] == "Chr" & colnames(df)[2] == "Start" & colnames(df)[3] == "End"){
  cat("featureCounts format table detected\n")
  featureCountsFormat <- T
  # clean names
  colnames(df) <- gsub("\\.bam$", "", colnames(df))
  colnames(df) <- gsub("sorted_reads.", "", colnames(df))
  colnames(df) <- gsub("mapped_reads.", "", colnames(df))
  # convert to int
  df[, 6:ncol(df)] <- sapply(df[, 6:ncol(df)], as.integer)
  # get cts
  cts <- as.matrix(df[, 6:ncol(df)])
} else if (dim(df)[2] == dim(meta.df)[1] + 1){
  cat("clean format table detected\n")
  cleanCountFormat<- T
  row.names(df) <- df[, 1]
  df[, 2:ncol(df)] <- sapply(df[, 2:ncol(df)], as.integer)
  cts <- as.matrix(df[, 2:ncol(df)])
} else if (colnames(df)[1] == "Gene" & colnames(df)[2] == "Name" & colnames(df)[3] == "Type") {
  cat("OSR COUNT.xlsx format detected\n")
  osrCountFormat <- T
  # clean names
  colnames(df) <- gsub(".COUNT", "", colnames(df))
  # convert to int
  df[, 4:ncol(df)] <- sapply(df[, 4:ncol(df)], as.integer)
  row.names(df) <- df$Gene
  # get cts
  cts <- as.matrix(df[, 4:ncol(df)])
} else {
  stop("COUNT file not in featureCounts output format, nor cleanCountFormat format, nor osrCount format. \nMay have different number of samples in meta-data and count-table if cleanCountFormat is provided\n")
}

# remove blackSamples
df <- df[, !(colnames(df) %in% blackSamples)]
cts <- cts[, !(colnames(cts) %in% blackSamples)]

# check number of samples
if ( dim(cts)[2] != dim(meta.df)[1] ){
  paste("# samples in count-table:", dim(cts)[2] )
  paste("# samples in meta-data provided:", dim(meta.df)[1] )
  stop("the number of samples in count-table does not match the number of samples in meta-data")
}

# check sample names
if (all(colnames(cts) %in% meta.df$SAMPLE_LABEL) & all(meta.df$SAMPLE_LABEL %in% colnames(cts)) ){
  cat("all sample-labels in count-table and meta-data matches")
} else {
  warning("count-table: ", colnames(cts), "\n")
  warning("meta-data: ", meta.df$SAMPLE_LABEL, "\n")
  stop("not all sample-labels in count-table and meta-data matches, abort!!!", colnames(cts), meta.df$SAMPLE_LABEL)
}



# Report Summary
print(paste("lib-size in millions:"))
print(format(colSums(cts/1e6), digits=2))
print(paste("Dim of input data:"))
print(dim(cts))
```


## Importing annotation
```{r import_annotation}
# from GitHub (must be raw, not zipped)
getAnnotation <- function(urlpath) {
  tmp <- tempfile()
  download.file(urlpath, destfile = tmp, method = 'auto')
  return(read.table(tmp, sep="\t", header = TRUE, quote=""))
}

paste("Getting data from", annoFile)

if (grepl('https://', annoFile)){
  print("downloading remote annoFile")
  anno <- getAnnotation(annoFile)
}else{
  print("reading local annoFile")
  anno <- read.table(annoFile, sep = "\t", header = T, quote="")
}
anno <- anno[!duplicated(anno[1]), ] # assume first column is Gene ID
print("Dimention of annotation table: ")
dim(anno)
head(anno)

if (sum(anno[, 1] %in% row.names(cts)) < 1){
  warning("!!! Annotation file and count file have no ID in common")
  warning("The results will be unannotated")
  warning("Please Double check Annotation file")
  print("count table ID:")
  print(row.names(cts)[1:2])
  print("anno table ID:")
  print(anno[1:2, 1])
  }
```

## COUNT output (without filtering)
```{r}
count <- data.frame(cts)
colnames(count) <- paste(colnames(count),"COUNT", sep = ":")
count_out <- merge(anno, count, by.x=1, by.y=0, all.y=T, sort=F)
head(count_out)
WriteXLS(x = count_out, 
         ExcelFileName = 'COUNT.xlsx', row.names = F, SheetNames = 'sheet1', na = 'NA')
print("saved in COUNT.xlsx")
```

## Filtering
```{r filtering}
min_rowsum <- 10
expression_filter <- rowSums(cts) >= min_rowsum  # default 10
min_colsum <- 2e4
sample_filter <- colSums(cts) >= min_colsum
if (sum(expression_filter) > 2 & sum(sample_filter) >= 4){
  cts <- cts[expression_filter, sample_filter]
  if (featureCountsFormat) {
    df <- df[expression_filter, c(rep(TRUE,5), sample_filter )]
  }else if (cleanCountFormat) {
    df <- df[expression_filter, sample_filter]
  } else if (osrCountFormat) {
    df <- df[expression_filter, c(rep(TRUE,3), sample_filter )]
  }  # order of cts and df kept the same
  cat(paste("Removed genes with less than ", min_rowsum, "reads/fragments across all samples", "\n"))
  cat(paste("Removed samples with less than ", min_colsum, "reads/fragments across all genes", "\n"))
  
  cat(paste("Data dim after filtering:"))
  dim(cts)
}else{
  print("Library size too small: too few genes/samples would be left after filtering, so skipped filtering.")
  print("Please interpret results with caution")
}
# print("Head of filtered data:")
# head(cts)
boxplot(log10(cts+1), las=2, main = "library size after filtering")
```

## COUNT calculation again (after filtering)
```{r}
count <- data.frame(cts)
colnames(count) <- paste(colnames(count),"COUNT", sep = ":")
count_out <- merge(anno, count, by.x=1, by.y=0, all.y=T, sort=F)
head(count_out)
```


## TPM calculation
```{r tpm}
calculateTPM <- function(counts,len) {
  # michael's version
  # https://support.bioconductor.org/p/91218/
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

if (featureCountsFormat ){
  tpm <- calculateTPM(cts, df$Length)
  tpm <- data.frame(tpm)
  colnames(tpm) <- paste(colnames(tpm),"TPM",  sep = ":")
  tpm_out <- merge(anno, tpm, by.x=1, by.y=0, all.y=T, sort=F)
  # head(tpm_out)
  # tail(tpm_out)
  WriteXLS(x = tpm_out, 
           ExcelFileName = 'TPM.xlsx', row.names = F, SheetNames = 'sheet1', na = 'NA')
  print("saved in TPM.xlsx")
} else {
  cat ("TPM not calculated if COUNT-table is not featureCount format\n")
}
```

## FPKM calculation
```{r fpkm}
calculateFPKM <- function(counts,len) {
  x <- counts
  x <- t(t(x)*1e6/colSums(x))
  return (x/len*1e3)
}

if (featureCountsFormat ){
  fpkm <- calculateFPKM(cts, df$Length)
  fpkm <- data.frame(fpkm)
  colnames(fpkm) <- paste(colnames(fpkm), "FPKM", sep = ":")
  fpkm_out <- merge(anno, fpkm, by.x=1, by.y=0, all.y=T, sort=F)
  # head(fpkm_out)
  # tail(fpkm_out)
  WriteXLS(x = fpkm_out, 
           ExcelFileName = 'FPKM.xlsx', row.names = F, SheetNames = 'sheet1', na = 'NA')
  print("saved in FPKM.xlsx")
  print("Recommend to use TPM, rather than FPKM")
} else {
  cat ("FPKM not calculated if COUNT-table is not featureCount format\n")
}
```

## Design matrix extracted from meta-data
```{r}
meta.df <- readExcel(metaFile)
meta.df <- meta.df[match(colnames(cts), meta.df$SAMPLE_LABEL), ]  # # todo: what if meta has less/more?

contrast.df <- readExcel(contrastFile)

sample <- factor(meta.df$SAMPLE_LABEL)
batch <- factor(meta.df$BATCH)
group <- factor(meta.df$GROUP_LABEL)

coldata <- data.frame(row.names=colnames(cts), 
                      sample,
                      group,
                      batch
                      )
coldata
```

## Model fitting
```{r}
if (length(levels(batch)) > 1){
  dds <- DESeqDataSetFromMatrix(countData = cts, 
                                colData = coldata, 
                                design = ~  0 + group + batch)
}else{
  dds <- DESeqDataSetFromMatrix(countData = cts, 
                                colData = coldata, 
                                design = ~  0 + group)  # converted to alph-order
}
dds
dds <-DESeq(dds)
resultsNames(dds)
#saveRDS(dds, file = 'deseq2.dds.rds')
```

## Save DESeq2 normalized counts
```{r deseq2_norm}
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- data.frame(normalized_counts)
colnames(normalized_counts) <- paste(colnames(normalized_counts),"DESeq2NormalizedCount",  sep = ":")
normalized_counts_out <- merge(anno, normalized_counts, by.x=1, by.y=0, all.y=T, sort=F)
WriteXLS(x = normalized_counts_out, 
         ExcelFileName = 'DESeq2NormalizedCounts.xlsx', row.names = F, SheetNames = 'sheet1', na = 'NA')
print("saved in DESeq2NormalizedCounts.xlsx")
```

## QC Plots

<!-- ## Data transformation -->
<!-- ```{r} -->
<!-- #vsd <- vst(dds, blind=FALSE) -->
<!-- rld <- rlog(dds, blind=FALSE) -->
<!-- counts <- counts(dds, normalized=0) -->
<!-- logCounts <- log10(counts +1 ) -->

<!-- normed <- counts(dds, normalized=1) -->
<!-- logNormed <- log10(normed+1) -->
<!-- ``` -->

### Histogram of Log10(Counts)
```{r histogram}
log1p_count <- log10(counts(dds) + 1)
hist(log1p_count, 
     main = 'Histogram of log10(count + 1)', 
     xlab = "log10(count+1)",
     100) # by default, use non-normalized data by counts function
```

### Dispersion plot
```{r dispersion_plot}
plotDispEsts(dds, main="Dispersion plot")
```


### Sample PCA plot
```{r pca}
plotQC_PCA <- function(dds) {
  vsd <- varianceStabilizingTransformation(dds) # fixed num < num(rowsum>5)
  
  pcaData <- plotPCA(vsd, intgroup = 'group', returnData=TRUE) # labeling fixed
  percentVar <- round(100 * attr(pcaData, 'percentVar'), 1)
  
  if (length(levels(batch)) > 1){
    ggplot(pcaData, aes(PC1, PC2, color = group, shape = batch)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    geom_label_repel(aes(label = sample),
                      box.padding = 0.35,
                      point.padding = 1,
                      segment.color = 'grey50',
                     segment.alpha = 0.5,
                      show.legend = FALSE) + # if TRUE, legend display might not be correct
    theme_classic()
  }else{
    ggplot(pcaData, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    geom_label_repel(aes(label = sample),
                      box.padding = 0.35,
                      point.padding = 1,
                      segment.color = 'grey50',
                     segment.alpha = 0.5,
                      show.legend = FALSE) + # if TRUE, legend display might not be correct
    theme_classic()
  }
}

plotQC_PCA_no_label <- function(dds) {
  vsd <- varianceStabilizingTransformation(dds) # fixed num < num(rowsum>5)
  
  
  pcaData <- plotPCA(vsd, intgroup = 'group', returnData=TRUE) # labeling fixed
  percentVar <- round(100 * attr(pcaData, 'percentVar'), 1)
  
    if (length(levels(batch)) > 1){
        ggplot(pcaData, aes(PC1, PC2, color = group, shape = batch)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_classic()
    }else{
        ggplot(pcaData, aes(PC1, PC2, color = group)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_classic()
    }
  
}

pdf("sample_PCA.labeled.pdf")
plotQC_PCA(dds)
dev.off()

pdf("sample_PCA.pdf")
plotQC_PCA_no_label(dds)
dev.off()

plotQC_PCA(dds)
plotQC_PCA_no_label(dds)
```






### Sample heatmap (Poisson Distance Based)
```{r poisson_heatmap}
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd ) 
rownames(samplePoisDistMatrix) <- coldata$sample
colnames(samplePoisDistMatrix) <- NULL 
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors, 
         clustering_method='complete', 
         legend = T, 
         filename = "sample_poisson_distance.pdf")

pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors, 
         clustering_method='complete',
         legend = T)
```


```{r, include=F}
### Sample Heatmap (Spearman Corr Based)
# corr <- cor(counts(dds), method="spearman")
# spearman_corrplot <- function(corr){
#   #corr <- round(corr, 2)
#   corrplot(corr, 
#            method="square", type="upper", add=F,  # add-on to previous plot
#            is.corr=T,  # flexible coloring if False
#            order="hclust", hclust.method="complete",
#            tl.pos="lt", diag = T)
#   
#   if (dim(corr)[1]<12){
#     corrplot(corr, 
#              method="number", type="lower", add=T,  # add-on to previous plot
#              is.corr=T,  # flexible coloring if False
#              order="hclust", hclust.method="complete",
#              tl.pos="n", diag = F, cl.pos = "n")
#   }else{
#     corrplot(corr, 
#              method="pie", type="lower", add=T,  # add-on to previous plot
#              is.corr=T,  # flexible coloring if False
#              order="hclust", hclust.method="complete",
#              tl.pos="n", diag = F, cl.pos = "n")
#   }
# }

# pdf(paste("sample_spearman_corr", "pdf", sep="."))
# spearman_corrplot(corr)
# dev.off()
# 
# spearman_corrplot(corr)
```

## DE analysis
```{r, warning=T, de}
gsub("group", "", resultsNames(dds))

for (i in 1:dim(contrast.df)[2]){
  name1 <- contrast.df[1,i]
  name2 <- contrast.df[2,i]
  name1 <- gsub(" ", "", name1)
  name2 <- gsub(" ", "", name2)
  name1 <- gsub(";$", "", name1)
  name2 <- gsub(";$", "", name2)
  name1s <- strsplit(name1, ";") [[1]]
  name2s <- strsplit(name2, ";") [[1]]
  name1 <- gsub(";", ".", name1)
  name2 <- gsub(";", ".", name2)
  
  name <- paste(name1, name2, sep = "_vs_")
  if (nchar(name) > 100) {name = paste0('contrast', i)}
  print(paste(">>>", i, name))
  
  # contrast <- c("group", name1, name2)
  
  poss <- match(name1s, gsub("group", "", resultsNames(dds)))
  negs <- match(name2s, gsub("group", "", resultsNames(dds)))
  contrast <- rep(0, length(resultsNames(dds)))
  contrast[poss] <- 1/length(poss)
  contrast[negs] <- -1/length(negs)
  print(data.frame(resNames=gsub("group", "", resultsNames(dds)), 
                   contrast=contrast))
  
  res <- lfcShrink(dds, contrast = contrast, type = 'ashr')
  res2 <- results(dds, contrast = contrast, 
                  independentFilter=indfilter, cooksCutoff=cookscutoff)
  if (featureCountsFormat){
    process_deseq_res(res = res, res2=res2, name = name, anno = anno, norm_exp = tpm)
  } else if (cleanCountFormat | osrCountFormat){
      process_deseq_res(res = res, res2=res2, name = name, anno = anno, norm_exp = normalized_counts)
  } else {
    stop("COUNT table format was wrong")
  }
} 
```



## Log
```{r log}
#save.image(file="../../DESeq2/RSession")
sessionInfo()
```

\name{build.plotter}
\docType{methods}
\alias{build.plotter}
\alias{build.plotter}
\alias{build.plotter.highlightSample.colorByLane}
\alias{build.plotter.highlightSample}
\alias{build.plotter.colorByLane}
\alias{build.plotter.colorByGroup}
\alias{build.plotter.colorBySample}
\alias{build.plotter.basic}
\alias{build.plotter.colorByX}
\alias{build.plotter.advanced}
\alias{QoRTs.default.plotting.params}
\alias{plotter}
\title{
   Generating plotters
}
\description{
   Generating QC_Plotter objects, which can be used in many of the QoRT utilities to organize samples in various ways to allow for easy comparison and detection of consistent biases and artifacts. 
}
\usage{
  build.plotter.basic(res, plotter.params = list())

  build.plotter.colorByGroup(res, plotter.params = list())
  
  build.plotter.colorByLane(res, plotter.params = list())
  
  build.plotter.colorBySample(res, plotter.params = list())
  
  build.plotter.highlightSample(curr.sample,
                                  res,
                                  plotter.params = list(),
                                  merge.offset.outgroup = TRUE)
                                  
  build.plotter.highlightSample.colorByLane(curr.sample,
                                  res,
                                  plotter.params = list(),
                                  merge.offset.outgroup = TRUE,
                                  lane.column.name = "lane.ID")

  build.plotter.colorByX(res, color.by.name,
                         color.by.title.name = color.by.name, 
                         plotter.params = list())

  build.plotter.advanced(res, 
                         colorBy = NULL,
                         color.title = "?",
                         highlightBy = NULL, 
                         highlight = "CURR", 
                         highlightTitle.singular = NULL,
                         highlightTitle.plural = highlightTitle.singular,
                         outgroup.title = "Other",
                         plotter.params = list())

}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  

  
  \item{curr.sample}{
    A character string. For the sample highlight summary plots,
    this should be the sample.ID of the sample that is to be
    highlighted.
  }
  \item{merge.offset.outgroup}{
    (For advanced users) A logical value. For the sample highlight plots, determines whether the all lanes that do not include the current sample should be treated as a single "outgroup".
  }
  \item{plotter.params}{
    (For advanced users) A named list. Allows you to specify 
    colors, offsets, and other similar patterns. By default 
    these will all be set to reasonable values, however, if 
    you want more control over colors, line-transparency, point 
    plotting characters, or similar, then you can specify a 
    named list.
    
    Any parameters that are not specified in the 
    \code{plotter.params} list will be left as default.
    
    Legal parameters are:
      \itemize{
        \item \code{"contrasting.colors"}: colors to use for contrast. By default these are set to a series of reasonably-contrasting colors. However, if you have too many different categories then it may be hard to tell some colors apart.
        \item \code{"contrasting.pch"}: point types to use for contrast (see \code{pch} in \link{graphical parameters}). By default this is set to the basic point types, then following through the upper and lower case letters.

        \item \code{"std.color"}: Color to use for the "highlighted" replicates.
        \item \code{"std.lines.lty"}: Line type to use for the "highlighted" replicates. (see \code{lty} in \link{graphical parameters})
        \item \code{"std.lines.lwd"}: Line width to use for the "highlighted" replicates. (see \code{lwd} in \link{graphical parameters})
        \item \code{"std.lines.alpha"}: Alpha transparency value to use on lines for the "highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"std.points.pch"}: Character to use for points for the "highlighted" replicates. (see \code{pch} in \link{graphical parameters})
        \item \code{"std.points.alpha"}: Alpha transparency value to use on points for the "highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"std.points.color"}: Color to use for the "highlighted" replicates.
        \item \code{"std.NVC.colors"}: A named list with elements named "A", "T", "C", and "G", with each element specifying a color. The colors used to indicate each base for the "highlighted" replicates in the nucleotide-rate-by-position plots.

        \item \code{"alt.color"}: Color to use for the "non-highlighted" replicates.
        \item \code{"alt.lines.lty"}: Line type to use for the "non-highlighted" replicates. (see \code{lty} in \link{graphical parameters})
        \item \code{"alt.lines.lwd"}: Line width to use for the "non-highlighted" replicates. (see \code{lwd} in \link{graphical parameters})
        \item \code{"alt.lines.alpha"}: Alpha transparency value to use on lines for the "non-highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"alt.points.pch"}: Character to use for points for the "non-highlighted" replicates. (see \code{pch} in \link{graphical parameters})
        \item \code{"alt.points.alpha"}: Alpha transparency value to use on points for the "non-highlighted" replicates. Numeric value between 0 and 255.
        \item \code{"alt.NVC.colors"}: A named list with elements named "A", "T","C", and "G", with each element specifying a color. The colors used to indicate each base for the "non-highlighted" replicates in the nucleotide-rate-by-position plots.
        
        \item \code{"show.legend"}: DEPRECIATED. Currently nonfunctional.

      }
  }
  \item{color.by.name}{
    (For advanced users) A character string. (TODO: document functionality)
  }
  \item{color.by.title.name}{
    (For advanced users) A character string. (TODO: document functionality)
  }
  \item{lane.column.name}{
    The name of the column in the decoder containing the "lane" names.
  }
  
  \item{colorBy}{
    A named character vector. Each unique colorBy string will be assigned a unique color. The names of colorBy must match \code{res@decoder$unique.ID}, and must be in the same order.
  }
  \item{color.title}{
    A character string. This is the title of the colorby category, used in the titles and figure legends.
  }
  \item{highlightBy}{
    A named character vector. Used to determine which replicates to highlight. The names of colorBy must match \code{res@decoder$unique.ID}, and must be in the same order.
  }
  \item{highlight}{
    A character string. Replicates where highlight equals highlightBy will be highlighted.
  }
  \item{highlightTitle.singular}{
    A character string. The singular form of the name of the category highlighted.
  }
  \item{highlightTitle.plural}{
    A character string. The plural form of the name of the category highlighted.
  }
  \item{outgroup.title}{
    A character string. The description of the non-highlighted category. Used in the figure legends.
  }
}
\value{
  A QoRT_Plotter reference object used to create QC summary plots. Depending on which plotter is used, samples/lane-bams can be organized by group, sample, lane, or any arbitrary variable found in the decoder.
}
\examples{
data(res,package="QoRTsExampleData");
plotter.basic <- build.plotter.basic(res);
makePlot.insert.size(plotter.basic);

plotter.colorByGroup <- build.plotter.colorByGroup(res);
makePlot.insert.size(plotter.colorByGroup);
makePlot.legend.over("topright",plotter.colorByGroup);

plotter.colorByLane <- build.plotter.colorByLane(res);
makePlot.insert.size(plotter.colorByLane);
makePlot.legend.over("topright",plotter.colorByLane);

plotter.colorBySample <- build.plotter.colorBySample(res);
makePlot.insert.size(plotter.colorBySample);
makePlot.legend.over("topright",plotter.colorBySample);

plotter.HS <- build.plotter.highlightSample("SAMP1",
                                            res);
makePlot.insert.size(plotter.HS);
makePlot.legend.over("topright",plotter.HS);

plotter.HSCBL <- build.plotter.highlightSample.colorByLane("SAMP1",
                                                           res);
makePlot.insert.size(plotter.HSCBL);
makePlot.legend.over("topright",plotter.HSCBL);


#FOR ADVANCED USERS:
#  With the build.plotter.advanced function, you can
#  set coloring and highlighting to match anything you
#  want.
#  The parameters are a little more complex...

#In order to control color, you must create a named
#  vector with names equal to the unique.ID's
#  in the decoder, and in the same order:
#  (this requirement is purely to prevent mistakes)

#For example: to color each sample differently:
colorBy <- res@decoder$sample.ID
names(colorBy) <- res@decoder$unique.ID;

plotter <- build.plotter.advanced(res, colorBy = colorBy);
makePlot.insert.size(plotter);
makePlot.legend.over("topright",plotter);

#Now, to highlight a subgroup of the dataset, you
#  must set the "highlightBy" parameter to a 
#  named vector with names equal to the decoder 
#  unique.ID's, and in the same order.
#  (this requirement is purely to prevent mistakes)
#Then you must tell the plotter which subgroup
#  you want to highlight using the "highlight" 
#  parameter.

#For example, to highlight all lanebams in lane L1:
highlightBy <- res@decoder$lane.ID
names(highlightBy) <- res@decoder$unique.ID;

plotter <- build.plotter.advanced(res, 
                                 highlightBy = highlightBy,
                                 highlight = "L1");
makePlot.insert.size(plotter);
makePlot.legend.over("topright",plotter);

#Other parameters are available to change the title
# and legends:
plotter <- build.plotter.advanced(res, 
                                 highlightBy = highlightBy,
                                 highlight = "L1",
                                 highlightTitle.singular = "Lane",
                                 highlightTitle.plural  = "Lanes",
                                 outgroup.title = "Other");
makePlot.insert.size(plotter);
makePlot.legend.over("topright",plotter);

#You can also color and highlight together.
#  If you do this, only the "highlighted" group will be
#  colored, all the others will be colored gray and will be 
#  drawn in the background. This can be useful for finding
#  biases that are restricted to a subset of the data.

plotter <- build.plotter.advanced(res, 
  colorBy = colorBy,
  highlightBy = highlightBy,
  highlight = "L1",
  color.title = "sample",
  highlightTitle.singular = "Lane",
  highlightTitle.plural = "Lanes",
  outgroup.title = "Other"
);
makePlot.insert.size(plotter);
makePlot.legend.over("topright",plotter);


#You can make multiplots using a given plotter object by
# using the "makeMultiPlot.withPlotter" function:

#makeMultiPlot.withPlotter(plotter);

}

\seealso{
  \code{\link{read.qc.results.data}} 
}
\name{makePlot.clipping}
\docType{methods}
\alias{makePlot.clipping}
\title{
   Plot Alignment Clipping
}
\description{
   Plots the rate at which the aligner soft-clips off portions of
   aligned reads.
}
\usage{
  makePlot.clipping(plotter, 
                 rate.per.million = FALSE,
                 use.readLength.denominator = TRUE, 
                 r2.buffer,
                 debugMode, 
                 singleEndMode, 
                 rasterize.plotting.area = FALSE, 
                 raster.height = 1000, 
                 raster.width = 1000, 
                 plot= TRUE, ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{rate.per.million}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
  }
  \item{use.readLength.denominator}{
    Logical. If TRUE, then use the total number of reads with at least the given position's length as the denominator of the rate. This is only relevant when the dataset has variable-length trimming prior to alignment.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: The read cycle (ie. the base-pair position in the read).
  
  y-axis: The rate at which the bases at the given read-cycle is
  clipped off.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.clipping(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{get.summary.table}
\docType{methods}
\alias{get.summary.table}
\alias{get.size.factors}

\title{
   Get summary data tables
}
\description{
   Retrieves and compiles a summary data table.
}
\usage{
  get.summary.table(res, outfile, debugMode);
  
  get.size.factors(res, 
       sf.method = c("DESeq2","DESeq2_GEO","TC",
                     "edgeR","edgeR_TMM","edgeR_UQ","edgeR_RLE"), 
       outfile, debugMode)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile}{
    Optional. A file name where the table should be written.
  }
  \item{sf.method}{
    The size factor method to use. Note that most of these methods (except "TC") are reliant
    on external packages not included with QoRTs. You will need to install DESeq2 and edgeR
    to use these methods.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
}
\details{
  Returns summary data in table form.
}
\examples{
  data(res,package="QoRTsExampleData");
  get.summary.table(res);
}

\seealso{
  \code{\link{read.qc.results.data}} 
}\name{makePlot.runTimePerformance}
\docType{methods}
\alias{makePlot.runTimePerformance}
\title{
   Chromosome type rate plot
}
\description{
   Plots the number or percent of read-pairs falling on each type of chromosome.
}
\usage{
  makePlot.runTimePerformance(plotter,
                 debugMode, singleEndMode,
                 plot = TRUE,
                 ...)

}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }

}
\details{
   For each sample run, indicates the amount of time spent running the QoRTs QC data processing tool
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.runTimePerformance(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}
\name{makePlot.readLengthDist}
\docType{methods}
\alias{makePlot.readLengthDist}
\title{
   Plot the distribution of read lengths.
}
\description{
   Plots the distribution of read lengths. Only useful for data with variable trimming (which is generally not recommended at least for RNA-Seq data).
}
\usage{
makePlot.readLengthDist(plotter,
                        plot.rates = TRUE, 
                        plot.means = TRUE, 
                        plot.medians = NULL,
                        include.full.length = FALSE, 
                        cumulative = TRUE,
                        singleEndMode,
                        rasterize.plotting.area = FALSE, 
                        raster.height = 1000, 
                        raster.width = 1000,
                        debugMode,
                        r2.buffer,
                        plot = TRUE,
                        ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether or not the X-axis should be 
    the raw number of nucleotides that are G/C, vs the rate G/C.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average GC content for each bam file at the bottom 
    of the plot.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average GC content for each bam file at the bottom 
    of the plot. Overrides \code{plot.means}.
  }
  \item{include.full.length}{
    Logical. If FALSE, omit the full-length read length from the x-axis of the plot.
  }
  \item{cumulative}{
    Logical. If TRUE, plot shows cumulative rates.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }

  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Read Length
  
  y-axis: Percentage of reads with length equal to the given length. If cumulative == TRUE, then it is the percentage of reads with length less than or equal to the given length.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.readLengthDist(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.strandedness.test}
\docType{methods}
\alias{makePlot.strandedness.test}
\title{
   Strandedness Test Plot
}
\description{
   Plots the apparent strandedness of the reads.
}
\usage{
makePlot.strandedness.test(plotter, plot.target.boxes = FALSE, 
                           debugMode, singleEndMode, plot = TRUE, \dots)

}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.target.boxes}{
    A logical value. If true, then green target boxes will be printed over the area in which all points should be expected to fall.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
    Logical. Determines whether the dataset consists of single-ended reads. 
    By default this is determined from the data. Thus, this parameter should 
    never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   TODO!
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.strandedness.test(plotter);
  
  #Add a legend:
  makePlot.legend.over("topright", plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.chrom.type.rates}
\docType{methods}
\alias{makePlot.chrom.type.rates}
\title{
   Chromosome type rate plot
}
\description{
   Plots the number or percent of read-pairs falling on each type of chromosome.
}
\usage{
makePlot.chrom.type.rates(plotter, 
                      plot.rates = TRUE,
                      chromosome.name.style = "UCSC",
                      exclude.autosomes = FALSE,
                      chrom.norm.factors = NULL,
                      custom.chromosome.style.def.function = NULL,
                      return.table = FALSE, 
                      debugMode, singleEndMode,
                      plot = TRUE,
                      ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether to plot percent of the total 
    (for each bam file), rather than read-counts.
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. There are 4 legal
    options:
    \itemize{
      \item "UCSC": The default. Chromosomes are named: "chr1,chr2,...,chrX,chrY,chrXY,chrM". There are 6 categories: autosome, X, Y, XY, mitochondrial, and other.
      \item "ENSEMBL": Chromosomes are named: "1,2,...,X,Y,XY,MT". There are 6 categories: autosome, X, Y, XY, mitochondrial, and other.
      \item "UCSC_WITH_ERCC": As UCSC, but there is an additional category, which contains all chromosomes that begin with the text "ERCC".
      \item "ENSEMBL_WITH_ERCC": As ENSEMBL, but there is an additional category, which contains all chromosomes that begin with the text "ERCC". 
    }
  }
  \item{chrom.norm.factors}{
    (Advanced users)
  }
  \item{exclude.autosomes}{
    A logical value indicating whether to exclude autosomes from the plot.
  }
  \item{custom.chromosome.style.def.function}{
    (For advanced users) If your chromosomes do not match any of the above styles, then you can set your own chromosome style by handing this option a function. The function must take one argument. When handed NULL, it must return a list of all legal categories. When handed a single chromosome name, it must return one of those categories. 
  }
  \item{return.table}{
    A Logical value. If TRUE, the function will return a table containing the plotted data.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
  
  
}
\details{
   For each sample-run, this function plots the number of read-pairs mapping to each category of chromosome. 
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.chrom.type.rates(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.dropped.rates}
\docType{methods}
\alias{makePlot.dropped.rates}
\title{
   Read Drop Plot
}
\description{
   Plots the rates at which reads are dropped from
   analysis for various causes.
}
\usage{
makePlot.dropped.rates(plotter, dropAlwaysZeroRows = FALSE, 
                       debugMode, singleEndMode,
                       plot = TRUE,
                       \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{dropAlwaysZeroRows}{
    Logical. If TRUE, drop-reasons that never occur in the dataset will not be plotted. This often cleans up the plot somewhat, since in many production pipelines reads that fail many of the filtering steps may have already been filtered out.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   For each bam file, this function plots the rates and reasons for reads being dropped from QC analysis. 
   
   Note that in the example dataset reads were never dropped. This is a consequence of the pre-processing steps in the example pipeline.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.dropped.rates(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.splice.junction.event.rates}
\docType{methods}
\alias{makePlot.splice.junction.event}
\alias{makePlot.splice.junction.event.ratesPerRead}
\alias{makePlot.splice.junction.event.counts}
\alias{makePlot.splice.junction.event.proportions}
\alias{makePlot.splice.junction.event.proportionsByType}
\title{
   Plot Splice Junction Event Rates
}
\description{
   Plots the rates at which splice junctions occur.
}
\usage{
   makePlot.splice.junction.event.counts(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, singleEndMode, 
                                    plot = TRUE,
                                    ...)
   makePlot.splice.junction.event.ratesPerRead(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, singleEndMode,
                                    plot = TRUE,
                                    ...)   
   makePlot.splice.junction.event.proportions(plotter,
                                    high.low.cutoff = 4, 
                                    debugMode, singleEndMode,
                                    plot = TRUE,
                                    ...) 
   makePlot.splice.junction.event.proportionsByType(plotter, 
                                    high.low.cutoff = 4, 
                                    debugMode, singleEndMode,
                                    plot = TRUE,
                                    ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{high.low.cutoff}{
    For Advanced users only! The cutoff between high and low expression splice junctions. Note that in order to function, this same cutoff MUST be used by the QoRTs jar utility that generates these counts.
  }
  \item{debugMode}{
    Activates debug mode, which causes more verbose reporting.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   These functions plot various metrics for the rate at which splice junction "events" occur. A splice junction "event" is an occurrance of a mapped read-pair bridging a splice junction. Some read-pairs may contain multiple splice junction events, and some read-pairs may contain none.
   
   Splice junctions are characterized into six categories:
      \itemize{
        \item Known: The splice junction locus is found in the supplied transcript annotation gtf file.
        \item Novel: The splice junction locus is NOT found in the supplied transcript annotation gtf file.
        \item Known, 1-3 reads: The locus is known, and is only covered by 1-3 read-pairs.
        \item Known, 4+ reads: The locus is known, and is covered by 4 or more read-pairs.
        \item Novel, 1-3 reads: The locus is novel, and is only covered by 1-3 read-pairs.
        \item Novel, 4+ reads: The locus is novel, and is covered by 4 or more read-pairs.
      }
    
   \code{makePlot.splice.junction.event.counts} plots the number (y-axis) of all splice junction events falling into each of six categories. Note that because different samples/runs may have different total read counts and/or library sizes, this function is generally not the best for comparing between samples. For most purposes, \code{makePlot.splice.junction.event.ratesPerRead} will be preferable.
   
   \code{makePlot.splice.junction.event.ratesPerRead} plots the rate (y-axis) at which each type of splice junction events appear, per read-pair. 
   
   \code{makePlot.splice.junction.event.proportions} plots the proportion (y-axis) of all splice junction events falling into the six categories.
   
   \code{makePlot.splice.junction.event.proportionsByType} plots the proportion (y-axis) at which splice junction events appear on known vs novel splice junction loci, the proportion of known splice junction events that occur on low-coverage junctions (1-3 read-pairs) vs high-coverage junctions (4 or more read-pairs), and the proportion of novel splice junction events that occur on low vs high coverage junctions.
   
   All of these plots are generally used to detect whether sample-specific or batch effects have a substantial or biased effect on splice junction appearance, either due to differences in the original RNA, or due to artifacts that alter the rate at which the aligner maps across splice junctions.
}
\examples{
data(res,package="QoRTsExampleData");
plotter <- build.plotter.colorByGroup(res);
makePlot.splice.junction.event.counts(plotter);  
makePlot.splice.junction.event.ratesPerRead(plotter);
makePlot.splice.junction.event.proportions(plotter);
makePlot.splice.junction.event.proportionsByType(plotter);
  
#Legend:
makePlot.legend.box(plotter);
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.splice.junction.loci.counts}}
}\name{makePlot.all.std}
\docType{methods}
\alias{makeMultiPlot.all}
\title{
   Generating all default plots
}
\description{
   Saves MANY compiled QC plots for the given dataset.
}
\usage{
makeMultiPlot.all(res, outfile.dir = "./",
                  plotter.params = list(), 
                  plot.device.name = "png", 
                  plotting.device.params = list(), 
                  debugMode, 
                  rasterize.large.plots, 
                  rasterize.medium.plots,
                  raster.height = 1050, 
                  raster.width = 1050,
                  exclude.autosomes.chrom.rate.plot = TRUE,
                  chromosome.name.style = "UCSC",
                  fig.res = 150, fig.base.height.inches = 7,
                  insertSize.plot.xlim ,
                  sequencing.type = c("RNA","Exome","Genome"),
                  maxColumns,
                  maxRows,
                  plotList,
                  labelPlots = TRUE,
                  ...)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile.dir}{
    A file prefix, used for all output files. Usually the directory to
    which you want all files to be written.
  }
  \item{plotter.params}{
    Additional (advanced) parameters used in creation of the Plotter 
    objects. See \link{build.plotter}. 
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Can be one of:
      \itemize{
        \item \code{"png"} for standard png compression,
        \item \code{"CairoPNG"} for png compression using package Cairo.
          Note that this requires the package Cairo.
      }
  }
  \item{plotting.device.params}{ 
    A named list of parameters to be passed to the graphics device.
    For example:
      \itemize{
        \item \code{"width = 2000"}
      }
    Reasonable values for height, width, and pointsize will be chosen
    by default.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.large.plots}{
    Logical. If TRUE, then if the currently selected plotting device
    is a vector device (or the "curr" device), then certain large plots
    will have their plotting areas rasterized. The axis labels, titles,
    and text will all remain vectorized, only the plotting areas will
    be flattened. \emph{Note that this requires the png package.}
    
    By default, this parameter will be set to TRUE when a vector
    device is selected.
  }
  \item{rasterize.medium.plots}{
    Logical. as rasterize.large.plots, but applies to moderately-large plots. By default, this parameter will be set to TRUE for pdf/CairoPDF output only.
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels. Double-pane plots
    will be twice this width.
  }
  \item{exclude.autosomes.chrom.rate.plot}{
    A logical value indicating whether to exclude autosomes from the plot.
    See \code{\link{makePlot.chrom.type.rates}}
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. See 
    \code{\link{makePlot.chrom.type.rates}}
  }
  \item{fig.res}{
    Numeric value. The number of pixels per "inch" (for raster devices only).
    For some plotting devices the figure height will be in pixels not inches, and
    will equal this value multiplied by the fig.base.height.inches value.
  }
  \item{fig.base.height.inches}{
    Numeric value. The base height, in inches, of each sub-figure in the plot. This
    will be equal to the height for vector devices, or used to calculate the height
    in pixels using the fig.res value (see above).
  }
  \item{insertSize.plot.xlim}{
    A numeric vector of length 2. The x-axis limits for the insert size plot. By default QoRTs will attempt to find reasonable values for this, 
    but there are always situations where the default behavior is not ideal. Using this parameter you can set
    it explicitly.
  }
  \item{sequencing.type}{
    The type of sequencing data being analyzed. This only changes the default plot set, which can be overriden with the \code{plotList} parameter.
  }
  
  \item{maxColumns}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxColumns columns. Extra rows will be added to fit all the plots, as needed.
  }
  \item{maxRows}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxRows rows. Extra columns will be added to fit all the plots as needed. To set the number of rows and columns manually, you can set both maxColumns and maxRows.
  }
  \item{plotList}{
    The list of desired plots.
  }
  \item{labelPlots}{
    Logical. If TRUE, then label each plot with a letter.
  }
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  Saves MANY compiled QC plots for the given dataset, calling all 
  standard variants of \code{\link{makeMultiPlot}}.
}
\examples{
  \dontrun{
  data(res,package="QoRTsExampleData");
    makeMultiPlot.all(res, outfile.dir = "./");
  }
}
\seealso{
  \code{\link{read.qc.results.data}}, \code{\link{build.plotter}}, \code{\link{makeMultiPlot}}
}\name{makePlot.biotype.rates}
\docType{methods}
\alias{makePlot.biotype.rates}
\title{
   Plot Biotype rates
}
\description{
   Plots counts for each gene biotype. This plot is only useful and informative when QoRTs is run on a GTF file that contains "gene_biotype" tags.
}
\usage{
  makePlot.biotype.rates(plotter, 
              plot.rates = TRUE,
              count.type = c("all","unambigOnly"),
              log.y = TRUE,
              return.table = FALSE, 
              debugMode = DEFAULTDEBUGMODE,  
              singleEndMode = plotter$res@singleEnd,
              showTypes,
              plot = TRUE,
              ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value. If \code{TRUE}, then reads/read-pairs 
    per million should be plotted, rather than raw counts. 
    Default is \code{TRUE}.
  }
  \item{count.type}{
    A logical value. If \code{TRUE}, then both ambiguous and unambiguous
    reads will be counted. Otherwise only unambiguous reads will be counted.
  }
  \item{log.y}{
    A logical value indicating that the y-axis should be log-scaled.
  }
  \item{return.table}{
    Logical. If \code{TRUE}, then return a data table.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{showTypes}{
      Character vector. An optional list of biotypes to include in the plot. By
      default all observed biotypes will be plotted.
  }
  \item{plot}{
      Logical. Default is \code{TRUE}. If \code{FALSE}, then do NOT plot any results, 
      instead just return a data table.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Gene "BioTypes", from the annotation GTF.
  
  y-axis: Total read counts across all genes in the BioType, 
  measured either in reads or reads-per-million.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.biotype.rates(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makeMultiPlot.highlightSample.all}
\docType{methods}
\alias{makeMultiPlot.highlightSample.all}
\alias{makeMultiPlot.highlightSample.colorByLane.all}
\title{
   Generating sample highlight multi-plots
}
\description{
   Generates multiple sample highlight summary plots, 
   one for every sample.
}
\usage{
  makeMultiPlot.highlightSample.all(res, outfile.dir = "./",
                plotter.params = list(), 
                plot.device.name = "png", 
                plotting.device.params = list(), 
                verbose = TRUE,
                debugMode , 
                rasterize.large.plots, 
                rasterize.medium.plots,
                raster.height = 1050, 
                raster.width = 1050,
                exclude.autosomes.chrom.rate.plot = TRUE,
                chromosome.name.style = "UCSC",
                fig.res = 150, fig.base.height.inches = 7,
                insertSize.plot.xlim,
                sequencing.type = c("RNA","Exome","Genome"),
                maxColumns,
                maxRows,
                plotList,
                labelPlots=TRUE,
                ...)

makeMultiPlot.highlightSample.colorByLane.all(res, 
                outfile.dir = "./",
                plotter.params = list(), 
                plot.device.name = "png", 
                plotting.device.params = list(), 
                verbose = TRUE, 
                debugMode ,
                rasterize.large.plots, 
                rasterize.medium.plots,
                raster.height = 1050, 
                raster.width = 1050,
                exclude.autosomes.chrom.rate.plot = TRUE,
                chromosome.name.style = "UCSC",
                fig.res = 150, fig.base.height.inches = 7,
                insertSize.plot.xlim,
                sequencing.type = c("RNA","Exome","Genome"),
                maxColumns,
                maxRows,
                plotList,
                labelPlots=TRUE,
                ...)
}
\arguments{
  \item{res}{
    A \code{QoRT_QC_Results} object, created by \code{\link{read.qc.results.data}}.
  }
  \item{outfile.dir}{
    A file prefix, used for all output files. Usually the directory to
    which you want all files to be written.
  }
  \item{plotter.params}{
    Additional parameters (advanced) used in creation of the Plotter 
    objects. See \link{build.plotter}. 
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Can be one of:
      \itemize{
        \item \code{"png"} for standard png compression,
        \item \code{"CairoPNG"} for png compression using package Cairo.
          Note that this requires the package Cairo.
      }
  }
  \item{plotting.device.params}{ 
    A named list of parameters to be passed to the graphics device.
    For example:
      \itemize{
        \item \code{"width = 2000"}
      }
    Reasonable values for height, width, and pointsize will be chosen
    by default.
  }
  \item{verbose}{
    Logical. If TRUE, more info will be printed to the console.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.large.plots}{
    Logical. If TRUE, then if the currently selected plotting device
    is a vector device (or the "curr" device), then certain large plots
    will have their plotting areas rasterized. The axis labels, titles,
    and text will all remain vectorized, only the plotting areas will
    be flattened. \emph{Note that this requires the png package.}
    
    By default, this parameter will be set to TRUE when a vector
    device is selected.
  }
  \item{rasterize.medium.plots}{
    Logical. as rasterize.large.plots, but applies to moderately-large plots. By default, this parameter will be set to TRUE for pdf/CairoPDF output only.
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels. Double-pane plots
    will be twice this width.
  }
  \item{exclude.autosomes.chrom.rate.plot}{
    A logical value indicating whether to exclude autosomes from the plot.
    See \code{\link{makePlot.chrom.type.rates}}
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. See 
    \code{\link{makePlot.chrom.type.rates}}
  }
  \item{fig.res}{
    Numeric value. The number of pixels per "inch" (for raster devices only).
    For some plotting devices the figure height will be in pixels not inches, and
    will equal this value multiplied by the fig.base.height.inches value.
  }
  \item{fig.base.height.inches}{
    Numeric value. The base height, in inches, of each sub-figure in the plot. This
    will be equal to the height for vector devices, or used to calculate the height
    in pixels using the fig.res value (see above).
  }
  \item{insertSize.plot.xlim}{
    A numeric vector of length 2. The x-axis limits for the insert size plot. By default QoRTs will attempt to find reasonable values for this, 
    but there are always situations where the default behavior is not ideal. Using this parameter you can set
    it explicitly.
  }
  \item{sequencing.type}{
      The type of sequencing data being analyzed. This only changes the default plot set, which can be overriden with the \code{plotList} parameter.
  }
  \item{maxColumns}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxColumns columns. Extra rows will be added to fit all the plots, as needed.
  }
  \item{maxRows}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxRows rows. Extra columns will be added to fit all the plots as needed. To set the number of rows and columns manually, you can set both maxColumns and maxRows.
  }
  \item{plotList}{
    The list of desired plots.
  }
  \item{labelPlots}{
    Logical. If TRUE, then label each plot with a letter.
  }
  
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  Generates a sample highlight plot for each sample in the dataset.
}

\examples{
  \dontrun{
  data(res,package="QoRTsExampleData");
    makeMultiPlot.highlightSample.all(res);
    makeMultiPlot.highlightSample.colorByLane.all(res);
  }
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makeMultiPlot}}
}\name{makePlot.qual.pair}
\docType{methods}
\alias{makePlot.qual.pair}
\title{
   Plot quality score by read cycle
}
\description{
   Plots the Phred quality score as a function of the read cycle
   for both reads.
}
\usage{
  makePlot.qual.pair(plotter, y.name, r2.buffer = NULL, 
                 debugMode, singleEndMode, 
                 rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                 plot = TRUE,
                 ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{y.name}{
    The name of the quality score metric to plot. Must be one of:
      \itemize{
        \item \code{"min"}
        \item \code{"lowerQuartile"} 
        \item \code{"median"} 
        \item \code{"upperQuartile"} 
        \item \code{"max"}
      }
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
These plots display information about the phred quality score 
(y-axis) as a function of the position in the read (x-axis). Five 
statistics can be plotted: minimum, maximum, upper and lower 
quartiles, and median. These statistics are calculated individually
for each bam file and each read position (ie, each plotted line 
corresponds to a bam file).

Note that the Phred score is always an integer, and as such these 
plots would normally be very difficult to read because lines 
would be plotted directly on top of one another. To reduce this 
problem, the lines are vertically offset from one another. Most
plotters offset each line by lane.ID.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.qual.pair(plotter,y.name="min");
  makePlot.qual.pair(plotter,y.name="median");
  makePlot.qual.pair(plotter,y.name="max");
  
  #Legend:
  makePlot.legend.box(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.NVC.lead.clip}
\docType{methods}
\alias{makePlot.NVC.lead.clip}
\alias{makePlot.NVC.tail.clip}
\title{
   Clipped sequence nucleotide-by-position plot
}
\description{
   Plots the nucleotide rates for clipped segments
}
\usage{
makePlot.NVC.lead.clip(plotter, clip.amt,  r2.buffer, 
                   points.highlighted = TRUE,
                   label.majority.bases = TRUE, 
                   label.majority.bases.threshold = 0.5, 
                   label.majority.bases.cex = 1, 
                   rasterize.plotting.area = FALSE, 
                   raster.height = 1000, 
                   raster.width = 1000,
                   show.base.legend = TRUE, 
                   debugMode, singleEndMode,
                   plot = TRUE,
                   \dots)


makePlot.NVC.tail.clip(plotter, clip.amt,  r2.buffer, 
                       points.highlighted = TRUE, 
                       label.majority.bases = TRUE, 
                       label.majority.bases.threshold = 0.5, 
                       label.majority.bases.cex = 1, 
                       rasterize.plotting.area = FALSE, 
                       raster.height = 1000, 
                       raster.width = 1000,
                       show.base.legend = TRUE, 
                       debugMode, singleEndMode,
                       plot = TRUE,
                       \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{clip.amt}{
    Numeric value. The number of bases clipped. These methods
    only plot the sequence for a single specific clip.amt. 
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  \item{points.highlighted}{
    A logical value. Determines whether to ever plot points on top of the lines. This
    can be useful for identifying outliers. If TRUE, then all highlighted lane-bams
    will be overlaid with points using their designated pch symbol. If the plotter does
    not highlight any lanebams, then points will be overlaid on ALL lines.
  }
  \item{label.majority.bases}{
    A logical value indicating whether to label the genotypes of
    read cycles in which the most common base has a frequency 
    exceeding label.majority.bases.threshold (see below).
  }
  \item{label.majority.bases.threshold}{
    A numeric value indicating the cutoff above which the most 
    common base will be labelled, if label.majority.bases is set
    TRUE (see above).
  }
  \item{label.majority.bases.cex}{
    The cex value fed to text() that is used to determine how
    big the labels are to be, if label.majority.bases is TRUE.
    (see \code{\link{par}} for information on cex).
  }
  \item{rasterize.plotting.area}{
    Logical. If TRUE, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  \item{show.base.legend}{
    A logical value indicating whether to print the base-color
    legend along the right edge of the plot.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   For general information on reading NVC plots, see \code{\link{makePlot.NVC}}.
   
   This function plots the nucleotide rates for the sections of reads that were soft-clipped by the aligner. Note that these will not function on reads that have been aligned using an aligner that does not practice soft clipping (such as TopHat2).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.NVC.lead.clip(plotter, clip.amt = 12);
  makePlot.NVC.tail.clip(plotter, clip.amt = 12);
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.NVC}}
}\name{read.qc.results.data}
\docType{methods}
\alias{read.qc.results.data}
\alias{completeAndCheckDecoder}
\title{
   Reading QC results data
}
\description{
   Creates a QoRT_QC_Results object using a set of QC result 
   data files.
}
\usage{
  read.qc.results.data(infile.dir, 
                       decoder, 
                       decoder.files, 
                       calc.DESeq2 = FALSE, 
                       calc.edgeR = FALSE, 
                       debugMode,
                       autodetectMissingSamples = FALSE,
                       skip.files = c())
                       
  completeAndCheckDecoder(decoder, decoder.files)
}
\arguments{
  \item{infile.dir}{
    REQUIRED. The base file directory where all the QC results data
    is stored.
  }
  \item{decoder}{
    A character vector or data.frame containing the decoder information. See details below.
  }
  \item{decoder.files}{
    Character vector. Either one or two character strings. Either decoder.files OR decoder must be set, never both. See details below.
  }
  \item{calc.DESeq2}{
    Logical. If TRUE, this function will attempt to load the DESeq2 package. If the DESeq2 package is found, it will then calculate DESeq2's geometric normalization factors (also known as "size factors") for each replicate and for each sample.
  }
  \item{calc.edgeR}{
    Logical. If TRUE, this function will attempt to load the edgeR package. If the edgeR package is found, it will then calculate all of edgeR's  normalization factors (also known as "size factors") for each replicate and for each sample.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{autodetectMissingSamples}{
    Logical. If TRUE, automatically drop replicates for which the QC files cannot be found. By default this function will throw an error.
  }
  \item{skip.files}{
    Character vector of QC data file titles to skip. This can be useful for performing a faster data load when only querying a subset of the available QC metrics. See examples below.
  }
}
\details{
  read.qc.results.data reads in a full QoRTs dataset of multiple QoRTs QC runs and compiles them into a QoRTs_Results object.
  
  completeAndCheckDecoder simply reads a decoder and "fills in" all missing parameters, returning a data.frame.
  
  The "decoder" is used to describe each replicate/sample. The standard decoder is a data frame that has one row per replicate, with the following columns:
  \itemize{
    \item unique.ID: The base identifier for the individual replicate. Must be unique. lanebam.ID is a synonym.
    \item lane.ID: (OPTIONAL) The identifier for the lane, run, or batch. The default is "UNKNOWN".
    \item group.ID: (OPTIONAL) The identifier for the biological condition for the given replicate. The default is "UNKNOWN".
    \item sample.ID: (OPTIONAL) The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis. By default, it is assumed that each replicate comes from a different sample, and sample.ID is set to equal unique.ID.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the unique.ID.
    \item input.read.pair.count: (OPTIONAL) This column contains the number of read-pairs (or just reads, for single-end data), before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of read-pairs (or just reads, for single-end data) that were multi-mapped. This must be included for multi-mapping rate to be calculated.
  }
  All the parameters except for unique.ID are optional. The decoder can even be supplied as a simple character vector, which is assumed to be the unique.ID's. All the other variables will be set to their default values.
  
  Alternatively, the decoder can be supplied as a file given by the decoder.files parameter.
  
  Dual Decoder: Optionally, two decoders can be supplied. In this case the first decoder should be the technical-replicate decoder and the second should be the biological-replicate decoder.
  The technical-replicate decoder should have one row per unique.ID, with the following columns:
  \itemize{
    \item unique.ID: The base identifier for the individual replicate. Must be unique!
    \item lane.ID: (OPTIONAL) The identifier for the lane, run, or batch.
    \item sample.ID: (OPTIONAL) The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "lanebam.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item qc.data.dir: (OPTIONAL) This column indicates the subdirectory in which the replicate's QC data was written. If this column is missing, it is assumed to be equal to the lanebam.ID. Must be unique!
    \item input.read.pair.count: (OPTIONAL) This column contains the number of input reads, before alignment. This is used later to calculate mapping rate.
    \item multi.mapped.read.pair.count: (OPTIONAL) This column contains the number of reads that were multi-mapped. This must be included for multi-mapping rate to be calculated.
  }
  
  The biological-replicate decoder should have one row per sample.ID, with the following columns:
  \itemize{
    \item sample.ID: The identifier for the specific biological replicate from which the replicate belongs. Note that this is distinct from "unique.ID" because in many RNA-Seq studies each "sample" can have multiple technical replicates, as multiple sequencing runs may be needed to acquire sufficient reads for analysis.
    \item group.ID (OPTIONAL): The identifier for the biological condition for the given replicate.
  }
  
  All decoders are allowed to contain other columns in addition to the ones listed here, so long as their names are distinct. Columns do not need to appear in any particular order, so long as they are named according to the specifications above.
}

\examples{
  #Load the decoder from the example dataset:
  directory <- paste0(system.file("extdata/", 
                      package="QoRTsExampleData", 
                      mustWork=TRUE),"/");
  decoder.file <- system.file("extdata/decoder.txt", 
                              package="QoRTsExampleData", 
                              mustWork=TRUE);
  decoder.data <- read.table(decoder.file,
                             header=TRUE,
                             stringsAsFactors=FALSE);
  print(decoder.data);

  #This command produces the example dataset used in all the other
  #  examples:
  res <- read.qc.results.data(directory, 
                              decoder = decoder.data,
                              calc.DESeq2 = TRUE, 
                              calc.edgeR = TRUE);
  #Note that DESeq2 and edgeR are required in order to 
  #  calculate the size factors using the options above.
  
  #You can also specify incomplete decoders, and use
  #  the following command to fill in the defaults:
  completeAndCheckDecoder(c("SAMP1","SAMP2",
                            "SAMP3","SAMP4",
                            "SAMP5","SAMP6"))

  #You don't actually have to use completeAndCheckDecoder,
  #You can just pass the incomplete decoder directly to QoRTs.
  #For example, to load a small subset of the example data 
  #(without phenotype data):
  res <- read.qc.results.data(paste0(directory,"/ex/"), 
                              decoder = c("SAMP1_RG1","SAMP2_RG1",
                                          "SAMP3_RG1","SAMP4_RG1"));

  #The list of available QC names for use with skip.files:
  names(res@qc.data);
  #Skip some of the files using a command like this:
  res.quick <- read.qc.results.data(directory, 
                              decoder = decoder.data,
                              skip.files=c(
                              "QC.NVC.raw.R1.txt.gz",
                              "QC.NVC.raw.R2.txt.gz",
                              "QC.NVC.lead.clip.R1.txt.gz",
                              "QC.NVC.lead.clip.R2.txt.gz",
                              "QC.NVC.tail.clip.R1.txt.gz",
                              "QC.NVC.tail.clip.R2.txt.gz"
                              ));
}

\seealso{
  \code{\link{build.plotter}}
}\name{makeMultiPlot}
\docType{methods}
\alias{makeMultiPlot}
\alias{makeMultiPlot.basic}
\alias{makeMultiPlot.colorByGroup}
\alias{makeMultiPlot.colorByLane}
\alias{makeMultiPlot.colorBySample}
\alias{makeMultiPlot.highlightSample}
\alias{makeMultiPlot.highlightSample.colorByLane}
\alias{makeMultiPlot.withPlotter}
\title{
   Generating summary multi-plots
}
\description{
   Creates a large multi-frame summary plot, or a report PDF file.
}
\usage{
  makeMultiPlot.basic(res,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = "plot-basic", 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots)
                      
  makeMultiPlot.colorByGroup(res,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = "plot-colorByGroup", 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      ...)
                      
  makeMultiPlot.colorByLane(res,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = "plot-colorByLane", 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots)               

  makeMultiPlot.colorBySample(res,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = "plot-colorByLane", 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots)

  makeMultiPlot.highlightSample(res, curr.sample,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = paste0("plot-sampleHL-",curr.sample), 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots) 

  makeMultiPlot.highlightSample.colorByLane(res, curr.sample,
                      outfile, 
                      outfile.dir = "./",
                      outfile.prefix = paste0("plot-sampleHL-coloredByLane-",curr.sample), 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots)  

  makeMultiPlot.withPlotter(plotter, 
                      res = plotter$res, 
                      outfile,  
                      outfile.dir = "./",
                      outfile.prefix = "plot-custom", 
                      outfile.ext, 
                      plotter.params = list(), 
                      plot.device.name = "curr", 
                      plotting.device.params = list(),
                      verbose = TRUE, 
                      debugMode, 
                      rasterize.large.plots, 
                      rasterize.medium.plots,
                      raster.height = 1050, 
                      raster.width = 1050,
                      separatePlots = FALSE,
                      exclude.autosomes.chrom.rate.plot = TRUE,
                      chromosome.name.style = "UCSC",
                      fig.res = 150, fig.base.height.inches = 7,
                      insertSize.plot.xlim,
                      sequencing.type = c("RNA","Exome","Genome"),
                      nvc.mark.points = TRUE,
                      maxColumns,
                      maxRows,
                      plotList,
                      labelPlots = TRUE,
                      plot = TRUE,
                      \dots)

}
\arguments{

  \item{res}{
    A \code{QoRT_QC_Results} object, created by 
    \code{\link{read.qc.results.data}}.
  }
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{curr.sample}{
    A character string. For the sample highlight summary plots,
    this should be the sample.ID of the sample that is to be
    highlighted.
  }
  \item{outfile}{
    The output file can be specified in two ways: first, by 
    setting the outfile parameter, which should be the full path
    for the output file. Alternatively, the file path will
    be the concatenation of the three parameters:
    outfile.dir, outfile.prefix, and outfile.ext.
    
    If plot.device.name == "curr", then all the outfile 
    parameters will be nonfunctional.
  }
  \item{outfile.dir}{
    A prefix to be appended to the output filename. Usually 
    the path to the destination directory. By default, this
    will be the current working directory.
  }
  \item{outfile.prefix}{
    A second prefix to be appended to the output filename, after
    outfile.dir. This is usually the name of the file, without
    the file extension. By default, a reasonable file name will be
    selected.
  }
  \item{outfile.ext}{
    The file extension. By default, this will be set to match
    the selected graphics device.
  }
  \item{plotter.params}{
    Additional parameters (advanced) used in creation of the Plotter 
    objects. See \link{build.plotter}.
  }
  \item{plot.device.name}{ 
    The method to use to save plots. Several formats and devices are supported:
      \itemize{
        \item \code{"curr"} (default) Do not create output files. Plot directly to the current or default device.
        \item \code{"png"} for standard png compression.
        \item \code{"CairoPNG"} for png compression using package Cairo. Note that this requires the package Cairo.
        \item \code{"tiff"} for the tiff device.
        \item \code{"jpeg"} for the jpeg device (\emph{not recommended!})
        \item \code{"svg"} for the vector svg device. Note that for all vector devices, by default, certain large plots will be rasterized if the \code{png} package is found. If not, the function will still work, but output files may be very large and may have trouble rendering and printing.
        \item \code{"pdf"} for a multi-page pdf report. 
        \item \code{"CairoPDF"} for a multi-page pdf report using package Cairo. Note that this requires the package Cairo. 
      }
  }
  \item{plotting.device.params}{ 
    A named list of parameters to be passed to the graphics device.
    For example:
      \itemize{
        \item \code{width = 2000}
        \item \code{height = 2000}
        \item \code{units = "px"}
      }
    Reasonable values for height, width, and pointsize will be chosen
    by default. Generally all raster plots will be set to 1000 by 1000 
    pixels, and all vector plots will be set to 7 inches by 7 inches.
  }
  \item{verbose}{
    Logical. If TRUE, more information will be printed to the console.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.large.plots}{
    Logical. If TRUE, then if the currently selected plotting device
    is a vector device (or the "curr" device), then certain large plots
    will have their plotting areas rasterized. The axis labels, titles,
    and text will all remain vectorized, only the plotting areas will
    be flattened. \emph{Note that this requires the png package.}
    
    By default, this parameter will be set to TRUE when a vector
    device is selected.
  }
  \item{rasterize.medium.plots}{
    Logical. as rasterize.large.plots, but applies to moderately-large plots. By default, this parameter will be set to TRUE for pdf/CairoPDF output only.
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels. Double-pane plots
    will be twice this width.
  }
  \item{separatePlots}{
    Logical. If TRUE, then separate image files will be saved for 
    each individual plot, rather than saving one large multi-pane 
    plot. Note that this is not compatible with the "curr" device.
    Also note: if this parameter is set to TRUE, then the output
    files will be saved using the outfile.dir, outfile.prefix
    and outfile.ext parameters. The "outfile" parameter cannot be
    set. The files will be saved to the contatenation of 
    outfile.dir, outfile.prefix, the name of the plot type, 
    and then outfile.ext.
  }
  \item{exclude.autosomes.chrom.rate.plot}{
    A logical value indicating whether to exclude autosomes from the plot.
    See \code{\link{makePlot.chrom.type.rates}}
  }
  \item{chromosome.name.style}{
    A string value indicating the style of the chromosome names, 
    and also how to split up the categories. See 
    \code{\link{makePlot.chrom.type.rates}}
  }
  \item{fig.res}{
    Numeric value. The number of pixels per "inch" (for raster devices only).
    For some plotting devices the figure height will be in pixels not inches, and
    will equal this value multiplied by the fig.base.height.inches value.
  }
  \item{fig.base.height.inches}{
    Numeric value. The base height, in inches, of each sub-figure in the plot. This
    will be equal to the height for vector devices, or used to calculate the height
    in pixels using the fig.res value (see above).
  }
  \item{insertSize.plot.xlim}{
    A numeric vector of length 2. The x-axis limits for the insert size plot. By default QoRTs will attempt to find reasonable values for this, 
    but there are always situations where the default behavior is not ideal. Using this parameter you can set
    it explicitly.
  }
  \item{nvc.mark.points}{
    Logical. If TRUE, then points should be marked with shapes on the NVC plots.
  }
  \item{sequencing.type}{
      The type of sequencing data being analyzed. This only changes the default plot set, which can be overriden with the \code{plotList} parameter.
  }
  \item{maxColumns}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxColumns columns. Extra rows will be added to fit all the plots, as needed.
  }
  \item{maxRows}{
    If set, QoRTs will attempt to create a multiplot that has (at most) maxRows rows. Extra columns will be added to fit all the plots as needed. To set the number of rows and columns manually, you can set both maxColumns and maxRows.
  }
  \item{plotList}{
    The list of desired plots.
  }
  \item{labelPlots}{
    Logical. If TRUE, then label each plot with a letter.
  }
  \item{plot}{
    Logical. If FALSE, do not create any plots. (In other words, do nothing but test the ability to create plots).
  } 
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
    Produces large, multi-frame summary plots, with a large number of different
    plots and graphs.
    
    For devices pdf and CairoPDF, this function will produce a multi-page document
    suitable for printing, with (by default) 6 panels on each page.
    
    Note that for all vector devices (svg, pdf, and CairoPDF), by default, 
    QoRTs will attempt to load the \code{png} package. If this package is 
    found, then certain large plots (the NVC plots and the cumulative 
    diversity plots) will be rasterized within the plotting area. The labels
    and text will still be vectorized. If the \code{png} package is not installed, 
    the function will still work, but output files may be very large and may 
    have trouble rendering and printing.
    
    Rasteration of large plots can be turned off via the rasterize.large.plots option.
}


\examples{
  \dontrun{
    data(res,package="QoRTsExampleData");
    makeMultiPlot.colorByGroup(res);
  }
}

\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.gene.assignment.rates}
\docType{methods}
\alias{makePlot.gene.assignment.rates}
\title{
   Gene assignment rate plot
}
\description{
   Plots the rate at which reads are assigned into various 
   categories.
}
\usage{
makePlot.gene.assignment.rates(plotter,  debugMode, singleEndMode, 
                               plot = TRUE, \dots)
                         
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  For each bam-file, this function plots the rate (y-axis) for which the bam-file's read-pairs are assigned to the given categories.
  
The categories are:
  \itemize{
    \item \emph{Unique Gene}: The read-pair overlaps with the exonic segments of one and only one gene. For many downstream analyses tools, such as DESeq, DESeq2, and EdgeR, only read-pairs in this category are used.
    \item \emph{Ambig Gene}: The read-pair overlaps with the exons of more than one gene.
    \item \emph{No Gene}: The read-pair does not overlap with the exons of any annotated gene.
    \item \emph{No Gene, Intronic}: The read-pair does not overlap with the exons of any annotated gene, but appears in a region that is bridged by an annotated splice junction.
    \item \emph{No Gene, 1kb from gene}: The read-pair does not overlap with the exons of any annotated gene, but is within 1 kilobase from the nearest annotated gene.
    \item \emph{No Gene, 10kb from gene}: The read-pair does not overlap with the exons of any annotated gene, but is within 10 kilobases from the nearest annotated gene.
    \item \emph{No Gene, middle of nowhere}: The read-pair does not overlap with the exons of any annotated gene, and is \emph{more} than 10 kilobases from the nearest annotated gene.
  }
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.gene.assignment.rates(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.cigarOp.byCycle}
\docType{methods}
\alias{makePlot.cigarOp.byCycle}
\title{
   Plot Cigar Operator Rate
}
\description{
   Plots the rate at which the given CIGAR operator occurs, by
   read cycle.
}
\usage{
  makePlot.cigarOp.byCycle(plotter, 
           op=c("SoftClip","Del","Ins","HardClip","Pad","Splice","Aln"), 
           r2.buffer = NULL, 
           rate.per.million = TRUE, 
           use.readLength.denominator = TRUE, 
           debugMode, 
           singleEndMode,
           rasterize.plotting.area = FALSE, 
           raster.height = 1000, 
           raster.width = 1000,
           plot = TRUE,
           ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{op}{
    A character string naming which cigar op to plot.
    Must be one of the following:
      \itemize{
        \item \code{"SoftClip"} Soft Clip (Cigar Op "S")
        \item \code{"HardClip"} Hard Clip (Cigar Op "H")
        \item \code{"Del"} Deletion from reference (Cigar Op "D")
        \item \code{"Ins"} Insertion from reference (Cigar Op "I")
        \item \code{"Pad"} Padding (Cigar Op "P")
        \item \code{"Splice"} Splice Junction (Cigar Op "N")
        \item \code{"Aln"} Aligned to reference (Cigar Op "M")
      }
    Note that makePlot.cigarOp.byCycle(plotter,"SoftClip") gives
    identical results as makePlot.clipping, although the default
    value for the \code{rate.per.million} argument is different.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  
  \item{rate.per.million}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
  }

  \item{use.readLength.denominator}{
    Logical. If TRUE, use the read-length counts as the denominator
    when calculating op rates. This is only relevant if operating on
    trimmed reads, where the read length is variable.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: The read cycle (ie. the base-pair position in the read).
  
  y-axis: The rate at which the bases at the given read-cycle is
  clipped off.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.cigarOp.byCycle(plotter, op = "Del");
  makePlot.cigarOp.byCycle(plotter, op = "Ins");
  makePlot.cigarOp.byCycle(plotter, op = "Splice");
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.onTarget}
\docType{methods}
\alias{makePlot.targetDistribution}
\alias{makePlot.onTarget.counts}
\alias{makePlot.onTarget.rates}

\title{
   On-Target Rates
}
\description{
   On-Target Rates.
}
\usage{
makePlot.targetDistribution(plotter,
       plot.rates = TRUE, 
       plot.hist = TRUE, 
       log.y = TRUE, 
       singleEndMode,
       byPair = !singleEndMode, 
       rasterize.plotting.area = FALSE, 
       raster.height = 1000, 
       raster.width = 1000,
       debugMode,
       plot = TRUE,
       ...)

makePlot.onTarget.counts(plotter, 
      y.counts.in.millions = TRUE, 
      debugMode, 
      singleEndMode,
      plot = TRUE, ...)

makePlot.onTarget.rates(plotter,  
      debugMode, 
      singleEndMode,
      plot = TRUE, ...)



}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{y.counts.in.millions}{
    Logical. If TRUE, y axis labels are labelled in millions.
  }
  \item{plot.rates}{
    Logical. If TRUE, plot rates, otherwise plot raw counts.
  }
  \item{plot.hist}{
    Plot a histogram.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis should be log-scaled.
  }
  \item{byPair}{
    Logical. If TRUE, count read-pairs rather than reads. In other words, if two reads overlap, the overlapping region is counted only once.
  }

  \item{rasterize.plotting.area}{
    Logical. If \code{TRUE}, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }

}
\details{
   For each sample run, indicates the amount of time spent running the QoRTs QC data processing tool
}
\value{
   By default, this function returns nothing. If the return.table parameter is TRUE, then it returns a data.frame with the data that was plotted.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.runTimePerformance(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}
\name{makePlot.genebody.coverage}
\docType{methods}
\alias{makePlot.genebody}
\alias{makePlot.genebody.coverage}
\alias{makePlot.genebody.coverage.UMQuartile}
\alias{makePlot.genebody.coverage.lowExpress}
\title{
   Plot Gene-Body coverage distribution
}
\description{
   Plots Gene-Body coverage distribution
}
\usage{
  makePlot.genebody(plotter, 
                  geneset = c("Overall","90-100","75-90","50-75","0-50"),
                  avgMethod = c("TotalCounts", "AvgPercentile"), 
                  plot.medians, 
                  plot.means = TRUE, 
                  debugMode, 
                  singleEndMode, 
                  rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                  plot = TRUE,
                  ... )

  #NOTE: The preferred method to access the below functions
  #  is to use makeplot.genebody and set:
  #  avgMethod = "TotalCounts"
  #DEPRECIATED:
  makePlot.genebody.coverage(plotter, plot.medians,
                       plot.means = TRUE,  
                       debugMode, singleEndMode, 
                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                       plot = TRUE,
                       ...)
  #DEPRECIATED:
  makePlot.genebody.coverage.UMQuartile(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, singleEndMode, 
                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                       plot = TRUE,
                       ...)
  #DEPRECIATED:
  makePlot.genebody.coverage.lowExpress(plotter, plot.medians, 
                       plot.means = TRUE, 
                       debugMode, singleEndMode,
                       rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                       plot = TRUE,
                       ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  
  \item{geneset}{
    The set of genes to plot the gene-body coverage over. Genes are grouped
    into four quantiles by their by their total read counts: the 90-100 quantile,
    the 75-90 quantile, the 50-75 quantile (or "upper-middle quartile"), and the
    0-50 quantile. By default it will plot the overall gene-body coverage across
    all genes.
  }
  \item{avgMethod}{
    The method used to generate average gene body coverage. 
    
    The "TotalCounts" (default) method simply takes the sum of each bin across all genes in the geneset,
    then for each replicate normalizes the coverage into frequencies.
    
    The "AvgPercentile" calculates the coverage distribution across each percentile bin 
    and normalizes the counts so they add up to 1, FOR EACH GENE, and 
    THEN averages those values across all genes in the geneset. 
    
    The latter variant is experimental. the idea was for it to reduce the impact of 
    highly-expressed genes and gain a better estimate of the degree of 5'/3' bias.
  }

  
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average for each bamfile at the bottom of the plot. 
    Overrides \code{plot.means}.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average for each bamfile at the bottom of the plot.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Gene-body quantile. By default this is broken up into
  40 quantiles containing 2.5% each.
  
  y-axis: Rate at which reads falls into the given quantile of the
  genes' bodies.
  
  \code{makePlot.genebody.coverage} plots the gene body coverage across all genes.
  \code{makePlot.genebody.coverage.UMQuartile} plots the gene body coverage across the genes that fall in the upper-middle quartile of total expression for each sample-run (excluding genes with zero reads).
  \code{makePlot.genebody.coverage.lowExpress} plots the gene body coverage across the genes that fall in the lower two quartiles of total expression for each sample-run (excluding genes with zero reads).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.genebody(plotter, geneset = "Overall");
  makePlot.genebody(plotter, geneset = "90-100");
  makePlot.genebody(plotter, geneset = "50-75");
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.missingness.rate}
\docType{methods}
\alias{makePlot.missingness.rate}
\title{
   Plot N-Rate by read cycle
}
\description{
   Plots rate by which the sequencer cannot determine the base, 
   by read cycle.
}
\usage{
  makePlot.missingness.rate(plotter, r2.buffer = NULL, 
                        debugMode, singleEndMode,
                        rasterize.plotting.area = FALSE, 
                        raster.height = 1000, 
                        raster.width = 1000,
                        log.y = FALSE, 
                        ylim = NULL,
                        plot = TRUE,
                        ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{log.y}{
      Logical. If TRUE, the y-axis will be log-scale.
  }
  \item{ylim}{
      Numeric of length 2. Sets the y-axis limits.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Read Cycle
  
  y-axis: Rate at which the sequencer assigns an 'N' base.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.missingness.rate(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.overlap}
\docType{methods}
\alias{makePlot.overlap}
\alias{makePlot.overlap.coverage}
\alias{makePlot.overlapMismatch.perRead}
\alias{makePlot.overlapMismatch.size}
\alias{makePlot.overlapMismatch.byQual.min}
\alias{makePlot.overlapMismatch.byQual.read}
\alias{makePlot.overlapMismatch.byBase}
\alias{makePlot.overlapMismatch.byBase.atScore}
\alias{makePlot.overlapMismatch.byCycle}
\title{
   Overlap Mismatch Rates
}
\description{
   For paired-end data only, plots the rate at which the two reads have point-mismatch on overlapping regions.
}
\usage{
  makePlot.overlap.coverage(plotter,
                            plot.rates = TRUE, 
                            singleEndMode,
                            rasterize.plotting.area = FALSE, 
                            raster.height = 1000, 
                            raster.width = 1000,
                            debugMode,
                            r2.buffer=NULL,
                            plot = TRUE,
                            ...)

  makePlot.overlapMismatch.size(plotter,
                                plot.rates = TRUE,
                                log.y = TRUE,
                                noIndel = TRUE,
                                draw.decade.lines = c(TRUE,TRUE),
                                xlim = NULL,
                                pct.cutoff=0.95,
                                singleEndMode,
                                debugMode,
                                rasterize.plotting.area = FALSE, 
                                raster.height = 1000, 
                                raster.width = 1000,
                                plot = TRUE,
                                ...)
                                
  makePlot.overlapMismatch.byCycle(plotter,
                                   plot.rates = TRUE, 
                                   noIndel = TRUE, 
                                   log.y = TRUE,
                                   singleEndMode = plotter$res@singleEnd,
                                   rasterize.plotting.area = FALSE, 
                                   raster.height = 1000, 
                                   raster.width = 1000,
                                   debugMode,
                                   r2.buffer,
                                   plot = TRUE,
                                   ...)
  
  makePlot.overlapMismatch.byQual.min(plotter,
                                      plot.rates = TRUE,
                                      log.y = TRUE,
                                      noIndel = TRUE,
                                      draw.decade.lines = TRUE,
                                      singleEndMode,
                                      debugMode,
                                      rasterize.plotting.area = FALSE, 
                                      raster.height = 1000, 
                                      raster.width = 1000,
                                      plot = TRUE,
                                      ...)

  makePlot.overlapMismatch.byQual.read(plotter,
                                       plot.rates = TRUE,
                                       noIndel = TRUE,
                                       log.y = TRUE,
                                       draw.decade.lines = TRUE,
                                       singleEndMode,
                                       debugMode,
                                       rasterize.plotting.area = FALSE, 
                                       raster.height = 1000, 
                                       raster.width = 1000,
                                       r2.buffer = NULL,
                                       plot = TRUE,
                                       ...)

  makePlot.overlapMismatch.byBase(plotter, 
                                  noIndel = TRUE, 
                                  plot.rates = TRUE,
                                  log.y = FALSE, 
                                  y.rate.per.kb = FALSE, 
                                  debugMode, 
                                  draw.vert.dotted.lines = TRUE,
                                  singleEndMode, 
                                  plot = TRUE,
                                  ...)
  
  makePlot.overlapMismatch.byBase.atScore(plotter,
                                     atScore = 41,
                                     overlapScoreMethod = c("pair","min"),
                                     plot.rates = TRUE,
                                     noIndel = TRUE,
                                     log.y = FALSE,
                                     singleEndMode,
                                     rasterize.plotting.area = FALSE,
                                     raster.height = 1000,
                                     raster.width = 1000,
                                     debugMode = DEFAULTDEBUGMODE,
                                     plot = TRUE,
                                     ...)
  
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{atScore}{
    Numeric integer. For makePlot.overlapMismatch.byBase.atScore, plot the overlap mismatch rates for each base-pair swap at the given phred score.
    
    If \code{overlapScoreMethod} is "min", then this will plot the rate of overlap mismatches where the minimum of the two reads' quality score is \code{atScore}. If \code{overlapScoreMethod} is "pair", then 
    this will plot the rate of overlap mismatches where both reads have the same quality score, equal to \code{atScore}.
  }
  \item{overlapScoreMethod}{
    Character string. See documentation for atScore, above.
  }
  \item{plot.rates}{
    A logical value indicating whether or not to plot mismatch rates or mismatch counts.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis will be log-scaled.
  }
  \item{noIndel}{
    Logical. If TRUE, only count overlapping reads if both reads have no aligned indels.
  }
  \item{draw.decade.lines}{
    Logical. If TRUE, draw mini tick lines at decade points on the y-axis when plotting in log-scale.
  }
  \item{xlim}{
    Numeric. The x-axis limits. If NULL, the x-limits will be autodetected using the pct.cutoff value.
  }
  \item{pct.cutoff}{
    Numeric. The percentile cutoff value for the x-axis upper limit.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{r2.buffer}{
      Numeric. The desired distance between the read1/read2 sub-plots.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
 
  
  \item{y.rate.per.kb}{
    Logical. If TRUE, the y axis should be labelled in overlap rate per kilobase.
  }

  \item{draw.vert.dotted.lines}{
    Logical. If TRUE, draw dotted lines at decade y axis points.
  }
  
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  These plotting functions create plots that summarize the overlap-mismatch between paired end reads. 
  When paired-end reads overlap, the sequence on the two reads can be compared in order to estimate 
  the sequencer error rate. This error rate can be visualized in a variety of ways to detect sequencer artifacts
  or errors.
  
  makePlot.overlap.coverage: The rates at which each position in read 1 and read 2 are covered by the other read.
     x-axis: Read position.
     y-axis: How often the given position on the given read is covered by the other read.
  
  makePlot.overlap.mismatchCountPerRead: A histogram of the amount of mismatch between the two reads.
     x-axis: Number of mismatches.
     y-axis: How often a read is observed with the given number of mismatches.
  
  makePlot.overlap.mismatch.byMinQual: The mismatch rate as a function of the minimum of the two base quality scores.
     x-axis: Minimum quality score
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byAvgQual: The mismatch rate as a function of the average of the two base quality scores.
     x-axis: Average quality score
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byReadQual: The mismatch rate as a function of each read's quality score.
     x-axis: Quality score for read 1 and read 2.
     y-axis: Mismatch rate.
  
  makePlot.overlap.mismatch.byBase: The rate at which overlap-mismatches occur for each possible base-pair swap.
     x-axis: Read1/read2 base-pair values (note: read2 base is complemented).
     y-axis: Mismatch rate.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.overlap.coverage(plotter)
  makePlot.overlapMismatch.size(plotter)
  makePlot.overlapMismatch.byCycle(plotter)
  makePlot.overlapMismatch.byQual.min(plotter)
  makePlot.overlapMismatch.byQual.read(plotter)
  makePlot.overlapMismatch.byBase(plotter)
  makePlot.overlapMismatch.byBase.atScore(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.insert.size}
\docType{methods}
\alias{makePlot.insert.size}
\title{
   Plot Insert Size Distribution.
}
\description{
   Plots the distribution of the size of the region covered by the
   two paired reads.
}
\usage{
  makePlot.insert.size(plotter, calc.rate = TRUE, pct.cutoff = 0.98,
                   plot.medians = TRUE, plot.means = NULL,
                   debugMode, singleEndMode, xlim = NULL,
                   rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                   plot = TRUE,
                   \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{calc.rate}{
    A logical value indicating whether or not to calculate and
    plot the rate at which each insert size occurs as the y-axis.
    If this is set to false, it will instead plot the total number
    of times each insert size occurs.
  }
  \item{pct.cutoff}{
    A numeric value between 0 and 1, indicating the quantile 
    within which to limit the x-axis. Generally the default value
    of 0.98 is perfectly usable.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average for each bamfile at the bottom of the plot.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average for each bamfile at the bottom of the plot. 
    Overrides \code{plot.means}.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{xlim}{
      Numeric vector of length 2. Sets the x-axis limits. By default
      QoRTs will attempt to autodetect reasonable values for this, but
      there are always cases where you want to set this explicitly.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: The insert size. This is the genomic distance from the 
  start of one read to the end of the other. In other words, the
  size of the full region covered by the paired reads.
  
  y-axis: The rate at which the given insert size occurs, for each
  bamfile.
  
  Note: The insert size is calculated by using the alignment  
  as well as the provided gene/splice-junction annotation.
  \enumerate{
    \item If the paired reads align to overlapping regions of the
    reference genome, then the insert size can be calculated
    exactly. This is NOT dependant on the annotation. If the two 
    reads align inconsistantly (for example, one read showing a 
    splice junction and the other not), then the read is ignored.
    \item If the paired reads do NOT overlap, then the annotation
    information is required. Using the full set of all known splice
    junctions that lie between the inner alignment endpoints for 
    the two reads, the shortest possible path from the two 
    endpoints is found. 
    In some rare cases this may underestimate the insert size, if
    the actual insert does not take the minimal path, but this is
    rare.
    In other (somewhat more common) cases this may overestimate 
    the insert size, when the region between the paired reads 
    bridges an unannotated splice junction.
  }
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.insert.size(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.norm.factors}
\docType{methods}
\alias{makePlot.norm.factors}
\alias{makePlot.norm.factors.vs.TC}
\title{
   Plot Alignment Clipping
}
\description{
   Plots the rate at which the aligner soft-clips off portions of
   aligned reads.
}
\usage{
  makePlot.norm.factors(plotter, by.sample = TRUE, 
                return.table = FALSE,
                debugMode, singleEndMode,
                plot = TRUE,
                ...)
  makePlot.norm.factors.vs.TC(plotter, 
                by.sample = TRUE, 
                return.table = FALSE,
                debugMode, singleEndMode,
                plot = TRUE,
                ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{by.sample}{
    Logical. Whether to combine all lanebams for each sample before 
    calculating normalization factors. By default, normalization factors
    for each lanebam will be calculated seperately.
  }
  \item{return.table}{
    Logical. Whether to return a data table containing the data that is plotted.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  makePlot.norm.factors plots the normalization factors for each sample/lanebam. 
  Note that unless DESeq2 and edgeR are installed, it will only plot total-count
  normalization by default. Also note that unless calc.DESeq2 = TRUE and/or
  calc.edgeR = TRUE in the execution of the read.qc.results.data that produced the
  QC results, only the total counts normalization factors will be calculated.
  
  makePlot.norm.factors.vs.TC plots the ratio of alternative normalization factors against
  the total count normalization.
}
\value{
  Usually returns nothing, unless return.table is TRUE, in which case it returns a data.frame containing the plotted data for each sample. 
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.norm.factors(plotter);
  makePlot.norm.factors.vs.TC(plotter);
  
  #Legend:
  makePlot.legend.box(plotter);
}
\seealso{
  \code{\link{build.plotter}}

}\name{makePlot.legend.box}
\docType{methods}
\alias{makePlot.legend.box}
\alias{makePlot.legend.over}
\title{
   Plot legend
}
\description{
   Plots the universal legend for a given plotter object.
}
\usage{
makePlot.legend.box(plotter, debugMode, singleEndMode, 
                    cex.legend, ncol, 
                    plot = TRUE, 
                    \dots)
makePlot.legend.over(position, 
                 plotter, 
                 debugMode ,
                 singleEndMode, 
                 ncol = NULL,
                 plot = TRUE, 
                 \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{position}{
    For makePlot.legend.over, a character string indicating the location 
    where you want the legend to appear. This is passed to 
    \code{\link{legend}}, and can be any keyword allowed by 
    \code{\link{xy.coords}}. For example: "top", 
    "topleft","bottomright", etc.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{cex.legend}{
      Numeric. The cex expansion factor passed to \code{legend}. By default
      the cex is autofitted to fill the available space.
  }
  \item{ncol}{
      Integer value. The number of columns in the legend. By default this will be
      selected automatically depending on the number of items in the legend.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   makePlot.legend.box creates a new plot (opening the next graphics
   frame), and writes a small description of the given plotter 
   type along with plotting a color-coded legend (if applicable).
   
   makePlot.legend.over adds a legend to the current plotting frame.
}
\examples{
data(res,package="QoRTsExampleData");
plotter <- build.plotter.colorByGroup(res);
#Add a legend to an existing plot:
makePlot.strandedness.test(plotter);
makePlot.legend.over("topright", plotter)

#Or make a legend as a separate plot:
makePlot.legend.box(plotter);
}


\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.cigarLength.distribution}
\docType{methods}
\alias{makePlot.cigarLength.distribution}
\title{
   Plot the length distribution of a given cigar operation
}
\description{
   Plots the length distribution of a given cigar operation
}
\usage{
  makePlot.cigarLength.distribution(plotter, op, 
                                r2.buffer = NULL,
                                perMillion = TRUE,
                                log.x = FALSE, 
                                log.y = FALSE, 
                                debugMode, singleEndMode,
                                rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
                                plot = TRUE,
                                \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{op}{
    A character string naming which cigar op to plot.
    Must be one of the following:
      \itemize{
        \item \code{"SoftClip"} Soft Clip (Cigar Op "S")
        \item \code{"HardClip"} Hard Clip (Cigar Op "H")
        \item \code{"Del"} Deletion from reference (Cigar Op "D")
        \item \code{"Ins"} Insertion from reference (Cigar Op "I")
        \item \code{"Pad"} Padding (Cigar Op "P")
        \item \code{"Splice"} Splice Junction (Cigar Op "N")
        \item \code{"Aln"} Aligned to reference (Cigar Op "M")
      }
    Note that makePlot.cigarOp.byCycle(plotter,"SoftClip") gives
    identical results as makePlot.clipping, although the default
    value for the \code{rate.per.million} argument is different.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1 and read 2.
    By default this will choose a reasonable value.
  }
  \item{perMillion}{
    A logical value indicating whether or not to scale the y 
    axis to rate-per-million-reads, rather than rate-per-read.
    Some people may find the results more readable this way, even
    though the plots themselves will appear the same.
  }
  \item{log.x}{
    A logical value indicating whether or not to log-scale the x 
    axis.
  }
  \item{log.y}{
    A logical value indicating whether or not to log-scale the y 
    axis.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Length of the cigar operation.
  
  y-axis: Frequency of the given length.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.cigarLength.distribution(plotter, op = "Del");
  makePlot.cigarLength.distribution(plotter, op = "Ins");
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.raw.NVC}
\docType{methods}
\alias{makePlot.NVC}
\alias{makePlot.raw.NVC}
\alias{makePlot.minus.clipping.NVC}
\title{
   Nucleotide-rate by read cycle plot
}
\description{
   Plots the nucleotide rate by read cycle
}
\usage{
makePlot.raw.NVC(plotter,  r2.buffer = NULL, 
                         points.highlighted = TRUE,
                         label.majority.bases = FALSE, 
                         label.majority.bases.threshold = 0.5, 
                         label.majority.bases.cex = 0.5, 
                         rasterize.plotting.area = FALSE, 
                         raster.height = 1000, 
                         raster.width = 2000,
                         show.base.legend = TRUE, 
                         debugMode, singleEndMode,
                         useFQ = FALSE,
                         plot = TRUE,
                         \dots)
                         
makePlot.minus.clipping.NVC(plotter,  r2.buffer = NULL, 
                         points.highlighted = TRUE,
                         label.majority.bases = FALSE, 
                         label.majority.bases.threshold = 0.5, 
                         label.majority.bases.cex = 0.5, 
                         rasterize.plotting.area = FALSE, 
                         raster.height = 1000, 
                         raster.width = 2000,
                         show.base.legend = TRUE, 
                         debugMode, singleEndMode,
                         plot = TRUE,
                         \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{r2.buffer}{
    Buffer space to place between the plotting of read 1
    and read 2. By default this will choose a reasonable value.
  }
  \item{points.highlighted}{
    A logical value. Determines whether to ever plot points on top of the lines. This
    can be useful for identifying outliers. If TRUE, then all highlighted lane-bams
    will be overlaid with points using their designated pch symbol. If the plotter does
    not highlight any lanebams, then points will be overlaid on ALL lines.
  }
  \item{label.majority.bases}{
    A logical value indicating whether to label the genotypes of
    read cycles in which the most common base has a frequency 
    exceeding label.majority.bases.threshold (see below).
  }
  \item{label.majority.bases.threshold}{
    A numeric value indicating the cutoff above which the most 
    common base will be labelled, if label.majority.bases is set
    TRUE (see above).
  }
  \item{label.majority.bases.cex}{
    The cex value fed to text() that is used to determine how
    big the labels are to be, if label.majority.bases is TRUE.
    (see \code{\link{par}} for information on cex).
  }
  \item{rasterize.plotting.area}{
    Logical. If \code{TRUE}, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  \item{show.base.legend}{
    A logical value indicating whether to print the base-color
    legend along the right edge of the plot.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{useFQ}{
    Logical. If TRUE, plot the G/C rate from the unaligned FASTQ data. This requires that the FASTQ data was supplied to QoRTs in the QoRTs QC step.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   This plot displays the nucleotide rates by read cycle. The color scheme for NVC plots is different from all other plots. Rather than being used for emphasis or to allow cross-comparisons by sample, biological-condition, or lane, the colors are used to indicate the four nucleotides: A (green), T (red), G (orange), or C (blue). Depending on the type of QoRT_Plotter being used, sample-runs will be marked and differentiated by marking the lines with shapes (R points). In many cases the points will be unreadable due to overplotting, but clear outliers that stray from the general trends can be readily identified.
   
   When used with a "sample.highlight" type QoRT_Plotter (see \code{\link{build.plotter}}), "highlighted" samples will be drawn with a deeper shade of the given color.
   
   \code{\link{makePlot.raw.NVC}} plots the nucleotides of all cycles of all aligned reads. \code{\link{makePlot.minus.clipping.NVC}} plots the nucleotides for all cycles that are NOT soft-clipped by the aligner. Note that for reads aligned with an aligner that does not perform soft-clipping (such as TopHat2), the two plots will be identical.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.minus.clipping.NVC(plotter);
  makePlot.raw.NVC(plotter);
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.reference}
\docType{methods}

\alias{makePlot.reference}
\alias{makePlot.referenceMismatch}
\alias{makePlot.referenceMismatch.byScore}
\alias{makePlot.referenceMismatch.perRead}
\alias{makePlot.referenceMismatch.byBase}
\alias{makePlot.referenceMismatch.byBase.atScore}
\alias{makePlot.referenceMismatch.byCycle}

\title{
   Plot Reference Mismatch Rates
}
\description{
   Plots various rates of single-base mismatches against the reference genome.
}
\usage{
makePlot.referenceMismatch.byScore(plotter,
             plot.rates = TRUE,
             draw.decade.lines = c(FALSE,TRUE),
             log.y = TRUE,
             singleEndMode = plotter$res@singleEnd,
             debugMode = DEFAULTDEBUGMODE,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             r2.buffer = NULL,
             plot = TRUE,
             ...)

makePlot.referenceMismatch.byBase(plotter, 
             y.rate.per.kb = TRUE, draw.vert.dotted.lines = TRUE,
             debugMode = DEFAULTDEBUGMODE, 
             singleEndMode = plotter$res@singleEnd, 
             plot = TRUE,
             ...)

makePlot.referenceMismatch.byBase.atScore(plotter, 
             atScore = 41,
             forRead = c("BOTH","R1","R2"),
             plot.rates = TRUE,
             log.y = FALSE,
             singleEndMode = plotter$res@singleEnd,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             debugMode = DEFAULTDEBUGMODE,
             plot = TRUE,...)

makePlot.referenceMismatch.byCycle(plotter,
             plot.rates = TRUE,
             log.y = TRUE,
             ylim = NULL,
             singleEndMode = plotter$res@singleEnd,
             debugMode = DEFAULTDEBUGMODE,
             rasterize.plotting.area = FALSE, 
             raster.height = 1000, 
             raster.width = 1000,
             r2.buffer = NULL,
             plot = TRUE,
             ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.rates}{
    A logical value indicating whether or not to plot mismatch rates or mismatch counts.
  }
  \item{log.y}{
    Logical. If TRUE, the y-axis will be log-scaled.
  }
  \item{draw.decade.lines}{
    Logical. If TRUE, draw mini tick lines at decade points on the y-axis when plotting in log-scale.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{r2.buffer}{
      Numeric. The desired distance between the read1/read2 sub-plots.
  }
  \item{y.rate.per.kb}{
    Logical. If TRUE, then y-axis should be labelled as the rate per kilobase.
  }
  \item{draw.vert.dotted.lines}{
    Logical. If TRUE, then draw dotted guide-lines at reasonable intervals.
  }
  \item{atScore}{
    Numeric integer. For \code{makePlot.referenceMismatch.byBase.atScore}, this sets the PHRED score.
  }
  \item{forRead}{
    Character string. For \code{makePlot.referenceMismatch.byBase.atScore}, this sets whether the plot should be for read 1, read 2, or both side by side.
  }
  \item{ylim}{
    Numeric vector of length 2. This sets the y-axis limits.
  }
  
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  These plotting functions create plots that summarize the rate of single-base-pair reference-mismatches of reads against the reference genome. A "reference mismatch"
  is defined as an aligned read base that does not match the genomic base.

  These reference mismatches may be caused by a number of different factors and may not necessarily indicate a serious quality issue. To begin with, many
  reference mismatches may be caused by real single nucleotide polymorphisms present in the subject. QoRTs produces several different plots that allow these 
  mismatch rates to be measured in a variety of ways.
  
  \emph{\code{makePlot.referenceMismatch.byCycle}}: Plots the overall reference mismatch rate as a function of read cycle.

  \emph{\code{makePlot.referenceMismatch.byBase}}: Plots the overall reference mismatch rate for each of the 12 possible X -> Y base-pair substitutions.
  
  \emph{\code{makePlot.referenceMismatch.byScore}}: Plots the reference mismatch rate as a function of the PHRED score.
  
  \emph{\code{makePlot.referenceMismatch.byBase.atScore}}: Plots the reference mismatch rate for each of the 12 possible X -> Y base-pair substitutions, for the
  subset of bases with PHRED score equal to \code{atScore}.

  \emph{NOTE:} Creation of the reference mismatch plots requires that the QoRTs QC run be executed with additional optional parameters. The BAM file must be sorted
  by coordinate, and a genome fasta file must be specified via the \code{--genomeFA} parameter.

}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.mapping.rates}
\docType{methods}
\alias{makePlot.mapping.rates}
\title{
   Plot mapping rates
}
\description{
   Plots the rates at which reads map to the genome.
}
\usage{
makePlot.mapping.rates(plotter, plot.mm = NULL, 
                   y.counts.in.millions = TRUE, 
                   debugMode, singleEndMode,
                   plot = TRUE,
                   \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.mm}{
    Plot multi-mapping rates. By default, it autodetects whether multimapping data was included in the original decoder.
  }
  \item{y.counts.in.millions}{
    Label/scale the y-axis in millions of reads.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
For each sample-run, this function plots the number of mapped reads and the rate at which reads map (if the total reads is provided).
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.mapping.rates(plotter);
}
\seealso{
  \code{\link{build.plotter}}
  \code{\link{read.qc.results.data}}
}\name{makePlot.NVC.clip.matchByClipPosition}
\docType{methods}
\alias{makePlot.NVC.clip.matchByClipPosition}
\alias{makePlot.NVC.lead.clip.matchByClipPosition}
\alias{makePlot.NVC.tail.clip.matchByClipPosition}
\title{
   Plot clipped nucleotide rates, aligned by the distance 
   from the point of clipping.
}
\description{
   WARNING: THESE FUNCTIONS ARE BETA, AND ARE NOT FULLY TESTED OR 
   READY FOR GENERAL USE.
}
\usage{
makePlot.NVC.lead.clip.matchByClipPosition(plotter, 
            clip.amt,  r2.buffer, 
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            rasterize.plotting.area = FALSE, 
            raster.height = 1000, 
            raster.width = 1000,
            show.base.legend = TRUE, 
            load.results = TRUE, 
            debugMode, singleEndMode, 
            plot = TRUE,
            ...)

makePlot.NVC.tail.clip.matchByClipPosition(plotter, 
            clip.amt, r2.buffer,
            label.majority.bases = TRUE, 
            label.majority.bases.threshold = 0.5, 
            label.majority.bases.cex = 1, 
            rasterize.plotting.area = FALSE, 
            raster.height = 1000, 
            raster.width = 1000,
            show.base.legend = TRUE, 
            debugMode, singleEndMode, 
            plot = TRUE,
            ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{clip.amt}{
    UNDOCUMENTED
  }
  \item{r2.buffer}{
    UNDOCUMENTED
  }
  \item{label.majority.bases}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.threshold}{
    UNDOCUMENTED
  }
  \item{label.majority.bases.cex}{
    UNDOCUMENTED
  }
  \item{show.base.legend}{
    UNDOCUMENTED
  }
  \item{rasterize.plotting.area}{
    Logical. If TRUE, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  \item{load.results}{
    UNDOCUMENTED
  }
  \item{debugMode}{
    UNDOCUMENTED
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   This function is currently BETA, and is not intended for general use. Documentation and testing is still pending.
}

\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.NVC}}
}\name{makePlot.gene.cdf}
\docType{methods}
\alias{makePlot.gene.cdf}
\title{
   Plot the cumulative gene diversity curve
}
\description{
   Plots the cumulative gene diversity curve
}
\usage{
  makePlot.gene.cdf(plotter, 
                sampleWise = FALSE,
                plot.intercepts = TRUE,
                label.intercepts = FALSE,
                debugMode, 
                rasterize.plotting.area = FALSE, 
                raster.height = 1000, 
                raster.width = 1000,
                singleEndMode,
                plot = TRUE,
                ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{sampleWise}{
    A logical value indicating whether to compile each sample
    together across all runs.
  }
  \item{plot.intercepts}{
    A logical value indicating whether or not to plot the
    intercepts with the round numbers on the x-axis, in 
    dotted lines.
  }
  \item{label.intercepts}{
    A logical value indicating whether or not to label the
    intercepts. No effect if plot.intercepts is FALSE.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{rasterize.plotting.area}{
    Logical. If TRUE, the plotting area will be written to a temp
    png file then drawn to the current device as a raster image.
    This option is generally preferred for vector devices, because
    NVC plots can be very large when drawn in vector format. 
    \emph{Note: this requires the \code{png} package!}
  }
  \item{raster.height}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    height of the plotting area raster image, in pixels.
  }
  \item{raster.width}{
    Numeric. If rasterize.plotting.area is TRUE, then this is the
    width of the plotting area raster image, in pixels.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  For each bam-file, this function plots the cumulative gene diversity. For each bam-file, the genes are sorted by read-count. Then, a cumulative function is calculated for the percent of the total proportion of reads as a function of the number of genes. Intercepts are plotted as well, showing the cumulative percent for 1 gene, 10 genes, 100 genes, 1000 genes, and 10000 genes. 

  So, for example, across all the bam-files, around 50 to 55 percent of the read-pairs were found to map to the top 1000 genes. Around 20 percent of the reads were found in the top 100 genes. And so on.

  This can be used as an indicator of whether a large proportion of the reads stem from of a small number of genes.

  Note that this is restricted to only the reads that map to a single unique gene. Reads that map to more than one gene, or that map to intronic or intergenic areas are ignored.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.gene.cdf(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.gc}
\docType{methods}
\alias{makePlot.gc}
\title{
   Plot GC content
}
\description{
   Plots GC content.
}
\usage{
  makePlot.gc(plotter, plot.medians = NULL, plot.means = TRUE, 
              plotRate = FALSE, byPair = FALSE, useFQ = FALSE,
              debugMode, singleEndMode, 
              rasterize.plotting.area = FALSE, raster.height = 1000, raster.width = 1000,
              plot = TRUE,
              ...)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{plot.medians}{
    A logical value indicating whether or not to plot the median 
    average GC content for each bam file at the bottom 
    of the plot. Overrides \code{plot.means}.
  }
  \item{plot.means}{
    A logical value indicating whether or not to plot the mean 
    average GC content for each bam file at the bottom 
    of the plot.
  }
  \item{plotRate}{
    A logical value indicating whether or not the X-axis should be 
    the raw number of nucleotides that are G/C, vs the rate G/C.
  }
  \item{byPair}{
    Logical. If FALSE, GC content rates will be calculated for
    all reads individually. If TRUE, GC content rates will be
    calculated for all read-pairs. Note that when the insert
    size is small, the paired plot can sometimes appear jagged,
    as completely-overlapped read-pairs will always have an even
    number of G/C nucleotides.
  }
  \item{useFQ}{
    Logical. If TRUE, plot the G/C rate from the unaligned FASTQ data. This requires that the FASTQ data was supplied to QoRTs in the QoRTs QC step.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
      Logical. Determines whether the dataset consists of single-ended reads. 
      By default this is determined from the data. Thus, this parameter should 
      never need to be set by the user.
  }
  \item{rasterize.plotting.area}{
      Logical. If \code{TRUE}, then "flatten" the plotting lines into a raster format. 
      This requires support for png file creation and the installation of the "png" 
      package. Only the plotting lines will be rasterized, the 
      axes and annotations will be vector format. Default is \code{FALSE}.
  }
  \item{raster.height}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the height of the rasterized plot, in pixels.
  }
  \item{raster.width}{
      Numeric integer. If \code{rasterize.plotting.area} is \code{TRUE}, then this
      will set the width of the rasterized plot, in pixels.
  }
  \item{plot}{
        Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
  x-axis: Percent of the read-pairs that is made up of G's or C's.
  
  y-axis: Rate at which the given percent occurs.
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.gc(plotter)
}
\seealso{
  \code{\link{build.plotter}}
}\name{makePlot.splice.junction.loci.counts}
\docType{methods}
\alias{makePlot.splice.junction.loci.counts}
\title{
   Splice Junction Loci Count Plot
}
\description{
   Plots the rate at which splice junction loci
   fall into various categories.
}
\usage{
makePlot.splice.junction.loci.counts(plotter, 
                                 calc.rate = FALSE, 
                                 high.low.cutoff = 4,
                                 debugMode, singleEndMode, 
                                 plot = TRUE, \dots)
}
\arguments{
  \item{plotter}{
    A \code{QoRT_Plotter} reference object. See \code{\link{build.plotter}}.
  }
  \item{calc.rate}{
    Logical. If TRUE, the proportion of loci in each category will be calculated and plotted, rather than the raw number.
  }
  \item{high.low.cutoff}{
    (For advanced users only!) The cutoff between "high" and "low" expression junction loci. Note that this must match the cutoff used by the jarfile QC execution.
  }
  \item{debugMode}{
    Logical. If TRUE, debugging data will be printed to the console.
  }
  \item{singleEndMode}{
    Logical. Determines whether the dataset consists of single-ended reads. 
    By default this is determined from the data. Thus, this parameter should 
    never need to be set by the user.
  }
  \item{plot}{
      Logical. If FALSE, suppress plotting and return \code{TRUE} if and only if the data is present in the dataset to allow plotting. Useful to automatically filter out missing data plots.
  }
  \item{\dots}{ 
    Arguments to be passed to methods, such as
    \link{graphical parameters} (see \code{\link{par}}).
  }
}
\details{
   This function plots the number (y-axis) of splice junction \emph{loci} of each type that appear in the bam-file's reads. Splice junctions are split into 4 groups, first by whether the splice junction appears in the transcript annotation gtf ("known" vs "novel"), and then by whether the splice junction has 4 or more reads covering it, or 1-3 reads ("Hi" vs "Lo").
}
\examples{
  data(res,package="QoRTsExampleData");
  plotter <- build.plotter.colorByGroup(res);
  makePlot.splice.junction.loci.counts(plotter);
  #Add a legend:
  makePlot.legend.over("topright", plotter)
}
\seealso{
  \code{\link{build.plotter}}, \code{\link{makePlot.splice.junction.event}}
}