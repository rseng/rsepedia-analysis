
# iwa-miRNA:  A Web-based Platform for Interactive Annotation of Plant MiRNAs

[![webserver](https://img.shields.io/badge/Web_server-ready-red.svg)](http://iwa-miRNA.omicstudio.cloud/) [![docker](https://img.shields.io/badge/Docker_image-ready-red.svg)](https://hub.docker.com/r/malab/iwa-mirna/) [![source](https://img.shields.io/badge/Source_code-support-blue.svg)](https://github.com/cma2015/iwa-miRNA/tree/master/Source_code) [![testdata](https://img.shields.io/badge/Test_data-support-blue.svg)](https://github.com/cma2015/iwa-miRNA/tree/master/Test_data)

- iwa-miRNA allows users to generate a comprehensive collection of miRNA candidates, and to interrogate miRNA annotation in a straightforward way, without the need for computational skills.

- iwa-miRNA Docker image is available at https://hub.docker.com/r/malab/iwa-mirna. Source codes and user manual are available at https://github.com/cma2015/iwa-miRNA. The web server of iwa-miRNA is accessible at http://iwa-miRNA.omicstudio.cloud/.

- The iwa-miRNA is composed with three functional modules: MiRNA Compilation, MiRNA Selection, and Manual Curation. It also provides users with some useful tools for downstream exploratory analysis. More details regarding these functional modules can be found [here](https://github.com/cma2015/iwa-miRNA/blob/master/Tutorials/Modules.md).

<img src="Tutorials/_images/Graphical_summary.png" alt="Graphical summary of iwa-miRNA" style="zoom:18%">

## How to use iwa-miRNA

- Tutorial for iwa-miRNA: http://iwa-miRNA.omicstudio.cloud/static/assets/html/index.html
- Test datasets referred in the tutorials for iwa-miRNA: https://github.com/cma2015/iwa-miRNA/tree/master/Test_data

## Abbreviations in iwa-miRNA

- SVM: support vector machine;
- HT criteria: high-throughput criteria;
- PsRNA: [plant small RNA genes](http://plantsmallrnagenes.science.psu.edu/);
- MFE: minimal free energy;
- AMFE: adjusted minimal free energy;
- sRNA-Seq: small RNA sequencing.

## Changelog

- 2020-11-30: Update the description of output results in the tutorial.
- 2020-06-24: A demo server of iwa-miRNA was released for users running small RNA sequencing datasets.
- 2020-06-20: Source codes and Docker image of iwa-miRNA were released for the first time.

## How to access help

- For any feedback and tool suggestions, please feel free to leave a message at Github [issues](https://github.com/cma2015/iwa-miRNA/issues). We will try our best to deal with all issues as soon as possible.
- In addition, if any suggestions are available, feel free to contact: **Ting Zhang** [zting135@gmail.com](mailto:zting135@gmail.com) or **Chuang Ma** [chuangma2006@gmail.com](mailto:chuangma2006@gmail.com)

## Citation

T. Zhang, J. Zhai, X. Zhang, L. Ling, M. Li, S. Xie, M. Song, C. Ma, Interactive Web-based Annotation of Plant MicroRNAs with iwa-miRNA, *Genomics, Proteomics & Bioinformatics* (2021), doi: https://doi.org/10.1016/j.gpb.2021.02.010

### Module I 

This module generates a comprehensive collection of miRNA candidates by aggregating already annotated miRNAs from four plant miRNA databases (i.e., miRBase, PmiREN, sRNAanno, and PsRNA) and predicted miRNAs from user-submitted sRNA-Seq data.

|       Tool       |                     **Input**                      |                            Output                            |                         Applications                         |
| :--------------: | :------------------------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| *miRNARetrival*  |        Name of species and miRNA databases         | Already annotated miRNAs<br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/Overview.html">Overview</a> | Aggregate annotated miRNAs provided by different miRNA databases |
|  *miRNAPredict*  |    SRA accession numbers or uploaded fastq file    | Predicted miRNAs<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNAPredict_output.txt">Prediction</a> |              Predict miRNAs from sRNA-Seq data               |
| *genomeRetrival* | Name of species or genome sequences and annotation |      Path of formatted genome sequences and annotation       |             Get genome sequences and annotation              |
| *miRNATranslate* |     Output from miRNARetrival and miRNAPredict     | miRNA and miRNA precursors with a uniform format<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNATranslate_output.txt">Aggregation</a><br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/Aggregation.html">Report</a> | Translate annotated and predicted miRNAs into the genomic coordinate system |

### Module II

This module selects a subset of miRNA candidates according to the high-throughput criteria and/or using an ML-based approach.

|    **Tools**     |         **Input**          |                          **Output**                          |           Applications            |
| :--------------: | :------------------------: | :----------------------------------------------------------: | :-------------------------------: |
| *miRNASelection* | Output from miRNATranslate | Selected miRNAs<br><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/miRNASelection_output.txt">Output</a> | Select promising miRNA candidates |

### Module III

This module provides the information for all miRNA candidates  for rapid curation of the quality of selected miRNAs.

|    **Tools**     |          **Input**          |                          **Output**                          |          Applications           |
| :--------------: | :-------------------------: | :----------------------------------------------------------: | :-----------------------------: |
| *manualCuration* | Output from MiRNA Selection | Summary and report pages<br/><a href="http://iwa-miRNA.omicstudio.cloud/static/assets/Test_results/manualCuration_output.html">Output</a> | Determine the quality of miRNAs |



> The table above briefly describes each function in this module. The sample data are displayed through a hyperlink. Users can preview the results by clicking it. Due to the limitation of file size, some results may not be displayed properly. For more detailed descriptions about inputs and outputs, please refer to [tutorial for iwa-miRNA](http://iwa-miRNA.omicstudio.cloud/static/assets/html/index.html) or the interpretation of each functional analysis interface in [web server](http://iwa-miRNA.omicstudio.cloud/).

# enaBrowserTools

enaBrowserTools is a set of scripts that interface with the ENA web services to download data from ENA easily, without any knowledge of scripting required.

*Important: Python 2 based scripts in the python/ folder are deprecated. Use the scripts in python3/

Please note that v1.5.5 is the last update that will support python 2.7. Due to the 2019 scheduled retirement of python 2.7, the next release (1.6.0) will be restructured to hold only python 3 scripts.*

# License

Copyright 2017 EMBL - European Bioinformatics Institute Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

# Installation and Setup

## Python installation

Both Python 2 and Python 3 scripts are available.  The Python 2 scripts can be found in the "python" folder, and the Python 3 scripts in the "python3" folder.   

To run these scripts you will need to have Python installed.  You can download either Python 2 or Python 3 from [here](https://www.python.org/downloads/). If you already have Python installed, you can find out which version when you start the python interpreter.  If using Python 2, we suggest you upgrade to the latest version if you don't already have it: 2.7.

Note that EBI now uses HTTPS servers. This can create a problem when using Python 3 on a Mac due to an oft-missed
installation step. Please run the "Install Certificates.command" command to ensure your Python 3 installation on
the Mac can correctly authenticate against the servers. To do this, run the following from a terminal window, updating
the Python version with the correct version of Python 3 that you have installed:
open "/Applications/Python 3.6/Install Certificates.command"

We have had a report from a user than when Python 3 was installed using homebrew, the following code needed to be run instead:
```
# install_certifi.py
#
# sample script to install or update a set of default Root Certificates
# for the ssl module.  Uses the certificates provided by the certifi package:
#       https://pypi.python.org/pypi/certifi

import os
import os.path
import ssl
import stat
import subprocess
import sys

STAT_0o775 = ( stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
             | stat.S_IRGRP | stat.S_IWGRP | stat.S_IXGRP
             | stat.S_IROTH |                stat.S_IXOTH )

openssl_dir, openssl_cafile = os.path.split(
    ssl.get_default_verify_paths().openssl_cafile)

print(" -- pip install --upgrade certifi")
subprocess.check_call([sys.executable,
    "-E", "-s", "-m", "pip", "install", "--upgrade", "certifi"])

import certifi

# change working directory to the default SSL directory
os.chdir(openssl_dir)
relpath_to_certifi_cafile = os.path.relpath(certifi.where())
print(" -- removing any existing file or link")
try:
    os.remove(openssl_cafile)
except FileNotFoundError:
    pass
print(" -- creating symlink to certifi certificate bundle")
os.symlink(relpath_to_certifi_cafile, openssl_cafile)
print(" -- setting permissions")
os.chmod(openssl_cafile, STAT_0o775)
print(" -- update complete")
```

## Installing and running the scripts

Download the [latest release](https://github.com/enasequence/enaBrowserTools/releases/latest) and extract it to the preferred location on your computer. You will now have the enaBrowserTools folder, containing both the python 2 and 3 option scripts.  If you are using a Unix/Linux or Mac computer, we suggest you add the following aliases to your .bashrc or .bash_profile file. Where INSTALLATION_DIR is the location where you have saved the enaBrowserTools to and PYTHON_CHOICE will depend on whether you are using the Python 2 or Python 3 scripts.

```
alias enaDataGet=INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaDataGet
alias enaGroupGet=INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaGroupGet
```

This will allow you to run the tools from any location on your computer.

You can run install and run these scripts on Windows as you would in Unix/Linux using [Cygwin](https://cygwin.com). If you have python installed on your Windows machine, you can run the python scripts directly without Cygwin. However the call is a bit more complicated.

For example, instead of calling ```enaDataGet```

you would need to call ```python INSTALLATION_DIR/enaBrowserTools/PYTHON_CHOICE/enaDataGet.py```

We will look more into the Windows equivalent of aliases to run batch files from the command line and hopefully be able to provide a better solution to Windows users shortly.

## Setting up for Aspera

*Important: there has been a change to using aspera in version 1.4.*

If you wish to use Aspera to download read or analysis files, you will need to use the aspera_settings.ini file.  Please save it to a static location on your local computer, and edit the file to include the location of your aspera binary (ASPERA_BIN) and the private key file (ASPERA_PRIVATE_KEY):

```
[aspera]
ASPERA_BIN = /path/to/ascp
ASPERA_PRIVATE_KEY = /path/to/aspera_dsa.openssh
ASPERA_OPTIONS =
ASPERA_SPEED = 100M
```

There are two command line flags/options available if you wish to use aspera for download. These are:
1. -a (or --aspera): use this flag if you'd like to download with aspera.
2. -as (or --aspera-settings): use this option if you'd like to specify the location of your aspera settings file. If you use this option, you don't need to use the --aspera flag.

If you use the --aspera-settings option, you don't need to also use the --aspera flag, e.g:

```
enaDataGet -f fastq -as /path/to/aspera_settings.ini ACCESSION
```

If you don't wish to specify the location of the aspera settings file each time you use the scripts, you have the option to either set an ENA_ASPERA_INIFILE environment variable to save the location:

```
export ENA_ASPERA_INIFILE="/path/to/aspera_settings.ini"
```

or you can use the default location for the file, this is the enaBrowserTools directory.  Note that if you use this option, you will have to be careful about how you update your scripts so that you don't overwrite your aspera settings file.

Using just the --aspera flag will result in the scripts looking first for the ENA_ASPERA_INIFILE environment variable, and second for the default file location.

```
enaDataGet -f fastq -a ACCESSION
```

Regardless of which option you have selected, if the aspera settings file cannot be found or the licence key file declared within your settings file does not exist, the scripts will default to using FTP for the download.

# Command line

There are two main tools for downloading data from ENA:  enaDataGet and enaGroupGet.  

## enaDataGet

This tool will download all data for a given sequence, assembly, read or analysis accession or WGS set.  Usage of this tool is described below.  Note that unless a destination directory is provided, the data will be downloaded to the directory from which you run the command. When using an assembly, run, experiment or analysis accession, a subdirectory will be created using that accession as its name.

Accepted WGS set accession formats are:
- AAAK03
- AAAK03000000
- AAAK
- AAAK00000000

```
usage: enaDataGet [-h] [-f {embl,fasta,submitted,fastq,sra}] [-d DEST] [-w]
                  [-m] [-i] [-a] [-as ASPERA_SETTINGS] [-v]
                  accession

Download data for a given accession

positional arguments:
  accession             Sequence, coding, assembly, run, experiment or
                        analysis accession or WGS prefix (LLLLVV) to download

optional arguments:
  -h, --help            show this help message and exit
  -f {embl,fasta,submitted,fastq,sra}, --format {embl,fasta,submitted,fastq,sra}
                        File format required. Format requested must be
                        permitted for data type selected. sequence, assembly
                        and wgs accessions: embl(default) and fasta formats.
                        read group: submitted, fastq and sra formats. analysis
                        group: submitted only.
  -d DEST, --dest DEST  Destination directory (default is current running
                        directory)
  -w, --wgs             Download WGS set for each assembly if available
                        (default is false)
  -e, --extract-wgs     Extract WGS scaffolds for each assembly if available
                        (default is false)
  -exp, --expanded      Expand CON scaffolds when downloading embl format
                        (default is false)
  -m, --meta            Download read or analysis XML in addition to data
                        files (default is false)
  -i, --index           Download CRAM index files with submitted CRAM files,
                        if any (default is false). This flag is ignored for
                        fastq and sra format options.
  -a, --aspera          Use the aspera command line client to download,
                        instead of FTP.
  -as ASPERA_SETTINGS, --aspera-settings ASPERA_SETTINGS
                        Use the provided settings file, will otherwise check
                        for environment variable or default settings file
                        location.
  -v, --version         show program's version number and exit
```

## enaGroupGet

This tool will allow you to download all data of a particular group (sequence, WGS, assembly, read or analysis) for a given sample or study accession. You can also download all sequence, WGS or assembly data for a given NCBI tax ID. When fetching data for a tax ID, the default is to only search for the specific tax ID, however you can use the subtree option to download the data associated with either the requested taxon or any of its subordinate taxa in the NCBI taxonomy tree.

Downloading read and analysis data by tax ID is currently disabled. We will be adding a data volume sanity check in place before we enable this functionality.

Usage of this tool is described below.  A new directory will be created using the provided accession as the name, and all data will be downloaded here. There will also be a separate subdirectory created for each assembly, run and analysis being fetched.  Note that unless a destination directory is provided, this group directory will be created in the directory from which you run the command.  

```
usage: enaGroupGet [-h] [-g {sequence,wgs,assembly,read,analysis}]
                   [-f {embl,fasta,submitted,fastq,sra}] [-d DEST] [-w] [-m]
                   [-i] [-a] [-as ASPERA_SETTINGS] [-t] [-v]
                   accession

Download data for a given study or sample, or (for sequence and assembly) taxon

positional arguments:
  accession             Study or sample accession or NCBI tax ID to fetch data
                        for

optional arguments:
  -h, --help            show this help message and exit
  -g {sequence,wgs,assembly,read,analysis}, --group {sequence,wgs,assembly,read,analysis}
                        Data group to be downloaded for this
                        study/sample/taxon (default is read)
  -f {embl,fasta,submitted,fastq,sra}, --format {embl,fasta,submitted,fastq,sra}
                        File format required. Format requested must be
                        permitted for data group selected. sequence, assembly
                        and wgs groups: embl and fasta formats. read group:
                        submitted, fastq and sra formats. analysis group:
                        submitted only.
  -d DEST, --dest DEST  Destination directory (default is current running
                        directory)
  -w, --wgs             Download WGS set for each assembly if available
                        (default is false)
  -e, --extract-wgs     Extract WGS scaffolds for each assembly if available
                        (default is false)
  -exp, --expanded      Expand CON scaffolds when downloading embl format
                        (default is false)
  -m, --meta            Download read or analysis XML in addition to data
                        files (default is false)
  -i, --index           Download CRAM index files with submitted CRAM files,
                        if any (default is false). This flag is ignored for
                        fastq and sra format options.
  -a, --aspera          Use the aspera command line client to download,
                        instead of FTP.
  -as ASPERA_SETTINGS, --aspera-settings ASPERA_SETTINGS
                        Use the provided settings file, will otherwise check
                        for environment variable or default settings file
                        location.
  -t, --subtree         Include subordinate taxa (taxon subtree) when querying
                        with NCBI tax ID (default is false)
  -v, --version         show program's version number and exit
```

# Tips

From version 1.4, when downloading read data if you use the default format (that is, don't use the format option), the scripts will look for available files in the following priority: submitted, sra, fastq.

A word of advice for read formats:
- submitted: only read data submitted to ENA have files available as submitted by the user.
- sra:  this is the NCBI SRA format, and is the format in which all NCBI/DDBJ data is mirrored to ENA.
- fastq:  not all submitted format files can be converted to FASTQ

# Problems

For any problems, please contact the ENA helpdesk (https://www.ebi.ac.uk/ena/browser/support) with 'enaBrowserTools' in your subject line.

We have had a couple of reports that the R2 FASTQ files are not always successfully downloading for paired runs. We have been unable to replicate this problem and have therefore exposed the error message to you should there be any download failure of read/analysis files via FTP or Aspera. If you get one of these errors, please copy the error into the support form.
**sRNA-Seq_test_data**

This file contains raw sequencing reads in FASTQ format, which can be compressed into **sRNA-Seq_test_data.zip** and used as input of **miRNAPredict** module. In detail, two Arabidopsis sRNA-seq data (SRR11347201 and SRR11829907) were selected, and the first 10 million rows were packaged and named as test1.fastq.gz and test2.fastq.gz.

**sample_information.txt**

This file contains the tissue name of corresponding samples involved in **sRNA-Seq_test_data**. It has two columns with tab-delimited format. The first column is sample names, and the second column is tissue names.

**gene_description.txt**

This file contains genes descriptions and is used as input of **manualCuration** function. It has two columns with tab-delimited format. The first column is the gene names and the second column is the gene description. Gene descriptions were extracted from [Arabidopsis GFF3 annotation file](ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.gff3.gz).

**genome_annotation.txt**

This file contains gene function annotation information, which can be uploaded to iwa-miRNA server and used as inputs of the **manualCuration** and **associationAnalysis** modules. It has seven columns with tab-delimited format, which represent chromosomes, start sites, end sites, ID, strand, gene types, and detailed types. These gene attributes were extracted from [Arabidopsis GFF3 annotation file](ftp://ftp.ensemblgenomes.org/pub/plants/release-47/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.47.gff3.gz). TEs (Transposed elements) were downloaded from [TAIR10_Transposable_Elements](https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt).

**SNPs_in_Arabidopsis.txt**

This file contains part of SNPs from The 1001 Genomes Project and is used as input of the **sequenceVariation** function. It has five columns with tab-delimited format, representing chromosomes, start sites, ID, reference, and alternative alleles, respectively. This information can be downloaded from [Ensembl Plants Variations](https://plants.ensembl.org/biomart/martview).---
title: "The characteristics of annotated miRNAs"
output:
    flexdashboard::flex_dashboard:
        orientation: rows
        vertical_layout: scroll
        social: menu
        source_code: embed
        theme: cosmo
        self_contained: no
---

Results {data-icon=fa-area-chart}
=====================================

Row
-------------------------------------

```{r  include=FALSE}
library(flexdashboard)
library(knitr)
library(ggplot2)
library(dplyr)
library(plotly)
library(patchwork)
library(ggsci)

options(stringsAsFactors = F)
TE_inter <- read.csv("miRNA_TE.txt", sep = "\t", header = F)
TE_inter_names <- unique(TE_inter[,4])
non_TE_inter <- read.csv("miRNA_non_TE.txt", sep = "\t", header = F)
non_TE_inter <- non_TE_inter[!non_TE_inter[,4]%in%TE_inter_names, ]
non_TE_PCG_names <- non_TE_inter[!grepl("[Pp][CcEe][Gg]", non_TE_inter[,14]), 4]
PCG_index <- grepl("[Pp][CcEe][Gg]", non_TE_inter[,14])
non_TE_PCG_index <-  !(PCG_index&non_TE_inter[,4]%in%non_TE_PCG_names)
non_TE_inter <- non_TE_inter[non_TE_PCG_index, ]

new_table <- rbind(TE_inter, non_TE_inter)
all_mirs <- read.csv("miRNA_list.txt", sep = "\t", header = F)
all_mirs <- all_mirs[!all_mirs[,4]%in%new_table[,4], ]
all_mirs[,(ncol(all_mirs)+1):ncol(new_table)] <- "."
data <- rbind(new_table, all_mirs)

location_names <- apply(data, 1, function(x){paste0(c(x[1:3], x[6]), collapse = ":")})
uni_index <- vector()
for(i in unique(location_names)){
  uni_index <- c(uni_index, which(location_names%in%i)[1])
}

data <- data[uni_index,]
colnames(data)[c(5, 7, 14, 15)] <- c("Len", "TPM", "Type", "Subtype")
data[data[,14]==".", 14] <- "Intergenic"
data[data[,15]==".", 15] <- "Intergenic"
data[data[,15]=="LTR_retrotransposon",15] <- "LTR"
data[data[,15]=="solo_LTR",15] <- "LTR"
data[data[,15]=="terminal_inverted_repeat_element",15] <- "TIR"

theme_Publication <- function(base_size=14, base_family="sans") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               # legend.position = "bottom",
               # legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               # plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
       ))
      
}
```

### miRNAs (miRNA precusors)

```{r}
valueBox(length(data[,4]),
         icon = "fa-search")
```

### Genomic features

```{r}
valueBox(length(unique(data[,14])),
         icon = 'fa-star')
```

Row
-------------------------------------

### Percentage of miRNAs with different lengths

```{r}
data$Len <- factor(data$Len, levels=sort(unique(data$Len)))
mycolors <- c("blue", "#FFC125", "darkgreen", "darkorange")
p1 <- data %>%
         group_by(Len) %>%
         summarise(count = n()) %>%
         plot_ly(labels = ~sort(Len),
                 values = ~count,
                 marker = list(colors = mycolors),
                 sort = FALSE) %>%
         add_pie(hole = 0.5) %>%
         layout(xaxis = list(zeroline = F,
                             showline = F,
                             showticklabels = F,
                             showgrid = F),
                yaxis = list(zeroline = F,
                             showline = F,
                             showticklabels=F,
                             showgrid=F),
                  legend = list(title = list(text="Length")))
p1 %>%
  config(
    toImageButtonOptions = list(
      filename = "Length count",
      format = "svg",
      width = 200,
      height = 100
    )
  )
```

### Percentage of miRNA abundance with different lengths

```{r}
if(!any(data[,7]=="-")){
  p2 <- data %>%
           group_by(Len) %>%
           summarise(count = sum(TPM)) %>%
           plot_ly(labels = ~Len,
                   values = ~count,
                   marker = list(colors = mycolors),
                   sort = FALSE) %>%
           add_pie(hole = 0.5) %>%
           layout(xaxis = list(zeroline = F,
                               showline = F,
                               showticklabels = F,
                               showgrid = F),
                  yaxis = list(zeroline = F,
                               showline = F,
                               showticklabels=F,
                               showgrid=F),
                  legend = list(title = list(text="Length")))
  p2 %>%
    config(
      toImageButtonOptions = list(
        filename = "Total TPM of different lengths",
        format = "svg",
        width = 200,
        height = 100
      )
    )
}
```

Row
-------------------------------------

### Percentage of the length and abundance of miRNAs {data-width=350}

```{r}
if(!any(data[,7]=="-")){
  raw_df <- matrix(ncol = 3)
  breaks = c(0,100,1000,10000,floor(max(data[,7])))
  for(i in 1:length(breaks)-1){
    index_tmp <- data[,7]>=breaks[i]&data[,7]<breaks[i+1]
    index_tmp <- table(data[index_tmp,5])
    for(j in names(index_tmp)){
       raw_df <- rbind(raw_df, c(paste0(breaks[i], "-", breaks[i+1]),
                             j, as.numeric(index_tmp[j])))
    }
  }

  raw_df <- as.data.frame(raw_df[-1,])
  raw_df[,3] <- as.numeric(raw_df[,3])
  raw_df[,1] <- factor(raw_df[,1],levels=unique(raw_df[,1]))
  colnames(raw_df) <- c("Region", "Length", "Percentage")

  ggplot(raw_df, aes(Region, Percentage, fill=Length)) +
    geom_bar(stat="identity", position="fill") +
    theme_Publication() +
    theme(axis.line = element_line(colour = "black", size = 0.5)) +
    guides(fill=guide_legend(title=NULL)) + coord_flip() +
    scale_fill_brewer(palette = "Set3")
 }
```

### Distribution of miRNA length among different genomic features {data-width=350}

```{r}
data_count <- table(data[,c(5,14)])
data_count <- reshape2::melt(data_count)
colnames(data_count) <- c("Length", "Regions", "Count") 

p3 <- ggplot()+geom_bar(data = data, aes(x=Type, fill=Type)) +
  theme_Publication() + theme(axis.ticks.x = element_blank(), 
                          axis.title.x = element_blank(),
                          axis.line.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.line.y = element_line(colour = "black", size = 0.5))+
  scale_y_continuous(breaks=seq(0, max(table(data[,14])), 20)) +
  scale_fill_jama() + ylab("Count")
  

p4 <- ggplot() + 
  geom_point(data = data_count, 
             aes(x=Regions, y=Length, 
                 size = Count), color="#4B6370") +
  scale_size(range = c(0, 8)) +
  theme_Publication() +
  theme(panel.grid.major = element_line(colour = NA)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(colour = "black", size = 0.5))+
  scale_y_continuous(breaks = seq(min(data_count$Length),
                                  max(data_count$Length),
                                  1))

p5 <- ggplot()+geom_bar(data = data, aes(x=Len), width=0.5) + 
  coord_flip() + theme_Publication() +
  theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(),
        axis.line.y = element_blank(),axis.text.y = element_blank(),
        axis.line.x = element_line(colour = "black", size = 0.5))  + ylab("Count")

design <- "A#
           BC"
wrap_plots(A = p3, B = p4, C = p5, 
           design = design, heights = c(1.5,1), guides = "collect")
```

### Distribution of miRNA length among different TEs {data-width=350}

```{r}
if("TE" %in% data$Type){
  TE_data <- data[data$Type=="TE", ]
if(nrow(TE_data)>0){
  TE_data_mat <- table(TE_data[,c(5,15)])
  TE_data_mat <- reshape2::melt(TE_data_mat)
  colnames(TE_data_mat) <- c("Length", "Subtype", "Count") 
  
  p6 <- ggplot()+geom_bar(data = TE_data, aes(x=Subtype, fill=Subtype)) +
    theme_Publication() + theme(axis.ticks.x = element_blank(), 
                            axis.title.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.y = element_line(colour = "black", size = 0.5))+
    scale_y_continuous(breaks=seq(0, max(table(data[,14])), 20)) +
    scale_fill_lancet() + ylab("Count")
    
  
  p7 <- ggplot() + 
    geom_point(data = TE_data_mat, 
               aes(x=Subtype, y=Length, 
                   size = Count), color="#4B6370") +
    scale_size(range = c(0, 8)) +
    theme_Publication() +
    theme(panel.grid.major = element_line(colour = NA)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.line = element_line(colour = "black", size = 0.5))+
      scale_y_continuous(breaks = seq(min(TE_data_mat$Length),
                                    max(TE_data_mat$Length),
                                    1))
  
  p8 <- ggplot()+geom_bar(data = TE_data, aes(x=Len), width=0.5) + 
    coord_flip() + theme_Publication() +
    theme(axis.ticks.y = element_blank(), axis.title.y = element_blank(),
          axis.line.y = element_blank(),axis.text.y = element_blank(),
          axis.line.x = element_line(colour = "black", size = 0.5))  + ylab("Count")
  
  # subplot(plot_spacer(), p1, p2, p3,
  #         nrows = 2, margin = 0.04, heights = c(0.7, 0.3))
  
  design <- "A#
             BC"
  wrap_plots(A = p6, B = p7, C = p8, 
             design = design, heights = c(1.5,1), guides = "collect")

 }
}
```
---
title: "Raw sequencing data preprocessing"
date: '(`r format(Sys.time(), format = "%Y-%m-%d")`)'
output:
  flexdashboard::flex_dashboard:
    vertical_layout: scroll
    social: menu
    source_code: embed
    theme: cosmo
    self_contained: no
---

```{r setup, include = FALSE}
library(DT)
library(ggplot2)
library(purrr)
library(highcharter)
options(stringsAsFactors = F)
knitr::opts_chunk$set(echo = FALSE, 
                      message = FALSE,
                      warning = FALSE,
                      fig.align='center',fig.pos='H',
                      fig.path = "output_plots/",
                      dev = "svg")
```

<div id="preloader"></div>

```{js, echo=FALSE}
$(function() {
  $(window).load(function() {
    $('#preloader').fadeOut('slow',function(){$(this).remove();});
  });
});
```

```{css, echo=FALSE}
.limitrow {
  width: 99.4% !important;
  flex: none !important;
}

.limitrow .chart-stage {
  overflow-x:scroll;
}

#row-4 .chart-stage,
#row-5 .chart-stage {
width:1475px !important;
overflow-x:scroll !important;
}
```

Summary  {data-orientation=rows data-icon="fa-table"}
=====================================

### Fetch large-scale small RNA-seq data from NCBI and filter the raw data.
Users can collect public data from NCBI (https://www.ncbi.nlm.nih.gov/). Samples with low data quality will be excluded. The remaining samples were analyzed and indicated in the following tables and figures.

Row {data-width=850}
--------------------------------------

### Summary table {.limitrow}
```{r}
inputMat <- read.table("00Table_Summary_of_sRNA-seq_data.txt", stringsAsFactors = F, sep = "\t", row.names = 1)
inputMat <- inputMat[order(rownames(inputMat)), ]
adapter <- read.table("Adapter.txt", row.names = 1)
inputMat <- cbind(adapter[rownames(inputMat),], inputMat)
colnames(inputMat) <- c("Adapter", 'Phred', 'Barcode', "Raw reads", "Clean reads", "Ratio(%)", "Length limited", "Aligned", "Remove t/r/sn/snoRNA", "Collapsed reads","Most abundant read", "Most abundant read count")

for(tmp_nn in c(4,5,7,8,9,10,12)){
  inputMat[,tmp_nn] <- formatC(inputMat[,tmp_nn], big.mark=",", format="d")
}

datatable(data = inputMat, extensions = 'Buttons', class="compact cell-border", options = list( dom = "Blfrtip", buttons = list("copy",list(extend = "collection", buttons = c("csv", "excel", "pdf"), text = "Download") ),lengthMenu = list( c(10, 20, -1), c(10, 20, "All") ), pageLength = 10, autoWidth = TRUE), rownames = TRUE)%>%
formatStyle(TRUE, `text-align` = 'right')
```

<p style="font-size:14px;"> This table displays read number, adaptors, and barcode information for each library. `Adapter`: The adapter sequence; `Phred`: Phred score; `Barcode`: Whether there is a barcode in each library, "-" means no; `Raw reads`: Number of reads generated by sequencing; `Clean reads`: Number of reads that are filtered by adapters and quality scores; `Ratio`: Clean reads/Raw reads; `Length-limited`: Number of reads satisfying the length criteria; `Aligned`: Number of reads in limited length and at most one mismatch during mapping (-v 1 --best --strata); `Remove t/r/sn/snoRNA`: Number of mapped reads after removing reads corresponding to tRNAs, rRNAs, snRNAs and snoRNAs; `Collapsed reads`: Number of unique reads after reads collapsing; `Most abundant read`: The sequence with the highest count; `Most abundant read count`: The highest count. </p>

Row {data-width=850}
--------------------------------------

### Length distribution of all small RNAs

```{r}
Length_count <- read.table("Length_count.txt", sep = "\t")
Length_count <- Length_count[Length_count[,1]=="All", ]
Length_count <- Length_count[order(Length_count[,3]), ]
ds <- map(unique(Length_count[,2]), function(x){
    dt <- Length_count[Length_count[,2] == x, 3:4]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>% 
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution of collapsed small RNAs")

```

> The figure displays an interactive line chart which contains the length distribution of all reads in each small RNA-seq sample. The sample information is available by hovering over the points.

### Length distribution of collapsed small RNAs

```{r}
Length_count <- read.table("Length_count.txt", sep = "\t")
Length_count <- Length_count[Length_count[,1]=="Unique", ]
Length_count <- Length_count[order(Length_count[,3]), ]
ds <- map(unique(Length_count[,2]), function(x){
    dt <- Length_count[Length_count[,2] == x, 3:4]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>% 
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution of collapsed small RNAs")

```

> The figure displays an interactive line chart which contains the length distribution of unique reads in each small RNA-seq sample. The sample information is available by hovering over the points.
---
title: "Aggregating already annotated miRNAs"
date: '(`r format(Sys.time(), format = "%Y-%m-%d")`)'
output:
  flexdashboard::flex_dashboard:
    vertical_layout: scroll
    social: menu
    source_code: embed
    theme: cosmo
    self_contained: no
---

```{css}
.correct:before {
  content: '\2611';
  color: #008100;
  font-style:normal;
  font-weight:bold;
  font-size: 150%;

}
.incorrect:before {
  content: '\2612';
  font-style:normal;    
  color: #b20610;
  font-size: 150%;
}

 .zoom {
  height:400px;
	display:block;
  left: 50%;
  margin: auto;
  border: 0
} 

```

```{r setup, include = FALSE}
library(DT)
library(dplyr)
library(tibble)
library(formattable)
knitr::opts_chunk$set(echo = FALSE, 
                      message = FALSE,
                      warning = FALSE)

options(stringsAsFactors = F)
options(digits=3)
premirTab <- read.table("Translate_out.txt", sep = "\t", header = T)
premirTab <- premirTab[, -ncol(premirTab)]
premirTab <- cbind(premirTab, Delete="<button class=\"btnDelete\">Delete</button>")
linkname <- gsub(":", "_", premirTab[,2])
premirTab[,1] <- paste0("<a href=\"javascript:void(0)\" onclick=\"mirnaplot('",linkname,"')\">", premirTab[,1],"</a>")
dbsum <- sum(colnames(premirTab)%in%c("miRBase", "PmiREN", "sRNAanno", "Psgenes", "sRNA_Seq"))
dbloc <- ncol(premirTab)-dbsum-1

for(numt in (dbloc+1):(ncol(premirTab)-1)){
  premirTab[,numt] <- factor(premirTab[,numt], levels = unique(premirTab[,numt]))
}

premirTab <- premirTab[,c(1,4,dbloc:ncol(premirTab))]


```

Results  {data-orientation=columns data-icon="fa-area-chart"}
================================

There are three sections in this HTML report. In **Overview and Summary Information**, we list the aggregated miRNAs that identified by four databases and small RNA-seq data. By clicking on the link in the table, the **RNAfold** and **CentroidFold** sections interactively show the RNA structure of RNAfold and CentroidFold prediction.

Column 1 {data-width=400}
-------------------------------------

### Overview and Summary Information 

```{r outdata}
premirTab <- premirTab %>% rowid_to_column("Row") %>% mutate(Row = "")
premirTab$Mature_arm <- factor(premirTab$Mature_arm, levels = unique(premirTab$Mature_arm))

datatable(
    premirTab,
    elementId = "linktable",
    class="compact cell-border",
    filter = 'top',
    extensions = c("Select", "Buttons", 'ColReorder'),
    options = list(
        autoWidth = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
        columnDefs = list(list(className = "select-checkbox",targets = 0, orderable = FALSE),
                          list(width = '2px', targets = "_all")),
        select = list(style = "multi", selector = "td:first-child"),
        searchHighlight = TRUE,
        colReorder = FALSE,
        scrollX = TRUE,
        fixedColumns = FALSE,
        extensions = 'Responsive',
        pageLength = 20),
    escape = FALSE,
    rownames = FALSE,
    callback = JS("$('#linktable tbody').on('click','.btnDelete',function(){table.row( $(this).parents('tr') ).remove().draw();});"))%>%
  formatStyle(TRUE, `text-align` = 'center')
```
<p style="font-size:14px;"> The table displays a collection of miRNAs from public databases and small RNA-Seq data. `Precursors`: The name of miRNA precursors; `pLen`: The length of miRNA precusors; `miRBase/PmiREN/sRNAanno/Psgenes`: Whether miRNAs are included in these databases. TRUE means being included, and NO means no; `sRNA_Seq`: Whether miRNAs are predicted in small RNA-Seq data (sRNA-Seq). </p>


Column 2 {data-width=200}
--------------------------------------------

### RNAfold

```{r RNAfold}
library(htmltools)
tags$img(id="content1",class="zoom", src=knitr::image_uri(paste0("png/", linkname[1], "_r.png")))
```

### CentroidFold

```{r CentroidFold}
library(htmltools)
tags$img(id="content2",class="zoom", src=knitr::image_uri(paste0("png/", linkname[1], "_c.png")))
```

```{js}
function mirnaplot(event)
{
document.getElementById("content1").src= "png/" + event + "_r.png";
document.getElementById("content2").src= "png/" + event + "_c.png";
}
```
---
title: "Report page"
output:
    flexdashboard::flex_dashboard:
        orientation: rows
        vertical_layout: scroll
        social: menu
        source_code: embed
        theme: cosmo
        self_contained: no
---

```{r setup, include=FALSE}
library(DT)
library(purrr)
library(readr)
library(dplyr)
library(cytoscape)
library(kableExtra)
library(highcharter)
library(svgPanZoom)
options(stringsAsFactors = F)
knitr::opts_chunk$set(echo = FALSE,
                      message = FALSE,
                      warning = FALSE)
```

```{r import, include=FALSE}
info_table_raw <- read.table("final_table.txt", sep = "\t", header = T)
info_table <- info_table_raw[,c(2,5:ncol(info_table_raw))]
info_table <- info_table[info_table$ID!="", ]

mature_exp <- read.table("miRNA_in_sample.txt", sep = "\t", stringsAsFactors = F, row.names = 1, header = T)
rownames(mature_exp) <- gsub("T", "U", rownames(mature_exp))

tmp_fullname <- info_table$Extended_stem_loop_loc
tmp_fullname <- gsub(":", "-", tmp_fullname)
tmp_index <- which(tmp_fullname%in%"xxxxxx")

tmp_name <- info_table[tmp_index, 2]
info_vec <- info_table[tmp_index, c(1:13,20:21)]

tmp_df <- data.frame("Location"=as.character(info_vec[c(2,5,7,10)]), 
                     "Sequence"=as.character(info_vec[c(3,6,8,11)]), 
                     "Length"=as.character(c(info_vec[4],  nchar(info_vec[6]), 
                                             info_vec[c(9,12)])))
rownames(tmp_df) <- c("miRNA precursors", "Stem loop", "5p", "3p")
mature_df <- t(mature_exp[c(info_vec$Seq5p, info_vec$Seq3p), ])
rownames(mature_df) <- colnames(mature_exp)
mature_df[is.na(mature_df)] <- 0
arm_ratio <- apply(mature_df, 1, function(x){
  if(max(x[1:2])<3){
    return(0)
  } else{
    if(info_vec$Mature_arm == '5p'){
      return(round(log2((x[1]+1)/(x[2]+1)), 2))
    } else{
      return(round(log2((x[2]+1)/(x[1]+1)), 2))
    }
  }
})

```

```{css}
p.seq {
  word-wrap:break-word;word-break:break-all;
  font-size: 15px;
}

p.one{
  margin-top:0.2em;
}

table{
  table-layout:fixed;
}

table td {
  word-wrap:break-word;
  
}
```

Browse {data-icon=fa-area-chart}
=====================================

Row 
-------------------------------------

### Location and sequences {data-width=250}

```{r table1, results='asis'}
eseq <- tmp_df[1,2]
eseqms <- regexpr(tmp_df[3,2], eseq)[1]
eseqss <- regexpr(tmp_df[4,2], eseq)[1]

tmp_df[1,2] <- paste0('<p class="seq">', substr(eseq, 1, eseqms-1),'<span style="background-color: #F88017">',
       tmp_df[3,2], "</span>", substr(eseq, eseqms+nchar(tmp_df[3,2]), eseqss-1),
       '<span style="background-color: #168EF7">', tmp_df[4,2], "</span>", substr(eseq, eseqss+nchar(tmp_df[4,2]), nchar(eseq)), '</p>')

tmp_df[2,2] <- paste0('<p class="seq"><span style="background-color: #F88017">',
       tmp_df[3,2], "</span>", substr(eseq, eseqms+nchar(tmp_df[3,2]), eseqss-1),
       '<span style="background-color: #168EF7">', tmp_df[4,2], "</span></p>")

knitr::kable(tmp_df, escape = FALSE) %>%kable_styling("striped", full_width = FALSE)
```

## {data-height=10}

<p class="one"></p>

Row
-------------------------------------

### RNAfold structure {data-height=100, data-width=350}

```{r pressure1, echo=FALSE, fig.cap="RNAfold", out.width = '33%'}
svgPanZoom::svgPanZoom( read_file(paste0("../miRNASelection/data/",  gsub(":", "_", tmp_name), '_r.svg')))
```

### Information {data-width=250}

```{r}
fold_info <- t(read.table(paste0("../miRNASelection/data/", tmp_name, '_fold.txt')))
colnames(fold_info) <- c("RNAfold", "Centriodfold")
rownames(fold_info) <- c("miRNAs", "meet criterion", "mismatch+bugleOne", "bugleOne", "bugleTwo")
knitr::kable(fold_info, escape = FALSE, align = "c") %>% 
  kable_styling("striped", full_width = TRUE)
```

### Centroidfold structure {data-width=350}

```{r pressure2, echo=FALSE, fig.cap="Centroidfold", out.width = '33%'}
svgPanZoom::svgPanZoom( read_file(paste0("../miRNASelection/data/",  gsub(":", "_", tmp_name), '_c.svg')))
```

## {data-height=10}

Row
-----------------------------------

### Overview of read-stacks on extended miRNA precursor  {data-height=150}

<div style="height:100%;width:100%;overflow:auto">
```{r readstacks, class="scroll-100", results='asis'}
htmltools::includeHTML(paste0("../miRNASelection/data/", tmp_name, '_map.html'))
```
</div>

## {data-height=10}

<p class="one"></p>

Row
-------------------------------------

### mature miRNA expression

```{r}
tissue_info <- read.table("sample_info.txt", sep = "\t", stringsAsFactors = F, row.names = 1)
tissue_name <- tissue_info[rownames(mature_df),]
exp_df <- data.frame("miRNAs"=rep(c(paste0(info_vec[1], "-5p"), paste0(info_vec[1], "-3p")), each=nrow(mature_df)), "TPM" = c(mature_df[,1], mature_df[,2]), "Tissues"= rep(tissue_name, 2), stringsAsFactors = F)
exp_df <- exp_df[exp_df[,3]!="-", ]
exp_df[,1] <- factor(exp_df[,1], levels = unique(exp_df[,1]))

hcboxplot(x = exp_df$TPM, var = exp_df$Tissues, var2 = exp_df$miRNAs, outliers = FALSE) %>%
  hc_chart(type = "column") %>%
  hc_xAxis(title = list(text = "")) %>%
  hc_yAxis(title  = list(text = "TPM")) %>%
  hc_exporting(enabled = TRUE, filename = "mature_miRNA_expression")
```

### Arm switch events

```{r}
mature_df_out <- data.frame('5p'=mature_df[,1], 
                            '3p'=mature_df[,2], "ratio"= arm_ratio, "grouplist"=tissue_info[rownames(mature_df),], stringsAsFactors = F)
mature_df_out <- mature_df_out[mature_df_out[,4]!="-", ]
mature_df_out <- mature_df_out[order(mature_df_out[["ratio"]], decreasing = T), ]
mature_df_out <- mature_df_out[order(mature_df_out[,3], decreasing=T), ]
mature_df_out <- mature_df_out[order(mature_df_out[,4]), ]

hchart(
  mature_df_out, 
  "scatter",
  hcaes(x=1:nrow(mature_df_out), y=ratio, group = grouplist, radius = 0), radius = 10
  ) %>% 
  hc_xAxis(title = list(text = "")) %>% 
  hc_yAxis(title  = list(text = "Log2(miR/miR*)")) %>% 
  hc_tooltip(
    pointFormat = "{series.name}: <b>{point.y}</b><br/>", 
    shared = TRUE,
    valueSuffix = "", 
    crosshairs = TRUE
  ) %>% 
  hc_add_theme(hc_theme_flat(chart = list(backgroundColor = "#FFF")))%>% 
  hc_exporting(enabled = TRUE, filename = "Arm_switch_events")
```

## {data-height=10}

<p class="one"></p>

Row {.tabset}
--------------------------------------------

### Network of miRNAs-target interactions {data-width=400}

```{r}
tmp_path <- paste0("../miRNASelection/data/", tmp_name, '.mti')
if(file.exists(tmp_path)&file.info(tmp_path)[1,1]>0){
    inputRaw <- read.table(tmp_path, sep = "\t", header = F, stringsAsFactors = F)
    if(nrow(inputRaw)>0){
      inputRaw[,2] <- gsub("\\.\\d", "", inputRaw[,2])
      if(file.exists("gene_description.txt")){
        gene_description <- read.table("gene_description.txt", row.names = 1, stringsAsFactors = F, sep = "\t")
        inputRaw[,4] <- gene_description[inputRaw[,2], 1]
      }else{
        inputRaw[,4] <- "-"
      }
      inputRaw[,1] <- gsub(tmp_name, info_vec$ID, inputRaw[,1])
      inputRaw[,1] <- gsub("\\+", "", inputRaw[,1])
      inputRaw <- inputRaw[, 1:12]
      colnames(inputRaw) <- c("miRNA_Acc.", "Target_Acc.", "Expectation", "Target description", "miRNA_start", "miRNA_end", "Target_start", "Target_end", "miRNA_aligned_fragment", "alignment", "Target_aligned_fragment", "Inhibition")
      inputRe <- unique(inputRaw[,1:2])
      nodes <- data.frame(id = unique(c(inputRe[,1],inputRe[,2])))%>%
            mutate(node_color = ifelse(id %in% inputRe[,1],
                                   "#DE3025",
                                   "#85C7E6"),
                   node_shape = ifelse(id %in% inputRe[,1],
                                   "rectangle",
                                   "ellipse"),
                   node_font_size = ifelse(id %in% inputRe[,1],
                                   10,5),
                   node_height = ifelse(id %in% inputRe[,1],
                                   20,15))
      edges <- data.frame(id = paste0(inputRe[,1], "_", inputRe[,2]), source = inputRe[,1], target = inputRe[,2])
      
      cytoscape(nodes = nodes, edges = edges, elementId="networkfig")%>% 
        layout(name = 'cose',
               directed = TRUE,
               padding = 4, 
               avoidOverlapPadding = 30) %>%
        node_style('height' = 'data(node_height)',
                   'width' = 'data(node_height)',
                   'font-size' = 'data(node_font_size)',
                   'background-fit' = 'cover',
                   'border-color' = '#000',
                   'border-width' = 1,
                   'border-opacity' = 0.5,
                   'background-color' = 'data(node_color)',
                   'shape' = 'data(node_shape)') %>%
        edge_style('curve-style' = 'bezier',
                   'width' = 1,
                   'arrow-scale' = 0.8,
                   'target-arrow-shape' = 'triangle',
                   'line-color' = '#88898D',
                   'target-arrow-color' = '#88898D')%>%
          panzoom()
    }
}

```

### Information of Network {data-width=600}

```{r}
if(file.exists(tmp_path)&file.info(tmp_path)[1,1]>0){
  table_options <- function() {
    list(dom = 'Bfrtip',
      pageLength = 8,
      buttons = list(c('copy', 'csv', 'excel', 'pdf')),
      searchHighlight = TRUE,
      colReorder = TRUE,
      scrollX = TRUE,
      fixedColumns = TRUE,
      extensions = 'Responsive',
      deferRender = TRUE,
      scroller = TRUE,
      lengthChange = FALSE
      )
  }
  
  datatable(
    inputRaw,
    rownames = FALSE,
    editable = TRUE,
    elementId = "linktable",
    class = 'cell-border',
    escape = FALSE,
    options = table_options(),
    extensions = c('Buttons', 'Select')) %>% 
    formatStyle(columns = colnames(inputRaw), fontSize = '80%')
}
```

## {data-height=10}

<p class="one"></p>

Row
-----------------------------------

### Secondary structures by strucVis {data-width=350}

```{r pressure3, echo=FALSE, fig.cap="strucVis", out.width = '33%'}
svgPanZoom::svgPanZoom( read_file(paste0("../miRNASelection/data/", tmp_name, '.svg')))
```

### The expression level of different types of isomiR {data-width=500}

```{r}
tmp_out <- paste0("../miRNASelection/data/", tmp_name, '.out')
if(file.exists(tmp_out)&file.info(tmp_out)[1,1]>0){
  iso_list <- read.table(tmp_out, sep = "\t", stringsAsFactors = F)
  iso_list <- iso_list[iso_list[,10]!="-", ]
  iso_list[,2] <- factor(iso_list[,2], levels = c("5p", "3p"))
  iso_list[,10] <- factor(iso_list[,10], levels = c("ref", "add5", "sub5", "add3", "sub3",
                                        "add5_add3", "add5_sub3", "sub5_add3", "sub5_sub3",
                                        "seed_snp", "tail_snp"))
  
  hchart(
    iso_list, 
    "column",
    hcaes(x=V10, y=V8, group = V2), color=c("#DE3025", "#257ADE")[1:length(unique(iso_list$V2))]
    ) %>% 
    hc_xAxis(title = list(text = "")) %>% 
    hc_yAxis(title  = list(text = "TPM")) %>% 
    hc_tooltip(
      pointFormat = "{series.name}: <b>{point.y}</b><br/>", 
      shared = TRUE,
      valueSuffix = "", 
      crosshairs = TRUE
    ) %>% 
    hc_add_theme(hc_theme_flat(chart = list(backgroundColor = "#FFF")))%>% 
    hc_exporting(enabled = TRUE, filename = "different_types_of_isomiR")
}
```

---
title: "Summary page"
date: '(`r format(Sys.time(), format = "%Y-%m-%d")`)'
output:
  flexdashboard::flex_dashboard:
    vertical_layout: scroll
    source_code: embed
    theme: cosmo
    self_contained: no
---

Quality control  {data-orientation=rows data-icon="fa-area-chart"}
================================
  
Row
-------------------------------------

### Overview and Summary Information

```{r setup, include = FALSE}
library(DT)
library(dplyr)
library(tibble)
library(formattable)
knitr::opts_chunk$set(echo = FALSE, 
                      message = FALSE,
                      warning = FALSE)
options(stringsAsFactors = F)
```

```{r }
options(stringsAsFactors = F)
info_table <- read.table("final_table_source.txt", sep = "\t", header = T)
info_table[info_table[, "TPM5p"]=="-", "TPM5p"] <- "0"
info_table[info_table[, "TPM3p"]=="-", "TPM3p"] <- "0"
info_table$Mean <- round(info_table$Mean, 2)
info_table$Max <- round(info_table$Max, 2)

colname_list <- c("ID", "HTcriteria", "One_class_SVM", "Genomic_source","Source","Stem_loop_loc", "Stem_loop_len", "Stem_loop_MFE", "Stem_loop_AMFE",
                     "Mature_arm", "Seq5p", "Len5p", "TPM5p", "The_number_of_sequences_in_pre.miRNAs", "Abundance_bias", "Strand_bias",
                     "RNAfold", "Centroidfold", "Mean", "Max", "Samples")
df <- info_table[, colname_list]

replaceCol <- c("Mature_arm", "Seq5p", "Len5p", "TPM5p")
for(i in 1:nrow(df)){
  if(info_table[i, "Mature_arm"]=="5p"){
    df[i, replaceCol] <- info_table[i, replaceCol]
  } else if(info_table[i, "Mature_arm"]=="3p"){
    df[i, replaceCol] <- info_table[i, c("Mature_arm", "Seq3p", "Len3p", "TPM3p")]
  } else {
    if(as.numeric(info_table[i, "TPM5p"])>=as.numeric(info_table[i, "TPM3p"])){
      df[i, replaceCol] <- c('5p', info_table[i, c("Seq5p", "Len5p", "TPM5p")])
    }else{
      df[i, replaceCol] <- c('3p', info_table[i, c("Seq3p", "Len3p", "TPM3p")])
    }
  }
}


# df <- read.table("summary_information.txt", sep = "\t", header = T)
#tissue_name <- colnames(info_table)[33:ncol(info_table)]
linkname <- gsub(":", "-", info_table$Extended_stem_loop_loc)
# for(i in 1:nrow(df)){
#  df_out <- data.frame("name" = colnames(df),  "value" = as.character(df[i,]))
#   writeLines(jsonlite::toJSON(df_out), paste0("data/", linkname[i], ".json"))
# }
df[,1] <- paste0("<a href=\"javascript:void(0)\" onclick=\"mirnaplot('",linkname,"')\">", df[,1],"</a>")

# set_list <- lapply(1:3,function(cc){
#   formatter("span", style = x ~ style(color = ifelse(x, "green", "red")),
#             x ~ icontext(ifelse(x, "ok", "remove"), ifelse(x, "Yes", "No")))
# })
# names(set_list) <- colnames(df)[c(2,16,17)]

df <- cbind(df, "remove" = '<input type="button" value="Delete" class="btnDelete">')
df <- df %>% rowid_to_column("Row") %>% mutate(Row = "")
# df <- formattable(df, set_list)

th_style = "padding: 5px;
color: #fff;
background-color:#517FB9;
text-align: center;
border-right-width: 1px; 
border-right-style: solid; 
border-right-color: white; 
border-bottom-width: 1px; 
border-bottom-style: solid; 
border-bottom-color: white;
word-wrap: break-word;      
overflow-wrap: break-word;"

table_frame <-
  function() {
    htmltools::withTags(
      table(class = 'display', style = "padding: 1px; font-size: 0.8em; 
            font-family: sans-serif; text-align: center; word-wrap: break-word; overflow-wrap: break-word;",
            thead(tr(
              th(rowspan = 2, style = th_style, 'Row'),
              th(rowspan = 2, style = th_style, 'ID'),
              th(rowspan = 2, style = th_style, 'HT criteria'),
              th(rowspan = 2, style = th_style, 'One class SVM'),
              th(rowspan = 2, style = th_style, 'Genomic source'),
              th(rowspan = 2, style = th_style, 'Source'),
              th(class = 'dt-center', style = th_style, colspan = 4, 'Stem loop'),
              th(class = 'dt-center', style = th_style, colspan = 4, 'Mature miRNA'),
              th(class = 'dt-center', style = th_style, colspan = 5, 'miRNA precursor'),
              th(class = 'dt-center', style = th_style, colspan = 3, 'Expressed samples'),
              th(rowspan = 2, style = th_style, 'Remove'),
              tr(lapply(c('Loc', 'Len', 'MFE', 'AMFE', 'Arm', 'Seq', 'Len', 'TPM',
                  'Seq count', 'Abundance bias', 'Strand bias',
                  'RNAfold', 'Centroidfold','Mean','Max','Sample(TPM>1)'), th, style = th_style)) #'Mean', 'Max','Sample (TPM>1)', tissue_name
            )
          )))
  }

#th(class = 'dt-center', style = th_style, colspan = 3, 'TPM in 1063 samples'), 
#th(class = 'dt-center', style = th_style, colspan = length(tissue_name), 'TPM in tissues'), 

table_options <- function() {
  list(dom = 'Bfrtip',
    pageLength = 20,
    buttons = list(c('copy', 'csv', 'excel', 'pdf', 'print')),
    columnDefs = list(list(className = "select-checkbox", targets = 0, orderable = FALSE)),
    select = list(style = "multi", selector = "td:first-child"),
    searchHighlight = TRUE,
    scrollX = TRUE,
    fixedColumns = TRUE,
    extensions = 'Responsive',
    deferRender = TRUE,
    scroller = TRUE,
    lengthChange = FALSE,
    initComplete = JS(
      "function(settings, json) {",
      "$(this.api().table().header()).css({'background-color': '#517fb9', 'color': '#fff'});",
      "}"
    )
      )
}

#brks <- quantile(df[,tissue_name], probs = seq(.05, .95, .05), na.rm = TRUE)
#clrs <- round(seq(255, 40, length.out = length(brks) + 1), 0) %>%
#{paste0("rgb(255,", ., ",", ., ")")}

df$One_class_SVM <- factor(df$One_class_SVM, levels = unique(df$One_class_SVM))
df$Source <- factor(df$Source, levels = unique(df$Source))
df$Genomic_source <- factor(df$Genomic_source, levels = unique(df$Genomic_source))

datatable(
  df,
  rownames = FALSE,
  editable = FALSE,
  elementId = "linktable",
  class="compact cell-border",
  filter = 'top',
  # class = 'cell-border',
  escape = FALSE,
  container = table_frame(),
  options = table_options(),
  extensions = c('Buttons', 'Select'),
  callback = JS("$('#linktable tbody').on('click','.btnDelete',function(){table.row( $(this).parents('tr') ).remove().draw();});")) %>% 
  formatStyle("Genomic_source", backgroundColor = styleEqual(unique(df$Genomic_source), c("lightblue", "lightgreen", "lightpink", "lightgrey", "lightred", "lightyellow")[1:length(unique(df$Genomic_source))]))
```

> The table displays the overview information of Precursor miRNAs. Each row presents one candidate pre-miRNAs with name, location in genome, strand information, abundance, and bias.

```{js}
function mirnaplot(event){
  window.open("miRNA_out/" + event +".html");
}
```
---
title: "The available miRNAs in several public databases"
output:
    flexdashboard::flex_dashboard:
        orientation: rows
        vertical_layout: scroll
        social: menu
        source_code: embed
        theme: cosmo
---

```{r setup, include=FALSE}
library(DT)
library(purrr)
library(ggplot2)
library(dplyr)
library(highcharter)
options(stringsAsFactors = F)
knitr::opts_chunk$set(echo = FALSE,
                      message = FALSE,
                      warning = FALSE)
```

```{css, echo=FALSE}
.limitrow {
  width: 99.4% !important;
  flex: none !important;
}

.limitrow .chart-stage {
  overflow-x:scroll;
}

p.seq {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    width: 200px;
    min-width: 0;
    margin: 0;
    font-size: 15px;
}

p.seq:hover {
    overflow: none;
    width: auto;
}
```

miRBase {data-icon=fa-bar-chart}
=====================================

Row
-------------------------------------

### Table: Detailed information of miRNAs {.limitrow}

```{r }
sRNA_data <- read.table("miRBase.txt", sep = "	", header = T)
sRNA_data <- sRNA_data[sRNA_data[,8]<30&sRNA_data[,11]<30, ]
sRNA_data[,3] <- paste0("<p class=\"seq\">", sRNA_data[,3],"</p>")
datatable(data = sRNA_data,
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA categories

```{r}
name_col_two <- as.data.frame(table(sRNA_data[,12]))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")

```

### miRNA family sizes and total number of family

```{r}
name_family <- gsub("\\D$", "", sRNA_data[,1])
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(sRNA_data[,12]), function(x){
    dt <- density(sRNA_data[,4][sRNA_data[,12] == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data)),
                      "Len" = c(sRNA_data[,8], sRNA_data[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Composition of the first base

```{r}
mat_char <- data.frame('Type'=rep(c("5p","3p"), each = nrow(sRNA_data)),
                      "Letter" = c(substr(sRNA_data[,7], 1, 1), substr(sRNA_data[,10], 1, 1)))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```

PmiREN {data-icon=fa-bar-chart}
=====================================

Row
-------------------------------------

### Table: Detailed information of miRNAs {.limitrow}

```{r }
sRNA_data2 <- read.table("PmiREN.txt", sep = "	", header = T)
sRNA_data2[,3] <- paste0("<p class=\"seq\">", sRNA_data2[,3],"</p>")
datatable(data = sRNA_data2,
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA family sizes and total number of family

```{r}
name_family <- gsub("\\D$", "", sRNA_data2[,1])
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

### miRNA categories

```{r}
def_col_two <- rep("Known", nrow(sRNA_data2))
def_col_two[grepl("MIRN", sRNA_data2[,1])] <- "Novel"
name_col_two <- as.data.frame(table(def_col_two))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = def_col_two, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(def_col_two), function(x){
    dt <- density(sRNA_data2[,4][def_col_two == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data2)),
                      "Len" = c(sRNA_data2[,8], sRNA_data2[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Length and distribution of mature miRNAs

```{r}
mat_len <- data.frame('Ter'=sRNA_data2[,5],
                      "Len" = sRNA_data2[,8])
mat_len <- mat_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length and distribution of mature miRNA")
```

sRNAanno {data-icon=fa-bar-chart}
=====================================

Row
-------------------------------------

### Table: Detailed information of miRNAs {.limitrow}

```{r }
sRNA_data3 <- read.table("sRNAanno.txt", sep = "	", header = T)
sRNA_data3[,3] <- paste0("<p class=\"seq\">", sRNA_data3[,3],"</p>")
datatable(data = sRNA_data3,
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA family sizes and total number of family

```{r}
name_col <- t(do.call("cbind", strsplit(sRNA_data3[,1], "-")))
name_col[,1] <- gsub("\\D$", "", name_col[,1] )
colnames(name_col) <- c("mirna", 'Type')

name_col_pre <- table(table(name_col[,1]))
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])

highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_col[,1]%in%names(table(name_col[,1])[table(name_col[,1])>1])
family_show <- as.data.frame(name_col[family_index, ])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

### miRNA categories

```{r}
name_col_two <- as.data.frame(table(name_col[,2]))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(name_col[,2]), function(x){
    dt <- density(sRNA_data3[,4][name_col[,2] == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())

ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data3)),
                      "Len" = c(sRNA_data3[,8], sRNA_data3[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Length and distribution of mature miRNA

```{r}
mat_len <- data.frame('Ter'=sRNA_data3[,5],
                      "Len" = sRNA_data3[,8])
mat_len <- mat_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length and distribution of mature miRNA")
```

### Composition of the first base

```{r}
mat_char <- data.frame('Type'=name_col[,2],
                      "Letter" = substr(sRNA_data3[,7], 1, 1))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```

Plant small RNA genes {data-icon=fa-bar-chart}
=====================================

Row
-------------------------------------

### Table: Detailed information of miRNAs {.limitrow}

```{r }
sRNA_data4 <- read.table("PlantsmallRNAgenes.txt", sep = "	", header = T)
data_index <- sRNA_data4$Mature_arm!='unmatchedRegion'&sRNA_data4$Length5p<30&sRNA_data4$Length3p<30
sRNA_data4 <- sRNA_data4[data_index, ]
sRNA_data4[,3] <- paste0("<p class=\"seq\">", sRNA_data4[,3],"</p>")
datatable(data = sRNA_data4,
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA categories

```{r}
name_col_two <- as.data.frame(table(sRNA_data4$Source))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

### miRNA family sizes and total number of family

```{r}
name_family <- gsub("\\D$", "", sRNA_data4$Precursors)
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
if( "Freq"%in%colnames(name_col_one) ){
  highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
}
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(sRNA_data4$Source), function(x){
    dt <- density(sRNA_data4$pLength[sRNA_data4$Source == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

db <- map(unique(sRNA_data4$Type), function(x){
    dt <- density(sRNA_data4$pLength[sRNA_data4$Type == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_series_list(db)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data4)),
                      "Len" = c(sRNA_data4$Length5p, sRNA_data4$Length3p))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Composition of the first base

```{r}
firstbase <- apply(sRNA_data4, 1, function(x){
  if(x[5] == "5p"){
  return(substr(x[8], 1, 1))
  }else{
  return(substr(x[11], 1, 1))
  }
})
mat_char <- data.frame('Type'=rep(c(sRNA_data4$Source, sRNA_data4$Type)),
                      "Letter" = c(firstbase, firstbase))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```
