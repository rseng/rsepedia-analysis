---
title: '  CPgeneProfiler: A lightweight R package to profile the Carbapenamase genes
  from genome assemblies'
tags:
- R
- Carbapenamase gene profile
- cocarriage
- Beta-lactamases
- Antimicrobial resistance (AMR)
date: "09 Oct 2020"
output:
  rticles::joss_article: default
  pdf_document: default
  word_document: default
authors:
- name: Prakki Sai Rama Sridatta
  orcid: 0000-0002-9254-2557
  affiliation: '1,2'
- name: Natascha M Thevasagayam
  affiliation: '1,2'
- name: Weizhen Xu
  affiliation: '1,2'
- name: Kalisvar Marimuthu
  affiliation: '1,2'
- name: Jeanette W P Teo
  affiliation: 3
- name: Indumathi Venkatachalam
  affiliation: 4
- name: Ng Oon Tek
  affiliation: '1,2,5'
year: 2020
bibliography: paper.bib
citation_author: Prakki et. al.
csl: apa.csl
journal: JOSS
affiliations:
- name: National Centre for Infectious Diseases, Singapore
  index: 1
- name: Tan Tock Seng Hospital, Singapore
  index: 2
- name: National University Hospital, Singapore
  index: 3
- name: Singapore General Hospital, Singapore
  index: 4
- name: Lee Kong Chian School of Medicine, Nanyang Technological University, Singapore
  index: 5
---

# Summary

“Carbapenems” are a specific subset of antibiotics considered to possess a higher spectrum of antimicrobial activity [@papp2011carbapenems] against Gram-positive and Gram-negative bacteria. Even so, there are pathogens which are resistant to carbapenems due to the presence of carbapenemase genes (CP genes) which have the ability to hydrolyze carbapenems. 

Studies show that those cases infected by carbapenem-resistant pathogens have a higher morbidity and mortality rate compared with those who are infected by non-carbapenem-resistant pathogens [@van2013carbapenem; @cai2017prevalence]. Therefore, early discerning of the CP genes and their resistance mechanisms are considered crucial to aid in infection control as well as lessen the likelihood of mortality, duration of hospitalization stay, and related medical costs [@van2013carbapenem; @nordmann2019epidemiology]. Further, it is understood that the cocarriage of genes encoding different classes of carbapenemases could confer higher resistance to carbapenem antibiotics, which may promote further spread of the disease [@wang2019cocarriage].

The detection of the resistance genes from various bacterial strains using techniques such as polymerase chain reaction (PCR) and microarrays in real-time is very time consuming and costly. With the advancement in whole-genome sequencing (WGS) technologies and decreased costs, this is more accessible and WGS provides an alternative method for detection of resistance genes, given that the relevant analysis tools are available.

To this end, several freely available bioinformatics tools such as ``ABRicate`` (https://github.com/tseemann/abricate), ``AMRPlusPlus`` [@doster2020megares], ``ARG-ANNOT`` [@gupta2014arg], ``ARIBA`` [@hunt2017ariba], ``Comprehensive Antibiotic Resistance Database – Resistance Gene Identifier`` (CARD-RGI) [@alcock2020card], ``NCBI AMRFinderPlus`` (https://ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/), ``KmerResistance`` [@clausen2016benchmarking; @clausen2018rapid], ``PointFinder`` [@zankari2017pointfinder], ``Resfinder`` [@bortolaia2020resfinder], ``sraX`` [@panunzi2020srax], and ``SRST2`` [@inouye2014srst2] assist in finding the antimicrobial resistance genes from the sequence data [@hendriksen2019using].  

# Statement of Need 

Undeniably, all the above-mentioned tools are focused around the antimicrobial-resistant genes, and tools such as ABRicate and CARD-RGI can even generate comparative tables across genomes, and sraX can help in visualization of comprehensive AMR gene complement. Nevertheless, they do not readily generate a genetic profile for the presence of CP genes, and extract and visualize the set intersections of cocarriage of CP genes. Achieving this currently necessitates a restructuring and transformation of the output from these tools. Furthermore, in the research settings where it is crucial to quickly examine the transmission of CP genes, it is useful to have a tool that is catered to the CP gene dataset that provides easily interpretable visualizations and statistics. Therefore, to address this need, we describe here a lightweight R package, **CPgeneProfiler**, that scans multiple bacterial genome assemblies to detect and visualize the presence of CP genes and their cocarriage using the R framework. Additionally, this package also allows one to assess the size of CP contigs to check if the CP genes are distributed on the particular sequence size by generating the contig length distribution plots. 

# Implementation

In order to detect CP genes from the genome assemblies, ``NCBI Bacterial Antimicrobial Resistance Reference Gene Database`` (2020-07-16.2) (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047) was used for generating a CP gene database. Only those resistant genes whose subclass is categorized as "CARBAPENEM" in the reference gene catalog were considered for database preparation. This excluded the possibility of having resistant gene variants which are beta-lactamases but do not show carbapenem-resistant activity. For example, although OXA-48 is a carbapenemase gene, other OXA variants such as OXA-163 and OXA-405 have been determined to be devoid of any carbapenemase activity [@dortet2017noncarbapenemase], and therefore their subclass is not categorized as "CARBAPENEM" in the NCBI Bacterial Antimicrobial Resistance Reference Gene Database. Therefore, OXA variants OXA-163 and OXA-405 were not included in the CP gene database.

The tool first uses the ``cpblast`` command, by which each fasta file is searched against the CP gene database using NCBI BLAST+ [@camacho2009blast] (version 2.9.0+), which is pre-installed in the local system as a dependency. The presence of a CP gene in an assembled genome is confirmed if the CP gene meets the identity and coverage thresholds (default: 100%) when aligned with the genome sequence. The genome sequences that meet the thresholds are extracted from the BLAST results using the ``filt_blast`` command.

Visualizing the presence of CP genes and their corresponding counts across all the genome assemblies in a simple heatmap enables one to find CP gene variants that are found across the samples and aids in exploring the pattern of CP gene presence with reference to species or sequence type (ST). In order to facilitate this, the ``cpprofile`` command generates a profile of CP genes (Figure 1A) across the genome assemblies, while the ``cocarriage`` command finds cocarriage of CP genes in the genome assemblies. In addition to this, the tool also generates plots to visualize the set intersections of CP genes across all the input genome assemblies using the command ``upsetR_plot`` (Figure 1B). It is understood that isolates that harbour multiple carbapenemase genes are considered to produce high resistant phenotypes, and running the commands ``cocarriage`` and ``upsetR_plot`` provides an overview of the CP genes as well as their cocarriages present in all the genomes.

Given a set of bacterial genomes that are of same species, it would be useful to explore if the CP genes are found on specific plasmids or scattered across multiple plasmids/chromosomes of different sequence lengths. This can be achieved by plotting the number of contigs across the contig length by using the ``plot_conlen`` command (Figure 2). 

Lastly, ``CPgeneProfiler`` can also generate the N50, N90, and assembly size statistics for each of the genome assemblies and also plots the assembly size against N50 and N90 using the ``assembly_stat`` command (Figure 3A, 3B). This helps in quickly assessing and comparing the quality of the assembled genomes provided as an input. All the generated output files from various commands of the package are arranged accordingly into their respective folders using the ``cp_summarize`` command.

# Availability and Implementation

The R package ``CPgeneProfiler`` (version 2.1.1) is supported on UNIX/Linux machines. The source code, guide and datasets are currently available on Github repository (https://github.com/ramadatta/CPgeneProfiler).
<br>  
<br>  

## Step 1: Download CP gene database using R console

```r
# Specify CP gene database URL 
url <- "https://raw.githubusercontent.com/ramadatta/CPgeneProfiler/
master/testData/db/NCBI_BARRGD_CPG_DB.fasta"

# Specify destination where CP gene database file should be saved 
path <- "/home/user/db" # Can be changed to preferred location
setwd(path)
destfile <- "NCBI_BARRGD_CPG_DB.fasta"

# Download the CP gene database file to the folder set in "path"
download.file(url, destfile)
```

## Step 2: Install CPgeneProfiler package

The R package ``CPgeneProfiler`` can be installed by typing the following in R:
```r
devtools::install_github("ramadatta/CPgeneProfiler")
```

# Figures

![A) CP Gene Profile B) Set Intersection plot of the available CP genes across genome assemblies](cpprof_upset.png)
<br>  
<br>  
<br> 
**Figure 1. (A)** CP gene profile obtained by ‘cpprofile’ command **(B)** Set intersection plot of the available CP genes across genome assemblies, obtained by the ‘upsetR_plot’ command
 
 
![CP gene contig length distribution plots](len1_2.png)
**Figure 2.** CP gene contig length distribution plots obtained by the ‘plot_conlen’ command.
<br>  
<br>
![A) Assembly Size vs N50 B) Assembly Size vs N90 plots](as1_2.png)
**Figure 3.** Plots generated by the ‘assembly_stat’ command **(A)**  Assembly size vs N50 **(B)** Assembly size vs N90 

# Conclusion

``CPgeneProfiler`` can be used to understand the CP gene profile of a set of bacterial genome assemblies. It generates a simple heatmap for visualization of the CP gene profile and reports details on cocarriage of CP genes within the genomes. The capability to identify and visualize the presence of CP genes across multiple genomes would have useful applications, for example, in a dataset of outbreak samples, and the CPgeneProfiler could aid researchers in obtaining an overview of the samples and their CP gene carriage.

# Acknowledgement
The authors would like to thank Victor Ong and Wang Liang De for generating the sequence data that was used for developing and testing the tool. 

# Funding
This work is supported by the Singapore Ministry of Health’s National Medical Research Council under its NMRC Collaborative Grant: Collaborative Solutions Targeting Antimicrobial Resistance Threats in Health Systems (CoSTAR-HS) (NMRC CGAug16C005) and NMRC Clinician Scientist Award (MOH-000276).  Any opinions, findings and conclusions or recommendations expressed in this material are those of the author(s) and do not reflect the views of MOH/NMRC.


# References  

## Contributor Code of Conduct

As the maintainers of ``CPgeneProfiler``, we pledge to respect all people who want to contribute to the package through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities. In addition, we pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.

Project maintainers will retain the right and responsibility to remove, edit, or reject comments, commits, code, issues, and other contributions that are not aligned to this Code of Conduct. Instances of unacceptable behavior may be reported by opening an issue, emailing the package developer (Prakki).

This Code of Conduct is adapted from the Contributor Covenant, version 2.0, available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla’s code of conduct enforcement ladder](https://github.com/mozilla/diversity).
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02473/status.svg)](https://doi.org/10.21105/joss.02473)

# CPgeneProfiler
Generate a profile of carbapenamase genes from the genome assemblies

## Content
  * [Synopsis](#synopsis)
  * [Input Requirements](#input-requirements)
  * [Installation](#installation)
  * [Example Usage](#example-usage)
  * [Working with testData](#working-with-testdata)
  * [Output Files](#output-files)
  * [Version](#version)
  * [A few example plots](#a-few-example-plots)
  * [Community Guidelines ](#community-guidelines)
  * [Bug Report and Support Request](#bug-report-and-support-request)


##### **Synopsis**

Current AMR detection tools generate comparative tables across genomes and can help in visualization of comprehensive AMR gene complement. Nevertheless, they do not readily generate a genetic profile for the presence of CP genes, extract and visualize the cocarriage of CP genes. Achieving this currently necessitates a restructuring and transformation of the output from these tools. 

To address this need, we describe here a lightweight R package, CPgeneProfiler that scans multiple bacterial genome assemblies to detect and visualize the presence of CP genes and their cocarriage using the R framework. Input required is a directory of FASTA file with genome assemblies.

Additionally, this package also allows to assess the size of CP contigs to check if the CP genes are distributed on the particular sequence size by generating the contig length distribution plots.

Other assembly statistics such as N50, N90, Assembly Size from the assembly are calculated and plots of length distribution of CP gene contigs from
 the list of assemblies are reported.  
 
Currently the package works only on Unix systems.
 
 
##### **Input Requirements**
* Path of a directory with multiple FASTA files (can be in multiple contigs) 
* Path of Carbapenemase Gene Database file (FASTA) directory

#### Conda (https://anaconda.org/bioconda/cpgeneprofiler) 

To install this package with conda run:
``` 
conda install -c bioconda cpgeneprofiler
```

##### **Requirements**

- **R packages (REQUIRED):**
	 tidyverse,
	 UpSetR,
	 scales,
	 ape,
	 BiocManager,
	 Biostrings,
	 reshape2,
	 gridExtra
	 png, 
	 tiff,
	 jpeg,
	 pdftools,
	 grid
	 
Install these packages using R (versions >=3.6):
	
``` r
install.packages(c("BiocManager", "tidyverse", "UpSetR", "scales", "ape", 
                    "reshape2", "gridExtra","png", "tiff", "jpeg", "pdftools", "grid"))

BiocManager::install("Biostrings")
```	 

- **External software (REQUIRED):** [NCBI BLAST 2.9.0+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
   
	- Go to page: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/

	- Save `ncbi-blast-2.9.0+-x64-linux.tar.gz` in the local directory

	- To install, simply extract the downloaded package using the following command

	```
	tar zxvpf ncbi-blast-2.9.0+-x64-linux.tar.gz
	```

	- Configure of NCBI BLAST+ executables using the following command. This would append the path to the new BLAST bin directory to the existing PATH setting.


	```
	export PATH=$PATH:$HOME/ncbi-blast-2.9.0+/bin
	```
	
	- Refer [here](https://www.ncbi.nlm.nih.gov/books/NBK52640/) for source documentation

	  Note 1: R package assumes `blastn` and `makeblastdb` files are in path
    
 	  Note 2: BLAST version 2.9.0+ was used for the present program although other BLAST+ similar to version 2.9.0+ parameters might also run without problems
    
- **CP gene Database Download (REQUIRED)**

 The CP gene database can be downloaded either of the following 3 ways:

- Using **R** console
	
``` r 	
# Specify CP gene database URL 
url <- "https://raw.githubusercontent.com/ramadatta/CPgeneProfiler/master/testData/db/NCBI_BARRGD_CPG_DB.fasta"

# Specify destination where CP gene database file should be saved 
path <- "/home/user/db" # Can change to prefarable location
setwd(path)
destfile <- "NCBI_BARRGD_CPG_DB.fasta"

# Download the CP gene database file to the folder set in "path"
download.file(url, destfile)
```

- Using **UNIX/Linux** command line terminal and download the `db` folder with SVN

```
svn export https://github.com/ramadatta/CPgeneProfiler/trunk/testData/db
```
- Else, simply [Click](https://downgit.github.io/#/home?url=https://github.com/ramadatta/CPgeneProfiler/tree/master/testData/db) to save database folder and uncompress the `db.zip` folder.    
    
#### **Installation**

##### **From Github**

The R package is available through github repository can be installed using devtools.

``` r 
install.packages("devtools")
devtools::install_github("ramadatta/CPgeneProfiler")
```

##### **Check installation**
Ensure you have installed the package properly:

``` r
library("CPgeneProfiler")
?CPgeneProfiler
```

### **Example Usage**
 
`CPgeneProfiler` package can be run using the following functions: `cpblast()`, `filt_blast()`, `cocarriage()`, `cpprofile()`, `upsetR_plot()`, `plot_conlen()`,`assemblystat()`, `cp_summarize()`, `db_summary()`. Below are the examples of `CPgeneProfiler` commands usage. 


- Generates NCBI BLAST Results by aligning input genome assemblies against Carbapenamase (CP) gene database

``` r
cpblast(fastalocation = "/home/user/CPgeneProfiler/testData/fasta", dblocation = "/home/user/CPgeneProfiler/testData/db")

```
- Filtering NCBI BLAST Results based on CP gene coverage and percent identity. By default, alignment coverage of CP gene (`cpgcov`) and percentage identity (`cpgpident`) are set to 100%.

``` r
filt_blast(cpgcov = 100, cpgpident = 100)
```

- Report cocarriage of CP genes across all the input genome assemblies (Default:`cpgcov = 100` and `cpgpident = 100`)

``` r
cocarriage(cpgcov = 100, cpgpident = 100)
```

- Generate CP gene Profile across all the input genome assemblies (Default:`cpgcov = 100` and `cpgpident = 100`)

``` r
cpprofile(xlab = "Carbapenamase Genes", ylab = "Assembly", title = "Carbapenamase Gene Profile Heatmap")
```

- Plot CP gene contig length distributions across all the input genome assemblies (Default:`cpgcov = 100` and `cpgpident = 100`)

``` r
plot_conlen(outputType = "png", xlab = "Contig Length", ylab = "Number of Contigs", title = " Contig Length Distribution", colorfill = "#F99245")
```
- Generate set intersection plot of CP genes across all the input genome assemblies (Default:`cpgcov = 100` and `cpgpident = 100`)

``` r
upsetR_plot(outputType="png", width = 2000, height = 2000, res = 250, xlab="Carbapenamase Gene Set Size", ylab="Number of genome assemblies",cpgcov=100, cpgpident=100, order.by = "degree",nsets = 40, number.angles = 0,point.size = 1.5, line.size = 1,sets.bar.color = "red")
```
- Generate basic assembly statistics such N50, N90 and Assembly Size and plots comparing Assembly Size with N50, N90 stats

``` r
assemblystat(fastalocation = "/home/user/CPgeneProfiler/testData/fasta")
```

- Summarize the plots and organize all the output files in specific folders

``` r
cp_summarize(outdir_loc = "/home/user/Desktop", outdir = "CPgeneProfiler_Output")
```

- Database information

``` r
db_summary()
```

### **Working with `testData`**

#### **Test Data Download**

To test the package with test input data, Go to the UNIX/Linux command line terminal and download the `fasta` folder with SVN

```
svn export https://github.com/ramadatta/CPgeneProfiler/trunk/testData/fasta
```
else simply [Click](https://downgit.github.io/#/home?url=https://github.com/ramadatta/CPgeneProfiler/tree/master/testData/fasta) to save fasta folder and uncompress the `fasta.zip` folder.


##### **1) A simple NCBI BLAST using `cpblast()` command**

As a first step, CPgeneProfiler generates NCBI BLAST Results by aligning input genome assemblies against Carbapenamase (CP) gene database. Now that you already have a directory with fasta files (should have extensions `.fasta` or `.fa`) in `fasta` folder and cp gene database sequence in `db` folder, you can specify the path of both directories as an input and run the package with `cpblast()` command.

``` r
cpblast(fastalocation = "/home/user/CPgeneProfiler/testData/fasta", dblocation = "/home/user/CPgeneProfiler/testData/db", num_threads = 4, evalue = "1e-3")
```
The users can adjust BLAST parameters `num_threads`, `evalue`, `word_size` and `max_target_seqs` accordingly. If not adjusted and command is simply executed with the file locations for `fasta` and `db`, then default parameters are used for the analysis.


##### **2) Filtering BLAST results using `filt_blast()` command**

`filt_blast()` then filters the output BLAST results obtained from `cpblast()` command. This filtering is to find the presence of CP genes given a particular CP gene coverage and Percentage Identity. Therefore, the BLAST hits are filtered based on CP gene coverage and Percentage Identity. By default, CP Gene Coverage and Percentage Identity are set at a threshold of 100% (cpgcov=100, cpgpident=100). This means that a CP gene should have 100% alignment length and 100% identity, without even a single mismatch with the input genome sequence. The default parameters can be adjusted.

``` r
filt_blast(cpgcov = 100, cpgpident = 100)
```

This should generate the following table:

| assemblyName     | qseqid                                    | sseqid  | qlen   | slen | qstart | qend   | length | pident | cov |
|------------------|-------------------------------------------|---------|--------|------|--------|--------|--------|--------|-----|
| genome_001.fasta | 4_length=71861_depth=1.95x_circular=true  | KPC-2   | 71861  | 918  | 3810   | 4727   | 918  | 100 | 100 |
| genome_001.fasta | 5_length=71851_depth=1.95x_circular=true  | KPC-2   | 71851  | 918  | 3810   | 4727   | 918  | 100 | 100 |
| genome_002.fasta | 5_length=51479_depth=1.31x_circular=true  | OXA-181 | 51479  | 998  | 31280  | 32277  | 998  | 100 | 100 |
| genome_003.fasta | 2_length=316292_depth=2.71x_circular=true | NDM-1   | 316292 | 1013 | 149582 | 150594 | 1013 | 100 | 100 |
| genome_003.fasta | 2_length=316292_depth=2.71x_circular=true | OXA-181 | 316292 | 998  | 49123  | 50120  | 998  | 100 | 100 |
| genome_004.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 49850  | 50847  | 998  | 100 | 100 |
| genome_004.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 43441  | 44438  | 998  | 100 | 100 |
| genome_004.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 37032  | 38029  | 998  | 100 | 100 |
| genome_005.fasta | 2_length=79441_depth=2.21x_circular=true  | KPC-2   | 79441  | 918  | 11390  | 12307  | 918  | 100 | 100 |
| genome_005.fasta | 2_length=79441_depth=2.21x_circular=true  | KPC-2   | 79441  | 918  | 3810   | 4727   | 918  | 100 | 100 |
| genome_007.fasta | 2_length=246497_depth=2.06x_circular=true | KPC-2   | 246497 | 918  | 178390 | 179307 | 918  | 100 | 100 |
| genome_007.fasta | 2_length=246497_depth=2.06x_circular=true | IMP-26  | 246497 | 861  | 145325 | 146185 | 861  | 100 | 100 |
| genome_008.fasta | 3_length=41186_depth=4.61x_circular=true  | NDM-1   | 41186  | 1013 | 23027  | 24039  | 1013 | 100 | 100 |
| genome_009.fasta | 4_length=41182_depth=4.10x_circular=true  | NDM-1   | 41182  | 1013 | 23023  | 24035  | 1013 | 100 | 100 |
| genome_010.fasta | 5_length=51479_depth=1.31x_circular=true  | OXA-181 | 51479  | 998  | 31280  | 32277  | 998  | 100 | 100 |
| genome_010.fasta | 3_length=41186_depth=4.61x_circular=true  | NDM-1   | 41186  | 1013 | 23027  | 24039  | 1013 | 100 | 100 |
| genome_012.fasta | 2_length=79441_depth=2.21x_circular=true  | KPC-2   | 79441  | 918  | 11390  | 12307  | 918  | 100 | 100 |
| genome_012.fasta | 2_length=79441_depth=2.21x_circular=true  | KPC-2   | 79441  | 918  | 3810   | 4727   | 918  | 100 | 100 |
| genome_013.fasta | 3_length=41186_depth=4.61x_circular=true  | NDM-1   | 41186  | 1013 | 23027  | 24039  | 1013 | 100 | 100 |
| genome_014.fasta | 4_length=41182_depth=4.10x_circular=true  | NDM-1   | 41182  | 1013 | 23023  | 24035  | 1013 | 100 | 100 |
| genome_015.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 49850  | 50847  | 998  | 100 | 100 |
| genome_015.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 43441  | 44438  | 998  | 100 | 100 |
| genome_015.fasta | 3_length=66727_depth=0.76x                | OXA-181 | 66727  | 998  | 37032  | 38029  | 998  | 100 | 100 |

##### **3) Finding cocarriage genes using `cocarriage()` command**

`cocarriage()` command finds if two or more CP genes exists in same contig or multiple contigs across all the input genome assemblies. This function can be used only after running `filt_blast()`. By default, parameters such as CP Gene Coverage and Percentage Identity are set to 100% (cpgcov=100, cpgpident=100) but can be adjusted as per requirement.

``` r
cocarriage(cpgcov = 100, cpgpident = 100)
```

##### **4) Finding CP gene profile using `cpprofile()` command**

`cpprofile()` creates a heatmap of CP gene profile across the input genome assemblies. By default, the command generates `png` image but user can also create image with other output formats (jpeg/tiff/pdf) and parameters such as width, height of image, label, titles and colors of the heatmap can be adjusted as per requirement.

``` r
cpprofile(outputType="png", width = 2000, height = 2000, res = 250, xlab="Carbapenamase Genes", ylab="Assembly", title="Carbapenamase Gene Profile Heatmap", titlesize=15, labelsize=12, colorcode_low = "#143D59", colorcode_high = "#F4B41A", cpgcov=100, cpgpident=100)
```
<p align="center">
<img src="https://user-images.githubusercontent.com/3212461/90124487-10cca700-dd93-11ea-9572-dfc44a190dd6.png" width="45%"></img> 
</p>

##### **5) Plot CP gene contig length distribution using `plot_conlen()` command**

`plot_conlen()` generates length distribution for all the CP gene contigs present across all the input genome assemblies. By default, the command generates `png` image but user can also create image with other output formats (jpeg/tiff/pdf) and parameters such as width, height of image, label, titles and colors of the heatmap can be adjusted as per requirement.

``` r
plot_conlen(outputType="tiff", width = 700, height = 700, res = 150, xlab="Contig Length", ylab="Number of Contigs", title=" Contig Length Distribution",element_text_angle=90,unit="KB", breaks=15, colorfill = "#F99245",cpgcov=100, cpgpident=100)
```
<p align="center">
<img src="https://user-images.githubusercontent.com/3212461/90124507-17f3b500-dd93-11ea-9f37-6167ec79b8a1.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90124510-188c4b80-dd93-11ea-8609-1d3fbfc3ef43.png" width="45%"></img>
</p>

##### **6) Generate assembly statistics using `assemblystat()` command**

`assemblystat()` generates basic assembly stats which includes N50 size, N90 size and Genome assembly size. This function also generates Assembly Size vs N50 & Assembly Size vs N90 plots. This function requires the location of fasta file directory. By default, the command generates `png` image plots.

``` r
assemblystat("/home/user/CPgeneProfiler/testData/fasta", outputType="png", width = 700, height = 700, res = 150, geom_point_size=3, n50colorfill = "#0072B2", n90colorfill = "#D55E00")
```
<p align="center">
<img src="https://user-images.githubusercontent.com/3212461/90954281-a825ae80-e4a5-11ea-863d-47678a00aba1.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90954280-a6f48180-e4a5-11ea-93ad-908e0e339674.png" width="45%"></img>
</p>

##### **7) Generate Set Intersection of CP genes using `upsetR_plot()` command**

`upsetR_plot()` generates set intersection plot of CP genes across all the input genome assemblies. By default, the command generates `png` image but user can change the output image type, width and height of image, label, titles and colors.

``` r
upsetR_plot(outputType="png", width = 2000, height = 2000, res = 250, xlab="Carbapenamase Gene Set Size", ylab="Number of genome assemblies",cpgcov=100, cpgpident=100, order.by = "degree",nsets = 40, number.angles = 0,point.size = 1.5, line.size = 1,sets.bar.color = "red")
```
<p align="center">
<img src="https://user-images.githubusercontent.com/3212461/90124536-1f1ac300-dd93-11ea-9acc-9435b79a7f1c.png" width="45%"></img> 
</p>

##### **8) Summarize all the results using `cp_summarize()` command**

`cp_summarize()` arranges all the output files generated from above commands into respective folders. This also creates a summary of all the plots from CPgeneProfiler output into a single PDF file. Users can specify the output directory name and summary pdf name and also can provide the location of where the output folder to be generated. Note: All the output image plots need to be in the same format i.e, either png/tiff/jpeg.

``` r
cp_summarize(outdir = "CPgeneProfiler_Output", report="Summary" , image = "png")
```

##### **9) Find Database summary details using `db_summary()` command**

`db_summary()` command displays the details of Database. This includes Database Name, Version, Total number of sequences in currently in CP gene Database, Last updated date, Reference web link from which database was downloaded.

``` r
db_summary()
[1] "DATABASE: NCBI Bacterial Antimicrobial Resistance Reference Gene Database"
[1] "VERSION: 2020-07-16.2"
[1] "SEQUENCES: 875"
[1] "DBTYPE: nucl"
[1] "DATE: 2020-Aug-23"
[1] "Reference Gene Catalog: ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.8/"
```

##### **Output Files**

* A folder "CPgeneProfiler_Output" with the following files in respective directories

Command | File | Description
--------|------|--------------
cpblast() | blastResults.txt | Blast Results of contigs against the CP genes
filt_blast() | blastResults.filt.txt | Filtered blast results with contains contigs matching CP genes (default: 100% identity and 100% coverage)
cocarriage() | Cocarriage_Report.txt | Information of Number of assemblies with the co-carriage broken down to category
cocarriage() | Cocarriage_combinedResults.txt | Combined cocarriage details across all the assemblies
cocarriage() | DiffCP_DiffContig.txt | Information of assemblies with different CP genes present in different contigs
cocarriage() | DiffCP_SameContig.txt | Information of assemblies with different CP genes present in same contigs 
cocarriage() | SameCP_DiffContig.txt | Information of assemblies with same CP genes present in different contigs
cocarriage() | SameCP_SameContig.txt | Information of assemblies with same CP genes present in same contigs
cpprofile() | CPgeneProfile.png | Carbapenamase Gene Profile (default:png)
plot_conlen() | CPContigSizeDist.txt | Contig Size distribution
plot_conlen() | "CPGene"_Contig_Dist.png | CP gene contig length distribution (for each CP gene a separate distribution plot is generated. Default image format: png)
upsetR_plot() | cp_presence-absence_matrix.csv | Presence-absence matrix of CP genes across assemblies
upsetR_plot() | upset_plot.pdf | Set intersection plot of CP genes across all the input genome assemblies
assemblystat() | N50_N90.pdf | Assembly Size vs N50, N90 plots
assemblystat() | assemblyStats.txt | A simple text file with N50, N90, Assembly Size for each assembly
cp_summarize() | SummaryPlots.pdf | All the plots in a single pdf file

## A few example plots 

# `CP Gene Profile HeatMap`
<img src="https://user-images.githubusercontent.com/3212461/90596495-a946a980-e221-11ea-9904-88d5b5c66cc7.png" width="45%"></img>
<img src="https://user-images.githubusercontent.com/3212461/90596501-aba90380-e221-11ea-9cbb-062031eb8618.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90596505-ad72c700-e221-11ea-8340-db09e1061353.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90596508-aea3f400-e221-11ea-89a2-562aa408cbdf.png" width="45%"></img> 

# `Assembly Size vs N50 and N90 plots`

<img src="https://user-images.githubusercontent.com/3212461/90954281-a825ae80-e4a5-11ea-863d-47678a00aba1.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90954280-a6f48180-e4a5-11ea-93ad-908e0e339674.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90953985-dce43680-e4a2-11ea-91ba-750b72e731dd.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90953986-de156380-e4a2-11ea-8755-337b95ec4662.png" width="45%"></img> 

# `UpsetR plot` (orderby: freq vs degree)
<img src="https://user-images.githubusercontent.com/3212461/90954026-52500700-e4a3-11ea-9b1c-ad87b04c845b.png" width="45%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90124536-1f1ac300-dd93-11ea-9acc-9435b79a7f1c.png" width="45%"></img>

# `Contig Length Distribution Plots`
<img src="https://user-images.githubusercontent.com/3212461/90124507-17f3b500-dd93-11ea-9f37-6167ec79b8a1.png" width="30%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90124510-188c4b80-dd93-11ea-8609-1d3fbfc3ef43.png" width="30%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90954140-69dbbf80-e4a4-11ea-95ec-03c398269f3a.png" width="30%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90953981-d6ee5580-e4a2-11ea-97a4-fe30ddabbd77.png" width="30%"></img> 
<img src="https://user-images.githubusercontent.com/3212461/90953983-d81f8280-e4a2-11ea-81af-84668ff530a5.png" width="30%"></img>
<img src="https://user-images.githubusercontent.com/3212461/90953984-d8b81900-e4a2-11ea-9793-0e0026e7ba38.png" width="30%"></img> 

##### **Version**

version 2.1.1

- Modularised single function to multiple functions
- Added more flexibility to each function
- Changed the database from ARG-ANNOT to NCBI Bacterial AMR Reference Gene Database

version 2.1.0

- Changed the heatmap colors
- User can decide to download the ARG-annot database to desired location
- Both fasta folder and db location are to be provided for the package

# Community Guidelines

## How to Contribute

In general, you can contribute to this project by creating [issues](https://github.com/ramadatta/CPgeneProfiler/issues).
You are also welcome to contribute to the source code directly by forking the project, modifying the code, and creating [pull requests](https://github.com/ramadatta/CPgeneProfiler/pulls).
If you are not familiar with pull requests, check out [this post](https://guides.github.com/activities/forking/).
Please use clear and organized descriptions when creating issues and pull requests.

Please note that ``CPgeneProfiler`` is released with a [Contributor Code of Conduct](https://github.com/ramadatta/CPgeneProfiler/blob/master/Code_of_Conduct.md). By contributing to this project, you agree to abide by its terms.

## Bug Report and Support Request

You can use [issues](https://github.com/ramadatta/CPgeneProfiler/issues) to report bugs and seek support.
Before creating any new issues, please check for similar ones in the issue list first. 

## Citation

If you publish the results of CPgeneProfiler, please also cite the following software and the database:

Sridatta et al., (2020). CPgeneProfiler: A lightweight R package to profile the Carbapenamase genes from genome assemblies. Journal of Open Source Software, 5(54), 2473, https://doi.org/10.21105/joss.02473

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10(1), 421.doi:10.1186/1471-2105-10-421

Jake R Conway, Alexander Lex, Nils Gehlenborg, UpSetR: an R package for the visualization of intersecting sets and their properties, Bioinformatics, Volume 33, Issue 18, 15 September 2017, Pages 2938–2940, https://doi.org/10.1093/bioinformatics/btx364

NCBI Bacterial Antimicrobial Resistance Reference Gene Database (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047)



