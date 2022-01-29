# nf-rnaSeqMetagen
[![GitHub license](https://img.shields.io/github/license/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/blob/master/LICENSE) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) [![GitHub stars](https://img.shields.io/github/stars/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/stargazers) [![GitHub issues](https://img.shields.io/github/issues/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/issues)

[biotools:nf-rnaseqmetagen](https://bio.tools/nf-rnaseqmetagen)

*nf-rnaSeqMetagen* is a [Nextflow](http://nextflow.io/)

For detailed documentation: https://phelelani.github.io/nf-rnaSeqMetagen/


<!-- To use the `nf-rnaSeqMetagen` pipeline, the following dependencies are required: -->
<!--    1. Installed softwares: -->
<!--       - [`Nextflow`](https://www.nextflow.io/) -->
<!--       - [`Singularity`](http://singularity.lbl.gov/) -->
<!--    2. Reference genome, annotation and indexes -->
<!--       - Reference genome (`.fa`/`.fasta`) and genome annotation (`.gtf`) files. -->
<!--       - Reference genome indexes (`bowtie2` & `STAR` - see *1.3.* below on how to generate the indexes). -->
 
<!-- --- -->

<!-- <p align="center"> -->
<!--   <img width="600" src="nf-rnaSeqMetagen.png"> -->
<!-- </p> -->

<!-- ## 1. Obtaining the `nf-rnaSeqMetagen` pipeline and preparing data -->
<!-- First, you need to clone the `nf-rnaSeqMetagen` repository onto you machine. You can either use `git` or `nextflow` (see the two methods below). I recommend using `nextflow` and creating you own `config` file (will explain later) for executing the workflow in the directory of your choosing. The rest of this documentation assumes that you have used `nextflow` to clone this workflow - If your're an expert and have used `git` to clone the workflow - you know what to do :) -->
<!-- ```bash -->
<!-- ## Using nextflow -->
<!-- nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen -->
<!-- ``` -->
<!-- Content of the repository (will be in "$HOME/.nextflow/assets/phelelani/nf-rnaSeqCount"): -->
<!-- ```bash -->
<!-- nf-rnaSeqMetagen -->
<!-- |--containers                       ## Folder for Singularity images and recipes (in case you want to build yourself). All downloaded images go here! -->
<!-- |  |--Singularity.kraken2           ## Singularity recipe file for -->
<!-- |  |--Singularity.multiQC           ## Singularity recipe file for -->
<!-- |  |--Singularity.star              ## Singularity recipe file for -->
<!-- |  |--Singularity.trinity           ## Singularity recipe file for -->
<!-- |  |--Singularity.upset             ## Singularity recipe file for -->
<!-- |--templates                        ## Folder for extra scripts for the pipeline. -->
<!-- |  |--create_matrix.R               ## Script for -->
<!-- |  |--get_taxons.sh                 ## Script for -->
<!-- |--LICENSE                          ## Duh! -->
<!-- |--main.config                      ## User configuration file! All inputs, outputs and options GO HERE!! ONLY file that SHOULD be modified by user! -->
<!-- |--main.nf                          ## Main nf-rnaSeqMetagen nextflow scripts. -->
<!-- |--nextflow.config                  ## Pipeline configuration file! DO NOT EDIT!!! -->
<!-- |--nf-rnaSeqMetagen.png             ## Pipeline flow diagram -->
<!-- |--README.md                        ## Duh! -->
<!-- ``` -->
<!-- To get the `help menu` for the workflow, execute the following from anywherre on your system aftercloning the repository: -->
<!-- ``` -->
<!-- nextflow run nf-rnaSeqMetagen --help -->
<!-- ``` -->
<!-- The command above will give you the following usage information and options for running the `nf-rnaSeqMetagen` workflow: -->
<!-- ``` -->
<!-- ==================================================================================================== -->
<!-- #####################################  nf-rnaSeqMetagen v0.2   ##################################### -->
<!-- ==================================================================================================== -->

<!-- USAGE: -->
<!-- nextflow run nf-rnaSeqMetagen -profile "slurm" --data "/path/to/data" --genome "/path/to/genome.fa" --genes "/path/to/genes.gtf" -->

<!-- HELP: -->
<!-- nextflow run nf-rnaSeqMetagen --help -->

<!-- MANDATORY ARGUEMENTS: -->
<!-- -profile     STRING    Executor to be used. Available options: -->
<!-- 				"standard"          : Local execution (no job scheduler). -->
<!--                 "slurm"             : SLURM scheduler. -->
<!-- --mode       STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqMetagen). -->
<!--                        Available options: -->
<!-- 				"prep.Containers"   : For downloading Singularity containers used in this workflow. -->
<!--                 "prep.STARIndex"    : For indexing your reference genome using STAR. -->
<!--                 "prep.BowtieIndex"  : For indexing your reference genome using Bowtie2. -->
<!--                 "prep.KrakenDB"     : For building the Kraken2 database. -->
<!--                 "run.FilterClassify": For performing metagenomics analysis, i.e., filtering and classification. -->
<!-- --data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files: -->
<!-- 				[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ] -->
<!-- --genome     FILE      The whole genome FASTA sequence. Supported FASTA files: -->
<!--     			[ fasta | fa | fna ] -->
<!-- --genes      FILE      The genome annotation GFT file. Supported GTF file: -->
<!-- 				[ gtf ] -->
<!-- --db         FOLDER    Path to where the Kraken2 database will be saved (or where it is located if already created). -->
<!--                        Default: $PWD/kraken2db -->

<!-- OPTIONAL ARGUEMENTS: -->
<!-- --help                 To show this menu. -->
<!-- --out        FOLDER    Path to where the output should be directed. -->
<!--                        Default: $PWD/results_nf-rnaSeqMetagen -->
<!-- --pairedEnd            If working with paired-end FASTQ files (default). -->
<!-- --singleEnd            If working with single-end FASTQ files. -->
<!-- --max_memory STRING    Maximum memory you have access to. -->
<!--                        Default: "200.GB" -->
<!-- --max_cpus   STRING    Maximum CPUs you have access to. -->
<!--                        Default: "24" -->
<!-- --max_time   STRING    Maximum time you have access to. -->
<!--                        Default: "24.h" -->
<!-- ==================================================================================================== -->
<!-- ``` -->

<!-- --- -->

<!-- ### 1.1. Download test datasets (optional) -->
<!-- We will now download the reference genome (along with its annotation file) from Ensembl. We will also download the FASTQ files from the H3ABioNet site, which we will analyse using the `nf-rnaSeqMetagen` workflow. *__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqMetagen` workflow!* -->

<!-- - [x] Download and decompress the mouse reference genome along with its annotation: -->
<!-- ``` -->
<!-- ## Make a directory for the reference genome: -->
<!-- mkdir reference -->

<!-- ## Download the reference genome (FASTA) and annotation file (GTF) files and put them into the newlly created directory: -->
<!-- wget -c -O reference/genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz -->
<!-- wget -c -O reference/genes.gtf.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz -->
<!-- gunzip reference/genome.fa.gz -->
<!-- gunzip reference/genes.gtf.gz -->
<!-- ``` -->

<!-- - [x] Download RNA-seq test dataset from H3ABioNet: -->
<!-- ``` -->
<!-- ## Make a directory for the data: -->
<!-- mkdir data -->

<!-- ## Download the data: -->
<!-- for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c -O data/$sample http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done -->
<!-- ``` -->
<!-- ### 1.2. Download the `Singularity` containers (required to execute the pipeline): -->
<!-- ```bash -->
<!-- nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.Containers -->
<!-- ``` -->

<!-- ### 1.3. Generating genome indexes. -->
<!-- To generate the `STAR` and `Bowtie2` genome indexes, run the following commands: -->
<!-- ```bash -->
<!-- ## Generate STAR indexes -->
<!-- nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.STARIndex --genome "$PWD/reference/genome.fa" --genes "$PWD/reference/genes.gtf" -->

<!-- ## Generate Bowtie2 indexes: -->
<!-- nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.BowtieIndex --genome "$PWD/reference/genome.fa" --genes "$PWD/reference/genes.gtf" -->
<!-- ``` -->

<!-- ### 1.4. Creating the Kraken2 database: -->
<!-- To create the Kraken2 database, run the following command: -->
<!-- ```bash -->
<!-- ## Create Kraken2 database -->
<!-- nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.KrakenDB --db $PWD/K2DB -->
<!-- ``` -->

<!-- We are now ready to execute the workflow! -->

<!-- --- -->

<!-- ## 2. Executing the main `nf-rnaSeqMetagen` pipeline -->
<!-- As seen on the `help menu` above, there are a couple of options that you can use with this workflow. It can become a bit tedious and confusing having to specify these commands everytime you have to execute the each section for the analysis. To make your life easier, we will create a configuration script that we will use in this tutorial (we will pass this using the `-c` option of `nextflow`). You can name it whatever you want, but for now, lets call it `myparams.config`. We will add the mandatory arguements for now, but as you become more farmiliar with the workflow - you can experiment with other options. You can use your favourite text editor to create the `myparams.config` file. Copy and paste the the parameters below: -->
<!-- ``` -->
<!-- params { -->
<!--     data    = "$PWD/data" -->
<!--     db      = "$PWD/K2DB" -->
<!--     genome  = "$PWD/reference/genome.fa" -->
<!--     genes   = "$PWD/reference/genes.fa" -->
<!-- } -->
<!-- ``` -->
<!-- Obviously - the above `myparams.config` assumes that you have been following this tutorial. If you have your data lying around somewhere in your system, you need to put the full path to where your the `data`, `genome` and `genes` files are. Since the `--mode` will keep changing, we will add this on the command as we do the analysis. Now that we have the mandatory arguements in our `myparams.config`, lets do some analysis -->

<!-- ### 2.1. Read Filtering and Classification: -->
<!-- To perform filtering of host reads and classification of exogeneous reads, use this command: -->
<!-- ```bash -->
<!-- nextflow run nf-rnaSeqMetagen -profile slurm --mode run.FilterClassify -c myparams.config -->
<!-- ``` -->

<!-- --- -->

<!-- ## 3. Explore `nf-rnaSeqMetagen` results -->

<!-- ``` -->
<!-- - [1] Sample analysis directories  =>    `<output_directory>/<sample_1> .. <sample_N>` -->
<!-- - [2] MultiQC                      =>    `<output_directory>/MultiQC` -->
<!-- - [3] Upset tool                   =>    `<output_directory>/upset` -->
<!-- - [4] Workflow tracing             =>    `<output_directory>/workflow-tracing -->
<!-- ``` -->
<!-- In addition to the directories created in the results directory, a directory `workflow-tracing` is created to monitor the resources used for filtering and classification. This directory will contain 4 files: -->
<!-- - `nf-rnaSeqMetagen_report.html` -->
<!-- - `nf-rnaSeqMetagen_timeline.html` -->
<!-- - `nf-rnaSeqMetagen_trace.txt` -->
<!-- - `nf-rnaSeqMetagen_flow.dot` -->

<!-- These files contain detailed information on the resources (CPU, MEMORY and TIME) usage of each of the process in the pipeline. The `<output_directory>` directory structure is summarized below: -->

<!-- ```bash -->
<!-- <output_directory> -->
<!-- |--<sample_1> -->
<!-- |  |--<sample_1>.fasta.html -->
<!-- |  |--<sample_1>.reads.html -->
<!-- |  |--<sample_1>_classified.fasta -->
<!-- |  |--<sample_1>_fasta.krak -->
<!-- |  |--<sample_1>_fasta.kron -->
<!-- |  |--<sample_1>_reads.kron -->
<!-- |  |--taxon_sequences -->
<!-- |  |  |--taxid<1>.fasta .. taxid<n>.fasta -->
<!-- |  |--trinity_<sample_1> -->
<!-- |  |  |--Trinity.fasta -->
<!-- .. -->
<!-- |--<sample_N> -->
<!-- |  |--<sample_N>.fasta.html -->
<!-- |  |--<sample_N>.reads.html -->
<!-- |  |--<sample_N>_classified.fasta -->
<!-- |  |--<sample_N>_fasta.krak -->
<!-- |  |--<sample_N>_fasta.kron -->
<!-- |  |--<sample_N>_reads.kron -->
<!-- |  |--taxon_sequences -->
<!-- |  |  |--taxid<1>.fasta .. taxid<n>.fasta -->
<!-- |  |--trinity_<sample_N> -->
<!-- |  |  |--Trinity.fasta -->
<!-- |--MultiQC -->
<!-- |  |--multiqc_data -->
<!-- |  |--multiqc_report.html -->
<!-- |--upset -->
<!-- |  |--data -->
<!-- |  |  |--nf-rnaSeqMetagen.csv -->
<!-- |  |  |--nf-rnaSeqMetagen.json -->
<!-- |--workflow-tracing -->
<!-- |  |--nf-rnaSeqMetagen_{report.html,timeline.html,trace.txt,flow.dot} -->
<!-- ``` -->
<!-- **NB:** I am working on further improving the pipleine and the associated documentation, feel free to share comments and suggestions! -->

<!-- --- -->
[![GitHub license](https://img.shields.io/github/license/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/blob/master/LICENSE) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu) [![GitHub stars](https://img.shields.io/github/stars/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/stargazers) [![GitHub issues](https://img.shields.io/github/issues/phelelani/nf-rnaSeqMetagen)](https://github.com/phelelani/nf-rnaSeqMetagen/issues)

[biotools:nf-rnaseqmetagen](https://bio.tools/nf-rnaseqmetagen)

`nf-rnaSeqMetagen` is a [Nextflow](http://nextflow.io/)
        
To use the `nf-rnaSeqMetagen` pipeline, the following are required:
1. Software dependencies:
   - [`Nextflow`](https://www.nextflow.io/)
   - [`Singularity`](http://singularity.lbl.gov/)
2. RNA-seq data (paired-end for now - support for single-ended reads to follow)
3. Reference genome (FASTA sequences) and its annotation file (GFT)

---

<p align="center">
  <img width="832px" src="assets/images/nf-rnaSeqMetagen.png">
</p>

---

## 1. Obtaining the `nf-rnaSeqMetagen` Pipeline and Preparing Data
First, you need to clone the `nf-rnaSeqMetagen` repository onto you machine. You can either use `git` or `nextflow`. I recommend you use `nextflow`. The rest of this documentation assumes that you have used `nextflow` to clone this workflow.
```bash
nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen
```
<script id="asciicast-308777" src="https://asciinema.org/a/308777.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="6" data-speed="3" data-loop="0"></script>

Content of the repository (located in `$HOME/.nextflow/assets/phelelani/nf-rnaSeqCount`):
<script id="asciicast-308808" src="https://asciinema.org/a/308808.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="12" data-speed="3" data-loop="0"></script>

To get the `help menu` for the workflow, execute the following command from anywherre on your system:
```
nextflow run nf-rnaSeqMetagen --help
```
<script id="asciicast-308804" src="https://asciinema.org/a/308804.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

---

### 1.1. Download test datasets (optional)
We will now download the reference genome (along with its annotation file) from Ensembl. We will also download the FASTQ files from the H3ABioNet site, which we will analyse using the `nf-rnaSeqMetagen` workflow. *__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqMetagen` workflow!* 

Make directories:
```
mkdir example
cd example
mkdir reference
mkdir data
```
<script id="asciicast-308945" src="https://asciinema.org/a/308945.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="8" data-speed="3" data-loop="0"></script>

Download and decompress the mouse reference genome along with its annotation:
```
wget -c -O reference/genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz
```
<script id="asciicast-308949" src="https://asciinema.org/a/308949.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="18" data-speed="3" data-loop="0"></script>


```
wget -c -O reference/genes.gtf.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
```
<script id="asciicast-308953" src="https://asciinema.org/a/308953.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="18" data-speed="3" data-loop="0"></script>

```
gunzip reference/genome.fa.gz
gunzip reference/genes.gtf.gz
```
<script id="asciicast-308955" src="https://asciinema.org/a/308955.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="8" data-speed="3" data-loop="0"></script>

Download RNA-seq test dataset from H3ABioNet: <a href="examples/data/get_data.sh" target="_blank">script</a>.
```
cd data
wget https://phelelani.github.io/nf-rnaSeqMetagen/examples/data/get_data.sh
```
<script id="asciicast-309156" src="https://asciinema.org/a/309156.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="14" data-speed="3" data-loop="0"></script>

```
sh get_data.sh
ls -l 
cd ..
```
<script id="asciicast-309177" src="https://asciinema.org/a/309177.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="12" data-speed="7" data-loop="0"></script>


### 1.2. Download the `singularity` containers (required to execute the pipeline):
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.Containers
```
<script id="asciicast-308816" src="https://asciinema.org/a/308816.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

### 1.3. Generating genome indexes.
To generate the `STAR` genome indexes, run the following commands:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.GenomeIndexes --genome "$PWD/reference/genome.fa" --genes "$PWD/reference/genes.gtf"
```
<script id="asciicast-309150" src="https://asciinema.org/a/309150.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>


### 1.4. Creating the Kraken2 database:
To create the Kraken2 database, run the following command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.KrakenDB --db $PWD/K2DB
```
<script id="asciicast-309125" src="https://asciinema.org/a/309125.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

We are now ready to execute the workflow!

---

## 2. Executing the Main `nf-rnaSeqMetagen` Pipeline
As seen on the `help menu` above, there are a couple of options that you can use with this workflow. It can become a bit tedious and confusing having to specify these commands everytime you have to execute the each section for the analysis. To make your life easier, we will create a configuration script that we will use in this tutorial (we will pass this using the `-c` option of `nextflow`). You can name it whatever you want, but for now, lets call it `myparams.config`. We will add the mandatory arguements for now, but as you become more farmiliar with the workflow - you can experiment with other options. You can use your favourite text editor to create the `myparams.config` file. Copy and paste the the parameters below:
```
params {
    data   = $PWD/data
    out    = $PWD/myresults
    genome = $PWD/reference/genome.fa
    genes  = $PWD/reference/gene.gtf
    db     = $PWD/K2DB
}

```

To perform filtering of host reads and classification of exogeneous reads, use this command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode run.FilterClassify -c myparams.config
```

---

## 3. Exploring `nf-rnaSeqMetagen` Results

```
- [1] Sample analysis directories  =>    `<output_directory>/<sample_1> .. <sample_N>`
- [2] MultiQC                      =>    `<output_directory>/MultiQC`
- [3] Upset tool                   =>    `<output_directory>/upset`
- [4] Workflow tracing             =>    `<output_directory>/workflow-tracing
```
### 3.1. MultiQC
View the full MultiQC report <a href="examples/output/MultiQC/multiqc_report.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/MultiQC/multiqc_report.html"></iframe>

### 3.2. Sample analysis directories
#### 3.2.1. Krona report: raw reads (SRR5074528)
View full Krona chart for raw reads <a href="examples/output/SRR5074528/SRR5074528_reads.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_reads.html"></iframe>

#### 3.2.2 Krona report: assembled reads (SRR5074528)
View full Krona chart for assembled reads <a href="examples/output/SRR5074528/SRR5074528_fasta.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_fasta.html"></iframe>

### 3.3. UpSet visualisation tool
View full UpSet plot <a href="examples/output/upset/index.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/upset/index.html"></iframe>

### 3.4. Workflow tracing
#### 3.4.1. Report
View full Nextflow report <a href="examples/output/workflow-tracing/nf-rnaSeqMetagen_report.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_report.html"></iframe>

#### 3.4.2. Timeline
View full timeline report <a href="examples/output/workflow-tracing/nf-rnaSeqMetagen_timeline.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_timeline.html"></iframe>

---
# UpSet

## About

UpSet is an interactive, web based visualization technique designed to analyze set-based data. UpSet visualizes both, set intersections and their properties, and the items (elements) in the dataset. Please see the project website at [http://vcg.github.io/upset/about](http://vcg.github.io/upset/about) for details about the technique, publications and videos.

## Demo

We are hosting a demo instance of UpSet at [http://vcg.github.io/upset](http://vcg.github.io/upset).

## R Package

An R package to generate UpSet plots is under development and available at [https://github.com/hms-dbmi/UpSetR](https://github.com/hms-dbmi/UpSetR).

## Local Deployment

1. Clone the repository using ```git clone``` or download and extract the [ZIP file](https://github.com/VCG/upset/archive/master.zip).
2. Launch the [Python SimpleHTTPServer](https://docs.python.org/2/library/simplehttpserver.html) in the project directory.
 
   ```
   $ python -m SimpleHTTPServer 8000
   ```

3. View UpSet in your browser at [localhost:8000](http://localhost:8000).

Alternatively you can also **run UpSet without a web server**. Chrome does not allow this by default, but Firefox works well. Simply open the index.html file in Firefox. 

## Configuring Datasets

See the project wiki for an [overview of the data definition file format](https://github.com/VCG/upset/wiki/Data-Import) used to describe tabular text files.


