
# ScanITD
[![Build Status](https://travis-ci.org/ylab-hi/ScanITD.svg?branch=master&status=passed)](https://travis-ci.org/ylab-hi/ScanITD)

Introduction
------------
ScanITD: detecting internal tandem duplication with robust variant allele frequency estimation

Getting Started
----------------
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

Prerequisites
----------------
You need Python 3.4 and later to run ScanITD.

### install necessary python packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) (python 3.7) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda pysam
conda install -c conda-forge scikit-bio
conda install -c anaconda numpy
conda install -c bioconda samtools  ## samtools/1.0 or newer is required
 ```
Usage
-------------------------
```
ScanITD.py -i input_bam_file -r indexed_refenence_genome_fasta -o output_vcf_filename_prefix [opts]
```

### Options:
```
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        BWA-MEM BAM file
  -r REF, --ref REF     reference genome in FASTA format
  -o OUTPUT, --output OUTPUT
                        output prefix
  -m MAPQ, --mapq MAPQ  minimal MAPQ in BAM for calling ITD (default: 15)
  -c AO, --ao AO        minimal observation count for ITD (default: 4)
  -d DP, --depth DP     minimal depth to call ITD (default: 10)
  -f VAF, --vaf VAF     minimal variant allele frequency (default: 0.1)
  -l ITD_LEN, --len ITD_LEN
                        minimal ITD length to report (default: 10)
  -n MISMATCH           maximum allowed mismatch bases of pairwise local
                        alignment (default: 3)
  -t TARGET, --target TARGET
                        Limit analysis to targets listed in the BED-format
                        file or a samtools region string
  -k, --keep            Kepp the ITD build BAM file
  -v, --version         show program's version number and exit
```  
#### Input:
```	
input_bam_file                    :input WES BAM file. (e.g., wes-seq.bam)
indexed_reference_genome_fasta    :specify reference genome in FASTA format (the reference genome should be indexed)
```
#### Output:
```	
output_vcf_filename_prefix        :specify the prefix of the output vcf file 
```
The name of the output VCF file will be __prefix__.itd.vcf.

The [__Test__](https://github.com/ylab-hi/ScanITD/tree/master/Test "Test directory") directory contains files and scripts for testing the installation.

License
----------------
This project is licensed under [MIT](https://opensource.org/licenses/MIT).

Contact
-----------------
Bug reports or feature requests can be submitted on the [ScanITD Github page](https://github.com/ylab-hi/ScanITD/issues).

Citation
------------------
Wang TY. and Yang R. [ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation](https://doi.org/10.1093/gigascience/giaa089 "ScanITD: Detecting internal tandem duplication with robust variant allele frequency estimation").
Sample data for testing
---
This directory contains one sample dataset for testing purpose to make sure you have installed ScanITD and its dependent packages sucessfully. __test.bam__ is a sliced WXS dataset (in BAM format) contains one FLT3-ITD. 

```
ScanITD.py -i test.bam -r /path/to/hg38.fa -o test 
```
The testing script will generate a VCF file named __test.itd.vcf__

The VCF file contains one ITD (POS=28034081; AB=0.12; SVLEN=102)
Step 1
---
Generate duplication positions using RSVSim
```
Rscript sim.R
```
tandemDuplications.csv will generated in this step.

Step 2
---
Prepare duplication file for svsim
```
prepare_svsim.py
```
tandemDuplications.txt will generated in this step.

Step 3
---
Generate rearranged genome in FASTA format
```
python create_indel_genome.py chr20.hg19.fa tandemDuplications.txt chr20.hg19.TDUP.fa
```
__chr20.hg19.TDUP.fa__ is the rearranged genome generated;
__chr20.hg19.fa__ is the unrranged genome.

Step 4
---
Generate simulated reads in various settings. (reads length, reads depth, variant allele frequency)
```
bash fastq_vaf10.sh # VAF=10%
bash fastq_vaf20.sh # VAF=20%
bash fastq_vaf50.sh # VAF=50%
```

Step 5
---
Map reads using BWA-MEM
```
bam_generator.py
```
Pay attention, you need to replace __/path/to/hg19.fa__ in bam_generator.py using a working hg19.fa
Instruction
----------------
If you start with raw reads (FASTQ) in application, you need do quality control and quality/adapter trimming for raw reads. 
Then you can use BWA-MEM or other soft-clipping aware NGS aligners to generate a BAM file feeding into ScanITD.

Prerequisites
----------------
You need Fastqc and Trim Galore to run QC.sh

### install necessary packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) (python 3.7) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda trim-galore
conda install -c bioconda fastqc
 ```
 
 Usage
-------------------------
```
trim_galore -o output_folder --paired --fastqc R1.fastq R2.fastq
```
#### Input:
```	
R1.fastq, R2.fastq   :input WES paired-end reads in FASTQ format file.
output_folder        :specify output folder name
```
