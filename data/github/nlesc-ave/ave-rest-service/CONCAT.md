# Data Processing

## getting vcfs and going from multiple (per accession) to single vcf file
generate compressed bcf file from vcf file:
```
bcftools view -O b 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz > 1001genomes_snp-short-indel_with_tair10_only_ACGTN.bcf

### Tomato

Download reference genome from ftp://ftp.solgenomics.net/tomato_genome/assembly/build_2.40/S_lycopersicum_chromosomes.2.40.fa.gz
Download all analysis (vcf) files from http://www.ebi.ac.uk/ena/data/view/PRJEB5235

gunzip: ./ERZ020502/RF_060_SZAXPI009336-14.vcf.gz: unexpected end of file
gunzip: ./ERZ020502/RF_060_SZAXPI009336-14.vcf.gz: uncompress failed
gunzip: ./ERZ020503/RF_062_SZAXPI009337-15.vcf.gz: unexpected end of file
gunzip: ./ERZ020503/RF_062_SZAXPI009337-15.vcf.gz: uncompress failed
```

~~concatenating vcf files:
`vcf-concat myvcfs/*.vcf.gz | gzip > out.vcf.gz`~~

files come from sra here:
ftp.sra.ebi.ac.uk
and the directory _vol1_
directories: ERZ020447 - ERZ020530


* ungzip vcf files
`gunzip sample.vcf.gz`
* each file compress with bgzip
`bgzip sample.vcf`
* index with tabix
`tabix -p vcf sample.vcf.gzip`
* merge into one variant file
`vcf-merge *.vcf.gz > all-snps.vcf`


## random access to fasta sequence
use pyfaidx
[GitHub - mdshw5/pyfaidx: Efficient pythonic random access to fasta subsequences](https://github.com/mdshw5/pyfaidx)
```py
from pyfaidx import Fasta
genome = Fasta('TAIR10.ga')

```

## random access to gff (gene annotatins)
* first need to sort the gff file by chrom/start
```sh
# doing it with bedtools
bedtools sort -i TAIR10_GFF3_genes.gff > TAIR10_GFF3_genes.sorted.gff
```

* then index with tabix
```sh
bgzip TAIR10_GFF3_genes.sorted.gff 
tabix -p gff TAIR10_GFF3_genes.sorted.gff.gz
```

* use `pysam.Tabixfile` to read the file, then iterate over it
```py
gff.fetch("Chr2", 1, 5000, parser=pysam.asGTF())
```

## random access to vcf
* fetch vcf or vcf.gz file
* sort by chromosome
```sh
 vcf-sort -t /mnt/disks/variant-store/tmp/ -c 1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz > variants_sorted.vcf
```
* block gzip
```sh
bgzip -c variants_sorted.vcf variants_sorted.vcf.gz
```
* index with tabix
```sh
tabix -p vcf variants_sorted.vcf.gz
```
* convert vcf to bcf
```sh
bcftools view -O b variants_sorted.vcf.gz > variants_sorted.bcf
```

## queries
1. region of 50 kB 
	* reference sequence
	* all variants
	* all annotations
2. find the gene of interest
there is no good way to do it based on gff file
one solution is to put genes names and their location in the sqlite db (will it be fast enough?)

# ave-rest-service

[![Build Status](https://travis-ci.org/nlesc-ave/ave-rest-service.svg?branch=master)](https://travis-ci.org/nlesc-ave/ave-rest-service)
[![SonarCloud Gate](https://sonarcloud.io/api/badges/gate?key=ave-rest-service)](https://sonarcloud.io/dashboard?id=ave-rest-service)
[![SonarCloud Coverage](https://sonarcloud.io/api/badges/measure?key=ave-rest-service&metric=coverage)](https://sonarcloud.io/component_measures/domain/Coverage?id=ave-rest-service)
[![Docker Automated build](https://img.shields.io/docker/automated/ave2/allelic-variation-explorer.svg)](https://hub.docker.com/r/ave2/allelic-variation-explorer/)
[![DOI](https://zenodo.org/badge/44877908.svg)](https://zenodo.org/badge/latestdoi/44877908)

The Allelic Variation Explorer (AVE) is a web application to visualize (clustered) single-nucleotide variants across genomes.
The Allelic Variation Explorer rest service clusters genomic variants and lists the available datasets.
Combined with [ave-app](https://github.com/nlesc-ave/ave-app) will visualize clustered genomic variants for a certain genomic range in a genome browser.

![Screenshot of Allelic Variation Explorer](https://github.com/nlesc-ave/ave-rest-service/raw/master/docs/screenshot.png)

This service is the back end for the [ave-app](https://github.com/nlesc-ave/ave-app) front end.
The front end runs in the users web browser and communicates with the back end running on a web server somewhere.
The front end is the user interface and the back end is the service serving the variant, annotation and genomic data.
The front end and back end communicate with each other according to the [Swagger specification](https://swagger.io/) in [swagger.yml](avedata/swagger.yml).

* [Architecture](#architecture)
* [Deployment](#deployment)
  * [Demo](#demo)
  * [Secure connection](#secure-connection)
  * [Shutting down](#shutting-down)
  * [Update image](#update-image)
* [Data pre processing](#data-pre-processing)
  * [Genome sequence](#genome-sequence)
  * [Variants](#variants)
  * [Genes](#genes)
  * [Genomic features annotations](#genomic-features-annotations)
* [Data registration](#data-registration)
  * [Register genome](#register-genome)
  * [Register variants](#register-variants)
  * [Register genes](#register-genes)
  * [Register feature annotations](#register-feature-annotations)
  * [Deregister](#deregister)
* [Develop](#develop)
  * [Setup](#setup)
  * [Configuration](#configuration)
  * [Run service](#run-service)
  * [Run commands](#run-commands)
  * [Build Docker image](#build-docker-image)

## Architecture

![Architecture](https://github.com/nlesc-ave/ave-rest-service/raw/master/docs/architecture.png)

The ave-rest-service and ave-app are wrapped up in a [Docker image](https://hub.docker.com/r/ave2/allelic-variation-explorer/).
The Docker image is used for deploying the Allelic Variation Explorer on a server.

The Allelic Variation Explorer consists of the following parts working together:
* a running ave rest service
* an extracted [ave-app](https://bintray.com/nlesc-ave/ave/ave-app/latest#files) build archive.
* 2bit (genome sequence), bcf (variants) and bigbed (genes and feature annotations) data files, green in diagram
* a directory with full text indices for genes and features in [Whoosh](https://whoosh.readthedocs.io) format, filled by [data registration commands](#data-registration), red in diagram
* an meta database file, contains list of available datasets inside the application, filled by [data registration commands](#data-registration), yellow in diagram
* a [NGINX web server](http://nginx.org/), for hosting app and data files and proxy-ing ave rest service behind a single port
* a Docker image combining all above, see `./Dockerfile` for the instructions used to install all the parts

## Deployment

A Docker image is available on [Docker Hub](https://hub.docker.com/r/ave2/allelic-variation-explorer/).

Any change to the master branch of this repo or the [ave-app](https://github.com/nlesc-ave/ave-app) will trigger an automatic build of the [Docker image](https://hub.docker.com/r/ave2/allelic-variation-explorer/).

The Docker image contains no data and when data is added then the data will be lost when the Docker container is stopped/started.
To get a deployment which is persists it's data we will use directories on the server and mount these as volumes in the Docker container.
It expects the following volumes:

* /data, location for 2bit, bcf and bigbed data files. Hosted as http://&lt;aveserver&gt;/data
* /whoosh, full text indices for genes and features
* /meta, directory in which ave meta database is stored

Run the service with
```bash
# Use sub directories in the current working directory to persist data
mkdir data
mkdir whoosh
mkdir meta
docker run -d \
  -v $PWD/data:/data -v $PWD/whoosh:/whoosh -v $PWD/meta:/meta \
  -p 80:80 \
  --name ave ave2/allelic-variation-explorer
```

Command above will run web server on port 80 of host machine.

After deployment the server is running, but contains no data, see [Data pre processing](#data-pre-processing) chapter how to prepare data followed by the [Data registration](#data-registration) chapter how to add data.

### Demo

A demo Docker image with a sample dataset is available at https://hub.docker.com/r/ave2/ave-demo/ .

### Secure connection

The Docker container uses http which is an unencryted connection.
To use a secure connection (https), a reverse proxy with Let's Encrypt certificate can be put in front of the Docker container.

The Docker container must be run a port which is not 443.

Configure a web server like NGINX to server https on port 443 and proxy all requests to the container.
Use [Certbot](https://certbot.eff.org/) to generate the certificate pair and configure the web server.

See example server conf in commented out block in `./nginx.conf` file.

### Shutting down

The Docker container can be stopped using
```bash
docker rm -f ave
```

### Update image

Make sure the Docker container is not running.

The Docker image can be updated using

```bash
docker pull ave2/allelic-variation-explorer
```

## Data pre processing

Before data can be registered it has to be converted in the right format. Below describes pre processing steps to convert common formats into the formats the application expects.

The tools used are available inside the [Docker container](#deployment) or can be installed in an [Anaconda environment](#setup).
To perform the pre processing inside the Docker container, copy the raw files to the `/data` Docker volume and login to the Docker container with `docker exec -ti ave bash`.

### Genome sequence

Genome sequence in [2bit](https://genome.ucsc.edu/goldenpath/help/twoBit.html) format is used.

When you have a genome sequence in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, where each chromosome is a sequence in the file.

The FASTA file can be converted to 2bit using:

```sh
faToTwoBit genome.fa genome.2bit
```

The sequence identifiers (chromosome/contig) should match the ones in corresponding variants (bcf) genes (bigbed) and genomic feature annotations (bigbed) files.

### Variants

Variants need to be provided in single [BCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) formatted file.

VCF files can be converted to a BCF file in the following way:

1. Merge VCF files with [SAMtools](http://www.htslib.org/doc/samtools.html) and [VCFtools](http://vcftools.sourceforge.net/perl_module.html)

```sh
# vcf-merge requires bgzipped and tabix indexed VCF files so do that first
for f in $(ls *.vcf)
do
bgzip $f
tabix -p vcf $f.gz
done
vcf-merge *.vcf.gz > variants.vcf
```

2. Sort by chromosome with [VCFtools](http://vcftools.sourceforge.net/perl_module.html)
```sh
vcf-sort -c variants.vcf > variants.sorted.vcf
```

3. Compress by [bgzip](http://www.htslib.org/doc/tabix.html) and index with [tabix](http://www.htslib.org/doc/tabix.html)
```sh
bgzip variants.sorted.vcf
tabix -p vcf variants.sorted.vcf.gz
```

4. Convert to [BCF](https://samtools.github.io/hts-specs/BCFv2_qref.pdf) with [bcftools](https://samtools.github.io/bcftools/bcftools.html)
```sh
bcftools view -O b variants.sorted.vcf.gz > variants.sorted.bcf
```

5. Index with [bcftools](http://www.htslib.org/doc/bcftools.html)
```sh
bcftools index variants.sorted.bcf
```

### Genes

The genes (or transcripts) are rendered in a gene track.

Genes must be provided as [bigBed](http://genome.ucsc.edu/goldenPath/help/bigBed.html) formatted.
A bigBed file can be converted from a [BED formatted](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file

The gene bed file is expected to have the following columns:
1. chrom, name of chromosome
2. chromStart, start of gene, zero-indexed
3. chromStart, end of gene, zero-indexed
4. name, transcript identifier
5. score
6. strand
7. thickStart, location of start codon
8. thinkEnd, location of stop codon
9. itemRgb
10. blockCount, number of exons
11. blockSizes
12. blockStarts
13. gene identifier
14. description of gene

So it should have exons and start/stop codons for one gene on a single line.

Some gene bed files are available at http://bioviz.org/quickload/

To convert a gene bed file to bigbed use:
```bash
# Fetch chrom sizes
twoBitInfo genome.2bit chrom.sizes
# the version for mac os is available also available http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/bedToBigBed
chmod +x bedToBigBed
# description field is too long for bedToBigBed so it must be trimmed
gunzip -c S_lycopersicum_May_2012.bed.gz | perl -n -e 'chomp;@F=split(/\t/);$F[13] = substr($F[13],0,255); print join("\t", @F),"\n";'  > S_lycopersicum_May_2012.bed.trimmed
bedToBigBed -tab -type=bed12+2 S_lycopersicum_May_2012.bed.trimmed chrom.sizes S_lycopersicum_May_2012.bb
```

A BigBed file with more than 256 chromosomes will not render, see [issue 37](https://github.com/nlesc-ave/ave-app/issues/37).

### Genomic features annotations

Each feature file will render a feature track. The feature track name is the same as the file name.

Feature annotations must be provided as [bigBed](http://genome.ucsc.edu/goldenPath/help/bigBed.html) formatted files.

To convert a gff feature file to bigbed use:

```bash
# Download a gff file
wget ftp://ftp.solgenomics.net/tomato_genome/microarrays_mapping/A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.gff
# Sort gff
bedtools sort -i A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.gff > A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.sorted.gff
# Convert gff to bed
gff2bed < A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.sorted.gff > A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.bed
# Fetch chrom sizes
twoBitInfo genome.2bit chrom.sizes
# Convert bed to bigbed
bedToBigBed -tab -type=bed6+4 -as=gff3.as A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.bed chrom.sizes A-AFFY-87_AffyGeneChipTomatoGenome.probes_ITAG2.3genome_mapping.bb
```

The [gff3.as](https://github.com/nlesc-ave/ave-rest-service/blob/master/gff3.as) (in Docker container available as `/app/gff3.as`) is used to describe the columns in the bed file.

A BigBed file with more than 256 chromosomes will not render, see [issue 37](https://github.com/nlesc-ave/ave-app/issues/37).

## Data registration

After deployment the server is running, but contains no data.
Data files must be registered so they show up in the web browser.

The following commands expect a running Docker container as described in the [Deployment](#deployment) chapter.

### Register genome

```bash
docker exec ave \
    avedata register \
    --species 'Solanum Lycopersicum' \
    --genome SL.2.40 \
    --datatype 2bit \
    /data/tomato/SL.2.40/genome.2bit
```

The last argument is the location of the genome 2bit file (`/data/tomato/SL.2.40/genome.2bit` in this example), it must be an absolute path which starts with `/data/` must be readable by anyone inside the Docker container.

The [ave-app](https://github.com/nlesc-ave/ave-app) front end will use this path as the relative http(s) path to fetch reference genome sequence in the selected region
and the AVE rest service will use it to determine the chromosome list and build haplotype sequence.

### Register variants

```bash
docker exec ave \
    avedata register \
    --species 'Solanum Lycopersicum' \
    --genome SL.2.40 \
    --datatype variants \
    /data/tomato/SL.2.40/tomato_snps.bcf
```

The `/data/tomato/SL.2.40/tomato_snps.bcf` file must be readable by anyone inside the Docker container.

To perform clustering a registered genome (2bit) file and corresponding variant (bcf) file is required.

### Register genes

One genome can have one gene track. The gene track shows exons, introns and untranslated regions.

```bash
docker exec ave \
    avedata register \
    --species 'Solanum Lycopersicum' \
    --genome SL.2.40 \
    --datatype genes \
    /data/tomato/SL.2.40/gene_models.bb
```

The last argument is the bigbed formatted file with genes (`/data/tomato/SL.2.40/gene_models.bb` in this example), it must be an absolute path which starts with `/data/` and must be readable by anyone inside the Docker container.

Registration can take some time because a [Whoosh](https://whoosh.readthedocs.io) full text index is build.

The [ave-app](https://github.com/nlesc-ave/ave-app) front end will use this path as the relative http(s) path to fetch the genes in the selected region.

### Register feature annotations

One genome can have one ore more feature annotation files registered.

```bash
docker exec ave \
    avedata register \
    --species 'Solanum Lycopersicum' \
    --genome SL.2.40 \
    --datatype features \
    /data/tomato/SL.2.40/A-AFFY-87.bb
```

The last argument is the bigbed formatted file with feature annotations (`/data/tomato/SL.2.40/A-AFFY-87.bb` in this example), it must be an absolute path which starts with `/data/` and must be readable by anyone inside the Docker container.

Registration can take some time because a [Whoosh](https://whoosh.readthedocs.io) full text index is build.

The [ave-app](https://github.com/nlesc-ave/ave-app) front end will use this path as the relative http(s) path to fetch the feature annotations in the selected region.
The basename of the file, in this case `A-AFFY-87`, will be used as the track label.

### Deregister

To deregister one of the files use the deregister command.

For example to deregister the `/data/tomato/SL.2.40/A-AFFY-87.bb` file use:
```bash
docker exec ave \
    avedata deregister \
    /data/tomato/SL.2.40/A-AFFY-87.bb
```

## Develop

Below are instructions how to get a development version of the service up and running.

Requirements:

* [Anaconda3](https://www.continuum.io/downloads) or [miniconda3](https://conda.io/miniconda.html)

### Setup

First clone the repository
```bash
git clone https://github.com/nlesc-ave/ave-rest-service.git
cd ave-rest-service
```

To create a new Anaconda environment with all the ave dependencies installed.
```bash
conda env create -f environment.yml
```
On osx use `enviroment.osx.yml` instead of `environment.yml`.

Activate the environment
```bash
source activate ave2
```

Install ave for development with
```bash
python setup.py develop
```

If dependencies are changed in `environment.yml` then update conda env by runnning
```
conda env update -f environment.yml
```

### Configuration

The service needs a configuration file called `setting.cfg`.

The repo contains an example config file called `settings.example.cfg`.
Copy the example config file to `settings.cfg` and edit it.

Make sure the `WHOOSH_BASE_DIR` directory exists.

### Run service

Change to directory with `settings.cfg` file.

```bash
gunicorn -w 4 --threads 2 -t 60 -b 127.0.0.1:8080 avedata.app:app
```

It will run the service on http://127.0.0.1:8080/ .
The api endpoint is at `http://127.0.0.1:8080/api/`.
The Swagger UI is at `http://127.0.0.1:8080/api/ui`.

This will only run the Python web service, the hosting of the application and data files is explained in the [deployment](#deployment) chapter.

### Run commands

The `avedata` command line tool expects to be run from the directory containing the `settings.cfg` file.

Available commands:
* deregister  Remove a filename or url from the database
* drop_db     Drops database
* init_db     Initializes database
* list        List of registered files/urls in the database
* register    Add file metadata information to the database
* run         Run as single threaded web service

### Build Docker image

Building the Docker image of the latest commit of the master branch is automatically build on Docker hub with the latest tag.

If you have local changes and want to test the Docker container locally then you can build the Docker image with
```bash
docker build -t ave2/allelic-variation-explorer .
```

The Docker image contains the latest version of [ave-app](https://github.com/nlesc-ave/ave-app).
If you want to run a different version of the app you will need to replace curl/tar command in `./Dockerfile` with commands to put the app you want in `/var/www/html` directory and build a Docker image.
