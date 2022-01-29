<p align="center">
   <img width="450" src="https://raw.githubusercontent.com/tobiasrausch/alfred/master/alfred.png">
   <h1></h1>
</p>

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/alfred/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/alfred/badges/downloads.svg)](https://anaconda.org/bioconda/alfred)
[![C/C++ CI](https://github.com/tobiasrausch/alfred/workflows/C/C++%20CI/badge.svg)](https://github.com/tobiasrausch/alfred/actions)
[![Docker CI](https://github.com/tobiasrausch/alfred/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/trausch/alfred/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/tobiasrausch/alfred/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/tobiasrausch/alfred.svg)](https://github.com/tobiasrausch/alfred/releases)

## Alfred: BAM alignment statistics, feature counting and feature annotation

Alfred is available as a [Bioconda package](https://anaconda.org/bioconda/alfred), as a pre-compiled statically linked binary from [Alfred's github release page](https://github.com/tobiasrausch/alfred/releases/), as a singularity container [SIF file](https://github.com/tobiasrausch/alfred/releases/) or as a minimal [Docker container](https://hub.docker.com/r/trausch/alfred/). Please have a look at [Alfred's documentation](https://www.gear-genomics.com/docs/alfred/) for any installation or usage questions.

[Source Code](https://github.com/tobiasrausch/alfred/)

[Web Application](https://www.gear-genomics.com/alfred/)

[Documentation](https://www.gear-genomics.com/docs/alfred/)

## Citation

Tobias Rausch, Markus Hsi-Yang Fritz, Jan O Korbel, Vladimir Benes.      
Alfred: interactive multi-sample BAM alignment statistics, feature counting and feature annotation for long- and short-read sequencing.
Bioinformatics. 2019 Jul 15;35(14):2489-2491.
[https://doi.org/10.1093/bioinformatics/bty1007](https://doi.org/10.1093/bioinformatics/bty1007)


License
-------
Alfred is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/tobiasrausch/alfred/blob/master/LICENSE) file for more details.
Contributing
------------
Thank you for considering contributing to Alfred.

Bugs
----
If you have noticed a bug in Alfred please search the issue tracker to see if someone else in the community had already created a ticket. If that is not the case please create one.

We try to fix all bugs as soon as possible but if you want to work on one yourself please consider using the fork and pull request mechanism of GitHub.

Usage questions
---------------
Extensive documentation is available at [https://www.gear-genomics.com/docs/alfred/](https://www.gear-genomics.com/docs/alfred/).

New Features
------------
If you want to suggest a new feature please go ahead and open an issue.
## Alfred: Bioconda installation test

Simple test makefile that pulls miniconda, adds the bioconda channel and installs alfred. It then displays the Alfred help message.

`make all`
You can build an [alfred](https://github.com/tobiasrausch/alfred) singularity container (SIF file) using

`sudo singularity build alfred.sif alfred.def`

Once you have built the container you can run analysis using

`singularity exec alfred.sif alfred qc -r ref.fa input.bam`
![Alfred logo](./images/alfred.png)

# Alfred: BAM Statistics, Feature Counting and Feature Annotation

[Alfred](https://github.com/tobiasrausch/alfred) is an efficient and versatile command-line application that computes multi-sample quality control metrics in a read-group aware manner. Alfred supports read counting, feature annotation and haplotype-resolved consensus computation using multiple sequence alignments. Alfred's [companion web application](https://www.gear-genomics.com/alfred/) enables interactive exploration of results. All code is open-source and hosted on [Alfred's GitHub page](https://github.com/tobiasrausch/alfred).

Contents:

1. [Installation](/installation/)
2. [Usage](/cli/)
3. [Web Application](/webapp/)
4. [FAQ](/faq/)

::: tip
For questions, help or feature requests please contact gear_genomics@embl.de
:::

Please cite Alfred's URL (https://www.gear-genomics.com/alfred) in publications.
# Usage

Alfred uses subcommands for [quality control](#alignment-quality-control) ([qc](#alignment-quality-control)), [feature counting](#bam-feature-counting) ([count_dna](#bam-read-counting-for-dna-seq), [count_rna](#bam-read-counting-for-rna-seq), [count_jct](#bam-read-counting-for-rna-seq)), [feature annotation](#bam-feature-annotation) ([annotate](#chip-seq-or-atac-seq-peak-annotation), [tracks](#browser-tracks)), [alignment](#pairwise-sequence-alignment) ([pwalign](#pairwise-sequence-alignment), [consensus](#bam-consensus-computation)) and [haplotype-resolved analysis](#haplotype-specific-bam-files) ([split](#haplotype-specific-bam-files), [ase](#allele-specific-expression)). The subcommands are explained below.

## Alignment Quality Control

Alfred supports a command-line interface to run alignment quality control and a [web application](https://www.gear-genomics.com) can be used to render all QC metrics.

## Command-line Interfact for BAM Quality Control

Alfred computes various alignment metrics and summary statistics by read group

```bash
alfred qc -r <ref.fa> -o qc.tsv.gz <align.bam>
```

Plotting alignment statistics

```bash
Rscript scripts/stats.R qc.tsv.gz
```

To convert all the alignment metrics from column format to row format for readability

```bash
zgrep ^ME qc.tsv.gz | cut -f 2- | datamash transpose | column -t
```

## Interactive Quality Control Browser

Quality control metrics can be browsed interactively using the [web front end of Alfred](https://www.gear-genomics.com/alfred).

```bash
alfred qc -r <ref.fa> -j qc.json.gz -o qc.tsv.gz <align.bam>
```

Then just upload the qc.json.gz file to the Alfred GUI [https://www.gear-genomics.com/alfred](https://www.gear-genomics.com/alfred). A convenient feature of the web-front end is that multiple samples can be uploaded and compared.


## BAM Alignment Quality Control for Targeted Sequencing

If target regions are provided, Alfred computes the average coverage for each target and the on-target rate.

```bash
alfred qc -r <ref.fa> -b <targets.bed.gz> -o <qc.tsv.gz> <align.bam>
```

For instance, for a human whole-exome data set.

```bash
cd maps/ && Rscript exon.R
alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -j qc.json.gz -o qc.tsv.gz <exome.bam>
Rscript scripts/stats.R qc.tsv.gz
```

Alternatively, one can use the [interactive GUI](https://www.gear-genomics.com/alfred) and upload the json file.

```bash
alfred qc -r <hg19.fa> -b maps/exonic.hg19.bed.gz -j qc.json.gz -o qc.tsv.gz <exome.bam>
```


## BAM Alignment Quality Control for ATAC-Seq

For ATAC-Seq data, the insert size distribution should reveal the DNA pitch and a clear nucleosome pattern with a peak for single nucleosomes and dimers. The transcription start site (TSS) enrichment should be >5 for a good ATAC-Seq library and ideally the duplicate rate is <20%, the alignment rate >70% and the standardized SD in coverage >0.3.

```bash
cd maps/ && Rscript promoter.R
alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -o qc.tsv.gz <atac.bam>
Rscript scripts/stats.R qc.tsv.gz
zgrep ^ME qc.tsv.gz | datamash transpose | egrep "^Dup|^MappedFraction|^SD|^Enrich"
```

ATAC-Seq libraries often have a large number of mitochondrial reads depending on the library preparation.

```bash
zgrep ^CM qc.tsv.gz | egrep "Mapped|chrM"
```

Alternatively, one can use the [interactive GUI](https://www.gear-genomics.com/alfred) and upload the json file.

```bash
alfred qc -r <hg19.fa> -b maps/hg19.promoter.bed.gz -j qc.json.gz -o qc.json.gz <atac.bam>
```

## BAM Feature Counting

Alfred supports counting reads in overlapping or non-overlapping windows, at predefined intervals in BED format,
or as gene and transcript counting for RNA-Seq in stranded or unstranded mode using a gtf or gff3 gene annotation
file. Expression values can be normalized as raw counts, FPKM, or FPKM-UQ values.

## BAM Read Counting for DNA-Seq

For DNA sequencing, Alfred can be used to calculate the coverage in overlapping or non-overlapping windows or in given set of intervals.

```bash
alfred count_dna -o <cov.gz> <align.GRCh37.bam>
```

To plot the whole-chromosome coverage profile for chr1-22 and chrX.

```bash
Rscript scripts/rd.R <cov.gz>
```

## BAM Read Counting for RNA-Seq

Alfred can also assign reads to gene annotation features from a GTF file such as counting reads by gene or transcript identifier.

```bash
cd gtf/ && ./downloadGTF.sh
alfred count_rna -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <align.GRCh37.bam>
```

An experimental feature of Alfred is to count splice junction supporting reads. This method generates exon-exon junction counts for intra-gene exon-exon junctions, inter-gene exon-exon junctions and completely novel (not annotated) intra-chromosomal junctions.

```bash
alfred count_jct -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <align.GRCh37.bam>
```

## BAM Feature Annotation

Alfred supports annotation of ChIP-Seq and ATAC-Seq peaks for neighboring genes or transcription factor binding sites (based on motif alignments). Additionally, browser tracks in UCSC bedgraph format can be computed with configurable resolution.

## ChIP-Seq or ATAC-Seq peak annotation

To annotate overlapping/neighboring genes up to a distance of 10,000bp:

```bash
alfred annotate -d 10000 -g gtf/Homo_sapiens.GRCh37.75.gtf.gz <peaks.bed>
```

The two output files summarize nearby genes by peak and vice versa (peaks by gene).


## Motif annotation

Motif annotation of peaks based on alignment scores of motif-specific position weight matrices can also be obtained. 

```bash
cd motif/ && ./downloadMotifs.sh
alfred annotate -r <hg19.fa> -m motif/jaspar.gz <peaks.bed>
```

Alfred further implements functionality to output the exact motif hits in each peak to perform for instance transcription factor footprinting.

```bash
alfred annotate -p hits.gz -r <hg19.fa> -m motif/jaspar.gz <peaks.bed>
```


## Joined motif hits

You can also obtain all motif hits across a reference genome. For instance, to identify all hits for recombination signal sequences (RSS) you can use:

```bash
cat <hg19.fa.fai>  | cut -f 1,2 | awk '{print $1"\t0\t"$2;}' > all.chrom.bed
alfred annotate -m motif/rss.gz -r <hg19.fa> -q 0.9 -p motif.hits.gz all.chrom.bed
```

To then join all motif hits of the conserved heptamer sequence (7bp), a spacer sequence (12bp or 23bp), and a conserved nonamer sequence (9bp) you can use:

```bash
alfred spaced_motif -l 11 -h 13 motif.hits.gz
```



## Browser Tracks

Normalized and file size optimized browser tracks are essential to compare peak profiles across samples and upload them quickly in online genome browsers such as the [UCSC Genome Browser](https://genome.ucsc.edu/). Alfred generates browser tracks in [UCSC bedgraph format](https://genome.ucsc.edu/goldenpath/help/bedgraph.html) with configurable resolution. Lower resolution values leed to coarse-grained windows at reduced file size. Contrarily, high resolution values leed to fine-grained windows at greater file size.

```bash
alfred tracks -r 0.2 -o track.bedGraph.gz <aligned.bam>
```

[IGV tools](https://software.broadinstitute.org/software/igv/igvtools) can be used to convert bedgraph files to IGV's proprietory tdf format.

```bash
igvtools totdf track.bedGraph.gz track.tdf hg19
```

This conversion enables comparing dozens of ATAC-Seq or ChIP-Seq samples at greatly improved speed in IGV compared to using the raw BAM files. By default the Alfred tracks command normalizes all input files to 30 million reads so peak heights are directly comparable across samples.

## Pairwise Sequence Alignment

Alfred supports global and local pairwise sequence alignments with configurable match scores and gap/mismatch penalties. Overlap and nested alignments can be achieved using penalty-free leading/trailing gaps. Affine gap penalties with separate gap open and gap extension penalties are supported.

```bash
alfred pwalign -a align.fa.gz <seq1.fasta> <seq2.fasta>
```

All computed pairwise alignments are "linear" alignments that is the order of sequence nucleotides is preserved. For more complex pairwise alignments involving inversions or duplications you can use [Maze](https://www.gear-genomics.com/maze/).

## BAM Consensus Computation

Alfred supports consensus computation of error-prone long reads using multiple sequence alignment principles. To compute a consensus sequence for a set of long reads at a given alignment position:

```bash
alfred consensus -f bam -t ont -p chr1:218992200 <ont_pacbio.bam>
```

The consensus method generates two output files. A simple fasta file with the consensus sequence and a FASTA align file that shows the multiple sequence alignment used for consensus generation in either horizontal or vertical format. The horizontal format is the classical FASTA align format. The vertical format transposes the horizontal alignment and in addition shows the consensus nucleotide for each alignment column after the vertical bar.

The consensus command is probably most useful by first separating haplotypes (Alfred's [split](#haplotype-specific-bam-files) subcommand) and then computing consensus sequences independently for each haplotype.

```bash
alfred split -r <ref.fa> -s NA12878 -p <haplotype1.bam> -q <haplotype2.bam> -v <phased.snps.bcf> <input.bam>
alfred consensus -c <hap1.fa.gz> -f bam -t ont -p chr1:chr4:500500 <haplotype1.bam>
alfred consensus -c <hap2.fa.gz> -f bam -t ont -p chr1:chr4:500500 <haplotype2.bam>
```

To identify variants, you can then compare the two locally assembled haplotypes using our online dotplot method [Maze](https://www.gear-genomics.com/maze) or align them pairwise against each other using the subcommand [pwalign](#pairwise-sequence-alignment).


## Haplotype-specific BAM files

Long reads, Hi-C or novel single-cell sequencing technologies such as Strand-Seq enable (local) haplotyping or phasing of variants. Once such a phased SNP scaffold has been obtained one can split a BAM file into the corresponding haplotypes.

```bash
alfred split -r <ref.fa> -s NA12878 -p <haplotype1.bam> -q <haplotype2.bam> -v <phased.snps.bcf> <input.bam>
```

## Allele-specific expression

Alfred implements methods to generate allele-specific expression tables. If the input SNPs are phased Alfred annotates the allele-specific expression haplotype-aware.

```bash
alfred ase -r <ref.fa> -s NA12878 -v <snps.bcf> -a <ase.tsv> <input.bam>
```

## Analysis of replication timing by NGS

Alfred implements a method to analyze replication timing using next-generation sequencing (Repli-Seq). The order of BAM files on the command-line must follow the cell-cycle.

```bash
alfred replication -r <ref.fa> -o outprefix <g1b.bam> <s1.bam> <s2.bam> <s3.bam> <s4.bam> <g2.bam>
```

There is a supporting script that plots the tag density for each cell-cycle fraction.

```bash
Rscript scripts/reppattern.R -f outprefix.profile.tsv -r chr12:24000000-26000000
```

There is also a script for plotting the replication time along a given chromosome.

```bash
Rscript scripts/reptime.R -f outprefix.reptime.tsv -r chr12
```
# Installation

Alfred is available as a [Bioconda package](https://anaconda.org/bioconda/alfred), as a pre-compiled statically linked binary from [Alfred's github release page](https://github.com/tobiasrausch/alfred/releases/) or as a minimal [Docker container](https://hub.docker.com/r/trausch/alfred/). All code is open-source and hosted on [Alfred's GitHub page](https://github.com/tobiasrausch/alfred).

## Installation from Source

To build Alfred from source you need some build essentials and the Boost libraries, i.e. for Ubuntu:

```bash
apt install \
    build-essential g++ \
    cmake \
    git-all \
    liblzma-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev
```

Once you have installed these system libraries you can compile and link Alfred.

```bash
git clone --recursive https://github.com/tobiasrausch/alfred.git
cd alfred/
make all
make install
./bin/alfred -h
```

## Custom Boost installation directory

Alfred requires Boost and you can install Boost locally using

```bash
git clone --recursive https://github.com/boostorg/boost.git
cd boost/
./bootstrap.sh --prefix=`pwd` --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time
./b2
./b2 headers
cd ..
```

You can then specify the custom Boost installation directory (i.e., `/opt/boost` below) using

```bash
# modify the following line accordingly
BOOSTROOT=/opt/boost
git clone --recursive https://github.com/tobiasrausch/alfred.git
cd alfred/
make CMDCXXFLAGS="-isystem $BOOSTROOT" CMDLDFLAGS="-L$BOOSTROOT/stage/lib -Wl,-rpath,$BOOSTROOT/stage/lib" all
make install
./bin/alfred -h
```
# Web application

Alfred's quality control JSON files can be interactively browsed with the
[companion web application](https://www.gear-genomics.com/alfred).
All charts support panning and zooming and can be downloaded as PNG images.
The summary QC table can be downloaded as a CSV file.

To generate a quality control file in JSON format run [Alfred's command-line tool](/cli/) as follows:

```bash
alfred qc -r <ref.fa> -f json -o qc.json.gz <align.bam>
```

The output file `qc.json.gz` can then be uploaded at
[https://www.gear-genomics.com/alfred/](https://www.gear-genomics.com/alfred/).

## Features

An overview of all available charts and the most important alignment statistics provided by Alfred is below.

| Alignment Metric               | DNA-Seq (WGS) | DNA-Seq (Capture) | RNA-Seq | ChIP-Seq/ATAC-Seq | Chart Type         |
| ------------------------------ | ------------- | ----------------- | ------- | ----------------- | ------------------ |
| Mapping Statistics             | ✔             | ✔                 | ✔       | ✔                 | Table              |
| Duplicate Statistics           | ✔             | ✔                 | ✔       | ✔                 | Table              |
| Sequencing Error Rates         | ✔             | ✔                 | ✔       | ✔                 | Table              |
| Base Content Distribution      | ✔             | ✔                 | ✔       | ✔                 | Grouped Line Chart |
| Read Length Distribution       | ✔             | ✔                 | ✔       | ✔                 | Line Chart         |
| Base Quality Distribution      | ✔             | ✔                 | ✔       | ✔                 | Line Chart         |
| Coverage Histogram             | ✔             | ✔                 | ✔       | ✔                 | Line Chart         |
| Insert Size Distribution       | ✔             | ✔                 | ✔       | ✔                 | Grouped Line Chart |
| InDel Size Distribution        | ✔             | ✔                 | ✔       | ✔                 | Grouped Line Chart |
| InDel Context                  | ✔             | ✔                 | ✔       | ✔                 | Bar Chart          |
| GC Content                     | ✔             | ✔                 | ✔       | ✔                 | Grouped Line Chart |
| On-Target Rate                 |               | ✔                 |         |                   | Line Chart         |
| Target Coverage Distribution   |               | ✔                 |         |                   | Line Chart         |
| TSS Enrichment                 |               |                   |         | ✔                 | Table              |
| DNA pitch / Nucleosome pattern |               |                   |         | ✔                 | Grouped Line Chart |

## Base content distribution

The base content distribution shows any base calling bias along the read.
For an ideal library the lines for A, C, G, and T should run in parallel.
For a whole-genome assay the GC-content of that genome should be reflected in the relative amounts of each base.
Some libraries are expected to show a biased leading base distribution such as many RNA-Seq libraries
because of random hexamer priming or restriction based assays.

## Read length distribution

Illumina sequencers produce reads of fixed read length but long read technologies usually have a median read length >1000bp
and a long tail of reads with read lengths >30,000bp. This plot is also insightful to understand adapter trimming results
or the removal of low quality bases at the start or end of a read.

## Mean base quality distribution

This plot shows the mean base quality along the read. A typical Illumina profile shows base qualities >Q30
before base 30 and then a gradual loss of base quality accuracy towards the end of the read.

## Mapping quality distribution

This plot shows the mapping quality distribution for all mapped reads. The reported quality scores are aligner-dependent.

## Coverage histogram

The coverage histogram shows how many bases of the sequenced genome are at a given coverage.
Please note that for targeted assays (capture assays) this plot is expected to show a large portion of the genome at coverage=0.
For targeted assays, we therefore recommend checking the on-target rate and the targets above coverage level plots.

## On-target rate and targets above a given coverage level

For targeted assays, the two major concerns are capture efficiency (on-target rate)
and how many of the targets are ready for downstream analysis
(targets above a pre-defined coverage threshold).
A standard whole-exome sequencing assay yields at least 70% of reads on-target
(+/-200bp target extension) and at least 70% of targets >20x coverage.

## Insert size histogram

The insert size plot shows the outer insert size distribution for all read pairs
stratified by read pair orientation. There are different nomenclatures around for
defining the different paired-end layouts. The default Illumina paired-end layout is R+
(or forward-reverse, FR), the default Illumina mate-pair layout is R- (or reverse-forward, RF).
For specific sequencing assays, the insert size distribution can serve as a key quality control metric.
For instance, ATAC-Seq libraries should show the characteristic nucleosome pattern and DNA pitch.

## InDel size distribution

Histogram of indel sizes collected from all mapped reads. This plot aggregates the length
of all Cigar `I` and `D` operations.

## InDel Homopolymer Context

The homopolymer plot shows for all InDels (Cigar I and D operations) if the preceding 3 bases are
all A, all C, all G, or all T. If at least 2 different nucleotides are observed the reported
homopolymer context is "None". For Illumina reads, almost 50% of all reported InDels occur in a
homopolymer context with greater counts for A and T compared to G and C.

## GC content

To estimate a GC bias curve even for low-coverage single-cell data, Alfred computes for each mapped read
the local GC-content and then compares the estimated sample GC content to the expected, genome-wide GC content.
If a targeted assay is analyzed, Alfred, in addition, computes the GC content of all target regions.

## GC-Content and Mapping Statistics by Chromosome

This table lists the size, the number of Ns, the GC-content, and the number of mapped reads for each chromosome
as well as the observed-to-expected ratio of mapped reads.

## Summary statistics

The summary tab aggregates quality control data in a simple table that can be downloaded in CSV format.
This table is ideal to compare QC metrics across samples and/or sequencing assays. Among many other statistics,
the table lists, for instance, the number of duplicate reads, the number of unmapped reads, the number of
secondary and supplementary alignments, base-pair exact error rates stratified by mismatch,
insertion and deletion errors, and the median coverage and insert size of the sequenced sample.
The table provides more detailed statistics for specialized assays, i.e.
for 10X Genomics it lists the number of MI tagged reads, the total number of UMIs,
the fraction of haplotype-tagged reads and the N50 phased block length.
For ATAC-Seq data, users can provide a BED file of promoter regions and then the `EnrichmentOverBed` column
corresponds to TSS enrichment whereas for WES data, the enrichment quantifies the capturing efficiency
if the BED file contains all target regions.

## Example Data Sets

The [web application](https://www.gear-genomics.com/alfred) hosts example data sets for a number of sequencing assays and sequencing technologies.

| Sequencing Assay  | Sequencing Technology |
| ----------------- | --------------------- |
| DNA-Seq (WGS)     | Illumina, PacBio, ONT |
| DNA-Seq (Capture) | Illumina              |
| RNA-Seq           | Illumina              |
| ATAC-Seq          | Illumina              |
| ChIP-Seq          | Illumina              |
# FAQ

::: tip
For questions, help or feature requests please contact gear_genomics@embl.de
:::

[[toc]]

## Does Alfred support the CRAM format?

Yes, Alfred uses [HTSlib](https://github.com/samtools/htslib) to read/write BAM/CRAM files.

## Is there an example data set to test my Alfred installation?

The github source code includes a minimal example to check that alfred compiled properly from source and that the web front end is working.

```bash
alfred qc -r example/E.coli.fa.gz -o example/stats.tsv.gz example/E.coli.cram
Rscript scripts/stats.R example/stats.tsv.gz
```

For the web front end.

```bash
alfred qc -r example/E.coli.fa.gz -f json -o ecoli.json.gz example/E.coli.cram
```

Please upload `ecoli.json.gz` to the [Alfred web application](https://www.gear-genomics.com/alfred).

## Is the feature counting paired-end aware?

Yes, Alfred counts fragments (read pairs) and not individual reads.

## Why are hard clipping statistics always zero?

Many aligners trim primary alignments using soft-clips and only secondary and supplementary alignments use hard clips. For long reads you may want to evaluate secondary and supplementary alignments using the `-s` and `-u` command-line flags.

```bash
alfred qc -su -r <genome.fa> <input.bam>
```

## Calculation of InDel rates and sequencing error rate?

The sequencing error rates are calculated over all aligned bases. The total number of deletions, character D in the Cigar string, is divided
by the total number of aligned bases (Cigar M). The same approach is used for insertion (Cigar I) and mismatches (Cigar M and mismatch to reference).
