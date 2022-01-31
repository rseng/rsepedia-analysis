## Quick start
```shell
    tandemquast.py --nano test_data/simulated_reads.fasta test_data/simulated_polished.fa -o simulated_res
```

## Introduction

**TandemTools** package includes **TandemQUAST** tool for evaluating and improving assemblies of extra-long tandem repeats (ETR) and **TandemMapper** tool for mapping long error-prone reads to ETRs.

Note: TandemTools is designed specifically for ETR (range in length from hundreds of thousands to millions of nucleotides). It is strongly not recommended to run TandemTools on shorter TRs.


## Installation

Requirements are listed in ```requirements.txt``` and can be installed through Conda as ```conda install --file requirements.txt```.

## Usage

```shell
tandemquast.py [options] --nano/--pacbio <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    --nano         PATH                 File with Oxford Nanopore reads used for ETR assembly 
      or
    --pacbio       PATH                 File with PacBio CLR reads used for ETR assembly
    -o             PATH                 Folder to store all result files

Optional arguments:    
    -t             INT                  Maximum number of threads [default: 4]
    -l             \"label,label,...\"  Human-readable names of assemblies to use in reports, comma-separated. If contain spaces, use quotes 
    -m             PATH                 FASTA file with monomer sequences
    --hifi         PATH                 File with accurate PacBio HiFi reads
    --only-polish                       Run tandemQUAST polishing module (metrics will not be calculated).
```

```shell
tandemmapper.py [options] --nano/--pacbio <reads_file> -o <output_dir> <assembly_file1> <assembly_file2>

Required arguments:
    --nano         PATH                 File with Oxford Nanopore reads used for ETR assembly 
      or
    --pacbio       PATH                 File with PacBio CLR reads used for ETR assembly
    -o             PATH                 Folder to store alignment results

Optional arguments:    
    -t             INT                  Maximum number of threads [default: 4]
```

## Output files

The following files are contained in `<output_dir>` directory (specified by `-o`) and include results 
for all input assemblies. See TandemQUAST paper for the detailed explanation of different metrics.

### TandemMapper output

`<output_dir>/*_alignment.bed` - TandemMapper alignments in BED format

`<output_dir>/*_alignment.sam` - TandemMapper alignments in SAM format

### TandemQUAST output

#### Basic metrics:

`<output_dir>/report/*_coverage.png` - plot of read coverage.

`<output_dir>/report/*_bp_analysis.png` - plot of breakpoint ratio values
 (see TandemQUAST paper for the details). The red peaks in this plot may correspond
 to large-scale assembly errors.

`<output_dir>/report/*_kmer_analysis.png` - distribution of different types of 
unique k-mers along the assembly. Each bar 
shows the number of different types of k-mers in a bin of length 20 kb. The blue bars
represent single-clump k-mers. The high number of single-clump k-mers suggests the good 
base-level quality of the assembly. The orange (multiple-clumps) and green (no-clumps) bars
suggest a low base-level quality in the region, caused by lack of polishing or an assembly error.

`<output_dir>/report/*_kmer_stats.txt` - distribution of different types of unique k-mers in TXT format.

#### Pairwise comparison:

`<output_dir>/report/*_vs_*.png` - a dot plot comparing mappings for two assemblies

`<output_dir>/report/discordance_*_vs_*.png` - a discordance plot showing coverage of two assemblies by discordant reads.
The peaks of coverage for one assembly suggest that this assembly is more "supported" by reads in this region than other assembly. 

#### Centromeric metrics:

`<output_dir>/report/*_monomer_lengths.html` - an interactive HTML-page showing 
monomer length distribution along the assembly.

`<output_dir>/report/*_units.txt` - file with a list of HOR units.

## Citation

Alla Mikheenko, Andrey V. Bzikadze, Alexey Gurevich, Karen H. Miga, and Pavel A. Pevzner. TandemMapper and TandemQUAST: mapping long reads and assessing/improving assembly quality in extra-long tandem repeats, 2019, bioRxiv


## Contacts

Please report any problems to the [issue tracker](https://github.com/ablab/tandemQUAST/issues). Alternatively, you can write directly to [a.mikheenko@spbu.ru](mailto:a.mikheenko@spbu.ru).Flye assembler
==============

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/flye.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/flye)

### Version: 2.5

Flye is a de novo assembler for single molecule sequencing reads,
such as those produced by PacBio and Oxford Nanopore Technologies.
It is designed for a wide range of datasets, from small bacterial projects
to large mammalian-scale assemblies. The package represents a complete
pipeline: it takes raw PB / ONT reads as input and outputs polished contigs.
Flye also includes a special mode for metagenome assembly.

Latest updates
--------------

### Flye 2.5 release (25 Jul 2019)
* Better ONT polishing for the latest basecallers (Guppy/flipflop)
* Improved consensus quality of repetitive regions
* More contiguous assemblies of real metagenomes
* Improvements for human genome assemblies
* Various bugfixes and performance optimizations
* Also check the new [FAQ section](docs/FAQ.md)

Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)
- [FAQ](docs/FAQ.md)


Repeat graph
------------

Flye is using repeat graph as a core data structure. 
In difference to de Bruijn graphs (which require exact k-mer matches),
repeat graphs are built using approximate sequence matches, and
can tolerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define
the junctions. Each edges is classified into unique or repetitive.
The genome traverses the graph (in an unknown way), so as each unique
edge appears exactly once in this traversal. Repeat graphs reveal the
repeat structure of the genome, which helps to reconstruct an optimal assembly.


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

Above is an example of the repeat graph of a bacterial assembly.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. Note that each edge is represented in 
two copies: forward and reverse complement (marked with +/- signs), 
therefore the entire genome is represented in two copies. This is necessary
because the orientation of input reads is unknown.

In this example, there are two unresolved repeats: (i) a red repeat of 
multiplicity two and length 35k and (ii) a green repeat cluster of multiplicity
three and length 34k - 36k. As the repeats remained unresolved, there are no reads
in the dataset that cover those repeats in full. Five unique edges 
will correspond to five contigs in the final assembly.

Repeat graphs produced by Flye could be visualized using
[AGB](https://github.com/almiheenko/AGB) or [Bandage](https://github.com/rrwick/Bandage).


Flye benchmarks
---------------

| Genome                   | Data           | Asm.Size  | NG50     | CPU time  | RAM    |
|--------------------------|----------------|-----------|----------|-----------|--------|
| [E.coli][ecoli]          | PB 50x         | 4.6 Mb    | 4.6 Mb   | 2 h       | 2 Gb   |
| [C.elegans][ce]          | PB 40x         | 102 Mb    | 2.9 Mb   | 100 h     | 31 Gb  |
| [A.thaliana][at]         | PB 75x         | 120 Mb    | 10.7 Mb  | 100 h     | 46 Gb  |
| [D.melanogaster][dm-ont] | ONT 30x        | 139 Mb    | 17.5 Mb  | 130 h     | 31 Gb  |     
| [D.melanogaster][dm-pb]  | PB 120x        | 142 Mb    | 17.5 Mb  | 150 h     | 75 Gb  |     
| [Human NA12878][na12878] | ONT 35x (rel6) | 2.9 Gb    | 22.6 Mb  | 2500 h    | 714 Gb |
| [Human CHM13 T2T][t2t]   | ONT 50x (rel2) | 2.9 Gb    | 57.9 Mb  | 3600 h    | 871 Gb |
| [Human HG002][hg002]     | PB CCS 30x     | 2.9 Gb    | 30.4 Mb  | 1400 h    | 272 Gb |
| [Human CHM1][chm1]       | PB 100x        | 2.8 Gb    | 18.8 Mb  | 2700 h    | 676 Gb |
| [HMP mock][hmp]          | PB meta 7 Gb   | 66 Mb     | 2.6 Mb   | 60 h      | 72 Gb  |
| [Zymo Even][zymo]        | ONT meta 14 Gb | 64 Mb     | 0.6 Mb   | 60 h      | 129 Gb |
| [Zymo Log][zymo]         | ONT meta 16 Gb | 23 Mb     | 1.3 Mb   | 100 h     | 76 Gb  |

[na12878]: https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[at]: https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/
[dm-pb]: https://github.com/PacificBiosciences/DevNet/wiki/Drosophila-sequence-and-assembly
[dm-ont]: https://www.ebi.ac.uk/ena/data/view/SRR6702603
[hg002]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/
[ecoli]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
[hmp]: https://github.com/PacificBiosciences/DevNet/wiki/Human_Microbiome_Project_MockB_Shotgun 
[chm1]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP044331
[t2t]: https://github.com/nanopore-wgs-consortium/CHM13
[zymo]: https://github.com/LomanLab/mockcommunity

The assemblies generated using Flye 2.5 could be downloaded from [Zenodo](https://zenodo.org/record/3353665).
All datasets were run with default parameters with the following exceptions:
CHM13 T2T was run with `--min-overlap 10000`; CHM1 was run with `--asm-overage 40`;
HG002 was run with maximum read error rate set to 1%.

Third-party
-----------

Flye package includes some third-party software:

* [libcuckoo](http://github.com/efficient/libcuckoo)
* [intervaltree](https://github.com/ekg/intervaltree)
* [lemon](http://lemon.cs.elte.hu/trac/lemon)
* [minimap2](https://github.com/lh3/minimap2)


License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Flye is developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Repeat graph and current package maintaining: Mikhail Kolmogorov
* Trestle module and original polisher code: Jeffrey Yuan
* Original contig extension code: Yu Lin
* Short plasmids recovery module: Evgeny Polevikov


Publications
------------
Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", Nature Biotechnology, 2019
[doi:10.1038/s41587-019-0072-8](https://doi.org/10.1038/s41587-019-0072-8)

Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS, 2016
[doi:10.1073/pnas.1604560113](https://www.doi.org/10.1073/pnas.1604560113)


How to get help
---------------
A preferred way report any problems or ask questions about Flye is the 
[issue tracker](https://github.com/fenderglass/Flye/issues). 
Before posting an issue/question, consider to look through the FAQ
and existing issues (opened and closed) - it is possble that your question
has already been answered.

If you reporting a problem, please include the `flye.log` file and provide some 
details about your dataset (if possible).

In case you prefer personal communication, please contact Mikhail at fenderglass@gmail.com.
## Table of Contents

- [Introduction & Installation](#intro)
- [Mapping Genomic Reads](#map-reads)
  * [Mapping long reads](#map-pb)
  * [Mapping Illumina paired-end reads](#map-sr)
  * [Evaluating mapping accuracy with simulated reads (for developers)](#mapeval)
- [Mapping Long RNA-seq Reads](#map-rna)
  * [Mapping Nanopore 2D cDNA reads](#map-ont-cdna-2d)
  * [Mapping Nanopore direct-RNA reads](#map-direct-rna)
  * [Mapping PacBio Iso-seq reads](#map-iso-seq)
- [Full-Genome Alignment](#genome-aln)
  * [Intra-species assembly alignment](#asm-to-ref)
  * [Cross-species full-genome alignment](#x-species)
  * [Eyeballing alignment](#view-aln)
  * [Calling variants from assembly-to-reference alignment](#asm-var)
  * [Constructing self-homology map](#hom-map)
  * [Lift Over (for developers)](#liftover)
- [Read Overlap](#read-overlap)
  * [Long-read overlap](#long-read-overlap)
  * [Evaluating overlap sensitivity (for developers)](#ov-eval)

## <a name="intro"></a>Introduction & Installation

This cookbook walks you through a variety of applications of minimap2 and its
companion script `paftools.js`. All data here are freely available from the
minimap2 release page at version tag [v2.10][v2.10]. Some examples only work
with v2.10 or later.

To acquire the data used in this cookbook and to install minimap2 and paftools,
please follow the command lines below:
```sh
# install minimap2 executables
curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar jxf -
cp minimap2-2.17_x64-linux/{minimap2,k8,paftools.js} .  # copy executables
export PATH="$PATH:"`pwd`                               # put the current directory on PATH
# download example datasets
curl -L https://github.com/lh3/minimap2/releases/download/v2.10/cookbook-data.tgz | tar zxf -
```

## <a name="map-reads"></a>Mapping Genomic Reads

### <a name="map-pb"></a>Mapping long reads
```sh
minimap2 -ax map-pb -t4 ecoli_ref.fa ecoli_p6_25x_canu.fa > mapped.sam
```
Alternatively, you can create a minimap2 index first and then map:
```sh
minimap2 -x map-pb -d ecoli-pb.mmi ecoli_ref.fa                      # create an index
minimap2 -ax map-pb ecoli-pb.mmi ecoli_p6_25x_canu.fa > mapped.sam
```
This will save you a couple of minutes when you map against the human genome.
**HOWEVER**, key algorithm parameters such as the k-mer length and window
size can't be changed after indexing. Minimap2 will give you a warning if
parameters used in a pre-built index doesn't match parameters on the command
line. **Please always make sure you are using an intended pre-built index.**

### <a name="map-sr"></a>Mapping Illumina paired-end reads:
```sh
minimap2 -ax sr -t4 ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq > mapped-sr.sam
```

### <a name="mapeval"></a>Evaluating mapping accuracy with simulated reads (for developers)
```sh
minimap2 -ax sr ecoli_ref.fa ecoli_mason_1.fq ecoli_mason_2.fq | paftools.js mapeval -
```
The output is:
```
Q       60      19712   0       0.000000000     19712
Q       0       282     219     0.010953286     19994
U       6
```
where a `U`-line gives the number of unmapped reads (for SAM input only); a
`Q`-line gives:

1. Mapping quality (mapQ) threshold
2. Number of mapped reads between this threshold and the previous mapQ threshold.
3. Number of wrong mappings in the same mapQ interval
4. Accumulative mapping error rate
5. Accumulative number of mappings

For `paftools.js mapeval` to work, you need to encode the true read positions
in read names in the right format. For [PBSIM][pbsim] and [mason2][mason2], we
provide scripts to generate the right format. Simulated reads in this cookbook
were created with the following command lines:
```sh
# in PBSIM source code directory:
src/pbsim ../ecoli_ref.fa --depth 1 --sample-fastq sample/sample.fastq
paftools.js pbsim2fq ../ecoli_ref.fa.fai sd_0001.maf > ../ecoli_pbsim.fa

# mason2 simulation
mason_simulator --illumina-prob-mismatch-scale 2.5 -ir ecoli_ref.fa -n 10000 -o tmp-l.fq -or tmp-r.fq -oa tmp.sam
paftools.js mason2fq tmp.sam | seqtk seq -1 > ecoli_mason_1.fq
paftools.js mason2fq tmp.sam | seqtk seq -2 > ecoli_mason_2.fq
```



## <a name="map-rna"></a>Mapping Long RNA-seq Reads

### <a name="map-ont-cdna-2d"></a>Mapping Nanopore 2D cDNA reads
```sh
minimap2 -ax splice SIRV_E2.fa SIRV_ont-cdna.fa > aln.sam
```
You can compare the alignment to the true annotations with:
```sh
paftools.js junceval SIRV_E2C.gtf aln.sam
```
It gives the percentage of introns found in the annotation. For SIRV data, it
is possible to achieve higher junction accuracy with
```sh
minimap2 -ax splice --splice-flank=no SIRV_E2.fa SIRV_ont-cdna.fa | paftools.js junceval SIRV_E2C.gtf
```
This is because minimap2 models one additional evolutionarily conserved base
around a canonical junction, but SIRV doesn't honor this signal. Option
`--splice-flank=no` asks minimap2 no to model this additional base.

In the output a tag `ts:A:+` indicates that the read strand is the same as the
transcript strand; `ts:A:-` indicates the read strand is opposite to the
transcript strand. This tag is inferred from the GT-AG signal and is thus only
available to spliced reads.

### <a name="map-direct-rna"></a>Mapping Nanopore direct-RNA reads
```sh
minimap2 -ax splice -k14 -uf SIRV_E2.fa SIRV_ont-drna.fa > aln.sam
```
Direct-RNA reads are noisier, so we use a shorter k-mer for improved
sensitivity. Here, option `-uf` forces minimap2 to map reads to the forward
transcript strand only because direct-RNA reads are stranded. Again, applying
`--splice-flank=no` helps junction accuracy for SIRV data.

### <a name="map-iso-seq"></a>Mapping PacBio Iso-seq reads
```sh
minimap2 -ax splice -uf -C5 SIRV_E2.fa SIRV_iso-seq.fq > aln.sam
```
Option `-C5` reduces the penalty on non-canonical splicing sites. It helps
to align such sites correctly for data with low error rate such as Iso-seq
reads and traditional cDNAs. On this example, minimap2 makes one junction
error. Applying `--splice-flank=no` fixes this alignment error.

Note that the command line above is optimized for the final Iso-seq reads.
PacBio's Iso-seq pipeline produces intermediate sequences at varying quality.
For example, some intermediate reads are not stranded. For these reads, option
`-uf` will lead to more errors. Please revise the minimap2 command line
accordingly.



## <a name="genome-aln"></a>Full-Genome Alignment

### <a name="asm-to-ref"></a>Intra-species assembly alignment
```sh
# option "--cs" is recommended as paftools.js may need it
minimap2 -cx asm5 --cs ecoli_ref.fa ecoli_canu.fa > ecoli_canu.paf
```
Here `ecoli_canu.fa` is the Canu assembly of `ecoli_p6_25x_canu.fa`. This
command line outputs alignments in the [PAF format][paf]. Use `-a` instead of
`-c` to get output in the SAM format.

### <a name="x-species"></a>Cross-species full-genome alignment
```sh
minimap2 -cx asm20 --cs ecoli_ref.fa ecoli_O104:H4.fa > ecoli_O104:H4.paf
sort -k6,6 -k8,8n ecoli_O104:H4.paf | paftools.js call -f ecoli_ref.fa -L10000 -l1000 - > out.vcf
```
Minimap2 has three presets for full-genome alignment: "asm5" for sequence
divergence below 1%, "asm10" for divergence around a couple of percent and
"asm20" for divergence not more than 10%. In theory, with the right setting,
minimap2 should work for sequence pairs with sequence divergence up to ~15%,
but this has not been carefully evaluated.

### <a name="view-aln"></a>Eyeballing alignment
```sh
# option "--cs" required; minimap2-r741 or higher required for the "asm20" preset
minimap2 -cx asm20 --cs ecoli_ref.fa ecoli_O104:H4.fa | paftools.js view - | less -S
```
This prints the alignment in a BLAST-like format.

### <a name="asm-var"></a>Calling variants from assembly-to-reference alignment
```sh
# don't forget the "--cs" option; otherwise it doesn't work
minimap2 -cx asm5 --cs ecoli_ref.fa ecoli_canu.fa \
  | sort -k6,6 -k8,8n \
  | paftools.js call -f ecoli_ref.fa - > out.vcf
```
Without option `-f`, `paftools.js call` outputs in a custom format. In this
format, lines starting with `R` give the regions covered by one contig only.
This information is not available in the VCF output.

### <a name="hom-map"></a>Constructing self-homology map
```sh
minimap2 -DP -k19 -w19 -m200 ecoli_ref.fa ecoli_ref.fa > out.paf
```
Option `-D` asks minimap2 to ignore anchors from perfect self match and `-P`
outputs all chains. For large nomes, we don't recommend to perform base-level
alignment (with `-c`, `-a` or `--cs`) when `-P` is applied. This is because
base-alignment is slow and occasionally gives wrong alignments close to the
diagonal of a dotter plot. For E. coli, though, base-alignment is still fast.

### <a name="liftover"></a>Lift over (for developers)
```sh
minimap2 -cx asm5 --cs ecoli_ref.fa ecoli_canu.fa > ecoli_canu.paf
echo -e 'tig00000001\t200000\t300000' | paftools.js liftover ecoli_canu.paf -
```
This lifts over a region on query sequences to one or multiple regions on
reference sequences. Note that this paftools.js command may not be efficient
enough to lift millions of regions.



## <a name="read-overlap"></a>Read Overlap

### <a name="long-read-overlap"></a>Long read overlap
```sh
# For pacbio reads:
minimap2 -x ava-pb ecoli_p6_25x_canu.fa ecoli_p6_25x_canu.fa > overlap.paf
# For Nanopore reads (ava-ont also works with PacBio but not as good):
minimap2 -x ava-ont -r 10000 ecoli_p6_25x_canu.fa ecoli_p6_25x_canu.fa > overlap.paf
# If you have miniasm installed:
miniasm -f ecoli_p6_25x_canu.fa overlap.paf > asm.gfa
```
Here we explicitly applied `-r 10000`. We are considering to set this as the
default for the `ava-ont` mode as this seems to improve the contiguity for
nanopore read assembly (Loman, personal communication).

*Minimap2 doesn't work well with short-read overlap.*

### <a name="ov-eval"></a>Evaluating overlap sensitivity (for developers)

```sh
# read to reference mapping
minimap2 -cx map-pb ecoli_ref.fa ecoli_p6_25x_canu.fa > to-ref.paf
# evaluate overlap sensitivity
sort -k6,6 -k8,8n to-ref.paf | paftools.js ov-eval - overlap.paf
```
You can see that for PacBio reads, minimap2 achieves higher overlap sensitivity
with `-x ava-pb` (99% vs 93% with `-x ava-ont`).



[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[v2.10]: https://github.com/lh3/minimap2/releases/tag/v2.10
[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/minimap2/total.svg?style=social&logo=github&label=Download)](https://github.com/lh3/minimap2/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/minimap2.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/minimap2)
[![PyPI](https://img.shields.io/pypi/v/mappy.svg?style=flat)](https://pypi.python.org/pypi/mappy)
[![Build Status](https://travis-ci.org/lh3/minimap2.svg?branch=master)](https://travis-ci.org/lh3/minimap2)
## <a name="started"></a>Getting Started
```sh
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# long sequences against a reference genome
./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam
# create an index first and then map
./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa
./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam
# use presets (no test data)
./minimap2 -ax map-pb ref.fa pacbio.fq.gz > aln.sam       # PacBio genomic reads
./minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads
./minimap2 -ax asm20 ref.fa pacbio-ccs.fq.gz > aln.sam    # PacBio CCS genomic reads
./minimap2 -ax sr ref.fa read1.fa read2.fa > aln.sam      # short genomic paired-end reads
./minimap2 -ax splice ref.fa rna-reads.fa > aln.sam       # spliced long reads (strand unknown)
./minimap2 -ax splice -uf -k14 ref.fa reads.fa > aln.sam  # noisy Nanopore Direct RNA-seq
./minimap2 -ax splice:hq -uf ref.fa query.fa > aln.sam    # Final PacBio Iso-seq or traditional cDNA
./minimap2 -cx asm5 asm1.fa asm2.fa > aln.paf             # intra-species asm-to-asm alignment
./minimap2 -x ava-pb reads.fa reads.fa > overlaps.paf     # PacBio read overlap
./minimap2 -x ava-ont reads.fa reads.fa > overlaps.paf    # Nanopore read overlap
# man page for detailed command line options
man ./minimap2.1
```
## Table of Contents

- [Getting Started](#started)
- [Users' Guide](#uguide)
  - [Installation](#install)
  - [General usage](#general)
  - [Use cases](#cases)
    - [Map long noisy genomic reads](#map-long-genomic)
    - [Map long mRNA/cDNA reads](#map-long-splice)
    - [Find overlaps between long reads](#long-overlap)
    - [Map short accurate genomic reads](#short-genomic)
    - [Full genome/assembly alignment](#full-genome)
  - [Advanced features](#advanced)
    - [Working with >65535 CIGAR operations](#long-cigar)
    - [The cs optional tag](#cs)
    - [Working with the PAF format](#paftools)
  - [Algorithm overview](#algo)
  - [Getting help](#help)
  - [Citing minimap2](#cite)
- [Developers' Guide](#dguide)
- [Limitations](#limit)

## <a name="uguide"></a>Users' Guide

Minimap2 is a versatile sequence alignment program that aligns DNA or mRNA
sequences against a large reference database. Typical use cases include: (1)
mapping PacBio or Oxford Nanopore genomic reads to the human genome; (2)
finding overlaps between long reads with error rate up to ~15%; (3)
splice-aware alignment of PacBio Iso-Seq or Nanopore cDNA or Direct RNA reads
against a reference genome; (4) aligning Illumina single- or paired-end reads;
(5) assembly-to-assembly alignment; (6) full-genome alignment between two
closely related species with divergence below ~15%.

For ~10kb noisy reads sequences, minimap2 is tens of times faster than
mainstream long-read mappers such as BLASR, BWA-MEM, NGMLR and GMAP. It is more
accurate on simulated long reads and produces biologically meaningful alignment
ready for downstream analyses. For >100bp Illumina short reads, minimap2 is
three times as fast as BWA-MEM and Bowtie2, and as accurate on simulated data.
Detailed evaluations are available from the [minimap2 paper][doi] or the
[preprint][preprint].

### <a name="install"></a>Installation

Minimap2 is optimized for x86-64 CPUs. You can acquire precompiled binaries from
the [release page][release] with:
```sh
curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar -jxvf -
./minimap2-2.17_x64-linux/minimap2
```
If you want to compile from the source, you need to have a C compiler, GNU make
and zlib development files installed. Then type `make` in the source code
directory to compile. If you see compilation errors, try `make sse2only=1`
to disable SSE4 code, which will make minimap2 slightly slower.

Minimap2 also works with ARM CPUs supporting the NEON instruction sets. To
compile for 32 bit ARM architectures (such as ARMv7), use `make arm_neon=1`. To compile for for 64 bit ARM architectures (such as ARMv8), use `make arm_neon=1 aarch64=1`.

### <a name="general"></a>General usage

Without any options, minimap2 takes a reference database and a query sequence
file as input and produce approximate mapping, without base-level alignment
(i.e. no CIGAR), in the [PAF format][paf]:
```sh
minimap2 ref.fa query.fq > approx-mapping.paf
```
You can ask minimap2 to generate CIGAR at the `cg` tag of PAF with:
```sh
minimap2 -c ref.fa query.fq > alignment.paf
```
or to output alignments in the [SAM format][sam]:
```sh
minimap2 -a ref.fa query.fq > alignment.sam
```
Minimap2 seamlessly works with gzip'd FASTA and FASTQ formats as input. You
don't need to convert between FASTA and FASTQ or decompress gzip'd files first.

For the human reference genome, minimap2 takes a few minutes to generate a
minimizer index for the reference before mapping. To reduce indexing time, you
can optionally save the index with option **-d** and replace the reference
sequence file with the index file on the minimap2 command line:
```sh
minimap2 -d ref.mmi ref.fa                     # indexing
minimap2 -a ref.mmi reads.fq > alignment.sam   # alignment
```
***Importantly***, it should be noted that once you build the index, indexing
parameters such as **-k**, **-w**, **-H** and **-I** can't be changed during
mapping. If you are running minimap2 for different data types, you will
probably need to keep multiple indexes generated with different parameters.
This makes minimap2 different from BWA which always uses the same index
regardless of query data types.

### <a name="cases"></a>Use cases

Minimap2 uses the same base algorithm for all applications. However, due to the
different data types it supports (e.g. short vs long reads; DNA vs mRNA reads),
minimap2 needs to be tuned for optimal performance and accuracy. It is usually
recommended to choose a preset with option **-x**, which sets multiple
parameters at the same time. The default setting is the same as `map-ont`.

#### <a name="map-long-genomic"></a>Map long noisy genomic reads

```sh
minimap2 -ax map-pb  ref.fa pacbio-reads.fq > aln.sam   # for PacBio subreads
minimap2 -ax map-ont ref.fa ont-reads.fq > aln.sam      # for Oxford Nanopore reads
```
The difference between `map-pb` and `map-ont` is that `map-pb` uses
homopolymer-compressed (HPC) minimizers as seeds, while `map-ont` uses ordinary
minimizers as seeds. Emperical evaluation suggests HPC minimizers improve
performance and sensitivity when aligning PacBio reads, but hurt when aligning
Nanopore reads.

#### <a name="map-long-splice"></a>Map long mRNA/cDNA reads

```sh
minimap2 -ax splice:hq -uf ref.fa iso-seq.fq > aln.sam       # PacBio Iso-seq/traditional cDNA
minimap2 -ax splice ref.fa nanopore-cdna.fa > aln.sam        # Nanopore 2D cDNA-seq
minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam  # Nanopore Direct RNA-seq
minimap2 -ax splice --splice-flank=no SIRV.fa SIRV-seq.fa    # mapping against SIRV control
```
There are different long-read RNA-seq technologies, including tranditional
full-length cDNA, EST, PacBio Iso-seq, Nanopore 2D cDNA-seq and Direct RNA-seq.
They produce data of varying quality and properties. By default, `-x splice`
assumes the read orientation relative to the transcript strand is unknown. It
tries two rounds of alignment to infer the orientation and write the strand to
the `ts` SAM/PAF tag if possible. For Iso-seq, Direct RNA-seq and tranditional
full-length cDNAs, it would be desired to apply `-u f` to force minimap2 to
consider the forward transcript strand only. This speeds up alignment with
slight improvement to accuracy. For noisy Nanopore Direct RNA-seq reads, it is
recommended to use a smaller k-mer size for increased sensitivity to the first
or the last exons.

Minimap2 rates an alignment by the score of the max-scoring sub-segment,
*excluding* introns, and marks the best alignment as primary in SAM. When a
spliced gene also has unspliced pseudogenes, minimap2 does not intentionally
prefer spliced alignment, though in practice it more often marks the spliced
alignment as the primary. By default, minimap2 outputs up to five secondary
alignments (i.e. likely pseudogenes in the context of RNA-seq mapping). This
can be tuned with option **-N**.

For long RNA-seq reads, minimap2 may produce chimeric alignments potentially
caused by gene fusions/structural variations or by an intron longer than the
max intron length **-G** (200k by default). For now, it is not recommended to
apply an excessively large **-G** as this slows down minimap2 and sometimes
leads to false alignments.

It is worth noting that by default `-x splice` prefers GT[A/G]..[C/T]AG
over GT[C/T]..[A/G]AG, and then over other splicing signals. Considering
one additional base improves the junction accuracy for noisy reads, but
reduces the accuracy when aligning against the widely used SIRV control data.
This is because SIRV does not honor the evolutionarily conservative splicing
signal. If you are studying SIRV, you may apply `--splice-flank=no` to let
minimap2 only model GT..AG, ignoring the additional base.

#### <a name="long-overlap"></a>Find overlaps between long reads

```sh
minimap2 -x ava-pb  reads.fq reads.fq > ovlp.paf    # PacBio read overlap
minimap2 -x ava-ont reads.fq reads.fq > ovlp.paf    # Oxford Nanopore read overlap
```
Similarly, `ava-pb` uses HPC minimizers while `ava-ont` uses ordinary
minimizers. It is usually not recommended to perform base-level alignment in
the overlapping mode because it is slow and may produce false positive
overlaps. However, if performance is not a concern, you may try to add `-a` or
`-c` anyway.

#### <a name="short-genomic"></a>Map short accurate genomic reads

```sh
minimap2 -ax sr ref.fa reads-se.fq > aln.sam           # single-end alignment
minimap2 -ax sr ref.fa read1.fq read2.fq > aln.sam     # paired-end alignment
minimap2 -ax sr ref.fa reads-interleaved.fq > aln.sam  # paired-end alignment
```
When two read files are specified, minimap2 reads from each file in turn and
merge them into an interleaved stream internally. Two reads are considered to
be paired if they are adjacent in the input stream and have the same name (with
the `/[0-9]` suffix trimmed if present). Single- and paired-end reads can be
mixed.

Minimap2 does not work well with short spliced reads. There are many capable
RNA-seq mappers for short reads.

#### <a name="full-genome"></a>Full genome/assembly alignment

```sh
minimap2 -ax asm5 ref.fa asm.fa > aln.sam       # assembly to assembly/ref alignment
```
For cross-species full-genome alignment, the scoring system needs to be tuned
according to the sequence divergence.

### <a name="advanced"></a>Advanced features

#### <a name="long-cigar"></a>Working with >65535 CIGAR operations

Due to a design flaw, BAM does not work with CIGAR strings with >65535
operations (SAM and CRAM work). However, for ultra-long nanopore reads minimap2
may align ~1% of read bases with long CIGARs beyond the capability of BAM. If
you convert such SAM/CRAM to BAM, Picard and recent samtools will throw an
error and abort. Older samtools and other tools may create corrupted BAM.

To avoid this issue, you can add option `-L` at the minimap2 command line.
This option moves a long CIGAR to the `CG` tag and leaves a fully clipped CIGAR
at the SAM CIGAR column. Current tools that don't read CIGAR (e.g. merging and
sorting) still work with such BAM records; tools that read CIGAR will
effectively ignore these records. It has been decided that future tools will
will seamlessly recognize long-cigar records generated by option `-L`.

**TL;DR**: if you work with ultra-long reads and use tools that only process
BAM files, please add option `-L`.

#### <a name="cs"></a>The cs optional tag

The `cs` SAM/PAF tag encodes bases at mismatches and INDELs. It matches regular
expression `/(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+)+/`. Like CIGAR, `cs`
consists of series of operations.  Each leading character specifies the
operation; the following sequence is the one involved in the operation.

The `cs` tag is enabled by command line option `--cs`. The following alignment,
for example:
```txt
CGATCGATAAATAGAGTAG---GAATAGCA
||||||   ||||||||||   |||| |||
CGATCG---AATAGAGTAGGTCGAATtGCA
```
is represented as `:6-ata:10+gtc:4*at:3`, where `:[0-9]+` represents an
identical block, `-ata` represents a deltion, `+gtc` an insertion and `*at`
indicates reference base `a` is substituted with a query base `t`. It is
similar to the `MD` SAM tag but is standalone and easier to parse.

If `--cs=long` is used, the `cs` string also contains identical sequences in
the alignment. The above example will become
`=CGATCG-ata=AATAGAGTAG+gtc=GAAT*at=GCA`. The long form of `cs` encodes both
reference and query sequences in one string. The `cs` tag also encodes intron
positions and splicing signals (see the [minimap2 manpage][manpage-cs] for
details).

#### <a name="paftools"></a>Working with the PAF format

Minimap2 also comes with a (java)script [paftools.js](misc/paftools.js) that
processes alignments in the PAF format. It calls variants from
assembly-to-reference alignment, lifts over BED files based on alignment,
converts between formats and provides utilities for various evaluations. For
details, please see [misc/README.md](misc/README.md).

### <a name="algo"></a>Algorithm overview

In the following, minimap2 command line options have a dash ahead and are
highlighted in bold. The description may help to tune minimap2 parameters.

1. Read **-I** [=*4G*] reference bases, extract (**-k**,**-w**)-minimizers and
   index them in a hash table.

2. Read **-K** [=*200M*] query bases. For each query sequence, do step 3
   through 7:

3. For each (**-k**,**-w**)-minimizer on the query, check against the reference
   index. If a reference minimizer is not among the top **-f** [=*2e-4*] most
   frequent, collect its the occurrences in the reference, which are called
   *seeds*.

4. Sort seeds by position in the reference. Chain them with dynamic
   programming. Each chain represents a potential mapping. For read
   overlapping, report all chains and then go to step 8. For reference mapping,
   do step 5 through 7:

5. Let *P* be the set of primary mappings, which is an empty set initially. For
   each chain from the best to the worst according to their chaining scores: if
   on the query, the chain overlaps with a chain in *P* by **--mask-level**
   [=*0.5*] or higher fraction of the shorter chain, mark the chain as
   *secondary* to the chain in *P*; otherwise, add the chain to *P*.

6. Retain all primary mappings. Also retain up to **-N** [=*5*] top secondary
   mappings if their chaining scores are higher than **-p** [=*0.8*] of their
   corresponding primary mappings.

7. If alignment is requested, filter out an internal seed if it potentially
   leads to both a long insertion and a long deletion. Extend from the
   left-most seed. Perform global alignments between internal seeds.  Split the
   chain if the accumulative score along the global alignment drops by **-z**
   [=*400*], disregarding long gaps. Extend from the right-most seed.  Output
   chains and their alignments.

8. If there are more query sequences in the input, go to step 2 until no more
   queries are left.

9. If there are more reference sequences, reopen the query file from the start
   and go to step 1; otherwise stop.

### <a name="help"></a>Getting help

Manpage [minimap2.1][manpage] provides detailed description of minimap2
command line options and optional tags. If you encounter bugs or have further
questions or requests, you can raise an issue at the [issue page][issue].
There is not a specific mailing list for the time being.

### <a name="cite"></a>Citing minimap2

If you use minimap2 in your work, please cite:

> Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
> *Bioinformatics*, **34**:3094-3100. [doi:10.1093/bioinformatics/bty191][doi]

## <a name="dguide"></a>Developers' Guide

Minimap2 is not only a command line tool, but also a programming library.
It provides C APIs to build/load index and to align sequences against the
index. File [example.c](example.c) demonstrates typical uses of C APIs. Header
file [minimap.h](minimap.h) gives more detailed API documentation. Minimap2
aims to keep APIs in this header stable. File [mmpriv.h](mmpriv.h) contains
additional private APIs which may be subjected to changes frequently.

This repository also provides Python bindings to a subset of C APIs. File
[python/README.rst](python/README.rst) gives the full documentation;
[python/minimap2.py](python/minimap2.py) shows an example. This Python
extension, mappy, is also [available from PyPI][mappypypi] via `pip install
mappy` or [from BioConda][mappyconda] via `conda install -c bioconda mappy`.

## <a name="limit"></a>Limitations

* Minimap2 may produce suboptimal alignments through long low-complexity
  regions where seed positions may be suboptimal. This should not be a big
  concern because even the optimal alignment may be wrong in such regions.

* Minimap2 requires SSE2 instructions on x86 CPUs or NEON on ARM CPUs. It is
  possible to add non-SIMD support, but it would make minimap2 slower by
  several times.

* Minimap2 does not work with a single query or database sequence ~2
  billion bases or longer (2,147,483,647 to be exact). The total length of all
  sequences can well exceed this threshold.

* Minimap2 often misses small exons.



[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[sam]: https://samtools.github.io/hts-specs/SAMv1.pdf
[minimap]: https://github.com/lh3/minimap
[smartdenovo]: https://github.com/ruanjue/smartdenovo
[longislnd]: https://www.ncbi.nlm.nih.gov/pubmed/27667791
[gaba]: https://github.com/ocxtal/libgaba
[ksw2]: https://github.com/lh3/ksw2
[preprint]: https://arxiv.org/abs/1708.01492
[release]: https://github.com/lh3/minimap2/releases
[mappypypi]: https://pypi.python.org/pypi/mappy
[mappyconda]: https://anaconda.org/bioconda/mappy
[issue]: https://github.com/lh3/minimap2/issues
[k8]: https://github.com/attractivechaos/k8
[manpage]: https://lh3.github.io/minimap2/minimap2.html
[manpage-cs]: https://lh3.github.io/minimap2/minimap2.html#10
[doi]: https://doi.org/10.1093/bioinformatics/bty191
Release 2.17-r941 (4 May 2019)
------------------------------

Changes since the last release:

 * Fixed flawed CIGARs like `5I6D7I` (#392).

 * Bugfix: TLEN should be 0 when either end is unmapped (#373 and #365).

 * Bugfix: mappy is unable to write index (#372).

 * Added option `--junc-bed` to load known gene annotations in the BED12
   format. Minimap2 prefers annotated junctions over novel junctions (#197 and
   #348). GTF can be converted to BED12 with `paftools.js gff2bed`.

 * Added option `--sam-hit-only` to suppress unmapped hits in SAM (#377).

 * Added preset `splice:hq` for high-quality CCS or mRNA sequences. It applies
   better scoring and improves the sensitivity to small exons. This preset may
   introduce false small introns, but the overall accuracy should be higher.

This version produces nearly identical alignments to v2.16, except for CIGARs
affected by the bug mentioned above.

(2.17: 5 May 2019, r941)



Release 2.16-r922 (28 February 2019)
------------------------------------

This release is 50% faster for mapping ultra-long nanopore reads at comparable
accuracy. For short-read mapping, long-read overlapping and ordinary long-read
mapping, the performance and accuracy remain similar. This speedup is achieved
with a new heuristic to limit the number of chaining iterations (#324). Users
can disable the heuristic by increasing a new option `--max-chain-iter` to a
huge number.

Other changes to minimap2:

 * Implemented option `--paf-no-hit` to output unmapped query sequences in PAF.
   The strand and reference name columns are both `*` at an unmapped line. The
   hidden option is available in earlier minimap2 but had a different 2-column
   output format instead of PAF.

 * Fixed a bug that leads to wrongly calculated `de` tags when ambiguous bases
   are involved (#309). This bug only affects v2.15.

 * Fixed a bug when parsing command-line option `--splice` (#344). This bug was
   introduced in v2.13.

 * Fixed two division-by-zero cases (#326). They don't affect final alignments
   because the results of the divisions are not used in both case.

 * Added an option `-o` to output alignments to a specified file. It is still
   recommended to use UNIX pipes for on-the-fly conversion or compression.

 * Output a new `rl` tag to give the length of query regions harboring
   repetitive seeds.

Changes to paftool.js:

 * Added a new option to convert the MD tag to the long form of the cs tag.

Changes to mappy:

 * Added the `mappy.Aligner.seq_names` method to return sequence names (#312).

For NA12878 ultra-long reads, this release changes the alignments of <0.1% of
reads in comparison to v2.15. All these reads have highly fragmented alignments
and are likely to be problematic anyway. For shorter or well aligned reads,
this release should produce mostly identical alignments to v2.15.

(2.16: 28 February 2019, r922)



Release 2.15-r905 (10 January 2019)
-----------------------------------

Changes to minimap2:

 * Fixed a rare segmentation fault when option -H is in use (#307). This may
   happen when there are very long homopolymers towards the 5'-end of a read.

 * Fixed wrong CIGARs when option --eqx is used (#266).

 * Fixed a typo in the base encoding table (#264). This should have no
   practical effect.

 * Fixed a typo in the example code (#265).

 * Improved the C++ compatibility by removing "register" (#261). However,
   minimap2 still can't be compiled in the pedantic C++ mode (#306).

 * Output a new "de" tag for gap-compressed sequence divergence.

Changes to paftools.js:

 * Added "asmgene" to evaluate the completeness of an assembly by measuring the
   uniquely mapped single-copy genes. This command learns the idea of BUSCO.

 * Added "vcfpair" to call a phased VCF from phased whole-genome assemblies. An
   earlier version of this script is used to produce the ground truth for the
   syndip benchmark [PMID:30013044].

This release produces identical alignment coordinates and CIGARs in comparison
to v2.14. Users are advised to upgrade due to the several bug fixes.

(2.15: 10 Janurary 2019, r905)



Release 2.14-r883 (5 November 2018)
-----------------------------------

Notable changes:

 * Fixed two minor bugs caused by typos (#254 and #266).

 * Fixed a bug that made minimap2 abort when --eqx was used together with --MD
   or --cs (#257).

 * Added --cap-sw-mem to cap the size of DP matrices (#259). Base alignment may
   take a lot of memory in the splicing mode. This may lead to issues when we
   run minimap2 on a cluster with a hard memory limit. The new option avoids
   unlimited memory usage at the cost of missing a few long introns.

 * Conforming to C99 and C11 when possible (#261).

 * Warn about malformatted FASTA or FASTQ (#252 and #255).

This release occasionally produces base alignments different from v2.13. The
overall alignment accuracy remain similar.

(2.14: 5 November 2018, r883)



Release 2.13-r850 (11 October 2018)
-----------------------------------

Changes to minimap2:

 * Fixed wrongly formatted SAM when -L is in use (#231 and #233).

 * Fixed an integer overflow in rare cases.

 * Added --hard-mask-level to fine control split alignments (#244).

 * Made --MD work with spliced alignment (#139).

 * Replaced musl's getopt with ketopt for portability.

 * Log peak memory usage on exit.

This release should produce alignments identical to v2.12 and v2.11.

(2.13: 11 October 2018, r850)



Release 2.12-r827 (6 August 2018)
---------------------------------

Changes to minimap2:

 * Added option --split-prefix to write proper alignments (correct mapping
   quality and clustered query sequences) given a multi-part index (#141 and
   #189; mostly by @hasindu2008).

 * Fixed a memory leak when option -y is in use.

Changes to mappy:

 * Support the MD/cs tag (#183 and #203).

 * Allow mappy to index a single sequence, to add extra flags and to change the
   scoring system.

Minimap2 should produce alignments identical to v2.11.

(2.12: 6 August 2018, r827)



Release 2.11-r797 (20 June 2018)
--------------------------------

Changes to minimap2:

 * Improved alignment accuracy in low-complexity regions for SV calling. Thank
   @armintoepfer for multiple offline examples.

 * Added option --eqx to encode sequence match/mismatch with the =/X CIGAR
   operators (#156, #157 and #175).

 * When compiled with VC++, minimap2 generated wrong alignments due to a
   comparison between a signed integer and an unsigned integer (#184). Also
   fixed warnings reported by "clang -Wextra".

 * Fixed incorrect anchor filtering due to a missing 64- to 32-bit cast.

 * Fixed incorrect mapping quality for inversions (#148).

 * Fixed incorrect alignment involving ambiguous bases (#155).

 * Fixed incorrect presets: option `-r 2000` is intended to be used with
   ava-ont, not ava-pb. The bug was introduced in 2.10.

 * Fixed a bug when --for-only/--rev-only is used together with --sr or
   --heap-sort=yes (#166).

 * Fixed option -Y that was not working in the previous releases.

 * Added option --lj-min-ratio to fine control the alignment of long gaps
   found by the "long-join" heuristic (#128).

 * Exposed `mm_idx_is_idx`, `mm_idx_load` and `mm_idx_dump` C APIs (#177).
   Also fixed a bug when indexing without reference names (this feature is not
   exposed to the command line).

Changes to mappy:

 * Added `__version__` (#165).

 * Exposed the maximum fragment length parameter to mappy (#174).

Changes to paftools:

 * Don't crash when there is no "cg" tag (#153).

 * Fixed wrong coverage report by "paftools.js call" (#145).

This version may produce slightly different base-level alignment. The overall
alignment statistics should remain similar.

(2.11: 20 June 2018, r797)



Release 2.10-r761 (27 March 2018)
---------------------------------

Changes to minimap2:

 * Optionally output the MD tag for compatibility with existing tools (#63,
   #118 and #137).

 * Use SSE compiler flags more precisely to prevent compiling errors on certain
   machines (#127).

 * Added option --min-occ-floor to set a minimum occurrence threshold. Presets
   intended for assembly-to-reference alignment set this option to 100. This
   option alleviates issues with regions having high copy numbers (#107).

 * Exit with non-zero code on file writing errors (e.g. disk full; #103 and
   #132).

 * Added option -y to copy FASTA/FASTQ comments in query sequences to the
   output (#136).

 * Added the asm20 preset for alignments between genomes at 5-10% sequence
   divergence.

 * Changed the band-width in the ava-ont preset from 500 to 2000. Oxford
   Nanopore reads may contain long deletion sequencing errors that break
   chaining.

Changes to mappy, the Python binding:

 * Fixed a typo in Align.seq() (#126).

Changes to paftools.js, the companion script:

 * Command sam2paf now converts the MD tag to cs.

 * Support VCF output for assembly-to-reference variant calling (#109).

This version should produce identical alignment for read overlapping, RNA-seq
read mapping, and genomic read mapping. We have also added a cook book to show
the variety uses of minimap2 on real datasets. Please see cookbook.md in the
minimap2 source code directory.

(2.10: 27 March 2017, r761)



Release 2.9-r720 (23 February 2018)
-----------------------------------

This release fixed multiple minor bugs.

* Fixed two bugs that lead to incorrect inversion alignment. Also improved the
  sensitivity to small inversions by using double Z-drop cutoff (#112).

* Fixed an issue that may cause the end of a query sequence unmapped (#104).

* Added a mappy API to retrieve sequences from the index (#126) and to reverse
  complement DNA sequences. Fixed a bug where the `best_n` parameter did not
  work (#117).

* Avoided segmentation fault given incorrect FASTQ input (#111).

* Combined all auxiliary javascripts to paftools.js. Fixed several bugs in
  these scripts at the same time.

(2.9: 24 February 2018, r720)



Release 2.8-r672 (1 February 2018)
----------------------------------

Notable changes in this release include:

 * Speed up short-read alignment by ~10%. The overall mapping accuracy stays
   the same, but the output alignments are not always identical to v2.7 due to
   unstable sorting employed during chaining. Long-read alignment is not
   affected by this change as the speedup is short-read specific.

 * Mappy now supports paired-end short-read alignment (#87). Please see
   python/README.rst for details.

 * Added option --for-only and --rev-only to perform alignment against the
   forward or the reverse strand of the reference genome only (#91).

 * Alleviated the issue with undesired diagonal alignment in the self mapping
   mode (#10). Even if the output is not ideal, it should not interfere with
   other alignments. Fully resolving the issue is intricate and may require
   additional heuristic thresholds.

 * Enhanced error checking against incorrect input (#92 and #96).

For long query sequences, minimap2 should output identical alignments to v2.7.

(2.8: 1 February 2018, r672)



Release 2.7-r654 (9 January 2018)
---------------------------------

This release fixed a bug in the splice mode and added a few minor features:

 * Fixed a bug that occasionally takes an intron as a long deletion in the
   splice mode. This was caused by wrong backtracking at the last CIGAR
   operator. The current fix eliminates the error, but it is not optimal in
   that it often produces a wrong junction when the last operator is an intron.
   A future version of minimap2 may improve upon this.

 * Support high-end ARM CPUs that implement the NEON instruction set (#81).
   This enables minimap2 to work on Raspberry Pi 3 and Odroid XU4.

 * Added a C API to construct a minimizer index from a set of C strings (#80).

 * Check scoring specified on the command line (#79). Due to the 8-bit limit,
   excessively large score penalties fail minimap2.

For genomic sequences, minimap2 should give identical alignments to v2.6.

(2.7: 9 January 2018, r654)



Release 2.6-r623 (12 December 2017)
-----------------------------------

This release adds several features and fixes two minor bugs:

 * Optionally build an index without sequences. This helps to reduce the
   peak memory for read overlapping and is automatically applied when
   base-level alignment is not requested.

 * Approximately estimate per-base sequence divergence (i.e. 1-identity)
   without performing base-level alignment, using a MashMap-like method. The
   estimate is written to a new dv:f tag.

 * Reduced the number of tiny terminal exons in RNA-seq alignment. The current
   setting is conservative. Increase --end-seed-pen to drop more such exons.

 * Reduced the peak memory when aligning long query sequences.

 * Fixed a bug that is caused by HPC minimizers longer than 256bp. This should
   have no effect in practice, but it is recommended to rebuild HPC indices if
   possible.

 * Fixed a bug when identifying identical hits (#71). This should only affect
   artifactual reference consisting of near identical sequences.

For genomic sequences, minimap2 should give nearly identical alignments to
v2.5, except the new dv:f tag.

(2.6: 12 December 2017, r623)



Release 2.5-r572 (11 November 2017)
-----------------------------------

This release fixes several bugs and brings a couple of minor improvements:

 * Fixed a severe bug that leads to incorrect mapping coordinates in rare
   corner cases.

 * Fixed underestimated mapping quality for chimeric alignments when the whole
   query sequence contain many repetitive minimizers, and for chimeric
   alignments caused by Z-drop.

 * Fixed two bugs in Python binding: incorrect strand field (#57) and incorrect
   sequence names for Python3 (#55).

 * Improved mapping accuracy for highly overlapping paired ends.

 * Added option -Y to use soft clipping for supplementary alignments (#56).

(2.5: 11 November 2017, r572)



Release 2.4-r555 (6 November 2017)
----------------------------------

As is planned, this release focuses on fine tuning the base algorithm. Notable
changes include

 * Changed the mapping quality scale to match the scale of BWA-MEM. This makes
   minimap2 and BWA-MEM achieve similar sensitivity-specificity balance on real
   short-read data.

 * Improved the accuracy of splice alignment by modeling one additional base
   close to the GT-AG signal. This model is used by default with `-x splice`.
   For SIRV control data, however, it is recommended to add `--splice-flank=no`
   to disable this feature as the SIRV splice signals are slightly different.

 * Tuned the parameters for Nanopore Direct RNA reads. The recommended command
   line is `-axsplice -k14 -uf` (#46).

 * Fixed a segmentation fault when aligning PacBio reads (#47 and #48). This
   bug is very rare but it affects all versions of minimap2. It is also
   recommended to re-index reference genomes created with `map-pb`. For human,
   two minimizers in an old index are wrong.

 * Changed option `-L` in sync with the final decision of hts-specs: a fake
   CIGAR takes the form of `<readLen>S<refLen>N`. Note that `-L` only enables
   future tools to recognize long CIGARs. It is not possible for older tools to
   work with such alignments in BAM (#43 and #51).

 * Fixed a tiny issue whereby minimap2 may waste 8 bytes per candidate
   alignment.

The minimap2 technical note hosted at arXiv has also been updated to reflect
recent changes.

(2.4: 6 November 2017, r555)



Release 2.3-r531 (22 October 2017)
----------------------------------

This release come with many improvements and bug fixes:

 * The **sr** preset now supports paired-end short-read alignment. Minimap2 is
   3-4 times as fast as BWA-MEM, but is slightly less accurate on simulated
   reads.

 * Meticulous improvements to assembly-to-assembly alignment (special thanks to
   Alexey Gurevich from the QUAST team): a) apply a small penalty to matches
   between ambiguous bases; b) reduce missing alignments due to spurious
   overlaps; c) introduce the short form of the `cs` tag, an improvement to the
   SAM MD tag.

 * Make sure gaps are always left-aligned.

 * Recognize `U` bases from Oxford Nanopore Direct RNA-seq (#33).

 * Fixed slightly wrong chaining score. Fixed slightly inaccurate coordinates
   for split alignment.

 * Fixed multiple reported bugs: 1) wrong reference name for inversion
   alignment (#30); 2) redundant SQ lines when multiple query files are
   specified (#39); 3) non-functioning option `-K` (#36).

This release has implemented all the major features I planned five months ago,
with the addition of spliced long-read alignment. The next couple of releases
will focus on fine tuning of the base algorithms.

(2.3: 22 October 2017, r531)



Release 2.2-r409 (17 September 2017)
------------------------------------

This is a feature release. It improves single-end short-read alignment and
comes with Python bindings. Detailed changes include:

 * Added the **sr** preset for single-end short-read alignment. In this mode,
   minimap2 runs faster than BWA-MEM, but is slightly less accurate on
   simulated data sets. Paired-end alignment is not supported as of now.

 * Improved mapping quality estimate with more accurate identification of
   repetitive hits. This mainly helps short-read alignment.

 * Implemented **mappy**, a Python binding for minimap2, which is available
   from PyPI and can be installed with `pip install --user mappy`. Python users
   can perform read alignment without the minimap2 executable.

 * Restructured the indexing APIs and documented key minimap2 APIs in the
   header file minimap.h. Updated example.c with the new APIs. Old APIs still
   work but may become deprecated in future.

This release may output alignments different from the previous version, though
the overall alignment statistics, such as the number of aligned bases and long
gaps, remain close.

(2.2: 17 September 2017, r409)



Release 2.1.1-r341 (6 September 2017)
-------------------------------------

This is a maintenance release that is expected to output identical alignment to
v2.1. Detailed changes include:

 * Support CPU dispatch. By default, minimap2 is compiled with both SSE2 and
   SSE4 based implementation of alignment and automatically chooses the right
   one at runtime. This avoids unexpected errors on older CPUs (#21).

 * Improved Windows support as is requested by Oxford Nanopore (#19). Minimap2
   now avoids variable-length stacked arrays, eliminates alloca(), ships with
   getopt_long() and provides timing functions implemented with Windows APIs.

 * Fixed a potential segmentation fault when specifying -k/-w/-H with
   multi-part index (#23).

 * Fixed two memory leaks in example.c

(2.1.1: 6 September 2017, r341)



Release 2.1-r311 (25 August 2017)
---------------------------------

This release adds spliced alignment for long noisy RNA-seq reads. On a SMRT
Iso-Seq and a Oxford Nanopore data sets, minimap2 appears to outperform
traditional mRNA aligners. For DNA alignment, this release gives almost
identical output to v2.0. Other changes include:

 * Added option `-R` to set the read group header line in SAM.

 * Optionally output the `cs:Z` tag in PAF to encode both the query and the
   reference sequences in the alignment.

 * Fixed an issue where DP alignment uses excessive memory.

The minimap2 technical report has been updated with more details and the
evaluation of spliced alignment:

 * Li, H. (2017). Minimap2: fast pairwise alignment for long nucleotide
   sequences. [arXiv:1708.01492v2](https://arxiv.org/abs/1708.01492v2).

(2.1: 25 August 2017, r311)



Release 2.0-r275 (8 August 2017)
--------------------------------

This release is identical to version 2.0rc1, except the version number. It is
described and evaluated in the following technical report:

 * Li, H. (2017). Minimap2: fast pairwise alignment for long DNA sequences.
   [arXiv:1708.01492v1](https://arxiv.org/abs/1708.01492v1).

(2.0: 8 August 2017, r275)



Release 2.0rc1-r232 (30 July 2017)
----------------------------------

This release improves the accuracy of long-read alignment and added several
minor features.

 * Improved mapping quality estimate for short alignments containing few seed
   hits.

 * Fixed a minor bug that affects the chaining accuracy towards the ends of a
   chain. Changed the gap cost for chaining to reduce false seeding.

 * Skip potentially wrong seeding and apply dynamic programming more frequently.
   This slightly increases run time, but greatly reduces false long gaps.

 * Perform local alignment at Z-drop break point to recover potential inversion
   alignment. Output the SA tag in the SAM format. Added scripts to evaluate
   mapping accuracy for reads simulated with pbsim.

This release completes features intended for v2.0. No major features will be
added to the master branch before the final v2.0.

(2.0rc1: 30 July 2017, r232)



Release r191 (19 July 2017)
---------------------------

This is the first public release of minimap2, an aligner for long reads and
assemblies. This release has a few issues and is generally not recommended for
production uses.

(19 July 2017, r191)
## <a name="started"></a>Getting Started

```sh
# install minimap2
git clone https://github.com/lh3/minimap2
cd minimap2 && make
# install the k8 javascript shell
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` k8              # or copy it to a directory on your $PATH
# export PATH="$PATH:`pwd`:`pwd`/misc"    # run this if k8, minimap2 or paftools.js not on your $PATH
minimap2 --cs test/MT-human.fa test/MT-orang.fa | paftools.js view -     # view alignment
minimap2 -c test/MT-human.fa test/MT-orang.fa | paftools.js stat -       # basic alignment statistics
minimap2 -c --cs test/MT-human.fa test/MT-orang.fa \
  | sort -k6,6 -k8,8n | paftools.js call -L15000 -     # calling variants from asm-to-ref alignment
minimap2 -c test/MT-human.fa test/MT-orang.fa \
  | paftools.js liftover -l10000 - <(echo -e "MT_orang\t2000\t5000")     # liftOver
# no test data for the following examples
paftools.js junceval -e anno.gtf splice.sam > out.txt  # compare splice junctions to annotations
paftools.js splice2bed anno.gtf > anno.bed             # convert GTF/GFF3 to BED12
```

## Table of Contents

- [Getting Started](#started)
- [Introduction](#intro)
- [Evaluation](#eval)
  - [Evaluating mapping accuracy with simulated reads](#mapeval)
  - [Evaluating read overlap sensitivity](#oveval)
- [Calling Variants from Assemblies](#asmvar)

## <a name="intro"></a>Introduction

paftools.js is a script that processes alignments in the [PAF format][paf],
such as converting between formats, evaluating mapping accuracy, lifting over
BED files based on alignment, and calling variants from assembly-to-assembly
alignment. This script *requires* the [k8 Javascript shell][k8] to run. On
Linux or Mac, you can download the precompiled k8 binary with:

```sh
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp k8-0.2.4/k8-`uname -s` $HOME/bin/k8  # assuming $HOME/bin in your $PATH
```

It is highly recommended to copy the executable `k8` to a directory on your
`$PATH` such as `/usr/bin/env` can find it. Like python scripts, once you
install `k8`, you can launch paftools.js in one of the two ways:

```sh
path/to/paftools.js             # only if k8 is on your $PATH
k8 path/to/paftools.js
```

In a nutshell, paftools.js has the following commands:

```
Usage: paftools.js <command> [arguments]
Commands:
  view       convert PAF to BLAST-like (for eyeballing) or MAF
  splice2bed convert spliced alignment in PAF/SAM to BED12
  sam2paf    convert SAM to PAF
  delta2paf  convert MUMmer's delta to PAF
  gff2bed    convert GTF/GFF3 to BED12

  stat       collect basic mapping information in PAF/SAM
  liftover   simplistic liftOver
  call       call variants from asm-to-ref alignment with the cs tag
  bedcov     compute the number of bases covered

  mapeval    evaluate mapping accuracy using mason2/PBSIM-simulated FASTQ
  mason2fq   convert mason2-simulated SAM to FASTQ
  pbsim2fq   convert PBSIM-simulated MAF to FASTQ
  junceval   evaluate splice junction consistency with known annotations
  ov-eval    evaluate read overlap sensitivity using read-to-ref mapping
```

paftools.js seamlessly reads both plain text files and gzip'd text files.

## <a name="eval"></a>Evaluation

### <a name="mapeval"></a>Evaluating mapping accuracy with simulated reads

The **pbsim2fq** command of paftools.js converts the MAF output of [pbsim][pbsim]
to FASTQ and encodes the true mapping position in the read name in a format like
`S1_33!chr1!225258409!225267761!-`. Similarly, the **mason2fq** command
converts [mason2][mason2] simulated SAM to FASTQ.

Command **mapeval** evaluates mapped SAM/PAF. Here is example output:

```
Q       60      32478   0       0.000000000     32478
Q       22      16      1       0.000030775     32494
Q       21      43      1       0.000061468     32537
Q       19      73      1       0.000091996     32610
Q       14      66      1       0.000122414     32676
Q       10      27      3       0.000214048     32703
Q       8       14      1       0.000244521     32717
Q       7       13      2       0.000305530     32730
Q       6       46      1       0.000335611     32776
Q       3       10      1       0.000366010     32786
Q       2       20      2       0.000426751     32806
Q       1       248     94      0.003267381     33054
Q       0       31      17      0.003778147     33085
U       3
```

where each Q-line gives the quality threshold, the number of reads mapped with
mapping quality equal to or greater than the threshold, number of wrong
mappings, accumulative mapping error rate and the accumulative number of
mapped reads. The U-line, if present, gives the number of unmapped reads if
they are present in the SAM file.

Suppose the reported mapping coordinate overlap with the true coordinate like
the following:

```
truth:   --------------------
mapper:           ----------------------
         |<- l1 ->|<-- o -->|<-- l2 -->|
```

Let `r=o/(l1+o+l2)`. The reported mapping is considered correct if `r>0.1` by
default.

### <a name="oveval"></a>Evaluating read overlap sensitivity

Command **ov-eval** takes *sorted* read-to-reference alignment and read
overlaps in PAF as input, and evaluates the sensitivity. For example:

```sh
minimap2 -cx map-pb ref.fa reads.fq.gz | sort -k6,6 -k8,8n > reads-to-ref.paf
minimap2 -x ava-pb reads.fq.gz reads.fq.gz > ovlp.paf
k8 ov-eval.js reads-to-ref.paf ovlp.paf
```

## <a name="asmvar"></a>Calling Variants from Haploid Assemblies

The **call** command of paftools.js calls variants from coordinate-sorted
assembly-to-reference alignment. It calls variants from the [cs tag][cs] and
identifies confident/callable regions as those covered by exactly one contig.
Here are example command lines:

```sh
minimap2 -cx asm5 -t8 --cs ref.fa asm.fa > asm.paf  # keeping this file is recommended; --cs required!
sort -k6,6 -k8,8n asm.paf > asm.srt.paf             # sort by reference start coordinate
k8 paftools.js call asm.srt.paf > asm.var.txt
```

Here is sample output:

```
V   chr1    2276040 2276041 1   60  c   g   LJII01000171.1  1217409 1217410 +
V   chr1    2280409 2280410 1   60  a   g   LJII01000171.1  1221778 1221779 +
V   chr1    2280504 2280505 1   60  a   g   LJII01000171.1  1221873 1221874 +
R   chr1    2325140 2436340
V   chr1    2325287 2325287 1   60  -   ct  LJII01000171.1  1272894 1272896 +
V   chr1    2325642 2325644 1   60  tt  -   LJII01000171.1  1273251 1273251 +
V   chr1    2326051 2326052 1   60  c   t   LJII01000171.1  1273658 1273659 +
V   chr1    2326287 2326288 1   60  c   t   LJII01000171.1  1273894 1273895 +
```

where a line starting with `R` gives regions covered by one query contig, and a
V-line encodes a variant in the following format: chr, start, end, query depth,
mapping quality, REF allele, ALT allele, query name, query start, end and the
query orientation. Generally, you should only look at variants where column 5
is one.

By default, when calling variants, "paftools.js call" ignores alignments 50kb
or shorter; when deriving callable regions, it ignores alignments 10kb or
shorter.  It uses two thresholds to avoid edge effects. These defaults are
designed for long-read assemblies. For short reads, both should be reduced.



[paf]: https://github.com/lh3/miniasm/blob/master/PAF.md
[cs]: https://github.com/lh3/minimap2#cs
[k8]: https://github.com/attractivechaos/k8
[maf]: https://genome.ucsc.edu/FAQ/FAQformat#format5
[pbsim]: https://github.com/pfaucon/PBSIM-PacBio-Simulator
[mason2]: https://github.com/seqan/seqan/tree/master/apps/mason2
# intervaltree

## Overview

An interval tree can be used to efficiently find a set of numeric intervals overlapping or containing another interval.

This library provides a basic implementation of an interval tree using C++ templates, allowing the insertion of arbitrary types into the tree.

## Usage

Add `#include "IntervalTree.h"` to the source files in which you will use the interval tree.

To make an IntervalTree to contain objects of class T, use:

```c++
vector<Interval<T> > intervals;
T a, b, c;
intervals.push_back(Interval<T>(2, 10, a));
intervals.push_back(Interval<T>(3, 4, b));
intervals.push_back(Interval<T>(20, 100, c));
IntervalTree<T> tree;
tree = IntervalTree<T>(intervals);
```

Now, it's possible to query the tree and obtain a set of intervals which are contained within the start and stop coordinates.

```c++
vector<Interval<T> > results;
tree.findContained(start, stop, results);
cout << "found " << results.size() << " overlapping intervals" << endl;
```

The function IntervalTree::findOverlapping provides a method to find all those intervals which are contained or partially overlap the interval (start, stop).

### Author: Erik Garrison <erik.garrison@gmail.com>

### License: MIT
Note to existing users: the iterator implementation has changed significantly
since we introduced the `locked_table` in [this
commit](https://github.com/efficient/libcuckoo/commit/2bedb3d0c811cd8b3adb3e78e2d2a28c66ba1d1d).
Please see the [`locked_table`
documentation](http://efficient.github.io/libcuckoo/classcuckoohash__map_1_1locked__table.html)
and [examples
directory](https://github.com/efficient/libcuckoo/tree/master/examples) for
information and examples of how to use iterators.

libcuckoo
=========

libcuckoo provides a high-performance, compact hash table that allows
multiple concurrent reader and writer threads.

The Doxygen-generated documentation is available at the
[project page](http://efficient.github.io/libcuckoo/).

Authors: Manu Goyal, Bin Fan, Xiaozhou Li, David G. Andersen, and Michael Kaminsky

For details about this algorithm and citations, please refer to
our papers in [NSDI 2013][1] and [EuroSys 2014][2]. Some of the details of the hashing
algorithm have been improved since that work (e.g., the previous algorithm
in [1] serializes all writer threads, while our current
implementation supports multiple concurrent writers), however, and this source
code is now the definitive reference.

   [1]: http://www.cs.cmu.edu/~dga/papers/memc3-nsdi2013.pdf "MemC3: Compact and Concurrent Memcache with Dumber Caching and Smarter Hashing"
   [2]: http://www.cs.princeton.edu/~mfreed/docs/cuckoo-eurosys14.pdf "Algorithmic Improvements for Fast Concurrent Cuckoo Hashing"

Requirements
================

This library has been tested on Mac OSX >= 10.8 and Ubuntu >= 12.04.

It compiles with clang++ >= 3.1 and g++ >= 4.7, however we strongly suggest
using the latest versions of both compilers, as they have greatly improved
support for atomic operations. Building the library requires the
autotools. Install them on Ubuntu

    $ sudo apt-get update && sudo apt-get install build-essential autoconf libtool

Building
==========

    $ autoreconf -fis
    $ ./configure
    $ make
    $ make install

Usage
==========

To build a program with the hash table, include
`libcuckoo/cuckoohash_map.hh` into your source file. If you want to
use CityHash, which we recommend, we have provided a wrapper
compatible with the `std::hash` type around it in the
`libcuckoo/city_hasher.hh` file. If compiling with CityHash, add the
`-lcityhash` flag. You must also enable C++11 features on your
compiler. Compiling the file `examples/count_freq.cpp` with g++
might look like this:

    $ g++ -std=c++11 examples/count_freq.cpp -lcityhash

The
[examples directory](https://github.com/efficient/libcuckoo/tree/master/examples)
contains some simple demonstrations of some of the basic features of the hash
table.

Tests
==========

The [tests directory](https://github.com/efficient/libcuckoo/tree/master/tests)
directory contains a number of tests and benchmarks of the hash table, which
also can serve as useful examples of how to use the table's various features.
After running `make all`, the entire test suite can be run with the `make check`
command. This will not run the benchmarks, which must be run individually. The
test executables, which have the suffix `.out`, can be run individually as well.

Issue Report
============

To let us know your questions or issues, we recommend you
[report an issue](https://github.com/efficient/libcuckoo/issues) on
github. You can also email us at
[libcuckoo-dev@googlegroups.com](mailto:libcuckoo-dev@googlegroups.com).

Licence
===========
Copyright (C) 2013, Carnegie Mellon University and Intel Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

---------------------------

CityHash (lib/city.h, lib/city.cc) is Copyright (c) Google, Inc. and
has its own license, as detailed in the source files.==============================
Mappy: Minimap2 Python Binding
==============================

Mappy provides a convenient interface to `minimap2
<https://github.com/lh3/minimap2>`_, a fast and accurate C program to align
genomic and transcribe nucleotide sequences.

Installation
------------

Mappy depends on `zlib <http://zlib.net>`_. It can be installed with `pip
<https://en.wikipedia.org/wiki/Pip_(package_manager)>`_:

.. code:: shell

	pip install --user mappy

or from the minimap2 github repo (`Cython <http://cython.org>`_ required):

.. code:: shell

	git clone https://github.com/lh3/minimap2
	cd minimap2
	python setup.py install

Usage
-----

The following Python script demonstrates the key functionality of mappy:

.. code:: python

	import mappy as mp
	a = mp.Aligner("test/MT-human.fa")  # load or build index
	if not a: raise Exception("ERROR: failed to load/build index")
	s = a.seq("MT_human", 100, 200)     # retrieve a subsequence from the index
	print(mp.revcomp(s))                # reverse complement
	for name, seq, qual in mp.fastx_read("test/MT-orang.fa"): # read a fasta/q sequence
		for hit in a.map(seq): # traverse alignments
			print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))

APIs
----

Mappy implements two classes and two global function.

Class mappy.Aligner
~~~~~~~~~~~~~~~~~~~

.. code:: python

	mappy.Aligner(fn_idx_in=None, preset=None, ...)

This constructor accepts the following arguments:

* **fn_idx_in**: index or sequence file name. Minimap2 automatically tests the
  file type. If a sequence file is provided, minimap2 builds an index. The
  sequence file can be optionally gzip'd. This option has no effect if **seq**
  is set.

* **seq**: a single sequence to index. The sequence name will be set to
  :code:`N/A`.

* **preset**: minimap2 preset. Currently, minimap2 supports the following
  presets: **sr** for single-end short reads; **map-pb** for PacBio
  read-to-reference mapping; **map-ont** for Oxford Nanopore read mapping;
  **splice** for long-read spliced alignment; **asm5** for assembly-to-assembly
  alignment; **asm10** for full genome alignment of closely related species. Note
  that the Python module does not support all-vs-all read overlapping.

* **k**: k-mer length, no larger than 28

* **w**: minimizer window size, no larger than 255

* **min_cnt**: mininum number of minimizers on a chain

* **min_chain_score**: minimum chaing score

* **bw**: chaining and alignment band width

* **best_n**: max number of alignments to return

* **n_threads**: number of indexing threads; 3 by default

* **extra_flags**: additional flags defined in minimap.h

* **fn_idx_out**: name of file to which the index is written. This parameter
  has no effect if **seq** is set.

* **scoring**: scoring system. It is a tuple/list consisting of 4, 6 or 7
  positive integers. The first 4 elements specify match scoring, mismatch
  penalty, gap open and gap extension penalty. The 5th and 6th elements, if
  present, set long-gap open and long-gap extension penalty. The 7th sets a
  mismatch penalty involving ambiguous bases.

.. code:: python

	mappy.Aligner.map(seq, seq2=None, cs=False, MD=False)

This method aligns :code:`seq` against the index. It is a generator, *yielding*
a series of :code:`mappy.Alignment` objects. If :code:`seq2` is present, mappy
performs paired-end alignment, assuming the two ends are in the FR orientation.
Alignments of the two ends can be distinguished by the :code:`read_num` field
(see Class mappy.Alignment below). Argument :code:`cs` asks mappy to generate
the :code:`cs` tag; :code:`MD` is similar. These two arguments might slightly
degrade performance and are not enabled by default.

.. code:: python

	mappy.Aligner.seq(name, start=0, end=0x7fffffff)

This method retrieves a (sub)sequence from the index and returns it as a Python
string. :code:`None` is returned if :code:`name` is not present in the index or
the start/end coordinates are invalid.

.. code:: python

	mappy.Aligner.seq_names

This property gives the array of sequence names in the index.

Class mappy.Alignment
~~~~~~~~~~~~~~~~~~~~~

This class describes an alignment. An object of this class has the following
properties:

* **ctg**: name of the reference sequence the query is mapped to

* **ctg_len**: total length of the reference sequence

* **r_st** and **r_en**: start and end positions on the reference

* **q_st** and **q_en**: start and end positions on the query

* **strand**: +1 if on the forward strand; -1 if on the reverse strand

* **mapq**: mapping quality

* **blen**: length of the alignment, including both alignment matches and gaps
  but excluding ambiguous bases.

* **mlen**: length of the matching bases in the alignment, excluding ambiguous
  base matches.

* **NM**: number of mismatches, gaps and ambiguous poistions in the alignment

* **trans_strand**: transcript strand. +1 if on the forward strand; -1 if on the
  reverse strand; 0 if unknown

* **is_primary**: if the alignment is primary (typically the best and the first
  to generate)

* **read_num**: read number that the alignment corresponds to; 1 for the first
  read and 2 for the second read

* **cigar_str**: CIGAR string

* **cigar**: CIGAR returned as an array of shape :code:`(n_cigar,2)`. The two
  numbers give the length and the operator of each CIGAR operation.

* **MD**: the :code:`MD` tag as in the SAM format. It is an empty string unless
  the :code:`MD` argument is applied when calling :code:`mappy.Aligner.map()`.

* **cs**: the :code:`cs` tag.

An :code:`Alignment` object can be converted to a string with :code:`str()` in
the following format:

::

	q_st  q_en  strand  ctg  ctg_len  r_st  r_en  mlen  blen  mapq  cg:Z:cigar_str

It is effectively the PAF format without the QueryName and QueryLength columns
(the first two columns in PAF).

Miscellaneous Functions
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

	mappy.fastx_read(fn, read_comment=False)

This generator function opens a FASTA/FASTQ file and *yields* a
:code:`(name,seq,qual)` tuple for each sequence entry. The input file may be
optionally gzip'd. If :code:`read_comment` is True, this generator yields
a :code:`(name,seq,qual,comment)` tuple instead.

.. code:: python

	mappy.revcomp(seq)

Return the reverse complement of DNA string :code:`seq`. This function
recognizes IUB code and preserves the letter cases. Uracil :code:`U` is
complemented to :code:`A`.
