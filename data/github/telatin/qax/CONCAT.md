# :package: Qiime2 Artifact eXtractor

[![Build Status](https://travis-ci.com/telatin/qax.svg?branch=main)](https://travis-ci.com/telatin/qax)
[![Repository Size](https://img.shields.io/github/languages/code-size/telatin/qax)](https://github.com/telatin/qax)
[![Latest release](https://img.shields.io/github/v/release/telatin/qax)](https://github.com/telatin/qax/releases)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/qax)](https://bioconda.github.io/recipes/qax/README.html)
[![Available via BioConda](https://img.shields.io/conda/vn/bioconda/qax)](https://bioconda.github.io/recipes/qax/README.html)
[![docs-badge](https://img.shields.io/badge/read-the%20docs-red)](https://telatin.github.io/qax)

## :book: Introduction

<a href="https://telatin.github.io/qax/" description="QAX documentation">
  <img alt="qax logo" align="right" width="333" height="153" src="https://raw.githubusercontent.com/telatin/qax/main/pages/qax.png">
</a>

Qiime2 is one of the most popular software used to analyze the output of metabarcoding experiment, and it introduced a unique data format in the bioinformatics scenario: the “_Qiime2 artifact_”.

Qiime2 artifacts are structured compressed archives containing a dataset (_e.g._, FASTQ reads, representative sequences in FASTA format, a phylogenetic tree in Newick format, etc.) and an exhaustive set of metadata (including the command that generated it, information on the execution environment, citations on the used software, and all the metadata of the artifacts used to produce it).

While artifacts can improve the shareability and reproducibility of Qiime workflows, they are less easily integrated with general bioinformatics pipelines, and even accessing metadata in the artifacts requires the full Qiime2 installation (not to mention that every release of Qiime2 will produce incompatible artifacts). Qiime Artifact Extractor (qxa) allows to easily interface with Qiime2 artifacts from the command line, without needing the full Qiime2 environment installed.


## Citation

If you use this tool, please cite

> Telatin A (2021) **Qiime Artifact eXtractor (qax): A Fast and Versatile Tool to Interact with Qiime2 Archives**. BioTech 10: 5. 
> Available: (doi.org/10.3390/biotech10010005)[http://dx.doi.org/10.3390/biotech10010005]

## :floppy_disk: Download and installation

Pre-compiled binaries are the fastest and easiest way to get _qax_. To get the latest version,
use the following command, otherwise check the [stable releases](https://github.com/telatin/qax/releases).  


```
# From linux
wget "https://github.com/telatin/qax/raw/main/bin/qax"
chmod +x qax

# From macOS
wget -O qax "https://github.com/telatin/qax/raw/main/bin/qax_mac"
chmod +x qax
```

Alternatively, you can install _qax_ from BioConda, if you have _conda_ installed:
```
conda install -c conda-forge -c bioconda qax
```

## :book: Usage

`qax` has four subprograms (general syntax is `qax [program] [program-arguments]`):

- **list** (default): list artifact(s) properties
- **citations**: extract citations in BibTeX format
- **extract**: extract artifact _data_ files
- **provenance**: describe artifact provenance, or generate its graph
- **view**: print the content of an artifact (eg. dna-sequences.fasta) to the terminal


### :page_facing_up: list


* See [**qax list** full documentation](pages/list.md)

This is the default module, and can be used to list the properties of one or more artifacts.

Some features:
* Supports multiple files at once
* 100X times faster than Qiime2
* Can be used to find an artifact given the ID

Example:
```
qax_mac -b -u input/*.*
┌───────────────────────────┬────────────────┬─────────────────────────┬─────────────────────────────┐
│ ID                        │ Basename       │ Type                    │ Format                      │
├───────────────────────────┼────────────────┼─────────────────────────┼─────────────────────────────┤
│ bb1b2e93-...-2afa2110b5fb │ rep-seqs.qza   │ FeatureData[Sequence]   │ DNASequencesDirectoryFormat │
│ 313a0cf3-...-befad4ebf2f3 │ table.qza      │ FeatureTable[Frequency] │ BIOMV210DirFmt              │
│ 35c32fe7-...-85ef27545f00 │ taxonomy.qzv   │ Visualization           │ HTML                        │
└───────────────────────────┴────────────────┴─────────────────────────┴─────────────────────────────┘
```

### :page_facing_up: extract


* See [**qax extract** full documentation](pages/extract.md)

This program extract the content of an artifact. By default, if a single file is present it will be extracted in the specified path. If multiple files are present, a directory containing them will be created instead.

_Example:_
```
# Extract representative sequences (will be called rep-seqs.fasta)
qax x -o ./ rep-seqs.qza

# Extract a visualization (a folder called "taxonomy" will be created)
qax x -o ./ taxonomy.qzv
```

### :page_facing_up: citations

* See [**qax citations** full documentation](pages/cite.md)

Each Qiime module provides the citations for the software and resources that it uses, storing the citations in BibTeX format inside the artifacts. The cite module allows to extract all the citations from a list of artifacts, removing the duplicates, thus effectively allowing to prepare the bibliography for a complete Qiime2 analysis.

_Example:_
```
qax c files/*.qza > bibliography.bib
```

### :page_facing_up: provenance

* See [**qax provenance** full documentation](pages/provenance.md)

This program allows to print the provenance of an artifact, or to produce a [publication grade graph](pages/qax-provenance.png) of the provenance.

_Example:_
```
# To view a summary
qax p taxonomy.qzv

# To save the plot
qax p -o graph.dot taxonomy.qza
```


### :page_facing_up: view

* See [**qax view** full documentation](pages/view.md)

This program allows to print the content of an artifact data file to the terminal.
If the artifact contains a single file, it will be printed. Otherwise the user can specify one or multiple files to be printed, and if none
is specified, a list of files will be printed.

```
# Example: count the number of representative sequences
qax view rep-seqs.qza | grep -c '>'
```

### :page_facing_up: make

To create a _visualization artifact_ from a folder with a website (index.html must
  be present).

```
qax make -o report.qza /path/to/report_dir/
```
# qax citations

This subprogram will print the BibTeX citations (deduplicated) from a list of files.
Can be abbreviated with **c** (_e. g._, `qax c ...`).

### Synopsis:

```
Usage: citations [options] [<inputfile> ...]

Options:
  -o, --output FILE      Save BibTeX citation to FILE
  -r, --recurse-parents  Retrieve BibTeX citations also from parents
  -f, --force-artifacts  Try to parse artifacts with non canonical extensions
  -v, --verbose          Verbose output
  -h, --help             Show this help
```

#### Recurse

By default only the citations for the given artifacts are printed. This means that feeding all the artifacts used for a specific
analysis will print the full bibliography:

```
qax citations  -o bibliography.bib /path/to/*.qz*
```

Since each artifact also contains the bibliography used to generate its parent artifacts, one can also use _citations_ in 
recursive mode to gather them:

```
qax citations  -r -o bibliography.bib /path/to/*.qz*
```

Since the citations are non redundant, there is no harm in using `-r` in any case.


### Performance


Using _hyperfine_ to compare the performance with Qiime2:

```
hyperfine --export-markdown docs/list_speed.md "qiime tools citations input/taxonomy.qzv" "./bin/qax citations input/taxonomy.qzv" 
```

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `./bin/qax citations input/taxonomy.qzv` | 20.3 ± 8.4 | 8.5 | 31.4 | 1.00 |
| `qiime tools citations  input/taxonomy.qzv` | 2616.7 ± 115.1 | 2476.7 | 3014.3 | 128.69 ± 53.76 |
# qax extract

This subprogram will extract data files from artifacts.
Can be abbreviated with **x** (_e. g._, `qax x ...`).

![qax extract](qax-extract.png)

### Synopsis:

```
Usage: extract [options] [<inputfile> ...]

Extract the artifact data. If multiple files are present, a new directory
will be created (artifact name), otherwise the artifact name will be used
to rename the single file (unless -k is specified).

Options:
  -o, --outdir DIR       Output directory [default: .]
  -k, --keepname         Keep original file names instead of using artifact's basename [default: false]
  -v, --verbose          Verbose output
  -h, --help             Show this help
```

### Behaviour

When a datafile contains a single file (feature table, representative sequences,...) it will be extracted to 
the output directory inheriting the name from the artifact's name.
Otherwise (visualizations, ...) it will be extracted to a subdirectory with the artifact's name.| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `./bin/qax_mac input/*` | 15.4 ± 6.6 | 8.4 | 29.3 | 1.00 |
| `for i in input/*; do qiime tools peek $i; done` | 4849.7 ± 306.8 | 4479.2 | 5297.2 | 315.62 ± 136.68 |
| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `qiime tools peek input/taxonomy.qzv` | 1.981 ± 0.056 | 1.925 | 2.087 | 109.58 ± 34.30 |
| `./bin/qax_mac input/taxonomy.qzv` | 0.018 ± 0.006 | 0.011 | 0.031 | 1.00 |
# qax view

This subprogram will print the BibTeX citations (deduplicated) from a list of files.
Can be abbreviated with **c** (_e. g._, `qax c ...`).

### Synopsis:

```
Usage: view [options] <inputfile> [<datafile>...]

View (a la `cat`) the content of data files in one artifact.

Options:
  -f, --force            Accept artifact with non caninical extension
  -v, --verbose          Verbose output
  -h, --help             Show this help
```

### Single file artifacts
Example:
```
qax view rep-seq.qza | grep -c '>'
```

### Artifact with multiple files

The file (or files) to extract have to be specified:

```
qax view taxonomy.qzv index.html | grep -c "body"
```# qax list

This is the default subprogram.

### Synopsis:

```
Usage: list [options] [<inputfile> ...]

Options:
  -a, --all              Show all fields [default: false]
  --abs                  Show absolute paths [default: false]
  -b, --basename         Use basename instead of path [default: false]
  -s, --sortby SORT      Column to sort (uuid, type, format, date) [default: type]
  -f, --force            Accept non standard extensions
  --verbose              Verbose output
  --help                 Show this help

Nice output (default):
  -d, --datetime         Show artifact's date time [default: false]
  -u, --uuid             Show uuid [default: false]

Raw table output:
  -r, --rawtable         Print a CSV table (-s) with all fields [default: false]   
  -h, --header           Print header [default: false]
  -z, --separator SEP    Separator when using --rawtable [default: tab]
```

This is the default module, so its name can be omitted. This means that the following two commands are equivalent:

```
qxa list directory/*.qza
qxa directory/*.qza
```

### File path

By default the file path will be printed as received by the program. Witht `-b` (`--basename`) the basename will be printed instead,
and wiht `--abs` the absolute path.

### Output

The program will prepare a table with the artifacts metadata, supporting a nice Unicode table (default) or a raw computer-friendly tabular
output (enabled with `-r`). 

#### Nice (Unicode) table

By default, only _path_, _type_ and _format_ of the artifact are printed:
```
┌────────────────────┬───────────────────────┬─────────────────────────────┐
│ Path               │ Type                  │ Format                      │
├────────────────────┼───────────────────────┼─────────────────────────────┤
│ input/rep-seqs.qza │ FeatureData[Sequence] │ DNASequencesDirectoryFormat │
│ input/taxonomy.qzv │ Visualization         │ HTML                        │
└────────────────────┴───────────────────────┴─────────────────────────────┘
```

With `-u` the UUID is added, with `-d` the generation timestamp is added. 
With `-a` all fields are printed:

```
┌──────────────────────────────────────┬──────────────┬───────────────────────┬─────────────────────────────┬───────────┬──────────────────┬───────┐
│ ID                                   │ Basename     │ Type                  │ Format                      │ Version   │ Date             │ Files │
├──────────────────────────────────────┼──────────────┼───────────────────────┼─────────────────────────────┼───────────┼──────────────────┼───────┤
│ bb1b2e93-0c45-4c8e-a140-2afa2110b5fb │ rep-seqs.qza │ FeatureData[Sequence] │ DNASequencesDirectoryFormat │ 2019.10.0 │ 2020-11-12;17:23 │ 1     │
│ 35c32fe7-3eb5-4b31-aa34-85ef27545f00 │ taxonomy.qzv │ Visualization         │ HTML                        │ 2019.10.0 │ 2020-11-12;17:23 │ 18    │
└──────────────────────────────────────┴──────────────┴───────────────────────┴─────────────────────────────┴───────────┴──────────────────┴───────┘
```

#### Raw tabular output

With the raw tabular output all the fields are printed, to allow easier parsing. By default the fields are separated by tabs, but a different
separator can be specified (`-s SEP`).

The fields are:
1. UUID (_e.g._ b1b2e93-0c45-4c8e-a140-2afa2110b5fb)
2. File path (_e.g._ rep-seqs.qza)
3. Artifact Type (_e.g._ FeatureData[Sequence])
4. Artifact Format (_e.g._ DNASequencesDirectoryFormat)
5. Qiime2 Version (_e.g._ 2019.10.0)
6. Timestamp (_e.g._ 2020-11-12;17:23 )
7. Data files contained (_e.g._ 1)


### Usage examples

To find an artifact given the ID:

```
find ${PATH_TO_SCAN} -name "*.qz?" | xargs qxa | grep ${UUID}
```

### Performance

Using _hyperfine_ to compare the performance with Qiime2:

```
hyperfine --export-markdown docs/list_speed.md "qiime tools peek input/taxonomy.qzv" "./bin/qax input/taxonomy.qzv" 
```


| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `qiime tools peek input/taxonomy.qzv` | 1.981 ± 0.056 | 1.925 | 2.087 | 109.58 ± 34.30 |
| `./bin/qax input/taxonomy.qzv` | 0.018 ± 0.006 | 0.011 | 0.031 | 1.00 |
| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `./bin/qax_mac citations input/taxonomy.qzv` | 20.3 ± 8.4 | 8.5 | 31.4 | 1.00 |
| `qiime tools citations  input/taxonomy.qzv` | 2616.7 ± 115.1 | 2476.7 | 3014.3 | 128.69 ± 53.76 |
# qax provenance

This subprogram will print the provenance of an artifact, optionally saved as graph.
Can be abbreviated with **p** (_e. g._, `qax p ...`).

![qax provenance](qax-provenance.png)

### Synopsis:

```
Usage: provenance [options] <inputfile>

Print ancestry of an artifact

Options:
  -o, --dotfile DOTFILE  Save graphviz graph (dot format)
  -p, --pdf              Generate PDF of the graph (requires -o and "dot" installed)
  -f, --font FONT        Font name for the PDF [default: Arial]
  -v, --verbose          Verbose output
  -h, --help             Show this help
```


### Text output

By default only a textual summary will be printed like: 
```
6e46b89a-0fcb-4a5f-860c-ec8ce26e7140    top     2020-01-24 09:01        reference_reads         [inputs: <none>]
ffdfda2e-be68-47f6-9bf5-07ccdd869689    top     2020-01-24 09:01        reference_taxonomy      [inputs: <none>]
71429e7c-6489-42d6-bbca-f4a6d9b85b48    |       2020-01-24 09:01        classifier              [inputs: 6e46b89a-0fcb-4a5f-860c-ec8ce26e7140;ffdfda2e-be68-47f6-9bf5-07ccdd869689;null]
a1ad1da7-8cc8-439b-bec5-c66a1125786f    top     2020-02-28 11:02        per_sample_sequences    [inputs: <none>]
3c984d76-82a7-4ff6-b64b-561834df9327    |       2020-02-28 11:02        demultiplexed_seqs      [inputs: a1ad1da7-8cc8-439b-bec5-c66a1125786f]
bb1b2e93-0c45-4c8e-a140-2afa2110b5fb    |       2020-02-28 13:02        reads                   [inputs: 3c984d76-82a7-4ff6-b64b-561834df9327]
9df28153-4e38-44aa-bbc4-58a37d699580    |       2020-02-28 13:02        input                   [inputs: bb1b2e93-0c45-4c8e-a140-2afa2110b5fb;71429e7c-6489-42d6-bbca-f4a6d9b85b48]
35c32fe7-3eb5-4b31-aa34-85ef27545f00    child   2020-02-28 13:02        visualization           [inputs: 9df28153-4e38-44aa-bbc4-58a37d699580:input.tsv]
```
### Publication grade graphs

To save the output as graph, the `-o FILE` option can be used. This will generate a dot file (graphviz), that can be rendered using the `dot`program (also [online](https://dreampuf.github.io/GraphvizOnline/#digraph%20G%20%7B%0A%0A%20%20subgraph%20cluster_0%20%7B%0A%20%20%20%20style%3Dfilled%3B%0A%20%20%20%20color%3Dlightgrey%3B%0A%20%20%20%20node%20%5Bstyle%3Dfilled%2Ccolor%3Dwhite%5D%3B%0A%20%20%20%20a0%20-%3E%20a1%20-%3E%20a2%20-%3E%20a3%3B%0A%20%20%20%20label%20%3D%20%22process%20%231%22%3B%0A%20%20%7D%0A%0A%20%20subgraph%20cluster_1%20%7B%0A%20%20%20%20node%20%5Bstyle%3Dfilled%5D%3B%0A%20%20%20%20b0%20-%3E%20b1%20-%3E%20b2%20-%3E%20b3%3B%0A%20%20%20%20label%20%3D%20%22process%20%232%22%3B%0A%20%20%20%20color%3Dblue%0A%20%20%7D%0A%20%20start%20-%3E%20a0%3B%0A%20%20start%20-%3E%20b0%3B%0A%20%20a1%20-%3E%20b3%3B%0A%20%20b2%20-%3E%20a3%3B%0A%20%20a3%20-%3E%20a0%3B%0A%20%20a3%20-%3E%20end%3B%0A%20%20b3%20-%3E%20end%3B%0A%0A%20%20start%20%5Bshape%3DMdiamond%5D%3B%0A%20%20end%20%5Bshape%3DMsquare%5D%3B%0A%7D)).

If the user has the `dot` program installed, a PDF can be automatically generated adding the `-p` switch.| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `qiime tools peek input/rep-seqs.qza` | 1.681 ± 0.458 | 1.391 | 2.961 | 90.79 ± 50.04 |
| `./bin/qax_mac input/rep-seqs.qza` | 0.019 ± 0.009 | 0.008 | 0.034 | 1.00 |
| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `qiime tools peek input/rep-seqs.qza; qiime tools peek input/taxonomy.qzv` | 3.325 ± 0.317 | 3.037 | 4.180 | 175.16 ± 81.37 |
| `./bin/qax_mac input/rep-seqs.qza input/taxonomy.qzv` | 0.019 ± 0.009 | 0.009 | 0.033 | 1.00 |
# `qax` source

![Commit](https://img.shields.io/github/last-commit/telatin/qax)
![Version 0.1.0](https://img.shields.io/badge/version-0.2.0-blue)
 | Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `vsearch --derep_fulllength ./input/filt.fa.gz --output /tmp/vsearch.fa` | 179.1 ± 9.6 | 170.1 | 198.0 | 1.06 ± 0.10 |
| `./bin/derep_Darwin ./input/filt.fa.gz > /tmp/derep.fa` | 169.2 ± 13.1 | 150.0 | 211.8 | 1.00 |
| `perl ./test/uniq.pl ./input/filt.fa.gz > /tmp/uniq.fa ` | 956.9 ± 64.5 | 866.4 | 1063.1 | 5.65 ± 0.58 |
| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `perl uniq.pl ./input/*.fa*` | 935.7 ± 38.4 | 886.1 | 1027.9 | 8.39 ± 0.36 |
| `derep_Darwin ./input/*.fa*` | 111.5 ± 1.6 | 108.9 | 115.1 | 1.00 |
# Tests

- [QAX Documentation](https://telatin.github.io/qax)


## Test script

The `all.sh` script tests the all the basic functionalities of **qax**.

```
[0] Synopsis
     OK Help for 'list': found
     OK Help for 'citations': found
     OK Help for 'provenance': found
     OK Help for 'extract': found
     OK Help for 'make': found
[1] List
     OK list (no params)
     OK list (-b)
     OK list (-b -u)
     OK list (-r)
[2] Extract
     OK Extract single artifact
     OK Extract multiple files artifact
     OK Extract multiple files artifact (subdirectory)
[3] Citations
     OK 1 citation found in taxonomy.qzv
     OK 8 _total_ citations found in taxonomy.qzv
[4] Provenance
     OK 3 parents for taxonomy.qzv
     OK 1 child for taxonomy.qzv
     OK output graph found
     OK output graph has 7 nodes
[5] View
     OK 769 seqs for rep-seqs.qza
[6] Make
     OK created visualization artifact
     OK correct UUID used and found 6 times
```
# QAX: the Qiime2 Artifacts eXtractor


[![Build Status](https://travis-ci.com/telatin/qax.svg?branch=main)](https://travis-ci.com/telatin/qax)
[![Repository Size](https://img.shields.io/github/languages/code-size/telatin/qax)](https://github.com/telatin/qax)
[![Latest release](https://img.shields.io/github/v/release/telatin/qax)](https://github.com/telatin/qax/releases)
[![Available via BioConda](https://img.shields.io/conda/vn/bioconda/qax)](https://bioconda.github.io/recipes/qax/README.html)
[![BioConda Downloads](https://img.shields.io/conda/dn/bioconda/qax)](https://bioconda.github.io/recipes/qax/README.html)

- Website: <https://telatin.github.io/qax/>
- Github: <https://github.com/telatin/qax>
- Paper: <https://doi.org/10.3390/biotech10010005>

## Introduction

<img alt="qax logo" align="right" width="333" height="153" src="https://raw.githubusercontent.com/telatin/qax/main/pages/qax.png">

Qiime2 is one of the most popular software used to analyze the output of metabarcoding experiment, and it introduced a unique data format in the bioinformatics scenario: the “_Qiime2 artifact_”.

Qiime2 artifacts are structured compressed archives containing a dataset (_e.g._, FASTQ reads, representative sequences in FASTA format, a phylogenetic tree in Newick format, etc.) and an exhaustive set of metadata (including the command that generated it, information on the execution environment, citations on the used software, and all the metadata of the artifacts used to produce it).

While artifacts can improve the shareability and reproducibility of Qiime workflows, they are less easily integrated with general bioinformatics pipelines, and even accessing metadata in the artifacts requires the full Qiime2 installation (not to mention that every release of Qiime2 will produce incompatible artifacts). Qiime Artifact Extractor (qxa) allows to easily interface with Qiime2 artifacts from the command line, without needing the full Qiime2 environment installed.

## Functions


`qax` has different subprograms (and the general syntax is `qax [program] [program-arguments]`):

- **list** (default): list artifact(s) properties
- **citations**: extract citations in BibTeX format
- **extract**: extract artifact _data_ files
- **provenance**: describe artifact provenance, or generate its graph
- **view**: print the content of an artifact (eg. dna-sequences.fasta) to the terminal
- **make**: create a visualization artifact from HTML

## Citation

Telatin, A. **Qiime Artifact eXtractor (qax): A Fast and Versatile Tool to Interact with Qiime2 Archives.** BioTech [doi.org/10.3390/biotech10010005](https://doi.org/10.3390/biotech10010005)
---
sort: 2
permalink: /usage
---
# General usage

`qax` is composed by five subprogram, and the general syntax is:

```
qax [program] parameters
```

The programs are:
- **list** (it's the default action and can be omitteed)
- **extract** or **x**
- **citations** or **c**
- **provenance** or **p**
- **view** or **v**

## list



This is the default module, and can be used to list the properties of one or more artifacts.

Some features:
* Supports multiple files at once
* 100X times faster than Qiime2
* Can be used to find an artifact given the ID

Example:
```
qax_mac -b -u input/*.*
┌───────────────────────────┬────────────────┬─────────────────────────┬─────────────────────────────┐
│ ID                        │ Basename       │ Type                    │ Format                      │
├───────────────────────────┼────────────────┼─────────────────────────┼─────────────────────────────┤
│ bb1b2e93-...-2afa2110b5fb │ rep-seqs.qza   │ FeatureData[Sequence]   │ DNASequencesDirectoryFormat │
│ 313a0cf3-...-befad4ebf2f3 │ table.qza      │ FeatureTable[Frequency] │ BIOMV210DirFmt              │
│ 35c32fe7-...-85ef27545f00 │ taxonomy.qzv   │ Visualization           │ HTML                        │
└───────────────────────────┴────────────────┴─────────────────────────┴─────────────────────────────┘
```

## extract



This program extract the content of an artifact. By default, if a single file is present it will be extracted in the specified path. If multiple files are present, a directory containing them will be created instead.

_Example:_
```
# Extract representative sequences (will be called rep-seqs.fasta)
qax x -o ./ rep-seqs.qza

# Extract a visualization (a folder called "taxonomy" will be created)
qax x -o ./ taxonomy.qzv
```

## citations

Each Qiime module provides the citations for the software and resources that it uses, storing the citations in BibTeX format inside the artifacts. The cite module allows to extract all the citations from a list of artifacts, removing the duplicates, thus effectively allowing to prepare the bibliography for a complete Qiime2 analysis.

_Example:_
```
qax c files/*.qza > bibliography.bib
```

## provenance

This program allows to print the provenance of an artifact, or to produce a [publication grade graph](docs/qax-provenance.png) of the provenance.

_Example:_
```
# To view a summary
qax p taxonomy.qzv

# To save the plot
qax p -o graph.dot taxonomy.qza
```


## view

This program allows to print the content of an artifact data file to the terminal.
If the artifact contains a single file, it will be printed. Otherwise the user can specify one or multiple files to be printed, and if none
is specified, a list of files will be printed.

```
# Example: count the number of representative sequences
qax view rep-seqs.qza | grep -c '>'
```

## make

This program converts a directory containing a website (index.html) into a
_visualization_ artifact.

```
qax make -o report.qzv /path/to/webpage/
```
---
sort: 1
permalink: /installation
---

# Installation

## Pre-compiled binaries

Pre-compiled binaries are the fastest and easiest way to get _qax_. To get the latest version,
use the following command, otherwise check the [stable releases](https://github.com/telatin/qax/releases).  


```
# From linux
wget "https://github.com/telatin/qax/raw/main/bin/qax"
chmod +x qax

# From macOS
wget -O qax "https://github.com/telatin/qax/raw/main/bin/qax_mac"
chmod +x qax
```

## Install via Miniconda

```note
Miniconda installation has been tested on MacOS and Linux, but being _qax_ a single binary, if the precompiled works for you we recommend it.
```

Alternatively, you can install _qax_ from BioConda, if you have _conda_ installed:

```
conda install -c conda-forge -c bioconda qax
```


---
sort: 3
permalink: /examples
---
# Usage examples

## Find artifact by UUID

If you know that an artifact has a specific UUID but do not know the exact location, you can
find it out combining `qax` with `find`.  

```
UUID="bb1b2e93-0c45-4c8e-a140-2afa2110b5fb"
find /path -name "*.qz?" | qax list -u  | grep $UUID
```

## Prepare the bibliography of a whole analysis

If we have a list of artifacts (_e. g._, `*.qz?`) we can extract the full bibliography and save it to file:
```
qax citations --output bibliography.bib *.qz?
```

## Check how many representative sequences have been produced

The _view_ subprogram will extract the main file and print it to the standard output, so we can combine 
it with [grep](https://linux.die.net/man/1/grep) to count how many lines have a ">" sign.t
```
qax view input/rep-seqs.qza | grep -c '>'
```
## v0.9.5

* `view` appears in the splashscreen (was available also before)
* `extract` can create missing directory
* Minor internal changes (can be compiled with `nimble build`)
---
sort: 2
---

# Releases


{% include list.liquid %}
