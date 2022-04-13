# gfastats

The swiss army knife for genome assembly.

**gfastats** is a single fast and exhaustive tool for **summary statistics** and simultaneous \*fa\* (fasta, fastq, gfa [.gz]) genome assembly file **manipulation**.
**gfastats** also allows seamless fasta<>fastq<>gfa[.gz] conversion.

Typical fast\* metrics include:

- scaffold, contig and gap size
- number of scaffolds, contigs and gaps
- total length of scaffolds, contigs and gaps
- scaffold, contig, gap N50 and statistics (full N\*/NG\* statistics with the `--nstar-report` flag)
- area under the curve (AuN/AuNG) values for scaffolds, contigs and gaps
- average scaffold, contig, gap size
- largest scaffold, contig and gap
- base composition and GC content
- soft-masked base counts (lower case bases)

Typical gfa metrics include:

- Number of nodes and edges
- Average degree
- Number of connected components, and length of the largets connected component
- Number of dead ends
- Number of disconnected components, and their total length
- Number of bubbles

Metrics for each scaffold/contig can be generated with the `--seq-report` flag.

`Bed` coordinates and sizes of scaffolds, contigs and gaps can be outputted with the options `--out-coord` and `--out-size`. By default, `--out-coord` produces a full representation of the assembly in `agp` format.

Additionally, input can be filtered using scaffold lists or `bed` coordinate files with the options `--include-bed` and `--exclude-bed`.

Importantly, the filtered input can be outputted in any \*fa\* (fasta, fastq, gfa [.gz]) format.

## Installation

Either download one of the releases or `git clone https://github.com/vgl-hub/gfastats.git` and `make` in `gfastats` folder.

## Usage

`gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]`

To check out all options and flags use `gfastats -h`.

You can test some typical usage with the files in the `testFiles` folder, e.g.:

```
gfastats testFiles/random1.fasta -o gfa // converts fasta to gfa
gfastats testFiles/random2.gfa2.gfa -o fa // converts gfa to fasta
```

## Assembly manipulation

**gfastats** allows extensive assembly manipulation at the sequence level. Manipulation is achieved using a set of _instructions_ provided as an ordered list in a file to the option `-k` / `--swiss-army-knife`:

```
gfastats testFiles/random1.fasta -k testFiles/random1.instructions.sak -o gfa // reads fasta applies a set of instructions and outputs gfa
```

The _instructions_ are sequentially processed to generate the final output. Examples of _instructions_ are:

```
JOIN contig1+ contig2+ 50 [gap1] [scaffold1] [this is a new scaffold] // introduces a new gap of 50 bp between scaffold1 and scaffold2 with optional id gap1, effectively joining the two sequences into a new sequences named scaffold1 with an optional comment
SPLIT contig1+ contig2+ // splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them
```

The _instructions_ directly provide the list of edits that were introduced. The _instructions_ could be from an automated tool or from manual annotation.

A prime example of manipulations using input from an automated tool is overlaying AGP coordinates on top of the graph to generate new scaffolds, which can be achieved with:
```
gfastats input.fasta|input.gfa -a input.agp -o output.fasta|output.gfa
```

See the <a href="instructions/">instruction wiki</a> for a full list of _instructions_.

## Description

Please refer to **gfastats** paper for a complete description. Briefly, **gfastats** reads and stores any fasta<>fastq<>gfa[.gz] in gfa format. **gfastats** then builds a bidirected graph representation of the assembly using adjacency lists, where each node is a segment, and each edge is a gap (see figure below). The original sequence can be directly manipulated from the graph. Finally, walking the graph allows to generate different kinds of outputs, including manipulated assemblies and feature coordinates.

<p align="center">
    <img src="images/graph.png" alt="alt gfastats assembly graph" width="70%" />
</p>

## How to cite

If you use **gfastats** in your work, please cite:

Gfastats: conversion, evaluation and manipulation of genome sequences using assembly graphs

Giulio Formenti, Linelle Abueg, Angelo Brajuka, Nadolina Brajuka, Cristo Gallardo, Alice Giani, Olivier Fedrigo, Erich D. Jarvis

doi: https://doi.org/10.1101/2022.03.24.485682
# Instruction wiki

## JOIN

The JOIN instruction introduces a new gap of 50 bp between scaffold1 and scaffold2 (two paths) with id gap1, effectively joining the two sequences into a new sequences named new_scaffold with an optional comment.

```
JOIN    scaffold1+    scaffold2+    50  gap1    new_scaffold
JOIN    scaffold1(1:100)+    scaffold2(1:100)+    50  gap1    new_scaffold // optional subsetting
```

## SPLIT

The SPLIT instruction splits the scaffold containing contig1 and contig2, effectively removing the existing gap between them. Two optional comments can be provided.

```
SPLIT   contig1+    contig2+    scaffold1   scaffold2 [this is a new scaffold1] [this is a new scaffold2]
```

## EXCISE

The EXCISE instruction removes contig1 from its scaffold, leading leaving it unplaced and adding a gap of the given size with optional id gap1 between the original sequences

```
EXCISE  contig1  50  [gap1] // new 50 bp gap
```

## REMOVE

The REMOVE instruction removes contig1 from the segment set. If it is part of a path, the path is also removed.

```
REMOVE  contig1
```

## ERASE

The ERASE instruction trims off the sequence range specified from the given segment.

```
ERASE   contig1:10-100 // deletes contig1 sequence between the coordinates provided (in bed format)
```

## RVCP

The RVCP instruction reverse-complements contig1 sequence in place

```
RVCP    contig1
```

## INVERT

The INVERT instruction inverts contig1 sequence in place

```
INVERT  contig1
```

## Yet to be implemented

```
ADD contig3 contig1+ 50 contig2+ 50 ACGT // introduces a new contig named contig3 with sequence ACGT between contig1 and contig2 leaving 50bp gaps on each side
REPLACE contig1:20-24 ACGT // replaces the sequence at coordinates contig1:20-24 with ACGT
```
