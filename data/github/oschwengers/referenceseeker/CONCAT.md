[![DOI](https://joss.theoj.org/papers/10.21105/joss.01994/status.svg)](https://doi.org/10.21105/joss.01994)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/referenceseeker/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/referenceseeker.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/referenceseeker.svg)
[![PyPI](https://img.shields.io/pypi/v/referenceseeker.svg)](https://pypi.org/project/referenceseeker)
![PyPI - Status](https://img.shields.io/pypi/status/referenceseeker.svg)
[![Conda](https://img.shields.io/conda/v/bioconda/referenceseeker.svg)](http://bioconda.github.io/recipes/referenceseeker/README.html)
![Python package](https://github.com/oschwengers/referenceseeker/workflows/Python%20package/badge.svg?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4415843.svg)](https://doi.org/10.5281/zenodo.4415843)

# ReferenceSeeker: rapid determination of appropriate reference genomes

## Contents

- [Description](#description)
- [Input & Output](#input-output)
- [Installation](#installation)
  - [BioConda](#bioconda)
  - [GitHub](#github)
- [Usage](#usage)
- [Examples](#examples)
- [Databases](#databases)
  - [RefSeq](#refseq-based)
  - [Custom](#custom-database)
- [Dependencies](#dependencies)
- [Citation](#citation)

## Description

ReferenceSeeker determines closely related reference genomes following a scalable hierarchical approach combining an fast kmer profile-based database lookup of candidate reference genomes and subsequent computation of specific average nucleotide identity (ANI) values for the rapid determination of suitable reference genomes.

ReferenceSeeker computes kmer-based genome distances between a query genome and potential reference genome candidates via Mash (Ondov et al. 2016). For resulting candidates ReferenceSeeker subsequently computes (bidirectional) ANI values picking genomes meeting community standard thresholds by default (ANI >= 95 % & conserved DNA >= 69 %) (Goris, Konstantinos et al. 2007) ranked by the product of ANI and conserved DNA values to take into account both genome coverage and identity.

Custom databases can be built with local genomes. For further convenience, we provide pre-built databases with sequences from RefSeq (<https://www.ncbi.nlm.nih.gov/refseq>), GTDB and PLSDB copmrising the following taxa:

- bacteria
- archaea
- fungi
- protozoa
- viruses

as well as *plasmids*.

The reasoning for subsequent calculations of both ANI and conserved DNA values is that Mash distance values correlate well with ANI values for closely related genomes, however the same is not true for conserved DNA values. A kmer fingerprint-based comparison alone cannot distinguish if a kmer is missing due to a SNP, for instance or a lack of the kmer-comprising subsequence. As DNA conservation (next to DNA identity) is very important for many kinds of analyses, *e.g.* reference based SNP detections, ranking potential reference genomes based on a mash distance alone is often not sufficient in order to select the most appropriate reference genomes. If desired, ANI and conserved DNA values can be computed bidirectionally.

![Mash D vs. ANI / conDNA](mash-ani-cdna.mini.png?raw=true)

## Input & Output

### Input

Path to a taxon database and a draft or finished genome in (zipped) fasta format:

```bash
$ referenceseeker ~/bacteria GCF_000013425.1.fna
```

### Output

Tab separated lines to STDOUT comprising the following columns:

Unidirectionally (query -> references):

- RefSeq Assembly ID
- Mash Distance
- ANI
- Conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism

```bash
#ID    Mash Distance    ANI    Con. DNA    Taxonomy ID    Assembly Status    Organism
GCF_000013425.1    0.00000    100.00    100.00    93061    complete    Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1    0.00002    100.00    99.89     46170    complete    Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1    0.00004    100.00    99.57     93061    complete    Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2    0.00016    100.00    99.28     1280     complete    Staphylococcus aureus FDAARGOS_10
GCF_003595465.1    0.00185    99.86     96.81     1280     complete    Staphylococcus aureus USA300-SUR6
GCF_003595385.1    0.00180    99.87     96.80     1280     complete    Staphylococcus aureus USA300-SUR2
GCF_003595365.1    0.00180    99.87     96.80     1280     complete    Staphylococcus aureus USA300-SUR1
GCF_001956815.1    0.00180    99.87     96.80     46170    complete    Staphylococcus aureus subsp. aureus USA300_SUR1
...
```

Bidirectionally (query -> references [QR] & references -> query [RQ]):

- RefSeq Assembly ID
- Mash Distance
- QR ANI
- QR Conserved DNA
- RQ ANI
- RQ Conserved DNA
- NCBI Taxonomy ID
- Assembly Status
- Organism

```bash
#ID    Mash Distance    QR ANI    QR Con. DNA    RQ ANI    RQ Con. DNA    Taxonomy ID    Assembly Status    Organism
GCF_000013425.1    0.00000    100.00    100.00    100.00    100.00    93061    complete    Staphylococcus aureus subsp. aureus NCTC 8325
GCF_001900185.1    0.00002    100.00    99.89     100.00    99.89     46170    complete    Staphylococcus aureus subsp. aureus HG001
GCF_900475245.1    0.00004    100.00    99.57     99.99     99.67     93061    complete    Staphylococcus aureus subsp. aureus NCTC 8325 NCTC8325
GCF_001018725.2    0.00016    100.00    99.28     99.95     98.88     1280     complete    Staphylococcus aureus FDAARGOS_10
GCF_001018915.2    0.00056    99.99     96.35     99.98     99.55     1280     complete    Staphylococcus aureus NRS133
GCF_001019415.2    0.00081    99.99     94.47     99.98     99.36     1280     complete    Staphylococcus aureus NRS146
GCF_001018735.2    0.00096    100.00    94.76     99.98     98.58     1280     complete    Staphylococcus aureus NRS137
GCF_003354885.1    0.00103    99.93     96.63     99.93     96.66     1280     complete    Staphylococcus aureus 164
...
```

## Installation

ReferenceSeeker can be installed via Conda and Git(Hub). In either case, a taxon database must be downloaded which we provide for download at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3562004.svg)](https://doi.org/10.5281/zenodo.3562004)
For more information have a look at [Databases](#databases).

### BioConda

The preferred way to install and run ReferenceSeeker is [Conda](https://conda.io/docs/install/quick.html) using the [Bioconda](https://bioconda.github.io/) channel:

```bash
$ conda install -c bioconda referenceseeker
$ referenceseeker --help
```

### GitHub

Alternatively, you can use this raw GitHub repository:

1. install necessary Python dependencies (if necessary)
2. clone the latest version of the repository
3. install necessary 3rd party executables (Mash, MUMmer4)

```bash
$ pip3 install --user biopython xopen
$ git clone https://github.com/oschwengers/referenceseeker.git
$ # install Mash & MUMmer
$ ./referenceseeker/bin/referenceseeker --help
```

### Test

To test your installation we prepared a tiny mock database comprising 4 `Salmonella spp` genomes and a query assembly (SRA: SRR498276) in the `tests` directory:

```bash
$ git clone https://github.com/oschwengers/referenceseeker.git

  # GitHub installation
$ ./referenceseeker/bin/referenceseeker referenceseeker/test/db referenceseeker/test/data/Salmonella_enterica_CFSAN000189.fasta

  # BioConda installation
$ referenceseeker referenceseeker/test/db referenceseeker/test/data/Salmonella_enterica_CFSAN000189.fasta
```

Expected output:

```bash
#ID    Mash Distance    ANI    Con. DNA    Taxonomy ID    Assembly Status    Organism
GCF_000439415.1    0.00003    100.00    99.55    1173427    complete    Salmonella enterica subsp. enterica serovar Bareilly str. CFSAN000189
GCF_900205275.1    0.01522    98.61     83.13    90370      complete    Salmonella enterica subsp. enterica serovar Typhi
```

## Usage

Usage:

```bash
usage: referenceseeker [--crg CRG] [--ani ANI] [--conserved-dna CONSERVED_DNA]
                       [--unfiltered] [--bidirectional] [--help] [--version]
                       [--verbose] [--threads THREADS]
                       <database> <genome>

Rapid determination of appropriate reference genomes.

positional arguments:
  <database>            ReferenceSeeker database path
  <genome>              target draft genome in fasta format

Filter options / thresholds:
  These options control the filtering and alignment workflow.

  --crg CRG, -r CRG     Max number of candidate reference genomes to pass kmer
                        prefilter (default = 100)
  --ani ANI, -a ANI     ANI threshold (default = 0.95)
  --conserved-dna CONSERVED_DNA, -c CONSERVED_DNA
                        Conserved DNA threshold (default = 0.69)
  --unfiltered, -u      Set kmer prefilter to extremely conservative values
                        and skip species level ANI cutoffs (ANI >= 0.95 and
                        conserved DNA >= 0.69
  --bidirectional, -b   Compute bidirectional ANI/conserved DNA values
                        (default = False)

Runtime & auxiliary options:
  --help, -h            Show this help message and exit
  --version, -V         show program's version number and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of used threads (default = number of available
                        CPU cores)
```

## Examples

Installation:

```bash
$ conda install -c bioconda referenceseeker
$ wget https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz
$ tar -xzf bacteria-refseq.tar.gz
$ rm bacteria-refseq.tar.gz
```

Simple:

```bash
$ # referenceseeker <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker bacteria-refseq/ genome.fasta
```

Expert: verbose output and increased output of candidate reference genomes using a defined number of threads:

```bash
$ # referenceseeker --crg 500 --verbose --threads 8 <REFERENCE_SEEKER_DB> <GENOME>
$ referenceseeker --crg 500 --verbose --threads 8 bacteria-refseq/ genome.fasta
```

## Databases

ReferenceSeeker depends on databases comprising taxonomic genome informations as well as kmer hash profiles for each entry.

### Pre-built

We provide pre-built databases based on public genome data hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4415843.svg)](https://doi.org/10.5281/zenodo.4415843) :

#### RefSeq

release: 205 (2021-04-01)

| Taxon | URL | # Genomes | Size |
| :---: | --- | ---: | :---: |
| bacteria | <https://zenodo.org/record/4415843/files/bacteria-refseq.tar.gz> | 30,941 | 40 Gb |
| archaea | <https://zenodo.org/record/4415843/files/archaea-refseq.tar.gz> | 606 | 553 Mb |
| fungi | <https://zenodo.org/record/4415843/files/fungi-refseq.tar.gz> | 347 | 3.3 Gb |
| protozoa | <https://zenodo.org/record/4415843/files/protozoa-refseq.tar.gz> | 88 | 1.1 Gb |
| viruses | <https://zenodo.org/record/4415843/files/viral-refseq.tar.gz> | 10,339 | 730 Mb |

#### GTDB

release: v95 (2021-01-06)

| Taxon | URL | # Genomes | Size |
| :---: | --- | ---: | :---: |
| bacteria | <https://zenodo.org/record/4415843/files/bacteria-gtdb.tar.gz> | 30,238 | 34 Gb |
| archaea | <https://zenodo.org/record/4415843/files/archaea-gtdb.tar.gz> | 1,672 | 1.1 Gb |

#### Plasmids

In addition to the genome based databases, we provide the following plasmid databases based on RefSeq and PLSDB:

| DB | URL | # Plasmids | Size |
| :---: | --- | ---: | :---: |
| RefSeq | <https://zenodo.org/record/4415843/files/plasmids-refseq.tar.gz> | 32,611 | 1.1 Gb |
| PLSDB | <https://zenodo.org/record/4415843/files/plasmids-plsdb.tar.gz> | 27,393 | 1.1 Gb |

### Custom database

If above mentiond RefSeq based databases do not contain sufficiently-close related genomes or are just too large, ReferenceSeeker provides auxiliary commands in order to either create databases from scratch or to expand existing ones. Therefore, a second executable `referenceseeker_db` accepts `init` and `import` subcommands:

Usage:

```bash
usage: referenceseeker_db [--help] [--version] {init,import} ...

Rapid determination of appropriate reference genomes.

positional arguments:
  {init,import}  sub-command help
    init         Initialize a new database
    import       Add a new genome to database

Runtime & auxiliary options:
  --help, -h     Show this help message and exit
  --version, -V  show program's version number and exit
```

If a new database should be created, use `referenceseeker_db init`:

```bash
usage: referenceseeker_db init [-h] [--output OUTPUT] --db DB

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        output directory (default = current working directory)
  --db DB, -d DB        Name of the new ReferenceSeeker database
```

This new database or an existing one can be used to import genomes in Fasta, GenBank or EMBL format:

```bash
usage: referenceseeker_db import [-h] --db DB --genome GENOME [--id ID]
                                 [--taxonomy TAXONOMY]
                                 [--status {complete,chromosome,scaffold,contig}]
                                 [--organism ORGANISM]

optional arguments:
  -h, --help            show this help message and exit
  --db DB, -d DB        ReferenceSeeker database path
  --genome GENOME, -g GENOME
                        Genome path [Fasta, GenBank, EMBL]
  --id ID, -i ID        Unique genome identifier (default sequence id of first
                        record)
  --taxonomy TAXONOMY, -t TAXONOMY
                        Taxonomy ID (default = 12908 [unclassified sequences])
  --status {complete,chromosome,scaffold,contig}, -s {complete,chromosome,scaffold,contig}
                        Assembly level (default = contig)
  --organism ORGANISM, -o ORGANISM
                        Organism name (default = "NA")
```

## Dependencies

ReferenceSeeker needs the following dependencies:

- Python (3.8, 3.9), Biopython (>=1.78), xopen(>=1.1.0)
- Mash (2.3) <https://github.com/marbl/Mash>
- MUMmer (4.0.0-beta2) <https://github.com/gmarcais/mummer>

ReferenceSeeker has been tested against aforementioned versions.

## Citation

> Schwengers et al., (2020). ReferenceSeeker: rapid determination of appropriate reference genomes. Journal of Open Source Software, 5(46), 1994, https://doi.org/10.21105/joss.01994

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
oliver.schwengers@cb.jlug.de.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
[https://www.contributor-covenant.org/version/2/0/code_of_conduct.html][v2.0].

Community Impact Guidelines were inspired by 
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available 
at [https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.0]: https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# Contribution Guidelines

## Setting up user information

Please, have **Git** set up with consistent user information before commiting. Preferably, provide your real name and a working email address.

**Example**:
```bash
git config --global user.name "Your Name Comes Here"
git config --global user.email you@yourdomain.example.com
```

## Make sure that your branch contains clean commits
- Follow the common sense guidelines for writing good commit messages (see below).
- Make separate commits for separate changes. If you cannot describe what the commit does in one sentence, it is probably a mix of changes and should be separated into several commits.
- Do not merge `master` into your branch. Instead use `git rebase`. If you need to resolve merge conflicts or include the latest changes.

## Check your coding style
- Make sure your contributions and changes follow the coding and indentation style of the code surrounding your changes.
- Do not commit commented-out code or files that are no longer needed. Remove the code or the files unless there is a good reason to keep it.

## Guidelines for good commit messages
1. Separate subject from body with a blank line
2. Use the imperative mood in the subject line ("Fix", "Add", "Change" instead of "Fixed", "Added", "Changed")
3. Limit the subject line to 50 characters
4. Reference an issue at the end of a subject line
5. Do not end the subject line with a period
6. Wrap the body at 72 characters
7. Use the body to explain what and why vs. how

Bad:
<pre>
some changes fixing this and that...
</pre>

Good:
<pre>
fix broken link in reports #123

As foo changed their internal data structure in the last release bar, we need to update our external links accordingly.
</pre>

# Helpful Sources
- https://www.atlassian.com/git/tutorials
---
title: 'ReferenceSeeker: rapid determination of appropriate reference genomes'
tags:
  - Python
  - Bioinformatics
  - WGS
  - NGS
  - Microbiology
authors:
  - name: Oliver Schwengers
    orcid: 0000-0003-4216-2721
    affiliation: "1, 2, 3"
  - name: Torsten Hain
    affiliation: "2, 3"
  - name: Trinad Chakraborty
    affiliation: "2, 3"
  - name: Alexander Goesmann
    orcid: 0000-0002-7086-2568
    affiliation: "1, 3"
affiliations:
 - name: Bioinformatics and Systems Biology, Justus Liebig University Giessen, Giessen, 35392, Germany
   index: 1
 - name: Institute of Medical Microbiology, Justus Liebig University Giessen, Giessen, 35392, Germany
   index: 2
 - name: German Centre for Infection Research (DZIF), partner site Giessen-Marburg-Langen, Giessen, Germany
   index: 3
date: 18 December 2019
bibliography: paper.bib
---

# Summary
The enormous success and ubiquitous application of next and third generation sequencing has
led to a large number of available high-quality draft and complete microbial genomes in the
public databases. Today, the NCBI RefSeq database release 90 alone contains 11,060
complete bacterial genomes ​[@Haft:2018​]. Concurrently, selection of appropriate reference
genomes (RGs) is increasingly important as it has enormous implications for routine in-silico
analyses, as for example in detection of single nucleotide polymorphisms, scaffolding of draft
assemblies, comparative genomics and metagenomic tasks. Therefore, a rigorously selected
RG is a prerequisite for the accurate and successful application of the aforementioned
bioinformatic analyses. In order to address this issue several new databases, methods and tools
have been published in recent years *e.g.* RefSeq, DNA-DNA hybridization [@Meier-Kolthoff:2013]​,
average nucleotide identity (ANI) as well as percentage of conserved DNA (conDNA) values [@Goris:2007] and Mash ​[@Ondov:2016]​.
Nevertheless, the sheer amount of currently available databases and potential RGs
contained therein, together with the plethora of tools available, often requires manual selection of
the most suitable RGs. To the best of the authors’ knowledge, there is currently no such tool
providing both an integrated, highly specific workflow and scalable and rapid implementation.
ReferenceSeeker was designed to overcome this bottleneck. As a novel command line tool, it
combines a fast kmer profile-based lookup of candidate reference genomes (CRGs) from high
quality databases with rapid computation of (mutual) highly specific ANI and conserved DNA values.

# Implementation
ReferenceSeeker is a linux command line tool implemented in Python 3. All necessary external
binaries are bundled with the software. The tool itself requires no external dependencies other
than Biopython for file input and output.

## Databases
ReferenceSeeker takes advantage of taxon-specific custom databases in order to reduce data
size and overall runtime. Pre-built databases for the taxonomic groups bacteria, archaea, fungi,
protozoa and viruses are provided. Each database integrates genomic as well as taxonomic
information comprising genome sequences of all RefSeq genomes with an assembly level
‘complete’ or whose RefSeq category is either denoted as ‘reference genome’ or ‘representative
genome’, as well as kmer profiles, related species names, NCBI Taxonomy identifiers and
RefSeq assembly identifiers. For convenient and fully automatic updates, we provide locally
executable scripts implemented in bash and Nextflow ​[@Di_Tommaso:2017]​. Non-public genomes can
be imported into existing or newly created databases by an auxiliary command line interface.

## Database Lookup of CRGs
To reduce the number of necessary ANI calculations a kmer profile-based lookup of CRGs
against custom databases is carried out. This step is implemented via Mash parameterized with
a Mash distance of 0.1, which was shown to correlate well with an ANI of roughly 90% ​
[@Ondov:2016] and thereby establishing a lower limit for reasonably related genomes.
The resulting set of CRGs is subsequently reduced to a configurable number of CRGs (default=100)
with the lowest Mash distances.

## Determination of RG
Mash distances used for the preliminary selection of CRGs were shown to correlate
well with ANI values capturing nucleotide-level sequence similarities. However,
Mash distances do not correlate well with the conDNA statistic, which captures the
query sequence coverage within a certain reference sequence (Figure 1).
In order to precisely calculate sequence similarities beyond the capability of kmer fingerprints
and to assure that RGs share an adequate portion of the query genome,
ReferenceSeeker calculates both ANI and conDNA to derive a highly specific measure
of microbial genome relationships [@Goris:2007].

![Figure 1: Scatter plots showing the correlation between Mash distance, ANI and conDNA
values. ANI and conserved DNA values are plotted against Mash distance values for 500
candidate reference genomes with the lowest Mash distance within the bacterial database for
10 randomly selected *Escherichia coli* genomes from the RefSeq database, each.](mash-ani-cdna.scatter.png)

Therefore, required sequence alignments are conducted via Nucmer of the
MUMmer package [@Marcais:2018] as it was recently shown that Nucmer based implementations (ANIn)
compare favourably to BLAST+ based implementations (ANIb) in terms of runtime.
Exact calculations of ANI and conDNA values were adopted from [@Goris:2007] and are carried out as follows.
Each query genome is split into consecutive 1,020 bp nucleotide fragments which are aligned to a reference genome via Nucmer.
The conDNA value is then calculated as the ratio between the sum of all aligned nucleotides within nucleotide fragments with an alignment with a sequence identity above 90% and the sum of nucleotides of all nucleotide fragments.
The ANI value is calculated as the mean sequence identity of all nucleotide fragments
with a sequence identity above 30% and an alignment length of at least 70% along the entire fragment length.

Given that compared genomes are closely related, *i.e.* they share an ANI of above 90%, it was also shown
that ANIn correlates well with ANIb [@Yoon:2017]. This requirement is ensured by the prior Mash-based
selection of CRGs. As an established threshold for species boundaries [@Goris:2007]​,
results are subsequently filtered by configurable ANI and conDNA values with a default of 95%
and 69%, respectively. Finally, CRGs are sorted according to the harmonic mean of ANI and
conDNA values in order to incorporate both the nucleotide identity and the genome coverage
between the query genome and resulting CRGs. In this manner, ReferenceSeeker ensures that the resulting RGs
sufficiently reflect the genomic landscape of a query genome. If desired by the user, this approach
can be extended to a bidirectional computation of aforementioned ANI and conDNA values.

# Application
ReferenceSeeker takes as input a microbial genome assembly in fasta format and the path to a
taxonomic database of choice. Results are returned as a tabular separated list comprising the
following information: RefSeq assembly identifier, ANI, conDNA, NCBI taxonomy identifier,
assembly status and organism name.
To illustrate the broad applicability at different scales we tested ReferenceSeeker with 12
microbial genomes from different taxonomic groups and measured overall runtimes on a
common consumer laptop providing 4 cores and a server providing 64 cores (Table \ref{Table 1}). For the
tested bacterial genomes, ReferenceSeeker limited the number of resulting RGs to a default
maximum of 100 genomes. Runtimes of archaeal and viral genomes are significantly shorter
due to a small number of available RGs in the database and overall smaller genome sizes,
respectively.

Table: Runtimes and numbers of resulting RG executed locally on a quad-core moderate consumer
laptop and a 64 core server machine. \label{Table 1}

------------------------------------------------------------------------------------
Genome                                           Genome Size   Laptop   Server  # RG
                                                        [kb]  [mm:ss]  [mm:ss]      
----------------------------------------------- ------------ -------- -------- -----
*Escherichia coli* str. K-12 substr. MG1665            4,641     3:24     0:30  100*
(GCF_000005845.2)                                                                   

*Pseudomonas aeruginosa* PAO1                          6,264     5:20     0:44  100*
(GCF_000006765.1)                                                                   

*Listeria monocytogenes* ​EGD-e                         2,944     2:52     0:24  100*
(GCF_000196035.1)                                                                   

*Staphylococcus aureus* subsp aureus NCTC 8325         2,821     2:31     0:21  100*
(GCF_000013425.1)                                                                   

*Halobacterium salinarum* NRC-1                        2,571     0:04     0:03     2
(GCF_000006805.1)                                                                   

*Methanococcus maripaludis* X1                         1,746     0:22     0:09     5
(GCF_000220645.1)                                                                   

*Aspergillus fumigatus* ​Af293                         29,384     3:11     2:07     1
(GCF_000002655.1)                                                                   

*Candida albicans* SC5314                             14,282     0:21     0:19     1
(GCF_000182965.3)                                                                   

*Entamoeba histolytica* HM-1:IMSS                     20,835     6:04     4:41     1
(GCF_000208925.1)                                                                   

*Plasmodium falciparum* ​3D7                           23,326     2:52     1:49     1
(GCF_000002765.4)                                                                   

*Influenza A* virus                                       13     0:03     0:02     1
(GCF_001343785.1)                                                                   

*Human coronavirus* NL63                                  27     0:04     0:02     1
(GCF_000853865.1)                                                                   
------------------------------------------------------------------------------------

# Availability
The source code is available on GitHub under a GPL3 license: https://github.com/oschwengers/referenceseeker​.
The software is packaged and publicly available via BioConda. Pre-built databases for
bacteria, archaea, fungi, protozoa and viruses are hosted at Zenodo: https://doi.org/10.5281/zenodo.3562005​.
All installation instructions, examples and download links are provided on GitHub.

# Funding
This work was supported by the German Center of Infection Research (DZIF) [DZIF grant 8000
 701–3 (HZI), TI06.001 and 8032808811 to T.C.]; the German Network for Bioinformatics
Infrastructure (de.NBI) [BMBF grant FKZ 031A533B to A.G.]; and the German Research
Foundation (DFG) [SFB-TR84 project A04 (TRR84/3 2018) to T.C., KFO309 Z1 (GO 2037/5-1)
to A.G., SFB-TR84 project B08 (TRR84/3 2018) to T.H., SFB1021 Z02 (SFB 1021/2 2017) to
T.H., KFO309 Z1 (HA 5225/1-1) to T.H.].

Authors declare that there are no conflicts of interest.

# Acknowledgement
The authors thank Karina Brinkrolf for valuable discussions, testing and bug reports.

# References
---
name: Bug report
about: Create a report to help us improving the ReferenceSeeker
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
To help users and fix any bugs and issues a concise description if beneficial and often necessary.

Therefore, please provide us with at least the following information:
- what exactly happened
- what exact command was executed: just copy-paste the command line
- what installation of ReferenceSeeker did you use: `BioConda`, `GitHub`, `Pip`
- what database was used: `public`, `custom`, which taxon
- which version of ReferenceSeeker was used
---
name: Feature request
about: Any suggestions and new ideas for this project are highly welcome!
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to an existing issue or bug?**
Briefly reference any existing issues this is related to.

**Is your new feature related to a general problem?**
Please provide us with a clear and concise description of what the problem is:
- what happened which was *not* expected
- what did *not* happen which was expected
- what would you like to happen or see in the output

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or information about your feature request here.
