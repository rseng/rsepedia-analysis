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
[![DOI:10.1099/mgen.0.000398](https://zenodo.org/badge/DOI/10.1099/mgen.0.000398.svg)](https://doi.org/10.1099/mgen.0.000398)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/oschwengers/platon/blob/master/LICENSE)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/cb-platon.svg)
![GitHub release](https://img.shields.io/github/release/oschwengers/platon.svg)
[![PyPI](https://img.shields.io/pypi/v/cb-platon.svg)](https://pypi.org/project/cb-platon)
![PyPI - Status](https://img.shields.io/pypi/status/cb-platon.svg)
[![Conda](https://img.shields.io/conda/v/bioconda/platon.svg)](https://bioconda.github.io/recipes/platon/README.html)

# Platon: identification and characterization of bacterial plasmid contigs from short-read draft assemblies

## Contents

- [Description](#description)
- [Input/Output](#inputoutput)
- [Installation](#installation)
  - [BioConda](#bioconda)
  - [GitHub/Pip](#githubpip)
- [Usage](#usage)
- [Examples](#examples)
- [Mode](#mode)
- [Database](#database)
- [Dependencies](#dependencies)
- [Citation](#citation)
- [Issues](#issues)

## Description

**TL;DR**
Platon detects plasmid-borne contigs within bacterial draft (meta) genomes assemblies. Therefore, Platon analyzes the distribution bias of protein-coding gene families among chromosomes and plasmids. This analysis is complemented by comprehensive contig characterizations follwoed by heuristic filters.

Platon conducts three analysis steps:

1. It predicts and searches protein sequences against a custom and pre-computed database comprising marker protein sequences (**MPS**) and related replicon distribution scores (**RDS**). These scores express the empirically measured bias of protein sequence family distributions among plasmids and chromosomes pre-computed on complete NCBI RefSeq replicons. Platon calculates the mean RDS for each contig and either classifies them as chromosome if the RDS is below a sensitivity cutoff determined to 95% sensitivity or as plasmid if the RDS is above a specificity cutoff determined to 99.9% specificity. Exact values for these thresholds have been computed based on Monte Carlo simulations of artifical replicon fragments created from complete RefSeq chromosome and plasmid sequences.
2. Contigs passing the sensitivity filter get comprehensivley characterized. Hereby, Platon tries to circularize the contig sequences, searches for rRNA, replication, mobilization and conjugation genes, oriT sequences, incompatibility group DNA probes and finally performs a BLAST+ search against the NCBI plasmid database.
3. Finally, to increase the overall sensitivity, Platon classifies all remaining contigs based on the gathered information by several heuristics.

| ![Replicon distribution and alignment hit frequencies of MPS](rds-ratio-counts.web.png?raw=true) |
| -- |
| *Fig: Replicon distribution and alignment hit frequencies of MPS. Shown are summed plasmid and chromosome alignment hit frequencies per MPS plotted against plasmid/chromosome hit count ratios scaled to [-1 (chromosome), 1 (plasmid)]; Hue: normalized RDS values (min=-100, max=100), hit count outliers below 10-4 and above 1 are discarded for the sake of readability.* |

## Input/Output

### Input

Platon accepts draft (meta) genome assemblies in fasta format. If contigs have been assembled with SPAdes, Platon is able to extract the coverage information from the contig names.

### Output

For each contig classified as plasmid sequence the following columns are printed to `STDOUT` as tab separated values:

- Contig ID
- Length
- Coverage
- \# ORFs
- RDS
- Circularity
- Incompatibility Type(s)
- \# Replication Genes
- \# Mobilization Genes
- \# OriT Sequences
- \# Conjugation Genes
- \# rRNA Genes
- \# Plasmid Database Hits

In addition, Platon writes the following files into the output directory:

- `<prefix>`.plasmid.fasta: contigs classified as plasmids or plasmodal origin
- `<prefix>`.chromosome.fasta: contigs classified as chromosomal origin
- `<prefix>`.tsv: dense information as printed to STDOUT (see above)
- `<prefix>`.json: comprehensive results and information on each single plasmid contig.
All files are prefixed (`<prefix>`) as the input genome fasta file.

## Installation

Platon can be installed via BioConda or Pip.
However, we encourage to use [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to automatically install all required 3rd party dependencies. In all cases a mandatory [database](#database_download) must be downloaded.

### BioConda

```bash
$ conda install -c conda-forge -c bioconda -c defaults platon
```

### Pip

```bash
$ python3 -m pip install --user cb-platon
```

Platon requires the following 3rd party executables which must be installed & executable:

- Prodigal (2.6.3) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648> <https://github.com/hyattpd/Prodigal>
- Diamond (2.0.6) <https://pubmed.ncbi.nlm.nih.gov/25402007> <http://www.diamondsearch.org>
- Blast+ (2.10.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
- MUMmer (4.0.0-beta2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC395750/> <https://github.com/gmarcais/mummer>
- HMMER (3.3.1) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3695513/> <http://hmmer.org/>
- INFERNAL (1.1.4) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854> <http://eddylab.org/infernal>

### Database download

Platon requires a mandatory database which is publicly hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4066768.svg)](https://doi.org/10.5281/zenodo.4066768)
Further information is provided in the [database](#database) section below.

```bash
$ wget https://zenodo.org/record/4066768/files/db.tar.gz
$ tar -xzf db.tar.gz
$ rm db.tar.gz
```

The db path can either be provided via parameter (`--db`) or environment variable (`PLATON_DB`):

```bash
$ platon --db <db-path> genome.fasta

$ export PLATON_DB=<db-path>
$ platon genome.fasta
```

Additionally, for a system-wide setup, the database can be copied to the Platon base directory:

```bash
$ cp -r db/ <platon-installation-dir>
```

## Usage

Usage:

```bash
usage: platon [--db DB] [--prefix PREFIX] [--output OUTPUT] [--mode {sensitivity,accuracy,specificity}] [--characterize] [--meta] [--help] [--verbose] [--threads THREADS] [--version] <genome>

Identification and characterization of bacterial plasmid contigs from short-read draft assemblies.

Input / Output:
  <genome>              draft genome in fasta format
  --db DB, -d DB        database path (default = <platon_path>/db)
  --prefix PREFIX, -p PREFIX
                        Prefix for output files
  --output OUTPUT, -o OUTPUT
                        Output directory (default = current working directory)

Workflow:
  --mode {sensitivity,accuracy,specificity}, -m {sensitivity,accuracy,specificity}
                        applied filter mode: sensitivity: RDS only (>= 95% sensitivity); specificity: RDS only (>=99.9% specificity); accuracy: RDS & characterization heuristics (highest accuracy) (default = accuracy)
  --characterize, -c    deactivate filters; characterize all contigs
  --meta                use metagenome gene prediction mode

General:
  --help, -h            Show this help message and exit
  --verbose, -v         Print verbose information
  --threads THREADS, -t THREADS
                        Number of threads to use (default = number of available CPUs)
  --version             show program's version number and exit
```

## Examples

Simple:

```bash
$ platon genome.fasta
```

Expert: writing results to `results` directory with verbose output using 8 threads:

```bash
$ platon --db ~/db --output results/ --verbose --threads 8 genome.fasta
```

## Mode

Platon provides 3 different modi controlling which filters will be used.
`Accuracy` mode is the preset default.

### Sensitivity

In the `sensitivity` mode Platon will classifiy all contigs with an `RDS` value *below* the sensitivity threshold as chromosomal and all remaining contigs as plasmid. This threshold was defined to account for 95% sensitivity and computed via Monte Carlo simulations of artifical contigs resulting in an RDS=-7.9.
-> use this mode to *exclude chromosomal* contigs.

### Specificity

In the `specificity` mode Platon will classifiy all contigs with an `RDS` value *above* the specificity threshold as plasmid and all remaining contigs as chromosomal. This threshold was defined to account for 99.9% specificity and computed via Monte Carlo simulations of artifical contigs resulting in an RDS=0.7.

### Accuracy (default)

In the `accuracy` mode Platon will classifiy all contigs with:

- an `RDS` value *below* the sensitivity threshold as chromosomal
- an `RDS` value *above* the specificity threshold as plasmid and in addition all contigs as plasmid for which one of the following is true: it
- can be circularized
- has an incompatibility group sequence
- has a replication or mobilization HMM hit
- has an oriT hit
- has an RDS above the conservative score (0.1), a RefSeq plasmid hit and *no* rRNA hit

## Database

Platon depends on a custom database based on MPS, RDS, RefSeq Plasmid database, PlasmidFinder db as well as manually curated MOB HMM models from MOBscan, custom conjugation and replication HMM models and oriT sequences from MOB-suite. This database based on UniProt UniRef90 release 202 can be downloaded here: (zipped 1.6 Gb, unzipped 2.8 Gb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4066768.svg)](https://doi.org/10.5281/zenodo.4066768)
[https://zenodo.org/record/4066768/files/db.tar.gz](https://zenodo.org/record/4066768/files/db.tar.gz)

*Please make sure that you use the latest Platon version along with the most recent database version! Older software versions are **not** compatible with the latest database version*

## Dependencies

Platon was developed and tested in Python 3.5 and depends on BioPython (>=1.71).

Additionally, it depends on the following 3rd party executables:

- Prodigal (2.6.3) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648> <https://github.com/hyattpd/Prodigal>
- Diamond (2.0.6) <https://pubmed.ncbi.nlm.nih.gov/25402007> <http://www.diamondsearch.org>
- Blast+ (2.10.1) <https://www.ncbi.nlm.nih.gov/pubmed/2231712> <https://blast.ncbi.nlm.nih.gov>
- MUMmer (4.0.0-beta2) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC395750/> <https://github.com/gmarcais/mummer>
- HMMER (3.3.1) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3695513/> <http://hmmer.org/>
- INFERNAL (1.1.4) <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3810854> <http://eddylab.org/infernal>

## Citation

> Schwengers O., Barth P., Falgenhauer L., Hain T., Chakraborty T., & Goesmann A. (2020). Platon: identification and characterization of bacterial plasmid contigs in short-read draft assemblies exploiting protein sequence-based replicon distribution scores. Microbial Genomics, 95, 295. https://doi.org/10.1099/mgen.0.000398

As Platon takes advantage of the inc groups, MOB HMMs and oriT sequences of the following databases, please also cite:

- > Carattoli A., Zankari E., Garcia-Fernandez A., Voldby Larsen M., Lund O., Villa L., Aarestrup F.M., Hasman H. (2014) PlasmidFinder and pMLST: in silico detection and typing of plasmids. Antimicrobial Agents and Chemotherapy, https://doi.org/10.1128/AAC.02412-14

- > Garcillán-Barcia M. P., Redondo-Salvo S., Vielva L., de la Cruz F. (2020) MOBscan: Automated Annotation of MOB Relaxases. Methods in Molecular Biology, https://doi.org/10.1007/978-1-4939-9877-7_21

- > Robertson J., Nash J. H. E. (2018) MOB-suite: Software Tools for Clustering, Reconstruction and Typing of Plasmids From Draft Assemblies. Microbial Genomics, https://doi.org/10.1099/mgen.0.000206

## Issues

If you run into any issues with Platon, we'd be happy to hear about it! Please, start the pipeline with `-v` (verbose) and do not hesitate to file an issue including as much of the following as possible:

- a detailed description of the issue
- the platon cmd line output
- the `<prefix>.json` file if possible
- A reproducible example of the issue with a small dataset that you can share (helps us identify whether the issue is specific to a particular computer, operating system, and/or dataset).

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
---
name: Bug report
about: Create a report to help us improving Platon
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
To help users and fix any bugs and issues a concise description if beneficial and often necessary.

Therefore, please provide us with at least the following information:
- what exactly happened
- what exact command was executed: just copy-paste the command line
- what installation of Platon did you use: `BioConda`, `GitHub`, `Pip`
- which version of Platon was used
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
