# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at [INSERT EMAIL ADDRESS]. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Opfi

[![Documentation Status](https://readthedocs.org/projects/opfi/badge/?version=latest)](https://opfi.readthedocs.io/en/latest/?badge=latest)
[![PyPI](http://img.shields.io/pypi/v/opfi.svg)](https://pypi.python.org/pypi/opfi/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/opfi/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

A python package for discovery, annotation, and analysis of gene clusters in genomics or metagenomics datasets.

## Installation

The recommended way to install Opfi is with [Bioconda](https://bioconda.github.io/), which requires the [conda](https://docs.conda.io/en/latest/) package manager. This will install Opfi and all of its dependencies (which you can read more about [here](https://opfi.readthedocs.io/en/latest/installation.html)).

Currently, Bioconda supports only 64-bit Linux and Mac OS. Windows users can still install Opfi with pip (see below); however, the complete installation procedure has not been fully tested on a Windows system. 

### Install with conda (Linux and Mac OS only)

First, set up conda and Bioconda following the [quickstart](https://bioconda.github.io/user/install.html) guide. Once this is done, run:

```
conda install -c bioconda opfi
```

And that's it! Note that this will install Opfi in the conda environment that is currently active. To create a fresh environment with Opfi installed, do:

```
conda create --name opfi-env -c bioconda opfi
conda activate opfi-env
```

### Install with pip 

This method does not automatically install non-Python dependencies, so they will need to be installed separately, following their individual installation instructions. A complete list of required software is available [here](https://opfi.readthedocs.io/en/latest/installation.html#dependencies). Once this step is complete, install Opfi with pip by running:

```
pip install opfi
```

For information about installing for development, check out the [documentation site](https://opfi.readthedocs.io/en/latest/installation.html).

## Gene Finder

Gene Finder iteratively executes homology searches to identify gene clusters of interest. Below is an example script that sets up a search for putative CRISPR-Cas systems in the Rippkaea orientalis PCC 8802 (cyanobacteria) genome. Data inputs are provided in the Opfi tutorial (`tutorials/tutorial.ipynb`).

```python
from gene_finder.pipeline import Pipeline
import os

genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"

p = Pipeline()
p.add_seed_step(db="cas1", name="cas1", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_filter_step(db="cas_all", name="cas_all", e_val=0.001, blast_type="PROT", num_threads=1)
p.add_crispr_step()

# use the input filename as the job id
job_id = os.path.basename(genomic_data)
results = p.run(job_id=job_id, data=genomic_data, min_prot_len=90, span=10000, gzip=True)
```

# Operon Analyzer

Operon Analyzer filters results from Gene Finder, and identifies promising candidate operons according to a given set of criteria. It also contains some tools for visualizing candidates and performing basic statistics.

Please note that the use of the word "operon" throughout this library is somewhat of an artifact from early development. At this time, Opfi does not predict whether a candidate system represents a true operon, that is, a set of genes under the control of a single promoter. Although a candidate gene cluster may certainly qualify as an operon, it is currently up to the user to make that distinction. 

## Analysis

The analysis module provides tools to identify operons that conform to certain rules, such as requiring that they contain a certain gene, or that two genes are within a given distance of each other (the full list is given below). CSV output is written to stdout, which identifies the outcome of the analysis for each putative operon.

Rules defined with the `RuleSet` determine whether an operon should be considered a candidate for further analysis. 
Filters defined with the `FilterSet` help define which features to consider when evaluating rules. You might, for example, want to exclude any operon containing a particular gene, but if a copy of that gene coincidentally exists 5 kb from the true operon, you might want to ignore it for the purposes of evaluating your rules. 

A sample script that performs this task is given here:

```python
import sys
from operon_analyzer.analyze import analyze
from operon_analyzer.rules import RuleSet, FilterSet


rs = RuleSet().require('transposase') \
              .exclude('cas3') \
              .at_most_n_bp_from_anything('transposase', 500) \
              .same_orientation()

fs = FilterSet().pick_overlapping_features_by_bit_score(0.9) \
                .must_be_within_n_bp_of_anything(1000)

if __name__ == '__main__':
    analyze(sys.stdin, rs, fs)
```

### List of available rules

  * `exclude(feature_name: str)`: Forbid the presence of a particular feature. 
  * `require(feature_name: str)`: Require the presence of a particular feature. 
  * `max_distance(feature1_name: str, feature2_name: str, distance_bp: int)`: The two given features must be no further than `distance_bp` base pairs apart. Requires exactly one of each feature to be present.
  * `at_least_n_bp_from_anything(feature_name: str, distance_bp: int)`: Requires that a feature be at least `distance_bp` base pairs away from any other feature.  This is mostly useful for eliminating overlapping features.
  * `at_most_n_bp_from_anything(feature_name: str, distance_bp: int)`: A given feature must be within `distance_bp` base pairs of another feature. Requires exactly one matching feature to be present. Returns `False` if the given feature is the only feature.
  * `same_orientation(exceptions: Optional[List[str]] = None)`: All features in the operon must have the same orientation.
  * `contains_any_set_of_features(sets: List[List[str]])`: Returns `True` if the operon contains features with all of the names in at least one of the lists. Useful for determining if an operon contains all of the essential genes for a particular system, for example.
  * `contains_exactly_one_of(feature1_name: str, feature2_name: str)`: An exclusive-or of the presence of two features.  That is, one of the features must be present and the other must not.
  * `contains_at_least_n_features(feature_names: List[str], feature_count: int, count_multiple_copies: bool = False)`: The operon must contain at least `feature_count` features in the list. By default, a matching feature that appears multiple times in the operon will only be counted once; to count multiple copies of the same feature, set `count_multiple_copies=True`.
  * `contains_group(self, feature_names: List[str], max_gap_distance_bp: int, require_same_orientation: bool)`: The operon must contain a contiguous set of features (in any order) separated by no more than max_gap_distance_bp. Optionally, the user may require that the features must all have the same orientation.
  * `maximum_size(self, feature_name: str, max_bp: int, all_matching_features_must_pass: bool = False, regex: bool = False)`: The operon must contain at least one feature with feature_name with a size (in base pairs) of max_bp or smaller. If all_matching_features_must_pass is True, every matching Feature must be at least max_bp long.
  * `minimum_size(self, feature_name: str, min_bp: int, all_matching_features_must_pass: bool = False, regex: bool = False)`: The operon must contain at least one feature with feature_name with a size (in base pairs) of min_bp or larger. If all_matching_features_must_pass is True, every matching Feature must be at least min_bp long. 
  * `custom(rule: 'Rule')`: Add a rule with a user-defined function. 

### List of available filters

  * `must_be_within_n_bp_of_anything(distance_bp: int)`: If a feature is very far away from anything it's probably not part of an operon.
  * `must_be_within_n_bp_of_feature(feature_name: str, distance_bp: int)`: There may be situations where two features always appear near each other in functional operons.  
  * `pick_overlapping_features_by_bit_score(minimum_overlap_threshold: float)`: If two features overlap by more than `minimum_overlap_threshold`, the one with the lower bit score is ignored.
  * `custom(filt: 'Filter')`: Add a filter with a user-defined function. 

### Analysis Output 

Each line of the CSV will contain an accession ID and the path to the file that contains it, the contig coordinates, and whether it passed or failed the given rules. If it passed, the last column will contain the word `pass` only. Otherwise it will start with `fail` followed by a comma-delimited list of the serialized rules that it failed to adhere to (with the name and parameters that were passed to the method).

## Visualization

Interesting operons can be visualized with a simple gene diagram. It is up to the user to decide how to define this, though this sample script below creates diagrams for all operons that passed all given rules:

```python
import sys
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.visualize import build_operon_dictionary, plot_operons

analysis_csv, pipeline_csv, image_directory = sys.argv[1:]
good_operons = []

with open(pipeline_csv) as f:
    operons = build_operon_dictionary(f)
with open(analysis_csv) as f:
    for contig, filename, start, end, result in load_analyzed_operons(f):
        if result[0] != 'pass':
            continue
        op = operons.get((contig, filename, start, end))
        if op is None:
            continue
        good_operons.append(op)
plot_operons(good_operons, image_directory)
```

## Overview Statistics

Some basic tools are provided to inspect the nature of operons that did not pass all given rules. The intent here is to help researchers determine if their filtering is too aggressive (or not aggressive enough), and to get an overall better feel for the data.

Simple bar plots can be produced as follows:

```python
import sys
import matplotlib.pyplot as plt
from operon_analyzer.analyze import load_analyzed_operons
from operon_analyzer.overview import load_counts


def plot_bar_chart(filename, title, data, rotate=True):
    fig, ax = plt.subplots()
    x = [str(d[0]).replace(":", "\n") for d in data]
    y = [d[1] for d in data]
    ax.bar(x, y, edgecolor='k')
    if rotate:
        plt.xticks(rotation=90)
    ax.set_ylabel("Count")
    ax.set_title(title)
    plt.savefig("%s.png" % filename, bbox_inches='tight')
    plt.close()


if __name__ == '__main__':
    unique_rule_violated, failed_rule_occurrences, rule_failure_counts = load_counts(sys.stdin)
    plot_bar_chart("sole-failure.png", "Number of times that each rule\nwas the only one that failed", sorted(unique_rule_violated.items()))
    plot_bar_chart("total-failures", "Total number of rule failures", sorted(failed_rule_occurrences.items()))
    plot_bar_chart("failures-at-each-contig", "Number of rules failed at each contig", sorted(rule_failure_counts.items()), rotate=False)
```
# Changelog

## Opfi 0.1.2

Update version number for JOSS publication and archive.

## Opfi 0.1.1

Pypi badge added to documentation.

## Opfi 0.1.0

Initial release.
# Contributing
Thank you for your interest in contributing to Opfi. Contribututions can take many forms, including:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

## Issues
Issues should be used to report problems or bugs with the library, to request new features, or to discuss potential changes before PRs are created.

Good bug reports have:

- A summary of the issue
- Specific steps to reproduce the problem, preferably in the form of code examples
- A description of expected behavior
- A description of actual behavior

## Pull Requests
In general, we use [Github Flow](https://guides.github.com/introduction/flow/index.html). This means that all changes happen through Pull Requests:

1. Fork the repository to your own Github account
2. Clone the project to your machine
3. Create a branch locally with a succinct but descriptive name
4. Commit changes to the branch
5. Push changes to your fork
6. Open a PR in our repository and include a description of the changes and a reference to the issue the PR addresses

## Coding Style
Contributors should attempt to adhere to the [PEP8](https://www.python.org/dev/peps/pep-0008/) style guide when possible, although this is not currently strictly enforced. At a minimum, new contributions should follow a style that is consistent with the preexisting code base. 

## Testing
We use the [pytest](https://docs.pytest.org/en/6.2.x/) framework for creating unit tests. Please write tests for any new code contributed to this project.

## License
By contributing, you agree that your contributions will be licensed under its MIT License.

## Code of Conduct
We take our open source community seriously and hold ourselves and other contributors to high standards of communication. By participating and contributing to this project, you agree to uphold our [Code of Conduct](https://github.com/wilkelab/Opfi/blob/contributing-guide/CODE-OF-CONDUCT.md).

## Attribution
This document was adapted from the [General Contributing Guidelines](https://github.com/extendr/rextendr/blob/main/CONTRIBUTING.md) of the rextendr project, and from the [contributing template](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62) developed by [briandk](https://gist.github.com/briandk).
---
title: 'Opfi: A Python package for identifying gene clusters in large genomics and metagenomics data sets'
tags:
  - Python
  - bioinformatics
  - metagenomics
  - gene cluster analysis
authors:
  - name: Alexis M. Hill^[co-first author, corresponding author]
    affiliation: 1
  - name: James R. Rybarski^[co-first author]
    affiliation: 2
  - name: Kuang Hu
    affiliation: "1,2"
  - name: Ilya J. Finkelstein
    affiliation: "2,3"
  - name: Claus O. Wilke
    affiliation: 1
affiliations:
 - name: Department of Integrative Biology, The University of Texas at Austin, Austin, Texas 78712, USA
   index: 1
 - name: Department of Molecular Biosciences, The University of Texas at Austin, Austin, Texas 78712, USA
   index: 2
 - name: Center for Systems and Synthetic Biology, The University of Texas at Austin, Austin, Texas, 78712, USA
   index: 3
date: 28 July 2021
bibliography: paper.bib

---

# Summary

Gene clusters are sets of co-localized, often contiguous genes that together perform specific functions, many of which are relevant to biotechnology. There is a need for software tools that can extract candidate gene clusters from vast amounts of available genomic data. Therefore, we developed Opfi: a modular pipeline for identification of arbitrary gene clusters in assembled genomic or metagenomic sequences. Opfi contains functions for annotation, de-deduplication, and visualization of putative gene clusters. It utilizes a customizable rule-based filtering approach for selection of candidate systems that adhere to user-defined criteria. Opfi is implemented in Python, and is available on the Python Package Index and on Bioconda [@Grüning:2018].  

# Statement of need

Gene clusters have been successfully repurposed for a number of biotechnical applications, including biofuel production, organic compound synthesis, and gene editing [@Fischbach:2010]. Despite the broad utility of known gene clusters, identification of novel gene clusters remains a challenging task. While there are many tools available for annotation of singular genes (or protein domains) in biological sequence data [@Camacho:2009; @Steinegger:2017; @Buchfink:2021], these programs do not identify whole gene clusters out of the box. In many cases, researchers must combine bioinformatics tools ad hoc, resulting in one-off pipelines that can be difficult to reproduce. Several software packages have been developed for the discovery of specific types of gene clusters [@Blin:2019; @Santos-Aberturas:2019; @vanHeel:2018], but these tools may not be sufficiently flexible to identify clusters of an arbitrary genomic composition. To address these gaps, we developed a modular pipeline that integrates multiple bioinformatics tools, providing a flexible, uniform computational framework for identification of arbitrary gene clusters. In a recent study, we used Opfi to uncover novel CRISPR-associated transposons (CASTs) in a large metagenomics dataset [@Rybarski:2021].

# Implementation

Opfi is implemented in Python, and uses several bioinformatics tools for feature annotation [@Camacho:2009; @Steinegger:2017; @Buchfink:2021; @Edgar:2007; @Shi:2019]. Users can install Opfi and all of its dependencies through Bioconda [@Grüning:2018]. Opfi consists of two major components: Gene Finder, for discovery of gene clusters, and Operon Analyzer, for rule-based filtering, deduplication, and visualization of gene clusters identified by Gene Finder. All modules generate output in a comma-separated (CSV) format that is common to the entire package.

## Example Gene Finder usage

The following example script searches for putative CRISPR-Cas loci in the genome of *Rippkaea orientalis PCC 8802*. Information about the biological significance of this example, as well as data inputs and descriptions, can be found in the `tutorials` directory in the project GitHub repository. The example illustrates the use of the `Pipeline` class for setting up a gene cluster search. First, `add_seed_step` specifies a step to annotate *cas1* genes, using protein BLAST (BLASTP) [@Camacho:2009] and a database of representative Cas1 protein sequences. 10,000 bp regions directly up- and downstream of each putative *cas1* gene are selected for further analysis, and all other regions are discarded. Next, `add_filter_step` adds a step to annotate candidate regions for additonal *cas* genes. Candidates that do not have at least one additional *cas* gene are discarded from the master list of putative systems. Finally, `add_crispr_step` adds a step to search remaining candidates for CRISPR arrays, i.e. regions of alternating ~30 bp direct repeat and variable sequences, using the PILER-CR repeat finding software [@Edgar:2007]. 

```python
from gene_finder.pipeline import Pipeline
import os

genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"
job_id = "r_orientalis"

p = Pipeline()
p.add_seed_step(db="cas1", name="cas1", e_val=0.001, blast_type="PROT")
p.add_filter_step(db="cas_all", name="cas", e_val=0.001, blast_type="PROT")
p.add_crispr_step()

p.run(job_id=job_id, data=genomic_data, span=10000, gzip=True)
```

Running this code creates the CSV file `r_orientalis_results.csv`, which contains information about each system identified; in this example, that is two canonical CRISPR-Cas systems, and one locus with weak homology to *cas* genes. Each line in the file represents a single putative feature in a candidate locus. Features from the same candidate are grouped together in the CSV. Detailed information about the output format can be found in the Opfi [documentation](https://opfi.readthedocs.io/).

## Example Operon Analyzer usage

In the previous example, passing systems must meet the relatively permissive criterion of having at least one *cas1* gene co-localized with one additional *cas* gene. This is sufficient to identify CRISPR-Cas loci, but may also capture regions that do not contain functional CRISPR-Cas systems, but rather consist of open reading frames (ORFs) with weak homology to *cas* genes. These improbable systems could be eliminated during the homology search by making the match acceptance threshold more restrictive (i.e., by decreasing the e-value), however, this could result in the loss of interesting, highly diverged systems. Therefore, we implemented a module that enables post-homology search filtering of candidate systems, using flexible rules that can be combined to create sophisticated elimination functions. This allows the user to first perform a broad homology search with permissive parameters, and then apply rules to cull unlikely candidates without losing interesting and/or novel systems. Additionally, rules may be useful for selecting candidates with a specific genomic composition for downstream analysis. It should be noted that the use of the term "operon" throughout this library is an artifact from early development of Opfi. At this time, Opfi does not predict whether a candidate system represents a true operon, that is, a set of genes under the control of a single promoter. Although a candidate gene cluster may certainly qualify as an operon, it is currently up to the user to make that distinction. 

Rule-based filtering is illustrated with the following example. The sample script takes the output generated by the previous example and reconstructs each system as an `Operon` object. Next, the `RuleSet` class is used to assess each candidate; here, passing systems must contain two cascade genes (*cas5* and *cas7*) no more than 1000 bp apart, and at least one *cas3* (effector) gene. For a complete list of rules, see the Opfi [documentation](https://opfi.readthedocs.io/). 

```python
from operon_analyzer import analyze, rules

rs = rules.RuleSet()
rs.contains_group(["cas5", "cas7"], max_gap_distance_bp = 1000)
rs.require("cas3"))

with open("r_orientalis_results.csv", "r") as input_csv:
    with open("filtered_output.csv", "w") as output_csv:
        analyze.evaluate_rules_and_reserialize(input_csv, rs, output_csv)
```

After running this code, the file `filtered_output.csv` contains only high-confidence type-I CRISPR-Cas systems (re-serialized to CSV format) that passed all rules in the rule set. 

## Candidate visualization

Opfi integrates the `DNAFeaturesViewer` package [@Zulkower:2020] to create gene diagrams of candidate systems. Each input system is visualized as a single PNG image. The sample script below reads in output from the previous example, and generates two gene diagram images, one for each CRISPR-Cas system present in *Rippkaea orientalis*. One image is provided for reference in \autoref{fig:operon}. 

```python
from operon_analyzer import load, rules, visualize

feature_colors = { "cas1": "lightblue",
                    "cas2": "seagreen",
                    "cas3": "gold",
                    "cas4": "springgreen",
                    "cas5": "darkred",
                    "cas6": "thistle",
                    "cas7": "coral",
                    "cas8": "red",
                    "cas9": "palegreen",
                    "cas10": "blue",
                    "cas11": "tan",
                    "cas12": "orange",
                    "cas13": "saddlebrown",
                    "CRISPR array": "purple"
                    }

fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
with open("filtered_output.csv", "r") as operon_data:
    operons = [operon for operon in load.load_operons(operon_data)]
    for operon in operons:
      fs.evaluate(operon)
    visualize.plot_operons(operons, output_directory=".", \
      plot_ignored=False, feature_colors=feature_colors)
```

The `FilterSet` class is used to resolve features with sequences that overlap by more than 90%. Specifically, only the overlapping feature with the highest bitscore value (a quantity that describes the overall quality of an alignment) is rendered when `pick_overlapping_features_by_bit_score` is applied. Note that is not a requirement for candidate visualization, but can improve gene diagram clarity.

![One of two type-I CRISPR-Cas systems present in the genome of *Rippkaea orientalis PCC 8802*. Note that the ORF beginning at position ~2500 has homology with both *cas1* and *cas4*. These alignments have identical bitscores (i.e., the goodness of alignments is quivalent, using this metric), so both annotations appear in the diagram, even though `pick_overlapping_features_by_bit_score` was applied.\label{fig:operon}](operon_diagram.png)

# Acknowledgements

The authors would like to thank the staff of the Texas Advanced Computing Center for providing computational resources, and members of the Finkelstein and Wilke labs for helpful discussions. This work was supported by an NIGMS grant R01GM124141 (to I.J.F.), the Welch Foundation grant F-1808 (to I.J.F.), NIGMS grant R01 GM088344 (to C.O.W.), and the College of Natural Sciences Catalyst Award for seed funding.

# References
Contributing
============

Thank you for your interest in contributing to Opfi. Contribututions can take many forms, including:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features

Issues
------

Issues should be used to report problems or bugs with the library, to request new features, or to discuss potential changes before PRs are created.

Good bug reports have:

- A summary of the issue
- Specific steps to reproduce the problem, preferably in the form of code examples
- A description of expected behavior
- A description of actual behavior

Pull Requests
-------------

In general, we use `Github Flow <https://guides.github.com/introduction/flow/index.html>`_. This means that all changes happen through Pull Requests:

1. Fork the repository to your own Github account
2. Clone the project to your machine
3. Create a branch locally with a succinct but descriptive name
4. Commit changes to the branch
5. Push changes to your fork
6. Open a PR in our repository and include a description of the changes and a reference to the issue the PR addresses

Coding Style
------------

Contributors should attempt to adhere to the `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_ style guide when possible, although this is not currently strictly enforced. At a minimum, new contributions should follow a style that is consistent with the preexisting code base. 

Testing
-------

We use the `pytest <https://docs.pytest.org/en/6.2.x/>`_ framework for creating unit tests. Please write tests for any new code contributed to this project.

License
-------

By contributing, you agree that your contributions will be licensed under its MIT License.

Code of Conduct
---------------

We take our open source community seriously and hold ourselves and other contributors to high standards of communication. By participating and contributing to this project, you agree to uphold our `Code of Conduct <https://github.com/wilkelab/Opfi/blob/contributing-guide/CODE-OF-CONDUCT.md>`_.

Attribution
-----------

This document was adapted from the `General Contributing Guidelines <https://github.com/extendr/rextendr/blob/main/CONTRIBUTING.md>`_ of the rextendr project, and from the `contributing template <https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62>`_ developed by `briandk <https://gist.github.com/briandk>`_.Inputs and Outputs
==================

.. _building-sequence-databases:

Building sequence databases
---------------------------

To search for gene clusters with Opfi, users must compile representative protein (or nucleic acid) sequences for any genes expected in target clusters (or for any non-essential accessory genes of interest). These may be from a pre-existing, private collection of sequences (perhaps from a previous bioinformatics analysis). Alternatively, users may download sequences from a publically available database such as `Uniprot <https://www.uniprot.org/>`_ (maintained by the `European Bioinformatics Institute <https://www.ebi.ac.uk/>`_ ) or one of the `databases <https://www.ncbi.nlm.nih.gov/>`_ provided by the National Center for Biotechnology Information. 

Once target sequences have been compiled, they must be converted to an application-specific database format. Opfi currently supports :program:`BLAST+`, :program:`mmseqs2`, and :program:`diamond` for homology searching:

* `Instructions for creating sequence databases for BLAST using makeblastdb <https://www.ncbi.nlm.nih.gov/books/NBK569841/>`_
* `Instructions for creating sequence databases for mmseqs2 using mmseqs createdb <https://github.com/soedinglab/mmseqs2/wiki#searching>`_
* `Diamond makedb command options <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options>`_

The FASTA file format
#####################

Both genomic input data and reference sequence data should be in `FASTA <https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp>`_ format. This is a simple flat text representation of biological sequence data, where individual sequences are delineated by the ``>`` greater than character. For example:

.. code-block:: 

    >UniRef50_Q02ML7 CRISPR-associated endonuclease Cas1 n=1700 RepID=CAS1_PSEAB
    MDDISPSELKTILHSKRANLYYLQHCRVLVNGGRVEYVTDEGRHSHYWNIPIANTTSLLL
    GTGTSITQAAMRELARAGVLVGFCGGGGTPLFSANEVDVEVSWLTPQSEYRPTEYLQRWV
    GFWFDEEKRLVAARHFQRARLERIRHSWLEDRVLRDAGFAVDATALAVAVEDSARALEQA
    PNHEHLLTEEARLSKRLFKLAAQATRYGEFVRAKRGSGGDPANRFLDHGNYLAYGLAATA
    TWVLGIPHGLAVLHGKTRRGGLVFDVADLIKDSLILPQAFLSAMRGDEEQDFRQACLDNL
    SRAQALDFMIDTLKDVAQRSTVSA
    >UniRef50_Q2RY21 CRISPR-associated endonuclease Cas1 1 n=1034 RepID=CAS1A_RHORT
    MADPAFVPLRPIAIKDRSSIVFLQRGQLDVVDGAFVLIDQEGVRVQIPVGGLACLMLEPG
    TRITHAAIVLCARVGCLVIWVGERGTRLYAAGQPGGARADRLLFQARNALDETARLNVVR
    EMYRRRFDDDPPARRSVDQLRGMEGVRVREIYRLLAKKYAVDWNARRYDHNDWDGADIPN
    RCLSAATACLYGLCEAAILAAGYAPAIGFLHRGKPQSFVYDVADLYKVETVVPTAFSIAA
    KIAAGKGDDSPPERQVRIACRDQFRKSGLLEKIIPDIEEILRAGGLEPPLDAPEAVDPVI
    PPEEPSGDDGHRG

The sequence definition (defline) comes directly after the ``>`` character, and should be on a separate line from the sequence (which can be on one or more subsequent lines). There is no specific defline format, however, Opfi requires that, for both genomic input and sequence data, each definition line contain a unique sequence identifer. This should be a single word/token immediately following the ``>`` character (i.e. spaces between the ``>`` character and the identifier are not allowed). Any additional text on the defline is parsed as a single string, and appears in the output CSV (see :ref:`opfi-output-format`).

.. tip::

    Biological sequences downloaded from most public databases will have an accession number/identifier by default.

.. _labeling-sequences:

Annotating sequence databases
#############################

To take full advantage of the rule-based filtering methods in :mod:`operon_analyzer.rules`, users are encouraged to annotate reference sequences with a name/label that is easily searched. Labels can be as broad or as specific as is necessary to provide meaningful annotation of target gene clusters.

Gene labels are parsed from sequence deflines; specifically, Opfi looks for the second word/token following the ``>`` character. For example, the following FASTA sequence has been annotated with the label "cas1":

.. code-block:: 

    >UniRef50_Q02ML7 cas1 CRISPR-associated endonuclease Cas1 n=1700 RepID=CAS1_PSEAB
    MDDISPSELKTILHSKRANLYYLQHCRVLVNGGRVEYVTDEGRHSHYWNIPIANTTSLLL
    GTGTSITQAAMRELARAGVLVGFCGGGGTPLFSANEVDVEVSWLTPQSEYRPTEYLQRWV
    GFWFDEEKRLVAARHFQRARLERIRHSWLEDRVLRDAGFAVDATALAVAVEDSARALEQA
    PNHEHLLTEEARLSKRLFKLAAQATRYGEFVRAKRGSGGDPANRFLDHGNYLAYGLAATA
    TWVLGIPHGLAVLHGKTRRGGLVFDVADLIKDSLILPQAFLSAMRGDEEQDFRQACLDNL
    SRAQALDFMIDTLKDVAQRSTVSA

After running :class:`gene_finder.pipeline.Pipeline`, users could select candidates with hits against this sequence using the following rule set:

.. code-block:: Python

    from operon_analyzer.rules import RuleSet

    rs = RuleSet.require("cas1")

In practice, a genomics search might use a reference database of hundreds (or even thousands) of representative protein sequences, in which case labeling each sequence individually would be tedious. It is recommended to organize sequences into groups of related proteins that can be given a single label. This script uses the Python package :program:`Biopython` to annotate sequences in a multi-sequence FASTA file:

.. code-block:: Python 

    from Bio import SeqIO
    import os, sys

    def annotate_reference(prot_ref_file, label):
        records = list(SeqIO.parse(ref_fasta, "fasta"))
            
        for record in records:
            des = record.description.split()
            prot_id = des.pop(0)
            des_with_label = "{} {} {}".format(prot_id, label, " ".join(des))
            record.description = des_with_label

        SeqIO.write(records, ref_fasta, "fasta")

    if __name__ == "__main__":
        ref_fasta = sys.argv[1]
        label = sys.argv[2]
        annotate_reference(ref_fasta, label)

It is possible to use the entire sequence description (i.e. all text following the sequence identifier) as the gene label. This is particularly useful when using a pre-built database like `nr <https://www.ncbi.nlm.nih.gov/refseq/about/nonredundantproteins/>`_, which contains representative protein sequences for many different protein families. When using sequence databases that haven't been annotated, users should set ``parse_descriptions=False`` for each :class:`gene_finder.pipeline.Pipeline` ``add_step()`` method call.

Converting sequence files to a sequence database
################################################

Once reference sequences have been compiled (and, optionally, labeled) they must be converted to a sequence database format that is specific to the homology search program used. Currently, Opfi supports :program:`BLAST`, :program:`mmseqs2`, and :program:`diamond`. Each software package is automatically installed with a companion utility program for generating sequence databases. The following example shows what a typical call to :program:`makeblastdb`, the BLAST+ database utility program, might look like:

.. code-block:: bash 

    makeblastdb -in "my_sequences.fasta" -out my_sequences/db -dbtype prot -title "my_sequences" -hash_index

The command takes a text/FASTA file ``my_sequences.fasta`` as input, and writes the resulting database files to the directory ``my_sequences``. Database files are prefixed with "db". ``-dbtype prot`` specifies that the input is amino acid sequences. We use ``-title`` to name the database (required by BLAST). ``-hash_index`` directs :program:`makeblastdb` to generate a hash index of protein sequences, which can speed up computation time.

.. tip::

    :program:`mmseqs2` and :program:`diamond` have similar database creation commands, see :ref:`building-sequence-databases`. 

BLAST advanced options
----------------------

BLAST+ programs have a number of tunable parameters that can, for example, be used to adjust the sensitivity of the search algorithm. We anticipate that application defaults will be sufficient for most users; nevertheless, it is possible to use non-default program options by passing them as keyword arguments to :class:`gene_finder.pipeline.Pipeline` ``add_step()`` methods. 

For example, when using :program:`blastp` on the command line, we could adjust the number of CPUs to four by passing the argument ``-num_threads 4`` to the program. When using Opfi, this would look like ``num_threads=4``. 

Flags (boolean arguments that generally do not precede additional data) are also possible. For example, the command line flag ``-use_sw_tback`` tells :program:`blastp` to compute locally optimal Smith-Waterman alignments. The correct way to specify this behavior via the :class:`gene_finder.pipeline.Pipeline` API would be to use the argument ``use_sw_tback=True``. 

Below is a list of options accepted by Opfi. Note that some BLAST+ options are not allowed, mainly those that modify BLAST output.

.. csv-table::
    :header: "Program", "Allowed Options"

    ":program:`blastp` and :program:`psiblast`", "dbsize word_size gapopen gapextend qcov_hsp_perc xdrop_ungap xdrop_gap xdrop_gap_final searchsp sum_stats seg soft_masking matrix threshold culling_limit window_size num_threads comp_based_stats gilist seqidlist negative_gilistdb_soft_mask db_hard_mask entrez_query max_hspsbest_hit_overhang best_hit_score_edge max_target_seqsimport_search_strategy export_search_strategy num_alignments"
    ":program:`blastp` only", "task"
    ":program:`psiblast` only", "gap_trigger num_iterations out_pssm out_ascii_pssm pseudocount inclusion_ethresh"
    ":program:`blastp` (flags)", "lcase_masking ungapped use_sw_tback remote"
    ":program:`psiblast` (flags)", "lcase_masking use_sw_tback save_pssm_after_last_round save_each_pssm remote"
    ":program:`blastn`", "filtering_algorithm sum_stats window_masker_db window_size template_type version parse_deflines min_raw_gapped_score string format max_hsps taxids negative_taxids num_alignments strand off_diagonal_range subject_besthit num_sequences no_greedy negative_taxidlist culling_limit xdrop_ungap open_penalty DUST_options sorthits xdrop_gap_final negative_gilist subject use_index bool_value filename seqidlist task_name sort_hits database_name lcase_masking query_loc subject_loc sort_hsps line_length boolean db_hard_mask negative_seqidlist template_length filtering_db filtering_database penalty searchsp ungapped type gapextend db_soft_mask dbsize qcov_hsp_perc sorthsps window_masker_taxid index_name export_search_strategy float_value soft_masking gilist entrez_query show_gis best_hit_score_edge gapopen subject_input_file range html word_size best_hit_overhang perc_identity input_file num_descriptions xdrop_gap dust taxidlist max_target_seqs num_threads task remote int_value extend_penalty reward import_search_strategy num_letters"

You can read more about BLAST+ options in the `BLAST+ appendices <https://www.ncbi.nlm.nih.gov/books/NBK279684/>`_. 

.. note::

    Using advanced options with :program:`mmseqs2` and :program:`diamond` is not supported at this time. 

.. _opfi-output-format:

Opfi output format
------------------

Results from :class:`gene_finder.pipeline.Pipeline` searches are written to a single CSV file. Below is an example from the tutorial (see :ref:`example-usage`):

.. csv-table::
    :file: csv/example_output.csv
    :header-rows: 0

The first two columns contain the input genome/contig sequence ID (sometimes called an accession number) and the coordinates of the candidate gene cluster, respectively. Since an input file can have multiple genomic sequences, these two fields together uniquely specify a candidate gene cluster. Each row represents a single annotated feature in the candidate locus. Features from the same candidate are always grouped together in the CSV. 

Descriptions of each output field are provided below. Alignment statistic naming conventions are from the BLAST documentation, see `BLAST+ appendices <https://www.ncbi.nlm.nih.gov/books/NBK279684/>`_ (specifically "outfmt" in table C1). This `glossary <https://www.ncbi.nlm.nih.gov/books/NBK62051/>`_ of common BLAST terms may also be useful in interpreting alignment statistic meaning. 

.. csv-table::
    :file: csv/fieldnames.csv
    :header-rows: 1
.. _example-usage:

Example Usage
=============

Example 1: Finding CRISPR-Cas systems in a cyanobacteria genome
---------------------------------------------------------------

In this example, we will annotate and visualize CRISPR-Cas systems in the cyanobacteria species Rippkaea orientalis. CRISPR-Cas is a widespread bacterial defense system, found in at least 50% of all known prokaryotic species. This system is significant in that it can be leveraged as a precision gene editing tool, an advancement that was awarded the 2020 Nobel Prize in Chemistry. The genome of R. orientalis harbors two complete CRISPR-Cas loci (one chromosomal, and one extrachromosomal/plasmid).

You can download the complete assembled genome `here <https://www.ncbi.nlm.nih.gov/assembly/GCF_000024045.1/>`_; it is also available at `<https://github.com/wilkelab/Opfi>`_ under ``tutorials``, along with the other data files necessary to run these examples, and an interactive jupyter notebook version of this tutorial. 

This tutorial assumes the user has already installed Opfi and all dependencies (if installing with conda, this is done automatically). Some familiarity with BLAST and the basic homology search algorithm may also be helpful, but is not required. 

1. Use the makeblastdb utility to convert a Cas protein database to BLAST format
################################################################################

We start by converting a Cas sequence database to a format that BLAST can recognize, using the command line utility :program:`makeblastdb`, which is part of the core NCBI BLAST+ distribution. A set of ~20,000 non-redundant Cas sequences, downloaded from `Uniprot <https://www.uniprot.org/uniref/>`_ is available as a tar archive ``tutorials/cas_database.tar.gz`` . We'll make a new directory, "blastdb", and extract sequences there:

.. code-block:: bash

    mkdir blastdb
    cd blastdb && tar -xzf cas_database.tar.gz && cd ..

Next, create two BLAST databases for the sequence data: one containing Cas1 sequences only, and another that contains the remaining Cas sequences.

.. code-block:: bash

    cd blastdb && cat cas1.fasta | makeblastdb -dbtype prot -title cas1 -hash_index -out cas1_db && cd ..
    cd blastdb && cat cas[2-9].fasta cas1[0-2].fasta casphi.fasta | makeblastdb -dbtype prot -title cas_all -hash_index -out cas_all_but_1_db && cd ..

``-dbtype prot`` simply tells :program:`makeblastdb` to expect amino acid sequences. We use ``-title`` and ``-out`` to name the database (required by BLAST) and to prefix the database files, respectively. ``-hash_index`` directs :program:`makeblastdb` to generate a hash index of protein sequences, which can speed up computation time.

2. Use Gene Finder to search for CRISPR-Cas loci
################################################

CRISPR-Cas systems are extremely diverse. The most recent `classification effort <https://www.nature.com/articles/s41579-019-0299-x>`_ identifies 6 major types, and over 40 subtypes, of compositionally destinct systems. Although there is sufficent sequence similarity between subtypes to infer the existence of a common ancestor, the only protein family present in the majority of CRISPR-cas subtypes is the conserved endonuclease Cas1. For our search, we will define candidate CRISPR-cas loci as having, minimally, a cas1 gene.

First, create another directory for output:

.. code-block:: bash

    mkdir example_1_output

The following bit of code uses Opfi's :mod:`gene_finder.pipeline` module to search for CRISPR-Cas systems:

.. code-block:: Python

    from gene_finder.pipeline import Pipeline
    import os

    genomic_data = "GCF_000024045.1_ASM2404v1_genomic.fna.gz"
    output_directory = "example_1_output"

    p = Pipeline()
    p.add_seed_step(db="blastdb/cas1_db", name="cas1", e_val=0.001, blast_type="PROT", num_threads=1)
    p.add_filter_step(db="blastdb/cas_all_but_1_db", name="cas_all", e_val=0.001, blast_type="PROT", num_threads=1)
    p.add_crispr_step()

    # use the input filename as the job id
    # results will be written to the file <job id>_results.csv
    job_id = os.path.basename(genomic_data)
    results = p.run(job_id=job_id, data=genomic_data, output_directory=output_directory, min_prot_len=90, span=10000, gzip=True)

First, we initialize a :class:`gene_finder.pipeline.Pipeline` object, which keeps track of all search parameters, as well as a running list of systems that meet search criteria. Next, we add three search steps to the pipeline:

1. :meth:`gene_finder.pipeline.Pipeline.add_seed_step` : BLAST is used to search the input genome against a database of Cas1 sequences. Regions around putative Cas1 hits become the intial candidates, and the rest of the genome is ignored.
2. :meth:`gene_finder.pipeline.Pipeline.add_filter_step` : Candidate regions are searched for any additional Cas genes. Candidates without at least one additional putative Cas gene are also discarded.
3. :meth:`gene_finder.pipeline.Pipeline.add_crispr_step` : Remaining candidates are annotated for CRISPR repeat sequences using PILER-CR. 

Finally, we run the pipeline, executing steps in the order they we added. ``min_prot_len`` sets the minimum length (in amino acid residues) of hits to keep (really short hits are unlikely real protein encoding genes). ``span`` is the region directly up- and downstream of initial hits. So, each candidate system will be about 20 kbp in length. Results are written to a single CSV file. Final candidate loci contain at least one putative Cas1 gene and one additional Cas gene. As we will see, this relatively permissive criteria captures some non-CRISPR-Cas loci. Opfi has additional modules for reducing unlikely systems after the gene finding stage.

3. Visualize annotated CRISPR-Cas gene clusters with Operon Analyzer
####################################################################

It is sometimes useful to visualize candidate systems, especially during the exploratory phase of a genomics survey. Opfi provides a few functions for visualizing candidate systems in :mod:`operon_analyzer.visualize`. We'll use these to visualize the CRISPR-Cas gene clusters in R. orientalis:

.. code-block:: Python

    import csv
    import sys
    from operon_analyzer import load, visualize

    feature_colors = { "cas1": "lightblue",
                        "cas2": "seagreen",
                        "cas3": "gold",
                        "cas4": "springgreen",
                        "cas5": "darkred",
                        "cas6": "thistle",
                        "cas7": "coral",
                        "cas8": "red",
                        "cas9": "palegreen",
                        "cas10": "yellow",
                        "cas11": "tan",
                        "cas12": "orange",
                        "cas13": "saddlebrown",
                        "casphi": "olive",
                        "CRISPR array": "purple"
                        }

    # read in the output from Gene Finder and create a gene diagram for each cluster (operon)
    with open("example_1_output/GCF_000024045.1_ASM2404v1_genomic.fna.gz_results.csv", "r") as operon_data:
        operons = load.load_operons(operon_data)
        visualize.plot_operons(operons=operons, output_directory="example_1_output", feature_colors=feature_colors, nucl_per_line=25000)

Running this script produces the following three gene diagrams, one for each system in the input CSV:

.. _fig-1:
.. figure:: img/operon_image_1.png
    
    A CRISPR-Cas system in the chromosome of R. orientalis.  

.. _fig-2:
.. figure:: img/operon_image_2.png

    A second CRISPR-Cas system in R. orientalis plasmid 1.  

.. _fig-3:
.. figure:: img/operon_image_3.png

    An R. orientalis locus with a putative CRISPR-Cas gene.

   
We can see that both CRISPR-Cas systems were identified (:numref:`fig-1` and :numref:`fig-2`). We also see some systems that don't resemble functional CRISPR-Cas operons (:numref:`fig-3`). Because we used a relatively permissive e-value threshhold of 0.001 when running BLAST, Opfi retained regions with very low sequence similarity to true CRISPR-Cas genes. In fact, these regions are likely not CRISPR-Cas loci at all. Using a lower e-value would likely eliminate these "false positive" systems, but :mod:`operon_analyzer.rules` exposes functions for filtering out unlikely candidates *after* the intial BLAST search. 

In general, we have found that using permissive BLAST parameters intially, and then filtering or eliminating candidates during the downstream analysis, is an effective way to search for gene clusters in large amounts of genomic/metagenomic data. In this toy example, we could re-run BLAST many times without significant cost. But on a more realistic dataset, needing to re-do the computationally expensive homology search could detrail a project. Since the optimal search parameters may not be known *a priori*, it can be better to do a permissive homology search initially, and then narrow down results later.

Finally, clean up the temporary directories, if desired:

.. code-block:: bash

    rm -r example_1_output blastdb

Example 2: Filter and classify CRISPR-Cas systems based on genomic composition
------------------------------------------------------------------------------

As discussed in the previous example, known CRISPR-Cas systems fall into 6 broad categories, based on the presence of particular "signature" genes, as well as overall composition and genomic architecture. In this example, we will use Opfi to search for and classify CRISPR-Cas systems in ~300 strains of fusobacteria. 

This dataset was chosen because it is more representative (in magnitude) of what would be encountered in a real genomics study. Additionally, the fusobacteria phylum contains a variety of CRISPR-Cas subtypes. Given that the homology search portion of the analysis takes several hours (using a single core) to complete, we have pre-run Gene Finder using the same setup as the previous example. 

1. Make another temporary directory for output:
###############################################

.. code-block:: bash

    mkdir example_2_output

2. Filter Gene Finder output and extract high-confidence CRISPR-Cas systems
###########################################################################

The following code reads in unfiltered output from :class:`gene_finder.pipeline.Pipeline` and applies a set of conditions ("rules") to accomplish two things:
1. Select (and bin) systems according to type, and,
2. Eliminate candidates that likely do not represent true CRISPR-Cas systems

To do this, we'll leverage the :mod:`operon_analyzer.rules` and :mod:`operon_analyzer.analyze` modules.

.. code-block:: Python

    from operon_analyzer import analyze, rules


    fs = rules.FilterSet().pick_overlapping_features_by_bit_score(0.9)
    cas_types = ["I", "II", "III", "V"]

    rulesets = []
    # type I rules
    rulesets.append(rules.RuleSet().contains_group(feature_names = ["cas5", "cas7"], max_gap_distance_bp = 1000, require_same_orientation = True) \
                                .require("cas3"))
    # type II rules
    rulesets.append(rules.RuleSet().contains_at_least_n_features(feature_names = ["cas1", "cas2", "cas9"], feature_count = 3) \
                                .minimum_size("cas9", 3000))
    # type III rules
    rulesets.append(rules.RuleSet().contains_group(feature_names = ["cas5", "cas7"], max_gap_distance_bp = 1000, require_same_orientation = True) \
                                .require("cas10"))
    # type V rules
    rulesets.append(rules.RuleSet().contains_at_least_n_features(feature_names = ["cas1", "cas2", "cas12"], feature_count = 3))

    for rs, cas_type in zip(rulesets, cas_types):
        with open("refseq_fusobacteria.csv", "r") as input_csv:
            with open(f"example_2_output/refseq_fuso_filtered_type{cas_type}.csv", "w") as output_csv:
                analyze.evaluate_rules_and_reserialize(input_csv, rs, fs, output_csv)

The rule sets are informed by an established CRISPR-Cas classification system, which you can learn more about in this `paper <https://www.nature.com/articles/s41579-019-0299-x>`_ . The most recent system recognizes 6 major CRISPR-Cas types, but since fusobacteria doesn't contain type IV or VI systems that can be identified with our protein dataset, we didn't define the corresponding rule sets.

3. Verify results with additional visualizations
################################################

Altogther, this analysis will identify several hundred systems. We won't look at each system individually (but you are free to do so!). For the sake of confirming that the code ran as expected, we'll create gene diagrams for just the type V systems, since there are only two:

.. code-block:: Python

    import csv
    import sys
    from operon_analyzer import load, visualize

    feature_colors = { "cas1": "lightblue",
                        "cas2": "seagreen",
                        "cas3": "gold",
                        "cas4": "springgreen",
                        "cas5": "darkred",
                        "cas6": "thistle",
                        "cas7": "coral",
                        "cas8": "red",
                        "cas9": "palegreen",
                        "cas10": "yellow",
                        "cas11": "tan",
                        "cas12": "orange",
                        "cas13": "saddlebrown",
                        "casphi": "olive",
                        "CRISPR array": "purple"
                        }

    # read in the output from Gene Finder and create a gene diagram for each cluster (operon)
    with open("example_2_output/refseq_fuso_filtered_typeV.csv", "r") as operon_data:
        operons = load.load_operons(operon_data)
        visualize.plot_operons(operons=operons, output_directory="example_2_output", feature_colors=feature_colors, nucl_per_line=25000)

The plotted systems should look like this:

.. figure:: img/operon_image_4.png

    A type V CRISPR-Cas system.

.. figure:: img/operon_image_5.png

    A second type V CRISPR-Cas system.


Finally, clean up the temporary output directory, if desired:

.. code-block:: bash

    rm -r example_2_output
Gene Finder
============

gene\_finder.pipeline
---------------------

.. automodule:: gene_finder.pipeline
   :members:
   :member-order: bysource
Getting Started
===============

.. _installation:

Installation
------------

The recommended way to install Opfi is with `Bioconda <https://bioconda.github.io/>`_, which requires the `conda <https://docs.conda.io/en/latest/>`_ package manager. This will install Opfi and all of its dependencies (which you can read more about below, see :ref:`dependencies`).

Currently, Bioconda supports only 64-bit Linux and Mac OS. Windows users can still install Opfi with pip (see below); however, the complete installation procedure has not been fully tested on a Windows system. 

.. _install-with-conda:

Install with conda (Linux and Mac OS only)
##########################################

First, set up conda and Bioconda following the `quickstart <https://bioconda.github.io/user/install.html>`_ guide. Once this is done, run:

.. code-block:: bash

    conda install -c bioconda opfi

And that's it! Note that this will install Opfi in the conda environment that is currently active. To create a fresh environment with Opfi installed, do:

.. code-block:: bash

    conda create --name opfi-env -c bioconda opfi
    conda activate opfi-env

.. _install-with-pip:

Install with pip
################

This method does not automatically install non-Python dependencies, so they will need to be installed separately, following their individual installation instructions. A complete list of required software is provided below, see :ref:`dependencies`. Once this step is complete, install Opfi with pip by running:

.. code-block:: bash

    pip install opfi

Install from source
###################

Finally, the latest development build may be installed directly from Github. First, non-Python :ref:`dependencies` will need to be installed in the working environment. An easy way to do this is to first install Opfi with conda using the :ref:`install-with-conda` method (we'll re-install the development version of the Opfi package in the next step). Alternatively, dependencies can be installed individually.

Once dependencies have been installed in the working environment, run the following code to download and install the development build:

.. code-block:: bash

    git clone https://github.com/wilkelab/Opfi.git
    cd Opfi
    pip install . # or pip install -e . for an editable version
    pip install -r requirements # if conda was used, this can be skipped

Testing the build
#################

Regardless of installation method, users can download and run Opfi's suite of unit tests to confirm that the build is working as expected. First download the tests from Github:

.. code-block:: bash

    git clone https://github.com/wilkelab/Opfi
    cd Opfi

And then run the test suite using pytest:

.. code-block:: bash

    pytest --runslow --runmmseqs --rundiamond

This may take a minute or so to complete. 

.. _dependencies:    

Dependencies
------------

Opfi uses the following bioinformatics software packages to find and annotate genomic features:

.. csv-table:: Software dependencies
   :header: "Application", "Description"

   "`NCBI BLAST+ <https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs>`_", "Protein and nucleic acid homology search tool"
   "`Diamond <https://github.com/bbuchfink/diamond>`_", "Alternative to BLAST+ for fast protein homology searches"
   "`MMseqs2 <https://github.com/soedinglab/MMseqs2>`_", "Alternative to BLAST+ for fast protein homology searches"
   "`PILER-CR <https://www.drive5.com/pilercr/>`_", "CRISPR repeat detection"
   "`Generic Repeat Finder <https://github.com/bioinfolabmu/GenericRepeatFinder>`_", "Transposon-associated repeat detection"

The first three (BLAST+, Diamond, and MMseqs2) are popular homology search applications, that is, programs that look for local similarities between input sequences (either protein or nucleic acid) and a target. These are used by Opfi in :class:`gene_finder.pipeline.Pipeline` for annotation of genes or non-coding regions of interest in the input genome/contig. The user specifies which homology search tool to use during pipeline setup (see :class:`gene_finder.pipeline.Pipeline` for details). Note that the BLAST+ distribution contains multiple programs for homology searching, three of which (blastp, blastn, and PSI-BLAST) are currently supported by Opfi. 

The following table summarizes the main difference between each homology search program. It may help users decide which application will best meet their needs. Note that performance tests are inherently hardware and context dependent, so this should be taken as a loose guide, rather than a definitive comparison. 

.. csv-table:: Comparison of homology search programs supported by Opfi
    :header: "Application", "Relative sensitivity", "Relative speed", "Requires a protein or nucleic acid sequence database?"

    "Diamond", `+`, `++++`, "protein"
    "MMseqs2", `++`, `+++`, "protein"
    "blastp", `+++`, `++`, "protein"
    "PSI-BLAST", `++++`, `+`, "protein"
    "blastn", "NA", "NA", "nucleic acid"

The last two software dependencies, PILER-CR and Generic Repeat Finder (GRF), deal with annotation of repetive sequences in DNA. PILER-CR identifies CRISPR arrays, regions of alternatating ~30 bp direct repeat and variable sequences that play a role in prokaryotic immunity. GRF identifies repeats associated with transposable elements, such as terminal inverted repeats (TIRs) and long terminal repeats (LTRs).
Operon Analyzer
===============

The following modules comprise the core Operon Analyzer functionality.

operon\_analyzer.genes
----------------------

.. automodule:: operon_analyzer.genes
   :members:
   :member-order: bysource

operon\_analyzer.rules
----------------------

.. automodule:: operon_analyzer.rules
   :members:
   :member-order: bysource

operon\_analyzer.analyze
------------------------

.. automodule:: operon_analyzer.analyze
   :members:
   :member-order: bysource

operon\_analyzer.visualize
--------------------------

.. automodule:: operon_analyzer.visualize
   :members:
   :member-order: bysource
   :exclude-members: build_image_filename, calculate_adjusted_operon_bounds, plot_operon_pair, make_operon_pairs, create_operon_figure, save_operon_figure, save_pair_figure, build_operon_dictionary

operon\_analyzer.overview
-------------------------

.. automodule:: operon_analyzer.overview
   :members:
   :member-order: bysource

operon\_analyzer.reannotation
-----------------------------

.. automodule:: operon_analyzer.reannotation
   :members:
   :member-order: bysource

operon\_analyzer.load
---------------------

.. automodule:: operon_analyzer.load
   :members:
   :member-order: bysource

operon\_analyzer.parse
----------------------

.. automodule:: operon_analyzer.parse
   :members:
   :exclude-members: parse_coordinates, read_pipeline_output

The next set of modules expose simple functions for dealing with CRISPR array sequences. Presumably, these would only be useful to researchers interested in CRISPR-Cas genomic systems. 

operon\_analyzer.piler\_parse
-----------------------------

.. automodule:: operon_analyzer.piler_parse
   :members:
   :member-order: bysource

operon\_analyzer.repeat\_finder
-------------------------------

.. automodule:: operon_analyzer.repeat_finder
   :members:
   :member-order: bysource
   :exclude-members: BufferedSequence

operon\_analyzer.spacers
------------------------

.. automodule:: operon_analyzer.spacers
   :members:
   :member-order: bysource
API Reference
=============

.. toctree::
   :maxdepth: 4

   gene_finder
   operon_analyzer.. Opfi documentation master file, created by
    sphinx-quickstart on Mon Aug 10 12:29:45 2020.
    You can adapt this file completely to your liking, but it should at least
    contain the root `toctree` directive.

Opfi
====

Welcome to the Opfi documentation site! Opfi is a modular, rule-based framework for creating gene cluster identification pipelines, particularly for large genomics or metagenomics datasets. 

Opfi is implemented entirely in Python, and can be downloaded with conda or the from the Python Package Index. It consists of two major modules: Gene Finder, for discovery of novel gene clusters, and Operon Analyzer, for rule-based filtering, deduplication, visualization, and re-annotation of systems identified by Gene Finder.

Contents
--------

.. toctree::
    :maxdepth: 2
    
    installation
    examples
    tips
    modules
    contributing
