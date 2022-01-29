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
