<p align="center">
  <img src="https://github.com/BioMeCIS-Lab/OpenOmics/raw/master/openomics_web/assets/openomics_logo.png" max-height="200">
</p>

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![pyOpenSci](https://tinyurl.com/y22nb8up)](https://github.com/pyOpenSci/software-review/issues/31)
[![status](https://joss.theoj.org/papers/aca43e3c2989a803b514faef72dd3294/status.svg)](https://joss.theoj.org/papers/aca43e3c2989a803b514faef72dd3294)
[![DOI](https://zenodo.org/badge/125549505.svg)](https://zenodo.org/badge/latestdoi/125549505)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics/branch/master/graph/badge.svg?token=6N1UZ27MPH)](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics)

This Python package provide a series of tools to integrate and explore the genomics, transcriptomics, proteomics, and
clinical data (aka multi-omics data). With interfaces to popular annotation databases and scalable data-frame manipulation tools, OpenOmics facilitates the common
data wrangling tasks when preparing data for RNA-seq bioinformatics analysis.

Documentation ([Latest](https://openomics.readthedocs.io/en/latest/)
| [Stable](https://openomics.readthedocs.io/en/stable/))
| [OpenOmics at a glance](https://openomics.readthedocs.io/en/latest/usage/getting-started.html)

## Features
OpenOmics assist in integration of heterogeneous multi-omics bioinformatics data. The library provides a Python API as well as an interactive Dash web interface.
It features support for:
- Genomics, Transcriptomics, Proteomics, and Clinical data.
- Harmonization with 20+ popular annotation, interaction, disease-association databases.

OpenOmics also has an efficient data pipeline that bridges the popular data manipulation Pandas library and Dask distributed processing to address the following use cases:

- Providing a standard pipeline for dataset indexing, table joining and querying, which are transparent and customizable
  for end-users.
- Providing Efficient disk storage for large multi-omics dataset with Parquet data structures.
- Integrating various data types including interactions and sequence data, then exporting to NetworkX graphs or data generators for down-stream machine learning.
- Accessible by both developers and scientists with a Python API that works seamlessly with an external Galaxy tool interface or the built-in Dash web interface (WIP).


## Installation via pip:

```
$ pip install openomics
```

## 

## Citations
The journal paper for this scientific package was reviewed by JOSS at <https://joss.theoj.org/papers/10.21105/joss.03249#>, and can be cited with:

    # BibTeX
    @article{Tran2021,
      doi = {10.21105/joss.03249},
      url = {https://doi.org/10.21105/joss.03249},
      year = {2021},
      publisher = {The Open Journal},
      volume = {6},
      number = {61},
      pages = {3249},
      author = {Nhat C. Tran and Jean X. Gao},
      title = {OpenOmics: A bioinformatics API to integrate multi-omics datasets and interface with public databases.},
      journal = {Journal of Open Source Software}
    }



## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [pyOpenSci/cookiecutter-pyopensci](https://github.com/pyOpenSci/cookiecutter-pyopensci) project template, based off [audreyr/cookiecutter-pypackage](https://github.com/audreyr/cookiecutter-pypackage).
* openOmics version:
* Python version:
* Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
---
title: 'OpenOmics: A bioinformatics API to integrate multi-omics datasets and interface with public databases.'
tags:
  - Python
  - bioinformatics
  - multiomics
  - data integration
  - big data
authors:
  - name: Nhat C. Tran^[corresponding author]
    orcid: 0000-0002-2575-9633
    affiliation: 1
  - name: Jean X. Gao
    affiliation: 1
affiliations:
  - name: Department of Computer Science and Engineering, The University of Texas at Arlington
    index: 1
date: 25 January 2021
bibliography: paper.bib
---

# Summary

Leveraging large-scale multi-omics data is emerging as the primary approach for systemic research of human diseases and
general biological processes. As data integration and feature engineering are the vital steps in these bioinformatics
projects, there currently lacks a tool for standardized preprocessing of heterogeneous multi-omics and annotation data
within the context of a clinical cohort. OpenOmics is a Python library for integrating heterogeneous multi-omics data
and interfacing with popular public annotation databases, e.g., GENCODE, Ensembl, BioGRID. The library is designed to be
highly flexible to allow the user to parameterize the construction of integrated datasets, interactive to assist complex
data exploratory analyses, and scalable to facilitate working with large datasets on standard machines. In this paper,
we demonstrate the software design choices to support the wide-ranging use cases of OpenOmics with the goal of
maximizing usability and reproducibility of the data integration framework.

# Statement of need

Recent advances in sequencing technology and computational methods have enabled the means to generate large-scale,
high-throughput multi-omics data [@lappalainen2013transcriptome], providing unprecedented research opportunities for
cancer and other diseases. These methods have already been applied to a number of problems within bioinformatics, and
indeed several integrative disease
studies [@zhang2014proteogenomic; @cancer2014comprehensive; @ren2016integration; @hassan2020integration]. In addition to
the genome-wide measurements of different genetic characterizations, the growing public knowledge-base of functional
annotations [@rnacentral2016rnacentral; @derrien2012gencode], experimentally-verified
interactions [@chou2015mirtarbase; @yuan2013npinter; @chou2017mirtarbase; @oughtred2019biogrid], and gene-disease
associations [@huang2018hmdd; @pinero2016disgenet; @chen2012lncrnadisease] also provides the prior-knowledge essential
for system-level analyses. Leveraging these data sources allow for a systematic investigation of disease mechanisms at
multiple molecular and regulatory layers; however, such task remains nontrivial due to the complexity of multi-omics
data.

While researchers have developed several mature tools to access or analyze a particular single omic data
type [@wolf2018scanpy; @stuart2019integrative], the current state of integrative data platforms for multi-omics data is
lacking due to three reasons. First, pipelines for data integration carry out a sequential tasks that does not process
multi-omics datasets holistically. Second, the vast size and heterogeneity of the data poses a challenge on the
necessary data storage and computational processing. And third, implementations of data pipelines are close-ended for
down-stream analysis or not conductive to data exploration use-cases. Additionally, there is currently a need for
increased transparency in the process of multi-omics data integration, and a standardized data preprocessing strategy is
important for the interpretation and exchange of bioinformatic projects. Currently, there exist very few systems that,
on the one hand, supports standardized handling of multi-omics datasets but also allows to query the integrated dataset
within the context of a clinical cohort.

# Related works

There are several existing platforms that aids in the integration of multi-omics data, such as Galaxy, Anduril, MixOmics
and O-Miner. First, Galaxy [@boekel2015multi] and Anduril [@cervera2019anduril] are mature platforms and has an
established workflow framework for genomic and transcriptomic data analysis. Galaxy contains hundreds of
state-of-the-art tools of these core domains for processing and assembling high-throughput sequencing data. Second,
MixOmics [@rohart2017mixomics] is an R library dedicated to the multivariate analysis of biological data sets with a
specific focus on data exploration, dimension reduction and visualisation. Third, O-Miner [@sangaralingam2019multi] is
web tool that provides a pipeline for analysis of both transcriptomic and genomic data starting from raw image files
through in-depth bioinformatics analysis. However, as large-scale multi-omic data analysis demands continue to grow, the
technologies and data analysis needs continually change to adapt with `big data`. For instance, the data manipulation
required for multi-omics integration requires a multitude of complex operations, but the point and click interface given
in existing Galaxy tools can be limiting or not computationally efficient. Although the MixOmics toolkit provides an R
programming interface, it doesn't yet leverage high-performance distributed storage or computing resources. Finally,
while O-Miner can perform end-to-end analysis in an integrated platform, its interim analysis results cannot be exported
elsewhere for down-stream analysis.

![Overall OpenOmics System Architecture, Data Flow, and Use Cases.\label{architecture}](figure.pdf)

# The OpenOmics library

OpenOmics consists of two core modules: multi-omics integration and annotation interface. An overview visualization of
the OpenOmics system architecture is provided in \autoref{architecture}.

## Multi-omics integration

Tabular data are everywhere in bioinformatics. To record expression quantifications, annotations, or variant calls, data
are typically stored in various tabular-like formats, such as BED, GTF, MAF, and VCF, which can be preprocessed and
normalized to row indexed formats. Given any processed single-omic dataset, the library generalizes the data as a
tabular structure where rows correspond to observation samples and columns correspond to measurements of different
biomolecules. The core functionality of the Multi-omics Integration module is to integrate the multiple single-omic
datasets for the overlapping samples. By generating multi-omics data for the same set of samples, our tool can provide
the necessary data structure to develop insights into the flow of biological information across multiple genome,
epigenome, transcriptome, proteome, metabolome and phenome levels. The user can import and integrate the following
supported omic types:

- Genomics: single nucleotide variants (SNV), copy number variation (CNV)
- Epigenomics: DNA methylation
- Transcriptomics: RNA-Seq, miRNA expression, lncRNA expression, microarrays
- Proteomics: reverse phase protein array (RPPA), iTRAQ

After importing each single omics data, OpenOmics stores a Pandas Dataframe [@mckinney-proc-scipy-2010] that is flexible
for a wide range of tabular operations. For instance, the user is presented with several functions for preprocessing of
the expression quantifications to normalize, filter outliers, or reduce noise.

Within a study cohort, the clinical characteristics are crucial for the study of a disease or biological phenomenon. The
user can characterize the set of samples using the Clinical Data structure, which is comprised of two levels: Patient
and Biospecimen. A Patient can have attribute fields on demographics, clinical diagnosis, disease progression, treatment
responses, and survival outcomes. Typically, multi-omics data observations are captured at the Biospecimen level and
each Patient can have multiple Biospecimens. OpenOmics tracks the ID's of biospecimens and the patient it belongs to, so
the multi-omics data are organized in a hierarchical order to enable aggregated operations.

## Annotation interface

After importing and integrating the multi-omic data, the user can supplement their dataset with various annotation
attributes from public data repositories such as GENCODE, Ensembl, and RNA Central. With just a few operations, the user
can easily download a data repository of choice, select relevant attributes, and efficiently join a variable number of
annotation columns to their genomics, transcriptomics, and proteomics data. The full list of databases and the
availability of annotation attributes is listed in Table 1.

For each public database, the Annotation Interface module provides a series of interfaces to perform specific importing,
preprocessing, and annotation tasks. At the import step, the module can either fetch the database files via a
file-transfer-protocol (ftp) URL or load a locally downloaded file. At this step, the user can specify the species,
genome build, and version of the database by providing a ftp URL of choice. To streamline this process, the module
automatically caches downloaded file to disk, uncompress them, and handle different file extensions, including FASTA,
GTF, VCF, and other tabular formats. Then, at the preprocessing step, the module selects only the relevant attribute
fields specified by the user and perform necessary data cleanings. Finally, the annotation data can be annotated to an
omics dataset by performing a SQL-like join operation on a user-specified index of the biomolecule name or ID. If the
user wishes to import an annotation database not yet included in OpenOmics, they can extend the Annotation Dataset API
to specify their own importing, preprocessing, and annotation tasks in an object-oriented manner.

An innovative feature of our integration module is the ability to cross-reference the gene ID's between different
annotation systems or data sources. When importing a dataset, the user can specify the level of genomic index, such as
at the gene, transcript, protein, or peptide level, and whether it is a gene name or gene ID. Since multiple
single-omics datasets can use different gene nomenclatures, the user is able to convert between the different gene
indexing methods by reindexing the annotation dataframe with a index column of choice. This not only allows the
Annotation Interface to select and join the annotation data to the correct index level, but also allow the user to
customize the selection and aggregation of biological measurements at different levels.

| Data Repository | Annotation Data Available                       | Index     | # entries  |
| --------------- | ----------------------------------------------- | --------- | ---------- |
| GENCODE         | Genomic annotations, primary sequence           | RNAs      | 60,660     |
| Ensembl         | Genomic annotations                             | Genes     | 232,186    |
| MiRBase         | MicroRNA sequences and annotatinos              | MicroRNAs | 38,589     |
| RNA Central     | ncRNA sequence and annotation collection        | ncRNAs    | 14,784,981 |
| NONCODE         | lncRNA sequences and annotations                | LncRNAs   | 173,112    |
| lncrnadb        | lncRNA functional annotations                   | LncRNAs   | 100        |
| Pfam            | Protein family annotation                       | Proteins  | 18,259     |
| Rfam            | RNA family annotations                          | ncRNAs    | 2,600      |
| Gene Ontology   | Functional, cellular, and molecular annotations | Genes     | 44,117     |
| KEGG            | High-level functional pathways                  | Genes     | 22,409     |
| DisGeNet        | gene-disease associations                       | Genes     | 1,134,942  |
| HMDD            | microRNA-disease associations                   | MicroRNAs | 35,547     |
| lncRNAdisease   | lncRNA-disease associations                     | LncRNAs   | 3,000      |
| OMIM            | Ontology of human diseases                      | Diseases  | 25,670     |

Table 1: Public annotation databases and availability of data in the Human genome.

# System design

This section describes the various implementation details behind the scalable processing and efficient data storage, and
the design choices in the development operations.

While the in-memory Pandas dataframes utilized in our data structures are fast, they have size and speed limitations
when the dataset size approaches the system memory limit. When this is an issue, the user can enable out-of-memory
distributed data processing on all OpenOmics operations, implemented by the Dask
framework [@matthew_rocklin-proc-scipy-2015]. When memory resources is limited, data in a Dask dataframe can be read
directly from disk and is only brought into memory when needed during computations (also called lazy evaluations). When
performing data query operations on Dask dataframes, a task graph containing each operation is built and is only
evaluated on command, in a process called lazy loading.

Operations on Dask dataframes are the same as Pandas dataframes, but can utilize multiple workers and can scale up to
clusters by connecting to a cluster client with minimal configuration. To enable this feature in OpenOmics, the user
simply needs to explicitly enable an option when importing an omics dataset, importing an annotation/interaction
database, or importing a MultiOmics file structure on disk.

## Software requirements

OpenOmics is distributed as a readily installable Python package from the Python Package Index (PyPI) repository. For
users to install OpenOmics in their own Python environment, several software dependencies are automatically downloaded
to reproduce the computing environment.

OpenOmics is compatible with Python 3.6 or higher, and is operational on both Linux and Windows operating systems. The
software requires as little as 4 GB of RAM and 2 CPU cores, and can computationally scale up to large-memory
multi-worker distributed systems such as a compute cluster. To take advantage of increased computational resource,
OpenOmics simply requires one line of code to activate parallel computing functionalities.

## Development operations

We developed OpenOmics following modern software best-practices and package publishing standards. For the version
control of our source-code, we utilized a public GitHub repository which contains two branches, master and develop. The
master branch contains stable and well-tested releases of the package, while the develop branch is used for building new
features or software refactoring. Before each version is released, we utilize Github Actions for continuous integration,
building, and testing for version and dependency compatibility. Our automated test suite covers essential functions of
the package and a reasonable range of inputs and conditions.

# Conclusion

A standardized data preprocessing strategy is essential for the interpretation and exchange of bioinformatics research.
OpenOmics provides researchers with the means to consistently describe the processing and analysis of their experimental
datasets. It equips the user, a bioinformatician, with the ability to preprocess, query, and analyze data with modern
and scalable software technology. As the wide array of tools and methods available in the public domain are largely
isolated, OpenOmics aims toward a uniform framework that can effectively process and analyze multi-omics data in an
end-to-end manner along with biologist-friendly visualization and interpretation.

# Acknowledgements

N/A.

# References
# Welcome to OpenOmics's documentation!

[![PyPI version](https://badge.fury.io/py/openomics.svg)](https://badge.fury.io/py/openomics)
[![Documentation Status](https://readthedocs.org/projects/openomics/badge/?version=latest)](https://openomics.readthedocs.io/en/latest/?badge=latest)
[![OpenOmics](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml/badge.svg?branch=master)](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
[![codecov](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics/branch/master/graph/badge.svg?token=WAN3PJwM42)](https://codecov.io/gh/BioMeCIS-Lab/OpenOmics)

This Python package provide a series of tools to integrate and query the genomics, transcriptomics, proteomics, and
clinical data (aka, the multi-omics data). With scalable data-frame manipulation tools, OpenOmics facilitates the common
coding tasks when preparing data for RNA-seq bioinformatics analysis.

## Features

- Provides a bioinformatics workflow to generate integrative results from multi-omics data.
- Facilitates integration of various bio-databases, multi-omics expression, genomics, and clinical data.
- Highly flexible to different data types and missing data.
- Provides researchers with means to consistently store and explore their experimental datasets.
- Enables scalable performance with parallel computing, while easily configurable to deploy on both single machine and a
  cluster.

## Table of Content

```{toctree}
:maxdepth: 2
:caption: Using OpenOmics Python API
:name: mastertoc
installation
usage/getting-started
usage/import-your-dataset
usage/annotate-external-databases
usage/network-data
usage/preprocess-downstream-analysis
```

```{toctree}
:maxdepth: 1
:caption: API Reference

modules/openomics.multiomics
modules/openomics.annotate
modules/openomics.database.annotation
modules/openomics.database.sequence
modules/openomics.database.interaction
modules/openomics.database.disease
modules/openomics.database.ontology
modules/openomics.utils
```


```{toctree}
:maxdepth: 1
:caption: MISC

misc/faq
```


```{toctree}
:maxdepth: 1
:caption: Contributing and releases

contributing
history
```
# Installation

## Stable release

To install OpenOmics, run this command in your terminal:

    $ pip install openomics

This is the preferred method to install OpenOmics, as it will always install the most recent stable release.

If you don't have `pip` installed, it is recommended to install
the [Anaconda Python distribution](https://www.anaconda.com/products/individual) to get started with Python 3.7, 3.8, or
3.9 on either Mac OS, Linux, or Windows.


## From sources

The sources for OpenOmics can be downloaded from the [Github repo](https://github.com/BioMeCIS-Lab/OpenOmics).

You can either clone the public repository:

    $ git clone https://github.com/BioMeCIS-Lab/OpenOmics

Once you have a copy of the source, you can install it with:

    $ cd OpenOmics
    $ pip install ./ -e

The versioned release history can be found at [OpenOmics/releases](https://github.com/BioMeCIS-Lab/OpenOmics/releases).
# Contributing

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

## Types of Contributions
- Web development with the Dash for the OpenOmics dashboard webserver.
- Documentation standards for the OpenOmics Python API.
- Organize the library of genomics, functional ontologies, interactions, and sequence databases for variety of
  biological studies.
- Implement general purpose utilities for importing various fasta, gtf and sequencing files.

## Report Bugs

Report bugs at [openomics/issues](https://github.com/BioMeCIS-Lab/openomics/issues).

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

## Fix Bugs

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

## Implement Features

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

## Write Documentation

OpenOmics could always use more documentation, whether as part of the
[official OpenOmics docs](https://openomics.readthedocs.io/), in docstrings within the API, or even on the web in blog posts,
articles, and such.

If you'd like to help write RTD documentations, note:
- Documentation pages are written in markdown using [myst-parser](https://myst-parser.readthedocs.io/en/latest/index.html)
- The Sphinx theme used is [furo](https://pradyunsg.me/furo/)
- The autodoc package used is [sphinx-automodapi](https://sphinx-automodapi.readthedocs.io/en/latest/)

## Submit Feedback

The best way to send feedback is to file an issue at [openomics/issues](https://github.com/BioMeCIS-Lab/openomics/issues).

If you are proposing a feature:

- Explain in detail how it would work.
- Keep the scope as narrow as possible, to make it easier to implement.
- Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

## Get Started!

Ready to contribute? Here's how to set up `openomics` for local development.

1. Fork the `openomics` repo on GitHub.
2. Clone your fork locally and work on the develop branch:
```
$ git clone git@github.com:your_name_here/openomics.git
$ git checkout develop
```

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development:
```
$ mkvirtualenv openomics
$ cd openomics/
$ python setup.py develop
```

4. Create a branch for local development:
```
$ git checkout -b name-of-your-bugfix-or-feature
```
   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox:
```
$ flake8 openomics tests
$ python setup.py test or py.test $ tox
```

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub:
```
$ git add . $ git commit -m "Your detailed description of your changes."
$ git push develop name-of-your-bugfix-or-feature
```
7. Submit a pull request through the GitHub website to the develop branch. Once major features are tested, we can create
   another pull-request to the **master** branch.

## Pull Request Guidelines

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests. Run tests by with `pytest ./` and make sure tests are 100% passing.
2. If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a
   docstring, and add the feature to the list in [docs/history.md](https://github.com/BioMeCIS-Lab/OpenOmics/blob/master/docs/history.md).
3. The pull request should work for Python 3.6 or higher, and for PyPi. Check
   [Github Actions Tests](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml)
   and make sure that the tests pass for all supported Python versions and operating systems.

## Tips

To run the automated tests locally, run this at the root directory:

    pytest ./

```{hint}
To run a subset of tests:

    $ py.test tests.test_openomics

```

To run tests targeting various operating systems and Python versions, make a pull-request to the **master** branch which
will run as [Github Actions Tests](https://github.com/BioMeCIS-Lab/OpenOmics/actions/workflows/python-package.yml).

## Deploying

A reminder for the maintainers on how to deploy. Make sure all your changes are committed (including an entry in
HISTORY.rst). Then run:

    $ bumpversion patch # possible: major / minor / patch
    $ git push --tags

Github Actions will then deploy to PyPI if tests pass.

## Code of Conduct
Please note that the OpenOmics project is released with a Contributor Code of Conduct. By contributing to this project you agree to abide by its terms.

[openomics/issues]: https://github.com/BioMeCIS-Lab/openomics/issues
# Release history

## 0.9.0 (Future)
- Build web-app dashboard interface for importing user multiomics data files

## 0.8.6 (Pending)
- Changed to Github Actions CI from Travis CI
- Revamped openomics.readthedocs.io
- Fixed bugs from pyOpenSci reviews.

## 0.8.5 (2020-02-19)
- Importing GENCODE gtf files using dask now works with gzip compressed files.
- Improved coverage and performance of automated tests.

## 0.8.4 (2020-01-07)
- Enabled the support for Dask dataframes for ExpressionData and Annotation Dataset's. To use this feature, simply use
  the npartitions argument when instantiating an ExpressionData or Annotation/Sequence/Interaction Dataset.

## 0.7.2 (2019-09-01)
- Added compatibility for Python 2.7
- Refactored ClinicalData
- Built working documentations with Sphinx on readthedocs
- Added pytests for MultiOmicsData
- First successful build on Travis CI on Python 3.4-3.7, 2.7

## 0.6.0 (2019-08-31)
- First release on PyPI.
# Loading a multi-omics dataset

Suppose you have your own -omics dataset(s) and you'd like to load them. One of OpenOmics's primary goal is to
encapsulate the data import process with one line of code along with a few parameters. Given any processed single-omic
dataset, the library loads the data as a tabular structure where rows correspond to observation samples and columns
correspond to measurements of different biomolecules.

Import TCGA LUAD data included in tests dataset (preprocessed from TCGA-Assembler). It is located at [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD).

```{code-block} python
folder_path = "tests/data/TCGA_LUAD/"
```

Load the multiomics: Gene Expression, MicroRNA expression lncRNA expression, Copy Number Variation, Somatic Mutation, DNA Methylation, and Protein Expression data

```{code-block} python
from openomics import MessengerRNA, MicroRNA, LncRNA, SomaticMutation, Protein

# Load each expression dataframe
mRNA = MessengerRNA(data=folder_path + "LUAD__geneExp.txt",
                    transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="gene_name")
miRNA = MicroRNA(data=folder_path + "LUAD__miRNAExp__RPM.txt",
                 transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="transcript_name")
lncRNA = LncRNA(data=folder_path + "TCGA-rnaexpr.tsv",
                transpose=True, usecols="Gene_ID|TCGA", gene_index="Gene_ID", gene_level="gene_id")
som = SomaticMutation(data=folder_path + "LUAD__somaticMutation_geneLevel.txt",
                      transpose=True, usecols="GeneSymbol|TCGA", gene_index="gene_name")
pro = Protein(data=folder_path + "protein_RPPA.txt",
              transpose=True, usecols="GeneSymbol|TCGA", gene_index="GeneSymbol", gene_level="protein_name")

# Create an integrated MultiOmics dataset
luad_data = MultiOmics(cohort_name="LUAD")
luad_data.add_clinical_data(
    clinical=folder_path + "nationwidechildrens.org_clinical_patient_luad.txt")

luad_data.add_omic(mRNA)
luad_data.add_omic(miRNA)
luad_data.add_omic(lncRNA)
luad_data.add_omic(som)
luad_data.add_omic(pro)

luad_data.build_samples()
```

Each data is stored as a Pandas DataFrame. Below are all the data imported for TCGA LUAD. For each, the first number represents the number of samples, the second number is the number of features.

> PATIENTS (522, 5)
  SAMPLES (1160, 6)
  DRUGS (461, 4)
  MessengerRNA (576, 20472)
  SomaticMutation (587, 21070)
  MicroRNA (494, 1870)
  LncRNA (546, 12727)
  Protein (364, 154)

You may notice that in this dataset, the samples index (e.g. TCGA-XX-XXXX) across different omics does not match. It may
be necessary to change them to be 12 characters in total.

```python
lncRNA.expressions.index = lncRNA.expressions.index.str.slice(-12, )
miRNA.expressions.index = miRNA.expressions.index.str.slice(0, 12)
mRNA.expressions.index = mRNA.expressions.index.str.slice(0, 12)
som.expressions.index = som.expressions.index.str.slice(0, 12)
pro.expressions.index = pro.expressions.index.str.slice(0, 12)

luad_data.build_samples()
luad_data.samples
```
> Index(['TCGA-05-4244', 'TCGA-05-4249', 'TCGA-05-4250', 'TCGA-05-4382',
    'TCGA-05-4384', 'TCGA-05-4389', 'TCGA-05-4390', 'TCGA-05-4395',
    'TCGA-05-4396', 'TCGA-05-4397', ...
    'TCGA-NJ-A4YG', 'TCGA-NJ-A4YI', 'TCGA-NJ-A4YP', 'TCGA-NJ-A4YQ',
    'TCGA-NJ-A55A', 'TCGA-NJ-A55O', 'TCGA-NJ-A55R', 'TCGA-NJ-A7XG',
    'TCGA-O1-A52J', 'TCGA-S2-AA1A'], dtype='object', length=952)

## Load single omics expressions for MessengerRNA, MicroRNA, LncRNA

We instantiate the MessengerRNA, MicroRNA and LncRNA -omics expression data from `gtex.data`. Since the gene expression
were not seperated by RNA type, we use GENCODE and Ensembl gene annotations to filter the list of mRNA, miRNA, and
lncRNAs.

```{code-block} python
from openomics import MessengerRNA, MicroRNA, LncRNA

# Gene Expression
messengerRNA_id = gtex_transcripts_gene_id & pd.Index(gencode.data[gencode.data["gene_type"] == "protein_coding"]["gene_id"].unique())

messengerRNA = MessengerRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(messengerRNA_id)],
                           transpose=True, gene_index="gene_name", usecols=None, npartitions=4)

# MicroRNA expression
microRNA_id = pd.Index(ensembl.data[ensembl.data["gene_biotype"] == "miRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id

microRNA = MicroRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(microRNA_id)],
                   gene_index="gene_id", transpose=True, usecols=None, )

# LncRNA expression
lncRNA_id = pd.Index(gencode.data[gencode.data["gene_type"] == "lncRNA"]["gene_id"].unique()) & gtex_transcripts_gene_id
lncRNA = LncRNA(gtex_transcripts[gtex_transcripts["gene_id"].isin(lncRNA_id)],
               gene_index="gene_id", transpose=True, usecols=None, )
```

## Create a MultiOmics dataset

Now, we create a MultiOmics dataset object by combining the messengerRNA, microRNA, and lncRNA.

```{code-block} python
   from openomics import MultiOmics

   gtex_data = MultiOmics(cohort_name="GTEx Tissue Avg Expressions")

   gtex_data.add_omic(messengerRNA)
   gtex_data.add_omic(microRNA)
   gtex_data.add_omic(lncRNA)

   gtex_data.build_samples()
```

## Accessing clinical data
Each multi-omics and clinical data can be accessed through luad_data.data[], like:

```{code-block} python
luad_data.data["PATIENTS"]
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gender</th>
      <th>race</th>
      <th>histologic_subtype</th>
      <th>pathologic_stage</th>
    </tr>
    <tr>
      <th>bcr_patient_barcode</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage IV</td>
    </tr>
    <tr>
      <th>TCGA-05-4245</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4249</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4250</th>
      <td>FEMALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma- Not Otherwise Specified (...</td>
      <td>Stage III</td>
    </tr>
    <tr>
      <th>TCGA-05-4382</th>
      <td>MALE</td>
      <td>NaN</td>
      <td>Lung Adenocarcinoma Mixed Subtype</td>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>522 rows × 4 columns</p>
</div>


```{code-block} python
luad_data.data["MessengerRNA"]
```
<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>gene_name</th>
      <th>A1BG</th>
      <th>A1BG-AS1</th>
      <th>A1CF</th>
      <th>A2M</th>
      <th>A2ML1</th>
      <th>A4GALT</th>
      <th>A4GNT</th>
      <th>AAAS</th>
      <th>AACS</th>
      <th>AACSP1</th>
      <th>...</th>
      <th>ZXDA</th>
      <th>ZXDB</th>
      <th>ZXDC</th>
      <th>ZYG11A</th>
      <th>ZYG11B</th>
      <th>ZYX</th>
      <th>ZZEF1</th>
      <th>ZZZ3</th>
      <th>psiTPTE22</th>
      <th>tAKR</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244-01A</th>
      <td>4.756500</td>
      <td>5.239211</td>
      <td>0.000000</td>
      <td>13.265291</td>
      <td>0.431997</td>
      <td>7.043317</td>
      <td>1.033652</td>
      <td>9.348765</td>
      <td>9.652057</td>
      <td>0.763921</td>
      <td>...</td>
      <td>5.350285</td>
      <td>8.197321</td>
      <td>9.907260</td>
      <td>0.763921</td>
      <td>10.088859</td>
      <td>11.471139</td>
      <td>9.768648</td>
      <td>9.170597</td>
      <td>2.932118</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4249-01A</th>
      <td>6.920471</td>
      <td>7.056843</td>
      <td>0.402722</td>
      <td>14.650247</td>
      <td>1.383939</td>
      <td>9.178805</td>
      <td>0.717123</td>
      <td>9.241537</td>
      <td>9.967223</td>
      <td>0.000000</td>
      <td>...</td>
      <td>5.980428</td>
      <td>8.950001</td>
      <td>10.204971</td>
      <td>4.411650</td>
      <td>9.622978</td>
      <td>11.199826</td>
      <td>10.153700</td>
      <td>9.433116</td>
      <td>7.499637</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4250-01A</th>
      <td>5.696542</td>
      <td>6.136327</td>
      <td>0.000000</td>
      <td>14.048541</td>
      <td>0.000000</td>
      <td>8.481646</td>
      <td>0.996244</td>
      <td>9.203535</td>
      <td>9.560412</td>
      <td>0.733962</td>
      <td>...</td>
      <td>5.931168</td>
      <td>8.517334</td>
      <td>9.722642</td>
      <td>4.782796</td>
      <td>8.895339</td>
      <td>12.408981</td>
      <td>10.194168</td>
      <td>9.060342</td>
      <td>2.867956</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>TCGA-05-4382-01A</th>
      <td>7.198727</td>
      <td>6.809804</td>
      <td>0.000000</td>
      <td>14.509730</td>
      <td>2.532591</td>
      <td>9.117559</td>
      <td>1.657045</td>
      <td>9.251035</td>
      <td>10.078124</td>
      <td>1.860883</td>
      <td>...</td>
      <td>5.373036</td>
      <td>8.441914</td>
      <td>9.888267</td>
      <td>6.041142</td>
      <td>9.828389</td>
      <td>12.725186</td>
      <td>10.192589</td>
      <td>9.376841</td>
      <td>5.177029</td>
      <td>0.000000</td>
    </tr>
  </tbody>
</table>
<p>576 rows × 20472 columns</p>
</div>

## To match samples accross different multi-omics, use
```{code-block} python
luad_data.match_samples(modalities=["MicroRNA", "MessengerRNA"])
```

    Index(['TCGA-05-4384-01A', 'TCGA-05-4390-01A', 'TCGA-05-4396-01A',
           'TCGA-05-4405-01A', 'TCGA-05-4410-01A', 'TCGA-05-4415-01A',
           'TCGA-05-4417-01A', 'TCGA-05-4424-01A', 'TCGA-05-4425-01A',
           'TCGA-05-4427-01A',
           ...
           'TCGA-NJ-A4YG-01A', 'TCGA-NJ-A4YI-01A', 'TCGA-NJ-A4YP-01A',
           'TCGA-NJ-A4YQ-01A', 'TCGA-NJ-A55A-01A', 'TCGA-NJ-A55O-01A',
           'TCGA-NJ-A55R-01A', 'TCGA-NJ-A7XG-01A', 'TCGA-O1-A52J-01A',
           'TCGA-S2-AA1A-01A'],
          dtype='object', length=465)

# Creating biological interaction networks
# Preparing data for downstream analyses

## To prepare the data for classification

```python
X_multiomics, y = luad_data.load_data(omics="all", target=["pathologic_stage"], remove_duplicates=True)

print(X_multiomics['MessengerRNA'].shape,
      X_multiomics['MicroRNA'].shape,
      X_multiomics['LncRNA'].shape,
      y.shape)
```

> (338, 20472) (338, 1870) (338, 12727) (338, 1)


```python
print(y)
```


<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>pathologic_stage</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4390-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4405-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4410-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4417-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-4424-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4427-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-4433-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-05-5423-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5425-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5428-01A</th>
      <td>Stage II</td>
    </tr>
    <tr>
      <th>TCGA-05-5715-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-4631-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-7271-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-38-A44F-01A</th>
      <td>Stage I</td>
    </tr>
    <tr>
      <th>TCGA-44-2655-11A</th>
      <td>Stage I</td>
    </tr>
  </tbody>
</table>
<p>336 rows × 1 columns</p>
</div>



## Log2 transform the mRNA, microRNA, and lncRNA expression values


```python
def expression_val_transform(x):
    return np.log2(x+1)
X_multiomics['MessengerRNA'] = X_multiomics['MessengerRNA'].applymap(expression_val_transform)
X_multiomics['MicroRNA'] = X_multiomics['MicroRNA'].applymap(expression_val_transform)
# X_multiomics['LncRNA'] = X_multiomics['LncRNA'].applymap(expression_val_transform)
```

## Classification of Cancer Stage


```python
from sklearn import preprocessing
from sklearn import metrics
from sklearn.svm import SVC, LinearSVC
import sklearn.linear_model
from sklearn.model_selection import train_test_split

```


```python
binarizer = preprocessing.LabelEncoder()
binarizer.fit(y)
binarizer.transform(y)
```


    array([0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0,
           0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
           0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1,
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1,
           0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1,
           0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0,
           0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0,
           0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0,
           1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1,
           1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1,
           1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0,
           0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0,
           1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0])




```python
for omic in ["MessengerRNA", "MicroRNA"]:
    print(omic)
    scaler = sklearn.preprocessing.StandardScaler(copy=True, with_mean=True, with_std=False)
    scaler.fit(X_multiomics[omic])

    X_train, X_test, Y_train, Y_test = \
        train_test_split(X_multiomics[omic], y, test_size=0.3, random_state=np.random.randint(0, 10000), stratify=y)
    print(X_train.shape, X_test.shape)


    X_train = scaler.transform(X_train)

    model = LinearSVC(C=1e-2, penalty='l1', class_weight='balanced', dual=False, multi_class="ovr")
#     model = sklearn.linear_model.LogisticRegression(C=1e-0, penalty='l1', fit_intercept=False, class_weight="balanced")
#     model = SVC(C=1e0, kernel="rbf", class_weight="balanced", decision_function_shape="ovo")

    model.fit(X=X_train, y=Y_train)
    print("NONZERO", len(np.nonzero(model.coef_)[0]))
    print("Training accuracy", metrics.accuracy_score(model.predict(X_train), Y_train))
    print(metrics.classification_report(y_pred=model.predict(X_test), y_true=Y_test))

```

    MessengerRNA
    (254, 20472) (109, 20472)
    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support

        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34

    avg / total       0.47      0.69      0.56       109

    MicroRNA
    (254, 1870) (109, 1870)
    NONZERO 0
    Training accuracy 0.6929133858267716
                 precision    recall  f1-score   support

        Stage I       0.69      1.00      0.82        75
       Stage II       0.00      0.00      0.00        34

    avg / total       0.47      0.69      0.56       109
# Getting started

Welcome! This tutorial highlights the OpenOmics API’s core features; for in-depth details and conceptual guides, see the links within, or the documentation index which has links to use cases, and API reference sections.

## Loading a single-omics dataframe

Suppose you have a single-omics dataset and would like to load them as a dataframe.

As an example, we use the `TGCA` Lung Adenocarcinoma dataset from [tests/data/TCGA_LUAD](https://github.com/BioMeCIS-Lab/OpenOmics/tree/master/tests/data/TCGA_LUAD). Data tables are tab-delimited and have the following format:

| GeneSymbol | EntrezID  | TCGA-05-4244-01A-01R-1107-07 | TCGA-05-4249-01A-01R-1107-07 | ...  |
| ---------- | --------- | ---------------------------- | ---------------------------- | ---- |
| A1BG       | 100133144 | 10.8123                      | 3.7927                       | ...  |
| ⋮ | ⋮ | ⋮ | ⋮ |

Depending on whether your data table is stored locally as a single file, splitted into multiple files, or was already a dataframe, you can load it using the class {class}`openomics.transcriptomics.Expression` or any of its subclasses.

````{tab} From a single file
If the dataset is a local file in a tabular format, OpenOmics can help you load them to Pandas dataframe.

```{code-block} python
from openomics.multiomics import MessengerRNA

mrna = MessengerRNA(
    data="https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/LUAD__geneExp.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA", # A regex that matches all column name with either "GeneSymbol" or "TCGA substring
    gene_index="GeneSymbol", # This column contains the gene index
    )
```

One thing to pay attention is that the raw data file given is column-oriented where columns corresponds to samples, so we have use the argument `transpose=True` to convert to row-oriented.
> MessengerRNA (576, 20472)
````

````{tab} From multiple files (glob)
If your dataset is large, it may be broken up into multiple files with a similar file name prefix/suffix. Assuming all the files have similar tabular format, OpenOmics can load all files and contruct an integrated data table using the memory-efficient Dask dataframe.

```python
from openomics.multiomics import MessengerRNA

mrna = MessengerRNA("TCGA_LUAD/LUAD__*", # Files must be stored locally
                    transpose=True,
                    usecols="GeneSymbol|TCGA",
                    gene_index="GeneSymbol")
```

> INFO: Files matched: ['LUAD__miRNAExp__RPM.txt', 'LUAD__protein_RPPA.txt', 'LUAD__geneExp.txt']
````

````{tab} From DataFrame
If your workflow already produced a dataframe, you can encapsulate it directly with {class}`openomics.transcriptomics.Expression`.

```python
import pandas as pd
import numpy as np
from openomics.multiomics import MessengerRNA

# A random dataframe of microRNA gene_id's.
df = pd.DataFrame(data={"ENSG00000194717": np.random.rand(5),
                        "ENSG00000198973": np.random.rand(5),
                        "ENSG00000198974": np.random.rand(5),
                        "ENSG00000198975": np.random.rand(5),
                        "ENSG00000198976": np.random.rand(5),
                        "ENSG00000198982": np.random.rand(5),
                        "ENSG00000198983": np.random.rand(5)},
                  index=range(5))
mrna = MessengerRNA(df, transpose=False, sample_level="sample_id")
```
````
---
To access the {class}`DataFrame`, simply use {obj}`mrna.expressions`:
```python
print(mrna.expressions)
```
<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>GeneSymbol</th>
      <th>A1BG</th>
      <th>A1BG-AS1</th>
      <th>A1CF</th>
      <th>A2M</th>
    </tr>
    <tr>
      <th>sample_index</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244-01A-01R-1107-07</th>
      <td>26.0302</td>
      <td>36.7711</td>
      <td>0.000</td>
      <td>9844.7858</td>
    </tr>
    <tr>
      <th>TCGA-05-4249-01A-01R-1107-07</th>
      <td>120.1349</td>
      <td>132.1439</td>
      <td>0.322</td>
      <td>25712.6617</td>
    </tr>
  </tbody>
</table>
</div>

<br/>

## Creating a multi-omics dataset

With multiple single-omics, each with different sets of genes and samples, you can use the {class}`openomics.MultiOmics` to integrate them.

```{code-block} python
from openomics.multiomics import MultiOmics, MessengerRNA, MicroRNA, LncRNA, SomaticMutation, Protein

path = "https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/"

# Load each expression dataframe
mRNA = MessengerRNA(path+"LUAD__geneExp.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
miRNA = MicroRNA(path+"LUAD__miRNAExp__RPM.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
lncRNA = LncRNA(path+"TCGA-rnaexpr.tsv",
    transpose=True,
    usecols="Gene_ID|TCGA",
    gene_index="Gene_ID")
som = SomaticMutation(path+"LUAD__somaticMutation_geneLevel.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")
pro = Protein(path+"protein_RPPA.txt",
    transpose=True,
    usecols="GeneSymbol|TCGA",
    gene_index="GeneSymbol")

# Create an integrated MultiOmics dataset
luad_data = MultiOmics(cohort_name="LUAD", omics_data=[mRNA, mRNA, lncRNA, som, pro])
# You can also add individual -omics one at a time `luad_data.add_omic(mRNA)`

luad_data.build_samples()
```
The `luad_data` is a {class}`MultiOmics` object builds the samples list from all the samples given in each -omics data.

> MessengerRNA (576, 20472)
> MicroRNA (494, 1870)
> LncRNA (546, 12727)
> SomaticMutation (587, 21070)
> Protein (364, 154)

To access individual -omics data within `luad_data`, such as the {obj}`mRNA`, simply use the `.` accessor with the class name {class}`MessengerRNA`:
```python
luad_data.MessengerRNA
# or
luad_data.data["MessengerRNA"]
```

<br/>

## Adding clinical data as sample attributes

When sample attributes are provided for the study cohort, load it as a data table with the {class}`openomics.clinical.ClinicalData`, then add it to the {class}`openomics.multiomics.MultiOmics` dataset to enable querying for subsets of samples across the multi-omics.

```python
from openomics import ClinicalData

clinical = ClinicalData(
    "https://raw.githubusercontent.com/BioMeCIS-Lab/OpenOmics/master/tests/data/TCGA_LUAD/nationwidechildrens.org_clinical_patient_luad.txt",
    patient_index="bcr_patient_barcode")

luad_data.add_clinical_data(clinical)

luad_data.clinical.patient
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>bcr_patient_uuid</th>
      <th>form_completion_date</th>
      <th>histologic_diagnosis</th>
      <th>prospective_collection</th>
      <th>retrospective_collection</th>
    </tr>
    <tr>
      <th>bcr_patient_barcode</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>TCGA-05-4244</th>
      <td>34040b83-7e8a-4264-a551-b16621843e28</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
    <tr>
      <th>TCGA-05-4245</th>
      <td>03d09c05-49ab-4ba6-a8d7-e7ccf71fafd2</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
    <tr>
      <th>TCGA-05-4249</th>
      <td>4addf05f-3668-4b3f-a17f-c0227329ca52</td>
      <td>2010-7-22</td>
      <td>Lung Adenocarcinoma</td>
      <td>NO</td>
      <td>YES</td>
    </tr>
  </tbody>
</table>
</div>

Note that in the clinical data table, `bcr_patient_barcode` is the column with `TCGA-XX-XXXX` patient IDs, which matches
that of the `sample_index` index column in the `mrna.expressions` dataframe.

````{note}
In our `TCGA_LUAD` example, mismatches in the `bcr_patient_barcode` sample index of clinical dataframe may happen because the `sample_index` in `mRNA` may have a longer form `TCGA-XX-XXXX-XXX-XXX-XXXX-XX` that contain the samples number and aliquot ID's. To make them match, you can modify the index strings on-the-fly using the [Pandas's extensible API](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.str.slice.html):
```python
mRNA.expressions.index = mRNA.expressions.index.str.slice(0, 12) # Selects only the first 12 characters
```
````

<br/>

## Import an external database

Next, we may want to annotate the genes list in our RNA-seq expression dataset with genomics annotation. To do so, we'd need to download annotations from the [GENCODE database](https://www.gencodegenes.org/), preprocess annotation files into a dataframe, and then match them with the genes in our dataset.

OpenOmics provides a simple, hassle-free API to download the GENCODE annotation files via FTP with these steps:
1. First, provide the base `path` of the FTP download server - usually found in the direct download link on GENCODE's website. Most of the time, selecting the right base `path` allows you to specify the specific species, genome assembly, and database version for your study.
2. Secondly, use the `file_resources` dict parameter to select the data files and the file paths required to construct the annotation dataframe. For each entry in the `file_resources`, the key is the alias of the file required, and the value is the filename with the FTP base `path`.

   For example, the entry `{"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz"}` indicates the GENCODE class to preprocess a `.gtf` file with the alias `"long_noncoding_RNAs.gtf"`, downloaded from the FTP path `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz`

   To see which file alias keys are required to construct a dataframe, refer to the docstring in {class}`openomics.database.sequence.GENCODE`.

```python
from openomics.database import GENCODE

gencode = GENCODE(
    path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
    file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                    "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                    "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz", # lncRNA sequences
                    "transcripts.fa": "gencode.v32.transcripts.fa.gz" # mRNA sequences
                    },
    npartitions=0, # if > 1, then use Dask partition the dataframe and leverage out-of-core multiprocessing
)
```
To access the attributes constructed from the combination of annotations `long_noncoding_RNAs.gtf` and `
basic.annotation.gtf`, use:

```python
gencode.data
```

<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }

</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>gene_id</th>
      <th>gene_name</th>
      <th>index</th>
      <th>seqname</th>
      <th>source</th>
      <th>feature</th>
      <th>start</th>
      <th>end</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>ENSG00000243485</td>
      <td>MIR1302-2HG</td>
      <td>0</td>
      <td>chr1</td>
      <td>HAVANA</td>
      <td>gene</td>
      <td>29554</td>
      <td>31109</td>
    </tr>
    <tr>
      <th>1</th>
      <td>ENSG00000243485</td>
      <td>MIR1302-2HG</td>
      <td>1</td>
      <td>chr1</td>
      <td>HAVANA</td>
      <td>transcript</td>
      <td>29554</td>
      <td>31097</td>
    </tr>
  </tbody>
</table>
</div>


<br/>

## Annotate your expression dataset with attributes
With the annotation database, you can perform a join operation to add gene attributes to your {class}`openomics.transcriptomics.Expression` dataset. To annotate attributes for the `gene_id` list `mRNA.expression`, you must first select the corresponding column in `gencode.data` with matching `gene_id` keys. The following are code snippets for a variety of database types.

````{tab} Genomics attributes
```python
luad_data.MessengerRNA.annotate_attributes(gencode,
    index="gene_id",
    columns=['gene_name', 'start', 'end', 'strand'] # Add these columns to the .annotations dataframe
    )
```

````

````{tab} Sequences
```python
luad_data.MessengerRNA.annotate_sequences(gencode,
    index="gene_name",
    agg_sequences="all", # Collect all sequences with the gene_name into a list
    )
```
````

````{tab} Disease Associations
```python
from openomics.database.disease import DisGeNet
disgenet = DisGeNet(path="https://www.disgenet.org/static/disgenet_ap1/files/downloads/", curated=True)

luad_data.MessengerRNA.annotate_diseases(disgenet, index="gene_name")
```
````

---
To view the resulting annotations dataframe, use:
```python
luad_data.MessengerRNA.annotations
```


For more detailed guide, refer to the [annotation interfaces API](../modules/openomics.annotate.md).
# External annotation databases

## Import GENCODE human release 32

Next, we can annotate the genes in our GTEx expression dataset with genomics annotation from GENCODE. In this example,
we use the URL path prefix "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/" which specifies the
species and release version. We also pass a dictionary `file_resources`, with key-value pairs where the key is name of
file and value is the suffix of the file download URL.

For example, file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz"} will download file
located at <ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.long_noncoding_RNAs.gtf.gz>
to process the `long_noncoding_RNAs.gtf` file.

Here, we loads both "long_noncoding_RNAs.gtf" and "basic.annotation.gtf" which builds a dataframe of combined
annotations for both lncRNAs and mRNAs. You can specify different annotation files options from GENCODE by modifying
the `file_resources` dict argument.

```{code-block} python
from openomics.database import GENCODE, EnsemblGenes

gencode = GENCODE(path="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/",
                 file_resources={"long_noncoding_RNAs.gtf": "gencode.v32.long_noncoding_RNAs.gtf.gz",
                                 "basic.annotation.gtf": "gencode.v32.basic.annotation.gtf.gz",
                                 "lncRNA_transcripts.fa": "gencode.v32.lncRNA_transcripts.fa.gz",
                                 "transcripts.fa": "gencode.v32.transcripts.fa.gz"},
                 remove_version_num=True)

# We also loads Ensembl genes to get list of miRNA gene IDs
ensembl = EnsemblGenes(biomart='hsapiens_gene_ensembl', npartitions=8, )
```

## Setting the cache download directory
The package `astropy` is used to automatically cache downloaded files. It defaults to saving the files at
`~/.astropy/cache/`, where the cached content is retrieved given the matching URL. To change the path for the cache download file, run:

```python
import openomics

openomics.set_cache_dir(path="PATH/OF/YOUR/CHOICE/")
```

```{note}
Note that this setting doesn't persist across different programming sessions. Ideally, the cache dir should be in one location to minimize automatic FTP downloads, which may cause unnecessary stress on the database server.
```

# Frequently Asked Questions (FAQ)

### How to change the download directory where database files are cached?

You can make the setting at `openomics.set_cache_dir(path="PATH/OF/YOUR/CHOICE/")` which will change the default
cache location. However, this setting doesn't persist across user sessions. You can also change the default cache
directory at a `~/.openomics/conf.json` to persist the setting for each user.


# Interaction databases

```{eval-rst}
.. automodapi:: openomics.database.interaction
```
# Genomic databases

```{eval-rst}
.. automodapi:: openomics.database.annotation
```
# Annotation interfaces

```{eval-rst}
.. automodapi:: openomics.database.base
    :no-heading:
    :no-inheritance-diagram:
    :skip: ABC, Dict, List, Tuple
```
# Sequence databases

```{eval-rst}
.. automodapi:: openomics.database.sequence
```
# MultiOmics data

```{eval-rst}
.. automodapi:: openomics
    :classes-only:
    :skip: set_backend, set_cache_dir
```
# Utilities

```{eval-rst}
.. automodapi:: openomics.utils
```
# Disease databases

```{eval-rst}
.. automodapi:: openomics.database.disease
```
# Ontology databases

```{eval-rst}
.. automodapi:: openomics.database.ontology
```
{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module Attributes

   .. autosummary::
      :toctree:                                          <-- add this line
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:                                          <-- add this line
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:                                          <-- add this line
      :template: custom-class-template.rst               <-- add this line
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:                                          <-- add this line
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst                 <-- add this line
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}
{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:                                    <-- add at least this line
   :show-inheritance:                           <-- plus I want to show inheritance...
   :inherited-members:                          <-- ...and inherited members too

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
MirBase
=======

.. currentmodule:: openomics.database.sequence

.. autoclass:: MirBase
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MirBase.get_sequences
      ~MirBase.load_dataframe
      ~MirBase.read_fasta

   .. rubric:: Methods Documentation

   .. automethod:: get_sequences
   .. automethod:: load_dataframe
   .. automethod:: read_fasta
LncRNA
======

.. currentmodule:: openomics

.. autoclass:: LncRNA
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~LncRNA.name

   .. rubric:: Methods Documentation

   .. automethod:: name
ProteinAtlas
============

.. currentmodule:: openomics.database.annotation

.. autoclass:: ProteinAtlas
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~ProteinAtlas.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~ProteinAtlas.get_expressions
      ~ProteinAtlas.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: get_expressions
   .. automethod:: load_dataframe
SomaticMutation
===============

.. currentmodule:: openomics

.. autoclass:: SomaticMutation
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~SomaticMutation.name

   .. rubric:: Methods Documentation

   .. automethod:: name
lncRNome
========

.. currentmodule:: openomics.database.interaction

.. autoclass:: lncRNome
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~lncRNome.load_dataframe
      ~lncRNome.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
   .. automethod:: load_network
write_taxonomy
==============

.. currentmodule:: openomics.database.ontology

.. autofunction:: write_taxonomy
GeneMania
=========

.. currentmodule:: openomics.database.interaction

.. autoclass:: GeneMania
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~GeneMania.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
MicroRNA
========

.. currentmodule:: openomics

.. autoclass:: MicroRNA
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MicroRNA.name

   .. rubric:: Methods Documentation

   .. automethod:: name
Interactions
============

.. currentmodule:: openomics.database.interaction

.. autoclass:: Interactions
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~Interactions.filter_values
      ~Interactions.get_interactions
      ~Interactions.info
      ~Interactions.load_network
      ~Interactions.name

   .. rubric:: Methods Documentation

   .. automethod:: filter_values
   .. automethod:: get_interactions
   .. automethod:: info
   .. automethod:: load_network
   .. automethod:: name
read_gtf
========

.. currentmodule:: openomics.utils

.. autofunction:: read_gtf
GeneOntology
============

.. currentmodule:: openomics.database.ontology

.. autoclass:: GeneOntology
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~GeneOntology.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~GeneOntology.add_predecessor_terms
      ~GeneOntology.filter_network
      ~GeneOntology.get_predecessor_terms
      ~GeneOntology.info
      ~GeneOntology.load_dataframe
      ~GeneOntology.load_network

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: add_predecessor_terms
   .. automethod:: filter_network
   .. automethod:: get_predecessor_terms
   .. automethod:: info
   .. automethod:: load_dataframe
   .. automethod:: load_network
LncBase
=======

.. currentmodule:: openomics.database.interaction

.. autoclass:: LncBase
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~LncBase.get_rename_dict
      ~LncBase.load_network

   .. rubric:: Methods Documentation

   .. automethod:: get_rename_dict
   .. automethod:: load_network
Database
========

.. currentmodule:: openomics.database.base

.. autoclass:: Database
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~Database.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~Database.close
      ~Database.get_annotations
      ~Database.get_expressions
      ~Database.info
      ~Database.list_databases
      ~Database.load_dataframe
      ~Database.name
      ~Database.validate_file_resources

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: close
   .. automethod:: get_annotations
   .. automethod:: get_expressions
   .. automethod:: info
   .. automethod:: list_databases
   .. automethod:: load_dataframe
   .. automethod:: name
   .. automethod:: validate_file_resources
MiRTarBase
==========

.. currentmodule:: openomics.database.interaction

.. autoclass:: MiRTarBase
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MiRTarBase.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
MalaCards
=========

.. currentmodule:: openomics.database.disease

.. autoclass:: MalaCards
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~MalaCards.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~MalaCards.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
ClinicalData
============

.. currentmodule:: openomics

.. autoclass:: ClinicalData
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~ClinicalData.pathologic_stage_map

   .. rubric:: Methods Summary

   .. autosummary::

      ~ClinicalData.add_biospecimen_data
      ~ClinicalData.add_drug_response_data
      ~ClinicalData.build_clinical_samples
      ~ClinicalData.get_patient_barcodes
      ~ClinicalData.get_sample_barcodes
      ~ClinicalData.name

   .. rubric:: Attributes Documentation

   .. autoattribute:: pathologic_stage_map

   .. rubric:: Methods Documentation

   .. automethod:: add_biospecimen_data
   .. automethod:: add_drug_response_data
   .. automethod:: build_clinical_samples
   .. automethod:: get_patient_barcodes
   .. automethod:: get_sample_barcodes
   .. automethod:: name
DNAMethylation
==============

.. currentmodule:: openomics

.. autoclass:: DNAMethylation
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~DNAMethylation.name

   .. rubric:: Methods Documentation

   .. automethod:: name
TargetScan
==========

.. currentmodule:: openomics.database.interaction

.. autoclass:: TargetScan
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~TargetScan.load_network
      ~TargetScan.process_interactions_table
      ~TargetScan.process_miR_family_info_table

   .. rubric:: Methods Documentation

   .. automethod:: load_network
   .. automethod:: process_interactions_table
   .. automethod:: process_miR_family_info_table
EnsemblGeneSequences
====================

.. currentmodule:: openomics.database.annotation

.. autoclass:: EnsemblGeneSequences
   :show-inheritance:
MultiOmics
==========

.. currentmodule:: openomics

.. autoclass:: MultiOmics
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MultiOmics.add_clinical_data
      ~MultiOmics.add_omic
      ~MultiOmics.annotate_patients
      ~MultiOmics.build_samples
      ~MultiOmics.get_omics_list
      ~MultiOmics.get_patients_clinical
      ~MultiOmics.load_data
      ~MultiOmics.match_samples
      ~MultiOmics.print_sample_sizes
      ~MultiOmics.remove_duplicate_genes

   .. rubric:: Methods Documentation

   .. automethod:: add_clinical_data
   .. automethod:: add_omic
   .. automethod:: annotate_patients
   .. automethod:: build_samples
   .. automethod:: get_omics_list
   .. automethod:: get_patients_clinical
   .. automethod:: load_data
   .. automethod:: match_samples
   .. automethod:: print_sample_sizes
   .. automethod:: remove_duplicate_genes
STRING
======

.. currentmodule:: openomics.database.interaction

.. autoclass:: STRING
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~STRING.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~STRING.get_sequences
      ~STRING.load_network

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: get_sequences
   .. automethod:: load_network
DisGeNet
========

.. currentmodule:: openomics.database.disease

.. autoclass:: DisGeNet
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~DisGeNet.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~DisGeNet.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
set_cache_dir
=============

.. currentmodule:: openomics

.. autofunction:: set_cache_dir
SequenceDatabase
================

.. currentmodule:: openomics.database.sequence

.. autoclass:: SequenceDatabase
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~SequenceDatabase.get_aggregator
      ~SequenceDatabase.get_sequences
      ~SequenceDatabase.read_fasta

   .. rubric:: Methods Documentation

   .. automethod:: get_aggregator
   .. automethod:: get_sequences
   .. automethod:: read_fasta
NPInter
=======

.. currentmodule:: openomics.database.interaction

.. autoclass:: NPInter
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~NPInter.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
LncReg
======

.. currentmodule:: openomics.database.interaction

.. autoclass:: LncReg
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~LncReg.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
RNAcentral
==========

.. currentmodule:: openomics.database.annotation

.. autoclass:: RNAcentral
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~RNAcentral.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~RNAcentral.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
lncRInter
=========

.. currentmodule:: openomics.database.interaction

.. autoclass:: lncRInter
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~lncRInter.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
flatten_list
============

.. currentmodule:: openomics.database.ontology

.. autofunction:: flatten_list
dfs_path
========

.. currentmodule:: openomics.database.ontology

.. autofunction:: dfs_path
GTEx
====

.. currentmodule:: openomics.database.annotation

.. autoclass:: GTEx
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~GTEx.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~GTEx.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
OMIM
====

.. currentmodule:: openomics.database.disease

.. autoclass:: OMIM
   :show-inheritance:
Annotatable
===========

.. currentmodule:: openomics.database.base

.. autoclass:: Annotatable
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~Annotatable.DISEASE_ASSOCIATIONS_COL
      ~Annotatable.SEQUENCE_COL_NAME

   .. rubric:: Methods Summary

   .. autosummary::

      ~Annotatable.annotate_attributes
      ~Annotatable.annotate_diseases
      ~Annotatable.annotate_expressions
      ~Annotatable.annotate_interactions
      ~Annotatable.annotate_sequences
      ~Annotatable.get_annotation_expressions
      ~Annotatable.get_annotations
      ~Annotatable.get_rename_dict
      ~Annotatable.initialize_annotations
      ~Annotatable.set_index

   .. rubric:: Attributes Documentation

   .. autoattribute:: DISEASE_ASSOCIATIONS_COL
   .. autoattribute:: SEQUENCE_COL_NAME

   .. rubric:: Methods Documentation

   .. automethod:: annotate_attributes
   .. automethod:: annotate_diseases
   .. automethod:: annotate_expressions
   .. automethod:: annotate_interactions
   .. automethod:: annotate_sequences
   .. automethod:: get_annotation_expressions
   .. automethod:: get_annotations
   .. automethod:: get_rename_dict
   .. automethod:: initialize_annotations
   .. automethod:: set_index
BioGRID
=======

.. currentmodule:: openomics.database.interaction

.. autoclass:: BioGRID
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~BioGRID.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
MessengerRNA
============

.. currentmodule:: openomics

.. autoclass:: MessengerRNA
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MessengerRNA.name

   .. rubric:: Methods Documentation

   .. automethod:: name
filter_dfs_paths
================

.. currentmodule:: openomics.database.ontology

.. autofunction:: filter_dfs_paths
GENCODE
=======

.. currentmodule:: openomics.database.sequence

.. autoclass:: GENCODE
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~GENCODE.get_rename_dict
      ~GENCODE.get_sequences
      ~GENCODE.load_dataframe
      ~GENCODE.read_fasta

   .. rubric:: Methods Documentation

   .. automethod:: get_rename_dict
   .. automethod:: get_sequences
   .. automethod:: load_dataframe
   .. automethod:: read_fasta
traverse_predecessors
=====================

.. currentmodule:: openomics.database.ontology

.. autofunction:: traverse_predecessors
Expression
==========

.. currentmodule:: openomics

.. autoclass:: Expression
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~Expression.features
      ~Expression.gene_index
      ~Expression.samples

   .. rubric:: Methods Summary

   .. autosummary::

      ~Expression.drop_genes
      ~Expression.drop_samples
      ~Expression.get_genes_list
      ~Expression.get_samples_list
      ~Expression.load_dataframe
      ~Expression.load_dataframe_glob
      ~Expression.name
      ~Expression.preprocess_table
      ~Expression.set_genes_index

   .. rubric:: Attributes Documentation

   .. autoattribute:: features
   .. autoattribute:: gene_index
   .. autoattribute:: samples

   .. rubric:: Methods Documentation

   .. automethod:: drop_genes
   .. automethod:: drop_samples
   .. automethod:: get_genes_list
   .. automethod:: get_samples_list
   .. automethod:: load_dataframe
   .. automethod:: load_dataframe_glob
   .. automethod:: name
   .. automethod:: preprocess_table
   .. automethod:: set_genes_index
set_backend
===========

.. currentmodule:: openomics

.. autofunction:: set_backend
flatten
=======

.. currentmodule:: openomics.database.ontology

.. autofunction:: flatten
Ontology
========

.. currentmodule:: openomics.database.ontology

.. autoclass:: Ontology
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~Ontology.DELIM

   .. rubric:: Methods Summary

   .. autosummary::

      ~Ontology.filter_annotation
      ~Ontology.filter_network
      ~Ontology.get_adjacency_matrix
      ~Ontology.get_child_nodes
      ~Ontology.get_dfs_paths
      ~Ontology.get_node_color
      ~Ontology.get_root_nodes
      ~Ontology.load_network
      ~Ontology.remove_predecessor_terms

   .. rubric:: Attributes Documentation

   .. autoattribute:: DELIM

   .. rubric:: Methods Documentation

   .. automethod:: filter_annotation
   .. automethod:: filter_network
   .. automethod:: get_adjacency_matrix
   .. automethod:: get_child_nodes
   .. automethod:: get_dfs_paths
   .. automethod:: get_node_color
   .. automethod:: get_root_nodes
   .. automethod:: load_network
   .. automethod:: remove_predecessor_terms
get_pkg_data_filename
=====================

.. currentmodule:: openomics.utils

.. autofunction:: get_pkg_data_filename
EnsemblGenes
============

.. currentmodule:: openomics.database.annotation

.. autoclass:: EnsemblGenes
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~EnsemblGenes.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~EnsemblGenes.get_functional_annotations
      ~EnsemblGenes.get_rename_dict
      ~EnsemblGenes.load_data

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: get_functional_annotations
   .. automethod:: get_rename_dict
   .. automethod:: load_data
CopyNumberVariation
===================

.. currentmodule:: openomics

.. autoclass:: CopyNumberVariation
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~CopyNumberVariation.name

   .. rubric:: Methods Documentation

   .. automethod:: name
BioMartManager
==============

.. currentmodule:: openomics.database.annotation

.. autoclass:: BioMartManager
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~BioMartManager.cache_dataset
      ~BioMartManager.query_biomart
      ~BioMartManager.retrieve_dataset

   .. rubric:: Methods Documentation

   .. automethod:: cache_dataset
   .. automethod:: query_biomart
   .. automethod:: retrieve_dataset
DiseaseAssociation
==================

.. currentmodule:: openomics.database.disease

.. autoclass:: DiseaseAssociation
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~DiseaseAssociation.get_disease_assocs

   .. rubric:: Methods Documentation

   .. automethod:: get_disease_assocs
EnsemblSomaticVariation
=======================

.. currentmodule:: openomics.database.annotation

.. autoclass:: EnsemblSomaticVariation
   :show-inheritance:
EnsemblTranscriptSequences
==========================

.. currentmodule:: openomics.database.annotation

.. autoclass:: EnsemblTranscriptSequences
   :show-inheritance:
HMDD
====

.. currentmodule:: openomics.database.disease

.. autoclass:: HMDD
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~HMDD.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~HMDD.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
LncRNADisease
=============

.. currentmodule:: openomics.database.disease

.. autoclass:: LncRNADisease
   :show-inheritance:

   .. rubric:: Attributes Summary

   .. autosummary::

      ~LncRNADisease.COLUMNS_RENAME_DICT

   .. rubric:: Methods Summary

   .. autosummary::

      ~LncRNADisease.load_dataframe

   .. rubric:: Attributes Documentation

   .. autoattribute:: COLUMNS_RENAME_DICT

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
Protein
=======

.. currentmodule:: openomics

.. autoclass:: Protein
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~Protein.name
      ~Protein.process_HPRD_PPI_network

   .. rubric:: Methods Documentation

   .. automethod:: name
   .. automethod:: process_HPRD_PPI_network
EnsemblSNP
==========

.. currentmodule:: openomics.database.annotation

.. autoclass:: EnsemblSNP
   :show-inheritance:
StarBase
========

.. currentmodule:: openomics.database.interaction

.. autoclass:: StarBase
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~StarBase.load_network

   .. rubric:: Methods Documentation

   .. automethod:: load_network
NONCODE
=======

.. currentmodule:: openomics.database.annotation

.. autoclass:: NONCODE
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~NONCODE.load_dataframe

   .. rubric:: Methods Documentation

   .. automethod:: load_dataframe
TANRIC
======

.. currentmodule:: openomics.database.annotation

.. autoclass:: TANRIC
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~TANRIC.get_expressions
      ~TANRIC.load_dataframe

   .. rubric:: Methods Documentation

   .. automethod:: get_expressions
   .. automethod:: load_dataframe
HumanPhenotypeOntology
======================

.. currentmodule:: openomics.database.ontology

.. autoclass:: HumanPhenotypeOntology
   :show-inheritance:
MultiOmics
==========

.. currentmodule:: openomics.multiomics

.. autoclass:: MultiOmics
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~MultiOmics.add_clinical_data
      ~MultiOmics.add_omic
      ~MultiOmics.annotate_samples
      ~MultiOmics.build_samples
      ~MultiOmics.get_omics_list
      ~MultiOmics.get_sample_attributes
      ~MultiOmics.load_data
      ~MultiOmics.match_samples
      ~MultiOmics.print_sample_sizes
      ~MultiOmics.remove_duplicate_genes

   .. rubric:: Methods Documentation

   .. automethod:: add_clinical_data
   .. automethod:: add_omic
   .. automethod:: annotate_samples
   .. automethod:: build_samples
   .. automethod:: get_omics_list
   .. automethod:: get_sample_attributes
   .. automethod:: load_data
   .. automethod:: match_samples
   .. automethod:: print_sample_sizes
   .. automethod:: remove_duplicate_genes
LncRNA2Target
=============

.. currentmodule:: openomics.database.interaction

.. autoclass:: LncRNA2Target
   :show-inheritance:

   .. rubric:: Methods Summary

   .. autosummary::

      ~LncRNA2Target.load_network
      ~LncRNA2Target.load_network_high_throughput
      ~LncRNA2Target.load_network_low_throughput

   .. rubric:: Methods Documentation

   .. automethod:: load_network
   .. automethod:: load_network_high_throughput
   .. automethod:: load_network_low_throughput
