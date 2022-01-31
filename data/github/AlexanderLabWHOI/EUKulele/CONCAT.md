# EUKulele

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master)](https://travis-ci.com/AlexanderLabWHOI/EUKulele)
[![Coverage Status](https://coveralls.io/repos/github/AlexanderLabWHOI/EUKulele/badge.svg?branch=master)](https://coveralls.io/github/AlexanderLabWHOI/EUKulele?branch=master)
[![Documentation Status](https://readthedocs.org/projects/eukulele/badge/?version=latest)](https://eukulele.readthedocs.io/en/latest/?badge=latest)
[![Read the Docs](https://img.shields.io/badge/read-the%20docs-green)](https://eukulele.readthedocs.io/en/latest/)

<p align="center">
  <img width = "500" src="/.infrastructure/eukulele-logo.png"/>
</p>

## Formalizing environmental eukaryotic taxonomic assignment

### About EUKulele

`EUKulele` is a Python program for taxonomic annotation of microbes in metatranscriptomic and metagenomic samples, with special emphasis on eukaryote discovery. `EUKulele` can be downloaded from [PyPI](https://pypi.org/), or it may be downloaded via `conda` and used as a command-line program. The software includes four major features:
- Database setup and formatting
- Database creation, alignment, and taxonomic estimation
- Assessment of the BUSCO completeness of subsets of contigs at each taxonomic level
- Assessment of taxonomic classification using only BUSCO-identified core eukaryotic genes

### Prerequisites for running EUKulele

In principle, there are two prerequisites for running the software:
1. Metagenomic or metatranscriptomic sample files (unless using the provided sample data)
2. A database to align the contigs from the metagenome/metatranscriptome to

Three databases are supported by default from within `EUKulele`, and may be downloaded and formatted automatically if the user chooses (or if another reference directory is not specified/does not exist):
- [PhyloDB](https://drive.google.com/drive/u/0/folders/0B-BsLZUMHrDQfldGeDRIUHNZMEREY0g3ekpEZFhrTDlQSjQtbm5heC1QX2V6TUxBeFlOejQ)
- [EukProt](https://figshare.com/articles/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881/2)
- [MMETSP](https://zenodo.org/record/1212585#.Xw3PoJNKhTZ)

### Basic usage

If installed either with pip or `conda`, `EUKulele` can be invoked via::

```
EUKulele <arguments>
```
    
Where the minimal command would be

```
EUKulele --mets_or_mags <choice of data type> --sample_dir <where samples are located>
```

See the [documentation](https://eukulele.readthedocs.io/en/latest/) for further details.

### Community guidelines 

#### How to contribute to `EUKulele`
If you are interested in modifying `EUKulele`, you may fork the project for your own use, as detailed in the [MIT License](https://github.com/AlexanderLabWHOI/EUKulele/blob/master/LICENSE) we have adopted for the project. In order to contribute, please contact the developers via Arianna Krinos (akrinos (at) mit (dot) edu) after making the desired changes, after which a pull request may be submitted. 

#### Submitting an issue
If you have any suggestions for feature additions or any problems with the software that you would like addressed with the development community, please submit an issue on the [Issues tab](https://github.com/AlexanderLabWHOI/EUKulele/issues) of the project `GitHub` repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

#### Asking for help
If you have questions about how to use `EUKulele`, or would like to seek out collaborations related to this project, you may contact Arianna Krinos at akrinos (at) mit (dot) edu. 

### Acknowledgments

Authors: Arianna Krinos, Sarah Hu, Natalie Cohen, and Harriet Alexander.
# Community Guidelines

## How to contribute to `EUKulele`

If you are interested in modifying `EUKulele`, you may fork the project for your own use, as detailed in the [MIT License](https://github.com/AlexanderLabWHOI/EUKulele/blob/master/LICENSE) we have adopted for the project. In order to contribute, please contact the developers after making the desired changes, after which a pull request may be submitted. 

## Submitting an issue / Asking for help

If you have any suggestions for feature additions or any problems with or questions about how to use `EUKulele` that you would like addressed by and with the development community, please submit an issue on the [Issues tab](https://github.com/AlexanderLabWHOI/EUKulele/issues) of the project `GitHub` repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

## Collaborations

If you would like to seek out collaborations related to this project, you may contact Arianna Krinos at akrinos (at) mit (dot) edu. 

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
reported by contacting the project team at [akrinos@mit.edu](mailto:akrinos@mit.edu). All
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
---
title: 'EUKulele: Taxonomic annotation of the unsung eukaryotic microbes'
tags:
  - Python
  - Taxonomy
  - Metagenomics
  - Metatranscriptomics
authors:
  - name: Arianna I. Krinos
    orcid: 0000-0001-9767-8392
    affiliation: "1,2"
  - name: Sarah K. Hu 
    affiliation: "3,4"
    orcid: 0000-0002-4439-1360
  - name: Natalie R. Cohen
    affiliation: "3"
    orcid: 0000-0001-6156-9186
  - name: Harriet Alexander^[Corresponding author]
    orcid: 0000-0003-1308-8008
    affiliation: "1" # (Multiple affiliations must be quoted)
affiliations:
 - name: Biology Department, Woods Hole Oceanographic Institution, Woods Hole, MA, USA
   index: 1
 - name: MIT-WHOI Joint Program in Oceanography, Cambridge and Woods Hole, MA, USA
   index: 2
 - name: Marine Chemistry and Geochemistry, Woods Hole Oceanographic Institution, Woods Hole, MA, USA
   index: 3
 - name: Center for Dark Energy Biosphere Investigations, University of Southern California, Los Angeles, CA, USA
   index: 4
   
date: 13 December 2020
bibliography: paper.bib

---

# Summary

The assessment of microbial species biodiversity is essential in ecology and evolutionary biology [@reaka1996biodiversity], but especially challenging for communities of microorganisms found in the environment [@das2006marine;@hillebrand2018climate]. Beyond providing a census of organisms in the ocean, assessing marine microbial biodiversity can reveal how microbes respond to environmental change [@salazar2017marine], clarify the ecological roles of community members [@hehemann2016adaptive], and lead to biotechnology discoveries [@das2006marine]. Computational approaches to characterize taxonomic diversity and phylogeny based on the quality of available data for environmental sequence datasets is fundamental for advancing our understanding of the role of these organisms in the environment. Even more pressing is the need for comprehensive and consistent methods to assign taxonomy to environmentally-relevant microbial eukaryotes. Here, we present `EUKulele`, an open-source software tool designed to assign taxonomy to microeukaryotes detected in meta-omic samples, and complement analysis approaches in other domains by accommodating assembly output and providing concrete metrics reporting the taxonomic completeness of each sample.

`EUKulele` is motivated by ongoing efforts in our community to create and curate databases of genetic and genomic information [@phylodb;@caron2017probing;@marineref;@richter2020eukprot]. For decades, it has been recognized that genetic and genomic techniques are key to understanding microbial diversity [@fell1992partial]. Genetic approaches are particularly useful in poorly-understood or difficult-to-access environmental systems, which may have a high degree of species diversity [@das2006marine;@mock2016bridging]. The most common approach for censusing microbial diversity is genetic barcoding, which targets the hyper-variable regions of highly conserved genes such as 16S or 18S rRNA  [@leray2016censusing]. Computational approaches to assess the origin of these barcode-based studies (or tag-sequencing) have been well established [@Bolyen2018;@schloss2009introducing], and enable biologists to compare microbial communities and estimate sequence phylogeny. The recent collation of reference databases, e.g. PR2 and EukRef, for ribosomal RNA in eukaryotes have enabled more accurate taxonomic assessment [@del2018eukref;@guillou2012protist]. However, barcoding approaches that focus on single marker genes or variable regions limit the field of view of microbes--especially protists, which have complex and highly variable genomes [@del2014others]--potentially limiting the organisms recovered and leaving the “true” diversity poorly constrained [@piganeau2011and;@caron2019we]. 

Shotgun sequencing approaches (e.g., metagenomics and metatranscriptomics) have become increasingly tractable, emerging as a  viable, untargeted means to simultaneously assess community diversity and function. Large-scale meta-omic surveys, such as the Tara Oceans project [@zhang2015tara], have presented opportunities to assemble and annotate full "genomes" from environmental metagenomic samples [@tully2018reconstruction;@Delmont2018] and assemble massive eukaryotic gene catalogs from environmental metatranscriptomic samples [@Carradec2018]. The interpretation of these meta-omic surveys hinges upon curated, culture-based reference material. Several curated databases that contain predicted proteins from a mixture of genomic and transcriptomic references from eukaryotes, as well as bacteria and archaea have been created [e.g., @phylodb;@eukzoo;@marineref;@richter2020eukprot]. Building upon the creation of high-quality reference databases, we sought to create a tool similar to MEGAN [@beier2017functional], CCMetagen [@marcelino2020ccmetagen], and MG-RAST [@keegan2016mg], but independent of NCBI databases and useful for both metagenomes and metatranscriptomes, as well as the study of environmental eukaryotes. Further, we sought to create a tool with a single function to download and format databases, which is necessary for computational tools to remain relevant and usable as reference databases grow.

![A flowchart describing the general workflow of the software as it relates to metatranscriptomes (METs) and metagenomes (MAGs). \label{fig:eukulelediagram}](eukulele_diagram_simplified.png){ height=50% }

## Implementation

We built a tool with default databases MMETSP [@caron2017probing], PhyloDB [@phylodb], EUKZoo [@eukzoo], MarineRefII [@marineref], and EukProt [@richter2020eukprot], for optimum compatibility with environmental eukaryotic sequences. In particular, the Marine Microbial Eukaryote Transcriptome Sequencing Project (MMETSP) database, which contains over 650 fully-assembled reference transcriptomes [@keeling2014marine;@johnson2019re], is among the largest single projects to create a unified reference. These databases are an invaluable resource yet, to our knowledge, no single integrated software tool currently exists to enable an end-user to harness databases in a consistent and reproducible manner. 

``EUKulele`` [@eukulele] (Figure \ref{fig:eukulelediagram}) is an open-source ``Python``-based package designed to simplify taxonomic identification of marine eukaryotes in meta-omic samples. The package is written in `Python`, but may be installed as a `Python` module via [PyPI](https://pypi.org/), as a standalone tool via `conda`, or through download of the `EUKulele` tarball through `GitHub`. User-provided metatranscriptomic or metagenomic samples are aligned against a database of the user's choosing, using a user-chosen aligner (``BLAST`` [@kent2002blat] or ``DIAMOND`` [@buchfink2015fast]). The "blastx" utility is used by default if metatranscriptomic samples are only provided in nucleotide format, while the "blastp" utility is used for samples available as translated protein sequences. Any consistently-formatted database may be used, but five microbial eukaryotic database options are provided by default: MMETSP [@keeling2014marine;@caron2017probing;@johnson2019re], PhyloDB [@phylodb], EukProt [@richter2020eukprot], EukZoo [@eukzoo], and a combination of the MMETSP and MarRef [@keeling2014marine;@caron2017probing;@johnson2019re;@klemetsen2018mar] (referred to as MarRef-MMETSP). This final database is the default database option, and allows the eukaryotic sequences to be compared against the expansive and high-quality MMETSP, while also distinguishing prokaryotic sequences that may be present in the sample. The package returns comma-separated files containing all of the contig matches from the metatranscriptome or metagenome, as well as the total number of transcripts that matched, at each taxonomic level, from domain or supergroup to species. If a quantification tool has been used to estimate the number of counts associated with each transcript ID, counts may also be returned. 

After a desired database is either specified by the user from a previous install, or downloaded by the program (with the MarRef-MMETSP database downloaded by default), the user-selected alignment tool will create a database from the reads in the database peptide file. That database is aligned against the sample metatranscriptomic or metagenomic reads, resulting in a transcript ID, a percentage identity, e-value, and bitscore, all of which are common metrics in bioinformatics for assessing the quality of an alignment comparison between sequences. The user can specify which metric should be used for filtering out low-quality matches.

The alignment output is compared to an accompanying phylogenetic reference specific to the database (which can be generated via a script included in the package). Taxonomy is estimated at eight levels of taxonomic resolution, labeled as they are defined in the MMETSP [@keeling2014marine] and MarRef [@klemetsen2018mar] from “species” to “domain"/"supergroup”. Additionally, the software returns barplots displaying the relative composition of each sample at each taxonomic level, according to the number of transcripts or number of estimated counts if provided from `Salmon` (an external transcript quantification tool [@patro2017salmon]), which enable users to get a quick sense of the diversity in their metagenomic or metatranscriptomic sample. For metagenomic samples, a consensus taxonomic annotation is assigned based on the majority assignment of the contigs in the metagenome-assembled genome (MAG). For the metatranscriptomic option, only the taxonomic breakdown of the mixed community detected in the assembly will be returned. 

`EUKulele` will assess the relative "completeness" of a given taxonomic group by taking a user-inputted list of names at some taxonomic level to determine BUSCO completeness and redundancy [@simao2015busco]. For example, if the user was interested whether there was a set of relatively complete contigs available for genus *Phaeocystis* within their metagenomic sample, they could pass *Phaeocystis*, along with its taxonomic level, "genus", to `EUKulele`. By default, `EUKulele` will assess the BUSCO completeness of the most commonly encountered classifications at each taxonomic level.  To do this, `BUSCO` [@simao2015busco] is used to identify the core eukaryotic genes present in each sample. Using the list of genes identified as "core", a secondary taxonomic estimation step (and consensus assignment step, for MAGs) is performed to compare the taxonomic assignment predicted using all of the genes in comparison to the assignment made using only the genes that would be expected to be found in most reference transcriptomes. This approach is particularly useful for MAGs, and offers a method for avoiding conflicting or spurious matches made due to strain-level inconsistencies. For metatranscriptome samples, BUSCO completeness can be used to estimate the completeness of taxonomic groups to better inform their downstream interpretation. 

## Statement of Need

A growing number of databases have been created to catalog eukaryotic and bacterial diversity, but even when the same database is used, taxonomic assessment is not always consistent and fully documented [@rasheed2012metagenomic;@menzel2016fast]. Databases often contain distinct compilations of organisms and custom databases are commonly compiled for only a particular study [@kranzler2019silicon;@geisen2015metatranscriptomic;@obiol2020metagenomic]. Database variability might influence interpretation by splitting taxonomic annotations between groups, and often impacts the proportion of contigs that are annotated [@price2017robust]. A software tool can bridge the gap between database availability and efficient taxonomic assessment, making environmental meta-omic analyses more reproducible. Further, such a tool can control and assess the quality of the annotation, enable inference for specific organisms or taxonomic groups, and provide more conservative annotation in the case of organisms with exceptional amounts of inter-strain variability. We have designed the `EUKulele` [@eukulele] package to enable efficient and consistent taxonomic annotation of metagenomes and metatranscriptomes, in particular for eukaryote-dominated samples.

### Future Outlook

As single species isolates continue to be sequenced, databases are growing and becoming more reliable for assigning taxonomy in diverse environmental communities. `EUKulele` provides a platform to enable the repeated and consistent linkage of these databases to metagenomic and metatranscriptomic analyses. Taxonomic annotation is not the only desired outcome of meta-omic datasets against available databases, hence we envision eventually integrating functional annotation into the `EUKulele` package.

# Acknowledgements

This software was developed with support from the Computational Science Graduate Fellowship (DOE; DE-SC0020347) awarded to AIK and from the Woods Hole Oceanographic Independent Research & Development grant awarded to HA. NRC was supported by grant 544236 from the Simons Foundation. The Center for Dark Energy Biosphere Investigations (C-DEBI; OCE-0939564) supported the participation of SKH through a C-DEBI Postdoctoral Fellowship. The High-Performance Computing cluster at Woods Hole Oceanographic Institution (Poseidon) was used to generate assemblies and run `EUKulele`.


# Author Contributions

HA and SKH conceived the original idea for the tool. HA wrote the initial code for the tool. HA and AIK refined ideas for `EUKulele` related to designing it as an installable package and adding the feature of classification via core gene taxonomy for metagenomic applications. AIK developed the Python package code, wrote the conda package, implemented multiple alignment tools, and the `BUSCO` integration. AIK, NC, SKH, and HA wrote tests and documentation. HA and AIK wrote the paper.

# References
---
title: EUKulele
summary: A quick introduction to the documentation
---

# Welcome to EUKulele's documentation!

## Table of contents
1. [Source Code](#source)
2. [Badging](#badging)
3. 

## Source Code <a name="source"></a>

You can find the [source](http://github.com/AlexanderLabWHOI/EUKulele) for EUKulele on GitHub.

## Badging <a name="badging"></a>

[![Coverage Status](https://coveralls.io/repos/github/AlexanderLabWHOI/EUKulele/badge.svg?branch=remodeling)](https://coveralls.io/github/AlexanderLabWHOI/EUKulele?branch=master)

[![Build Status](https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master)](https://travis-ci.com/AlexanderLabWHOI/EUKulele)


There are a few ways to install `EUKulele`, but the easiest is:

```
pip install EUKulele
```

There may some dependencies that are not immediately downloaded when `EUKulele` is pip-installed. If this is the case, the offending dependencies may be installed individually... EUKulele documentation master file, created by
   sphinx-quickstart on Fri Aug  7 07:51:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EUKulele's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   quickstart
   about
   running-eukulele
   install
   databaseandconfig
   outputstructure
   tutorial
   community
   
.. toctree::
   :maxdepth: 2
   :caption: Documentation
   
   api
   parameters
   auxiliaryscripts
   
.. toctree::
   :maxdepth: 2
   :caption: Example Analyses  
   
   EUKulele_read_breakdown
   available-databases-R
   
Source Code
-----------

You can find the `source <http://github.com/AlexanderLabWHOI/EUKulele>`_ for EUKulele on GitHub.

Badging
-------

.. image:: https://coveralls.io/repos/github/AlexanderLabWHOI/EUKulele/badge.svg?branch=remodeling
   :target: https://coveralls.io/github/AlexanderLabWHOI/EUKulele?branch=remodeling

.. image:: https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master
   :target: https://travis-ci.com/AlexanderLabWHOI/EUKulele
   
Flowchart
---------

.. image:: eukulele_flowchart.png
  :width: 400
  :alt: Flowchart of EUKulele overview
  
Caveats
-------

See :ref:`Using EUKulele<usingeukulele>` for caveats about the use of predicted proteins, rather than nucleotide contigs, for metagenomic sequences. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _GitHub: http://github.com/AlexanderLabWHOI/EUKulele
.. _Git Issues: http://github.com/AlexanderLabWHOI/EUKulele/issues
Auxiliary Scripts
=================

A few scripts are included with ``EUKulele`` that can be used outside of the software, or before invoking the program to obtain the required configuration.

Create Protein Table
--------------------

``create_protein_table.py`` is used to generate a taxonomy table and protein map JSON file from a provided taxonomy table tab-delimited file, and a FASTA file containing database protein sequences. The file is invoked using the following arguments:

- ``--infile_peptide``: the FASTA file containing the database peptide sequences
- ``--infile_taxonomy``: the tab-separated file containing the taxonomy of each peptide sequence in the database
- ``--outfile_json``: the protein map JSON file to be written
- ``--output``: the formatted taxonomy table to be written
- ``--delim``: the delimiter that separates tokens in the FASTA headers in the peptide sequence file
- ``--col_source_id``: the column in your tab-separated taxonomy file containing the name of the strain
- ``--taxonomy_col_id``: the column containing the strain taxonomy in the tab-separated taxonomy file
- ``--column``: either a character name for the token in the FASTA headers of the peptide sequence file containing the strain name (matching what is in the taxonomy file), or a numeric value indicating the order of the token
- ``--reformat_tax``: if included, indicates that the taxonomy table should be reformatted instead of left as-is
- ``--euk-prot``: if included, means we are using input from the EukProt database, which has different formatting features.

Download Database
-----------------

Used within ``EUKulele``, ``download_database.sh`` may also be used independently of ``EUKulele`` to download one of the available databases provided with ``EUKulele`` (see Section :ref:`databases`). Invoked via::

    download_database.sh <DATABASE> <REF_FASTA> <REF_TABLE> <REF_FASTA_URL> <REF_TABLE_URL> <REFERENCE_DIR>
    
Where <DATABASE> is the name of the database, <REF_FASTA> is the name of the reference FASTA file to be generated/downloaded, REF_TABLE is the same but for the tab-delimited taxonomy table, <REF_FASTA_URL> and <REF_TABLE_URL> are the URLs to download said databases from (provided in ``reference_url.yaml``), and <REFERENCE_DIR> is where to store the resulting database.
.. _parameters:
====================================
Parameters
====================================

A full list of parameters can be found in the table at the bottom of this page. However, in practice, only a few parameters will be relevant for most users of :code:`EUKulele`. These are the required ones:

- mets_or_mags: Whether the user intends to run the analysis for metatranscriptomic samples ("mets") or metagenomic samples ("mags")


.. list-table:: Full list of ``EUKulele`` parameters
   :widths: 25 25 50
   :header-rows: 1

   * - Flag
     - Configuration File Entry
     - Meaning
   * - ``--config``
     - N/A 
     - The path to the configuration file which should be used to retrieve the equivalent of command-line arguments.
   * - ``-m/--mets_or_mags`` 
     - ``mets_or_mags`` 
     - A required flag to indicate whether metatranscriptomic ("mets") or metagenomic ("mags") samples are being used as input.
   * - ``-s/--sample_dir`` 
     - samples 
     - A required flag to indicate where the samples (metagenomic or metatranscriptomic, depending on "mets_or_mags" flag) are located. 
   * - ``-o/--out_dir`` 
     - output 
     - The path to the directory where output will be stored. Defaults to a folder called ``output`` in the present working directory.
   * - ``--reference_dir`` 
     - reference 
     - A flag to indicate where the reference FASTA is stored, or a keyword argument for the dataset to be downloaded and used. Only used if not downloading automatically.
   * - ``--ref_fasta`` 
     - ref_fasta 
     - The name of the reference FASTA file in ``reference_dir``; defaults to reference.pep.fa if not specified, or is set according to the downloaded file if using a keyword argument.
   * - ``--database`` 
     - database 
     - An optional additional argument for specifying the database name. If the database specified is one of the supported databases (currently, "mmetsp", "eukprot", or "phylodb", it will be downloaded automatically. Otherwise, MMETSP is used as a default. 
   * - ``--run_transdecoder``
     - run_transdecoder (set to 0 or 1)
     - An argument for the user to specify whether or not TransDecoder should be used to translate input nucleotide sequences, prior to ``blastp`` being used (i.e., the equivalent protein-protein alignment with the tool of choice). If included in command line or set to 1 in configuration file, ``TransDecoder`` is run. Otherwise, ``blastp`` is run if protein files are found (according to files in the sample directory ending in ``--p_ext`` (below), or ``blastx`` is run if only nucleotide format files are found. 
   * - ``--nucleotide_extension/--n_ext`` 
     - nucleotide_extension 
     - The file extension for samples in nucleotide format (metatranscriptomes). Defaults to .fasta.
   * - ``--protein_extension/--p_ext`` 
     - protein_extension 
     - The file extension for samples in protein format (metatranscriptomes). Defaults to .faa.
   * - ``-f/--force_rerun`` 
     - force_rerun 
     - If included in a command line argument or set to 1 in a configuration file, this argument forces all steps to be re-run, regardless of whether output is already present.
   * - ``--use_salmon_counts`` 
     - use_salmon_counts 
     - If included in a command line argument or set to 1 in a configuration file, this argument causes classifications to be made based both on number of classified transcripts and by counts.
   * - ``--salmon_dir`` 
     - salmon_dir 
     - If ``--use_salmon_counts`` is true, this must be specified, which is the directory location of the ``salmon`` output/quantification files.
   * - ``--names_to_reads`` 
     - names_to_reads 
     - A file that creates a correspondence between each transcript name and the number of ``salmon``-quantified reads. Can be generated manually via the ``names_to_reads.py`` script, or will be generated automatically if it does not exist. \
   * - ``--transdecoder_orfsize`` 
     - transdecoder_orfsize 
     - The minimum cutoff size for an open reading frame (ORF) detected by ``TransDecoder``. Only relevant if ``--use_transdecoder`` is specified.
   * - ``--alignment_choice`` 
     - alignment_choice 
     - A choice of aligner to use, currently ``BLAST`` or ``DIAMOND``.
   * - ``--cutoff_file`` 
     - cutoff_file 
     - A ``YAML`` file, provided in ``src/EUKulele/static/``, that contains the percent identity cutoffs for various taxonomic classifications. Any path may be provided here to a user-specified file.
   * - ``--filter_metric`` 
     - filter_metric 
     - Either evalue, pid, or bitscore (default evalue) - the metric to be used to filter hits based on their quality prior to taxonomic estimation. 
   * - ``--consensus_cutoff`` 
     - consensus_cutoff 
     - The value to be used to decide whether enough of the taxonomic matches are identical to overlook a discrepancy in classification based on hits associated with a contig. Defaults to 0.75 (75%). 
   * - ``--busco_file`` 
     - busco_file 
     - Overrides specific organism and taxonomy parameters (next two entries below) in favor of a tab-separated file containing each organism/group of interest and the taxonomic level of the query. \
   * - ``--organisms``
     - organisms
     - A list of organisms/groups to test the BUSCO completeness of matching contigs for.
   * - ``--taxonomy_organisms`` 
     - taxonomy_organisms 
     - The taxonomic level of the groupings indicated in the list of ``--organisms``; also a list.
   * - ``--individual_or_summary / -i``
     - individual_or_summary 
     - Defaults to summary. Whether BUSCO assessment should just be performed for the top organism matches, or whether the list of organisms + their taxonomies or BUSCO file (above parameters) should be used (individual). When ``-i`` is specified, individual mode is chosen.
   * - ``--busco_threshold``
     - busco_threshold 
     - The threshold for BUSCO completeness for a set of contigs to be considered reasonably BUSCO-complete.
   * - ``--tax_table`` 
     - tax_table 
     - The name of the formatted taxonomy table; defaults to "tax-table.txt.". If this file is not found, it can be generated from the reference FASTA and original taxonomy file using the provided script ``create_protein_file.py``, or the database specified will be automatically downloaded, if it is one of the supported databases.
   * - ``--protein_map``
     - protein_map 
     - The name of the JSON file containing protein correspondences; defaults to "protein-map.json". If this file is not found, it can be generated from the reference FASTA and original taxonomy file using the provided script ``create_protein_file.py``, or the database specified will be automatically downloaded, if it is one of the supported databases.
     ``EUKulele`` Quick Start
====================================

We recommend consulting the full documentation to explore all of the ``EUKulele`` features, capabilities, and intended usage. A quick start approach for annotating eukaryotic metagenome-assembled genomes and metatranscriptomes are outlined below using ``EUKulele``-provided databases and default parameters. 

Installation
------------

The conda installation is the most straightforward approach and the environment will contain all dependencies needed to run ``EUKulele``. If you do not already have conda installed on your machine, see the conda installation documentation `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_. Other installation options are included under Installation and Invocation::

    conda create -n EUKulele
    conda activate EUKulele
    conda install -c akrinos -c bioconda -c conda-forge EUKulele

Generalized flow of EUKulele annotation
---------------------------------------

.. image:: eukulele_diagram_simplified.jpg
  :width: 400
  :alt: Flowchart of EUKulele overview
  
 
Metagenome-assembled genome (MAG) annotation
--------------------------------------------

``EUKulele`` can determine the taxonomic identity of binned MAGs using their consensus contig annotations. It run using the following command::

    EUKulele --sample_dir/output_directory -m mags

where ``output_directory`` contains one or more assembly fasta files with the extension ``.faa`` (it is recommended that MAG files are provided in protein format; see :ref:`Using EUKulele<usingeukulele>`). See :ref:`Parameters<parameters>` for other file extension accommodations. This will annotate the assemblies using the MMETSP database (default) and DIAMOND aligner (default). 

Useful output is located in:

``output/max-level-mag`` - the majority level annotation by level (supergroup - species) for each mag and the proportion of proteins that share that same annotation

``output/levels_mags`` - consensus annotations based on the majority of annotated contigs at each classification level

``output/core_taxonomy_estimation`` - contains contig annotations for core (shared) genes among MAGs as determined by BUSCO. Useful for identifying core genes among strains/variants during pangenomic analyses

``output/busco_assessment`` -  BUSCO assessment for the most abundant taxonomic annotation at each classification level. Provides information on estimated MAG completion based on conserved eukaryotic genes expected to be present in a full genome.

Metatranscriptomic (MET) annotation
-----------------------------------

``EUKulele`` can be run on metatranscriptomic assemblies using the following command:

``EUKulele --sample_dir /output_directory -m mets``

where ``output_directory`` contains one or more assembly fasta files with the extension ``.fasta``. See Parameters for other file extension accommodations. This will annotate the input assemblies using the MMETSP database (default) and DIAMOND aligner (default).

Useful output files will be located in:

``output/taxonomy_estimation`` - contains each contig’s annotation and the least common ancestor (LCA) classification threshold determined

``output/taxonomy_visualization`` - shows the breakdown of contig annotations at the seven classification levels (Supergroup, Division, Class, Order, Family, Genus, Species) for the mixed community. 


If you have any issues using ``EUKulele`` or suggestions to improve the software, please submit an `issue on GitHub <https://github.com/AlexanderLabWHOI/EUKulele/issues>`_.
Expected Output of ``EUKulele``
================================

Below is what you should expect to see when you run ``EUKulele``. ``output-folder-name`` is the top level directory specified by the ``-o/--out_dir`` parameter, which defaults to a folder called "output" in the current working directory.

| output-folder-name
| ├── busco
| │   ├── build
| │   ├── make.bat
| │   ├── Makefile
| │   └── source
| ├── busco_assessment
| ├── core_taxonomy_counts
| ├── core_taxonomy_estimation
| ├── core_taxonomy_visualization
| ├── mets (if running ``TransDecoder``, or protein files provided)
| │   ├── *.faa files*
| │   └── transdecoder
| │   |   └── *other TransDecoder outputs*
| ├── mets_full / mags_full (depending on ``--mets_or_mags`` parameter)
| │   └── diamond (or blast)
| ├── mets_core / mags_core (depending on ``--mets_or_mags`` parameter)
| │   └── diamond (or blast)
| ├── taxonomy_counts
| ├── taxonomy_estimation
| └── taxonomy_visualization
|


Taxonomy Estimation Folders
---------------------------

Inside each of the taxonomy estimation folders (``core_taxonomy_estimation``, for exclusively transcripts annotated as core genes, and ``taxonomy_estimation``), there are files labeled ``<sample_name>-estimated-taxonomy.out``. Each of these files has the following columns:

- transcript_name
    - The name of the matched transcript/contig from this sample file
- classification_level
    - The most specific taxonomic level that the transcript/contig was classified at
- full_classification
    - The full taxonomic classification for the transcript/contig, up to and including the most specific classification
- classification
    - Just the most specific classification obtained, at the taxonomic level specified by classification_level
- max_pid
    - The maximum percentage identity reported by the alignment program for the match
- counts
    - 1 if not using ``Salmon``, or the number of counts reported by the quantification program for that transcript/contig if using ``Salmon``
- ambiguous
    - 1 if there were multiple disagreeing matches for this transcript/contig, resolved by either consensus annotation using the user-provided cutoff (defaults to 75%) or by last common ancestor (LCA) approaches, otherwise 0

Taxonomy Counts Folders
-----------------------

Inside each of the taxonomy counts folders (``core_taxonomy_counts`` and ``taxonomy_counts``), there is a file with the aggregated counts belonging to each annotation at the possible taxonomic levels. These files are derived diectly from the taxonomy estimation and used for the visualization steps. Inside this folder are comma-separated files with the following headings:

- *Taxonomic Level*
    - The taxonomic level specified by the file name
- GroupedTranscripts
    - A semicolon-separated list of the identification names for transcripts that matched to this taxonomic label
- NumTranscripts
    - The total number of transcripts (the length of the semicolon-separated list of transcripts)
- Counts
    - Provided if ``Salmon`` counts are used - the matching summed number of counts for each label
- Sample
    - The original metagenomic/metatranscriptomic sample that this count is from (a separate row would be provided if the match is found in multiple samples)
    
The taxonomic count files are named according to the convention ``<output-folder-name>_all_<taxonomic-level>_counts.csv``.

Taxonomy Visualization Folders
------------------------------

Inside each of the taxonomy visualization folders (``core_taxonomy_visualization``, for exclusively transcripts annotated as core genes, and ``taxonomy_visualization``), there are auto-generated barplots that show:

- x-axis: samples
- y-axis, left subplot (if using counts): relative number of transcripts
- bars, left subplot (if using counts): each of the top represented taxonomic groups (must represent >= 5% of total number of transcripts)
- y-axis, right subplot (if using counts): relative number of counts
- bars, right subplot (if using counts): each of the top represented taxonomic groups (must represent >= 5% of total counts)

The right subplot is only generated if counts from a quantification tool (namely, ``Salmon``) are provided.
.. _databases:
   =====================================================
Installing Databases and Creating Configuration Files
=====================================================

Default Databases
-----------------

Four databases can be downloaded and formatted automatically when invoking ``EUKulele``. Currently the supported databases are:

- `PhyloDB <https://drive.google.com/drive/u/0/folders/0B-BsLZUMHrDQfldGeDRIUHNZMEREY0g3ekpEZFhrTDlQSjQtbm5heC1QX2V6TUxBeFlOejQ>`_
- `EukProt <https://figshare.com/articles/EukProt_a_database_of_genome-scale_predicted_proteins_across_the_diversity_of_eukaryotic_life/12417881/2>`_
- `MMETSP <https://zenodo.org/record/1212585#.Xw3PoJNKhTZ>`_ 
- `MMETSP <https://zenodo.org/record/1212585#.Xw3PoJNKhTZ>`_ and `MMETSP <https://mmp.sfb.uit.no/databases/marref/#/>`_ *Default*
- `EukZoo <https://github.com/zxl124/EukZoo-database>`_

Note that the MMETSP database is generated using cleaned MMETSP assemblies originally derived from, but not identical to, the assemblies stored in full at the link above. In order to download the cleaned assemblies used to create the ``EUKulele`` MMETSP database, please follow the instructions recorded in `this Github repository <https://github.com/shu251/download-cleaned-mmetsp>`_.

To use these databases, all you need to do is specify ``--database phylodb``, ``--database eukprot``, or ``--database mmetsp``, respectively, when invoking ``EUKulele``. 

A database (for example ``phylodb``) can be setup prior to running by using::

    EUKulele download --database phylodb

If a database is not found automatically by ``EUKuele`` it will automatically download the database specified by the flag. If you downloaded a database previously you can specify the ``--reference_dir`` flag indicating the path to the previously downloaded database. If no reference database is specified with ``--reference_dir``, EUKulele will automatically download and use the MMETSP database. You can also (1) download the other databases and use the flag ``reference_dir`` to point EUKulele to the location of already downloaded databases or (2) use your own databases.

Composition of Default Databases
--------------------------------

Several databases are automatically provided and formatted for ``EUKulele`` in order to let the user make the best decision. Often database choice depends on the relevant research question. Therefore, we also provide information on the contents of each database below.

Note that databases provided through ``EUKulele`` include the taxonomic structure used for the creation of the database. Therefore, the structure of the taxonomic assignment varies database to database. For instance, PhyloDB includes other microbial domains classified under "Supergroup", while the other databases only include eukaryotic references. Users should take note when they compiled downstream results from different databases. 

.. list-table:: Broad overview of each database.
      :widths: 25 12 25 13 25
   :header-rows: 1

   * - Database
     - Total entries
     - Domains
     - Unique species (strains) [*Note: not all databases have both species and strain level*]
     - Taxonomy Assignment Structure
   * - MMETSP
     - 678
     - Eukaryote
     - 316 (405)
     - Supergroup, Division, Class, Order, Family, Genus, Species, Strain
   * - PhyloDB
     - 25,992
     - Eukaryote, bacteria, archaea, virus
     - 25,992
     - Supergroup, Division, Class, Order, Family, Genus_Species
   * - EukZoo
     - 739
     - Eukaryote
     - 361 (441)
     - Supergroup, Phylum, Class, Order, Family, Genus, Species
   * - EukProt
     - 742
     - Eukaryote
     - (621)
     - Supergroup_UniEuk, Taxogroup_UniEuk, Epithet_UniEuk, Genus_UniEuk, Strain
     
.. image:: mmetsp-doughnut.png
     :width: 400
  :alt: Breakdown of MMETSP database composition.
  
.. image:: phylodb-doughnut.png
     :width: 400
  :alt: Breakdown of PhyloDB database composition.
  
.. image:: eukprot-doughnut.png
     :width: 400
  :alt: Breakdown of EukProt database composition.
  
.. image:: eukzoo-doughnut.png
     :width: 400
  :alt: Breakdown of EukZoo database composition.
  
Recommendations for Database Usage
----------------------------------

``EUKulele`` was initially designed for use with the MMETSP database. As this is the most complete resource for reference transcriptomes of eukaryotic species, it is the recommended database. 

A highly recommended approach is to run EUKulele with both MMETSP and PhyloDB. Since PhyloDB includes non-eukaryotic domains, a best hit eukaryotic reference from both MMETSP and PhyloDB will provide higher confidence. 


Using Other Databases
---------------------

The basic requirements for using a database with ``EUKulele`` are:

- A singular protein FASTA file containing the sequences
- A taxonomy table file which contains, for each transcriptome sample in the protein FASTA database:
    - Source_ID: what identifier/organism the transcript came from, which typically should be specified in the header within the FASTA file
    - Supergroup 
    - Division
    - Class
    - Order
    - Family
    - Genus
    - Species
- A JSON file containing a list of dictionary correspondences between each Source ID and transcript ID 
    - If you have a separate correspondence between transcript IDs and the organism each transcript ID came from, this prevents you from having to have the Source ID in the transcript header
    - Example: ``{"CAMPEP_0174983734": "MMETSP0004", "CAMPEP_0174982176": "MMETSP0004", "CAMPEP_0184404416": "MMETSP0007"}`` for a database of three transcripts coming from two different Source IDs (``MMETSP0004`` and ``MMETSP0007``)
    
These taxonomy table and JSON file can be generated using the ``create_protein_file`` script provided with ``EUKulele``. This script is invoked via::

    create-protein-table.py --infile_peptide <peptide fasta file> --infile_taxonomy <taxonomy file> --outfile_json <name of protein map JSON file> --output <name of taxonomy file> [--delim <delimiter> --column <column>] 
    
when ``EUKulele`` is installed. 

- ``--infile_peptide``
    - The peptide FASTA file for the database
- ``--infile_taxonomy``
    - The original taxonomy file
- ``--col_source_id``
    - Optional; defaults to "Source_ID"; the column in the taxonomy file that corresponds to the Source ID in the database
- ``--reformat_tax``
    - If this tag is included, the taxonomy will be split according to the contents of the column labeled with the ``taxonomy_col_id`` that is specified by the tag below (instead of 7 different columns corresponding to each taxonomic level as in the listing above)
- ``--taxonomy_col_id``
    - Only relevant if ``--reformat_tax`` is specified. The column (e.g. "taxonomy" as in the default) that contains a semicolon-separated list of the taxonomic levels to be separated into columns
- ``--outfile_json``
    - The name of the output protein map file to be created. To use the output most easily with ``EUKulele``, this file should be called ``prot-map.json`` (as is the default) and placed in the same )nce directory with the reference protein FASTA file, which ideally would be named ``reference.pep.fa`` to facilitate working with the defaults. Then, just specify this output folder as ``--reference_dir`` when invoking ``EUKulele``
- ``--output``
    - The name of the output taxonomy table file to be created. To use the output most easily with ``EUKulele``, this file should be called ``tax-table.txt`` (as is the default) and placed in the same reference directory with the reference protein FASTA file, which ideally would be named ``reference.pep.fa`` to facilitate working with the defaults. Then, just specify this output folder as ``--reference_dir`` when invoking ``EUKulele``
- ``--delim``
    - What to split the FASTA headers on in the protein database file, typically ``\t``
- ``--column``
    - The label to be used for the Source_ID parsed from the reference peptide FASTA headers. This is such that the protein map JSON file can be created from the transcript IDs. So if your transcripts include a tab-separated list of entries that includes ``SOURCE_ID=XXXXX``, as in the MMETSP, include a string here for the label before the equals sign. If instead the Source ID occurs at a predictable position in the parsed FASTA headers, a number can be included for this parameter
- ``--euk-prot``
    - Should only be used if you are specifically creating a table and protein map for the EukProt database, which has a few particular features to take into account
    
Customizing the Taxonomic Identification Cutoffs
------------------------------------------------

By default, ``EUKulele`` uses the following percent identity cutoffs to determine taxonomic matches::

    species: 95
    genus: 80
    family: 65
    order: 50
    class: 30
    
To change these cutoffs, simply create a YAML file containing these entries exactly as written above, and provide this cutoff file as input to ``EUKulele`` via ``--cutoff_file <name of YAML file you created>``. A YAML file is a Markdown document that can be used to quickly parse and deliver new variables to a script. In order to produce this YAML file, you would create a text file containing exactly the text above (i.e., line 1 would be "species: 95", or whatever you desire as a cutoff, and so on), and then save the file with the YAML extension.

Community Guidelines 
====================================

How to contribute to :code:`EUKulele`
-------------------------------------

If you are interested in modifying :code:`EUKulele`, you may fork the project for your own use, as detailed in the `MIT License
<https://github.com/AlexanderLabWHOI/EUKulele/blob/master/LICENSE>`_ we have adopted for the project. In order to contribute, please contact the developers after making the desired changes, after which a pull request may be submitted. 

Submitting an issue / Asking for help
-------------------------------------

If you have any suggestions for feature additions or any problems with or questions about how to use :code:`EUKulele` that you would like addressed by and with the development community, please submit an issue on the `Issues tab
<https://github.com/AlexanderLabWHOI/EUKulele/issues>`_ of the project :code:`GitHub` repository. You may want to search the existing issues before submitting, to avoid asking a question or requesting a feature that has already been discussed.

Collaborations
---------------

If you would like to seek out collaborations related to this project, you may contact Arianna Krinos at akrinos (at) mit (dot) edu. 
Trying out ``EUKulele``
=======================

A minimal working example for ``EUKulele`` can be tested using the instructions on this page.

You can obtain sample metagenomic data with simulated MAGs from the following Dropbox link, using the following commands on a Linux system::

    wget https://www.dropbox.com/s/l4kvbpqftdad5ib/sample_eukulele.tar.gz?dl=1
    tar xzf sample_eukulele.tar.gz?dl=1

Then ``cd`` into ``sample_EUKulele``. We have noticed a problem with tarring across systems and through Dropbox. So you may be required to run:: 

    rm samples_MAGs/._sample_*

Now, you should have a clean folder containing only metagenomic samples that you can use for test driving ``EUKulele``!

First, let's create an environment to run ``EUKulele`` in.  Create a ``conda`` environment and activate it, using::

    conda env create -f EUKulele-env.yaml
    conda activate EUKulele

This must also be done inside the directory you installed via Dropbox, which contains a ``conda`` configuration file that will make sure that every dependency of ``EUKulele`` that cannot be installed via pip is available on your system.

You should now download ``EUKulele`` via ``pip`` using::

    python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps EUKulele

If any dependency is not satisfied, you can install it manually using ``pip install <requirement> --user``. If you install ``EUKulele`` via ``conda`` instead of via ``pip``, all of the dependencies are installed for you automatically.

Once everything is setup and all dependencies satisfied, from within the ``sample_EUKulele`` directory, run::
    
    EUKulele --config curr_config.yaml

    
Where ``curr_config.yaml`` is a configuration file in the directory you downloaded which contains all of the flags needed for a basic ``EUKulele`` run. Feel free to open this file if you're curious, and compare it to the full list of parameters available to customize ``EUKulele``. 

Then check the folder ``test_out_23July`` in the current directory.
.. EUKulele documentation master file, created by
   sphinx-quickstart on Fri Aug  7 07:51:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EUKulele's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents
   
   quickstart
   about
   running-eukulele
   install
   databaseandconfig
   outputstructure
   tutorial
   community
   
.. toctree::
   :maxdepth: 2
   :caption: Documentation
   
   api
   parameters
   auxiliaryscripts
   
.. toctree::
   :maxdepth: 2
   :caption: Example Analyses  
   
   EUKulele_read_breakdown
   EUKulele-tuneR
   
Source Code
-----------

You can find the `source <http://github.com/AlexanderLabWHOI/EUKulele>`_ for EUKulele on GitHub.

Badging
-------

.. image:: https://coveralls.io/repos/github/AlexanderLabWHOI/EUKulele/badge.svg?branch=remodeling
   :target: https://coveralls.io/github/AlexanderLabWHOI/EUKulele?branch=remodeling

.. image:: https://travis-ci.com/AlexanderLabWHOI/EUKulele.svg?branch=master
   :target: https://travis-ci.com/AlexanderLabWHOI/EUKulele
   
Flowchart
---------

.. image:: eukulele_diagram_simplified.png
  :width: 400
  :alt: Flowchart of EUKulele overview
 
Caveats
-------

See :ref:`Using EUKulele<usingeukulele>` for caveats about the use of predicted proteins, rather than nucleotide contigs, for metagenomic sequences. 


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _GitHub: http://github.com/AlexanderLabWHOI/EUKulele
.. _Git Issues: http://github.com/AlexanderLabWHOI/EUKulele/issues
====================================
About ``EUKulele``
====================================

``EUKulele`` is a python package that enables the rapid taxonomic annotation of metagenomic and metatranscriptomic data using a Last Common Ancestor (LCA) approach to align samples to a reference database. Assignments are chosen using either (a) the consensus annotation above a user-specified threshold, by default 75% identical predicted taxonomic annotations, or (b) the last conserved taxonomic level between the identified matches. ``EUKulele`` has been designed to be flexible and reproducible. The software is flexible in that it allows any sequences the user is interested in to be added as a reference, and also provides the user the freedom to choose their own levels of taxonomic specificity and confidence in alignment matches. The pipeline is reproducible because the specific parameters used to run ``EUKulele`` and the output can be saved, containing information about which database was run using which version, allowing users to share metadata and easily defend their taxonomic annotation results. In addition to estimating a taxonomic annotation for each putative transcript (“contig”), `EUKulele` can predict the taxonomic identity of eukaryotic Metagenome Assembled Genomes (MAGs) using the consensus of the contigs in each MAG. `EUKulele` also offers users a secondary taxonomic annotation approach limited only to sequenced reads annotated as core eukaryotic genes. 

A variety of curated protein databases are available to use with `EUKulele`, which are downloaded and formatted automatically. A custom local database can also be created by the user, using user-provided sequences and accompanying taxonomy. 

Functionality
====================================

``EUKulele`` :cite:`eukulele` is an open-source ``Python``-based package designed to simplify the process of taxonomic identification of marine eukaryotes in meta-omic samples. User-provided metatranscriptomic or metagenomic samples are aligned against a database of the user's choosing, with an aligner of the user's choice (``BLAST`` :cite:`kent2002blat` or ``DIAMOND`` :cite:`buchfink2015fast`). The "blastx" utility is used by default if metatranscriptomic samples are only provided in nucleotide format, while the "blastp" utility is used for metagenomic samples and metatranscriptomic samples available as translated protein sequences. Optionally, the user may indicate a preference to translate nucleotide input sequences using the ``TransDecoder`` software :cite:`haastransdecoder`, with the output provided to "blastp". Any consistently-formatted database may be used, but three published microbial eukaryotic database options are provided by default: MMETSP :cite:`keeling2014marine,caron2017probing`, PhyloDB :cite:`phylodb`, and EukProt :cite:`richter2020eukprot`. The package returns comma-separated files containing all of the contig matches from the metatranscriptome or metagenome, as well as the total number of transcripts that matched, at each taxonomic level, from supergroup to species. If a quantification tool has been used to estimate the number of counts associated with each transcript ID, counts may also be returned. Additionally, the software returns barplots displaying the relative composition of each sample at each taxonomic level, according to the number of transcripts or number of estimated counts if provided from ``Salmon`` (an external transcript quantification tool :cite:`patro2017salmon`).

``EUKulele`` will assess the relative 'completeness' of a given taxonomic group by taking a user-inputted list of names at some taxonomic level to determine BUSCO completeness and redundancy :cite:`simao2015busco`. For example, if the user was interested whether there was a set of relatively complete contigs available for genus *Phaeocystis* within their metagenomic sample, they could pass *Phaeocystis*, along with its taxonomic level, "genus", to ``EUKulele``. By default, ``EUKulele`` will assess the BUSCO completeness of the most commonly encountered classifications at each taxonomic level. 

Usage and Dependencies
====================================

The package is written in ``Python``, but may be installed as a ``Python`` module via `PyPI <https://pypi.org/>`_, as a standalone tool via ``conda``, or through download of the ``EUKulele`` tarball through ``GitHub``.

After a desired database is either specified by the user from a previous install, or downloaded by the program, the user-selected alignment tool will create a database from the reads in the database peptide file. That database is aligned against the sample metatranscriptomic or metagenomic reads, resulting in a transcript ID, a percentage identity, e-value, and bitscore, all of which are common metrics in bioinformatics for assessing the quality of an alignment comparison between sequences. The user can specify which metric should be used for filtering out low-quality matches.

The alignment output is compared to an accompanying phylogenetic reference specific to the database (which can be generated via a script included in the package). Taxonomy is estimated at six levels of taxonomic resolution, labeled as they are defined in the MMETSP :cite:`keeling2014marine` from “species” to “supergroup”, and preliminary visualizations are provided as part of the package output, which enable users to get a quick sense of the diversity in their metagenomic or metatranscriptomic sample. For metagenomic samples, a consensus taxonomic annotation is assigned based on the majority assignment of the contigs in the metagenome-assembled genome (MAG). For the metatranscriptomic option, only the taxonomic breakdown of the mixed community detected in the assembly will be returned. If counts from ``Salmon`` :cite:`patro2017salmon` are provided, ``EUKulele`` also provides and visualizes the counts associated with each taxonomic classification.

Subsequently, ``BUSCO`` :cite:`simao2015busco` is used to identify the core eukaryotic genes present in each sample. Using the list of genes identified as "core", a secondary taxonomic estimation step (and consensus assignment step, for MAGs) is performed to compare the taxonomic assignment predicted using all of the genes in comparison to the assignment made using only the genes that would be expected to be found in most reference transcriptomes. This approach is particularly useful for MAGs, and offers a method for avoiding conflicting or spurious matches made due to strain-level inconsistencies. For metatranscriptome samples, BUSCO completeness can be used to estimate the completeness of taxonomic groups to better inform their downstream interpretation. 

.. bibliography:: refs.bib
   :cited:
Installation
============

Installing with ``conda``
-------------------------

:code:`EUKulele` may also be downloaded as a :code:`conda` package, which will eventually become the easiest option, as `conda` automatically installs all dependencies for the user. The package can currently be downloaded via::

    conda install -c akrinos -c bioconda -c conda-forge EUKulele
    
Eventually, the :code:`-c akrinos` designation will be replaced, when ``EUKulele`` becomes a ``bioconda`` package, at which point it may be installed from that channel (and this documentation will be updated). 

To ensure you have the most recent version of ``EUKulele``, you can run ``conda update EUKulele``, or you can check the ``EUKulele`` landing page on Anaconda Cloud to check the most recent version number.

If you find that the ``conda`` install is running slowly, you may also, alternatively, follow the steps below to install with ``mamba``.

Installing with ``mamba``
------------------------

:code:`EUKulele` can be installed more rapidly remotely using ``mamba``, a parallelized fast installer that is obtainable via ``conda``. In order to install via ``mamba``, the following steps are necessary. 

- Create a ``conda`` environment
- Install ``mamba``, using::
    
    conda install mamba -c conda-forge

- Install the package as before, but using ``mamba`` as the installer, rather than ``conda``::

    mamba install -c akrinos -c bioconda -c conda-forge EUKulele

This should provide a considerable speedup as compared to ``conda``.

Installing with Pip
-------------------

There are few different ways to install :code:`EUKulele`, but the easiest is::

    pip install EUKulele
    
There may some Python dependencies that are not immediately downloaded when :code:`EUKulele` is pip-installed. If this is the case, the offending dependencies may be installed individually.

The external dependencies can be installed individual, or with conda using::
   
    conda create -n EUKulele
    conda activate EUKulele
    conda install -c bioconda -c conda-forge blast busco=4.0.6 diamond transdecoder
    
As the dependencies external to PyPI are:

- BLAST
- BUSCO
- Diamond
- TransDecoder (if using metatranscriptome samples and the ``--use_transdecoder`` flag)

Cloning the Development Code from GitHub
----------------------------------------

In addition, you can clone :code:`EUKulele` from GitHub (the current development version) using::

    git clone https://github.com/AlexanderLabWHOI/EUKulele
    
And then invoke :code:`EUKulele` either by executing the script :code:`bin/EUKulele` directly, or by installing :code:`EUKulele` as a local package, by calling the following from the :code:`EUKulele` directory::

    python3 -m pip install -e . --user
    python setup.py install  --user
    
Again, in this case, external dependencies may be installed via::
   
    conda create -n EUKulele
    conda activate EUKulele
    conda install -c bioconda -c conda-forge blast busco=4.0.6 diamond transdecoder
    
Or individually for each software.

Invocation
==========

If installed either with pip or ``conda``, ``EUKulele`` can be invoked via::

    EUKulele <arguments>
    
Where the minimal command would be::

    EUKulele --mets_or_mags <choice of data type> --sample_dir <where samples are located>
    
In which case ``EUKulele`` would be run with mostly parameter defaults and using the MMETSP database, by default.

Less common use cases
---------------------

``EUKulele`` may also be run as a module within Python. Include the phrase ``import EUKulele`` in the header of a Python file. Then, you may execute ``EUKulele`` using::

    EUKulele.eukulele(config=config_file)

in the case that you have a configuration file to specify (replace the ``config_file`` variable with this path), or with::

    EUKulele.eukulele(string_arguments=string_of_arguments)

where ``string_of_arguments`` is a string containing the non-default `EUKulele` options you wish to specify.
.. _usingeukulele:
====================================
Using EUKulele
====================================

``EUKulele`` has been designed to provide taxonomic annotation for two primary data types: 1) contigs derived from metatranscriptomes (METs) and 2) metagenome assembled genomes (MAGs). The focus of ``EUKulele`` is on the annotation of microbial eukaryotes; however, it is conceivable to use any database as the foundation of your analyses (see Databases).

Metatranscriptomes (METs)
=========================
In the first case, metatranscriptomes (shortened in ``EUKulele`` to ``mets``), are assumed to be contigs generated from shotgun-style sequencing and assembly of metatranscriptomic data (RNA) from a mixed community. These contigs can be provide to ``EUKulele`` as either nucleotide sequences (such as those output by `Trinity <https://github.com/trinityrnaseq/trinityrnaseq/wiki>`_) or predicted protein sequences from these contigs (such as those output by `Transdecoder <https://github.com/transdecoder>`_). 

The most basic running of ``EUKulele`` on metatranscriptome samples  would be::

    EUKulele --sample_dir path/to/metatranscriptome/samples -m mets

This command will automatically download (if not downloaded already) the MMETSP database and use ``DIAMOND`` to align sequences to the database. ``EUKulele`` will automatically inspect the directory for files ending with ``.fna`` (nucleotide files) or ``.faa`` (amino acid files). If no files are found within the directory with those extensions an error will be raised. If nucleotide files are found (``.fna``) a blastx-style search will automatically performed; alternatively, if amino acid files are found a blastp-style search will be automatically performed.  The results of the ``DIAMOND`` searches and estimated taxonomy for each of the contigs within the input datasets as well as various figures will be output in a directory called ``output``. More on output files in the :ref:`Output Documentation Section<documentation>`! 

``EUKulele`` is highly customizable and can be easily adapted to work with your dataset. For example, if you wanted to run ``BLAST`` against the ``PhyloDB`` database on your protein files that end in ``.protein`` this can be accomplished with the following command::

    EUKulele --sample_dir path/to/metatranscriptome/samples -m mets --protein_extension .protein --database phylodb --alignment_choice BLAST

A full list of parameters and customizations can be found here :ref:`Parameters Section<parameters>`.  Users might find it simpler  to modify a config file such as that provided `here <https://github.com/AlexanderLabWHOI/EUKulele/blob/master/config.yaml>`_). A config file can be passed to ``EUKulele`` as follows. A config file will be automatically generated after the first run of ``EUKulele`` in you working directory to facilitate re-running:: 

    EUKulele --config config.yaml

.. note::
    It is feasible to run ``EUKulele`` on metagenome derived contigs (not MAGs, those are discussed in detail below). If you wish to analyze metagenomic contigs as is described above for metatranscriptomes, we **STRONGLY RECOMMEND** that you provide predicted proteins from your metagenome rather than the nucleotide sequences from your metagenomic assembly. 
    Metagenomic contigs often consist of many proteins. ``EUKulele`` can predict proteins from metatranscriptomes with ``TransDecoder``, but this is **NOT** advised for metagenomic contigs. Additionally, a ``blastx``-style search will no be optimal for metagenomic contigs. It is up to the user to provide predicted proteins from their metagenomic contigs, as this can be a complex process (particularly for eukaryotic metagenomes) and is not within the scope of ``EUKulele``.
    Note that this is a problem in particular for eukaryotic sequences, where protein calling is more complicated due to the presence of introns.

Metagenome Assembled Genomes (MAGs)
===================================
In the second case, metagenome assembled genomes (shortened in ``EUKulele`` to ``mags``), are assumed to be the predicted proteins from a MAG derived from the binning of like contigs from a metagenome assembly based on tetranucleotide frequency, abundance, etc. as done by programs like `CONCOCT <https://github.com/BinPro/CONCOCT>`_ or `metaBAT <https://bitbucket.org/berkeleylab/metabat>`_. Please note that we only recommend using predicted proteins for the taxonomic annotation of MAGs and not raw nucleotide sequences (see above note for some justification). 

The most basic invocation of a MAG based analysis with EUKulele is the following::

    EUKulele --sample_dir path/to/MAGs -m mags

This will align all MAG proteins against the ``MMETSP`` with ``DIAMOND``. As detailed above in the section on metatranscriptomes, the parameters (see :ref:`Parameters Section<parameters>`) or config file can be used to adjust the running of ``EUKulele`` to suit the user. 

Of particular interest to the user might be the ``--tax_cutoffs`` parameter and associated file (``tax-cutoffs.yaml``). The default ``tax-cutoffs.yaml`` file details the percent id cutoffs that are used to assign taxonomy at various levels::

    species: 95
    genus: 80
    family: 65
    order: 50
    class: 30

These cutoffs are based on commonly used cutoffs in the literature, but can be modified by the user to suit their purposes. Using the above example, proteins that align with 95% identity to a reference protein will be assigned at the level of the species, those aligning with 65% identity at the level of family, etc. 

As with the metatranscriptome analysis, the taxonomic estimation of every protein will be output in the ``taxonomy_estimation`` folder within the output directory. Additionally, however, the taxonomic consensus across all proteins within a MAG will also be returned for MAGs. 

.. list-table:: Example ``max-level-mag`` output file
   :widths: 25 25 25
   :header-rows: 1

   * - 
     - max_taxa
     - proportion_id
   * - supergroup 
     - Eukaryota      
     -  0.9998
   * - division 
     - Archaeplastida      
     - 0.9898
   * - class   
     - Chlorophyta  
     - 0.9878
   * - family   
     - Mamiellales  
     - 0.9764
   * - genus   
     - Mamiellaceae  
     - 0.9382
   * - species   
     - Micromonas sp. NEPCC29
     - 0.4002

The ``max-level-mag`` files detail the relative proportion of all the proteins within a MAG that agree at a particular level. The file reports each of the six level considered by default in ``EUKulele`` (more on this in the section on databases :ref:`here<databases>`). For each level, the taxa which recruited that most reads within that level is reported. For example, the majority of proteins in the division level were annotated as Archaeplastida. The proportion of proteins that are annotated as that max level are also reported. 

So, in the above example 99.98% of the proteins in the dataset have a best hit to the supergroup level Eukaryota, meaning that the vast majority of the proteins had the same annotation at the supergroup level. This is largely true, where all proteins are annotated consistently (>90%) from supergroup to genus. However, only 40% of the proteins annotated consistently at the species level. It is up to the user to decide where and how they want to make a final taxonomic annotation for their MAG. In the above example, one might choose to annotate with confidence to the level of genus given the universally high consensus across proteins.

LCA Algorithm
=============

In some cases, multiple hits from alignment via ``blast`` or ``diamond`` will be reported and will meet the threshold specified by the user (see the :ref:`Parameters Section<parameters>`). In this case, the hits available at each taxonomic level will be evaluated using a simple Last Common Ancestor (LCA) algorithm. This simple implementation of the algorithm accepts input from the user (detailed in the :ref:`Parameters Section<parameters>`; parameter is ``--consensus_cutoff`` and has default of 0.75/75%) on what percentage of alignment-derived annotations need to be identical in order for the annotation to be adopted. If, for instance, only 50% of alignment hits match at the species level, less specific taxonomic levels are assessed until a 75% consensus is reached. For example, if two of four hits have the same species annotation, but all four hits have the same genus annotation, the genus annotation would be used, even if all hits meet the defined percentage identity threshold for the species level. 

LCA, while a robust annotation approach, is not the only means of predicting taxonomic level. We are currently exploring adding a phylogenetic estimate of eukaryotic taxonomy, particularly for the taxonomic placement of MAGs.
API 
===

def eukulele(string_arguments="", config=""):

    :param string_arguments: a space-separated string containing the ``EUKulele`` parameters and appropriate flags
    :param config: the path to a configuration file
    :type string_arguments: str
    :type config: str