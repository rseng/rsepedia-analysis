---
title: 'Parent-map: analysis of parental contributions to evolved or engineered protein or DNA sequences'
tags:
  - Python
  - Directed evolution
  - AAV
  - Sequence analysis
authors:
  - name: Damien Marsic
    orcid: 0000-0003-0847-8095
    affiliation: 1
affiliations:
 - name: Porton Biologics, 388 Xinping Street, Suzhou Industrial Park, Jiangsu 215021, China
   index: 1
date: 18 November 2020
bibliography: paper.bib
---

# Summary

Parent-map analyzes protein or DNA sequences which are derived from one or multiple parent sequences, and shows parental contributions as well as differences from relevant parents. Originally developed to analyze capsid protein sequences obtained by directed evolution, parent-map can be used in any case where variant sequences are to be compared to parent sequences from which they are derived. Parent-map detects sequence shuffling as well as substitutions, insertions and deletions, and displays results in user-friendly formats. Parent-map is an open-source, platform-independent Python 3 script, available as a Bioconda package as well as a Windows program.

Source code: <https://github.com/damienmarsic/Parent-map>

Python package: <https://pypi.org/project/parent-map/>

Bioconda recipe and package: <http://bioconda.github.io/recipes/parent-map/README.html>

Windows installer: <https://sourceforge.net/projects/parent-map/>

Documentation: <https://parent-map.readthedocs.io/>

# Statement of need

Adeno-associated virus (AAV) capsid directed evolution projects typically generate multiple enriched variant sequences after 2 to 5 rounds of selection starting from complex capsid libraries. For libraries developed from a single parental serotype, through random peptide insertion at a specific position or surface loop diversification in well-defined variable regions for example, a single multiple alignment of all enriched variant sequences against the parent sequence conveniently shows how each variant differs from the parent. However, when more than one parental sequence is involved, such as when different libraries are mixed together, or when a library design involves DNA shuffling from several parents, such alignments can quickly become illegible, particularly when the complete capsid gene is sequenced. In such cases, in the absence of appropriate software tools, each variant needs to be separately aligned against all possible parents, a time-consuming and cumbersome process. An added difficulty in the case of shuffled libraries is that, because of high sequence homology between parents, multiple regions will share sequence identities with more than one parent, complicating attempts at comprehensively defining the variant sequences in terms of parental contributions. To date, SALANTO [@herrmann_robust_2019] seems to be the only relevant publicly available software. However, it only applies to shuffled libraries, and its user-friendliness is limited as it requires the user to perform a multiple sequence alignment beforehand, and to further process the data manually after analysis. The software described in this article, parent-map, provides a user-friendly and comprehensive solution. It can be used with sequences derived from any type of library, or even with naturally-occurring mutants or rationally engineered variants. It is not limited to protein sequences. It only requires one file containing the variant sequences to be analyzed, and one file containing parental sequences, without any prior manipulation. It generates a set of five files covering most end-users’ needs, in directly usable formats. Finally, although it was developed to address a need in the field of AAV capsid directed evolution, parent-map can be used whenever protein or DNA sequences, whether originating from natural evolution, directed evolution or rational design, are to be compared with one or more possible parental sequences.

# Methods

Parent-map was written under Python 3.7 as both a command-line interface (CLI) and a graphical user interface (GUI) application, by allowing parser modules [argparse](https://docs.python.org/3/library/argparse.html) and [Gooey](https://github.com/chriskiehl/Gooey) to coexist within a single file (the GUI will start if no argument is present, while any argument will cause parent-map to start in CLI mode). A parent-map Python package was created and uploaded to the Python Package Index (PyPI) according to [packaging instructions](https://packaging.python.org/tutorials/packaging-projects/). A parent-map Bioconda [@gruning_bioconda_2018] recipe based on the PyPI package was written and submitted according to [instructions](https://bioconda.github.io/contributor/index.html). A stand-alone Windows executable and its installation program were created using respectively [PyInstaller](https://www.pyinstaller.org/) and [Inno Setup](https://jrsoftware.org/isinfo.php). The documentation was written using [Sphinx](https://www.sphinx-doc.org/en/master/).

# Implementation

Parent-map is a platform-independent Python script that generates a set of five output files from two input files. Input file names and options can be entered as arguments at launch time, resulting in parent-map running in CLI mode, or within the GUI, which starts if parent-map is launched without arguments. This flexibility allows parent-map to be deployed in a variety of settings, as a simple desktop application or even as a bioinformatics pipeline component. The first input file contains the variant sequences, typically the most frequent or the most enriched sequences obtained at the completion of a directed evolution experiment. The other input file is a set of potential parental sequences to the variant sequences. The most useful files generated by parent-map, particularly in the case of variants derived from DNA shuffling, are parental contribution maps (file names ending in –par.txt and –par.html, the latter being a colorized version of the former). Instead of all possible combinations, the simplest map that can accurately describe the variant is shown, using as few parents and as few fragments as possible. Other output files include a statistics file summarizing the variant sequences main features, a sequence definition file comprehensively defining each variant in terms of its parents, and an alignment file showing how variants differ from their common parent.

Parent-map can be tested using the provided [variant](https://github.com/damienmarsic/Parent-map/blob/master/example_variants.fasta) and [parent](https://github.com/damienmarsic/Parent-map/blob/master/example_parents.fasta) sample files, based on available literature describing evolved and rationally designed AAV capsid variants. Variants AAV-DJ [@grimm_vitro_2008], AAV2.5T [@excoffon_directed_2009], NP84 [@paulk_bioengineered_2018] and OLIG001 [@powell_characterization_2016] are derived from shuffled DNA libraries. Variants AAV-F [@hanlon_selection_2019], AAV-PHP.B [@deverman_cre-dependent_2016], 7m8 [@dalkara_vivo-directed_2013] and rAAV2-retro [@tervo_designer_2016] are derived from peptide insertion libraries. Variants SCH2, SCH9 [@ojala_vivo_2018], LI-A and LI-C [@marsic_vector_2014] are derived from more complex rationally designed libraries. Variants AAV2i8 [@asokan_reengineering_2010] and AAV2-sept-Y-F [@petrs-silva_novel_2011] were rationally designed. Using default settings, parent-map correctly identifies single parental contributions from AAV9 for variants AAV-F and AAV-PHP.B, single parental contributions from AAV2 for variants 7m8, rAAV2-retro, LI-A, LI-C, AAV2-sept-Y-F, and multiple parental contributions from AAV2, AAV8 and AAV9 for AAV-DJ, from AAV2 and AAV5 for AAV2.5T, from AAV2, AAV3B and AAV6 for NP84, from AAV2, AAV6, AAV8 and AAV9 for OLIG001, SCH2 and SCH9, and from AAV2 and AAV8 for AAV2i8. Parent-map also correctly detects peptide insertions FVVGQSY for AAV-F and TLAVPFK for AAV-PHP.B, both at position 588, and peptide insertions LALGETTRPA for 7m8 and LADQDYTKTA for rAAV2-retro, both at position 587. Finally, parent-map correctly identifies substitutions A to T at position 457 for AAV-DJ and at position 582 for AAV2.5T, substitutions K to E at 532 and R to G at 585 for NP84, E to K substitution at 532 and unmatched H at 726 for OLIG001, substitutions I to T at 240 and V to I at 718 for 7m8, substitutions N to D at 382 and V to I at 718 for rAAV2-retro, the 14 and 4 substitutions for LI-A and LI-C respectively, as well as the 7 Y to F substitutions at 252, 272, 444, 500, 700, 704 and 730 for AAV2-sept-Y-F.   

A comprehensive description of parent-map is provided in the [documentation](https://parent-map.rtfd.io).

# Acknowledgements

We thank Yan Chen and Oleksandr Kondratov for testing parent-map and providing valuable feedback.

# References
# Parent-map

Analyze parental contributions to evolved or engineered protein or DNA sequences

[Read the JOSS article](https://joss.theoj.org/papers/10.21105/joss.02864)

Run parent-map without arguments to start in GUI mode.
Run parent-map with arguments to start in console mode.
Run parent-map.py -h for help and list of arguments.

Install from Bioconda:
```
conda install -c bioconda parent-map
```

Install from PyPI:
````
pip install parent-map
````
Note that you may need to individually install some dependencies. If you get an error message about missing modules when trying to run parent-map after installing it with pip, install them using pip install. Example: pip install gooey


Run:
```
python -m parent-map
```
Windows installer: https://sourceforge.net/projects/parent-map/

Full documentation: https://parent-map.readthedocs.io/

GitHub issues (bug reports, feature requests): https://github.com/damienmarsic/Parent-map/issues/new/choose

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/parent-map/README.html)

[![Download parent-map](https://a.fsdn.com/con/app/sf-download-button)](https://sourceforge.net/projects/parent-map/files/latest/download)

[![Parent-map documentation](https://img.shields.io/badge/Parent--map-Documentation-yellow)](https://parent-map.readthedocs.io/)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02864/status.svg)](https://doi.org/10.21105/joss.02864)

<div itemscope itemtype="https://schema.org/Person"><a itemprop="sameAs" content="https://orcid.org/0000-0003-0847-8095" href="https://orcid.org/0000-0003-0847-8095" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon">https://orcid.org/0000-0003-0847-8095</a></div>
# Contributing to parent-map

## How can I contribute?

### Reporting Bugs

When you are creating a bug report, please include as many details as possible. Fill out [the required template](https://github.com/damienmarsic/Parent-map/blob/master/.github/ISSUE_TEMPLATE/bug_report.md), the information it asks for helps us resolve issues faster.

> **Note:** If you find a **Closed** issue that seems like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one.

#### How do I submit a bug report?

Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue and provide the following information by filling in [the template](https://github.com/damienmarsic/Parent-map/blob/master/.github/ISSUE_TEMPLATE/bug_report.md).

Explain the problem and include additional details to help maintainers reproduce the problem.


### Feature requests

This section guides you through submitting improvements to parent-map, including completely new features and minor improvements to existing functionality. Following these guidelines helps maintainers and the community understand your suggestion :pencil: and find related suggestions :mag_right:.

#### How do I submit a feature request?

Feature requests are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue and provide the following information:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Include screenshots and animated GIFs** which help you demonstrate the steps or point out the part of parent-map which the suggestion is related to. You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast) or [this tool](https://github.com/GNOME/byzanz) on Linux.
* **Explain why this enhancement would be useful** to parent-map users
* **Specify which version of parent-map you are using.** You can get the exact version by running `parent-map -v` in your terminal.
* **Specify the name and version of the OS you're using.**

---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
