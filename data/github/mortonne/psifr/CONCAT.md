# Psifr
[![PyPI version](https://badge.fury.io/py/psifr.svg)](https://badge.fury.io/py/psifr)
[![Documentation Status](https://readthedocs.org/projects/psifr/badge/?version=latest)](https://psifr.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.com/mortonne/psifr.svg?branch=master)](https://travis-ci.com/mortonne/psifr)
[![codecov](https://codecov.io/gh/mortonne/psifr/branch/master/graph/badge.svg)](https://codecov.io/gh/mortonne/psifr)
[![status](https://joss.theoj.org/papers/712d4452e465229d61d0e281d3d6f299/status.svg)](https://joss.theoj.org/papers/712d4452e465229d61d0e281d3d6f299)
[![DOI](https://zenodo.org/badge/248593723.svg)](https://zenodo.org/badge/latestdoi/248593723)

Advanced analysis and visualization of free recall data in Python.

In free recall, participants study a list of items and then name all of the items they can remember in any order they choose. Many sophisticated analyses have been developed to analyze data from free recall experiments, but these analyses are often complicated and difficult to implement. 

Psifr leverages the Pandas data analysis package to make precise and flexible analysis of free recall data faster and easier. The name Psifr is pronounced "cipher". It's taken from Psi, in reference to the field of psychology, and FR for free recall.

<div align="center">
  <div style="max-width:500px; margin:0 20px;">
    <img src="https://github.com/mortonne/psifr/blob/master/images/raster.png" alt="free recall visualization" width="500px">
    <div style="text-align:left; padding:10px 0;">
      Raster plot showing each recall in a free recall experiment. Color indicates serial position; yellow items were presented late in the list, while purple items were presented at the beginning. Magenta squares indicate intrusions of items that were not presented during the study list.
    </div>
  </div>
</div>

## Installation

You can install the latest stable version of Psifr using pip:

```bash
pip install psifr
```

You can also install the development version directly from the code
repository on GitHub:

```bash
pip install git+git://github.com/mortonne/psifr
```

## Quickstart

To plot a serial position curve for a sample dataset:

```python
from psifr import fr
df = fr.sample_data('Morton2013')
data = fr.merge_free_recall(df)
recall = fr.spc(data)
g = fr.plot_spc(recall)
```

See the [user guide](https://psifr.readthedocs.io/en/latest/guide/overview.html) for detailed documentation on importing and analyzing free recall datasets.

Also see the Jupyter notebooks for more analysis examples:
* [Recall performance](https://github.com/mortonne/psifr-notebooks/blob/master/demo_recall.ipynb)
* [Temporal clustering](https://github.com/mortonne/psifr-notebooks/blob/master/demo_lag_crp.ipynb)

## Importing data

Generally the best way to get your data into shape for analysis in Psifr is to create a CSV file with one row for each event in the experiment, including study events (i.e., item presentations) and all recall attempts (including repeats and intrusions). See [importing data](https://psifr.readthedocs.io/en/latest/guide/import.html) for details.

If you have data in the standard [EMBAM](https://github.com/vucml/EMBAM) format, use `scripts/frdata2table.m` to convert your data struct to a table with standard format. Then use the Matlab function `writetable` to write a CSV file which can then be read into Python for analysis.

## Related projects

### EMBAM
Analyses supported by Psifr are based on analyses implemented in the Matlab toolbox [EMBAM](https://github.com/vucml/EMBAM).

### pybeh
[pybeh](https://github.com/pennmem/pybeh) is a direct Python port of EMBAM that supports a wide range of analyses.

### Quail
[Quail](https://github.com/ContextLab/quail) runs automatic scoring of free recall data, supports calculation and plotting of some common free recall measures, and has tools for measuring the "memory fingerprint" of individuals.

## Contributing to Psifr

Contributions are welcome to suggest new features, add documentation, and identify bugs. See the [contributing guidelines](.github/CONTRIBUTING.md) for an overview. 
---
title: 'Psifr: Analysis and visualization of free recall data'
tags:
  - Python
  - psychology
  - memory
  - research
authors:
  - name: Neal W Morton
    orcid: 0000-0002-2631-2710
    affiliation: 1
affiliations:
 - name: Research Fellow, Center for Learning and Memory, The University of Texas at Austin
   index: 1
date: 24 August 2020
bibliography: paper.bib
---

# Summary

Research on human memory has been strongly influenced by data from free 
recall experiments, wherein participants study a list of items (such as 
words) and then freely recall them in any order they wish [@Murdock:1962]. 
Free recall provides an extremely rich dataset that not only reflects 
which items were recalled but also the order in which they were recalled. 
However, analysis of free recall data is difficult, as many different influences
on recall must be taken into account. 
For example, one influential analysis, conditional response probability 
as a function of lag (lag-CRP), has been used to measure the tendency of 
participants to successively recall items that were originally presented 
near to each other in time [@Kahana:1996].
This analysis requires taking into account the items that are still
available for recall at each transition between recalled items. 
This analysis may need to be made conditional on other factors, such as
the category of the items being recalled [@Polyn:2011; @Morton:2017], 
thus complicating the analysis further.

`Psifr` (pronounced like "cypher") was developed to consolidate a number 
of free recall analysis methods (often implemented in MATLAB) within a flexible 
Python package. 
The `Psifr` package includes core utilities that simplify
and standardize a number of different analyses of recall sequences,
including analyses focused on serial position [@Murdock:1962],
temporal order [@Kahana:1996; @Polyn:2011], 
stimulus category [@Polyn:2009; @Morton:2016], and the semantic meaning 
of presented items [@Howard:2002]. 
The core utilities are also designed to facilitate implementation of 
extensions to tailor analyses for specific experiments.

![Example of a serial position curve showing the probability of recalling 
an item based on its position in the list. Plots may be flexibly divided by 
condition using grouping semantics supported by Seaborn. In this case,
different list types (mixed-category or pure-category) are plotted as separate
curves.\label{fig:spc}](spc_list_type.pdf)

# Statement of Need

Existing packages for analysis of free recall data include `EMBAM`
and `Quail`. 
[EMBAM](https://github.com/vucml/EMBAM), which is a fork of the 
[Behavioral Toolbox](http://memory.psych.upenn.edu/Behavioral_toolbox), 
is implemented in MATLAB, making it difficult to use with the extensive 
data science ecosystem in Python. 
The [pybeh](https://github.com/pennmem/pybeh) package is a Python port 
of `EMBAM` written using `numpy`.
As it is a fairly direct port of EMBAM, `pybeh` does not make use of some of 
the advantages  of Python, such as the advanced data science packages of 
`Pandas` and `Seaborn` [@Reback:2020; @Waskom:2020].
`Quail`, a Python package, provides some similar functionality to `Psifr`,
including analysis of recall order [@Heusser:2017]. 
However, while `Quail` uses a custom data structure to store free 
recall sequences, `Psifr` uses `Pandas` `DataFrame` objects. 
This design makes it possible for the user to make full use 
of the split-apply-combine operations of `Pandas` to quickly run complex analyses.

![Serial position curve split by list type, with a separate panel for each
participant in an experiment.\label{fig:spc_subject}](spc_subject.pdf)

Similarly, `Psifr` makes available the full power of the `Seaborn` 
visualization package to provide expressive visualization capabilities. 
The plotting functions in `Psifr` allow the user to easily view analysis 
results in different ways; for example, an analysis of recall by serial 
position can be visualized either as a single plot with error bars 
(\autoref{fig:spc}) or as a grid of individual plots for each 
participant in the experiment (\autoref{fig:spc_subject}).
`Psifr` also supports creation of raster plots 
(\autoref{fig:raster}), a method for visualizing whole 
free recall datasets to facilitate quick discovery of patterns in
the order of recalls [@Romani:2016].

![Raster plot displaying the order of every recall for one participant.
Each marker indicates one recall, and the color of the marker reflects
the serial position of the recalled item. Red markers indicate intrusions
of items not on the studied list.\label{fig:raster}](raster.pdf)

`Psifr` was designed to be used by memory researchers and students.
It is currently being used in two ongoing projects that require advanced
analysis and visualization. 
The interface is designed to simplify common tasks while also allowing 
for substantial customization to facilitate analysis of specific episodic
memory experiments.
Advanced visualization further helps to support better understanding of 
complex datasets. 
The source code for `Psifr` has  been archived to Zenodo with the linked DOI: 
[10.5281/zenodo.4086188](https://doi.org/10.5281/zenodo.4086187).

# Acknowledgements

Thanks to Sean Polyn for helpful discussions.
`Psifr` is inspired by functions initially developed for `EMBAM`,
which was developed by Sean Polyn, Neal Morton, Richard Lawrence,
Joshua McCluey, Karl Healey, Lynn Lohnas, and Jin Jeon. 
`Psifr` also takes inspiration from the `Quail` package developed 
by Andrew Heusser, Paxton Fitzpatrick, Campbell Field, Kirsten Ziman, 
and Jeremy Manning.

# References# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
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

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at mortonne@gmail.com. All
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
# Contributing to Psifr

Whether you are a novice or experienced software developer, 
all contributions and suggestions are welcome!

## Getting Started

If you are looking to contribute to the *Psifr* codebase, the 
best place to start is the [GitHub "issues" tab](https://github.com/mortonne/psifr/issues). 
This is also a great place for filing bug reports and making suggestions for 
ways in which the code and documentation can be improved.

## Filing Issues

If you notice a bug in the code or documentation, or have suggestions for 
how we can improve either, feel free to create an issue on the 
[GitHub "issues" tab](https://github.com/mortonne/psifr/issues) using 
[GitHub's "issue" form](https://github.com/mortonne/psifr/issues/new). 
The form contains some questions to help clarify the issue.

## Contributing to the Codebase

The code is hosted on [GitHub](https://www.github.com/mortonne/psifr), so 
you will need to use [Git](https://git-scm.com/) to clone the project and 
make changes to the codebase. 

If you do not have edit access to the original repository (mortonne/psifr), 
then first make a fork of the repository. Next, create a branch (either within 
the original repository or your fork) with a name indicating what you'll be 
working on. Your changes from the main branch should all be committed to 
this branch.

Once you have created a branch on GitHub (and a fork if necessary), 
you should create a development environment using Conda or virtualenv 
that is separate from your existing Python environment so that you can make 
and test changes without compromising your own work environment. To
install code locally for development, clone the project, checkout your
branch, and run `pip install -e .` from the main project directory. That
will install Psifr in editable mode so changes to the codebase will be
reflected the next time you import (or reload) Psifr.

Unit tests will be run automatically when you open a pull request. 
You can also run the test suite manually on your machine by installing
[pytest](https://docs.pytest.org/en/stable/), then running `pytest` from
the main project directory.

Please run [black](https://black.readthedocs.io/en/stable/) on any modules 
you edit so that your code will have standard formatting. For example, run 
`black -S src/psifr/fr.py` to reformat the `fr.py` module. The `-S` flag
indicates that strings should be left alone, rather than being reformatted
to use double quotes.

When making commits to your branch, please format your commit messages
following the
[conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) 
standard. This make it easier to see what kinds of changes are being
proposed (e.g., `feat` for an added feature or `fix` for a proposed fix).

Once your changes are ready to be submitted, push your changes 
to GitHub and then create a pull request from your branch. Your code will 
be reviewed, and you may be asked to make changes before your code is 
approved to be merged into the master branch. 
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
1. Load/create test data (can use `raw = fr.sample_data('Morton2013')` to load project sample data)
2. Run command on test data
3. See error/unexpected output

**Expected behavior**
A clear and concise description of what you expected to happen.

**Error message**
If applicable, copy the full stack trace from the error message.

**Environment (please complete the following information):**
 - Python version [e.g. 3.8]
 - Psifr version [e.g. 0.5.1]

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
A clear and concise description of what the problem is.

**Describe the solution you'd like**
A clear and concise description of what you want to happen. If proposing to add an analysis or visualization, provide a reference to a publication if applicable.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context about the feature request here.
