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
.. Psifr documentation master file, created by
   sphinx-quickstart on Sun Mar 29 17:31:18 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Psifr documentation
===================

In free recall, participants study a list of items and then name all of the
items they can remember in any order they choose. Many sophisticated analyses
have been developed to analyze data from free recall experiments, but these
analyses are often complicated and difficult to implement.

Psifr leverages the Pandas data analysis package to make precise and flexible
analysis of free recall data faster and easier.

See the code repository for version `release notes`_.

.. toctree::
   :maxdepth: 2

   /usage/installation
   /guide/overview
   /tutorials/overview
   /api/overview
   /development/overview

.. _release notes: https://github.com/mortonne/psifr/releases

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   import matplotlib as mpl
   import matplotlib.pyplot as plt

   plt.style.use('default')
   mpl.rcParams['axes.labelsize'] = 'large'
   mpl.rcParams['savefig.bbox'] = 'tight'
   mpl.rcParams['savefig.pad_inches'] = 0.1

   pd.options.display.max_rows = 15

Comparing conditions
====================

When analyzing a dataset, it's often important to compare different
experimental conditions. Psifr is built on the Pandas DataFrame, which
has powerful ways of splitting data and applying operations to it.
This makes it possible to analyze and plot different conditions using
very little code.

.. _custom-columns:

Working with custom columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, load some sample data and create a merged DataFrame:

.. ipython:: python

    from psifr import fr
    df = fr.sample_data('Morton2013')
    data = fr.merge_free_recall(
        df, study_keys=['category'], list_keys=['list_type']
    )
    data.head()

The :py:func:`~psifr.fr.merge_free_recall` function only includes columns from the
raw data if they are one of the standard columns or if they've explictly been
included using :code:`study_keys`, :code:`recall_keys`, or :code:`list_keys`.
:code:`list_keys` apply to all events in a list, while :code:`study_keys` and
:code:`recall_keys` are relevant only for study and recall events, respectively.

We've included a list key here, to indicate that the :code:`list_type`
field should be included for all study and recall events in each list, even
intrusions. The :code:`category` field will be included for all study events
and all valid recalls. Intrusions will have an undefined category.

Analysis by condition
~~~~~~~~~~~~~~~~~~~~~

Now we can run any analysis separately for the different conditions. We'll
use the serial position curve analysis as an example.

.. ipython:: python

    spc = data.groupby('list_type').apply(fr.spc)
    spc.head()

The :code:`spc` DataFrame has separate groups with the results for each
:code:`list_type`.

.. warning::

    When using :code:`groupby` with order-based analyses like
    :py:func:`~psifr.fr.lag_crp`, make sure all recalls in all recall
    sequences for a given list have the same label. Otherwise, you will
    be breaking up recall sequences, which could result in an invalid
    analysis.

Plotting by condition
~~~~~~~~~~~~~~~~~~~~~

We can then plot a separate curve for each condition. All plotting functions
take optional :code:`hue`, :code:`col`, :code:`col_wrap`, and :code:`row`
inputs that can be used to divide up data when plotting. See the
`Seaborn documentation <https://seaborn.pydata.org/generated/seaborn.relplot.html>`_
for details. Most inputs to :py:func:`seaborn.relplot` are supported.

For example, we can plot two curves for the different list types:

.. ipython:: python

    @savefig spc_list_type.svg
    g = fr.plot_spc(spc, hue='list_type').add_legend()

We can also plot the curves in different axes using the :code:`col` option:

.. ipython:: python

    @savefig spc_list_type_col.svg
    g = fr.plot_spc(spc, col='list_type')

We can also plot all combinations of two conditions:

.. ipython:: python

    spc_split = data.groupby(['list_type', 'category']).apply(fr.spc)
    @savefig spc_split.svg
    g = fr.plot_spc(spc_split, col='list_type', row='category')

Plotting by subject
~~~~~~~~~~~~~~~~~~~

All analyses can be plotted separately by subject. A nice way to do this is
using the :code:`col` and :code:`col_wrap` optional inputs, to make a grid
of plots with 6 columns per row:

.. ipython:: python

    @savefig spc_subject.svg
    g = fr.plot_spc(
        spc, hue='list_type', col='subject', col_wrap=6, height=2
    ).add_legend()

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   import matplotlib as mpl
   import matplotlib.pyplot as plt

   plt.style.use('default')
   mpl.rcParams['axes.labelsize'] = 'large'
   mpl.rcParams['savefig.bbox'] = 'tight'
   mpl.rcParams['savefig.pad_inches'] = 0.1

   pd.options.display.max_rows = 15

============
Recall order
============

A key advantage of free recall is that it provides information not only about
what items are recalled, but also the order in which they are recalled. A
number of analyses have been developed to charactize different influences on
recall order, such as the temporal order in which the items were presented at
study, the category of the items themselves, or the semantic similarity between
pairs of items.

Each conditional response probability (CRP) analysis involves calculating the
probability of some type of transition event. For the lag-CRP analysis,
transition events of interest are the different lags between serial positions
of items recalled adjacent to one another. Similar analyses focus not on
the serial position in which items are presented, but the properties of the
items themselves. A semantic-CRP analysis calculates the probability of
transitions between items in different semantic relatedness bins. A special
case of this analysis is when item pairs are placed into one of two bins,
depending on whether they are in the same stimulus category or not. In Psifr,
this is referred to as a category-CRP analysis.

Lag-CRP
~~~~~~~

In all CRP analyses, transition probabilities are calculated conditional
on a given transition being available. For example, in a six-item list,
if the items 6, 1, and 4 have been recalled, then possible items that could
have been recalled next are 2, 3, or 5; therefore, possible lags at
that point in the recall sequence are -2, -1, or +1. The number of actual
transitions observed for each lag is divided by the number of times that
lag was possible, to obtain the CRP for each lag.

First, load some sample data and create a merged DataFrame:

.. ipython:: python

    from psifr import fr
    df = fr.sample_data('Morton2013')
    data = fr.merge_free_recall(df, study_keys=['category'])

Next, call :py:func:`~psifr.fr.lag_crp` to calculate conditional response
probability as a function of lag.

.. ipython:: python

    crp = fr.lag_crp(data)
    crp

The results show the count of times a given transition actually happened
in the observed recall sequences (:code:`actual`) and the number of times a
transition could have occurred (:code:`possible`). Finally, the :code:`prob` column
gives the estimated probability of a given transition occurring, calculated
by dividing the actual count by the possible count.

Use :py:func:`~psifr.fr.plot_lag_crp` to display the results:

.. ipython:: python

   @savefig lag_crp.svg
   g = fr.plot_lag_crp(crp)

The peaks at small lags (e.g., +1 and -1) indicate that the recall sequences
show evidence of a temporal contiguity effect; that is, items presented near
to one another in the list are more likely to be recalled successively than
items that are distant from one another in the list.

Lag rank
~~~~~~~~

We can summarize the tendency to group together nearby items using a lag
rank analysis. For each recall, this determines the absolute lag of all
remaining items available for recall and then calculates their percentile
rank. Then the rank of the actual transition made is taken, scaled to vary
between 0 (furthest item chosen) and 1 (nearest item chosen). Chance
clustering will be 0.5; clustering above that value is evidence of a
temporal contiguity effect.

.. ipython:: python

    ranks = fr.lag_rank(data)
    ranks
    ranks.agg(['mean', 'sem'])

Category CRP
~~~~~~~~~~~~

If there are multiple categories or conditions of trials in a list, we
can test whether participants tend to successively recall items from the
same category. The category-CRP estimates the probability of successively
recalling two items from the same category.

.. ipython:: python

    cat_crp = fr.category_crp(data, category_key='category')
    cat_crp
    cat_crp[['prob']].agg(['mean', 'sem'])

The expected probability due to chance depends on the number of
categories in the list. In this case, there are three categories, so
a category CRP of 0.33 would be predicted if recalls were sampled
randomly from the list.

Restricting analysis to specific items
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you may want to focus an analysis on a subset of recalls. For
example, in order to exclude the period of high clustering commonly
observed at the start of recall, lag-CRP analyses are sometimes
restricted to transitions after the first three output positions.

You can restrict the recalls included in a transition analysis using
the optional :code:`item_query` argument. This is built on the Pandas
query/eval system, which makes it possible to select rows of a
:code:`DataFrame` using a query string. This string can refer to any
column in the data. Any items for which the expression evaluates to
:code:`True` will be included in the analysis.

For example, we can use the :code:`item_query` argument to exclude any
items recalled in the first three output positions from analysis. Note
that, because non-recalled items have no output position, we need to
include them explicitly using :code:`output > 3 or not recall`.

.. ipython:: python

    crp_op3 = fr.lag_crp(data, item_query='output > 3 or not recall')
    @savefig lag_crp_op3.svg
    g = fr.plot_lag_crp(crp_op3)

Restricting analysis to specific transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In other cases, you may want to focus an analysis on a subset of
transitions based on some criteria. For example, if a list contains
items from different categories, it is a good idea to take this into
account when measuring temporal clustering using a lag-CRP analysis.
One approach is to separately analyze within- and across-category
transitions.

Transitions can be selected for inclusion using the optional
:code:`test_key` and :code:`test` inputs. The :code:`test_key`
indicates a column of the data to use for testing transitions; for
example, here we will use the :code:`category` column. The
:code:`test` input should be a function that takes in the test value
of the previous recall and the current recall and returns True or False
to indicate whether that transition should be included. Here, we will
use a lambda (anonymous) function to define the test.

.. ipython:: python

    crp_within = fr.lag_crp(data, test_key='category', test=lambda x, y: x == y)
    crp_across = fr.lag_crp(data, test_key='category', test=lambda x, y: x != y)
    crp_combined = pd.concat([crp_within, crp_across], keys=['within', 'across'], axis=0)
    crp_combined.index.set_names('transition', level=0, inplace=True)
    @savefig lag_crp_cat.svg
    g = fr.plot_lag_crp(crp_combined, hue='transition').add_legend()

The :code:`within` curve shows the lag-CRP for transitions between
items of the same category, while the :code:`across` curve shows
transitions between items of different categories.

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   import matplotlib as mpl
   import matplotlib.pyplot as plt

   plt.style.use('default')
   mpl.rcParams['axes.labelsize'] = 'large'
   mpl.rcParams['savefig.bbox'] = 'tight'
   mpl.rcParams['savefig.pad_inches'] = 0.1

   pd.options.display.max_rows = 15

==================
Recall performance
==================

First, load some sample data and create a merged DataFrame:

.. ipython:: python

    from psifr import fr
    df = fr.sample_data('Morton2013')
    data = fr.merge_free_recall(df)

Raster plot
~~~~~~~~~~~

Raster plots can give you a quick overview of a whole dataset. We'll look at
all of the first subject's recalls. This will plot every individual recall,
colored by the serial position of the recalled item in the list. Items near
the end of the list are shown in yellow, and items near the beginning of the
list are shown in purple. Intrusions of items not on the list are shown in red.

.. ipython:: python

    subj = fr.filter_data(data, 1)
    @savefig raster_subject.svg
    g = fr.plot_raster(subj).add_legend()

Serial position curve
~~~~~~~~~~~~~~~~~~~~~

We can calculate average recall for each serial position
using :py:func:`~psifr.fr.spc` and plot using :py:func:`~psifr.fr.plot_spc`.

.. ipython:: python

    recall = fr.spc(data)
    @savefig spc.svg
    g = fr.plot_spc(recall)

Using the same plotting function, we can plot the curve for each
individual subject:

.. ipython:: python

    @savefig spc_indiv.svg
    g = fr.plot_spc(recall, col='subject', col_wrap=5)

Probability of Nth recall
~~~~~~~~~~~~~~~~~~~~~~~~~

We can also split up recalls, to test for example how likely participants
were to initiate recall with the last item on the list.

.. ipython:: python

    prob = fr.pnr(data)
    prob

This gives us the probability of recall by output position (:code:`'output'`)
and serial or input position (:code:`'input'`). This is a lot to look at all
at once, so it may be useful to plot just the first three output positions.
We can plot the curves using :py:func:`~psifr.fr.plot_spc`, which takes an
optional :code:`hue` input to specify a variable to use to split the data
into curves of different colors.

.. ipython:: python

    pfr = prob.query('output <= 3')
    @savefig pnr.svg
    g = fr.plot_spc(pfr, hue='output').add_legend()

This plot shows what items tend to be recalled early in the recall sequence.
Importing data
==============

In Psifr, free recall data are imported in the form of a "long" format
table. Each row corresponds to one *study* or *recall* event. Study
events include any time an item was presented to the participant.
Recall events correspond to any recall attempt; this includes *repeats*
of items there were already recalled and *intrusions* of items that
were not present in the study list.

This type of information is well represented in a CSV spreadsheet,
though any file format supported by pandas may be used for input. To
import from a CSV, use pandas. For example:

.. code-block:: python

    import pandas as pd
    data = pd.read_csv("my_data.csv")

Trial information
-----------------

The basic information that must be included for each event is the
following:

subject
    Some code (numeric or string) indicating individual participants.
    Must be unique for a given experiment. For example, ``sub-101``.

list
    Numeric code indicating individual lists. Must be unique within
    subject.

trial_type
    String indicating whether each event is a ``study`` event or a
    ``recall`` event.

position
    Integer indicating position within a given phase of the list. For
    ``study`` events, this corresponds to *input position* (also
    referred to as *serial position*). For ``recall`` events, this
    corresponds to *output position*.

item
    Individual thing being recalled, such as a word. May be specified
    with text (e.g., ``pumpkin``, ``Jack Nicholson``) or a numeric code
    (``682``, ``121``). Either way, the text or number must be unique
    to that item. Text is easier to read and does not require any
    additional information for interpretation and is therefore
    preferred if available.

Example
-------

.. csv-table:: Sample data
    :header: "subject", "list", "trial_type", "position", "item"
    :widths: 8, 8, 8, 8, 8

    1, 1, "study", 1, "absence"
    1, 1, "study", 2, "hollow"
    1, 1, "study", 3, "pupil"
    1, 1, "recall", 1, "pupil"
    1, 1, "recall", 2, "absence"

Additional information
----------------------

Additional fields may be included in the data to indicate other
aspects of the experiment, such as presentation time, stimulus
category, experimental session, distraction length, etc. All of
these fields can then be used for analysis in Psifr.
Scoring data
============

After :doc:`importing free recall data</guide/import>`, we have a DataFrame with
a row for each study event and a row for each recall event. Next, we need to
score the data by matching study events with recall events.

Scoring list recall
-------------------

First, let's create a simple sample dataset with two lists:

.. ipython:: python

    import pandas as pd
    data = pd.DataFrame({
        'subject': [
            1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1,
        ],
        'list': [
            1, 1, 1, 1, 1, 1,
            2, 2, 2, 2, 2, 2,
        ],
        'trial_type': [
            'study', 'study', 'study', 'recall', 'recall', 'recall',
            'study', 'study', 'study', 'recall', 'recall', 'recall',
        ],
        'position': [
            1, 2, 3, 1, 2, 3,
            1, 2, 3, 1, 2, 3,
        ],
        'item': [
            'absence', 'hollow', 'pupil', 'pupil', 'absence', 'empty',
            'fountain', 'piano', 'pillow', 'pillow', 'fountain', 'pillow',
        ],
    })
    data

Next, we'll merge together the study and recall events by matching up
corresponding events:

.. ipython:: python

    from psifr import fr
    merged = fr.merge_free_recall(data)
    merged

For each item, there is one row for each unique combination of input and
output position. For example, if an item is presented once in the list, but
is recalled multiple times, there is one row for each of the recall attempts.
Repeated recalls are indicated by the `repeat` column, which is greater than
zero for recalls of an item after the first. Unique study events are indicated
by the `study` column; this excludes intrusions and repeated recalls.

Items that were not recalled have the `recall` column set to `False`. Because
they were not recalled, they have no defined output position, so `output` is
set to `NaN`. Finally, intrusions have an output position but no input position
because they did not appear in the list. There is an `intrusion` field for
convenience to label these recall attempts.

:py:func:`~psifr.fr.merge_free_recall` can also handle additional attributes beyond
the standard ones, such as codes indicating stimulus category or list condition.
See :ref:`custom-columns` for details.

Filtering and sorting
---------------------

Now that we have a merged `DataFrame`, we can use `pandas` methods to quickly
get different views of the data. For some analyses, we may want to organize in
terms of the study list by removing repeats and intrusions. Because our data
are in a `DataFrame`, we can use the `DataFrame.query` method:

.. ipython:: python

    merged.query('study')

Alternatively, we may also want to get just the recall events, sorted by
output position instead of input position:

.. ipython:: python

    merged.query('recall').sort_values(['list', 'output'])

Note that we first sort by list, then output position, to keep the
lists together.
User guide
==========

.. toctree::
   :maxdepth: 2

   /guide/import
   /guide/score
   /guide/recall
   /guide/order
   /guide/grouping
===========
Transitions
===========

Psifr has a core set of tools for analyzing transitions in free recall data.
These tools focus on measuring what transitions actually occurred, and which
transitions were possible given the order in which participants recalled items.

Actual and possible transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calculating a conditional response probability involves two parts: the frequency
at which a given event actually occurred in the data and frequency at which a
given event could have occurred. The frequency of possible events is
calculated conditional on the recalls that have been made leading up to each
transition. For example, a transition between item :math:`i` and item :math:`j`
is not considered "possible" in a CRP analysis if item :math:`i` was never
recalled. The transition is also not considered "possible" if, when item
:math:`i` is recalled, item :math:`j` has already been recalled previously.

Repeated recall events are typically excluded from the counts of both actual
and possible transition events. That is, the transition event frequencies are
conditional on the transition not being either to or from a repeated item.

Calculating a CRP measure involves tallying how many transitions of a given
type were made during a free recall test. For example, one common measure is
the serial position lag between items. For a list of length :math:`N`, possible
lags are in the range :math:`[-N+1, N-1]`. Because repeats are excluded, a lag
of zero is never possible. The count of actual and possible transitions for
each lag is calculated first, and then the CRP for each lag is calculated as
the actual count divided by the possible count.

The transitions masker
~~~~~~~~~~~~~~~~~~~~~~

The :py:func:`psifr.transitions.transitions_masker` is a generator that makes
it simple to iterate over transitions while "masking" out events such as
intrusions of items not on the list and repeats of items that have already
been recalled.

On each step of the iterator, the previous, current, and possible items are
yielded. The *previous*
item is the item being transitioned from. The *current* item is the item being
transitioned to. The *possible* items includes an array of all items that
were valid to be recalled next, given the recall sequence up to that point (not
including the current item).

.. ipython::

    In [1]: from psifr.transitions import transitions_masker

    In [2]: pool = [1, 2, 3, 4, 5, 6]

    In [3]: recs = [6, 2, 3, 6, 1, 4]

    In [4]: masker = transitions_masker(pool_items=pool, recall_items=recs,
       ...:                             pool_output=pool, recall_output=recs)

    In [5]: for prev, curr, poss in masker:
       ...:     print(prev, curr, poss)
       ...:

Only valid transitions are yielded, so the code
for a specific analysis only needs to calculate the transition measure of
interest and count the number of actual and possible transitions in each bin
of interest.

Four inputs are required:

`pool_items`
    List of identifiers for all items available for recall. Identifiers
    can be anything that is unique to each item in the list (e.g., serial
    position, a string representation of the item, an index in the stimulus
    pool).

`recall_items`
    List of identifiers for the sequence of recalls, in order. Valid recalls
    must match an item in `pool_items`. Other items are considered intrusions.

`pool_output`
    Output codes for each item in the pool. This should be whatever you need to
    calculate your transition measure.

`recall_output`
    Output codes for each recall in the sequence of recalls.

By using different values for these four inputs and defining different
transition measures, a wide range of analyses can be implemented.
===========
Development
===========

.. toctree::
   :maxdepth: 2

   /development/transitions
Tutorials
=========

See the psifr-notebooks_ project for a set of Jupyter notebooks with
sample code. These examples go more in depth into the options available
for each analysis and how they can be used for advanced analyses such as
conditionalizing CRP analysis on specific transitions.

.. _psifr-notebooks: https://github.com/mortonne/psifr-notebooks
Installation
============

You can install the latest stable version of Psifr using pip:

.. code-block:: bash

    pip install psifr

You can also install the development version directly from the code
repository on GitHub:

.. code-block:: bash

    pip install git+git://github.com/mortonne/psifr
===========
Transitions
===========

.. currentmodule:: psifr.transitions

Counting transitions
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    count_lags
    count_category
    count_distance

Ranking transitions
~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    percentile_rank
    rank_lags
    rank_distance

Iterating over transitions
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    transitions_masker
========
Measures
========

.. currentmodule:: psifr.measures

Transition measure base class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    TransitionMeasure
    TransitionMeasure.split_lists
    TransitionMeasure.analyze
    TransitionMeasure.analyze_subject

Transition measures
~~~~~~~~~~~~~~~~~~~

.. autosummary::

    TransitionOutputs
    TransitionLag
    TransitionLagRank
    TransitionCategory
    TransitionDistance
    TransitionDistanceRank
====================
Free recall analysis
====================

.. currentmodule:: psifr.fr

Managing data
~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    merge_free_recall
    merge_lists
    filter_data
    reset_list
    split_lists

Recall probability
~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    spc
    pnr

Transition probability
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    lag_crp
    category_crp
    distance_crp

Transition rank
~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    lag_rank
    distance_rank

Plotting
~~~~~~~~

.. autosummary::
    :toctree: api/

    plot_raster
    plot_spc
    plot_lag_crp
    plot_distance_crp
    plot_swarm_error
API reference
=============

.. toctree::
    :maxdepth: 2

    /api/fr
    /api/measures
    /api/transitions
    /api/outputs
=======
Outputs
=======

.. currentmodule:: psifr.outputs

Counting recalls by serial position and output position
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    count_outputs

Iterating over output positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
    :toctree: api/

    outputs_masker
