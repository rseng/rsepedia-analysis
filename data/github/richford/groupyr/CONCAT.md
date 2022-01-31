![groupyr logo](https://raw.githubusercontent.com/richford/groupyr/main/doc/_static/groupyr-logo-large.svg)

# _Groupyr_: Sparse Group Lasso in Python

[![Build Status](https://github.com/richford/groupyr/workflows/Build/badge.svg)](https://github.com/richford/groupyr/workflows/Build/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/richford/groupyr/badge.svg?branch=main&service=github)](https://coveralls.io/github/richford/groupyr?branch=main&service=github)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
<br>
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03024/status.svg)](https://doi.org/10.21105/joss.03024)
[![DOI](https://zenodo.org/badge/300933639.svg)](https://zenodo.org/badge/latestdoi/300933639)

_Groupyr_ is a Python library for penalized regression of grouped covariates.
This is the _groupyr_ development site. You can view the source code, file new issues, and contribute to _groupyr_'s development. If you just want to learn how to install and use _groupyr_, please look at the [_groupyr_ documentation][link_groupyr_docs].

## Contributing

We love contributions! _Groupyr_ is open source, built on open source,
and we'd love to have you hang out in our community.

We have developed some [guidelines](.github/CONTRIBUTING.md) for contributing to
_groupyr_.

## Citing _groupyr_

If you use _groupyr_ in a scientific publication, please see cite us:

Richie-Halford et al., (2021). Groupyr: Sparse Group Lasso in Python. Journal of Open Source Software, 6(58), 3024, https://doi.org/10.21105/joss.03024

```
@article{richie-halford-groupyr,
    doi = {10.21105/joss.03024},
    url = {https://doi.org/10.21105/joss.03024},
    year = {2021},
    publisher = {The Open Journal},
    volume = {6},
    number = {58},
    pages = {3024},
    author = {Adam {R}ichie-{H}alford and Manjari Narayan and Noah Simon and Jason Yeatman and Ariel Rokem},
    title = {{G}roupyr: {S}parse {G}roup {L}asso in {P}ython},
    journal = {Journal of Open Source Software}
}
```

## Acknowledgements

_Groupyr_ development is supported through a grant from the [Gordon
and Betty Moore Foundation](https://www.moore.org/) and from the
[Alfred P. Sloan Foundation](https://sloan.org/) to the [University of
Washington eScience Institute](http://escience.washington.edu/), as
well as
[NIMH BRAIN Initiative grant 1RF1MH121868-01](https://projectreporter.nih.gov/project_info_details.cfm?aid=9886761&icde=46874320&ddparam=&ddvalue=&ddsub=&cr=2&csb=default&cs=ASC&pball=)
to Ariel Rokem (University of Washington).

The API design of _groupyr_ was facilitated by the [scikit-learn project
template](https://github.com/scikit-learn-contrib/project-template) and it
therefore borrows heavily from
[scikit-learn](https://scikit-learn.org/stable/index.html). _Groupyr_ relies
on the [copt optimization library](http://openo.pt/copt/index.html) for its
solver. The _groupyr_ logo is a flipped silhouette of an [image from J. E.
Randall](https://commons.wikimedia.org/wiki/File:Epinephelus_amblycephalus,_banded_grouper.jpg)
and is licensed [CC BY-SA](https://creativecommons.org/licenses/by-sa/3.0).

[link_groupyr_docs]: https://richford.github.io/groupyr/
---
title: 'Groupyr: Sparse Group Lasso in Python'
tags:
  - Python
  - group lasso
  - penalized regression
  - classification
authors:
  - name: Adam Richie-Halford
    orcid: 0000-0001-9276-9084
    affiliation: 1
  - name: Manjari Narayan
    orcid: 0000-0001-5348-270X
    affiliation: 2
  - name: Noah Simon
    orcid: 0000-0002-8985-2474
    affiliation: 4
  - name: Jason Yeatman
    orcid: 0000-0002-2686-1293
    affiliation: 5
  - name: Ariel Rokem
    orcid: 0000-0003-0679-1985
    affiliation: 3
affiliations:
  - name: eScience Institute, University of Washington
    index: 1
  - name: Department of Psychiatry and Behavioral Sciences, Stanford University
    index: 2
  - name: Department of Psychology, University of Washington
    index: 3
  - name: Department of Biostatistics, University of Washington
    index: 4
  - name: Graduate School of Education and Division of Developmental and Behavioral Pediatrics, Stanford University
    index: 5
date: 25 Dec 2021
bibliography: paper.bib
---

## Summary

For high-dimensional supervised learning, it is often beneficial to use
domain-specific knowledge to improve the performance of statistical learning
models. When the problem contains covariates which form groups, researchers
can include this grouping information to find parsimonious representations
of the relationship between covariates and targets. These groups may arise
artificially, as from the polynomial expansion of a smaller feature space, or
naturally, as from the anatomical grouping of different brain regions or the
geographical grouping of different cities. When the number of features is
large compared to the number of observations, one seeks a subset of the
features which is sparse at both the group and global level.

The sparse group lasso [@simon2013sparse] is a penalized regression technique
designed for exactly these situations. It combines the original lasso
[@tibshirani1996regression], which induces global sparsity, with the group
lasso [@yuan2006model], which induces group-level sparsity. It estimates a target variable $\hat{y}$ from a
feature matrix $\mathbf{X}$, using

$$
\hat{y} = \mathbf{X} \hat{\beta},
$$

as depicted in \autoref{fig:sgl_model}, with color encoding the group
structure of the covariates in $\mathbf{X}$. The coefficients in
$\hat{\beta}$ characterize the relationship between the features and the
target and must satisfy [@simon2013sparse]

$$
\hat{\beta} = \min_{\beta} \frac{1}{2}
|| y - \sum_{\ell = 1}^{G} \mathbf{X}^{(\ell)} \beta^{(\ell)} ||_2^2
+ (1 - \lambda) \alpha \sum_{\ell = 1}^{G} \sqrt{p_{\ell}} ||\beta^{(\ell)}||_2
+ \lambda \alpha ||\beta||_1,
$$
where $G$ is the total number of groups, $\mathbf{X}^{(\ell)}$ is the
submatrix of $\mathbf{X}$ with columns belonging to group $\ell$,
$\beta^{(\ell)}$ is the coefficient vector of group $\ell$, and $p_{\ell}$ is
the length of $\beta^{(\ell)}$. The model hyperparameter $\lambda$ controls
the combination of the group-lasso and the lasso, with $\lambda=0$ giving the
group lasso fit and $\lambda=1$ yielding the lasso fit. The hyperparameter
$\alpha$ controls the overall strength of the regularization.

![A linear model, $y = \mathbf{X} \cdot \beta$, with grouped covariates. The feature matrix $\mathbf{X}$ is color-coded to reveal a group structure. The coefficients in $\beta$ follow the same grouping. \label{fig:sgl_model}](groupyr_linear_model.pdf)

## Statement of need

*Groupyr* is a Python library that implements the sparse group lasso
as scikit-learn [@sklearn; @sklearn_api] compatible estimators.
It satisfies the need for grouped penalized regression models that
can be used interoperably in researcher's real-world scikit-learn
workflows. Some pre-existing Python libraries come close to satisfying
this need. [*Lightning*](http://contrib.scikit-learn.org/lightning/) [@lightning-2016]
is a Python library for large-scale linear classification and
regression. It supports many solvers with a combination of the
L1 and L2 penalties. However, it does not allow the user to
specify groups of covariates (see, for example, [this GitHub
issue](https://github.com/scikit-learn-contrib/lightning/issues/39)).
Another Python package,
[*group_lasso*](https://group-lasso.readthedocs.io/en/latest/#) [@group-lasso], is a
well-designed and well-documented implementation of the sparse group lasso.
It meets the basic API requirements of scikit-learn compatible estimators.
However, we found that our implementation in *groupyr*, which relies on the
*copt* optimization library [@copt], was between two and ten times faster
for the problem sizes that we encounter in our research (see the
repository's examples directory for a performance comparison).
Additionally, we needed estimators with built-in cross-validation
support using both grid search and sequential model based optimization
strategies. For example, the speed and cross-validation enhancements
were crucial to using *groupyr* in *AFQ-Insight*, a neuroinformatics
research library [@richiehalford2019multidimensional].

## Usage

*Groupyr* is available on the Python Package Index (PyPI) and can be installed
with

```shell
pip install groupyr
```

*Groupyr* is compatible with the scikit-learn API and its estimators offer the
same instantiate, ``fit``, ``predict`` workflow that will be familiar to
scikit-learn users. See the online documentation for a detailed description of the
API and examples in both classification and regression settings. Here, we describe
only the key differences necessary for scikit-learn users to get started with *groupyr*.

For syntactic parallelism with the scikit-learn ``ElasticNet`` estimator, we use the
keyword ``l1_ratio`` to refer to SGL's $\lambda$ hyperparameter. In addition
to keyword parameters shared with scikit-learn's ``ElasticNet``,
``ElasticNetCV``, ``LogisticRegression``, and ``LogisticRegressionCV``
estimators, users must specify the group assignments for the columns of the
feature matrix ``X``. This is done during estimator instantiation using the
``groups`` parameter, which accepts a list of numpy arrays, where the $i$-th
array specifies the feature indices of the $i$-th group. If no grouping
information is provided, the default behavior assigns all features to one
group.

*Groupyr* also offers cross-validation estimators that automatically select
the best values of the hyperparameters $\alpha$ and $\lambda$ using either an
exhaustive grid search (with ``tuning_strategy="grid"``) or sequential model
based optimization (SMBO) using the scikit-optimize library (with
``tuning_strategy="bayes"``). For the grid search strategy, our
implementation is more efficient than using the base estimator with
scikit-learn's ``GridSearchCV`` because it makes use of warm-starting, where
the model is fit along a pre-defined regularization path and the solution
from the previous fit is used as the initial guess for the current
hyperparameter value. The randomness associated with SMBO complicates the use
of a warm start strategy; it can be difficult to determine which of the
previously attempted hyperparameter combinations should provide the initial
guess for the current evaluation. However, even without warm-starting, we
find that the SMBO strategy usually outperforms grid search because far fewer
evaluations are needed to arrive at the optimal hyperparameters. We provide
examples of both strategies (grid search for a classification example and
SMBO for a regression example) in the online documentation.

## Author statements and acknowledgments

The first author (referred to as A.R.H. below) is the lead and corresponding
author. The last author (referred to as A.R.) is the primary supervisor and
is responsible for funding acquisition. All other authors are listed in
alphabetical order by surname. We describe contributions to the paper using
the CRediT taxonomy [@credit].
Writing – Original Draft: A.R.H.;
Writing – Review & Editing: A.R.H., N.S., J.Y., and A.R.;
Conceptualization and methodology: A.R.H., N.S., and A.R.;
Software and data curation: A.R.H., M.N., and A.R.;
Validation: A.R.H. and M.N.;
Resources: A.R.H. and A.R;
Visualization: A.R.H.;
Supervision: N.S., J.Y., and A.R.;
Project Administration: A.R.H;
Funding Acquisition: A.R.;

Groupyr development was supported through a grant from the Gordon and
Betty Moore Foundation and from the Alfred P. Sloan Foundation to the
University of Washington eScience Institute, as well as NIMH BRAIN
Initiative grant 1RF1MH121868-01 to Ariel Rokem at the University of
Washington and through cloud credits from the Google Cloud Platform.

## References
# Contributing to _groupyr_

Welcome to the _groupyr_ repository! We're excited you're here and want to
contribute.

**Imposter's syndrome disclaimer**[^1]: We want your help. No, really.

There may be a little voice inside your head that is telling you that
you're not ready to be an open-source contributor; that your skills
aren't nearly good enough to contribute. What could you possibly offer a
project like this one?

We assure you - the little voice in your head is wrong. If you can
write code at all, you can contribute code to open-source. Contributing
to open-source projects is a fantastic way to advance one's coding
skills. Writing perfect code isn't the measure of a good developer (that
would disqualify all of us!); it's trying to create something, making
mistakes, and learning from those mistakes. That's how we all improve,
and we are happy to help others learn.

Being an open-source contributor doesn't just mean writing code, either.
You can help out by writing documentation, tests, or even giving
feedback about the project (and yes - that includes giving feedback
about the contribution process). Some of these contributions may be the
most valuable to the project as a whole, because you're coming to the
project with fresh eyes, so you can see the errors and assumptions that
seasoned contributors have glossed over.

## Practical guide to submitting your contribution

These guidelines are designed to make it as easy as possible to get involved.
If you have any questions that aren't discussed below, please let us know by
opening an [issue][link_issues]!

Before you start, you'll need to set up a free [GitHub][link_github] account
and sign in. Here are some [instructions][link_signupinstructions].

Already know what you're looking for in this guide? Jump to the following sections:

- [Joining the conversation](#joining-the-conversation)
- [Contributing through Github](#contributing-through-github)
- [Understanding issues](#understanding-issues)
- [Making a change](#making-a-change)
- [Structuring contributions](#groupyr-coding-style-guide)
- [Licensing](#licensing)
- [Recognizing contributors](#recognizing-contributions)

## Joining the conversation

_Groupyr_ is primarily maintained by a [collaborative research group][link_autofq].
But we maintain this software as an open project. This means that we welcome
contributions from people outside our group and we make sure to give
contributors from outside our group credit in presentations of the work. In
other words, we're excited to have you join! Most of our discussions will
take place on open [issues][link_issues]. We actively monitor this space and
look forward to hearing from you!

## Contributing through GitHub

[git][link_git] is a really useful tool for version control.
[GitHub][link_github] sits on top of git and supports collaborative and distributed working.

If you're not yet familiar with `git`, there are lots of great resources to
help you _git_ started!
Some of our favorites include the [git Handbook][link_handbook] and
the [Software Carpentry introduction to git][link_swc_intro].

On GitHub, You'll use [Markdown][link_markdown] to chat in issues and pull
requests. You can think of Markdown as a few little symbols around your text
that will allow GitHub to render the text with a little bit of formatting.
For example, you could write words as bold (`**bold**`), or in italics
(`*italics*`), or as a [link][link_rick_roll]
(`[link](https://youtu.be/dQw4w9WgXcQ)`) to another webpage.

GitHub has a really helpful page for getting started with [writing and
formatting Markdown on GitHub][link_writing_formatting_github].

## Understanding issues

Every project on GitHub uses [issues][link_issues] slightly differently.

The following outlines how the _groupyr_ developers think about these tools.

- **Issues** are individual pieces of work that need to be completed to move the project forward.
  A general guideline: if you find yourself tempted to write a great big issue that
  is difficult to be described as one unit of work, please consider splitting it into two or more issues.

      Issues are assigned [labels](#issue-labels) which explain how they relate to the overall project's
      goals and immediate next steps.

### Issue Labels

The current list of issue labels are [here][link_labels] and include:

- [![Good first issue](https://img.shields.io/github/labels/richford/groupyr/good%20first%20issue)][link_firstissue] _These issues contain a task that is amenable to new contributors because it doesn't entail a steep learning curve._

  If you feel that you can contribute to one of these issues,
  we especially encourage you to do so!

- [![Bug](https://img.shields.io/github/labels/richford/groupyr/bug)][link_bugs] _These issues point to problems in the project._

  If you find new a bug, please give as much detail as possible in your issue,
  including steps to recreate the error.
  If you experience the same bug as one already listed,
  please add any additional information that you have as a comment.

- [![Enhancement](https://img.shields.io/github/labels/richford/groupyr/enhancement)][link_enhancement] _These issues are asking for new features and improvements to be considered by the project._

  Please try to make sure that your requested feature is distinct from any others
  that have already been requested or implemented.
  If you find one that's similar but there are subtle differences,
  please reference the other request in your issue.

In order to define priorities and directions in the development roadmap,
we have two sets of special labels:

| Label                                                                                                                                                                                                                                                                                 | Description                                                                               |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20high) <br> ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20medium) <br> ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20low) | Estimation of the downstream impact the proposed feature/bugfix will have.                |
| ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20high) <br> ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20medium) <br> ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20low) | Estimation of effort required to implement the requested feature or fix the reported bug. |

These labels help triage and set priorities to the development tasks.
For instance, one bug regression that has been reported to affect most of the users after
a release with an easy fix because it is a known old problem that came back.
Such an issue will typically be assigned the following labels ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/bug) ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20high) ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20low), and its priority will be maximal since addressing low-effort high-impact issues delivers the maximum turnout without increasing the churn by much.

Of course, the implementation of long-term goals may include the scheduling of ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20medium) ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20high).
Finally, ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/impact%3A%20low) ![GitHub labels](https://img.shields.io/github/labels/richford/groupyr/effort%3A%20high) issues are less likely to be addressed.

## Making a change

We appreciate all contributions to _groupyr_,
but those accepted fastest will follow a workflow similar to the following:

1. **Comment on an existing issue or open a new issue referencing your addition.**<br />
   This allows other members of the _groupyr_ development team to confirm that
   you aren't overlapping with work that's currently underway and that
   everyone is on the same page with the goal of the work you're going to
   carry out.<br /> [This blog][link_pushpullblog] is a nice explanation of
   why putting this work in up front is so useful to everyone involved.

1. **[Fork][link_fork] the [groupyr repository][link_groupyr] to your profile.**<br />
   This is now your own unique copy of _groupyr_.
   Changes here won't effect anyone else's work, so it's a safe space to
   explore edits to the code!
   On your own fork of the repository, select Settings -> Actions-> "Disable
   Actions for this repository" to avoid flooding your inbox with warnings
   from our continuous integration suite.

1. **[Clone][link_clone] your forked _groupyr_ repository to your machine/computer.**<br />
   While you can edit files [directly on github][link_githubedit], sometimes
   the changes you want to make will be complex and you will want to use a
   [text editor][link_texteditor] that you have installed on your local
   machine/computer. (One great text editor is [vscode][link_vscode]).<br />
   In order to work on the code locally, you must clone your forked repository.<br />
   To keep up with changes in the _groupyr_ repository,
   add the ["upstream" _groupyr_ repository as a remote][link_addremote]
   to your locally cloned repository.

   ```Shell
   git remote add upstream https://github.com/richford/groupyr.git
   ```

   Make sure to [keep your fork up to date][link_updateupstreamwiki] with the upstream repository.<br />
   For example, to update your main branch on your local cloned repository:

   ```Shell
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

1. **Install a development version of _groupyr_ so that your local changes are reflected in your local tests**<br />
   You can install a development version of _groupyr_ by navigating to the root of your _groupyr_ repository and then typing

   ```Shell
   pip install -e .[dev]
   ```

   [or][link_pythonmpip_1] [better][link_pythonmpip_2] [still][link_pythonmpip_3]

   ```Shell
   python -m pip install -e .[dev]
   ```

1. **Create a [new branch][link_branches] to develop and maintain the proposed code changes.**<br />
   For example:

   ```Shell
   git fetch upstream  # Always start with an updated upstream
   git checkout -b fix/bug-1222 upstream/main
   ```

   Please consider using appropriate branch names as those listed below:

   | Branch name             | Use case                       |
   | ----------------------- | ------------------------------ |
   | `fix/<some-identifier>` | for bugfixes                   |
   | `enh/<feature-name>`    | for new features               |
   | `doc/<some-identifier>` | for documentation improvements |

   You should name all your documentation branches with the prefix `doc/`
   as that will preempt triggering the full battery of continuous integration tests.

1. **Make the changes you've discussed, following the [_groupyr_ coding style guide](#groupyr-coding-style-guide).**<br />
   Try to keep the changes focused: it is generally easy to review changes that address one feature or bug at a time.
   Assuming you installed a development environment above, you can test your local changes using

   ```Shell
   make lint
   make test
   ```

   Once you are satisfied with your local changes, [add/commit/push them][link_add_commit_push]
   to the branch on your forked repository.

1. **Submit a [pull request][link_pullrequest].**<br />
   A member of the development team will review your changes to confirm
   that they can be merged into the main code base.<br />
   Pull request titles should begin with a descriptive prefix
   (for example, `ENH: Adding another estimator class`):

   - `ENH`: enhancements or new features ([example][ex_enh])
   - `FIX`: bug fixes ([example][ex_fix])
   - `TST`: new or updated tests ([example][ex_tst])
   - `DOC`: new or updated documentation ([example][ex_doc])
   - `STY`: style changes ([example][ex_sty])
   - `REF`: refactoring existing code ([example][ex_ref])
   - `CI`: updates to continous integration infrastructure ([example][ex_ci])
   - `MAINT`: general maintenance ([example][ex_maint])
   - For works-in-progress, add the `WIP` tag in addition to the descriptive prefix.
     Pull-requests tagged with `WIP:` will not be merged until the tag is removed.

1. **Have your PR reviewed by the development team, and update your changes accordingly in your branch.**<br />
   The reviewers will take special care in assisting you to address their
   comments, as well as dealing with conflicts and other tricky situations
   that could emerge from distributed development. And if you don't make the
   requested changes, we might ask [@bedevere-bot][link_bedevere] to [poke
   you with soft cushions][link_bedevere_video]!

## _Groupyr_ coding style guide

We use the [Black code formatter][link_black] for format our code
contributions to a common style. All pull requests will automatically be
checked for compliance using `flake8` and the Black code formatter.
We recommend you activate the pre-commit formatting hook by typing

```shell
pre-commit install
```

Afterward, all of your local changes will be automatically formatted before
you commit them. This is the easiest way to ensure that your much appreciated contribution is not delayed due to formatting or style compliance.

### Documentation

We use [Sphinx][link_sphinx] to generate documentation from files stored in the
`docs/source` folder. To generate proper documentation of functions, we use the
[numpy docstring standard][link_np_docstring] when documenting code inline in
docstrings.

## Licensing

_Groupyr_ is licensed under the BSD license. By contributing to _groupyr_, you
acknowledge that any contributions will be licensed under the same terms.

### Reminder note for maintainers

_Groupyr_ pushes a development version to
[Test-PyPI](https://test.pypi.org/) on every pull request merged into
the main branch. To release a new version of _groupyr_, use the `publish_release.sh` script from the root directory, i.e.:

```Shell
.maintenance/publish_release.sh <version_number>
```

For releases, use the following format for <version*number>:
"v<major>.<minor>.<micro>".
When executed, this will ask you if you want to customize the
`CHANGES.rst` document or the release notes. After that, *groupyr\*'s
GitHub actions will take care of publishing the new release on PyPI and
creating a release on GitHub.

[^1]:
    The imposter syndrome disclaimer was originally written by
    [Adrienne Lowe](https://github.com/adriennefriend) for a
    [PyCon talk](https://www.youtube.com/watch?v=6Uj746j9Heo), and was
    adapted based on its use in the README file for the
    [MetPy project](https://github.com/Unidata/MetPy).

[ex_ci]: https://github.com/richford/groupyr/pull/8
[ex_doc]: https://github.com/richford/groupyr/pull/10
[ex_enh]: https://github.com/richford/groupyr/pull/23
[ex_fix]: https://github.com/richford/groupyr/pull/16
[ex_maint]: https://github.com/richford/groupyr/pull/17
[ex_sty]: https://github.com/richford/groupyr/pull/21
[ex_tst]: https://github.com/richford/groupyr/pull/11
[link_add_commit_push]: https://help.github.com/articles/adding-a-file-to-a-repository-using-the-command-line
[link_addremote]: https://help.github.com/articles/configuring-a-remote-for-a-fork
[link_autofq]: https://autofq.org/
[link_bedevere]: https://github.com/search?q=commenter%3Abedevere-bot+soft+cushions
[link_bedevere_video]: https://youtu.be/XnS49c9KZw8?t=1m7s
[link_black]: https://black.readthedocs.io/en/stable/
[link_branches]: https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/
[link_bugs]: https://github.com/richford/groupyr/labels/bug
[link_clone]: https://help.github.com/articles/cloning-a-repository
[link_discussingissues]: https://help.github.com/articles/discussing-projects-in-issues-and-pull-requests
[link_enhancement]: https://github.com/richford/groupyr/labels/enhancement
[link_firstissue]: https://github.com/richford/groupyr/labels/good%20first%20issue
[link_fork]: https://help.github.com/articles/fork-a-repo/
[link_git]: https://git-scm.com/
[link_github]: https://github.com/
[link_githubedit]: https://help.github.com/articles/editing-files-in-your-repository
[link_groupyr]: https://github.com/richford/groupyr
[link_handbook]: https://guides.github.com/introduction/git-handbook/
[link_issues]: https://github.com/richford/groupyr/issues
[link_labels]: https://github.com/richford/groupyr/labels
[link_markdown]: https://daringfireball.net/projects/markdown
[link_np_docstring]: https://numpydoc.readthedocs.io/en/latest/format.html
[link_pullrequest]: https://help.github.com/articles/creating-a-pull-request-from-a-fork
[link_pushpullblog]: https://www.igvita.com/2011/12/19/dont-push-your-pull-requests/
[link_pythonmpip_1]: https://adamj.eu/tech/2020/02/25/use-python-m-pip-everywhere/
[link_pythonmpip_2]: https://snarky.ca/why-you-should-use-python-m-pip/
[link_pythonmpip_3]: https://github.com/pypa/pip/issues/3164#issue-109993120
[link_rick_roll]: https://www.youtube.com/watch?v=dQw4w9WgXcQ
[link_signupinstructions]: https://help.github.com/articles/signing-up-for-a-new-github-account
[link_sphinx]: http://www.sphinx-doc.org/en/master/
[link_swc_intro]: http://swcarpentry.github.io/git-novice/
[link_texteditor]: https://en.wikipedia.org/wiki/Text_editor
[link_updateupstreamwiki]: https://help.github.com/articles/syncing-a-fork/
[link_vscode]: https://code.visualstudio.com/
[link_writing_formatting_github]: https://help.github.com/articles/getting-started-with-writing-and-formatting-on-github
---
name: Bug report
about: Create a report to help us reproduce and correct the bug
title: ''
labels: 'bug'
assignees: ''
---

<!--
Before submitting a bug, please make sure the issue hasn't been already
addressed by searching through the past issues.
-->

#### Describe the bug

<!--
A clear and concise description of what the bug is.
-->

#### Steps/Code to Reproduce

<!--
Please add a minimal example to reproduce the error by running the code.
Please, be as succinct as possible and do not depend on external data. In
short, we would like to be able to copy-paste your code and get the same
result as you.

If the code is too long, feel free to put it in a public gist and link it in
the issue: https://gist.github.com
-->

    Sample code to reproduce the problem

#### Expected Results

<!-- Example: No error is thrown. Please paste or describe the expected results.-->

#### Actual Results

<!-- Please paste or specifically describe the actual output or traceback. -->

#### Versions

<!--
Please run the following snippet and paste the output below.

import groupyr as gpr
print(gpr.__version__)
-->

<!-- Thanks for contributing! -->
---
name: Other
about: For all other issues. 
title: ''
labels: ''
assignees: ''

---
---
name: Documentation improvement
about: Create a report to help us improve the documentation. Alternatively you can just open a pull request with the suggested change.
title: ''
labels: 'documentation'
assignees: ''
---

#### Describe the issue linked to the documentation

<!--
Tell us about the confusion introduced in the documentation.
-->

#### Suggest a potential alternative/fix

<!--
Tell us how we could improve the documentation in this regard.
-->
---
name: Feature request
about: Suggest a new algorithm, enhancement to an existing algorithm, etc.
title: ''
labels: 'enhancement'
assignees: ''

---

<!--
Thanks for suggesting an improvement to groupyr.
-->

#### Describe the workflow you want to enable

#### Describe your proposed solution

#### Describe alternatives you've considered, if relevant

#### Additional context
v0.2.5 (July 27, 2021)
======================
  * ENH: Add GroupResampler (#62)
  * ENH: Add select_intersection kwarg to transformers (#61)
  * ENH: Add GroupAggregator, tests, and doc API (#59)

v0.2.4 (June 22, 2021)
======================
  * Add sgl_path example to the documentation (#58)
  * Add GroupPCA, and supervised PCA variants (#55)

v0.2.3 (March 17, 2021)
=======================
  * ENH: Add GroupFPCA (#48)
  * ENH: Add group transformers (#51)
  * DOC: Update README.md with JOSS article info (#53)
  * DOC: One typo (#52)
  * CI: Only publish docs to GitHub pages for one Python version (#49)

v0.2.2 (February 24, 2021)
==========================
  * Micro release to accompany publication of JOSS paper

v0.2.1 (February 21, 2021)
==========================
  * DEP: Loosen dependency requirements (#46)
  * DOC: Add groupyr/group-lasso comparison example (#44)
  * Move matplotlib dependency to dev option (#45)
  * DOC: Add author contributions using the CRediT taxonomy (#42)
  * DOC: Add help target to makefile to make is self-documenting (#43)
  * ENH: Add lightning and group-lasso citations to paper (#41)
  * DEP: remove ipywidgets from setup.cfg dependencies (#40)
  * Remove redundant paper reference (#38)

v0.1.10 (December 10, 2020)
===========================
  * FIX: Assign error_score in BayesSearchCV (#34)

v0.1.9 (December 09, 2020)
==========================
  * ENH: Use sgl_scoring_path instead of sklearn's _path_residuals (#33)
  * FIX: Sets bayes_optimizer_ to None when "grid" strategy is used (#32)

v0.1.8 (December 05, 2020)
==========================
  * ENH: Add BayesSearchCV option to SGLCV and LogisticSGLCV (#31)

v0.1.7 (October 26, 2020)
=========================
  * ENH: Use joblib Parallel instead of custom _ProgressParallel wrapper (#28)


v0.1.6 (October 22, 2020)
=========================
  * DOC: Fix mathjax rendering in documentation (#27)


v0.1.5 (October 15, 2020)
=========================
  * ENH: Speed up `SparseGroupL1.prox()` (#23)


v0.1.4 (October 09, 2020)
=========================
  * DOC: Update citation instructions (#22)
  * STY: Prefer `python -m pip` over `pip`. Also use taxicab random_state instead of 42 (#21)
  * FIX: Use classifier=True in check_cv for LogisticSGLCV (#20)
  * MAINT: Automatically update zenodo file as part of release script (#17)
  * FIX: Shuffle groups make_group_regression. Use `generator.choice` in make_group_classification (#16)
  * TST: Add tests for `SGL`, `SGLCV`, etc. (#13)


v0.1.3 (October 05, 2020)
=========================
  * TST: Test sparsity masks and chosen groups/features in _base.py (#12)
  * TST: Test check_groups in utils.py (#11)
  * DOC: Add pull request examples to CONTRIBUTING.md (#10)
  * CI: Use pydocstyle in github actions (#8)
  * DOC: Replaces CRCNS grant with our BRAINI grant in README (#7)


v0.1.2 (October 03, 2020)
=========================

- Bump version to confirm GitHub action behavior.


v0.1.1 (October 03, 2020)
=========================

- Fix automatic documentation building.


v0.1.0 (October 03, 2020)
=========================

- Initial release

v0.2.5 (July 27, 2021)
======================
  * ENH: Add GroupResampler (#62)
  * ENH: Add select_intersection kwarg to transformers (#61)
  * ENH: Add GroupAggregator, tests, and doc API (#59)

#########################
Contributing to *groupyr*
#########################

*Groupyr* is an open-source software project. This means that you are welcome
to use it (see `LICENSE
<https://github.com/richford/groupyr/blob/main/LICENSE>`_ for details).

It also means that we welcome contributions from developers outside of our
research collaboration. We love contributions and we will credit them
appropriately!

We've written a more detailed contribution guide `here
<https://github.com/richford/groupyr/blob/main/.github/CONTRIBUTING.md>`_.
############################
Getting help using *groupyr*
############################

If you run into any bugs or issues using *groupyr*, please let us know by
`posting an issue <https://github.com/richford/groupyr/issues/new/choose>`_.
on our GitHub repository. You can browse existing issues `here
<https://github.com/richford/groupyr/issues>`_. Please provide all pertinent
information about the issue you are facing: we often need to know what
version of the software you are running and how you installed it, as well as
the operating system that you are using. If you run into issues with your
data, we might ask you to share a small sample of your data with us, so that
we can examine it. Please make sure that your study ethics procedure allows
you to share this data before sending it on to us.
.. _faq-label:

Frequently Asked Questions
==========================

Here we'll maintain a list of frequently asked questions about *groupyr*. Do
you have a question that isn't addressed here? If so, please see our `getting
help page <getting_help.html>`_ for information about how to file a new
issue.

.. dropdown:: Why did we create *groupyr* and how does it compare to other similar packages?

    We created *groupyr* to be useful in our own research and we hope it is
    useful in yours. There are other packages for penalized regression in
    python. The `lightning <http://contrib.scikit-learn.org/lightning/>`_
    package has a lasso penalty but does not allow the user to specify
    groups. The `group_lasso
    <https://group-lasso.readthedocs.io/en/latest/#>`_ package is well
    designed and documented, but we found that *groupyr*'s execution time was
    faster for most problems. We also wanted estimators with built-in
    cross-validation using both grid search and the ``BayesSearchCV``
    sequential model based optimization.

    In the future, we hope that *groupyr* will include other methods for
    statistical learning with grouped covariates (e.g. unsupervised learning
    methods). These would also be out of scope for the aforementioned
    libraries. However, we encourage you to try many tools. If you find that
    another one is better suited to your problem, please `leave us some
    feedback <https://github.com/richford/groupyr/issues/new/choose>`_, go
    forth, and do good work. Happy coding!
*Groupyr*: Sparse Group Lasso in Python
=======================================

*Groupyr* is a scikit-learn compatible implementation of the sparse group lasso
linear model. It is intended for high-dimensional supervised learning
problems where related covariates can be assigned to predefined groups.

The Sparse Group Lasso
----------------------

The sparse group lasso [1]_ is a penalized regression approach that combines the
group lasso with the normal lasso penalty to promote both global sparsity and
group-wise sparsity. It estimates a target variable :math:`\hat{y}` from a
feature matrix :math:`\mathbf{X}`, using

.. math::

    \hat{y} = \mathbf{X} \hat{\beta},

where the coefficients in :math:`\hat{\beta}` characterize the relationship
between the features and the target and must satisfy [1]_

.. math::

    \hat{\beta} = \min_{\beta} \frac{1}{2}
    || y - \sum_{\ell = 1}^{G} \mathbf{X}^{(\ell)} \beta^{(\ell)} ||_2^2
    + (1 - \alpha) \lambda \sum_{\ell = 1}^{G} \sqrt{p_{\ell}} ||\beta^{(\ell)}||_2
    + \alpha \lambda ||\beta||_1,
   
where :math:`G` is the total number of groups, :math:`\mathbf{X}^{(\ell)}` is
the submatrix of :math:`\mathbf{X}` with columns belonging to group
:math:`\ell`, :math:`\beta^{(\ell)}` is the coefficient vector of group
:math:`\ell`, and :math:`p_{\ell}` is the length of :math:`\beta^{(\ell)}`.
The model hyperparameter :math:`\alpha` controls the combination of the
group-lasso and the lasso, with :math:`\alpha=0` giving the group lasso fit
and :math:`\alpha=1` yielding the lasso fit. The hyperparameter
:math:`\lambda` controls the strength of the regularization.

.. toctree::
   :hidden:
   :titlesonly:

   Home <self>


.. toctree::
   :maxdepth: 3
   :hidden:

   install
   auto_examples/index
   getting_help
   api
   FAQ <faq>
   contributing
   Groupyr on GitHub <https://github.com/richford/groupyr>

`Installation <install.html>`_
------------------------------

See the `installation guide <install.html>`_ for installation instructions.

Usage
-----

*Groupyr* is compatible with the scikit-learn API and its estimators offer the
same instantiate, ``fit``, ``predict`` workflow that will be familiar to
scikit-learn users. See the `API <api.html>`_ and `examples
<auto_examples/index.html>`_ for full details. Here, we describe only the key
differences necessary for scikit-learn users to get started with *groupyr*.

For syntactic parallelism with the scikit-learn ``ElasticNet`` estimator, we
use the keyword ``l1_ratio`` to refer to SGL's :math:`\alpha` hyperparameter
above that controls the mixture of group lasso and lasso penalties. In
addition to keyword parameters shared with scikit-learn's ``ElasticNet``,
``ElasticNetCV``, ``LogisticRegression``, and ``LogisticRegressionCV``
estimators, users must specify the group assignments for the columns of the
feature matrix ``X``. This is done during estimator instantiation using the
``groups`` parameter, which accepts a list of numpy arrays, where the
:math:`i`-th array specifies the feature indices of the :math:`i`-th group.
If no grouping information is provided, the default behavior assigns all
features to one group.

*Groupyr* also offers cross-validation estimators that automatically select
the best values of the hyperparameters :math:`\alpha` and :math:`\lambda`
using either an exhaustive grid search (with ``tuning_strategy="grid"``) or
sequential model based optimization (SMBO) using the scikit-optimize library
(with ``tuning_strategy="bayes"``). For the grid search strategy, our
implementation is more efficient than using the base estimator with
scikit-learn's ``GridSearchCV`` because it makes use of warm-starting, where
the model is fit along a pre-defined regularization path and the solution
from the previous fit is used as the initial guess for the current
hyperparameter value. The randomness associated with SMBO complicates the use
of a warm start strategy; it can be difficult to determine which of the
previously attempted hyperparameter combinations should provide the initial
guess for the current evaluation. However, even without warm-starting, we
find that the SMBO strategy usually outperforms grid search because far fewer
evaluations are needed to arrive at the optimal hyperparameters. We provide
`examples <auto_examples/index.html>`_ of both strategies.

`API Documentation <api.html>`_
-------------------------------

See the `API Documentation <api.html>`_ for detailed documentation of the API.

`Examples <auto_examples/index.html>`_
--------------------------------------

And look at the `example gallery <auto_examples/index.html>`_ for a set of introductory examples.

Citing groupyr
--------------

If you use *groupyr* in a scientific publication, we would appreciate
citations. Please see our `citation instructions
<https://github.com/richford/groupyr#citing-groupyr>`_ for the latest
reference and a bibtex entry.

Acknowledgements
----------------

*Groupyr* development is supported through a grant from the `Gordon and Betty
Moore Foundation <https://www.moore.org/>`_ and from the `Alfred P. Sloan
Foundation <https://sloan.org/>`_ to the `University of Washington eScience
Institute <http://escience.washington.edu/>`_, as well as `NIMH BRAIN
Initiative grant 1RF1MH121868-01
<https://projectreporter.nih.gov/project_info_details.cfm?aid=9886761&icde=46874320&ddparam=&ddvalue=&ddsub=&cr=2&csb=default&cs=ASC&pball=)>`_
to Ariel Rokem (University of Washington).

The API design of *groupyr* was facilitated by the `scikit-learn project
template`_ and it therefore borrows heavily from `scikit-learn`_ [2]_.
*Groupyr* relies on the copt optimization library [3]_ for its solver. The
*groupyr* logo is a flipped silhouette of an `image from J. E. Randall`_ and is
licensed `CC BY-SA`_.

.. _scikit-learn project template: https://github.com/scikit-learn-contrib/project-template
.. _scikit-learn: https://scikit-learn.org/stable/index.html
.. _image from J. E. Randall: https://commons.wikimedia.org/wiki/File:Epinephelus_amblycephalus,_banded_grouper.jpg
.. _CC BY-SA: https://creativecommons.org/licenses/by-sa/3.0

References
----------
.. [1] Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013).
    A sparse-group lasso. Journal of Computational and Graphical
    Statistics, 22(2), 231-245.
.. [2] Pedregosa et al. (2011). `Scikit-learn: Machine Learning in Python`_.
    Journal of Machine Learning Research, 12, 2825-2830;
    Buitnick et al. (2013). `API design for machine learning software:
    experiences from the scikit-learn project`_. ECML PKDD Workshop: Languages
    for Data Mining and Machine Learning, 108-122.
.. [3] Pedregosa et al. (2020). `copt: composite optimization in Python`__.
    DOI:10.5281/zenodo.1283339.
    
.. _Scikit-learn\: Machine Learning in Python: http://jmlr.csail.mit.edu/papers/v12/pedregosa11a.html
.. _API design for machine learning software\: experiences from the scikit-learn project: https://arxiv.org/abs/1309.0238
.. __: http://openopt.github.io/copt/############
Installation
############

*Groupyr* requires Python 3.6, 3.7, or 3.8 and depends on

    copt
    numpy
    scikit-learn
    scipy
    scikit-optimize
    tqdm

Installing the release version
------------------------------

The recommended way to install *groupyr* is from PyPI,

.. code-block:: console

    $ pip install groupyr

This will install *groupyr* and all of its dependencies.

Installing the development version
----------------------------------

The development version is less stable but may include new features.
You can install the development version using ``pip``:

.. code-block:: console

    pip install git+https://github.com/richford/groupyr.git 

Alternatively, you can clone the source code from the `github repository
<https://github.com/richford/groupyr>`_:

.. code-block:: console

    $ git clone git@github.com:richford/groupyr.git
    $ cd groupyr
    $ pip install .

If you would like to contribute to *groupyr*, see the `contributing guidelines
<contributing.html>`_.

Next, go to the `user guide <user_guide.html>`_ or see the `example gallery
<auto_examples/index.html>`_ for further information on how to use *groupyr*.

Dependencies
------------

Installing *groupyr* using either of the methods above will install all of
its dependencies: copt, numpy, scikit-learn, scipy, scikit-optimize, and
tqdm.#############
API Reference
#############

*Groupyr* contains estimator classes that are fully compliant
with the `scikit-learn <https://scikit-learn.org>`_ ecosystem. Consequently,
their initialization, ``fit``, ``predict``, ``transform``, and ``score``
methods will be familiar to ``sklearn`` users.

.. currentmodule:: groupyr

Sparse Groups Lasso Estimators
==============================

These are *groupyr*'s canonical estimators. ``SGL`` is intended for regression
problems while ``LogisticSGL`` is intended for classification problems.

.. autoclass:: SGL

.. autoclass:: LogisticSGL

Cross-validation Estimators
===========================

These estimators have built-in `cross-validation
<https://scikit-learn.org/stable/modules/cross_validation.html#cross-validation-evaluating-estimator-performance>`_
capabilities to find the best values of the hyperparameters ``alpha`` and
``l1_ratio``. These are more efficient than using the canonical estimators
with grid search because they make use of warm-starting. Alternatively, you
can specify ``tuning_strategy = "bayes"`` to use `Bayesian optimization over
the hyperparameters
<https://scikit-optimize.github.io/stable/modules/generated/skopt.BayesSearchCV.html>`_
instead of a grid search.

.. autoclass:: SGLCV

.. autoclass:: LogisticSGLCV

Dataset Generation
==================

Use these functions to generate synthetic sparse grouped data.

.. currentmodule:: groupyr.datasets

.. autofunction:: make_group_classification

.. autofunction:: make_group_regression

Regularization Paths
====================

Use these functions to compute regression coefficients along a regularization path.

.. currentmodule:: groupyr

.. autofunction:: sgl_path

.. currentmodule:: groupyr.logistic

.. autofunction:: logistic_sgl_path

Group Transformers
==================

These classes perform group-wise transformations on their inputs.

.. currentmodule:: groupyr.transform

.. autoclass:: GroupExtractor

.. autoclass:: GroupRemover

.. autoclass:: GroupShuffler

.. autoclass:: GroupAggregator

.. autoclass:: GroupResampler:mod:`{{module}}`.{{objname}}
{{ underline }}====================

.. currentmodule:: {{ module }}

.. autofunction:: {{ objname }}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div style='clear:both'></div>
:mod:`{{module}}`.{{objname}}
{{ underline }}==============

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   .. automethod:: __init__
   {% endblock %}

.. include:: {{module}}.{{objname}}.examples

.. raw:: html

    <div style='clear:both'></div>
