# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as contributors and maintainers pledge to making participation in our project and our community a harassment-free experience for everyone, regardless of age, body size, disability, ethnicity, gender identity and expression, level of experience, nationality, personal appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable behavior and are expected to take appropriate and fair corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, or to ban temporarily or permanently any contributor for other behaviors that they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces when an individual is representing the project or its community. Examples of representing a project or community include using an official project e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event. Representation of a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at lukas.heinrich@cern.ch. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Contributing to pyhf

We are happy to accept contributions to `pyhf` via Pull Requests to the GitHub repository and welcome Issues.
To get started fork the repo.

## Issues

Making Issues is very helpful to the project &mdash; they help the dev team form the development roadmap and are where most important discussion takes place.
If you have suggestions, questions that you can't find answers to on the [documentation website](https://scikit-hep.org/pyhf/) or on the [GitHub Discussions](https://github.com/scikit-hep/pyhf/discussions), or have found a bug please [open an Issue](https://github.com/scikit-hep/pyhf/issues/new/choose)!

## Pull Requests

## Opening an Issue to Discuss

Unless your Pull Request is an obvious 1 line fix, please first [open an Issue](https://github.com/scikit-hep/pyhf/issues/new/choose) to discuss your PR with the dev team.
The Issue allows for discussion on the usefulness and scope of the PR to be publicly discussed and also allows for the PR to then be focused on the code review.
The `pyhf` dev team wants to encourage contributions and community involvement in the project and also avoid low quality PRs that don't actually fix an issue or contribute towards the current roadmap.
PRs that don't follow these guidelines might be rejected.

### Good Examples

If you're looking for some examples of high quality contributed pull requests we recommend you take a look at these:

- PR [#902](https://github.com/scikit-hep/pyhf/pull/902) by Nikolai Hartmann ([@nikoladze](https://github.com/nikoladze))

Many thanks goes out to our contributors!

### Drafts

Unless you are making a single commit pull request please create a draft pull request. Outline the work that will be done in this ongoing pull request. When you are close to being done please tag someone with Approver permissions to follow the pull request.

### Pull Request Procedure

If you would like to make a pull request please:

1. Make a fork of the project.
2. Open an Issue to discuss the planned PR with the project maintainers.
3. Commit your changes to a feature branch on your fork and push to your branch.
4. Start a pull request to let the project maintainers know you're working on it.
5. Test your changes with `pytest`.
6. Update your fork to make sure your changes don't conflict with the current state of the master branch.
7. Make sure that you've added your name to `docs/contributors.rst`.
If you haven't **please** do so by simply appending your name to the bottom of the list.
We are thankful for and value your contributions to `pyhf`, not matter the size.
8. Request your PR be reviewed by the project maintainers.

## Bug Reports

If you have found a bug please report it by filling out the [bug report template](https://github.com/scikit-hep/pyhf/issues/new?template=Bug-Report.md&labels=bug&title=Bug+Report+:+Title+Here).

## Installing the development environment

We recommend first reading the "[Developing](https://scikit-hep.org/pyhf/development.html)" page on the pyhf website and the coming back here.

You can install the development environment (which includes a number of extra) libraries and all others needed to run the tests via `pip`:

```console
python -m pip install --upgrade --editable .[complete]
```

To make the PR process much smoother we also strongly recommend that you setup the Git pre-commit hook for [Black](https://github.com/psf/black) by running

```console
pre-commit install
```

This will run `black` over your code each time you attempt to make a commit and warn you if there is an error, canceling the commit.

## Running the tests

You can run the unit tests (which should be fast!) via the following command.

```console
pytest --ignore=tests/test_notebooks.py
```

Note: This ignores the notebook tests (which are run via [papermill](https://github.com/nteract/papermill) which run somewhat slow.
Make sure to run the complete suite before submitting a PR

```console
pytest
```

## Making a pull request

We try to follow [Conventional Commit](https://www.conventionalcommits.org/) for commit messages and PR titles. Since we merge PR's using squash commits, it's fine if the final commit messages (proposed in the PR body) follow this convention.

## Generating Reference Visuals

New baseline visuals can be generated using this command:

```console
pytest tests/contrib/test_viz.py --mpl-generate-path=tests/contrib/baseline
```
# Description

Please first read [CONTRIBUTING.md](https://github.com/scikit-hep/pyhf/tree/master/CONTRIBUTING.md).

Please describe the purpose of this pull request in some detail. Reference and link to any relevant issues or pull requests.

# Checklist Before Requesting Reviewer

- [ ] Tests are passing
- [ ] "WIP" removed from the title of the pull request
- [ ] Selected an Assignee for the PR to be responsible for the log summary

# Before Merging

For the PR Assignees:

- [ ] Summarize commit messages into a comprehensive review of the PR
# Pull Request Description

Please first read [CONTRIBUTING.md](https://github.com/scikit-hep/pyhf/tree/master/CONTRIBUTING.md).

Please describe the purpose of this pull request in some detail and what the specific feature being added will do. Reference and link to any relevant issues or pull requests (such as the issue in which this feature was first suggested).

# Checklist Before Requesting Reviewer

- [ ] Tests are passing
- [ ] "WIP" removed from the title of the pull request
- [ ] Selected an Assignee for the PR to be responsible for the log summary

# Before Merging

For the PR Assignees:

- [ ] Summarize commit messages into a comprehensive review of the PR
# Pull Request Description

Please first read [CONTRIBUTING.md](https://github.com/scikit-hep/pyhf/tree/master/CONTRIBUTING.md).

Please describe the purpose of this pull request in some detail and what bug it fixes. Reference and link to any relevant issues or pull requests (such as the issue in which this bug was first discussed).

# Checklist Before Requesting Reviewer

- [ ] Tests are passing
- [ ] "WIP" removed from the title of the pull request
- [ ] Selected an Assignee for the PR to be responsible for the log summary

# Before Merging

For the PR Assignees:

- [ ] Summarize commit messages into a comprehensive review of the PR
---
name: ✅  Release Checklist (Maintainers Only)
about: Checklist for core developers to complete as part of making a release

---
# Release Checklist

## Before Release

* [ ] Migrate any unresolved Issues or PRs from the [release GitHub project board](https://github.com/scikit-hep/pyhf/projects/) to a new project board.
* [ ] Verify that there is a release notes file for the release under [``docs/release-notes``](https://github.com/scikit-hep/pyhf/tree/master/docs/release-notes).
* [ ] Verify that the release notes files correctly summarize all development changes since the last release.
* [ ] Draft email to [``pyhf-announcements`` mailing list](https://groups.google.com/group/pyhf-announcements/subscribe) that summarizes the main points of the release notes and circulate it for development team approval.
* [ ] Update the checklist Issue template in the [``.github/ISSUE_TEMPLATE``](https://github.com/scikit-hep/pyhf/tree/master/.github/ISSUE_TEMPLATE) directory if there are revisions.
* [ ] Make a release to [TestPyPI][TestPyPI_pyhf] using the [workflow dispatch event trigger](https://github.com/scikit-hep/pyhf/actions/workflows/publish-package.yml).
* [ ] Verify that the project README is displaying correctly on [TestPyPI][TestPyPI_pyhf].
* [ ] Add any new use citations or published statistical models to the [Use and Citations page][citations_page].
* [ ] Verify that the citations on the [Use and Citations page][citations_page] are up to date with their current [INSPIRE](https://inspirehep.net/) record.
* [ ] Update the [pypa/gh-action-pypi-publish](https://github.com/pypa/gh-action-pypi-publish) GitHub Action used for deployment to TestPyPI and PyPI to the latest stable release.
* [ ] Update the ``codemeta.json`` file in the release PR if its requirements have updated.

[TestPyPI_pyhf]: https://test.pypi.org/project/pyhf/
[citations_page]: https://scikit-hep.org/pyhf/citations.html

## Once Release PR is Merged

* [ ] Watch the CI to ensure that the deployment to [PyPI](https://pypi.org/project/pyhf/) is successful.
* [ ] Create a [GitHub release](https://github.com/scikit-hep/pyhf/releases) from the generated PR tag and copy the release notes published to the GitHub release page. The creation of the GitHub release triggers all other release related activities.
   - [ ] Before pasting in the release notes copy the changes that the GitHub bot has already queued up and pasted into the tag and place them in the "Changes" section of the release notes. If the release notes are published before these are copied then they will be overwritten and you'll have to add them back in by hand.
* [ ] Verify there is a new [Zenodo DOI](https://doi.org/10.5281/zenodo.1169739) minted for the release.
   - [ ] Verify that the new release archive metadata on Zenodo matches is being picked up as expected from [`CITATION.cff`](https://github.com/scikit-hep/pyhf/blob/master/CITATION.cff).
* [ ] Verify that a Binder has properly built for the new release.
* [ ] Watch for a GitHub notification that there is an automatic PR to the
  [Conda-forge feedstock](https://github.com/conda-forge/pyhf-feedstock). This may take multiple hours to happen. If there are any changes needed to the Conda-forge release make them **from a personal account** and not from an organization account to have workflows properly trigger.
   - [ ] Verify the requirements in the [Conda-forge feedstock](https://github.com/conda-forge/pyhf-feedstock) recipe `meta.yaml` match those in `setup.cfg` and `pyproject.toml`.

## After Release

* [ ] Verify that the release is installable from both [PyPI](https://pypi.org/project/pyhf/) and [Conda-forge](https://github.com/conda-forge/pyhf-feedstock).
* [ ] Send the drafted ``pyhf-announcements`` email out from the ``pyhf-announcements`` account email.
* [ ] Tweet the release out on both personal and team Twitter accounts.
* [ ] Announce the release on the [Scikit-HEP community Gitter](https://gitter.im/Scikit-HEP/community).
* [ ] Make a release for the [`pyhf` tutorial](https://github.com/pyhf/pyhf-tutorial/releases) corresponding to the **previous release** number. This release represents the last version of the tutorial that is guaranteed to work with previous release API.
* [ ] Update the [tutorial](https://github.com/pyhf/pyhf-tutorial) to use the new release number and API.
* [ ] Make a PR to use the new release in the [CUDA enabled Docker images](https://github.com/pyhf/cuda-images).
* [ ] If the release is a **major** or **minor** release, open a [GitHub Release Radar](https://github.com/github/release-radar) Issue for the release to potentially get featured on GitHub's [Release Radar blog](https://github.blog/?s=release+radar).
* [ ] Close the [release GitHub Project board](https://github.com/scikit-hep/pyhf/projects/).
---
title: 'pyhf: pure-Python implementation of HistFactory statistical models'
tags:
  - Python
  - physics
  - high energy physics
  - statistical modeling
  - fitting
  - auto-differentiation
authors:
  - name: Lukas Heinrich
    orcid: 0000-0002-4048-7584
    affiliation: 1
  - name: Matthew Feickert^[Corresponding author.]
    orcid: 0000-0003-4124-7862
    affiliation: 2
  - name: Giordon Stark
    orcid: 0000-0001-6616-3433
    affiliation: 3
  - name: Kyle Cranmer
    orcid: 0000-0002-5769-7094
    affiliation: 4
affiliations:
 - name: CERN
   index: 1
 - name: University of Illinois at Urbana-Champaign
   index: 2
 - name: SCIPP, University of California, Santa Cruz
   index: 3
 - name: New York University
   index: 4
date: 5 October 2020
bibliography: paper.bib
---

# Summary

Statistical analysis of High Energy Physics (HEP) data relies on quantifying the compatibility of observed collision events with theoretical predictions.
The relationship between them is often formalised in a statistical model $f(\mathbf{x}|\mathbf{\phi})$ describing the probability of data $\mathbf{x}$ given model parameters $\mathbf{\phi}$.
Given observed data, the likelihood $\mathcal{L}(\mathbf{\phi})$ then serves as the basis for inference on the parameters $\mathbf{\phi}$.
For measurements based on binned data (histograms), the `HistFactory` family of statistical models [@Cranmer:1456844] has been widely used in both Standard Model measurements [@HIGG-2013-02] as well as searches for new physics [@ATLAS-CONF-2018-041].
`pyhf` is a pure-Python implementation of the `HistFactory` model specification and implements a declarative, plain-text format for describing `HistFactory`-based likelihoods that is targeted for reinterpretation and long-term preservation in analysis data repositories such as HEPData [@Maguire:2017ypu].
The source code for `pyhf` has been archived on Zenodo with the linked DOI: [@pyhf_zenodo].
At the time of writing this paper, the most recent release of `pyhf` is [`v0.5.4`](https://doi.org/10.5281/zenodo.4318533).

# Statement of Need

Through adoption of open source "tensor" computational Python libraries, `pyhf` decreases the abstractions between a physicist performing an analysis and the statistical modeling without sacrificing computational speed.
By taking advantage of tensor calculations, `pyhf` outperforms the traditional `C++` implementation of `HistFactory` on data from real LHC analyses.
`pyhf`'s default computational backend is built from NumPy and SciPy, and supports TensorFlow, PyTorch, and JAX as alternative backend choices.
These alternative backends support hardware acceleration on GPUs, and in the case of JAX JIT compilation, as well as auto-differentiation allowing for calculating the full gradient of the likelihood function &mdash; all contributing to speeding up fits.

## Impact on Physics

In addition to enabling the first publication of full likelihoods by an LHC experiment [@ATL-PHYS-PUB-2019-029], `pyhf` has been used by the `SModelS` library to improve the reinterpretation of results of searches for new physics at LHC experiments [@Abdallah:2020pec; @Khosa:2020zar; @Alguero:2020grj].

## Future work

Future development aims to provide support for limit setting through pseudoexperiment generation in the regimes in which asymptotic approximations [@Cowan:2010js] are no longer valid.
Further improvements to the performance of the library as well as API refinement are also planned.

# Acknowledgements

We would like to thank everyone who has made contributions to `pyhf` and thank our fellow developers in the Scikit-HEP community for their continued support and feedback.
Matthew Feickert and Kyle Cranmer have received support for work related to `pyhf` provided by NSF cooperative agreement OAC-1836650 (IRIS-HEP) and grant OAC-1450377 (DIANA/HEP).

# References
.. image:: https://raw.githubusercontent.com/scikit-hep/pyhf/master/docs/_static/img/pyhf-logo-small.png
   :alt: pyhf logo
   :width: 320
   :align: center

pure-python fitting/limit-setting/interval estimation HistFactory-style
=======================================================================

|GitHub Project| |DOI| |JOSS DOI| |Scikit-HEP| |NSF Award Number|

|Docs from latest| |Docs from master| |Binder|

|PyPI version| |Conda-forge version| |Supported Python versions| |Docker Hub pyhf| |Docker Hub pyhf CUDA|

|Code Coverage| |CodeFactor| |pre-commit.ci Status| |Code style: black|

|GitHub Actions Status: CI| |GitHub Actions Status: Docs| |GitHub Actions Status: Publish|
|GitHub Actions Status: Docker|

The HistFactory p.d.f. template
[`CERN-OPEN-2012-016 <https://cds.cern.ch/record/1456844>`__] is per-se
independent of its implementation in ROOT and sometimes, it’s useful to
be able to run statistical analysis outside of ROOT, RooFit, RooStats
framework.

This repo is a pure-python implementation of that statistical model for
multi-bin histogram-based analysis and its interval estimation is based
on the asymptotic formulas of “Asymptotic formulae for likelihood-based
tests of new physics”
[`arXiv:1007.1727 <https://arxiv.org/abs/1007.1727>`__]. The aim is also
to support modern computational graph libraries such as PyTorch and
TensorFlow in order to make use of features such as autodifferentiation
and GPU acceleration.

Hello World
-----------

This is how you use the ``pyhf`` Python API to build a statistical model and run basic inference:

.. code:: pycon

   >>> import pyhf
   >>> pyhf.set_backend("numpy")
   >>> model = pyhf.simplemodels.uncorrelated_background(
   ...     signal=[12.0, 11.0], bkg=[50.0, 52.0], bkg_uncertainty=[3.0, 7.0]
   ... )
   >>> data = [51, 48] + model.config.auxdata
   >>> test_mu = 1.0
   >>> CLs_obs, CLs_exp = pyhf.infer.hypotest(
   ...     test_mu, data, model, test_stat="qtilde", return_expected=True
   ... )
   >>> print(f"Observed: {CLs_obs:.8f}, Expected: {CLs_exp:.8f}")
   Observed: 0.05251497, Expected: 0.06445321

Alternatively the statistical model and observational data can be read from its serialized JSON representation (see next section).

.. code:: pycon

   >>> import pyhf
   >>> import requests
   >>> pyhf.set_backend("numpy")
   >>> wspace = pyhf.Workspace(requests.get("https://git.io/JJYDE").json())
   >>> model = wspace.model()
   >>> data = wspace.data(model)
   >>> test_mu = 1.0
   >>> CLs_obs, CLs_exp = pyhf.infer.hypotest(
   ...     test_mu, data, model, test_stat="qtilde", return_expected=True
   ... )
   >>> print(f"Observed: {CLs_obs:.8f}, Expected: {CLs_exp:.8f}")
   Observed: 0.35998409, Expected: 0.35998409


Finally, you can also use the command line interface that ``pyhf`` provides

.. code:: bash

   $ cat << EOF  | tee likelihood.json | pyhf cls
   {
       "channels": [
           { "name": "singlechannel",
             "samples": [
               { "name": "signal",
                 "data": [12.0, 11.0],
                 "modifiers": [ { "name": "mu", "type": "normfactor", "data": null} ]
               },
               { "name": "background",
                 "data": [50.0, 52.0],
                 "modifiers": [ {"name": "uncorr_bkguncrt", "type": "shapesys", "data": [3.0, 7.0]} ]
               }
             ]
           }
       ],
       "observations": [
           { "name": "singlechannel", "data": [51.0, 48.0] }
       ],
       "measurements": [
           { "name": "Measurement", "config": {"poi": "mu", "parameters": []} }
       ],
       "version": "1.0.0"
   }
   EOF

which should produce the following JSON output:

.. code:: json

   {
      "CLs_exp": [
         0.0026062609501074576,
         0.01382005356161206,
         0.06445320535890459,
         0.23525643861460702,
         0.573036205919389
      ],
      "CLs_obs": 0.05251497423736956
   }

What does it support
--------------------

Implemented variations:
  - ☑ HistoSys
  - ☑ OverallSys
  - ☑ ShapeSys
  - ☑ NormFactor
  - ☑ Multiple Channels
  - ☑ Import from XML + ROOT via `uproot <https://github.com/scikit-hep/uproot4>`__
  - ☑ ShapeFactor
  - ☑ StatError
  - ☑ Lumi Uncertainty
  - ☑ Non-asymptotic calculators

Computational Backends:
  - ☑ NumPy
  - ☑ PyTorch
  - ☑ TensorFlow
  - ☑ JAX

Optimizers:
  - ☑ SciPy (``scipy.optimize``)
  - ☑ MINUIT (``iminuit``)

All backends can be used in combination with all optimizers.
Custom user backends and optimizers can be used as well.

Todo
----

-  ☐ StatConfig

results obtained from this package are validated against output computed
from HistFactory workspaces

A one bin example
-----------------

.. code:: python

   import pyhf
   import numpy as np
   import matplotlib.pyplot as plt
   from pyhf.contrib.viz import brazil

   pyhf.set_backend("numpy")
   model = pyhf.simplemodels.uncorrelated_background(
       signal=[10.0], bkg=[50.0], bkg_uncertainty=[7.0]
   )
   data = [55.0] + model.config.auxdata

   poi_vals = np.linspace(0, 5, 41)
   results = [
       pyhf.infer.hypotest(
           test_poi, data, model, test_stat="qtilde", return_expected_set=True
       )
       for test_poi in poi_vals
   ]

   fig, ax = plt.subplots()
   fig.set_size_inches(7, 5)
   brazil.plot_results(poi_vals, results, ax=ax)
   fig.show()

**pyhf**

.. image:: https://raw.githubusercontent.com/scikit-hep/pyhf/master/docs/_static/img/README_1bin_example.png
   :alt: manual
   :width: 500
   :align: center

**ROOT**

.. image:: https://raw.githubusercontent.com/scikit-hep/pyhf/master/docs/_static/img/hfh_1bin_55_50_7.png
   :alt: manual
   :width: 500
   :align: center

A two bin example
-----------------

.. code:: python

   import pyhf
   import numpy as np
   import matplotlib.pyplot as plt
   from pyhf.contrib.viz import brazil

   pyhf.set_backend("numpy")
   model = pyhf.simplemodels.uncorrelated_background(
       signal=[30.0, 45.0], bkg=[100.0, 150.0], bkg_uncertainty=[15.0, 20.0]
   )
   data = [100.0, 145.0] + model.config.auxdata

   poi_vals = np.linspace(0, 5, 41)
   results = [
       pyhf.infer.hypotest(
           test_poi, data, model, test_stat="qtilde", return_expected_set=True
       )
       for test_poi in poi_vals
   ]

   fig, ax = plt.subplots()
   fig.set_size_inches(7, 5)
   brazil.plot_results(poi_vals, results, ax=ax)
   fig.show()


**pyhf**

.. image:: https://raw.githubusercontent.com/scikit-hep/pyhf/master/docs/_static/img/README_2bin_example.png
   :alt: manual
   :width: 500
   :align: center

**ROOT**

.. image:: https://raw.githubusercontent.com/scikit-hep/pyhf/master/docs/_static/img/hfh_2_bin_100.0_145.0_100.0_150.0_15.0_20.0_30.0_45.0.png
   :alt: manual
   :width: 500
   :align: center

Installation
------------

To install ``pyhf`` from PyPI with the NumPy backend run

.. code:: bash

   python -m pip install pyhf

and to install ``pyhf`` with all additional backends run

.. code:: bash

   python -m pip install pyhf[backends]

or a subset of the options.

To uninstall run

.. code:: bash

   python -m pip uninstall pyhf

Questions
---------

If you have a question about the use of ``pyhf`` not covered in `the
documentation <https://pyhf.readthedocs.io/>`__, please ask a question
on the `GitHub Discussions <https://github.com/scikit-hep/pyhf/discussions>`__.

If you believe you have found a bug in ``pyhf``, please report it in the
`GitHub
Issues <https://github.com/scikit-hep/pyhf/issues/new?template=Bug-Report.md&labels=bug&title=Bug+Report+:+Title+Here>`__.
If you're interested in getting updates from the ``pyhf`` dev team and release
announcements you can join the |pyhf-announcements mailing list|_.

.. |pyhf-announcements mailing list| replace:: ``pyhf-announcements`` mailing list
.. _pyhf-announcements mailing list: https://groups.google.com/group/pyhf-announcements/subscribe

Citation
--------

As noted in `Use and Citations <https://scikit-hep.org/pyhf/citations.html>`__,
the preferred BibTeX entry for citation of ``pyhf`` includes both the
`Zenodo <https://zenodo.org/>`__ archive and the
`JOSS <https://joss.theoj.org/>`__ paper:

.. code:: bibtex

   @software{pyhf,
     author = {Lukas Heinrich and Matthew Feickert and Giordon Stark},
     title = "{pyhf: v0.6.3}",
     version = {0.6.3},
     doi = {10.5281/zenodo.1169739},
     url = {https://doi.org/10.5281/zenodo.1169739},
     note = {https://github.com/scikit-hep/pyhf/releases/tag/v0.6.3}
   }

   @article{pyhf_joss,
     doi = {10.21105/joss.02823},
     url = {https://doi.org/10.21105/joss.02823},
     year = {2021},
     publisher = {The Open Journal},
     volume = {6},
     number = {58},
     pages = {2823},
     author = {Lukas Heinrich and Matthew Feickert and Giordon Stark and Kyle Cranmer},
     title = {pyhf: pure-Python implementation of HistFactory statistical models},
     journal = {Journal of Open Source Software}
   }

Authors
-------

``pyhf`` is openly developed by Lukas Heinrich, Matthew Feickert, and Giordon Stark.

Please check the `contribution statistics for a list of
contributors <https://github.com/scikit-hep/pyhf/graphs/contributors>`__.

Milestones
----------

- 2021-12-09: 1000 commits to the project. (See PR `#1710 <https://github.com/scikit-hep/pyhf/pull/1710>`__)
- 2020-07-28: 1000 GitHub issues and pull requests. (See PR `#1000 <https://github.com/scikit-hep/pyhf/pull/1000>`__)

Acknowledgements
----------------

Matthew Feickert has received support to work on ``pyhf`` provided by NSF
cooperative agreement `OAC-1836650 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1836650>`__ (IRIS-HEP)
and grant `OAC-1450377 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1450377>`__ (DIANA/HEP).

.. |GitHub Project| image:: https://img.shields.io/badge/GitHub--blue?style=social&logo=GitHub
   :target: https://github.com/scikit-hep/pyhf
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1169739.svg
   :target: https://doi.org/10.5281/zenodo.1169739
.. |JOSS DOI| image:: https://joss.theoj.org/papers/10.21105/joss.02823/status.svg
   :target: https://doi.org/10.21105/joss.02823
.. |Scikit-HEP| image:: https://scikit-hep.org/assets/images/Scikit--HEP-Project-blue.svg
   :target: https://scikit-hep.org/
.. |NSF Award Number| image:: https://img.shields.io/badge/NSF-1836650-blue.svg
   :target: https://nsf.gov/awardsearch/showAward?AWD_ID=1836650
.. |Docs from latest| image:: https://img.shields.io/badge/docs-v0.6.3-blue.svg
   :target: https://pyhf.readthedocs.io/
.. |Docs from master| image:: https://img.shields.io/badge/docs-master-blue.svg
   :target: https://scikit-hep.github.io/pyhf
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/scikit-hep/pyhf/master?filepath=docs%2Fexamples%2Fnotebooks%2Fbinderexample%2FStatisticalAnalysis.ipynb

.. |PyPI version| image:: https://badge.fury.io/py/pyhf.svg
   :target: https://badge.fury.io/py/pyhf
.. |Conda-forge version| image:: https://img.shields.io/conda/vn/conda-forge/pyhf.svg
   :target: https://github.com/conda-forge/pyhf-feedstock
.. |Supported Python versions| image:: https://img.shields.io/pypi/pyversions/pyhf.svg
   :target: https://pypi.org/project/pyhf/
.. |Docker Hub pyhf| image:: https://img.shields.io/badge/pyhf-v0.6.3-blue?logo=Docker
   :target: https://hub.docker.com/r/pyhf/pyhf/tags
.. |Docker Hub pyhf CUDA| image:: https://img.shields.io/badge/pyhf-CUDA-blue?logo=Docker
   :target: https://hub.docker.com/r/pyhf/cuda/tags

.. |Code Coverage| image:: https://codecov.io/gh/scikit-hep/pyhf/graph/badge.svg?branch=master
   :target: https://codecov.io/gh/scikit-hep/pyhf?branch=master
.. |CodeFactor| image:: https://www.codefactor.io/repository/github/scikit-hep/pyhf/badge
   :target: https://www.codefactor.io/repository/github/scikit-hep/pyhf
.. |pre-commit.ci Status| image:: https://results.pre-commit.ci/badge/github/scikit-hep/pyhf/master.svg
  :target: https://results.pre-commit.ci/latest/github/scikit-hep/pyhf/master
  :alt: pre-commit.ci status
.. |Code style: black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black

.. |GitHub Actions Status: CI| image:: https://github.com/scikit-hep/pyhf/workflows/CI/CD/badge.svg?branch=master
   :target: https://github.com/scikit-hep/pyhf/actions?query=workflow%3ACI%2FCD+branch%3Amaster
.. |GitHub Actions Status: Docs| image:: https://github.com/scikit-hep/pyhf/workflows/Docs/badge.svg?branch=master
   :target: https://github.com/scikit-hep/pyhf/actions?query=workflow%3ADocs+branch%3Amaster
.. |GitHub Actions Status: Publish| image:: https://github.com/scikit-hep/pyhf/workflows/publish%20distributions/badge.svg?branch=master
   :target: https://github.com/scikit-hep/pyhf/actions?query=workflow%3A%22publish+distributions%22+branch%3Amaster
.. |GitHub Actions Status: Docker| image:: https://github.com/scikit-hep/pyhf/actions/workflows/docker.yml/badge.svg?branch=master
   :target: https://github.com/scikit-hep/pyhf/actions/workflows/docker.yml?query=branch%3Amaster
==========
Developing
==========

Developer Environment
---------------------

To develop, we suggest using `virtual environments <https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments>`__ together with ``pip`` or using `pipenv <https://pipenv.readthedocs.io/en/latest/>`__. Once the environment is activated, clone the repo from GitHub

.. code-block:: console

    git clone https://github.com/scikit-hep/pyhf.git

and install all necessary packages for development

.. code-block:: console

    python -m pip install --upgrade --editable .[complete]

Then setup the Git pre-commit hook for `Black <https://github.com/psf/black>`__  by running

.. code-block:: console

    pre-commit install

as the ``rev`` gets updated through time to track changes of different hooks,
simply run

.. code-block:: console

    pre-commit autoupdate

to have pre-commit install the new version.

Testing
-------

Data Files
~~~~~~~~~~

A function-scoped fixture called ``datadir`` exists for a given test module
which will automatically copy files from the associated test modules data
directory into a temporary directory for the given test execution. That is, for
example, if a test was defined in ``test_schema.py``, then data files located
in ``test_schema/`` will be copied to a temporary directory whose path is made
available by the ``datadir`` fixture. Therefore, one can do:

.. code-block:: python

    def test_patchset(datadir):
        data_file = open(datadir.join("test.txt"))
        ...

which will load the copy of ``text.txt`` in the temporary directory. This also
works for parameterizations as this will effectively sandbox the file
modifications made.

TestPyPI
~~~~~~~~

``pyhf`` tests packaging and distributing by publishing in advance of releases
to `TestPyPI <https://test.pypi.org/project/pyhf/>`__.
Installation of the latest test release from TestPyPI can be tested
by first installing ``pyhf`` normally, to ensure all dependencies are installed
from PyPI, and then upgrading ``pyhf`` to a test release from TestPyPI

.. code-block:: bash

  python -m pip install pyhf
  python -m pip install --upgrade --extra-index-url https://test.pypi.org/simple/ --pre pyhf

.. note::

  This adds TestPyPI as `an additional package index to search <https://pip.pypa.io/en/stable/reference/pip_install/#cmdoption-extra-index-url>`__
  when installing.
  PyPI will still be the default package index ``pip`` will attempt to install
  from for all dependencies, but if a package has a release on TestPyPI that
  is a more recent release then the package will be installed from TestPyPI instead.
  Note that dev releases are considered pre-releases, so ``0.1.2`` is a "newer"
  release than ``0.1.2.dev3``.

Publishing
----------

Publishing to `PyPI <https://pypi.org/project/pyhf/>`__ and `TestPyPI <https://test.pypi.org/project/pyhf/>`__
is automated through the `PyPA's PyPI publish GitHub Action <https://github.com/pypa/gh-action-pypi-publish>`__
and the ``pyhf`` `Tag Creator GitHub Actions workflow <https://github.com/scikit-hep/pyhf/blob/master/.github/workflows/tag.yml>`__.
A release can be created from any PR created by a core developer by adding a
``bumpversion`` tag to it that corresponds to the release type:
`major <https://github.com/scikit-hep/pyhf/labels/bumpversion%2Fmajor>`__,
`minor <https://github.com/scikit-hep/pyhf/labels/bumpversion%2Fminor>`__,
`patch <https://github.com/scikit-hep/pyhf/labels/bumpversion%2Fpatch>`__.
Once the PR is tagged with the label, the GitHub Actions bot will post a comment
with information on the actions it will take once the PR is merged. When the PR
has been reviewed, approved, and merged, the Tag Creator workflow will automatically
create a new release with ``bump2version`` and then deploy the release to PyPI.

Context Files and Archive Metadata
----------------------------------

The ``.zenodo.json`` and ``codemeta.json`` files have the version number
automatically updated through ``bump2version``, though their additional metadata
should be checked periodically by the dev team (probably every release).
The ``codemeta.json`` file can be generated automatically **from a PyPI install**
of ``pyhf`` using ``codemetapy``

.. code-block:: bash

  codemetapy --no-extras pyhf > codemeta.json

though the ``author`` metadata will still need to be checked and revised by hand.
The ``.zenodo.json`` is currently generated by hand, so it is worth using
``codemeta.json`` as a guide to edit it.

Release Checklist
-----------------

As part of the release process a checklist is required to be completed to make
sure steps aren't missed.
There is a GitHub Issue template for this that the developer in charge of the
release should step through and update if needed.
Use and Citations
=================

.. raw:: html

   <p id="dev-version"><strong>Warning:</strong> This is a development version and should not be cited. To find the specific version to cite, please go to <a href="https://pyhf.readthedocs.io/">ReadTheDocs</a>.</p>

Citation
--------

The preferred BibTeX entry for citation of ``pyhf`` includes both the  `Zenodo <https://zenodo.org/>`__
archive and the `JOSS <https://joss.theoj.org/>`__ paper:

.. literalinclude:: bib/preferred.bib
   :language: bibtex

Use in Publications
-------------------

The following is an updating list of citations and use cases of :code:`pyhf`.
There is an incomplete but automatically updated `list of citations on INSPIRE
<https://inspirehep.net/literature/1845084>`__ as well.

Use Citations
~~~~~~~~~~~~~

.. bibliography:: bib/use_citations.bib
   :list: enumerated
   :all:
   :style: unsrt

General Citations
~~~~~~~~~~~~~~~~~

.. bibliography:: bib/general_citations.bib
   :list: enumerated
   :all:
   :style: unsrt

Published Statistical Models
----------------------------

Updating list of HEPData entries for publications using ``HistFactory`` JSON statistical models:

.. bibliography:: bib/HEPData_likelihoods.bib
   :list: enumerated
   :all:
   :style: unsrt

.. note::

   ATLAS maintains a public listing of all published statistical models on the `ATLAS public results
   page <https://twiki.cern.ch/twiki/bin/view/AtlasPublic>`__ which can be found by filtering all
   public results by the "Likelihood available" analysis characteristics keyword.
.. _sec:likelihood:

Likelihood Specification
========================

The structure of the JSON specification of models follows closely the
original XML-based specification :cite:`likelihood-Cranmer:1456844`.

Workspace
---------

.. literalinclude:: ../src/pyhf/schemas/1.0.0/workspace.json
   :language: json

The overall document in the above code snippet describes a *workspace*, which includes

* **channels**: The channels in the model, which include a description of the samples
  within each channel and their possible parametrised modifiers.
* **measurements**: A set of measurements, which define among others the parameters of
  interest for a given statistical analysis objective.
* **observations**: The observed data, with which a likelihood can be constructed from the model.

A workspace consists of the channels, one set of observed data, but can
include multiple measurements. If provided a JSON file, one can quickly
check that it conforms to the provided workspace specification as follows:

.. code:: python

   import json, requests, jsonschema

   workspace = json.load(open("/path/to/analysis_workspace.json"))
   # if no exception is raised, it found and parsed the schema
   schema = requests.get("https://scikit-hep.org/pyhf/schemas/1.0.0/workspace.json").json()
   # If no exception is raised by validate(), the instance is valid.
   jsonschema.validate(instance=workspace, schema=schema)


.. _ssec:channel:

Channel
-------

A channel is defined by a channel name and a list of samples :cite:`likelihood-schema_defs`.

.. code:: json

    {
        "channel": {
            "type": "object",
            "properties": {
                "name": { "type": "string" },
                "samples": { "type": "array", "items": {"$ref": "#/definitions/sample"}, "minItems": 1 }
            },
            "required": ["name", "samples"],
            "additionalProperties": false
        },
    }

The Channel specification consists of a list of channel descriptions.
Each channel, an analysis region encompassing one or more measurement
bins, consists of a ``name`` field and a ``samples`` field (see :ref:`ssec:channel`), which
holds a list of sample definitions (see :ref:`ssec:sample`). Each sample definition in
turn has a ``name`` field, a ``data`` field for the nominal event rates
for all bins in the channel, and a ``modifiers`` field of the list of
modifiers for the sample.

.. _ssec:sample:

Sample
------

A sample is defined by a sample name, the sample event rate, and a list of modifiers :cite:`likelihood-schema_defs`.

.. _lst:schema:sample:

.. code:: json

    {
        "sample": {
            "type": "object",
            "properties": {
                "name": { "type": "string" },
                "data": { "type": "array", "items": {"type": "number"}, "minItems": 1 },
                "modifiers": {
                    "type": "array",
                    "items": {
                        "anyOf": [
                            { "$ref": "#/definitions/modifier/histosys" },
                            { "$ref": "#/definitions/modifier/lumi" },
                            { "$ref": "#/definitions/modifier/normfactor" },
                            { "$ref": "#/definitions/modifier/normsys" },
                            { "$ref": "#/definitions/modifier/shapefactor" },
                            { "$ref": "#/definitions/modifier/shapesys" },
                            { "$ref": "#/definitions/modifier/staterror" }
                        ]
                    }
                }
            },
            "required": ["name", "data", "modifiers"],
            "additionalProperties": false
        },
    }

Modifiers
---------

The modifiers that are applicable for a given sample are encoded as a
list of JSON objects with three fields. A name field, a type field
denoting the class of the modifier, and a data field which provides the
necessary input data as denoted in :ref:`tab:modifiers_and_constraints`.

Based on the declared modifiers, the set of parameters and their
constraint terms are derived implicitly as each type of modifier
unambiguously defines the constraint terms it requires. Correlated shape
modifiers and normalisation uncertainties have compatible constraint
terms and thus modifiers can be declared that *share* parameters by
re-using a name [1]_ for multiple modifiers. That is, a variation of a
single parameter causes a shift within sample rates due to both shape
and normalisation variations.

We review the structure of each modifier type below.

Uncorrelated Shape (shapesys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To construct the constraint term, the relative uncertainties
:math:`\sigma_b` are necessary for each bin. Therefore, we record the
absolute uncertainty as an array of floats, which combined with the
nominal sample data yield the desired :math:`\sigma_b`.
An example of an uncorrelated shape modifier with three absolute uncertainty
terms for a 3-bin channel is shown below:

.. code:: json

   { "name": "mod_name", "type": "shapesys", "data": [1.0, 1.5, 2.0] }

.. warning::

   For bins in the model where:

     * the samples nominal expected rate is zero, or
     * the absolute uncertainty is zero.

   nuisance parameters will be allocated, but will be fixed to ``1`` in the
   calculation (as shapesys is a multiplicative modifier this results in
   multiplying by ``1``).

   These values are, in the context of uncorrelated shape uncertainties,
   unphysical. If this situation occurs, one needs to go back and understand
   the inputs as this is undefined behavior in HistFactory.

The previous example will allocate three nuisance parameters for ``mod_name``.
The following example will also allocate three nuisance parameters for a 3-bin
channel, with the second nuisance parameter fixed to ``1``:

.. code:: json

   { "name": "mod_name", "type": "shapesys", "data": [1.0, 0.0, 2.0] }

Correlated Shape (histosys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This modifier represents the same source of uncertainty which has a
different effect on the various sample shapes, hence a correlated shape.
To implement an interpolation between sample distribution shapes, the
distributions with a "downward variation" ("lo") associated with
:math:`\alpha=-1` and an "upward variation" ("hi") associated with
:math:`\alpha=+1` are provided as arrays of floats.
An example of a correlated shape modifier with absolute shape variations
for a 2-bin channel is shown below:

.. code:: json

   { "name": "mod_name", "type": "histosys", "data": {"hi_data": [20,15], "lo_data": [10, 10]} }

Normalisation Uncertainty (normsys)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The normalisation uncertainty modifies the sample rate by a overall
factor :math:`\kappa(\alpha)` constructed as the interpolation between
downward ("lo") and upward ("hi") as well as the nominal setting, i.e.
:math:`\kappa(-1) = \kappa_{\alpha=-1}`, :math:`\kappa(0) = 1` and
:math:`\kappa(+1) = \kappa_{\alpha=+1}`. In the modifier definition we record
:math:`\kappa_{\alpha=+1}` and :math:`\kappa_{\alpha=-1}` as floats.
An example of a normalisation uncertainty modifier with scale factors recorded
for the up/down variations of an :math:`n`-bin channel is shown below:

.. code:: json

   { "name": "mod_name", "type": "normsys", "data": {"hi": 1.1, "lo": 0.9} }

MC Statistical Uncertainty (staterror)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As the sample counts are often derived from Monte Carlo (MC) datasets, they
necessarily carry an uncertainty due to the finite sample size of the datasets.
As explained in detail in :cite:`likelihood-Cranmer:1456844`, adding uncertainties for
each sample would yield a very large number of nuisance parameters with limited
utility. Therefore a set of bin-wise scale factors :math:`\gamma_b` is
introduced to model the overall uncertainty in the bin due to MC statistics.
The constrained term is constructed as a set of Gaussian constraints with a
central value equal to unity for each bin in the channel. The scales
:math:`\sigma_b` of the constraint are computed from the individual
uncertainties of samples defined within the channel relative to the total event
rate of all samples: :math:`\delta_{csb} = \sigma_{csb}/\sum_s \nu^0_{scb}`. As
not all samples are within a channel are estimated from MC simulations, only
the samples with a declared statistical uncertainty modifier enter the sum.
An example of a statistical uncertainty modifier for a single bin channel is
shown below:

.. code:: json

   { "name": "mod_name", "type": "staterror", "data": [0.1] }

.. warning::

   For bins in the model where:

     * the samples nominal expected rate is zero, or
     * the scale factor is zero.

   nuisance parameters will be allocated, but will be fixed to ``1`` in the
   calculation (as staterror is a multiplicative modifier this results in
   multiplying by ``1``).

Luminosity (lumi)
~~~~~~~~~~~~~~~~~

Sample rates derived from theory calculations, as opposed to data-driven
estimates, are scaled to the integrated luminosity corresponding to the
observed data. As the luminosity measurement is itself subject to an
uncertainty, it must be reflected in the rate estimates of such samples.  As
this modifier is of global nature, no additional per-sample information is
required and thus the data field is nulled. This uncertainty is relevant, in
particular, when the parameter of interest is a signal cross-section. The
luminosity uncertainty :math:`\sigma_\lambda` is provided as part of the
parameter configuration included in the measurement specification discussed
in :ref:`ssec:measurements`.
An example of a luminosity modifier is shown below:

.. code:: json

   { "name": "mod_name", "type": "lumi", "data": null }

Unconstrained Normalisation (normfactor)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The unconstrained normalisation modifier scales the event rates of a
sample by a free parameter :math:`\mu`. Common use cases are the signal
rate of a possible BSM signal or simultaneous in-situ measurements of
background samples. Such parameters are frequently the parameters of
interest of a given measurement. No additional per-sample data is
required.
An example of a normalisation modifier is shown below:

.. code:: json

   { "name": "mod_name", "type": "normfactor", "data": null }

Data-driven Shape (shapefactor)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to support data-driven estimation of sample rates (e.g. for
multijet backgrounds), the data-driven shape modifier adds free,
bin-wise multiplicative parameters. Similarly to the normalisation
factors, no additional data is required as no constraint is defined.
An example of an uncorrelated shape modifier is shown below:

.. code:: json

   { "name": "mod_name", "type": "shapefactor", "data": null }

Data
----

The data provided by the analysis are the observed data for each channel
(or region). This data is provided as a mapping from channel name to an
array of floats, which provide the observed rates in each bin of the
channel. The auxiliary data is not included as it is an input to the
likelihood that does not need to be archived and can be determined
automatically from the specification.
An example of channel data is shown below:

.. _lst:example:data:

.. code:: json

   { "chan_name_one": [10, 20], "chan_name_two": [4, 0]}

.. _ssec:measurements:

Measurements
------------

Given the data and the model definitions, a measurement can be defined.
In the current schema, the measurements defines the name of the
parameter of interest as well as parameter set configurations.  [2]_
Here, the remaining information not covered through the channel
definition is provided, e.g. for the luminosity parameter. For all
modifiers, the default settings can be overridden where possible:

* **inits**: Initial value of the parameter.
* **bounds**: Interval bounds of the parameter.
* **auxdata**: Auxiliary data for the associated constraint term.
* **sigmas**: Associated uncertainty of the parameter.

An example of a measurement is shown below:

.. code:: json

   {
       "name": "MyMeasurement",
       "config": {
           "poi": "SignalCrossSection", "parameters": [
               { "name":"lumi", "auxdata":[1.0],"sigmas":[0.017], "bounds":[[0.915,1.085]],"inits":[1.0] },
               { "name":"mu_ttbar", "bounds":[[0, 5]] },
               { "name":"rw_1CR", "fixed":true }
           ]
       }
   }

This measurement, which scans over the parameter of interest ``SignalCrossSection``, is setting configurations for the luminosity modifier, changing the default bounds for the normfactor modifier named ``mu_ttbar``, and specifying that the modifier ``rw_1CR`` is held constant (``fixed``).

.. _ssec:observations:

Observations
------------

This is what we evaluate the hypothesis testing against, to determine the
compatibility of signal+background hypothesis to the background-only
hypothesis. This is specified as a list of objects, with each object structured
as

* **name**: the channel for which the observations are recorded
* **data**: the bin-by-bin observations for the named channel

An example of an observation for a 2-bin channel ``channel1``, with values
``110.0`` and ``120.0`` is shown below:

.. code:: json

   {
       "name": "channel1", "data": [110.0, 120.0]
   }

Toy Example
-----------

.. # N.B. If the following literalinclude is changed test_examples.py must be changed accordingly
.. literalinclude:: ./examples/json/2-bin_1-channel.json
   :language: json

In the above example, we demonstrate a simple measurement of a
single two-bin channel with two samples: a signal sample and a background
sample. The signal sample has an unconstrained normalisation factor
:math:`\mu`, while the background sample carries an uncorrelated shape
systematic controlled by parameters :math:`\gamma_1` and :math:`\gamma_2`. The
background uncertainty for the bins is 10% and 20% respectively.

Additional Material
-------------------

Footnotes
~~~~~~~~~

.. [1]
   The name of a modifier specifies the parameter set it is controlled
   by. Modifiers with the same name share parameter sets.

.. [2]
   In this context a parameter set corresponds to a named
   lower-dimensional subspace of the full parameters :math:`\fullset`.
   In many cases these are one-dimensional subspaces, e.g. a specific
   interpolation parameter :math:`\alpha` or the luminosity parameter
   :math:`\lambda`. For multi-bin channels, however, e.g. all bin-wise
   nuisance parameters of the uncorrelated shape modifiers are grouped
   under a single name. Therefore in general a parameter set definition
   provides arrays of initial values, bounds, etc.

Bibliography
~~~~~~~~~~~~

.. bibliography:: bib/docs.bib
   :filter: docname in docnames
   :style: plain
   :keyprefix: likelihood-
   :labelprefix: likelihood-
Contributors
============

``pyhf`` is openly developed and benefits from the contributions and feedback
from its users.
The ``pyhf`` dev team would like to thank all contributors to the project for
their support and help.
Thank you!

Contributors include:

- Jessica Forde
- Ruggero Turra
- Tadej Novak
- Frank Sauerburger
- Lars Nielsen
- Kanishk Kalra
- Nikolai Hartmann
- Alexander Held
- Karthikeyan Singaravelan
- Marco Gorelli
- Pradyumna Rahul K
- Eric Schanet
- Henry Schreiner
- Saransh Chopra
- Sviatoslav Sydorenko
- Mason Proffitt
- Lars Henkelmann
- Aryan Roy
Examples
========

Try out in Binder! |Binder|

.. |Binder| image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/scikit-hep/pyhf/master?filepath=docs%2Fexamples%2Fnotebooks%2Fbinderexample%2FStatisticalAnalysis.ipynb

Notebooks:

.. toctree::
   :maxdepth: 2
   :glob:

   examples/notebooks/*
   examples/notebooks/learn/*
   examples/notebooks/binderexample/*
=============
Release Notes
=============

.. include:: release-notes/v0.6.3.rst
.. include:: release-notes/v0.6.2.rst
.. include:: release-notes/v0.6.1.rst
.. include:: release-notes/v0.6.0.rst
.. include:: release-notes/v0.5.4.rst
.. include:: release-notes/v0.5.3.rst
Command Line API
================

.. click:: pyhf.cli.cli:pyhf
   :prog: pyhf
   :show-nested:
Outreach
========

We are always interested in talking about :code:`pyhf`. See the abstract and a list of previously given presentations and feel free to invite us to your next conference/workshop/meeting!

Abstract
--------

    The HistFactory p.d.f. template `[CERN-OPEN-2012-016]
    <https://cds.cern.ch/record/1456844>`_ is per-se independent of its
    implementation in ROOT and it is useful to be able to run statistical
    analysis outside of the ROOT, RooFit, RooStats framework. pyhf is a
    pure-python implementation of that statistical model for multi-bin
    histogram-based analysis and its interval estimation is based on the
    asymptotic formulas of "Asymptotic formulae for likelihood-based tests of
    new physics" :xref:`arXiv:1007.1727`.
    pyhf supports modern computational graph libraries such as TensorFlow,
    PyTorch, and JAX in order to make use of features such as
    auto-differentiation and GPU acceleration.


    .. code-block:: latex

        The HistFactory p.d.f. template
        \href{https://cds.cern.ch/record/1456844}{[CERN-OPEN-2012-016]} is
        per-se independent of its implementation in ROOT and it is useful to be
        able to run statistical analysis outside of the ROOT, RooFit, RooStats
        framework. pyhf is a pure-python implementation of that statistical
        model for multi-bin histogram-based analysis and its interval
        estimation is based on the asymptotic formulas of "Asymptotic formulae
        for likelihood-based tests of new physics"
        \href{https://arxiv.org/abs/1007.1727}{[arXiv:1007.1727]}. pyhf
        supports modern computational graph libraries such as TensorFlow,
        PyTorch, and JAX in order to make use of features such as
        auto-differentiation and GPU acceleration.


Presentations
-------------

This list will be updated with talks given on :code:`pyhf`:

.. bibliography:: bib/talks.bib
   :list: bullet
   :all:
   :style: unsrt

Tutorials
---------

This list will be updated with tutorials and schools given on :code:`pyhf`:

.. bibliography:: bib/tutorials.bib
   :list: bullet
   :all:
   :style: unsrt


Posters
-------

This list will be updated with posters presented on :code:`pyhf`:

.. bibliography:: bib/posters.bib
   :list: bullet
   :all:
   :style: unsrt

In the Media
------------

This list will be updated with media publications featuring :code:`pyhf`:

.. bibliography:: bib/media.bib
   :list: bullet
   :all:
   :style: unsrt
Fundamentals
============

Notebooks:

.. toctree::
   :maxdepth: 2
   :glob:

   examples/notebooks/learn/*
..  _installation:

Installation
============

To install, we suggest first setting up a `virtual environment <https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments>`__

.. code-block:: console

    # Python3
    python3 -m venv pyhf

and activating it

.. code-block:: console

    source pyhf/bin/activate


Install latest stable release from `PyPI <https://pypi.org/project/pyhf/>`__...
-------------------------------------------------------------------------------

... with NumPy backend
++++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf

... with TensorFlow backend
+++++++++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf[tensorflow]

... with PyTorch backend
++++++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf[torch]

... with JAX backend
++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf[jax]

... with all backends
+++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf[backends]


... with xml import/export functionality
++++++++++++++++++++++++++++++++++++++++

.. code-block:: console

    python -m pip install pyhf[xmlio]


Install latest development version from `GitHub <https://github.com/scikit-hep/pyhf>`__...
------------------------------------------------------------------------------------------

... with NumPy backend
++++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf"

... with TensorFlow backend
+++++++++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf[tensorflow]"

... with PyTorch backend
++++++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf[torch]"

... with JAX backend
++++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf[jax]"

... with all backends
+++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf[backends]"


... with xml import/export functionality
++++++++++++++++++++++++++++++++++++++++

.. code-block:: console

    python -m pip install --upgrade "git+https://github.com/scikit-hep/pyhf.git#egg=pyhf[xmlio]"


Updating :code:`pyhf`
---------------------

Rerun the installation command. As the upgrade flag (:code:`-U`, :code:`--upgrade`) is used then the libraries will be updated.
.. _sec:faq:

FAQ
===

Frequently Asked Questions about :code:`pyhf` and its use.

Questions
---------

Where can I ask questions about ``pyhf`` use?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you have a question about the use of ``pyhf`` not covered in the `documentation <https://pyhf.readthedocs.io/>`__, please ask a question on the `GitHub Discussions <https://github.com/scikit-hep/pyhf/discussions>`__.

If you believe you have found a bug in ``pyhf``, please report it in the `GitHub Issues <https://github.com/scikit-hep/pyhf/issues/new?template=Bug-Report.md&labels=bug&title=Bug+Report+:+Title+Here>`__.

How can I get updates on ``pyhf``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you're interested in getting updates from the ``pyhf`` dev team and release
announcements you can join the |pyhf-announcements mailing list|_.

.. |pyhf-announcements mailing list| replace:: ``pyhf-announcements`` mailing list
.. _pyhf-announcements mailing list: https://groups.google.com/group/pyhf-announcements/subscribe

Is it possible to set the backend from the CLI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Yes.
Use the :code:`--backend` option for :code:`pyhf cls` to specify a tensor backend.
The default backend is NumPy.
For more information see :code:`pyhf cls --help`.

Does ``pyhf`` support Python 2?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No.
Like the rest of the Python community, as of January 2020 the latest releases of ``pyhf`` no longer support Python 2.
The last release of ``pyhf`` that was compatible with Python 2.7 is `v0.3.4 <https://pypi.org/project/pyhf/0.3.4/>`__, which can be installed with

    .. code-block:: console

        python -m pip install pyhf~=0.3

I only have access to Python 2. How can I use ``pyhf``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended that ``pyhf`` is used as a standalone step in any analysis, and its environment need not be the same as the rest of the analysis.
As Python 2 is not supported it is suggested that you setup a Python 3 runtime on whatever machine you're using.
If you're using a cluster, talk with your system administrators to get their help in doing so.
If you are unable to get a Python 3 runtime, versioned Docker images of ``pyhf`` are distributed through `Docker Hub <https://hub.docker.com/r/pyhf/pyhf>`__.

Once you have Python 3 installed, see the :ref:`installation` page to get started.

I validated my workspace by comparing ``pyhf`` and ``HistFactory``, and while the expected CLs matches, the observed CLs is different. Why is this?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you're using the right test statistic (:math:`q` or :math:`\tilde{q}`) in both situations.
In ``HistFactory``, the asymptotics calculator, for example, will do something more involved for the observed CLs if you choose a different test statistic.

I ran validation to compare ``HistFitter`` and ``pyhf``, but they don't match exactly. Why not?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``pyhf`` is validated against ``HistFactory``.
``HistFitter`` makes some particular implementation choices that ``pyhf`` doesn't reproduce.
Instead of trying to compare ``pyhf`` and ``HistFitter`` you should try to validate them both against ``HistFactory``.

How is ``pyhf`` typeset?
~~~~~~~~~~~~~~~~~~~~~~~~

As you may have guessed from this page, ``pyhf`` is typeset in all lowercase.
This is largely historical, as the core developers had just always typed it that way and it seemed a bit too short of a library name to write as ``PyHF``.
When typesetting in LaTeX the developers recommend introducing the command

    .. code-block:: latex

        \newcommand{\pyhf}{\texttt{pyhf}}

If the journal you are publishing in requires you to use ``textsc`` for software names it is okay to instead use

    .. code-block:: latex

        \newcommand{\pyhf}{\textsc{pyhf}}

Why use Python?
~~~~~~~~~~~~~~~

As of the late 2010's Python is widely considered the lingua franca of machine learning
libraries, and is sufficiently high-level and expressive for physicists of various computational
skill backgrounds to use.
Using Python as the language for development allows for the distribution of the software
--- as both source files and binary distributions --- through the Python Package Index (PyPI)
and Conda-forge, which significantly lowers the barrier for use as compared to ``C++``.
Additionally, a 2017 `DIANA/HEP <https://diana-hep.org/>`_ study :cite:`faq-feickert-diana-fellowship-report`
demonstrated the graph structure and automatic differentiation abilities of machine learning
frameworks allowed them to be quite effective tools for statistical fits.
As the frameworks considered in this study (TensorFlow, PyTorch, MXNet) all provided
low-level Python APIs to the libraries this made Python an obvious choice for a common
high-level control language.
Given all these considerations, Python was chosen as the development language.

How did ``pyhf`` get started?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In 2017 Lukas Heinrich was discussing with colleauge Holger Schulz how it would be convienent
to share and produce statistical results from LHC experiements if they were able to be
created with tools that didn't require the large ``C++`` dependencies and tooling expertise as
:math:`\HiFa{}`.
Around the same time that Lukas began thinking on these ideas, Matthew Feickert was working on
`a DIANA/HEP fellowship <https://twitter.com/SMUPhysics/status/861584474638766080>`_ with
Kyle Cranmer (co-author of :math:`\HiFa{}`) to study if the graph structure and automatic
differentiation abilities of machine learning frameworks would allow them to be effective
tools for statistical fits.
Lukas would give helpful friendly advice on Matthew's project and one night [1]_ over dinner
in CERN's R1 cafeteria the two were discussing the idea of implementing :math:`\HiFa{}`
in Python using machine learning libraries to drive the computation.
Continuing the discussion in Lukas's office, Lukas showed Matthew that the core statistical
machinery could be implemented rather succinctly, and that night
`proceeded to do so <https://github.com/scikit-hep/pyhf/commit/fd32503fb760f070a4047cb867757458b1687599>`_
and |dubbed the project pyhf|_.

Matthew joined him on the project to begin development and by April 2018 Giordon Stark had
learned about the project and began making contributions, quickly becoming
`the third core developer <https://twitter.com/KyleCranmer/status/1052186117452259328>`_.
The first physics paper to use ``pyhf`` followed closely in October 2018
:cite:`faq-Heinrich:2018nip`, making Lukas and Holger's original conversations a reality.
``pyhf`` was founded on the ideas of open contributions and community software and continues
in that mission today as a `Scikit-HEP project <https://scikit-hep.org/>`_, with an open
invitation for community contributions and new developers.

Troubleshooting
---------------

- :code:`import torch` or :code:`import pyhf` causes a :code:`Segmentation fault (core dumped)`

    This is may be the result of a conflict with the NVIDIA drivers that you
    have installed on your machine.  Try uninstalling and completely removing
    all of them from your machine

    .. code-block:: console

        # On Ubuntu/Debian
        sudo apt-get purge nvidia*

    and then installing the latest versions.

Footnotes
~~~~~~~~~

.. [1]
   24 January, 2018

Bibliography
~~~~~~~~~~~~

.. bibliography:: bib/docs.bib bib/use_citations.bib
   :filter: docname in docnames
   :style: plain
   :keyprefix: faq-
   :labelprefix: faq-

.. |dubbed the project pyhf| replace:: dubbed the project ``pyhf``
.. _`dubbed the project pyhf`: https://twitter.com/lukasheinrich_/status/956809112674885632
.. pyhf documentation master file, created by
   sphinx-quickstart on Fri Feb  9 11:58:49 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:

   intro
   likelihood
   learn
   examples
   outreach
   installation
   development
   faq
   babel
   cli
   api
   citations
   governance/ROADMAP
   release-notes
   contributors

.. raw:: html

   <a class="github-fork-ribbon right-top fixed" href="https://github.com/scikit-hep/pyhf/" data-ribbon="View me on GitHub" title="View me on GitHub">View me on GitHub</a>


.. raw:: html

   <p id="dev-version"><strong>Warning:</strong> This is a development version. The latest stable version is at <a href="https://pyhf.readthedocs.io/">ReadTheDocs</a>.</p>

.. include:: ../README.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Translations
============
One key goal of ``pyhf`` is to provide seamless translations between other statistical frameworks and ``pyhf``.
This page details the various ways to translate from a tool you might already be using as part of an existing analysis to ``pyhf``.
Many of these solutions involve extracting out the ``HistFactory`` workspace and then running `pyhf xml2json <cli.html#pyhf-xml2json>`_ which provides a single JSON workspace that can be loaded directly into ``pyhf``.

HistFitter
----------

In order to go from ``HistFitter`` to ``pyhf``, the first step is to extract out the ``HistFactory`` workspaces. Assuming you have an existing configuration file, ``config.py``, you likely run an exclusion fit like so:

.. code:: bash

  HistFitter.py -f -D "before,after,corrMatrix" -F excl config.py

The name of output workspace files depends on four parameters you define in your ``config.py``:

- ``analysisName`` is from ``configMgr.analysisName``
- ``prefix`` is defined in ``configMgr.addFitConfig({prefix})``
- ``measurementName`` is the first measurement you define via ``fitConfig.addMeasurement(name={measurementName},...)``
- ``channelName`` are the names of channels you define via ``fitConfig.addChannel("cuts", [{channelName}], ...)``
- ``cachePath`` is where ``HistFitter`` stores the cached histograms, defined by ``configMgr.histCacheFile`` which defaults to ``data/{analysisName}.root``

To dump the HistFactory workspace, you will modify the above to skip the fit ``-f`` and plotting ``-D`` so you end up with

.. code:: bash

  HistFitter.py -wx -F excl config.py

The ``-w`` flag tells ``HistFitter`` to (re)create the ``HistFactory`` workspace stored in ``results/{analysisName}/{prefix}_combined_{measurementName}.root``.
The ``-x`` flag tells ``HistFitter`` to dump the XML files into ``config/{analysisName}/``, with the top-level file being ``{prefix}.xml`` and all other files being ``{prefix}_{channelName}_cuts.xml``.

Typically, ``prefix = 'FitConfig'`` and ``measurementName = 'NormalMeasurement'``. For example, if the following exists in your ``config.py``

.. code:: python

  from configManager import configMgr

  # ...
  configMgr.analysisName = "3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg"
  configMgr.histCacheFile = f"cache/{configMgr.analysisName:s}.root"
  # ...
  fitConfig = configMgr.addFitConfig("Excl")
  # ...
  channel = fitConfig.addChannel("cuts", ["SR_0L"], 1, 0.5, 1.5)
  # ...
  meas1 = fitConfig.addMeasurement(name="DefaultMeasurement", lumi=1.0, lumiErr=0.029)
  meas1.addPOI("mu_SIG1")
  # ...
  meas2 = fitConfig.addMeasurement(name="DefaultMeasurement", lumi=1.0, lumiErr=0.029)
  meas2.addPOI("mu_SIG2")

Then, you expect the following files to be made:

- ``config/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg/Excl.xml``
- ``config/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg/Excl_SR_0L_cuts.xml``
- ``cache/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg.root``
- ``results/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg/Excl_combined_DefaultMeasurement.root``

These are all the files you need in order to use `pyhf xml2json <cli.html#pyhf-xml2json>`_. At this point, you could run

.. code:: bash

    pyhf xml2json config/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg/Excl.xml

which will read all of the XML files and load the histogram data from the histogram cache.

The ``HistFactory`` workspace in ``results/`` contains all of the information necessary to rebuild the XML files again. For debugging purposes, the ``pyhf`` developers will often ask for your workspace file, which means ``results/3b_tag21.2.27-1_RW_ExpSyst_36100_multibin_bkg/Excl_combined_DefaultMeasurement.root``. If you want to generate the XML, you can open this file in ``ROOT`` and run ``DefaultMeasurement->PrintXML()`` which puts all of the XML files into the current directory you are in.


TRExFitter
----------

.. note::

    For more details on this section, please refer to the ATLAS-internal `TRExFitter documentation <https://trexfitter-docs.web.cern.ch/trexfitter-docs/advanced_topics/pyhf/>`_.

In order to go from ``TRExFitter`` to ``pyhf``, the good news is that the ``RooWorkspace`` files (``XML`` and ``ROOT``) are already made for you. For a given configuration which looks like

.. code:: yaml

    Job: "pyhf_example"
    Label: "..."

You can expect some files to be made after the ``n``/``h`` and ``w`` steps:

- ``pyhf_example/RooStats/pyhf_example.xml``
- ``pyhf_example/RooStats/pyhf_example_Signal_region.xml``
- ``pyhf_example/Histograms/pyhf_example_Signal_region_histos.root``

These are all the files you need in order to use `pyhf xml2json <cli.html#pyhf-xml2json>`_. At this point, you could run

.. code:: bash

    pyhf xml2json pyhf_example/RooStats/pyhf_example.xml

which will read all of the XML files and load the histogram data from the histogram cache.

.. warning::

    There are a few caveats one needs to be aware of with this conversion:

    - Uncorrelated shape systematics cannot be pruned, see Issue :issue:`662`.
    - Custom expressions for normalization factors cannot be used, see Issue :issue:`850`.
.. _sec:intro:

Introduction
============

Measurements in High Energy Physics (HEP) rely on determining the
compatibility of observed collision events with theoretical predictions.
The relationship between them is often formalised in a statistical *model*
:math:`f(\bm{x}|\fullset)` describing the probability of data
:math:`\bm{x}` given model parameters :math:`\fullset`. Given observed
data, the *likelihood* :math:`\mathcal{L}(\fullset)` then serves as the basis to test
hypotheses on the parameters \ :math:`\fullset`. For measurements based
on binned data (*histograms*), the :math:`\HiFa{}` family of statistical models has been widely used
in both Standard Model measurements :cite:`intro-HIGG-2013-02` as
well as searches for new
physics :cite:`intro-ATLAS-CONF-2018-041`. In this package, a
declarative, plain-text format for describing :math:`\HiFa{}`-based likelihoods is
presented that is targeted for reinterpretation and long-term
preservation in analysis data repositories such as
HEPData :cite:`intro-Maguire:2017ypu`.

HistFactory
-----------

Statistical models described using :math:`\HiFa{}` :cite:`intro-Cranmer:1456844`
center around the simultaneous measurement of disjoint binned
distributions (*channels*) observed as event counts :math:`\channelcounts`. For
each channel, the overall expected event rate [1]_ is the sum over a
number of physics processes (*samples*). The sample rates may be subject to
parametrised variations, both to express the effect of *free parameters*
:math:`\freeset` [2]_ and to account for systematic uncertainties as a
function of *constrained parameters* :math:`\constrset`. The degree to which the latter can cause
a deviation of the expected event rates from the nominal rates is
limited by *constraint terms*. In a frequentist framework these constraint terms can be
viewed as *auxiliary measurements* with additional global observable data :math:`\auxdata`, which
paired with the channel data :math:`\channelcounts` completes the
observation :math:`\bm{x} =
(\channelcounts,\auxdata)`. In addition to the partition of the full
parameter set into free and constrained parameters :math:`\fullset =
(\freeset,\constrset)`, a separate partition :math:`\fullset =
(\poiset,\nuisset)` will be useful in the context of hypothesis testing,
where a subset of the parameters are declared *parameters of interest* :math:`\poiset` and the
remaining ones as *nuisance parameters* :math:`\nuisset`.

.. math::
    :label: eqn:parameters_partitions

    f(\bm{x}|\fullset) = f(\bm{x}|\overbrace{\freeset}^{\llap{\text{free}}},\underbrace{\constrset}_{\llap{\text{constrained}}}) = f(\bm{x}|\overbrace{\poiset}^{\rlap{\text{parameters of interest}}},\underbrace{\nuisset}_{\rlap{\text{nuisance parameters}}})

Thus, the overall structure of a :math:`\HiFa{}` probability model is a product of the
analysis-specific model term describing the measurements of the channels
and the analysis-independent set of constraint terms:

.. math::
    :label: eqn:hifa_template

    f(\channelcounts, \auxdata \,|\,\freeset,\constrset) = \underbrace{\color{blue}{\prod_{c\in\mathrm{\,channels}} \prod_{b \in \mathrm{\,bins}_c}\textrm{Pois}\left(n_{cb} \,\middle|\, \nu_{cb}\left(\freeset,\constrset\right)\right)}}_{\substack{\text{Simultaneous measurement}\\%
      \text{of multiple channels}}} \underbrace{\color{red}{\prod_{\singleconstr \in \constrset} c_{\singleconstr}(a_{\singleconstr} |\, \singleconstr)}}_{\substack{\text{constraint terms}\\%
      \text{for }\unicode{x201C}\text{auxiliary measurements}\unicode{x201D}}},

where within a certain integrated luminosity we observe :math:`n_{cb}`
events given the expected rate of events
:math:`\nu_{cb}(\freeset,\constrset)` as a function of unconstrained
parameters :math:`\freeset` and constrained parameters
:math:`\constrset`. The latter has corresponding one-dimensional
constraint terms
:math:`c_\singleconstr(a_\singleconstr|\,\singleconstr)` with auxiliary
data :math:`a_\singleconstr` constraining the parameter
:math:`\singleconstr`. The event rates :math:`\nu_{cb}` are defined as

.. math::
    :label: eqn:sample_rates

    \nu_{cb}\left(\fullset\right) = \sum_{s\in\mathrm{\,samples}} \nu_{scb}\left(\freeset,\constrset\right) = \sum_{s\in\mathrm{\,samples}}\underbrace{\left(\prod_{\kappa\in\,\bm{\kappa}} \kappa_{scb}\left(\freeset,\constrset\right)\right)}_{\text{multiplicative modifiers}}\, \Bigg(\nu_{scb}^0\left(\freeset, \constrset\right) + \underbrace{\sum_{\Delta\in\bm{\Delta}} \Delta_{scb}\left(\freeset,\constrset\right)}_{\text{additive modifiers}}\Bigg)\,.

The total rates are the sum over sample rates :math:`\nu_{csb}`, each
determined from a *nominal rate* :math:`\nu_{scb}^0` and a set of multiplicative and
additive denoted *rate modifiers* :math:`\bm{\kappa}(\fullset)` and
:math:`\bm{\Delta}(\fullset)`. These modifiers are functions of (usually
a single) model parameters. Starting from constant nominal rates, one
can derive the per-bin event rate modification by iterating over all
sample rate modifications as shown in :eq:`eqn:sample_rates`.

As summarised in :ref:`tab:modifiers_and_constraints`, rate modifications
are defined in :math:`\HiFa{}` for bin :math:`b`, sample :math:`s`, channel
:math:`c`.  Each modifier is represented by a parameter :math:`\phi \in
\{\gamma, \alpha, \lambda, \mu\}`.  By convention bin-wise parameters are
denoted with :math:`\gamma` and interpolation parameters with :math:`\alpha`.
The luminosity :math:`\lambda` and scale factors :math:`\mu` affect all bins
equally.  For constrained modifiers, the implied constraint term is given as
well as the necessary input data required to construct it.  :math:`\sigma_b`
corresponds to the relative uncertainty of the event rate, whereas
:math:`\delta_b` is the event rate uncertainty of the sample relative to the
total event rate :math:`\nu_b = \sum_s \nu^0_{sb}`.

Modifiers implementing uncertainties are paired with
a corresponding default constraint term on the parameter limiting the
rate modification. The available modifiers may affect only the total
number of expected events of a sample within a given channel, i.e. only
change its normalisation, while holding the distribution of events
across the bins of a channel, i.e. its “shape”, invariant.
Alternatively, modifiers may change the sample shapes. Here :math:`\HiFa{}` supports
correlated an uncorrelated bin-by-bin shape modifications. In the
former, a single nuisance parameter affects the expected sample rates
within the bins of a given channel, while the latter introduces one
nuisance parameter for each bin, each with their own constraint term.
For the correlated shape and normalisation uncertainties, :math:`\HiFa{}` makes use of
interpolating functions, :math:`f_p` and :math:`g_p`, constructed from a
small number of evaluations of the expected rate at fixed values of the
parameter :math:`\alpha` [3]_. For the remaining modifiers, the
parameter directly affects the rate.

.. _tab:modifiers_and_constraints:

.. table:: Modifiers and Constraints

    ==================== ============================================================================================================= ===================================================================================================== ================================
    Description          Modification                                                                                                  Constraint Term :math:`c_\singleconstr`                                                               Input
    ==================== ============================================================================================================= ===================================================================================================== ================================
    Uncorrelated Shape   :math:`\kappa_{scb}(\gamma_b) = \gamma_b`                                                                     :math:`\prod_b \mathrm{Pois}\left(r_b = \sigma_b^{-2}\middle|\,\rho_b = \sigma_b^{-2}\gamma_b\right)` :math:`\sigma_{b}`
    Correlated Shape     :math:`\Delta_{scb}(\alpha) = f_p\left(\alpha\middle|\,\Delta_{scb,\alpha=-1},\Delta_{scb,\alpha = 1}\right)` :math:`\displaystyle\mathrm{Gaus}\left(a = 0\middle|\,\alpha,\sigma = 1\right)`                       :math:`\Delta_{scb,\alpha=\pm1}`
    Normalisation Unc.   :math:`\kappa_{scb}(\alpha) = g_p\left(\alpha\middle|\,\kappa_{scb,\alpha=-1},\kappa_{scb,\alpha=1}\right)`   :math:`\displaystyle\mathrm{Gaus}\left(a = 0\middle|\,\alpha,\sigma = 1\right)`                       :math:`\kappa_{scb,\alpha=\pm1}`
    MC Stat. Uncertainty :math:`\kappa_{scb}(\gamma_b) = \gamma_b`                                                                     :math:`\prod_b \mathrm{Gaus}\left(a_{\gamma_b} = 1\middle|\,\gamma_b,\delta_b\right)`                 :math:`\delta_b^2 = \sum_s\delta^2_{sb}`
    Luminosity           :math:`\kappa_{scb}(\lambda) = \lambda`                                                                       :math:`\displaystyle\mathrm{Gaus}\left(l = \lambda_0\middle|\,\lambda,\sigma_\lambda\right)`          :math:`\lambda_0,\sigma_\lambda`
    Normalisation        :math:`\kappa_{scb}(\mu_b) = \mu_b`
    Data-driven Shape    :math:`\kappa_{scb}(\gamma_b) = \gamma_b`
    ==================== ============================================================================================================= ===================================================================================================== ================================

Given the likelihood :math:`\mathcal{L}(\fullset)`, constructed from
observed data in all channels and the implied auxiliary data, *measurements* in the
form of point and interval estimates can be defined. The majority of the
parameters are *nuisance parameters* — parameters that are not the main target of the
measurement but are necessary to correctly model the data. A small
subset of the unconstrained parameters may be declared as *parameters of interest* for which
measurements hypothesis tests are performed, e.g. profile likelihood
methods :cite:`intro-Cowan:2010js`. The :ref:`tab:symbol_summary` table provides a summary of all the
notation introduced in this documentation.

.. _tab:symbol_summary:

.. table:: Symbol Notation

    =================================================================== ===============================================================
    Symbol                                                              Name
    =================================================================== ===============================================================
    :math:`f(\bm{x} | \fullset)`                                        model
    :math:`\mathcal{L}(\fullset)`                                       likelihood
    :math:`\bm{x} = \{\channelcounts, \auxdata\}`                       full dataset (including auxiliary data)
    :math:`\channelcounts`                                              channel data (or event counts)
    :math:`\auxdata`                                                    auxiliary data
    :math:`\nu(\fullset)`                                               calculated event rates
    :math:`\fullset = \{\freeset, \constrset\} = \{\poiset, \nuisset\}` all parameters
    :math:`\freeset`                                                    free parameters
    :math:`\constrset`                                                  constrained parameters
    :math:`\poiset`                                                     parameters of interest
    :math:`\nuisset`                                                    nuisance parameters
    :math:`\bm{\kappa}(\fullset)`                                       multiplicative rate modifier
    :math:`\bm{\Delta}(\fullset)`                                       additive rate modifier
    :math:`c_\singleconstr(a_\singleconstr | \singleconstr)`            constraint term for constrained parameter :math:`\singleconstr`
    :math:`\sigma_\singleconstr`                                        relative uncertainty in the constrained parameter
    =================================================================== ===============================================================

Declarative Formats
-------------------

While flexible enough to describe a wide range of LHC measurements, the
design of the :math:`\HiFa{}` specification is sufficiently simple to admit a *declarative format* that fully
encodes the statistical model of the analysis. This format defines the
channels, all associated samples, their parameterised rate modifiers and
implied constraint terms as well as the measurements. Additionally, the
format represents the mathematical model, leaving the implementation of
the likelihood minimisation to be analysis-dependent and/or
language-dependent. Originally XML was chosen as a specification
language to define the structure of the model while introducing a
dependence on :math:`\Root{}` to encode the nominal rates and required input data of the
constraint terms :cite:`intro-Cranmer:1456844`. Using this
specification, a model can be constructed and evaluated within the
:math:`\RooFit{}` framework.

This package introduces an updated form of the specification based on
the ubiquitous plain-text JSON format and its schema-language *JSON Schema*.
Described in more detail in :ref:`sec:likelihood`, this schema fully specifies both structure
and necessary constrained data in a single document and thus is
implementation independent.

Additional Material
-------------------

Footnotes
~~~~~~~~~

.. [1]
   Here rate refers to the number of events expected to be observed
   within a given data-taking interval defined through its integrated
   luminosity. It often appears as the input parameter to the Poisson
   distribution, hence the name “rate”.

.. [2]
   These *free parameters* frequently include the of a given process, i.e. its cross-section
   normalised to a particular reference cross-section such as that expected
   from the Standard Model or a given BSM scenario.

.. [3]
   This is usually constructed from the nominal rate and measurements of the
   event rate at :math:`\alpha=\pm1`, where the value of the modifier at
   :math:`\alpha=\pm1` must be provided and the value at :math:`\alpha=0`
   corresponds to the corresponding identity operation of the modifier, i.e.
   :math:`f_{p}(\alpha=0) = 0` and :math:`g_{p}(\alpha = 0)=1` for additive and
   multiplicative modifiers respectively. See Section 4.1
   in :cite:`intro-Cranmer:1456844`.

Bibliography
~~~~~~~~~~~~

.. bibliography:: bib/docs.bib
   :filter: docname in docnames
   :style: plain
   :keyprefix: intro-
   :labelprefix: intro-
Python API
==========

Top-Level
---------

.. currentmodule:: pyhf

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   tensorlib
   optimizer
   get_backend
   set_backend
   readxml
   writexml
   compat

Probability Distribution Functions (PDFs)
-----------------------------------------

.. currentmodule:: pyhf.probability

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   Normal
   Poisson
   Independent
   Simultaneous

Making Models from PDFs
-----------------------

.. currentmodule:: pyhf

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   ~pdf.Model
   ~pdf._ModelConfig
   ~workspace.Workspace
   ~patchset.PatchSet
   ~patchset.Patch
   simplemodels.uncorrelated_background
   simplemodels.correlated_background

Backends
--------

The computational backends that :code:`pyhf` provides interfacing for the vector-based calculations.

.. currentmodule:: pyhf.tensor

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   numpy_backend.numpy_backend
   pytorch_backend.pytorch_backend
   tensorflow_backend.tensorflow_backend
   jax_backend.jax_backend

Optimizers
----------

.. currentmodule:: pyhf.optimize

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   mixins.OptimizerMixin
   opt_scipy.scipy_optimizer
   opt_minuit.minuit_optimizer

Modifiers
---------

.. currentmodule:: pyhf.modifiers

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   histosys
   normfactor
   normsys
   shapefactor
   shapesys
   staterror

Interpolators
-------------

.. currentmodule:: pyhf.interpolators

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   code0
   code1
   code2
   code4
   code4p

Inference
---------

.. currentmodule:: pyhf.infer


Test Statistics
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   test_statistics.q0
   test_statistics.qmu
   test_statistics.qmu_tilde
   test_statistics.tmu
   test_statistics.tmu_tilde
   utils.get_test_stat

Calculators
~~~~~~~~~~~

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   calculators.generate_asimov_data
   calculators.HypoTestFitResults
   calculators.AsymptoticTestStatDistribution
   calculators.EmpiricalDistribution
   calculators.AsymptoticCalculator
   calculators.ToyCalculator
   utils.create_calculator

Fits and Tests
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   mle.twice_nll
   mle.fit
   mle.fixed_poi_fit
   hypotest
   intervals.upperlimit
   utils.all_pois_floating

Exceptions
----------

Various exceptions, apart from standard python exceptions, that are raised from using the :code:`pyhf` API.

.. currentmodule:: pyhf.exceptions

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   InvalidMeasurement
   InvalidNameReuse
   InvalidSpecification
   InvalidPatchSet
   InvalidPatchLookup
   PatchSetVerificationError
   InvalidWorkspaceOperation
   InvalidModel
   InvalidModifier
   InvalidInterpCode
   ImportBackendError
   InvalidBackend
   InvalidOptimizer
   InvalidPdfParameters
   InvalidPdfData

Utilities
---------

.. currentmodule:: pyhf.utils

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   load_schema
   validate
   options_from_eqdelimstring
   digest
   citation

Contrib
-------

.. currentmodule:: pyhf.contrib

.. autosummary::
   :toctree: _generated/
   :nosignatures:

   viz.brazil
   utils.download
{{ fullname | escape | underline }}

.. rubric:: Description

.. automodule:: {{ fullname }}

.. currentmodule:: {{ fullname }}

{% if classes %}
.. rubric:: Classes

.. autosummary::
    :toctree: .
    {% for class in classes %}
    {{ class }}
    {% endfor %}

{% endif %}

{% if functions %}
.. rubric:: Functions

.. autosummary::
    :toctree: .
    {% for function in functions %}
    {{ function }}
    {% endfor %}

{% endif %}
{{ name | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :show-inheritance:

   .. automethod:: __init__

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   {% for item in attributes %}
   .. autoattribute:: {{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   {% for item in members %}
   {% if item not in attributes and item not in inherited_members and not item.startswith('__') %}
   .. automethod:: {{ name }}.{{ item }}
   {% endif %}
   {%- endfor %}

   {% endif %}
   {% endblock %}
|release v0.5.4|_
=================

This is a patch release from ``v0.5.3`` → ``v0.5.4``.

Fixes
-----

* Require ``uproot3`` instead of ``uproot`` ``v3.X`` releases to avoid conflicts when
  ``uproot4`` is installed in an environment with ``uproot`` ``v3.X`` installed and
  namespace conflicts with ``uproot-methods``.
  Adoption of ``uproot3`` in ``v0.5.4`` will ensure ``v0.5.4`` works far into the future
  if XML and ROOT I/O through uproot is required.

**Example:**

Without the ``v0.5.4`` patch release there is a regression in using ``uproot`` ``v3.X``
and ``uproot4`` in the same environment (which was swiftly identified and patched by the
fantastic ``uproot`` team)

.. code-block:: shell

   $ python -m pip install "pyhf[xmlio]<0.5.4"
   $ python -m pip list | grep "pyhf\|uproot"
   pyhf           0.5.3
   uproot         3.13.1
   uproot-methods 0.8.0
   $ python -m pip install uproot4
   $ python -m pip list | grep "pyhf\|uproot"
   pyhf           0.5.3
   uproot         4.0.0
   uproot-methods 0.8.0
   uproot4        4.0.0

this is resolved in ``v0.5.4`` with the requirement of ``uproot3``

.. code-block:: shell

   $ python -m pip install "pyhf[xmlio]>=0.5.4"
   $ python -m pip list | grep "pyhf\|uproot"
   pyhf            0.5.4
   uproot3         3.14.1
   uproot3-methods 0.10.0
   $ python -m pip install uproot4 # or uproot
   $ python -m pip list | grep "pyhf\|uproot"
   pyhf            0.5.4
   uproot          4.0.0
   uproot3         3.14.1
   uproot3-methods 0.10.0
   uproot4         4.0.0

.. |release v0.5.4| replace:: ``v0.5.4``
.. _`release v0.5.4`: https://github.com/scikit-hep/pyhf/releases/tag/v0.5.4
|release v0.6.2|_
=================

This is a patch release from ``v0.6.1`` → ``v0.6.2``.

Important Notes
---------------

* The :func:`pyhf.simplemodels.hepdata_like` API has been deprecated in favor of
  :func:`pyhf.simplemodels.uncorrelated_background`.
  The :func:`pyhf.simplemodels.hepdata_like` API will be removed in ``pyhf`` ``v0.7.0``.
  (PR :pr:`1438`)
* There is a small breaking API change for :func:`pyhf.contrib.viz.brazil.plot_results`.
  See the Python API changes section for more information.
* The :class:`pyhf.patchset.PatchSet` schema now allows string types for patch values in patchsets.
  (PR :pr:`1488`)
* Only lower bounds on core dependencies are now set.
  This allows for greater developer freedom and reduces the risk of breaking
  user's applications by unnecessarily constraining libraries.
  This also means that users will be responsible for ensuring that their
  installed dependencies do not conflict with or break ``pyhf``.
  c.f. Hynek Schlawack's blog post `Semantic Versioning Will Not Save You
  <https://hynek.me/articles/semver-will-not-save-you/>`_ for more in-depth coverage
  on this topic.
  For most users nothing should change.
  This mainly affects developers of other libraries in which ``pyhf`` is a dependency.
  (PR :pr:`1382`)
* Calling ``dir()`` on any ``pyhf`` module or trying to tab complete an API will
  now provide a more helpfully restricted view of the available APIs.
  This should help provide better exploration of the ``pyhf`` API.
  (PR :pr:`1403`)
* Docker images of releases are now published to both `Docker Hub
  <https://hub.docker.com/r/pyhf/pyhf/tags>`_ and to the `GitHub Container
  Registry <https://github.com/scikit-hep/pyhf/pkgs/container/pyhf>`_.
  (PR :pr:`1444`)
* CUDA enabled Docker images are now available for release ``v0.6.1`` and later
  on `Docker Hub <https://hub.docker.com/r/pyhf/cuda>`__ and the `GitHub
  Container Registry <https://github.com/pyhf/cuda-images/pkgs/container/cuda-images>`__.
  Visit `github.com/pyhf/cuda-images <https://github.com/pyhf/cuda-images>`_ for more
  information.

Fixes
-----

* Allow for precision to be properly set for the tensorlib ``ones`` and ``zeros``
  methods through a ``dtype`` argument.
  This allows for precision to be properly set through the :func:`pyhf.set_backend`
  ``precision`` argument.
  (PR :pr:`1369`)
* The default precision for all backends is now ``64b``.
  (PR :pr:`1400`)
* Add check to ensure that POIs are not fixed during a fit.
  (PR :pr:`1409`)
* Parameter name strings are now normalized to remove trailing spaces.
  (PR :pr:`1436`)
* The logging level is now not automatically set in :class:`pyhf.contrib.utils`.
  (PR :pr:`1460`)

Features
--------

Python API
~~~~~~~~~~

* The :func:`pyhf.simplemodels.hepdata_like` API has been deprecated in favor of
  :func:`pyhf.simplemodels.uncorrelated_background`.
  The :func:`pyhf.simplemodels.hepdata_like` API will be removed in ``pyhf`` ``v0.7.0``.
  (PR :pr:`1438`)
* The :func:`pyhf.simplemodels.correlated_background` API has been added to
  provide an example model with a single channel with a correlated background
  uncertainty.
  (PR :pr:`1435`)
* Add CLs component plotting kwargs to :func:`pyhf.contrib.viz.brazil.plot_results`.
  This allows CLs+b and CLb components of the CLs ratio to be plotted as well.
  To be more consistent with the ``matplotlib`` API,
  :func:`pyhf.contrib.viz.brazil.plot_results` now returns a lists of the artists
  drawn on the axis and moves the ``ax`` arguments to the to the last argument.
  (PR :pr:`1377`)
* The ``pyhf.compat`` module has been added to aid in translating to and from ROOT
  names.
  (PR :pr:`1439`)

CLI API
~~~~~~~

* The CLI API now supports a ``patchset inspect`` API to list the individual
  patches in a ``PatchSet``.
  (PR :pr:`1412`)

.. code-block:: shell

  pyhf patchset inspect [OPTIONS] [PATCHSET]

Contributors
------------

``v0.6.2`` benefited from contributions from:

* Alexander Held

.. |release v0.6.2| replace:: ``v0.6.2``
.. _`release v0.6.2`: https://github.com/scikit-hep/pyhf/releases/tag/v0.6.2
|release v0.6.0|_
=================

This is a minor release from ``v0.5.4`` → ``v0.6.0``.

Important Notes
---------------

* Please note this release has **API breaking changes** and carefully read these
  notes while updating your code to the ``v0.6.0`` API.
  Perhaps most relevant is the changes to the :func:`pyhf.infer.hypotest` API, which now
  uses a ``calctype`` argument to differentiate between using an asymptotic calculator
  or a toy calculator, and a ``test_stat`` kwarg to specify which test statistic
  the calculator should use, with ``'qtilde'``, corresponding to
  :func:`pyhf.infer.test_statistics.qmu_tilde`, now the default option.
  It also relies more heavily on using kwargs to pass options through to the optimizer.
* Following the recommendations of |NEP 29|_ ``pyhf`` ``v0.6.0`` drops support for
  Python 3.6.
  |PEP 494|_ also notes that Python 3.6 will be end of life in December 2021, so
  ``pyhf`` is moving forward with a minimum required runtime of Python 3.7.
* Support for the discovery test statistic, :math:`q_{0}`, has now been added through
  the :func:`pyhf.infer.test_statistics.q0` API.
* Support for pseudoexperiments (toys) has been added through the
  :func:`pyhf.infer.calculators.ToyCalculator` API.
  Please see the corresponding `example notebook`_ for more detailed exploration
  of the API.
* The ``minuit`` extra, ``python -m pip install pyhf[minuit]``, now uses and requires
  the |iminuit docs|_ ``v2.X`` release series and API.
  Note that ``iminuit`` ``v2.X`` can result in slight differences in minimization
  results from ``iminuit`` ``v1.X``.
* The documentation will now be versioned with releases on ReadTheDocs.
  Please use `pyhf.readthedocs.io`_ to access the documentation for the latest
  stable release of ``pyhf``.
* ``pyhf`` is transtioning away from Stack Overflow to `GitHub Discussions`_ for
  resolving user questions not covered in the documentation.
  Please check the `GitHub Discussions`_ page to search for discussions addressing
  your questions and to open up a new discussion if your question is not covered.
* ``pyhf`` has published a paper in the Journal of Open Source Software. |JOSS DOI|
  Please make sure to include the paper reference in all citations of ``pyhf``, as
  documented in the `Use and Citations`_ section of the documentation.

Fixes
-----

* Fix bug where all extras triggered warning for installation of the ``contrib`` extra.
* ``float``-like values are used in division for :func:`pyhf.writexml`.
* ``Model.spec`` now supports building new models from existing models.
* :math:`p`-values are now reported based on their quantiles, instead of interpolating
  test statistics and converting to :math:`p`-values.
* Namespace collisions between ``uproot3`` and ``uproot``/``uproot4`` have been fixed
  for the ``xmlio`` extra.
* The ``normsys`` modifier now uses the :mod:`pyhf.interpolators.code4` interpolation
  method by default.
* The ``histosys`` modifier now uses the :mod:`pyhf.interpolators.code4p` interpolation
  method by default.

Features
--------

Python API
~~~~~~~~~~

* The ``tensorlib`` API now supports a ``tensorlib.to_numpy`` and
  ``tensorlib.ravel`` API.
* The :func:`pyhf.infer.calculators.ToyCalculator` API has been added to support
  pseudoexperiments (toys).
* The empirical test statistic distribution API has been added to help support the
  ``ToyCalculator`` API.
* Add a ``tolerance`` kwarg to the optimizer API to set a ``float`` value as a
  tolerance for termination of the fit.
* The :func:`pyhf.optimize.opt_minuit.minuit_optimizer` optimizer now can return
  correlations of the fitted parameters through use of the ``return_correlation``
  Boolean kwarg.
* Add the ``pyhf.utils.citation`` API to get a ``str`` of the preferred BibTeX entry
  for citation of the version of ``pyhf`` installed.
  See the example for the CLI API for more information.
* The :func:`pyhf.infer.hypotest` API now uses a ``calctype`` argument to differentiate
  between using an asymptotic calculator or a toy calculator, and a ``test_stat`` kwarg
  to specify which test statistic to use.
  It also relies more heavily on using kwargs to pass options through to the optimizer.
* The default ``test_stat`` kwarg for :func:`pyhf.infer.hypotest` and the calculator
  APIs is ``'qtilde'``, which corresponds to the alternative test statistic
  :func:`pyhf.infer.test_statistics.qmu_tilde`.
* The return type of :math:`p`-value like functions is now a 0-dimensional ``tensor``
  (with shape ``()``) instead of a ``float``.
  This is required to support end-to-end automatic differentiation in future releases.

CLI API
~~~~~~~

* The CLI API now supports a ``--citation`` or ``--cite`` option to print the
  preferred BibTeX entry for citation of the version of ``pyhf`` installed.

.. code-block:: shell

   $ pyhf --citation
   @software{pyhf,
     author = {Lukas Heinrich and Matthew Feickert and Giordon Stark},
     title = "{pyhf: v0.6.0}",
     version = {0.6.0},
     doi = {10.5281/zenodo.1169739},
     url = {https://doi.org/10.5281/zenodo.1169739},
     note = {https://github.com/scikit-hep/pyhf/releases/tag/v0.6.0}
   }

   @article{pyhf_joss,
     doi = {10.21105/joss.02823},
     url = {https://doi.org/10.21105/joss.02823},
     year = {2021},
     publisher = {The Open Journal},
     volume = {6},
     number = {58},
     pages = {2823},
     author = {Lukas Heinrich and Matthew Feickert and Giordon Stark and Kyle Cranmer},
     title = {pyhf: pure-Python implementation of HistFactory statistical models},
     journal = {Journal of Open Source Software}
   }

Contributors
------------

``v0.6.0`` benefited from contributions from:

* Alexander Held
* Marco Gorelli
* Pradyumna Rahul K
* Eric Schanet
* Henry Schreiner

.. |release v0.6.0| replace:: ``v0.6.0``
.. _`release v0.6.0`: https://github.com/scikit-hep/pyhf/releases/tag/v0.6.0

.. |NEP 29| replace:: NEP 29 — Recommend Python and NumPy version support as a community policy standard
.. _`NEP 29`: https://numpy.org/neps/nep-0029-deprecation_policy.html

.. |PEP 494| replace:: PEP 494 -- Python 3.6 Release Schedule
.. _`PEP 494`: https://www.python.org/dev/peps/pep-0494/

.. _`example notebook`: https://pyhf.readthedocs.io/en/latest/examples/notebooks/toys.html

.. |iminuit docs| replace:: ``iminuit``
.. _`iminuit docs`: https://iminuit.readthedocs.io/

.. _`pyhf.readthedocs.io`: https://pyhf.readthedocs.io/

.. _`GitHub Discussions`: https://github.com/scikit-hep/pyhf/discussions

.. |JOSS DOI| image:: https://joss.theoj.org/papers/10.21105/joss.02823/status.svg
   :target: https://doi.org/10.21105/joss.02823

.. _`Use and Citations`: https://pyhf.readthedocs.io/en/latest/citations.html
|release v0.6.3|_
=================

This is a patch release from ``v0.6.2`` → ``v0.6.3``.

Important Notes
---------------

* With the addition of writing ROOT files in |uproot v4.1.0 release|_ the
  ``xmlio`` extra no longer requires ``uproot3`` and all dependencies on
  ``uproot3`` and ``uproot3-methods`` have been dropped.
  (PR :pr:`1567`)
  ``uproot4`` additionally brings large speedups to writing, which results in an
  order of magnitude faster conversion time for most workspace conversions from
  JSON back to XML + ROOT with ``pyhf json2xml``.
* All backends are now fully compatible and tested with
  `Python 3.9 <https://www.python.org/dev/peps/pep-0596/>`_.
  (PR :pr:`1574`)
* The TensorFlow backend now supports compatibility with TensorFlow ``v2.2.1``
  and later and TensorFlow Probability ``v0.10.1`` and later.
  (PR :pr:`1001`)
* The :func:`pyhf.workspace.Workspace.data` ``with_aux`` keyword arg has been
  renamed to ``include_auxdata`` to improve API consistency.
  (PR :pr:`1562`)

.. |uproot v4.1.0 release| replace:: ``uproot`` ``v4.1.0``
.. _`uproot v4.1.0 release`: https://github.com/scikit-hep/uproot4/releases/tag/4.1.0

Fixes
-----

* The weakref bug with Click ``v8.0+`` was resolved.
  ``pyhf`` is now fully compatible with Click ``v7`` and ``v8`` releases.
  (PR :pr:`1530`)

Features
--------

Python API
~~~~~~~~~~

* Model parameter names are now propagated to optimizers through addition of the
  :func:`pyhf.pdf._ModelConfig.par_names` API.
  :func:`pyhf.pdf._ModelConfig.par_names` also handles non-scalar modifiers with
  1 parameter.
  (PRs :pr:`1536`, :pr:`1560`)

  .. code:: pycon

      >>> import pyhf
      >>> model = pyhf.simplemodels.uncorrelated_background(
      ...     signal=[12.0, 11.0], bkg=[50.0, 52.0], bkg_uncertainty=[3.0, 7.0]
      ... )
      >>> model.config.parameters
      ['mu', 'uncorr_bkguncrt']
      >>> model.config.npars
      3
      >>> model.config.par_names()
      ['mu', 'uncorr_bkguncrt[0]', 'uncorr_bkguncrt[1]']

* The :class:`pyhf.pdf._ModelConfig` ``channel_nbins`` dict is now sorted by
  keys to match the order of the ``channels`` list.
  (PR :pr:`1546`)

* The :func:`pyhf.workspace.Workspace.data` ``with_aux`` keyword arg has been
  renamed to ``include_auxdata`` to improve API consistency.
  (PR :pr:`1562`)

.. |release v0.6.3| replace:: ``v0.6.3``
.. _`release v0.6.3`: https://github.com/scikit-hep/pyhf/releases/tag/v0.6.3
|release v0.5.3|_
=================

This is a patch release from ``v0.5.2`` → ``v0.5.3``.

Fixes
-----

* Workspaces are now immutable
* ShapeFactor support added to XML reading and writing
* An error is raised if a fit initialization parameter is outside of its bounds
  (preventing hypotest with POI outside of bounds)

Features
--------

Python API
~~~~~~~~~~

* Inverting hypothesis tests to get upper limits now has an API with
  ``pyhf.infer.intervals.upperlimit``
* Building workspaces from a model and data added with ``pyhf.workspace.build``

CLI API
~~~~~~~

* Added CLI API for ``pyhf.infer.fit``: ``pyhf fit``
* pyhf combine now allows for merging channels: ``pyhf combine --merge-channels --join <join option>``
* Added utility to download archived pyhf pallets (workspaces + patchsets) to contrib module: ``pyhf contrib download``

Contributors
------------

``v0.5.3`` benefited from contributions from:

* Karthikeyan Singaravelan

.. |release v0.5.3| replace:: ``v0.5.3``
.. _`release v0.5.3`: https://github.com/scikit-hep/pyhf/releases/tag/v0.5.3
|release v0.6.1|_
=================

This is a patch release from ``v0.6.0`` → ``v0.6.1``.

Important Notes
---------------

* As a result of changes to the default behavior of ``torch.distributions`` in
  PyTorch ``v1.8.0``, accommodating changes have been made in the underlying
  implementations for :func:`pyhf.tensor.pytorch_backend.pytorch_backend`.
  These changes require a new lower bound of ``torch`` ``v1.8.0`` for use of the
  PyTorch backend.

Fixes
-----

* In the PyTorch backend the ``validate_args`` kwarg is used with
  ``torch.distributions`` to ensure a continuous approximation of the Poisson
  distribution in ``torch`` ``v1.8.0+``.

Features
--------

Python API
~~~~~~~~~~

* The ``solver_options`` kwarg can be passed to the
  :func:`pyhf.optimize.opt_scipy.scipy_optimizer` optimizer for additional
  configuration of the minimization.
  See :func:`scipy.optimize.show_options` for additional options of optimization
  solvers.
* The ``torch`` API is now used to provide the implementations of the ``ravel``,
  ``tile``, and ``outer`` tensorlib methods for the PyTorch backend.

.. |release v0.6.1| replace:: ``v0.6.1``
.. _`release v0.6.1`: https://github.com/scikit-hep/pyhf/releases/tag/v0.6.1
Roadmap (2019-2020)
===================

This is the pyhf 2019 into 2020 Roadmap (Issue
`#561 <https://github.com/scikit-hep/pyhf/issues/561>`__).

Overview and Goals
------------------

We will follow loosely Seibert’s `Heirarchy of
Needs <https://twitter.com/FRoscheck/status/1159158552298229763>`__

|Seibert Hierarchy of Needs SciPy 2019| (`Stan
Seibert <https://github.com/seibert>`__, SciPy 2019)

As a general overview that will include:

-  Improvements to docs

   -  Add lots of examples
   -  Add at least 5 well documented case studies

-  Issue cleanup
-  Adding core feature support
-  "pyhf evolution": integration with columnar data analysis systems
-  GPU support and testing
-  Publications

   -  Submit pyhf to JOSS
   -  Submit pyhf to pyOpenSci
   -  Start pyhf paper in 2020

-  Align with IRIS-HEP Analysis Systems NSF milestones

Time scale
----------

The roadmap will be executed over mostly Quarter 3 of 2019 through
Quarter 1 of 2020, with some projects continuing into Quarter 2 of 2020

-  2019-Q3
-  2019-Q4
-  2020-Q1
-  (2020-Q2)

Roadmap
-------

1. **Documentation and Deployment**

   -  |uncheck| Add docstrings to all functions and classes (Issues #38, #349)
      [2019-Q3]
   -  |uncheck| `Greatly revise and expand
      examples <https://github.com/scikit-hep/pyhf/issues?q=is%3Aopen+is%3Aissue+label%3Adocs>`__
      (Issues #168, #202, #212, #325, #342, #349, #367) [2019-Q3 →
      2019-Q4]

      -  |uncheck| Add small case studies with published sbottom likelihood from
         HEPData

   -  |check| Move to `scikit-hep <https://github.com/scikit-hep>`__ GitHub
      organization [2019-Q3]
   -  |uncheck| Develop a release schedule/criteria [2019-Q4]
   -  |check| Automate deployment with [STRIKEOUT:Azure pipeline (talk with
      Henry Schreiner) (Issue #517)] GitHub Actions (Issue #508)
      [2019-Q3]
   -  |uncheck| Finalize logo and add it to website (Issue #453) [2019-Q3 →
      2019-Q4]
   -  |check| Write submission to `JOSS <https://joss.theoj.org/>`__ (Issue
      #502) and write submission to
      `pyOpenSci <https://www.pyopensci.org/>`__ [2019-Q4 → 2020-Q2]
   -  |check| Contribute to `IRIS-HEP Analysis Systems
      Milestones <https://docs.google.com/spreadsheets/d/1VKpHlQWXu_p8AUv5E5H_BzqF_i7hh2Z-Id0XPwNHu8o/edit#gid=1864915304>`__
      "`Initial roadmap for ecosystem
      coherency <https://github.com/iris-hep/project-milestones/issues/8>`__"
      and "`Initial roadmap for high-level cyberinfrastructure
      components of analysis
      system <https://github.com/iris-hep/project-milestones/issues/11>`__"
      [2019-Q4 → 2020-Q2]

2. **Revision and Maintenance**

   -  |check| Add tests using HEPData published sbottom likelihoods (Issue
      #518) [2019-Q3]
   -  |check| Add CI with GitHub Actions and Azure Pipelines (PR #527, Issue
      #517) [2019-Q3]
   -  |check| Investigate rewrite of pytest fixtures to use modern pytest
      (Issue #370) [2019-Q3 → 2019-Q4]
   -  |check| Factorize out the statistical fitting portion into
      ``pyhf.infer`` (PR #531) [2019-Q3 → 2019-Q4]
   -  |uncheck| Bug squashing at large [2019-Q3 → 2020-Q2]

      -  |uncheck| Unexpected use cases (Issues #324, #325, #529)
      -  |uncheck| Computational edge cases (Issues #332, #445)

   -  |uncheck| Make sure that all backends reproduce sbottom results [2019-Q4 →
      2020-Q2]

3. **Development**

   -  |check| Batch support (PR #503) [2019-Q3]
   -  |check| Add ParamViewer support (PR #519) [2019-Q3]
   -  |check| Add setting of NPs constant/fixed (PR #653) [2019-Q3]
   -  |check| Implement pdf as subclass of distributions (PR #551) [2019-Q3]
   -  |check| Add sampling with toys (PR #558) [2019-Q3]
   -  |uncheck| Make general modeling choices (e.g., Issue #293) [2019-Q4 →
      2020-Q1]
   -  |check| Add "discovery" test stats (p0) (PR #520) [2019-Q4 → 2020-Q1]
   -  |uncheck| Add better Model creation [2019-Q4 → 2020-Q1]
   -  |check| Add background model support (Issues #514, #946) [2019-Q4 → 2020-Q1]
   -  |check| Develop interface for the optimizers similar to tensor/backend
      (Issue #754, PR #951) [2019-Q4 → 2020-Q1]
   -  |check| Migrate to TensorFlow v2.0 (PR #541) [2019-Q4]
   -  |check| Drop Python 2.7 support at end of 2019 (Issue #469) [2019-Q4
      (last week of December 2019)]
   -  |uncheck| Finalize public API [2020-Q1]
   -  |uncheck| Integrate pyfitcore/Statisfactory API [2020-Q1]

4. **Research**

   -  |uncheck| Add pyfitcore/Statisfactory integrations (Issue #344, `zfit
      Issue 120 <https://github.com/zfit/zfit/issues/120>`__) [2019-Q4]
   -  |uncheck| Hardware acceleration scaling studies (Issues #93, #301)
      [2019-Q4 → 2020-Q1]
   -  |uncheck| Speedup through Numba (Issue #364) [2019-Q3 → 2019-Q4]
   -  |uncheck| Dask backend (Issue #259) [2019-Q3 → 2020-Q1]
   -  |uncheck| Attempt to use pyhf as fitting tool for full Analysis Systems
      pipeline test in early 2020 [2019-Q4 → 2020-Q1]
   -  |uncheck| pyhf should satisfy `IRIS-HEP Analysis Systems
      Milestone <https://docs.google.com/spreadsheets/d/1VKpHlQWXu_p8AUv5E5H_BzqF_i7hh2Z-Id0XPwNHu8o/edit#gid=1864915304>`__
      "`GPU/accelerator-based implementation of statistical and other
      appropriate
      components <https://github.com/iris-hep/project-milestones/issues/15>`__"
      [2020-Q1 → 2020-Q2] and contributes to "`Benchmarking and
      assessment of prototype analysis system
      components <https://github.com/iris-hep/project-milestones/issues/17>`__"
      [2020-Q3 → 2020-Q4]

Roadmap as Gantt Chart
~~~~~~~~~~~~~~~~~~~~~~

.. figure:: https://user-images.githubusercontent.com/5142394/64583069-53049180-d355-11e9-8b39-8b2a4599e21e.png
   :alt: pyhf_AS_gantt


Presentations During Roadmap Timeline
-------------------------------------

-  |check| `Talk at IRIS-HEP Institute
   Retreat <https://indico.cern.ch/event/840472/contributions/3564386/>`__
   (September 12-13th, 2019)
-  |check| Talk at `PyHEP 2019 <https://indico.cern.ch/event/833895/>`__
   (October 16-18th, 2019)
-  |check| `Talk at CHEP
   2019 <https://indico.cern.ch/event/773049/contributions/3476143/>`__
   (November 4-8th, 2019)
-  |check| `Poster at CHEP
   2019 <https://indico.cern.ch/event/773049/contributions/3476180/>`__
   (November 4-8th, 2019)

.. |Seibert Hierarchy of Needs SciPy 2019| image:: https://pbs.twimg.com/media/EBYojw8XUAERJhZ?format=png

.. |check| raw:: html

    <input checked=""  type="checkbox" disabled="true">

.. |uncheck| raw:: html

    <input type="checkbox" disabled="true">
