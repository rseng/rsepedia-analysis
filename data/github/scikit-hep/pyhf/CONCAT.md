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
