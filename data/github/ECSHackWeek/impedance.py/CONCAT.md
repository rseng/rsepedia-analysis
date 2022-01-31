[![DOI](https://zenodo.org/badge/136110609.svg)](https://zenodo.org/badge/latestdoi/136110609)  ![GitHub release](https://img.shields.io/github/release/ECSHackWeek/impedance.py)

![PyPI - Downloads](https://img.shields.io/pypi/dm/impedance?style=flat-square)  [![All Contributors](https://img.shields.io/badge/all_contributors-11-orange.svg?style=flat-square)](#contributors)

[![Build Status](https://travis-ci.org/ECSHackWeek/impedance.py.svg?branch=master&kill_cache=1)](https://travis-ci.org/ECSHackWeek/impedance.py)  [![Documentation Status](https://readthedocs.org/projects/impedancepy/badge/?version=latest&kill_cache=1)](https://impedancepy.readthedocs.io/en/latest/?badge=latest) [![Coverage Status](https://coveralls.io/repos/github/ECSHackWeek/impedance.py/badge.svg?branch=master&kill_cache=1)](https://coveralls.io/github/ECSHackWeek/impedance.py?branch=master)

impedance.py
------------

`impedance.py` is a Python package for making electrochemical impedance spectroscopy (EIS) analysis reproducible and easy-to-use.

Aiming to create a consistent, [scikit-learn-like API](https://arxiv.org/abs/1309.0238) for impedance analysis, impedance.py contains modules for data preprocessing, validation, model fitting, and visualization.

For a little more in-depth discussion of the package background and capabilities, check out our [Journal of Open Source Software paper](https://joss.theoj.org/papers/10.21105/joss.02349).

If you have a feature request or find a bug, please [file an issue](https://github.com/ECSHackWeek/impedance.py/issues) or, better yet, make the code improvements and [submit a pull request](https://help.github.com/articles/creating-a-pull-request-from-a-fork/)! The goal is to build an open-source tool that the entire impedance community can improve and use!

### Installation

The easiest way to install impedance.py is from [PyPI](https://pypi.org/project/impedance/) using pip.

```bash
pip install impedance
```

See [Getting started with impedance.py](https://impedancepy.readthedocs.io/en/latest/getting-started.html) for instructions on getting started from scratch.

#### Dependencies

impedance.py requires:

-   Python (>=3.6)
-   SciPy (>=1.0)
-   NumPy (>=1.14)
-   Matplotlib (>=3.0)
-   Altair (>=3.0)

Several example notebooks are provided in the `docs/source/examples/` directory. Opening these will require Jupyter notebook or Jupyter lab.

#### Examples and Documentation

Several examples can be found in the `docs/source/examples/` directory (the [Fitting impedance spectra notebook](https://impedancepy.readthedocs.io/en/latest/examples/fitting_example.html) is a great place to start) and the documentation can be found at [impedancepy.readthedocs.io](https://impedancepy.readthedocs.io/en/latest/).

## Citing impedance.py

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02349/status.svg)](https://doi.org/10.21105/joss.02349)

If you use impedance.py in published work, please consider citing https://joss.theoj.org/papers/10.21105/joss.02349 as

```bash
@article{Murbach2020,
  doi = {10.21105/joss.02349},
  url = {https://doi.org/10.21105/joss.02349},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {52},
  pages = {2349},
  author = {Matthew D. Murbach and Brian Gerwe and Neal Dawson-Elli and Lok-kun Tsui},
  title = {impedance.py: A Python package for electrochemical impedance analysis},
  journal = {Journal of Open Source Software}
}
```

## Contributors ✨

This project started at the [2018 Electrochemical Society (ECS) Hack Week in Seattle](https://www.electrochem.org/233/hack-week) and has benefited from a community of users and contributors since. Thanks goes to these wonderful people ([emoji key](https://allcontributors.org/docs/en/emoji-key)):

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/lktsui"><img src="https://avatars0.githubusercontent.com/u/22246069?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Lok-kun Tsui</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=lktsui" title="Code">💻</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=lktsui" title="Tests">⚠️</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=lktsui" title="Documentation">📖</a></td>
    <td align="center"><a href="https://github.com/BGerwe"><img src="https://avatars3.githubusercontent.com/u/38819321?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Brian Gerwe</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=BGerwe" title="Code">💻</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=BGerwe" title="Tests">⚠️</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=BGerwe" title="Documentation">📖</a> <a href="https://github.com/ECSHackWeek/impedance.py/pulls?q=is%3Apr+reviewed-by%3ABGerwe" title="Reviewed Pull Requests">👀</a></td>
    <td align="center"><a href="https://github.com/nealde"><img src="https://avatars2.githubusercontent.com/u/25877868?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Neal</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/pulls?q=is%3Apr+reviewed-by%3Anealde" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=nealde" title="Code">💻</a></td>
    <td align="center"><a href="http://mattmurbach.com"><img src="https://avatars3.githubusercontent.com/u/9369020?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matt Murbach</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=mdmurbach" title="Documentation">📖</a> <a href="https://github.com/ECSHackWeek/impedance.py/pulls?q=is%3Apr+reviewed-by%3Amdmurbach" title="Reviewed Pull Requests">👀</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=mdmurbach" title="Tests">⚠️</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=mdmurbach" title="Code">💻</a></td>
    <td align="center"><a href="https://kennyvh.com"><img src="https://avatars2.githubusercontent.com/u/29909203?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Kenny Huynh</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3Ahkennyv" title="Bug reports">🐛</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=hkennyv" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/lawrencerenna"><img src="https://avatars0.githubusercontent.com/u/49174337?v=4?s=100" width="100px;" alt=""/><br /><sub><b>lawrencerenna</b></sub></a><br /><a href="#ideas-lawrencerenna" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Rowin"><img src="https://avatars3.githubusercontent.com/u/1727478?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Rowin</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3ARowin" title="Bug reports">🐛</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=Rowin" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/michaelplews"><img src="https://avatars2.githubusercontent.com/u/14098929?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michael Plews</b></sub></a><br /><a href="#ideas-michaelplews" title="Ideas, Planning, & Feedback">🤔</a></td>
    <td align="center"><a href="https://github.com/Chebuskin"><img src="https://avatars0.githubusercontent.com/u/33787723?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Chebuskin</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3AChebuskin" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/environmat"><img src="https://avatars0.githubusercontent.com/u/9309353?v=4?s=100" width="100px;" alt=""/><br /><sub><b>environmat</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3Aenvironmat" title="Bug reports">🐛</a></td>
    <td align="center"><a href="http://www.abdullahsumbal.com"><img src="https://avatars2.githubusercontent.com/u/12946947?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Abdullah Sumbal</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3Aabdullahsumbal" title="Bug reports">🐛</a></td>
    <td align="center"><a href="https://github.com/nobkat"><img src="https://avatars3.githubusercontent.com/u/29077445?v=4?s=100" width="100px;" alt=""/><br /><sub><b>nobkat</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=nobkat" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/nickbrady"><img src="https://avatars1.githubusercontent.com/u/7471367?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Nick</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/issues?q=author%3Anickbrady" title="Bug reports">🐛</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=nickbrady" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/aokomorowski"><img src="https://avatars.githubusercontent.com/u/43665474?v=4?s=100" width="100px;" alt=""/><br /><sub><b>aokomorowski</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=aokomorowski" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://petermattia.com"><img src="https://avatars.githubusercontent.com/u/29551858?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Peter Attia</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=petermattia" title="Code">💻</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=petermattia" title="Tests">⚠️</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=petermattia" title="Documentation">📖</a></td>
    <td align="center"><a href="http://sdkang.org"><img src="https://avatars.githubusercontent.com/u/55116501?v=4?s=100" width="100px;" alt=""/><br /><sub><b>sdkang</b></sub></a><br /><a href="https://github.com/ECSHackWeek/impedance.py/commits?author=stephendkang" title="Tests">⚠️</a> <a href="https://github.com/ECSHackWeek/impedance.py/commits?author=stephendkang" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

This project follows the [all-contributors](https://github.com/all-contributors/all-contributors) specification. Contributions of any kind welcome!
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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at mmurbach AT uw DOT edu. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
:tada: Welcome to the open-source impedance analysis community! :tada:

Everyone is welcome to contribute.
We value all forms of contributions including code reviews, bug fixes, new features, examples, documentation, community participation, etc.

This document outlines the guidelines for contributing to the various aspects of the project.

# Contributing

If you find a bug in the code or a mistake in the [documentation](https://impedancepy.readthedocs.io/en/latest/?badge=latest) or want a new feature, you can help us by creating [an issue in our repository](https://github.com/ECSHackWeek/impedance.py/issues), or even submit a pull request.

# Development Guide

## Repository Setup

1.  To work on the impedance.py package, you should first fork the repository on GitHub using the button on the top right of the ECSHackWeek/impedance.py repository.

2.  You can then clone the fork to your computer

```bash
git clone https://github.com/<GitHubUsername>/impedance.py
```

3.  Make your changes and commit them to your fork (for an introduction to git, checkout the [tutorial from the ECS Hack Week](https://github.com/ECSHackWeek/ECSHackWeek_Dallas/blob/master/Version_Control.pptx))

For example,
```bash
git add changedfiles
git commit
git push
```

4.  [Submit a Pull Request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/about-pull-requests) (make sure to write a good message so the reviewer can understand what you're adding!) via GitHub.

5.  Add yourself to the list of collaborators (you can use the [all-contributors bot](https://allcontributors.org/docs/en/bot/usage))! You rock!

## Continuous Integration

`impedance.py` uses [Travis CI](https://travis-ci.org/ECSHackWeek/impedance.py) for Continuous Integration testing. This means that every time you submit a pull request, a series of tests will be run to make sure the changes don’t accidentally introduce any bugs :bug:. *Your PR will not be accepted until it passes all of these tests.* While you can certainly wait for the results of these tests after submitting a PR, you can also run them locally to speed up the review process.

We use flake8 to test for PEP 8 conformance. [PEP 8](https://www.python.org/dev/peps/pep-0008/) is a series of style guides for Python that provide suggestions for everything from variable naming to indentation. 
To run the flake8 (PEP8) code style checker:

```
conda install flake8
cd impedance.py/
flake8
```
:warning: if there is any output here, fix the errors and try running flake8 again.

To run the tests using py.test:

```
conda install pytest
cd impedance.py/
pytest
```
:warning: you should see all tests pass, if not try fixing the error or file an issue.

### Unit Tests

`impedance.py` aims to have complete test coverage of our package code. If you're adding a new feature, consider writing the test first and then the code to ensure it passes. PRs which decrease code coverage will need to add tests before they can be merged.
---
title: 'impedance.py: A Python package for electrochemical impedance analysis'
tags:
  - Python
  - electrochemistry
  - impedance
  - lithium-ion batteries
  - fuel cells
  - corrosion
authors:
  - name: Matthew D. Murbach
    email: matt@hivebattery.com
    orcid: 0000-0002-6583-5995
    affiliation: 1
  - name: Brian Gerwe
    email: brian.s.gerwe@gmail.com
    orcid: 0000-0002-1184-8483
    affiliation: 2
  - name: Neal Dawson-Elli
    affiliation: 3
    email: ndawsonelli@gmail.com
    orcid: 0000-0003-4559-7309
  - name: Lok-kun Tsui
    email: lktsui@unm.edu
    orcid: 0000-0001-7381-0686
    affiliation: 4
affiliations:
 - name: Hive Battery Labs, Seattle, WA, USA
   index: 1
 - name: University of Washington, Department of Chemical Engineering, Seattle, WA, USA
   index: 2
 - name: PayScale, Seattle, WA, USA
   index: 3
 - name: University of New Mexico, Center for MicroEngineered Materials, Albuquerque, NM, USA
   index: 4
date: 14 March 2020
bibliography: paper.bib
---

`impedance.py` is a community-driven Python package for making the analysis of electrochemical impedance spectroscopy (EIS)  data easier and more reproducible. `impedance.py` currently provides several useful features commonly used in the typical  impedance analysis workflow:

- _preprocessing_: functions for importing data from a variety of instruments and file types
- _data validation_: easy-to-use methods for checking measurement validity
- _model fitting_: a simple and powerful interface for quickly fitting models to analyze data
- _model selection_: parameter error estimates and model confidence bounds
- _visualization_: interactive and publication-ready Nyquist and Bode plots

# Background

Electrochemical impedance spectroscopy (EIS) is a powerful technique for noninvasively probing the physicochemical processes governing complex electrochemical systems by measuring the relationship between voltage and current (i.e. the impedance) as a function of frequency [@orazem_electrochemical_2008]. Analysis of EIS spectra often involves the construction and fitting of models to represent the system response. Equivalent circuit models, a popular method for interpreting impedance spectra, represent processes in an electrochemical system with elements in a circuit allowing quantitative characterization of the data. EIS inherently requires the system response is linear, causal and stable over the measurement time. As such, researchers must also ensure their data complies with these assumptions – often through some implementation of Kramers-Kronig analysis. Following the approach outlined above, workers have applied EIS to problems ranging from  measuring the dielectric properties of biological systems [@semenov_dielectric_2001] to corrosion rates of coated metals [@bonora_corrosion_1996]. More recently, EIS has been widely adopted for studying processes governing batteries [@rodrigues_batteries_2000; @troltzsch_lib_2006]  and fuel cells [@springer_pemfc_1996; @wagner_sofc_1998].

# Statement of Need
To date, typical impedance analysis solutions have relied on either instrument-specific, proprietary software or ad hoc, lab-specific code written for internal use. These proprietary or home-grown solutions lack cross-platform support, and users of Linux or MacOS are unable to perform analysis on their computers. In addition to making reproducible analysis difficult, these solutions can be restrictive to defining custom circuit elements and processing large datasets. By providing an open-source, cross-platform (Windows, MacOS, and Linux), community-driven package for the full impedance analysis pipeline from data management to parameter extraction and publication-ready figures, `impedance.py` seeks to encourage reproducible, easy-to-use, and transparent analysis. In addition to decades of electrochemical research, many new methods for validating [@schonleber_method_2014] and analyzing [@murbach_analysis_2018; @buteau_analysis_2019] impedance spectra have been developed by researchers. By lowering the barrier to use tried-and-true methods along side cutting-edge analytical techniques via a consistent interface, `impedance.py` also serves to grow as a community repository of best-practices while facilitating the adoption of new techniques.

# Current `impedance.py` functionality


## Preprocessing
Once you've collected your data, the last thing you want to do is worry about how to turn the files generated by the instrument into results you can learn from. `impedance.py` offers easy-to-use functions for importing data from BioLogic, CH Instruments, Gamry, PARSTAT, VersaStudio, and ZView files. Don't see your instrument yet? [Create an issue with a sample file](https://github.com/ECSHackWeek/impedance.py/issues/new?assignees=&labels=&template=data-file-support-request.md&title=%5BDATA%5D) and help contribute to the project!

## Data validation
Ensuring that the system under measurement is linear, stable, and causal is an important, but often overlooked part of impedance analysis. `impedance.py` provides several methods  (measurement models, Lin-KK algorithm [@schonleber_method_2014]) for data validation as a part of the same easy-to-use package.

## Equivalent circuit fitting
`impedance.py` equivalent circuit fitting combines an extremely flexible circuit definition scheme with a range of available circuit elements -- from simple L, R, C, and constant-phase elements to electrochemistry specific elements such as Warburg, Gerischer [@boukamp_gerischer_element_2003], and transmission-line models [@paasch_theory_1993]. Custom elements can be easily added to a single users installation or to the overall project by [creating an issue](https://github.com/ECSHackWeek/impedance.py/issues/new?assignees=&labels=&template=equivalent-circuit-element-request.md&title=%5BElement%5D). Additional model fitting features include parameter bounds, holding circuit elements constant, weighted fitting, and even saving/loading models to a human readable .json file. Fitting is performed by non-linear least squares regression of the circuit model to impedance data via `curve_fit` from the `scipy.optimize` package [@2020SciPy-NMeth].

## Model selection
One of the biggest challenges in extracting parameters from EIS spectra is quantitatively evaluating the fit of different models. `impedance.py` returns estimated parameter error bars with every fit and comes with a built in method for calculating the confidence bounds of a fit model.

## Data visualization
EIS is often qualitatively interpreted by examing the shape of a spectrum which relies on properly defined axes. `impedance.py` eliminates the tedium of manually formatting plots by providing built-in methods to generate clean, publication-ready figures including interactive plots via Altair [@vanderplas_altair_2018] and customizable Nyquist or Bode plots via matplotlib [@hunter_matplotlib_2007].


# An example of the simple model API
The documentation for `impedance.py` contains [a guide on getting started](https://impedancepy.readthedocs.io/en/latest/getting-started.html) and several examples of what a typical analysis workflow might look like using the package. Here we show how importing data, defining and fitting an equivalent circuit  model, and visualizing the results can be done with just a handful of lines in `impedance.py`:

```python
# 1. loading in data:
from impedance.preprocessing import readFile
f, Z = readFile('exampleData.csv')

# 2. remove positive imaginary impedance data:
from impedance.preprocessing import ignoreBelowX
f, Z = ignoreBelowX(f, Z)

# 3. importing and initializing a circuit:
from impedance.models.circuits import CustomCircuit
initial_guess = [1e-8, .01, .005, .1, .9, .005, .1, 200, .1, .9]
circuit = CustomCircuit('L_0-R_0-p(R_1,CPE_1)-p(R_2-Wo_1,CPE_2)',
                        initial_guess=initial_guess)

# 4. fitting the circuit to the data:
circuit.fit(f, Z)
print(circuit)
```

```text
Circuit string: L_0-R_0-p(R_1,CPE_1)-p(R_2-Wo_1,CPE_2)
Fit: True

Initial guesses:
      L_0 = 1.00e-08 [H]
      R_0 = 1.00e-02 [Ohm]
      R_1 = 5.00e-03 [Ohm]
  CPE_1_0 = 1.00e-01 [Ohm^-1 sec^a]
  CPE_1_1 = 9.00e-01 []
      R_2 = 5.00e-03 [Ohm]
   Wo_1_0 = 1.00e-01 [Ohm]
   Wo_1_1 = 2.00e+02 [sec]
  CPE_2_0 = 1.00e-01 [Ohm^-1 sec^a]
  CPE_2_1 = 9.00e-01 []

Fit parameters:
      L_0 = 1.28e-07  (+/- 4.92e-08) [H]
      R_0 = 1.49e-02  (+/- 8.31e-04) [Ohm]
      R_1 = 8.31e-03  (+/- 2.24e-03) [Ohm]
  CPE_1_0 = 1.13e+00  (+/- 8.07e-01) [Ohm^-1 sec^a]
  CPE_1_1 = 6.85e-01  (+/- 1.27e-01) []
      R_2 = 7.92e-03  (+/- 1.51e-03) [Ohm]
   Wo_1_0 = 1.37e-01  (+/- 9.03e-02) [Ohm]
   Wo_1_1 = 1.23e+03  (+/- 1.61e+03) [sec]
  CPE_2_0 = 3.69e+00  (+/- 3.09e-01) [Ohm^-1 sec^a]
  CPE_2_1 = 9.79e-01  (+/- 6.02e-02) []
```

```python
# 5. visualize the results:
circuit.plot(f_data=f, Z_data=Z)
```

![Interactive impedance plots are as easy as `.plot()`!](./example.png)

# Acknowledgements

We thank participants on the 2018 Electrochemical Society (ECS) Hack Week team in Seattle, WA as well as Dan Schwartz and David Beck for their guidance. An up-to-date [list of contributors can be found on GitHub](https://github.com/ECSHackWeek/impedance.py#contributors-). Example data is from [@murbach_nonlinear_2018].

# References
---
name: Data file support request
about: Request a new file type for preprocessing.readFile()
title: "[DATA]"
labels: ''
assignees: ''

---

* What instrument/software/file type are you requesting support for? *

[Please attach a sample file by dragging and dropping an example file (this file will be public)]
---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Code which generates the error '...'
2. Full error info

**Expected behavior**
A clear and concise description of what you expected to happen.

**Additional context**
Add any other context about the problem here.
---
name: Equivalent circuit element request
about: Request a new equivalent circuit element to be added
title: "[Element]"
labels: ''
assignees: ''

---

Please submit a reference and the equation for a proposed element! Attaching an example data set is also usually helpful :)
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
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
Circuit Elements
================

.. automodule:: impedance.models.circuits.elements
   :members:
Circuits
========

.. automodule:: impedance.models.circuits.circuits
   :members:
=========================================
Getting started with :code:`impedance.py`
=========================================

:code:`impedance.py` is a Python package for analyzing electrochemical impedance
spectroscopy (EIS) data. The following steps will show you how to get started
analyzing your own data using :code:`impedance.py` in a Jupyter notebook.

.. hint::
  If you get stuck or believe you have found a bug, please feel free to open an
  `issue on GitHub <https://github.com/ECSHackWeek/impedance.py/issues>`_.

Step 1: Installation
====================

If you already are familiar with managing Python packages, feel free to skip
straight to `Installing packages`_.
Otherwise, what follows is a quick introduction to the Python package ecosystem:

Installing Miniconda
--------------------
One of the easiest ways to get started with Python is using Miniconda.
Installation instructions for your OS can be found at
https://conda.io/miniconda.html.

After you have installed conda, you can run the following commands in your
Terminal/command prompt/Git BASH to update and test your installation:

1. Update conda’s listing of packages for your system: :code:`conda update conda`

2. Test your installation: :code:`conda list`

  For a successful installation, a list of installed packages appears.

3. Test that Python 3 is your default Python: :code:`python -V`

  You should see something like Python 3.x.x :: Anaconda, Inc.

You can interact with Python at this point, by simply typing :code:`python`.

Setting up a conda environment
------------------------------
*(Optional)* It is recommended that you use virtual environments to keep track
of the packages you've installed for a particular project. Much more info on
how conda makes this straightforward is given `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#viewing-a-list-of-your-environments>`_.

We will start by creating an environment called :code:`impedance-analysis`
which contains all the Python base distribution:

.. code-block:: bash

   conda create -n impedance-analysis python=3

After conda creates this environment, we need to activate it before we can
install anything into it by using:

.. code-block:: bash

   conda activate impedance-analysis

We've now activated our conda environment and are ready to install
:code:`impedance.py`!

Installing packages
-------------------

The easiest way to install :code:`impedance.py` and it's dependencies
(:code:`scipy`, :code:`numpy`, and :code:`matplotlib`) is from
`PyPI <https://pypi.org/project/impedance/>`_ using pip:

.. code-block:: bash

   pip install impedance

For this example we will also need Jupyter Lab which we can install with:

.. code-block:: bash

   conda install jupyter jupyterlab

We've now got everything in place to start analyzing our EIS data!

.. note::
  The next time you want to use this same environment, all you have to do is
  open your terminal and type :code:`conda activate impedance-analysis`.

Open Jupyter Lab
----------------

*(Optional)* Create a directory in your documents folder for this example:

.. code-block:: bash

   mkdir ~/Documents/impedance-example

   cd ~/Documents/impedance-example

Next, we will launch an instance of Jupyter Lab:

.. code-block:: bash

  jupyter lab

which should open a new tab in your browser. A tutorial on Jupyter Lab from the
Electrochemical Society HackWeek can be found `here <https://ecshackweek.github.io/tutorial/getting-started-with-jupyter/>`_.

.. tip::
  The code below can be found in the `getting-started.ipynb <_static/getting-started.ipynb>`_ notebook

Step 2: Import your data
========================

This example will assume the following dataset is located in your current working directory (feel free to replace it with your data): :download:`exampleData.csv <_static/exampleData.csv>`

For this dataset, importing the data looks something like:

.. code-block:: python

  from impedance import preprocessing

  # Load data from the example EIS result
  frequencies, Z = preprocessing.readCSV('./exampleData.csv')

  # keep only the impedance data in the first quandrant
  frequencies, Z = preprocessing.ignoreBelowX(frequencies, Z)

.. tip::
  Functions for reading in files from a variety of vendors (ZPlot, Gamry, Parstat, Autolab, ...) can be found in the `preprocessing module <preprocessing.html>`_!

Step 3: Define your impedance model
===================================

Next we want to define our impedance model. In order to enable a wide variety
of researchers to use the tool, :code:`impedance.py` allows you to define a
custom circuit with any combination of `circuit elements <circuit-elements.html>`_.

The circuit is defined as a string (i.e. using :code:`''` in Python), where elements in
series are separated by a dash (:code:`-`), and elements in parallel are wrapped in
a :code:`p( , )`. Each element is defined by the function (in `circuit-elements.py <circuit-elements.html>`_) followed by a single digit identifier.

For example, the circuit below:

.. image:: _static/two_time_constants.png

would be defined as :code:`R0-p(R1,C1)-p(R2-Wo1,C2)`.

Each circuit, we want to fit also needs to have an initial guess for each
of the parameters. These inital guesses are passed in as a list in order the
parameters are defined in the circuit string. For example, a good guess for this
battery data is :code:`initial_guess = [.01, .01, 100, .01, .05, 100, 1]`.

We create the circuit by importing the CustomCircuit object and passing it our
circuit string and initial guesses.

.. code-block:: python

  from impedance.models.circuits import CustomCircuit

  circuit = 'R0-p(R1,C1)-p(R2-Wo1,C2)'
  initial_guess = [.01, .01, 100, .01, .05, 100, 1]

  circuit = CustomCircuit(circuit, initial_guess=initial_guess)

Step 4: Fit the impedance model to data
=======================================

Once we've defined our circuit, fitting it to impedance data is as easy as
calling the `.fit()` method and passing it our experimental data!

.. code-block:: python

  circuit.fit(frequencies, Z)

We can access the fit parameters with :code:`circuit.parameters_` or by
printing the circuit object itself, :code:`print(circuit)`.

Step 5: Analyze/Visualize the results
=====================================

For this dataset, the resulting fit parameters are

================ ========
   Parameter      Value
---------------- --------
:math:`R_0`      1.65e-02
:math:`R_1`      8.68e-03
:math:`C_1`      3.32e+00
:math:`R_2`      5.39e-03
:math:`Wo_{1,0}` 6.31e-02
:math:`Wo_{1,1}` 2.33e+02
:math:`C_2`      2.20e-01
================ ========

We can get the resulting fit impedance by passing a list of frequencies to the :code:`.predict()` method.

.. code-block:: python

  Z_fit = circuit.predict(frequencies)

To easily visualize the fit, the :code:`plot_nyquist()` function can be handy.

.. code-block:: python

  import matplotlib.pyplot as plt
  from impedance.visualization import plot_nyquist

  fig, ax = plt.subplots()
  plot_nyquist(ax, Z, fmt='o')
  plot_nyquist(ax, Z_fit, fmt='-')

  plt.legend(['Data', 'Fit'])
  plt.show()

.. image:: _static/example_fit_fig.png

.. important::
  🎉 Congratulations! You're now up and running with impedance.py 🎉
Examples
========


.. toctree::
    :maxdepth: 1
    :glob:

    examples/fitting_example.ipynb
    examples/plotting_example.ipynb
    examples/model_io_example.ipynb
    examples/validation_example.ipynb
    examples/looping_files_example.ipynb
Validation
==========

Interpreting EIS data fundamentally relies on the the system conforming to
conditions of causality, linearity, and stability. For an example of how
the adherence to the Kramers-Kronig relations, see the `Validation Example Jupyter Notebook
<examples/validation_example.ipynb>`_

Lin-KK method
-------------

Validating your data with the lin-KK model requires fitting an optimal number
of RC-elements and analysis of the residual errors.

.. automodule:: impedance.validation
   :members:
Preprocessing
=============

.. automodule:: impedance.preprocessing
   :members:
Fitting
=======

.. automodule:: impedance.models.circuits.fitting
   :members:
Frequently Asked Questions
==========================

What method does impedance.py use for fitting equivalent circuit models?
------------------------------------------------------------------------
By default, fitting is performed by non-linear least squares regression of
the circuit model to impedance data via
`curve_fit <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_
from the `scipy.optimize` package.[1]
Real and imaginary components are fit simultaneously with uniform
weighting, i.e. the objective function to minimize is,

.. math::
    \chi^2 = \sum_{n=0}^{N} [Z^\prime_{data}(\omega_n) - Z^\prime_{model}(\omega_n)]^2 +
                   [Z^{\prime\prime}_{data}(\omega_n) - Z^{\prime\prime}_{model}(\omega_n)]^2

where N is the number of frequencies and :math:`Z^\prime` and
:math:`Z^{\prime\prime}` are the real and imaginary components of
the impedance, respectively.
The default optimization method is the
Levenberg-Marquardt algorithm (:code:`method='lm'`) for unconstrained
problems and the Trust Region Reflective algorithm
(:code:`method='trf'`) if bounds are provided. See `the SciPy documentation
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_
for more details and options.

While the default method converges quickly and often yields acceptable fits,
the results may be sensitive to the initial conditions.
EIS fitting can be prone to this issue given the high dimensionality
of typical equivalent circuit models.
`Global optimization algorithms <https://en.wikipedia.org/wiki/Global_optimization>`_
attempt to search the entire parameter landscape to minimize the error.
By setting :code:`global_opt=True` in :code:`circuit_fit`, :code:`impedance.py` will use the
`basinhopping <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html>`_
global optimization algorithm (also from the `scipy.optimize` package[1]) instead of :code:`curve_fit`.
Note that the computational time may increase.

[1] Virtanen, P., Gommers, R., Oliphant, T.E. et al.
SciPy 1.0: fundamental algorithms for scientific computing in Python.
Nat Methods 17, 261–272 (2020). `doi: 10.1038/s41592-019-0686-2 <https://doi.org/10.1038/s41592-019-0686-2>`_

How do I cite impedance.py?
---------------------------

.. image:: https://joss.theoj.org/papers/10.21105/joss.02349/status.svg
    :target: https://doi.org/10.21105/joss.02349

If you use impedance.py in published work, please consider citing https://joss.theoj.org/papers/10.21105/joss.02349 as

.. code:: text

    @article{Murbach2020,
        doi = {10.21105/joss.02349},
        url = {https://doi.org/10.21105/joss.02349},
        year = {2020},
        publisher = {The Open Journal},
        volume = {5},
        number = {52},
        pages = {2349},
        author = {Matthew D. Murbach and Brian Gerwe and Neal Dawson-Elli and Lok-kun Tsui},
        title = {impedance.py: A Python package for electrochemical impedance analysis},
        journal = {Journal of Open Source Software}
    }

How can I contribute to impedance.py?
-------------------------------------

First off, thank you for your interest in contributing to the
open-source electrochemical community! We're excited to welcome all
contributions including suggestions for new features, bug reports/fixes,
examples/documentation, and additional impedance analysis functionality.

Feature requests and bug reports
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to make a suggestion for a new feature, please `make an
issue <https://github.com/ECSHackWeek/impedance.py/issues/new/choose>`_
including as much detail as possible. If you're requesting a
new circuit element or data file type, there are special issue templates
that you can select and use.

Contributing code
~~~~~~~~~~~~~~~~~

The prefered method for contributing code to impedance.py is to fork
the repository on GitHub and submit a "pull request" (PR).
More detailed information on how to get started developing impedance.py
can be found in
`CONTRBUTING.md <https://github.com/ECSHackWeek/impedance.py/blob/master/CONTRIBUTING.md>`_.

Feel free to reach out via GitHub issues with any questions!

.. eisfit documentation master file, created by
   sphinx-quickstart on Wed May 16 16:54:07 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

impedance.py
=============

:code:`impedance.py` is a Python package for making 
electrochemical impedance spectroscopy (EIS) analysis
reproducible and easy-to-use.

Aiming to create a consistent,
`scikit-learn-like API <https://arxiv.org/abs/1309.0238>`_
for impedance analysis, :code:`impedance.py` contains
modules for data preprocessing, validation, model fitting,
and visualization.


If you have a feature request or find a bug, please 
`file an issue <https://github.com/ECSHackWeek/impedance.py/issues>`_
or, better yet, make the code improvements and 
`submit a pull request <https://help.github.com/articles/creating-a-pull-request-from-a-fork/>`_!
The goal is to build an open-source tool that the
entire impedance community can improve and use!

Installation
------------

The easiest way to install :code:`impedance.py` is
from `PyPI <https://pypi.org/project/impedance/>`_ 
using pip:

.. code-block:: bash

   pip install impedance

See :doc:`./getting-started` for instructions
on getting started from scratch.

Dependencies
~~~~~~~~~~~~

impedance.py requires:

-   Python (>=3.7)
-   SciPy (>=1.0)
-   NumPy (>=1.14)
-   Matplotlib (>=3.0)
-   Altair (>=3.0)

Several example notebooks are provided in the examples/ directory.
Opening these will require Jupyter notebook or Jupyter lab.

Examples and Documentation
---------------------------
:doc:`./getting-started` contains a detailed walk
through of how to get started from scratch. If you're already familiar with
Jupyter/Python, several examples can be found in the :code:`examples/` directory
(:doc:`./examples/fitting_example` is a great place to start). 
The documentation can be found at 
`impedancepy.readthedocs.io <https://impedancepy.readthedocs.io/en/latest/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting-started
   examples
   preprocessing
   validation
   circuits
   circuit-elements
   fitting
   faq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
