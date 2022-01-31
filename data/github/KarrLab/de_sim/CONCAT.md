[![PyPI package](https://img.shields.io/pypi/v/de_sim.svg)](https://pypi.python.org/pypi/de_sim)
[![Documentation](https://readthedocs.org/projects/de-sim/badge/?version=latest)](https://docs.karrlab.org/de_sim)
[![Test results](https://circleci.com/gh/KarrLab/de_sim.svg?style=shield)](https://circleci.com/gh/KarrLab/de_sim)
[![Test coverage](https://coveralls.io/repos/github/KarrLab/de_sim/badge.svg)](https://coveralls.io/github/KarrLab/de_sim)
[![Code analysis](https://api.codeclimate.com/v1/badges/2fa3ece22f571fd36b12/maintainability)](https://codeclimate.com/github/KarrLab/de_sim)
[![License](https://img.shields.io/github/license/KarrLab/de_sim.svg)](LICENSE)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02685/status.svg)](https://doi.org/10.21105/joss.02685)
![Analytics](https://ga-beacon.appspot.com/UA-86759801-1/de_sim/README.md?pixel)

# *DE-Sim*: a Python-based object-oriented discrete-event simulator for modeling complex systems

*DE-Sim* is an open-source, Python-based object-oriented discrete-event simulation (DES) tool that makes it easy to use large, heterogeneous datasets and high-level data science tools such as [NumPy](https://numpy.org/), [Scipy](https://scipy.org/scipylib/index.html), [pandas](https://pandas.pydata.org/), and [SQLAlchemy](https://www.sqlalchemy.org/) to build and simulate complex computational models. Similar to [Simula](http://www.simula67.info/), *DE-Sim* models are implemented by defining logical process objects which read the values of a set of variables and schedule events to modify their values at discrete instants in time.

To help users build and simulate complex, data-driven models, *DE-Sim* provides the following features:

* **High-level, object-oriented modeling:** *DE-Sim* makes it easy for users to use object-oriented Python programming to build models. This makes it easy to use large, heterogeneous datasets and high-level data science packages such as NumPy, pandas, SciPy, and SQLAlchemy to build complex models.
* **Stop conditions:** DE-Sim makes it easy to terminate simulations when specific criteria are reached. Researchers can specify stop conditions as functions that return true when a simulation should conclude.
* **Results checkpointing:** DE-Sim makes it easy to record the results of simulations by using a configurable checkpointing module.
* **Reproducible simulations:** To help researchers debug simulations, repeated executions of the same simulation with the same configuration and same random number generator seed produce the same results.
* **Space-time visualizations:** DE-Sim generates space-time visualizations of simulation trajectories. These diagrams can help researchers understand and debug simulations.

## Projects that use *DE-Sim*
*DE-Sim* has been used to develop [WC-Sim](https://github.com/KarrLab/wc_sim), a multi-algorithmic simulator for [whole-cell models](https://www.wholecell.org).

## Examples
* [Minimal simulation](de_sim/examples/minimal_simulation.py): a minimal example of a simulation
* [Random walk](de_sim/examples/random_walk.py): a random one-dimensional walk which increments or decrements a variable with equal probability at each event
* [Parallel hold (PHOLD)](de_sim/examples/phold.py): model developed by Richard Fujimoto for benchmarking parallel DES simulators
* [Epidemic](https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/sirs.py): an SIR model of an epidemic of an infectious disease

## Tutorial
Please see [sandbox.karrlab.org](https://sandbox.karrlab.org/tree/de_sim) for interactive tutorials on creating and executing models with *DE-Sim*.

## Template for models and simulations
[`de_sim/examples/minimal_simulation.py`](de_sim/examples/minimal_simulation.py) contains a template for implementing and simulating a model with *DE-Sim*.

## Installation
1. Install dependencies
    
    * Python >= 3.7
    * pip >= 19

2. Install this package using one of these methods

    * Install the latest release from PyPI
      ```
      pip install de_sim
      ```

    * Install a Docker image with the latest release from DockerHub
      ```
      docker pull karrlab/de_sim
      ```

    * Install the latest version from GitHub
      ```
      pip install git+https://github.com/KarrLab/de_sim.git#egg=de_sim
      ```

## API documentation
Please see the [API documentation](https://docs.karrlab.org/de_sim/source/de_sim.html).

## Performance
Please see the [*DE-Sim* article](joss_paper/paper.md) for information about the performance of *DE-Sim*.

## Strengths and weaknesses compared to other DES tools
Please see the [*DE-Sim* article](joss_paper/paper.md) for a comparison of *DE-Sim* with other DES tools.

## License
The package is released under the [MIT license](LICENSE).

## Citing *DE-Sim*
Please use the following reference to cite *DE-Sim*:

Arthur P. Goldberg & Jonathan Karr. (2020). [DE-Sim: an object-oriented, discrete-event simulation tool for data-intensive modeling of complex systems in Python. Journal of Open Source Software, 5(55), 2685.](https://doi.org/10.21105/joss.02685)

## Contributing to *DE-Sim*
We enthusiastically welcome contributions to *DE-Sim*! Please see the [guide to contributing](CONTRIBUTING.md) and the [developer's code of conduct](CODE_OF_CONDUCT.md).

## Development team
This package was developed by the [Karr Lab](https://www.karrlab.org) at the Icahn School of Medicine at Mount Sinai in New York, USA by the following individuals:

* [Arthur Goldberg](https://www.mountsinai.org/profiles/arthur-p-goldberg)
* [Jonathan Karr](https://www.karrlab.org)

## Acknowledgements
This work was supported by National Science Foundation award 1649014, National Institutes of Health award R35GM119771, and the Icahn Institute for Data Science and Genomic Technology.

## Questions and comments
Please submit questions and issues to [GitHub](https://github.com/KarrLab/de_sim/issues) or contact the [Karr Lab](mailto:info@karrlab.org).
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
[info@karrlab.org](mailto:info@karrlab.org).
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org/version/2/0/code_of_conduct.html).

The Community Impact Guidelines were inspired by [Mozilla's Code of Conduct Enforcement ladder](https://github.com/mozilla/diversity).

Answers to common questions about this code of conduct are available at
https://www.contributor-covenant.org/faq. Translations of this code of conduct are available at
https://www.contributor-covenant.org/translations.
# Contributing to `DE-Sim`

We enthusiastically welcome contributions to `DE-Sim`!

## Coordinating contributions

Before getting started, please contact the lead developers at [info@karrlab.org](mailto:info@karrlab.org) and [Arthur Goldberg](mailto:Arthur.Goldberg@mssm.edu) to coordinate your planned contributions with other ongoing efforts. Please also use GitHub issues to announce your plans to the community so that other developers can provide input into your plans and coordinate their own work. As the development community grows, we will institute additional infrastructure as needed such as a leadership committee and regular online meetings.

## Repository organization

`DE-Sim` follows standard Python conventions:

* `README.md`: Overview of `DE-Sim`
* `de_sim/`: Source code
* `docs/`: Documentation
* `tests/`: Unit tests
* `pytest.ini`: pytest configuration
* `setup.py`: pip installation script
* `setup.cfg`: Configuration for the pip installation
* `requirements.txt`: Dependencies
* `requirements.optional.txt`: Optional dependencies
* `manifest.in`: List of files to include in package
* `codemeta.json`: Package metadata
* `LICENSE`: License
* `CONTRIBUTING.md`: Guide to contributing to `DE-Sim` (this document)
* `CODE_OF_CONDUCT.md`: Code of conduct for developers

## Coding convention

`DE-Sim` follows standard Python style conventions:

* Class names: `UpperCamelCase`
* Function names: `lower_snake_case`
* Variable names: `lower_snake_case`

## Testing and continuous integration

We strive to have complete test coverage of `DE-Sim`. As such, all contributions to `DE-Sim` should be tested. 

The tests depend on additional Python packages. These can be installed by running one of the following commands:
```
pip install /path/to/de_sim[tests]
pip install -r /path/to/de_sim/tests/requirements.txt
```

The tests are located in the `tests`  directory. The tests can be executed by running the following commands.
Note that the current directory must be on the `PYTHONPATH` environment variable.
```
git clone https://github.com/KarrLab/de_sim.git
cd de_sim
pip install pytest
python -m pytest tests
```

The coverage of the tests can be evaluated by running the following commands and then opening `/path/to/de_sim/htmlcov/index.html` with your browser.
```
pip install pytest pytest-cov coverage
python -m pytest tests --cov de_sim
coverage html
```

Upon each push to GitHub, GitHub will trigger CircleCI to execute all of the tests. All test results are available at [CircleCI](https://circleci.com/gh/KarrLab/de_sim). And all coverage reports are available at [Coveralls](https://coveralls.io/github/KarrLab/de_sim).

## Documentation convention

`DE-Sim` is documented using [reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html) and the [napoleon Sphinx plugin](https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html). The documentation can be compiled with [Sphinx](https://www.sphinx-doc.org/) by running `sphinx-build docs docs/_build/html`. The compiled documentation is available at [https://docs.karrlab.org/de_sim/](https://docs.karrlab.org/de_sim/).

## Submitting changes

Please use GitHub pull requests to submit changes. Each request should include a brief description of the new and/or modified features.

## Releasing and deploying new versions

Contact [info@karrlab.org](mailto:info@karrlab.org) to request release and deployment of new changes. 

## Reporting issues

Please use [GitHub issues](https://github.com/KarrLab/de_sim/issues) to report any issues to the development community.

## Getting help

Please use [GitHub issues](https://github.com/KarrLab/de_sim/issues) to post questions or contact the lead developers at [info@karrlab.org](mailto:info@karrlab.org) and [Arthur Goldberg](mailto:Arthur.Goldberg@mssm.edu).
---
title: 'DE-Sim: an object-oriented, discrete-event simulation tool for data-intensive modeling of complex systems in Python'
tags:
  - dynamical modeling
  - simulation
  - discrete-event simulation
  - object-oriented simulation
  - parallel discrete-event simulation
  - biochemical modeling
  - whole-cell modeling
  - Python
authors:
  - name: Arthur P. Goldberg
    orcid: 0000-0003-2772-1484
    affiliation: "1"
  - name: Jonathan R. Karr
    orcid: 0000-0002-2605-5080
    affiliation: "1"
affiliations:
 - name: Icahn Institute for Data Science and Genomic Technology and Department of Genetics and Genomic Sciences, Icahn School of Medicine at Mount Sinai, New York, NY 10029, USA
   index: 1
date: 14 November 2020
bibliography: paper.bib
---

# Summary

Recent advances in data collection, storage, and sharing have created unprecedented opportunities to gain insights into complex systems such as the biochemical networks that generate cellular behavior. 
Understanding the behavior of such systems will likely require larger and more comprehensive dynamical models that are based on a combination of first principles and empirical data.
These models will likely represent each component and interaction using mechanistic approximations that are derived from first principles and calibrated with data. For example, dynamical models of biochemical networks often represent the interactions among molecules as chemical reactions whose rates are determined by combining approximations of chemical kinetics and empirically-observed reaction rates.
Furthermore, complex models that represent multiple types of components and their interactions will require diverse approximations and large, heterogeneous datasets. New tools are needed to build and simulate such data-intensive models.

One of the most promising methods for building and simulating data-intensive models is discrete-event simulation (DES). DES represents the dynamics of a system as a sequence of instantaneous events [@fishman2013discrete]. DES is used for a wide range of research, such as studying the dynamics of biochemical networks, characterizing the performance of computer networks, evaluating potential war strategies, and forecasting epidemics [@banks2005discrete].
Although multiple DES tools exist, it remains difficult to build and simulate data-intensive models.  First, it is cumbersome to create complex models with the low-level languages supported by many of the existing tools. Second, most of the existing tools are siloed from the ecosystems of data science tools that are exploding around Python and R.

To address this problem, we developed DE-Sim ([https://github.com/KarrLab/de_sim](https://github.com/KarrLab/de_sim)), an open-source, object-oriented (OO), Python-based DES tool.
DE-Sim helps researchers model complex systems by enabling them to use Python's powerful OO features to manage multiple types of components and multiple types of interactions.
By building upon Python, DE-Sim also makes it easy for researchers to employ Python's powerful data science tools, such as pandas [@mckinney2010data] and SciPy [@virtanen2020scipy], to use large, heterogeneous datasets to build comprehensive and detailed models.
We anticipate that DE-Sim will enable a new generation of models that capture systems with unprecedented breadth and depth.
For example, we are using DE-Sim to develop WC-Sim [@goldberg2020wc_sim], a multi-algorithmic simulation tool for whole-cell models [@karr2015principles; @goldberg2018emerging; @karr2012whole; @goldberg2016toward] that predict phenotype from genotype by capturing all of the biochemical activity in a cell.

Here, we describe the need for new tools for building and simulating more comprehensive and more detailed models, and outline how DE-Sim addresses this need. In addition, we summarize the strengths of DE-Sim over existing DES tools, and we report the simulation performance of DE-Sim. Finally, we outline our plans to increase the performance of simulations executed by DE-Sim. 
A tutorial that describes how to build and simulate models with DE-Sim, examples, and documentation are available online, as described in the 'Availability of DE-Sim' section below.

# Statement of Need

Many scientific fields can now collect detailed data about the components of complex systems and their interactions. For example, deep sequencing has dramatically increased the availability of molecular data about biochemical networks. Combined with advances in computing, we believe that it is now possible to use this data and first principles to create comprehensive and detailed models that can provide new insights into complex systems. For example, deep sequencing and other molecular data can be used to build whole-cell models.

Achieving such comprehensive and detailed models will likely require integrating disparate principles and diverse data. While there are several DES tools, such as SimEvents [@clune2006discrete] and SimPy [@matloff2008introduction], and numerous tools for working with large, heterogeneous datasets, such as pandas and SQLAlchemy [@bayer2020sqlalchemy], it is difficult to use these tools in combination. As a result, despite having most of the major ingredients, it remains difficult to build and simulate data-intensive models.

# DE-Sim provides critical features for building and simulating data-intensive models

DE-Sim simplifies the construction and simulation of *discrete-event models* through several features. First, DE-Sim structures discrete-event models as OO programs [@zeigler1987hierarchical]. This structure enables researchers to use *simulation object* classes to encapsulate the complex logic required to represent *model components*, and use *event message* classes to encapsulate the logic required to describe the *interactions* among model components. With DE-Sim, users define simulation object classes by creating subclasses of DE-Sim's simulation object class. DE-Sim simulation object classes can exploit all the features of Python classes. For example, users can encode relationships between the types of components in a model into hierarchies of subclasses of simulation objects. As a concrete example, a model of the biochemistry of RNA transcription and protein translation could be implemented using a superclass that captures the behavior of polymers and three subclasses that represent the specific properties of DNAs, RNAs, and proteins.
By representing model components as Python simulation objects, DE-Sim makes it easy to model complex systems that contain multiple types of components by defining multiple classes of simulation objects.
Users can then model arbitrarily many instances of each type of component by creating multiple instances of the corresponding simulation object class.

Second, by building on top of Python, DE-Sim enables researchers to conveniently use Python's extensive suite of data science tools to build models from heterogeneous, multidimensional datasets. For example, researchers can use tools such as ObjTables [@karr2020objtables], H5py, requests, SQLAlchemy, and pandas to access diverse data in spreadsheets, HDF5 files, REST APIs, databases, and other sources; use tools such as NumPy [@oliphant2006guide] to integrate this data into a unified model; and use tools such as SciPy and NumPy to perform calculations during simulations of models. DE-Sim also facilitates the use of Python tools to analyze simulation results.

In addition, DE-Sim provides several features to help users execute, analyze, and debug simulations:

* **Stop conditions:** DE-Sim makes it easy to terminate simulations when specific criteria are reached. Researchers can specify stop conditions as functions that return true when a simulation should conclude.
* **Results checkpointing:** The results of a simulation can be conveniently recorded by configuring periodic checkpoints of specified parts of the simulation's state.
* **Reproducible simulations:** To help researchers debug simulations, repeated executions of the same simulation with the same configuration and same random number generator seed produce the same results.
* **Space-time visualizations:** DE-Sim generates space-time visualizations of simulation trajectories (\autoref{fig:phold_space_time_plot}). These diagrams can help researchers understand and debug simulations.

![**DE-Sim can generate space-time visualizations of simulation trajectories.** 
This figure illustrates a space-time visualization of all of the events and messages in a simulation of the parallel hold (PHOLD) DES benchmark model [@fujimoto1990performance] with three simulation objects. The timeline (black line) for each object shows its events (grey dots). The blue and purple arrows illustrate events scheduled by simulation objects for themselves and other objects, respectively. The code for this simulation is available in the DE-Sim Git repository. 
\label{fig:phold_space_time_plot}](figures/phold_space_time_plot.png)

We believe that these features can simplify and accelerate the development of complex, data-intensive models.

# Comparison of DE-Sim with existing discrete-event simulation tools

Although multiple DES tools already exist, we believe that DE-Sim uniquely facilitates data-intensive modeling through a novel combination of OO modeling and support for numerous high-level data science tools. \autoref{fig:comparison} compares the features and characteristics of DE-Sim with some of the most popular DES tools.

![**Comparison of DE-Sim with some of the most popular DES tools.**
DE-Sim is the only open-source, OO DES tool based on Python.
This combination of features makes DE-Sim uniquely suitable for creating and simulating complex, data-intensive models. 
\label{fig:comparison}](figures/comparison.png)

SimPy is an open-source DES tool that enables users to write functions that describe simulation processes (SimPy's analog to DE-Sim's simulation objects). As another Python-based tool, SymPy also makes it easy for researchers to leverage the Python ecosystem to build models. However, we believe that DE-Sim makes it easier for researchers to build complex models by enabling them to implement models as collections of classes rather than collections of functions.
DE-Sim thereby enables modelers to use a Python object to encapsulate the state of a model component together with operations on the state, and use inheritance to share state and operations among related types of model components.
In addition, we believe that DE-Sim is simpler to use because DE-Sim supports a uniform approach for scheduling events, whereas SimPy simulation processes must use two different approaches: one to schedule events for themselves, and another to schedule events for other processes.

SimEvents is a library for DES within the MATLAB/Simulink environment. While SimEvents' graphical interface makes it easy to create simple models, we believe that DE-Sim makes it easier to implement more complex models. By facilitating use of the many Python-based data science tools, DE-Sim makes it easier to use data to create models than SimEvents, which builds on a smaller ecosystem of data science tools.

SystemC is a `C++`-based OO DES tool that is frequently used to model digital systems [@mueller2001simulation]. While SystemC provides many of the same core features as DE-Sim, we believe that DE-Sim is more accessible to researchers than SystemC because DE-Sim builds upon Python, which is more popular than `C++` in many fields of research.

SIMSCRIPT III [@rice2005simscript] and SIMUL8 [@concannon2003dynamic] are commercial DES tools that define proprietary languages which researchers can use to implement models. SIMSCRIPT III is a general-purpose simulation language designed for modeling decision support systems in domains such as war-gaming, transportation, and manufacturing.
We believe that DE-Sim is more powerful than SIMSCRIPT III for most scientific and engineering problems because it leverages Python's robust data science ecosystem.

SIMUL8 models business processes as workflows.
It provides a powerful GUI for describing the flow of *work items* through a network of queues and servers, and includes tools to analyze and visualize the potentially stochastic behavior of a process.
DE-Sim is more suitable than SIMUL8 for modeling scientific or engineering systems because modelers can use DE-Sim to describe processes that cannot be easily structured as workflows.

# Performance of DE-Sim

\autoref{fig:performance} illustrates the performance of DE-Sim simulating a model of a cyclic messaging network over a range of network sizes. A messaging network consists of a ring of nodes.
When a node handles an event, it schedules the same type of event for its forward neighbor with a one time-unit delay.
Each simulation is initialized by sending a message to each node at the first time-unit. 
The code for this performance test is available in the DE-Sim Git repository, and in a Jupyter notebook at [https://sandbox.karrlab.org/tree/de_sim](https://sandbox.karrlab.org/tree/de_sim/4.%20DE-Sim%20performance%20test.ipynb).

![**Performance of DE-Sim simulating a cyclic messaging network over a range of sizes.**
Each simulation was executed for 100 time-units. Each statistic represents the average of three simulation runs in a Docker container running on a 2.9 GHz Intel Core i5 processor.
\label{fig:performance}](figures/performance.png)

# Conclusion

In summary, DE-Sim is an open-source, object-oriented, discrete-event simulation tool implemented in Python that makes it easier for modelers to create and simulate complex, data-intensive models. First, DE-Sim enables researchers to conveniently use Python's OO features to manage multiple types of model components and interactions among them. Second, DE-Sim enables researchers to directly use Python data science tools, such as pandas and SciPy, and large, heterogeneous datasets to construct models. Together, we anticipate that DE-Sim will accelerate the construction and simulation of unprecedented models of complex systems, leading to new scientific discoveries and engineering innovations.

To further advance the simulation of data-intensive models, we aim to improve the simulation performance of DE-Sim. One potential direction is to use DE-Sim as a specification language for a parallel DES system such as ROSS [@carothers2000ross]. This combination of DE-Sim and ROSS would enable modelers to both create models with DE-Sim's high-level model specification semantics and quickly simulate models with ROSS.

# Availability of DE-Sim

DE-Sim is freely and openly available under the MIT license at the locations below.

* Source code repository: [GitHub: KarrLab/de_sim](https://github.com/KarrLab/de_sim/)
* Jupyter notebook tutorials: [https://sandbox.karrlab.org/tree/de_sim](https://sandbox.karrlab.org/tree/de_sim)
* Documentation: [docs.karrlab.org/de_sim](https://docs.karrlab.org/de_sim/)

DE-Sim requires [Python](https://www.python.org/) 3.7 or higher and [pip](https://pip.pypa.io/). This article discusses version 1.0.5 of DE-Sim.

# Acknowledgements

We thank Yin Hoon Chew for her helpful feedback. This work was supported by the National Science Foundation [1649014 to JRK], the National
Institutes of Health [R35GM119771 to JRK], and the Icahn Institute for Data Science and Genomic Technology.

# References

Getting started
===============

The following examples and tutorials illustrate how to use *DE-Sim* to build and simulate models.

-----------------------------------
Examples
-----------------------------------

* `Minimal simulation <https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/minimal_simulation.py>`_: a minimal example of a simulation
* `Random walk <https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/random_walk.py>`_: a one-dimensional random walk model, with random times between steps
* `Parallel hold (PHOLD) <https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/phold.py>`_: a model developed by Richard Fujimoto to benchmark parallel discrete-event simulators
* `Epidemic <https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/sirs.py>`_: two SIR models of an infectious disease epidemic

These examples have corresponding unit tests which run them in the *DE-Sim*'s `directory of unit tests of examples <https://github.com/KarrLab/de_sim/tree/master/tests/examples>`_.

-----------------------------------
Interactive tutorials
-----------------------------------

Please see `sandbox.karrlab.org <https://sandbox.karrlab.org/tree/de_sim>`_ for interactive Jupyter notebook tutorials about designing, building and simulating models with *DE-Sim*.
It includes tutorials that use the random walk, PHOLD, and epidemic models listed above.

-----------------------------------
Template for models and simulations
-----------------------------------

The minimal simulation, located at `de_sim/examples/minimal_simulation.py <https://github.com/KarrLab/de_sim/blob/master/de_sim/examples/minimal_simulation.py>`_, can be used as a template for implementing and simulating a model with *DE-Sim*.
Performance
===========

Please see Arthur P. Goldberg & Jonathan Karr. (2020). `DE-Sim: an object-oriented, discrete-event simulation tool for data-intensive modeling of complex systems in Python. Journal of Open Source Software, 5(55), 2685. <https://doi.org/10.21105/joss.02685>`_ for information about the performance of *DE-Sim*.Comparison to other DES tools
=============================

Please see Arthur P. Goldberg & Jonathan Karr. (2020). `DE-Sim: an object-oriented, discrete-event simulation tool for data-intensive modeling of complex systems in Python. Journal of Open Source Software, 5(55), 2685. <https://doi.org/10.21105/joss.02685>`_ for a comparison of *DE-Sim* with other DES tools.Installation
============

Prerequisites
--------------------------

* Python >= 3.7
* Pip >= 19

Latest release from PyPI
---------------------------
Run the following command to install the latest release from PyPI::

    pip install de_sim

Latest release from DockerHub
-----------------------------
Run the following command to install a Docker image with the latest release from DockerHub::

    docker pull karrlab/de_sim

Latest revision from GitHub
---------------------------
Run the following command to install the latest version from GitHub::

    pip install git+https://github.com/KarrLab/de_sim.git#egg=de_sim
*DE-Sim* documentation
======================

*DE-Sim* is an open-source, Python-based, object-oriented discrete-event simulation tool that helps modelers model complex systems.
First, *DE-Sim* enables them to use Python's powerful object-oriented features to manage multiple types of components in a complex system and multiple types of interactions between these components.
Second, by building upon Python, DE-Sim makes it easy for modelers to use Python's powerful data science tools, such as `NumPy <https://numpy.org/>`_, `Scipy <https://scipy.org/scipylib/index.html>`_, `pandas <https://pandas.pydata.org/>`_, and `SQLAlchemy <https://www.sqlalchemy.org/>`_,
to incorporate large, heterogeneous datasets into comprehensive and detailed models.
We anticipate that DE-Sim will enable a new generation of models that capture systems with unprecedented breadth and depth.

*DE-Sim* provides the following features to help users build and simulate complex, data-intensive models:

* **High-level, object-oriented modeling:** *DE-Sim* facilitates model designs that use classes of *simulation objects* to encapsulate the complex logic required to represent each *model component*, and use classes of *event messages* to encapsulate the logic required to describe the *interactions* between components.
* **Powerful stop conditions:** *DE-Sim* makes it easy to terminate simulations when specific criteria are reached. Modelers can specify stop conditions as functions that return true when the simulation should conclude.
* **Results checkpointing:** Models that use *DE-Sim* can record the results of simulations, and metadata such as the start and run time of each simulation, by simply configuring a checkpointing module.
* **Space-time visualizations:** *DE-Sim* can generate space-time visualizations of simulation objects and the event messages that they exchange. These diagrams can help modelers understand and debug simulations.
* **Reproducible simulations:** To help modelers debug simulations and analyze their results, repeated executions of a simulation with the same configuration and random number generator seed produce the same results.

We have used *DE-Sim* to develop `WC-Sim <https://github.com/KarrLab/wc_sim>`_, a multi-algorithmic simulator for `whole-cell models <https://www.wholecell.org>`_.

For more information, see the `interactive DE-Sim Jupyter notebooks <https://sandbox.karrlab.org/tree/de_sim>`_ that contain a *DE-Sim* tutorial and several example *DE-Sim* models.

Contents
--------

.. toctree::
   :maxdepth: 3
   :numbered:

   getting-started.rst
   installation.rst
   API documentation <source/de_sim.rst>
   performance.rst
   comparison.rst
   about.rst
About
=====

----------------------
License
----------------------

The software is released under the MIT license:

.. literalinclude:: ../LICENSE
    :language: text

---------------
Citing *DE-Sim*
---------------

Please use the following reference to cite *DE-Sim*:

Arthur P. Goldberg & Jonathan Karr. (2020). `DE-Sim: an object-oriented, discrete-event simulation tool for data-intensive modeling of complex systems in Python. Journal of Open Source Software, 5(55), 2685. <https://doi.org/10.21105/joss.02685>`_

------------------------
Contributing to *DE-Sim*
------------------------

We enthusiastically welcome contributions to *DE-Sim*! Please see the `guide to contributing <https://github.com/KarrLab/de_sim/blob/master/CONTRIBUTING.md>`_ and the `developer's code of conduct <https://github.com/KarrLab/de_sim/blob/master/CODE_OF_CONDUCT.md>`_.

----------------------
Development team
----------------------

This package was developed by the `Karr Lab <https://www.karrlab.org/>`_ at the Icahn School of Medicine at Mount Sinai in New York, USA by the following individuals: 

* `Arthur Goldberg <https://www.mountsinai.org/profiles/arthur-p-goldberg>`_
* `Jonathan Karr <https://www.karrlab.org>`_

----------------------
Acknowledgements
----------------------

This work was supported by National Science Foundation award 1649014, National Institutes of Health award R35GM119771, and the Icahn Institute for Data Science and Genomic Technology.

----------------------
Questions and comments
----------------------

Please submit questions and issues to `GitHub <https://github.com/KarrLab/de_sim/issues>`_ or contact the `Karr Lab <mailto:info@karrlab.org>`_.
References
==========

.. bibliography:: references.bib
    :encoding: latin
    :style: unsrt
