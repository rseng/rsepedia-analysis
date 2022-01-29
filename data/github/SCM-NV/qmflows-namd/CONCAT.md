# Change Log

# 0.12.2 (17/11/2021)

## Changed
* Only distribute source files on pypi.


# 0.12.1 (17/11/2021)

## Fixed
* Fixed a number of pypi classifiers.


# 0.12.0 (17/11/2021)

## New
* Support Python 3.9 and 3.10
* Allow to compute the spectrum of multiple stack geometries (324)

## Changed
* Check for duplicate keys when loading .yaml files.
* Make cell parameters optional.

## Fixed
* Various fixes.


# 0.11.0 (04/12/2020)
## New
* Print CP2K err/out files if the calculation fails (#150)

## Changed
* HDF5 storage layout (#300)


# 0.10.4 (26/10/2020)
## New
* Allow to compute both alphas/betas derivative couplings simultaneusly (#275)
* Allow to compute mutipole matrices for unrestricted calculation (#297, #299)

## Changed
* Do not remove the CP2K log files by default
* Do not remove the point folder wher ethe CP2K orbitals are stored

## Fixed
* Unrestricted Hamiltitonians name (#286)
* Hamiltonian units (#290)
* Schema error (#292)
* Distribute Slurm to retrieve Hamiltonians (#296)


# 0.10.3 (09/10/2020)
## New
* Template to create B3LYP computations (#269)
* Add support for derivative couplings for system with more than one spin state (#275)

## Fixed
* Distribution error (#272)
* Molecular orbital error [in qmflows](https://github.com/SCM-NV/qmflows/pull/213) (#270)
* Multivalue settings issue (#260)
* CP2K executable (#264)

# 0.10.1
## New
* Keywords to print eigenvalues and eigenvectors (#248)

## Fixed
* SLURM free-format specification (#245)
* Touch HDF5 file if it doesn't exist (#246)
* Create a separate folder for each distributed chunk (#247)
* Error creating the scratch path and the HDF5 file (#255)

# 0.10.0
## Changed
* Rename package to **nano-qmflows**

# 0.9.0
## Added
* [Support for Mac](https://github.com/SCM-NV/nano-qmflows/issues/231)
* [Allow to specify the CP2K executable](https://github.com/SCM-NV/nano-qmflows/issues/226)
* [MyPy static checking](https://github.com/SCM-NV/nano-qmflows/issues/237)

## Changed
* [Use New QMFlows API](https://github.com/SCM-NV/nano-qmflows/issues/227)
* [Use Libint==2.6.0](https://github.com/SCM-NV/nano-qmflows/issues/234)
* [Allow the user to enter her own slurm script](https://github.com/SCM-NV/nano-qmflows/issues/225)

# 0.8.3

## Changed
* Add the global run_type  keyword in the templates

# 0.8.2 [31/01/20]

## Changed

* Replace `qmflows.utils` with [more-itertools](https://more-itertools.readthedocs.io/en/stable/index.html)

## Added
* [smiles_calculator](https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/qmflows/smiles_calculator.py) script to compute molecular properties from smiles.

# 0.8.1 [17/10/19]

## Changed
* Use [f-strings](https://docs.python.org/3/reference/lexical_analysis.html#f-strings)
* Clean [C++ interface](https://cgithub.com/SCM-NV/nano-qmflows/blob/master/libint/compute_integrals.cc) to [libint](https://github.com/evaleev/libint)

## Removed
* Unused code to compile the [C++ interface](https://cgithub.com/SCM-NV/nano-qmflows/blob/master/libint/compute_integrals.cc) to [libint](https://github.com/evaleev/libint)

# 0.8.0

### Fixed

* Passed to libint2 the Cartesian coordinates in Angstrom instead *atomic units*.


# 0.7.0

### New

 * Allow to compute charge molecules in the *C2Pk* input.
 * Compute the multipole integrals in the center of mass.
 * A new variable called ``aux_fix`` has been introduced to change the quality of the auxiliar basis set
   for hybrid calculation. The possible values are: "low", "medium", "good", "verygood" and "excellent".
   The default value is: verygood.
 * Return a ``input_parameters.yml`` file containing the input after all the preprocessing steps.

## Change

 * The ``path_basis`` variable in the yaml input, points to the folder where all the CP2K basis are located.
   By Default this variable points to <Installation>/nac/basis where there are some default basis.

### Deleted

* The ``path_potential`` variable has been removed since it is superseded by the ``path_basis``.


# 0.6.0

### New
 * Compute the overlap integrals to calculate the derivative coupling and the multipole integrals using [libint2](https://github.com/evaleev/)
 * Used `openmp` to compute the integrals in all the available cores
 * New dependencies: [eigen](http://eigen.tuxfamily.org/dox/), [highfive](https://github.com/BlueBrain/HighFive/tree/master/include/highfive), [libint2](https://github.com/evaleev/libint/wiki) and [pybind11](https://pybind11.readthedocs.io/en/master/)

### Deleted

 * Python/Cython implementation of the overlap integrals
 * Unused functionality replaced by [libint2](https://github.com/evaleev/)

# 0.5.0

### New

* The user only need to provide an **active_space** and both the `mo_index_range` and `nHOMO`  keywords are computed automatically.

* Added fast test to [compute the couplings](https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_coupling.py)

### Deleted

* Removed all the Fourier transform for orbitals.

* Removed unused electron transfer functionality.

### Changed

* The `nHOMO` and the `kinds` for the *CP2K* input are computed using the [valence_electrons](https://github.com/SCM-NV/nano-qmflows/blob/master/nac/basisSet/valence_electrons.json) from the basis/pseudpotential combination.

* Use a configuration dictionary to around the initial input instead of many arguments functions.

* Import only the most important functions.


# 0.4.1

### Deleted

* Removed all the MPI unused functionality

### Changed

* Refactor the [distribute_jobs.py](https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/distribution/distribute_jobs.py) script to split the derivative coupling calculations.

# 0.4.0

### Deleted

* removed `workflow_oscillator_strength`. Use `workflow_stddft` instead

### Changed

* Moved `nHomo` keyword to `general_setting`
* Renamed the `ci_range` keyword and replaced it by the *CP2K* keyword `mo_index_range`

### New
* Templates to call functionals **pbe** and **pbe0** to compute the Molecular orbitals


## 0.3.1

### Changed

* Replace the `json schemas` with the [schemas](https://github.com/keleshev/schema) library


## 0.3.0

### Added

The following actions were performed:
* Removed nose and pandas dependencies
* Use pytest for testing
* Replace MIT license by Apache-2.0
* Allow only fast tests in Travis
* Added changelog
* made general mergeHDF5 script
* Added Runners: MPI and Multiprocessing(default)
* Introduce new input file (yaml)
* Validate input files with json schemas
* Refactor the workflow API
* Used [noodles==0.3.1](https://github.com/NLeSC/noodles) and [qmflows==0.3.0](https://github.com/SCM-NV/qmflows)


### Removed

* Dead code from `workflow_cube`
 Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

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
reported by contacting the project team at f.zapata@esciencecenter.nl. All
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
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

## You have a question

1. use the search functionality [here](https://github.com/SCM-NV/nano-qmflows/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality [here](https://github.com/SCM-NV/nano-qmflows/issues) to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the [SHA hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas) of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community _before you start working_. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
1. make sure the existing tests still work by running ``pytest test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. [push](http://rogerdudler.github.io/git-guide/) your feature branch to (your fork of) the [nano-qmflows](https://github.com/SCM-NV/nano-qmflows) repository on GitHub;
1. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
---
name: Bug report
about: Something doesn't work like I expected.
title: ''
labels: bug
assignees: ''

---

Use your best judgment to provide a useful level of information. Depending on the nature of the issue, consider including, e.g.

- Which version of the software you're running
- Any logging or error output you see
---
name: Help wanted
about: I need some help with the code
title: ''
labels: help wanted
assignees: ''

---
---
name: Feature request
about: I would like a new feature to be included in the library
title: ''
labels: enhancement
assignees: ''

---

Tell us what would be a nice feature to make the **nano-qmflows** library even better!
### All Submissions:

* [ ] Have you followed the guidelines in our Contributing document?
* [ ] Have you checked to ensure there aren't other open [Pull Requests](https://github.com/SCM-NV/nano-qmflows/pulls) for the same update/change?


### New Feature Submissions:

1. [ ] Does your submission pass tests?
2. [ ] Have you check your code quality (e.g. using [flake8](http://flake8.pycqa.org/en/latest/)) prior to submission?

### Changes to Core Features:

* [ ] Have you added an explanation of what your changes do and why you'd like us to include them?
* [ ] Have you written new tests for your core changes, as applicable?
* [ ] Have you successfully ran tests with your changes locally?
