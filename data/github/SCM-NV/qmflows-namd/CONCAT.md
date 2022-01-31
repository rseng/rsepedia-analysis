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
.. image:: https://readthedocs.org/projects/qmflows-namd/badge/?version=latest
   :target: https://qmflows-namd.readthedocs.io/en/latest/?badge=latest
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3661036.svg
   :target: https://doi.org/10.5281/zenodo.3661036
.. image:: https://github.com/SCM-NV/nano-qmflows/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/SCM-NV/nano-qmflows/actions
.. image:: https://codecov.io/gh/SCM-NV/nano-qmflows/branch/master/graph/badge.svg?token=L1W0fPrSUn
   :target: https://codecov.io/gh/SCM-NV/nano-qmflows

====================
nano-qmflows
====================

Nano-QMFlows is a generic python library for computing (numerically) electronic properties for nanomaterials like the non-adiabatic coupling vectors (NACV) using several quantum chemical (QM) packages.

One of the main problems to calculate (numerically) NACVs by standard QM software is the computation of the overlap matrices between two electronically excited states at two consecutive time-steps that are needed in the numerical differentiation to evaluate the coupling. This happens because most of these softwares are inherently static, i.e. properties are computed for a given structural configuration, and the computation of the overlap matrices at different times requires complicated scripting tools to handle input/outputs of several QM packages.

For further information on the theory behind nano-qmflows and how to use the program see the documentation_.

Installation
------------

In order to install the **nano-qmflows** library you need to install *Miniconda* as detailed here_.

.. _here: https://docs.conda.io/en/latest/miniconda.html

Then,  to install the **nano-qmflows** library type the following commands inside the conda environment:
  - ``conda create -n qmflows -c conda-forge boost eigen libint==2.6.0 highfive``
  - ``conda activate qmflows``
  - ``pip install nano-qmflows --upgrade``

.. note::
   For GCC <7 one has to pass ``eigen=3.3``.

Advantages and Limitations
--------------------------
nano-qmflows is based on the approximation that all excited states are represented by singly excited-state determinants. This means that the computation of the NACVs boils down to the computation of molecular orbitals (MOs) coefficients at given points of time using an electronic structure code and an overlap matrix S(t,t+dt) in atomic orbital basis (AO) computed between two consecutive time step. nano-qmflows main advantage is to use an internal module to compute efficiently the atomic overlap matrix S(t, t+dt) by employing the same basis-set used in the electronic structure calculation. In this way the QM codes are only needed to retrieve the MOs coefficients at time t and t+dt. This approach is very useful because the interfacing nano-qmflows to a QM code is reduced to writing a simple module that reads the MOs coefficients in the specific code format. At this moment, nano-qmflows handles output formats generated by CP2K, Orca, and Gamess, but, as said, it can be easily extended to other codes.

Finally, nano-qmflows can be also used in benchmarks studies to test new code developments in the field of excited state dynamics by providing a platform that uses all the functionalities of QMFlows, which automatizes the input preparation and execution of thousands of QM calculations.

In the near future, nano-qmflows is expected to offer new functionalities.


Interface to Pyxaid
-------------------

nano-qmflows has been designed mostly to be integrated with Pyxaid, a python program that performs non-adiabatic molecular dynamic (NAMD) simulations using the classical path approximation (CPA). The CPA is based on the assumption that nuclear dynamics of the system remains unaffected by the dynamics of the electronic degrees of freedom. Hence, the electronic dynamics remains driven by the ground state nuclear dynamics. CPA is usually valid for extended materials or cluster materials of nanometric size.

In this framework, nano-qmflows requires as input the coordinates of a pre-computed trajectory (at a lower level or at the same level of theory) in xyz format and the input parameters of the SCF code (HF and DFT). nano-qmflows will then calculate the overlap matrix between different MOs by correcting their phase and will also track the nature of each state at the crossing seam using a min-cost algorithm . The NACVs are computed using the Hammes-Schiffer-Tully (HST) 2-point approximation and the recent Meek-Levine approach. The NACVs are then written in Pyxaid format for subsequent NAMD simulations.


Overview
--------
 The Library contains a **C++** interface to the libint2_ library to compute the integrals and several numerical functions in Numpy_. While the scripts are set of workflows to compute different properties using different approximations that can be tuned by the user.

.. _libint2: https://github.com/evaleev/libint/wiki
.. _Numpy: http://www.numpy.org

Worflow to calculate Hamiltonians for nonadiabatic molecular simulations
************************************************************************
The figure represents schematically a Worflow to compute the **Hamiltonians** that described the behavior and coupling between the excited state of a molecular system. These **Hamiltonians** are used by thy PYXAID_ simulation package to carry out nonadiabatic molecular dynamics.

.. image:: docs/_images/nac_worflow.png

.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _documentation: https://qmflows-namd.readthedocs.io/en/latest/
Using coordination_ldos.py
--------------------------

The script prints local PDOS projected on subsets of atoms given through lists.
These lists are obtained using the nano-CAT module ``nanoCAT.recipes.coordination_number`` (see the relative documentation_)
that returns a nested dictionary

``{'Cd': {4: [0, 1, 2, 3, 4, ...], ...}, ...}``

with atomic symbol (*e.g.* ``'Cd'``) and coordination number (*e.g.* ``4``) as keys.


You thus have to install the nano-CAT package in your conda environment according to the installation instructions, reported here_.

.. _documentation: https://cat.readthedocs.io/en/latest/12_5_recipes.html
.. _here: https://github.com/nlesc-nano/nano-CAT
Derivative Couplings
--------------------
.. automodule:: nanoqm.integrals.nonAdiabaticCouplingCP2K Interface
--------------
.. automodule:: nanoqm.schedule.scheduleCP2KIntegrals
---------
.. automodule:: nanoqm.integrals.multipole_matricesCommand line interface
----------------------
Running a workflow
##################
.. automodule:: nanoqm.workflows.run_workflow

Workflows distribution
######################
.. automodule:: nanoqm.workflows.distribute_jobsWorkflows
---------

The following workflows are available:

.. autofunction:: nanoqm.workflows.workflow_coop.workflow_crystal_orbital_overlap_population
.. autofunction:: nanoqm.workflows.workflow_coupling.workflow_derivative_couplings
.. autofunction:: nanoqm.workflows.workflow_ipr.workflow_ipr
.. autofunction:: nanoqm.workflows.workflow_single_points.workflow_single_points
.. autofunction:: nanoqm.workflows.workflow_stddft_spectrum.workflow_stddft
Crystal Orbital Overlap Population (COOP) calculation
=====================================================

The workflow coop_calculation allows to compute the crystal orbital overlap population between two selected elements.

Preparing the input
-------------------

The following is an example of input file to perform the COOP calculation between Cd and Se for the Cd33Se33 system.

.. code-block:: yaml

    workflow:
      coop_calculation

    project_name: Cd33Se33
    active_space: [50, 50]
    path_hdf5: "Cd33Se33.hdf5"
    path_traj_xyz: "Cd33Se33.xyz"
    scratch_path: "/tmp/COOP"

    coop_elements: ["Cd", "Se"]

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 28.0
      periodic: none
      executable: cp2k.popt
      
      cp2k_settings_main:
        specific:
          template: pbe_main
          cp2k:
            force_eval:
              dft:
                scf:
                  eps_scf: 1e-6
 
      cp2k_settings_guess:
        specific:
          template:
            pbe_guess


In your working directory, copy the previous input into an *input_test_coop.yml* file. 
Also copy locally the file containing the coordinates of the relaxed Cd33Se33 system, Cd33Se33.xyz_.

Your *input_test_coop.yml* input file now contains all settings to perform the coop calculations and needs to be edited according to your system and preferences.
Please note that this input is very similar to the basic example of single point calculation provided in a previous tutorial_ (please refer to it for a more extensive description of the above options)
except for the following options: **workflow**, **coop_elements**.

- **workflow**: The workflow you need for your calculations, in this case set to coop_calculation is this case.
- **coop_elements**: List of the two elements to calculate the COOP for, here Cd and Se.

In the cp2k_general_settings, you can customize the settings used to generate the cp2k input. To help you creating your custom input requirements, please consult the cp2k manual_ and the templates_ available in nano-qmflows.

.. _Cd33Se33.xyz: https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_files/Cd33Se33.xyz
.. _tutorial: https://qmflows-namd.readthedocs.io/en/latest/single_points.html
.. _manual: https://manual.cp2k.org/
.. _templates: https://github.com/SCM-NV/nano-qmflows/blob/master/nanoqm/workflows/templates.py

Setting up the calculation 
---------------------------

Once all settings of your yml input have been customized, can to launch your coop calculation.

- First, activate the conda environment with QMFlows:

  ``conda activate qmflows``
  
- Then, load the module with your version of cp2k, for example:

  ``module load CP2K/7.1.0``
  
- Finally, use the command run_workflow.py to submit your calculation.

  ``run_workflow.py -i input_test_coop.yml``

Results 
-------

Once your calculation has finished successfully, you will find a *COOP.txt* file in your working directory.
The two columns of this file contain, respectively, the orbitals’ energies and the corresponding COOP values for the selected atoms pair.
Theory
==========

Nonadiabatic coupling matrix
-----------------------------

The current implementation of the nonadiabatic coupling is based on:
Plasser, F.; Granucci, G.; Pittner, j.; Barbatti, M.; Persico, M.;
Lischka. *Surface hopping dynamics using a locally diabatic formalism:
Charge transfer in the ethylene dimer cation and excited state dynamics
in the 2-pyridone dimer*. **J. Chem. Phys. 2012, 137, 22A514.**

The total time-dependent wave function :math:`\Psi(\mathbf{R}, t)` can be
expressed in terms of a linear combination of ``N`` adiabatic electronic
eigenstates :math:`\phi_{i}(\mathbf{R}(t))`,

.. math::
   \Psi(\mathbf{R}, t) = \sum^{N}_{i=1} c_i(t)\phi_{i}(\mathbf{R}(t)) \quad \mathbf(1)

The time-dependent coefficients are propagated according to

.. math::
   
   \frac{dc_j(t)}{dt} = -i\hbar^2 c_j(t) E_j(t) - \sum^{N}_{i=1}c_i(t)\sigma_{ji}(t) \quad \mathbf(2)

where :math:`E_j(t)` is the energy of the jth adiabatic state and :math:`\sigma_{ji}(t)` the nonadiabatic matrix, which elements are given by the expression

.. math::
  \sigma_{ji}(t) = \langle \phi_{j}(\mathbf{R}(t)) \mid \frac{\partial}{\partial t} \mid \phi_{i}(\mathbf{R}(t)) \rangle \quad \mathbf(3)

that can be approximate using three consecutive molecular geometries

.. math::
   \sigma_{ji}(t) \approx \frac{1}{4 \Delta t} (3\mathbf{S}{ji}(t) - 3\mathbf{S}{ij}(t) - \mathbf{S}{ji}(t-\Delta t) + \mathbf{S}{ij}(t-\Delta t)) \quad \mathbf(4)

where :math:`\mathbf{S}_{ji}(t)` is the overlap matrix between two consecutive time steps

.. math::
   \mathbf{S}{ij}(t) = \langle \phi{j}(\mathbf{R}(t-\Delta t)) \mid \phi_{i}(\mathbf{R}(t)) \rangle \quad \mathbf(5)

and the overlap matrix is calculated in terms of atomic orbitals

.. math::
   \mathbf{S}{ji}(t) = \sum{\mu} C^{*}{\mu i}(t) \sum{\nu} C_{\nu j}(t - \Delta t) \mathbf{S}_{\mu \nu}(t) \quad \mathbf(6)

Where :math:C_{\mu i} are the Molecular orbital coefficients and :math:`\mathbf{S}_{\mu \nu}` The atomic orbitals overlaps.

.. math::
   \mathbf{S}{\mu \nu}(\mathbf{R}(t), \mathbf{R}(t - \Delta t)) = \langle \chi{\mu}(\mathbf{R}(t)) \mid \chi_{\nu}(\mathbf{R}(t - \Delta t)\rangle \quad \mathbf(7)


Nonadiabatic coupling algorithm implementation
----------------------------------------------

The  figure belows shows schematically the workflow for calculating the Nonadiabatic 
coupling matrices from a molecular dynamic trajectory. The uppermost node represent
a molecular dynamics
trajectory that is subsequently divided in its components andfor each geometry the molecular
orbitals are computed. These molecular orbitals are stored in a HDF5_.
binary file and subsequently calculations retrieve sets of three molecular orbitals that are
use to calculate the nonadiabatic coupling matrix using equations **4** to **7**.
These coupling matrices are them feed to the PYXAID_ package to carry out nonadiabatic molecular dynamics.

The Overlap between primitives are calculated using the Obara-Saika recursive scheme and has been implemented using the C++ libint2_ library for efficiency reasons.
The libint2_ library uses either OpenMP_ or C++ threads to distribute the integrals among the available CPUs.
Also, all the heavy numerical processing is carried out by the highly optimized functions in NumPy_.

 The **nonadiabaticCoupling** package relies on *QMWorks* to run the Quantum mechanical simulations using the [CP2K](https://www.cp2k.org/) package.  Also, the noodles_ is used
 to schedule expensive numerical computations that are required to calculate the nonadiabatic coupling matrix.


.. _OpenMP: https://www.openmp.org/
.. _libint2: https://github.com/evaleev/libint/wiki
.. _HDF5: http://www.h5py.org/
.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _Cython: http://cython.org
.. _multiprocessing: https://docs.python.org/3.6/library/multiprocessing.html
.. _NumPy: http://www.numpy.org
.. _noodles: http://nlesc.github.io/noodles/

For a more detailed description of **nano-qmflows** read the documentation

.. toctree::
   docs_command_line
   docs_cp2k_interface
   docs_derivative_coupling
   docs_molecular_orbitals
   docs_integrals
   docs_workflowsInverse Participation Ratio (IPR) calculation
=============================================

The workflow ipr_calculation returns the inverse participation ratio for the selected orbitals. 
For finite systems, the IPR is defined as the inverse of number of atoms that contribute to a given electronic state i. 
It assumes its maximum value, 1, in the case of a state localized to a single atom (1/1) and tends to 0 (1/*N*, where *N* is the total number of atoms in the system) when the wave function is distributed equally over all atoms.

Preparing the input
-------------------

The following is an example of input file to perform the IPR calculation for the Cd33Se33 system.

.. code-block:: yaml

    workflow:
       ipr_calculation

    project_name: Cd33Se33
    active_space: [50, 50]
    path_hdf5: "Cd33Se33.hdf5"
    path_traj_xyz: "Cd33Se33.xyz"
    scratch_path: "/tmp/IPR"

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 28.0
      periodic: none
      executable: cp2k.popt
    
      cp2k_settings_main:
        specific:
          template: pbe_main
          cp2k:
            force_eval:
              dft:
                scf:
                  eps_scf: 1e-6
    
      cp2k_settings_guess:
        specific:
          template:
            pbe_guess
            

In your working directory, copy the previous input into an *input_test_ipr.yml* file. 
Also copy locally the file containing the coordinates of the relaxed Cd33Se33 system, Cd33Se33.xyz_.

Your *input_test_ipr.yml* input file now contains all settings to perform the coop calculations and needs to be edited according to your system and preferences.
Please note that this input is very similar to the basic example of single point calculation provided in a previous tutorial_ (please refer to it for a more extensive description of the above options)
except for the **workflow** option, set in this case to *ipr_calculation*.

Here again you can customize the settings used to generate the cp2k input in the cp2k_general_settings. To help you creating your custom input requirements, please consult the cp2k manual_ and the templates_ available in nano-qmflows.

.. _Cd33Se33.xyz: https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_files/Cd33Se33.xyz
.. _tutorial: https://qmflows-namd.readthedocs.io/en/latest/single_points.html
.. _manual: https://manual.cp2k.org/
.. _templates: https://github.com/SCM-NV/nano-qmflows/blob/master/nanoqm/workflows/templates.py

Setting up the calculation 
---------------------------

Once all settings of your yml input have been customized, can to launch your ipr calculation.

- First, activate the conda environment with QMFlows:

  ``conda activate qmflows``
  
- Then, load the module with your version of cp2k, for example:

  ``module load CP2K/7.1.0``
  
- Finally, use the command run_workflow.py to submit your calculation.

  ``run_workflow.py -i input_test_ipr.yml``

Results 
-------

Once your calculation has finished successfully, you will find a *IPR.txt* file in your working directory.
The two columns of this file contain, respectively, the orbitals’ energies and the corresponding IPR values.
Molecular Orbitals
------------------
.. automodule:: nanoqm.schedule.componentsDerivative coupling calculation
===============================

These tutorials focus on how to compute non-adiabatic coupling vectors between molecular orbitals belonging at two different time steps, t and t+dt, of a pre-computed molecular dynamics trajectory. What this program does is to compute at each point of the trajectory, the electronic structure using DFT, and then the overlap integrals :math:`\langle \psi_{i}(t) \mid \psi_{j}(t+dt)>`. These integrals are stored and finally used to compute numerically the non-adiabatic couplings. These and the orbital energies are written in a format readable by PYXAID to perform surface hopping dynamics. 
When using this tutorial, ensure you have the latest version of QMFlows and nano-qmflows installed.

Preparing the input
--------------------
The following is an example of the inputfile for the calculation of derivative couplings for the Cd33Se33 system. The calculations are carried out with the CP2k package, which you need to have pre-installed. The level of theory is DFT/PBE. 

.. code-block:: yaml

    workflow: distribute_derivative_couplings
    project_name: Cd33Se33
    dt: 1
    active_space: [10, 10]
    algorithm: "levine"
    tracking: False
    path_hdf5: "test/test_files/Cd33Se33.hdf5"
    path_traj_xyz: "test/test_files/Cd33Se33_fivePoints.xyz" 
    scratch_path: "/tmp/namd"
    workdir: "."
    blocks: 2

    job_scheduler:
      free_format: "
      #! /bin/bash\n
      #SBATCH --job-name=Cd33Se33\n
      #SBATCH -N 1\n
      #SBATCH -t 00:15:00\n
      #SBATCH -p short\n
      source activate qmflows\n
      module load cp2k/3.0\n\n"
    
    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 28.0
      file_cell_parameters: "test/test_files/file_distribute_cell_parameters.txt"
      periodic: none
      executable: cp2k.popt

      cp2k_settings_main:
        specific:
          template: pbe_main

      cp2k_settings_guess:
        specific:
          template: pbe_guess


The previous input can be found at input_test_distribute_derivative_couplings.yml_. Copy this file to a folder where you want start the QMFlows calculations. 

The *input_test_distribute_derivative_couplings.yml* file contains all settings to perform the calculations and needs to be edited according to your system and preferences. Pay attention to the following parameters: *project_name, dt, active_space, algorithm, tracking, path_hdf5, path_traj_xyz, scratch_path, workdir, blocks*. 

- **project_name**: Project name for the calculations. You can choose what you wish. 
- **dt**: The size of the timestep used in your MD simulations. 
- **active_space**: Range of `(occupied, virtual)` molecular orbitals to computed the derivate couplings. For example, if 50 occupied and 100 virtual should be considered in your calculations, the active space should be set to [50, 100]. 
- **algorithm**: Algorithm to calculate derivative couplings can be set to ‘levine’ or ‘3points’.
- **tracking**: If required, you can track each state over the whole trajectory. You can also disable this option.  
- **path_hdf5**: Path where the hdf5 should be created / can be found. The hdf5 is the format used to store the molecular orbitals and other information. 
- **path_traj_xyz**: Path to the pre-computed MD trajectory. It should be provided in xyz format. 
- **scratch_path**: A scratch path is required to perform the calculations. For large systems, the .hdf5 files can become quite large (hundredths of GBs) and calculations are instead performed in the scratch workspace. The final results will also be stored here.
- **workdir**: This is the location where the logfile and the results will be written. Default setting is current directory.
- **blocks**: The number of blocks (chunks) is related to how the MD trajectory is split up. As typical trajectories are quite large (+- 5000 structures), it is convenient to split the trajectory up into multiple chunks so that several calculations can be performed simultaneously. Generally around 4-5 blocks is sufficient, depending on the length of the trajectory and the size of the system. 
- **write_overlaps**: The overlap integrals are stored locally. This option is usually activated for debugging.
- **overlaps_deph**: The overlap integrals are computed between t=0 and all othe times: <psi_i (t=0) | psi_j (t + dt)>. This option is of interest to understand how long it takes to a molecular orbital to dephase from its starting configuration. This option is disabled by default. 

The **job_scheduler** can be found below these parameters. Customize these settings according to the system and environment you are using to perform the calculations. 

In the **cp2k_general_settings**, you can customize the settings used to generate the cp2k input. You can use the cp2k manual_ to create your custom input requirements. Remember to provide a path to the folder with the cp2k basis set anc potential files.

.. _manual: https://manual.cp2k.org/
.. _input_test_distribute_derivative_couplings.yml: https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_files/input_test_distribute_derivative_couplings.yml

Setting up the calculation 
---------------------------

Once all settings in *input_test_distribute_derivative_couplings.yml* have been customized, you will need to create the different chunks. 
  
- First, activate QMFlows:

  ``conda activate qmflows``  

- Use the command *distribute_jobs.py* to create the different chunks:

  ``distribute_jobs.py -i input_test_distribute_derivative_couplings.yml``

A number of new folders are created. In each folder you will find a launch.sh file, a chunk_xyz file and an input.yml file. In the input.yml file, you will find all your settings. Check for any possible manual errors.

- If you are satisfied with the inputs, submit each of your jobs for calculation.

You can keep track of the calculations by going to your scratch path. The location where all points of the chunks are calculated is your assigned scratch path plus project name plus a number. 

The overlaps and couplings between each state will be calculated once the single point calculations are finished. The progress can be tracked with the .log file in your working directory folders. The calculated couplings are meaningless at this point and need to be removed and recalculated, more on that later.  

Merging the chunks and recalculating the couplings 
---------------------------------------------------

Once the overlaps and couplings are all calculated, you need to merge the different chunks into a single chunk, as the overlaps between the different chunks still need to be calculated. For this you will use the *mergeHDF5.py* command that you will have if you have installed QMFlows correctly. 

You are free to choose your own HDF5 file name but for this tutorial we will use *chunk_01.hdf5* as an example. 

- Merge the different chunk into a single file using the *mergeHDF5.py* script:

  ``mergeHDF5.py -i chunk_0.hdf5 chunk_1.hdf5 -o chunk_01.hdf5``

Follow -i with the names of different chunks you want to merge and follow -o the name of the merged HDF5 file.  

- Remove the couplings from the chunk_01.hdf5 using the *removeHDF5folders.py* script. To run the script, use: 

  ``removeHDF5folders.py -hdf5 chunk_01.hdf5``

Using the script in this manner will only allow the couplings to be removed. 

.. Note::
   If required, you can remove all overlaps by by adding -o at the end of the previous command:

  ``removeHDF5folders.py -hdf5 chunk_01.hdf5 –o``


- Create a new subfolder in your original working directory and copy the *input.yml* file that was created for chunk 0 (when running the *distribute_jobs.py* script) to this folder. 

- Edit the *input.yml* file to include the path to the merged .hdf5, the full MD trajectory, and a new scratch path for the merged hdf5 calculations.

- Relaunch the calculation.

Once the remaining overlaps and the couplings have been calculated successfully, the hdf5 files and hamiltonians will be written to both the working directory as well as the scratch folder in a format suitable for PYXAID to run the non-adiabatic excited state molecular dynamics. If requested, also the overlap integrals can be found in the working directory.

.. note::
   There are several way to declare the parameters of the unit cell, you can passed to the cell_parameters
   variable either a number, a list or a list or list. A single number represent a cubic box, while a list
   represent a parallelepiped and finally a list of list contains the ABC vectors describing the unit cell.
   Alternatively, you can pass the angles of the cell using the cell_angles variable.

Restarting a Job
----------------

Both the *molecular orbitals* and the *derivative couplings* for a given molecular dynamic trajectory are stored in a HDF5_. The library check wether the *MO* orbitals or the coupling under consideration are already present in the HDF5_ file, otherwise compute it. Therefore  if the workflow computation fails due to a recoverable issue like:

  * Cancelation due to time limit.
  * Manual suspension or cancelation for another reasons.

Then, in order to restart the job you need to perform the following actions:

  * **Do Not remove** the file called ``cache.db`` from the current work  directory.


Reporting a bug or requesting a feature
---------------------------------------
To report an issue or request a new feature you can use the github issues_ tracker.

.. _HDF5: http://www.h5py.org/
.. _issues: https://github.com/SCM-NV/nano-qmflows/issues
.. _QMflows: https://github.com/SCM-NV/qmflows
.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _YAML: https://pyyaml.org/wiki/PyYAML


Single points calculation
=========================

The single_points workflow performs single point calculations and can be used, for example, to compute the eigenvalues and energies of triplet orbitals on a singlet geometry or viceversa.

Preparing the input
--------------------

Basic Example
^^^^^^^^^^^^^

We start with the very basic example of an input file to perform a single point calculation on the relaxed geometry of a Cd33Se33 system.

.. code-block:: yaml

    workflow:
      single_points

    project_name: Cd33Se33
    active_space: [50, 50]
    compute_orbitals: True
    path_hdf5: "Cd33Se33.hdf5"
    path_traj_xyz: "Cd33Se33.xyz"
    scratch_path: "/tmp/singlepoints_basic"

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 28.0
      periodic: none
      executable: cp2k.popt

      cp2k_settings_main:
        specific:
          template: pbe_main

      cp2k_settings_guess:
        specific:
          template:
            pbe_guess

In your working directory, create an *input_test_single_points_basic.yml* file and copy the previous input inside it, by respecting the indentation.
Also copy locally the file containing the coordinates of the Cd33Se33 system in an xyz format, Cd33Se33.xyz_.

Your *input_test_single_points_basic.yml* now contains all settings to perform the calculations and needs to be edited according to your system and preferences. 
Pay attention to the following parameters that are common to all input files: **workflow**, **project_name**, **active_space**, **path_hdf5**, **path_traj_xyz**, **scratch_path**.

- **workflow**: The workflow you need for your calculations, in this case single_points for a single point calculation.
- **project_name**: Project name for the calculations. You can choose what you wish.
- **active_space**: Range of (occupied, virtual) molecular orbitals to be computed. For example, if 50 occupied and 100 virtual should be considered in your calculations, the active space should be set to [50, 100].
- **compute_orbitals**: Specify if the energy and eigenvalues of the selected orbitals are to be computed. The default is set to True.
- **path_hdf5**: Path where the hdf5 should be created / can be found. The hdf5 is the format used to store the molecular orbitals and other information.
- **path_traj_xyz**: Path to the pre-optimized geometry of your system. It should be provided in xyz format.
- **scratch_path**: A scratch path is required to perform the calculations. For large systems, the .hdf5 files can become quite large (hundredths of GBs) and calculations are instead performed in the scratch workspace. The final results will also be stored here.

You can find the complete list of all the general options (common to all workflows) in this dictionary_.

In the cp2k_general_settings, you can customize the settings used to generate the cp2k input (see available options in schema_cp2k_general_settings_).

Here you can specify the level of theory you want to use in your cp2k calculation (basis set and potential) as well as the main characteristics of your system (cell parameters and angles, periodicity, charge, multiplicity, …). 

Note that the (fast) SCF routine in cp2k is based on the Orbital Transformation (OT) method, which works on the occupied orbital subspace. To obtain the full spectrum of molecular orbitals, one should perform a full diagonalization of the Fock matrix. For this reason, to obtain and store both occupied and unoccupied orbitals, defined using the active_space keyword, we have to follow a 2-step procedure: in the first step, which in the yaml input we define as cp2k_settings_guess, we perform a single point calculation using the fast OT approach; then in the second step, defined as cp2k_settings_main, we use the converged orbitals in the first step to start a full diagonalization calculation using the DIIS procedure.

In the cp2k_settings_guess and cp2k_settings_main subsections you can provide more detailed information about the cp2k input settings to be used to compute the guess wavefunction and to perform the main calculation respectively.
In this example, we have used one of the available templates_, specifically customized for calculations with a PBE exchange-correlation functional. 
You can use the cp2k manual_ to further personalize your input requirements.

.. _Cd33Se33.xyz: https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_files/Cd33Se33.xyz
.. _dictionary: https://github.com/SCM-NV/nano-qmflows/blob/e176ade9783677962d5146d8e6bc5dd6bb4f9102/nanoqm/workflows/schemas.py#L116
.. _schema_cp2k_general_settings: https://github.com/SCM-NV/nano-qmflows/blob/e176ade9783677962d5146d8e6bc5dd6bb4f9102/nanoqm/workflows/schemas.py#L55
.. _templates: https://github.com/SCM-NV/nano-qmflows/blob/master/nanoqm/workflows/templates.py
.. _manual: https://manual.cp2k.org/


Advanced Example
^^^^^^^^^^^^^^^^

We are now ready to move to a more advanced example in which we want to compute the orbitals' energies and eigenvalues for each point of a pre-computed MD trajectory for our Cd33Se33 system. The input file will look like that:

.. code-block:: yaml

    workflow:
      single_points

    project_name: Cd33Se33
    active_space: [50, 50]
    dt: 1
    path_hdf5: "Cd33Se33.hdf5"
    path_traj_xyz: "Cd33Se33_fivePoints.xyz"
    scratch_path: "/tmp/singlepoints_advanced"
    calculate_guesses: "first"

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 28.0
      periodic: none
      executable: cp2k.popt

      cp2k_settings_main:
        specific:
          template: pbe_main
          cp2k:
            force_eval:
              dft:
                scf:
                  eps_scf: 1e-6

      cp2k_settings_guess:
        specific:
          template:
            pbe_guess
            

In your working directory, create an *input_test_single_points_advanced.yml* file and copy the previous input inside it (remember to respect the indentation). 
Also copy locally the small pre-computed MD trajectory of the Cd33Se33 system, Cd33Se33_fivePoints.xyz_.

In the input file, pay particular attention to the following parameters that have been added/modified with respect to the previous example:

- **dt**: The size of the timestep used in your MD simulations (in fs).
- **path_traj_xyz**: Path to the pre-computed MD trajectory. It should be provided in xyz format.
- **calculate_guesses**: Specify whether to calculate the guess wave function only in the first point of the trajectory ("first") or in all ("all). Here, we keep the default value, first.

In this example, we also show how to further personalize the cp2k_general_settings. In particular, a cp2k subsection is added to overwrite some parameters of the pbe template_ and tighten the scf convergence criterion to 1e-6 (the default value in the pbe_main template is 5e-4). Please note that a specific indentation is used to reproduce the  structure of a typical cp2k input file. By using this approach, you can easily personalize your input requirements by referring to the cp2k manual_.

A more elaborate example would have involved the computation of the eigenvalues and energies of orbitals in the **triplet** state for each point of this **singlet** trajectory. This would have been done by simply adding ``multiplicity: 3`` under the cp2k_general_settings block.

.. _Cd33Se33_fivePoints.xyz: https://github.com/SCM-NV/nano-qmflows/blob/master/test/test_files/Cd33Se33_fivePoints.xyz
.. _template: https://github.com/SCM-NV/nano-qmflows/blob/master/nanoqm/workflows/templates.py
.. _manual: https://manual.cp2k.org/

Setting up the calculation 
---------------------------

Once all settings of your yml input have been customized, you are ready to launch your single point calculation.

- First, activate the conda environment with QMFlows:

  ``conda activate qmflows``
  
- Then, load the module with your version of cp2k, for example:

  ``module load CP2K/7.1.0``
  
- Finally, use the command run_workflow.py to submit your calculation:

  ``run_workflow.py -i input_test_single_points_basic.yml``
  
  or 

  ``run_workflow.py -i input_test_single_points_advanced.yml``
  
  for the advanced example.
Distribute Absorption Spectrum 
==============================

This workflow computes the absorption spectra for a given molecular system and returns a set of files in TXT format. The principle of distribution workflow is dividing the work in multiple, separated, instances (chunks), in order to be able to split time-consuming jobs into smaller, quicker ones.

In this tutorial, we want to compute the excited states at **each point** of a pre-computed MD trajectory for the guanine system. Please note that this trajectory has already been used in the advanced example of the absorption_spectrum tutorial_, where the spectrum analysis was performed on 1 out of 4 points only. Here we take advantage of the distribution workflow to increase by four times the accuracy of sampling with no substantial variation of the computational cost in terms of time by dividing the job in five chunks (each taking charge of 4 points out of 20).

Preparing the input
--------------------

The input is described in YAML format as showed in the following example:

.. code-block:: yaml 

    workflow:
      distribute_absorption_spectrum
    
    project_name: guanine_distribution
    active_space: [20, 20]
    dt: 1
    path_hdf5: "guanine.hdf5"
    path_traj_xyz: "guanine_twentyPoints.xyz"
    scratch_path: "/tmp/distribute_absorption_spectrum"
    calculate_guesses: "first"
    
    xc_dft: pbe
    tddft: stda
    stride: 1
    
    blocks: 5
    workdir: "."
    
    job_scheduler:
      free_format: "
      #! /bin/bash\n
      #SBATCH --job-name=guanine_distribution\n
      #SBATCH -N 1\n
      #SBATCH -t 1:00:00\n
      #SBATCH -p short\n
      source activate qmflows\n
      module load CP2K/7.1.0\n\n"
    
    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 25.0
      periodic: none
      executable: cp2k.popt
    
      cp2k_settings_main:
        specific:
          template: pbe_main
    
      cp2k_settings_guess:
        specific:
          template: pbe_guess


In your working directory, create an *input_test_distribute_absorption_spectrum.yml* file and copy the previous input inside it (remember to respect the indentation). 
Also copy locally the small pre-computed MD trajectory of the guanine system, guanine_twentyPoints.xyz.

In the input file, pay particular attention to the following parameters that have been added/modified with respect to the previous tutorial_ (advanced example):

- **stride**: Controls the accuracy of sampling of geometries contained in the MD trajectory of reference. Here, a value of stride: 1 indicates that the spectrum analysis will be performed on each point in the reference trajectory. Two important things have to be pointed out:

  #. The workflow will perform SCF calculations for each point in the trajectory; only afterwards it will sample the number of structures on which the spectrum analysis will be performed

  #. Down-sampling issues might arise from the number of points that are actually printed during the MD calculations. Some programs, indeed, offer the possibility to print (in the output file) only one point out of ten (or more) calculated. For example, applying a stride: 10 would in practice mean that you are sampling 1 point out of 100 points in the trajectory.

- **blocks**: Indicates into how many blocks has the job to be split. This will generate as many chunks’ folders in your working directory.

- **workdir**: Path to the chunks' folders.

The **job_scheduler** can also be found below these parameters. Customize these settings according to the system and environment you are using to perform the calculations.

.. _tutorial: https://qmflows-namd.readthedocs.io/en/latest/absorption_spectrum.html
.. _tutorial: https://qmflows-namd.readthedocs.io/en/latest/absorption_spectrum.html

Setting up the calculation 
---------------------------

Once all settings in *input_test_distribute_absorption_spectrum.yml* have been customized, you will need to create the different chunks. 
  
- First, activate QMFlows:

  ``conda activate qmflows``  

- Use the command *distribute_jobs.py* to create the different chunks:

  ``distribute_jobs.py -i input_test_distribute_absorption_spectrum.yml``

A number of new folders are created. In each folder you will find a submission file, launch.sh, a sub-trajectory file (containing the points assigned to that chunk), chunk_xyz, and an input.yml file. In the input.yml file, you will find all your settings. Check for any possible manual errors.

- If you are satisfied with the inputs, submit each of your jobs for calculation.

You can keep track of the calculations by going to your scratch path. The location where all points of the chunks are calculated is your assigned scratch path plus project name plus a number.

Results 
-------

Once the calculations are completed, you will find multiple *output_n_stda.txt* files in your scratch directories (with *n* being the index of the geometry at which the spectrum analysis has been performed). The first two lines of the file *output_0_stda.txt*, found in /tmp/distribute_absorption_spectrum/scratch_chunk_0/ are reported below.

::

    # state    energy       f      t_dip_x    t_dip_y    t_dip_y    weight   from   energy  to     energy     delta_E
        1      4.566    0.03832   -0.51792   -0.25870    0.08573    0.50158  20     -5.175  21     -1.261      3.914

For each excited state (line), the first six columns contain, from left to right:

- *# state*: Assigned index, in ascending order of energy. Here, the lowest excitation is reported and corresponds to # state 1.
- *energy*: Transition energy, in eV.
- *f*: Oscillator strength, dimensionless.
- *t_dip_x*, *t_dip_y*, *t_dip_z*: Transition dipole moment components along x, y and z.

The next six columns report some useful information about the dominant single orbital transition for the excited state under examination:

- *weight*: Weight in the overall transition. Always 1.0000 in the Single Orbital approximation.
- *from*: Index of the initial occupied orbital in the active space.
- *energy*: Energy of the initial occupied orbital.
- *to*: Index of the final virtual orbital in the active space.
- *energy*: Energy of the final virtual orbital.
- *delta_E*:Energy of the dominant single orbital transition. Corresponds to the excited state energy in the Single Orbital approximation.

Copy all the output files to your working directory and plot the absorption spectrum (averaged over all sampled structures) using the script convolution.py_:

  ``convolution.py -nm True``
  
To plot the absorption spectrum of a specific sample, for example our point 0, use the -n option.

  ``convolution.py -n 0 -nm True``

.. _convolution.py: https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/qmflows/convolution.py

Reporting a bug or requesting a feature
---------------------------------------
To report an issue or request a new feature you can use the github issues_ tracker.

.. _HDF5: http://www.h5py.org/
.. _issues: https://github.com/SCM-NV/nano-qmflows/issues
.. _QMflows: https://github.com/SCM-NV/qmflows
.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _YAML: https://pyyaml.org/wiki/PyYAML



Welcome to nano-qmflows's documentation!
========================================

Contents:

.. toctree::
   :maxdepth: 2
   :caption: Introduction 

   includereadme
   theory

.. toctree::
   :maxdepth: 2
   :caption: Tutorials 

   intro
   single_points
   coop
   ipr
   derivative_couplings
   absorption_spectrum
   distribute_absorption_spectrum

.. toctree::
   :maxdepth: 2
   :caption: Library Documentation

   documentation



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: ../README.rst
Introduction to the Tutorials
=============================

The *nano-qmflows* packages offers the following set of workflows to compute different properties:
 * single_points
 * coop_calculation
 * ipr_calculation
 * derivative_coupling
 * absorption_spectrum
 * distribute_absorption_spectrum

Known Issues
------------

Distribution of the workflow over multiple nodes
################################################

`CP2K` can uses multiple nodes to perform the computation of the molecular orbitals using the **MPI** protocol. Unfortunately, the `MPI` implementation for the computation of the *derivative coupling matrix* is experimental and unestable. The practical consequences of the aforemention issues, is that **the calculation of the coupling matrices are carried out in only 1 computational node**. It means that if you want ask for more than 1 node to compute the molecular orbitals with `CP2K`, once the workflow starts to compute the *derivative couplings* only 1 node will be used at a time and the rest will remain idle wating computational resources.


Reporting a bug or requesting a feature
---------------------------------------
To report an issue or request a new feature you can use the github issues_ tracker.

.. _HDF5: http://www.h5py.org/
.. _issues: https://github.com/SCM-NV/nano-qmflows/issues
.. _QMflows: https://github.com/SCM-NV/qmflows
.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _YAML: https://pyyaml.org/wiki/PyYAML


Absorption Spectrum 
===================

This other workflow computes the excited states energies, transition dipole moments and oscillator strength using the STDDFT approach.

Preparing the input
--------------------

Basic Example
^^^^^^^^^^^^^

Below is a basic example of input file for the computation of the first 400 (20*20, as setted in the active_space) lowest lying excited states of a guanine molecule at the sTDA level of approximation.

.. code-block:: yaml

    workflow:
      absorption_spectrum

    project_name: guanine
    active_space: [20, 20]
    compute_orbitals: True
    path_hdf5: "guanine.hdf5"
    path_traj_xyz: "guanine.xyz"
    scratch_path: "/tmp/absorption_spectrum_basic"

    xc_dft: pbe
    tddft: stda

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 25.0
      periodic: none
      executable: cp2k.popt

      cp2k_settings_main:
        specific:
          template: pbe_main

      cp2k_settings_guess:
        specific:
          template: pbe_guess

In your working directory, create an *input_test_absorption_spectrum_basic.yml* file and copy the previous input inside it, by paying attention to preserve the correct indentation.
Also copy locally the file containing the coordinates of the relaxed geometry of the guanine in an xyz format, guanine.xyz.

At this point, your *input_test_absorption_spectrum_basic.yml* contains all the settings to perform the excited states calculation and needs to be edited according to your system and preferences. First, let’s recall some parameters that are common to all input files: **workflow**, **project_name**, **active_space**, **path_hdf5**, **path_traj_xyz**, **scratch_path**.

- **workflow**: The workflow you need for your calculations, in this case absorption_spectrum to compute excited states properties.
- **project_name**: Project name for the calculations. You can choose what you wish.
- **active_space**: Range of (occupied, virtual) molecular orbitals to be computed. For example, if 50 occupied and 100 virtual should be considered in your calculations, the active space should be set to [50, 100].
- **compute_orbitals**: Specify if the energy and eigenvalues of the selected orbitals are to be computed. The default is set to True so we will not consider it in the advanced examples.
- **path_hdf5**: Path where the hdf5 should be created / can be found. The hdf5 is the format used to store the molecular orbitals and other information.
- **path_traj_xyz**: Path to the pre-optimized geometry of your system. It should be provided in xyz format.
- **scratch_path**: A scratch path is required to perform the calculations. For large systems, the .hdf5 files can become quite large (hundredths of GBs) and calculations are instead performed in the scratch workspace. The final results will also be stored here.

You can find the complete list of these general options in this dictionary_.

Also pay particular attention to the following parameters, specific to the absorption_spectrum workflow:

- **xc_dft**: Type of exchange-correlation functional used in your DFT calculations.
- **tddft**:  Type of approximation used in the excited states calculations. The Single Orbital (sing_orb), sTDA (stda) and sTDDFT (stddft) approximations are available.

In the cp2k_general_settings, you can customize the settings used to generate the cp2k input of the initial single point calculations (from which Molecular Orbital energies and coefficients are retrieved). For more details about this section please refer to the available tutorial_ on single point calculations. To further personalize the input requirements, also consult the cp2k manual_ and the templates_ available in nano-qmflows.

.. _dictionary: https://github.com/SCM-NV/nano-qmflows/blob/e176ade9783677962d5146d8e6bc5dd6bb4f9102/nanoqm/workflows/schemas.py#L116
.. _schema_cp2k_general_settings: https://github.com/SCM-NV/nano-qmflows/blob/e176ade9783677962d5146d8e6bc5dd6bb4f9102/nanoqm/workflows/schemas.py#L55
.. _templates: https://github.com/SCM-NV/nano-qmflows/blob/master/nanoqm/workflows/templates.py
.. _manual: https://manual.cp2k.org/
.. _tutorial: https://github.com/SCM-NV/nano-qmflows/blob/master/docs/single_points.rst


Advanced Example
^^^^^^^^^^^^^^^^

We are now ready to move to a more advanced example in which we want to compute the excited states of our guanine molecule starting from a pre-computed MD trajectory rather than a single geometry. The input file will look like that:

.. code-block:: yaml

    workflow:
      absorption_spectrum

    project_name: guanine
    active_space: [20, 20]
    dt: 1
    path_hdf5: "guanine.hdf5"
    path_traj_xyz: "guanine_twentyPoints.xyz"
    scratch_path: "/tmp/absorption_spectrum_advanced1"
    calculate_guesses: "first"

    xc_dft: pbe
    tddft: stda 
    stride: 4

    cp2k_general_settings:
      basis:  "DZVP-MOLOPT-SR-GTH"
      potential: "GTH-PBE"
      cell_parameters: 25.0
      periodic: none
      executable: cp2k.popt

      cp2k_settings_main:
        specific:
          template: pbe_main

      cp2k_settings_guess:
        specific:
          template: pbe_guess

In your working directory, create an *input_test_absorption_spectrum_advanced.yml* file and copy the previous input inside it (remember to respect the indentation). 
Also copy locally the small pre-computed MD trajectory of the guanine system, guanine_twentyPoints.xyz.

In the input file, pay particular attention to the following parameters that have been added/modified with respect to the previous example:

- **dt**: The size of the timestep used in your MD simulations (in fs).
- **path_traj_xyz**: Path to the pre-computed MD trajectory. It should be provided in xyz format.
- **calculate_guesses**: Specify whether to calculate the guess wave function only in the first point of the trajectory ("first") or in all ("all). Here, we keep the default value, first.
- **stride**: Controls the accuracy of sampling of geometries contained in the MD trajectory of reference. For example, our value of stride: 4 indicates that the spectrum analysis will be performed on 1 out of 4 points in the reference trajectory. Two important things have to be pointed out:

  #. The workflow will perform SCF calculations for each point in the trajectory (twenty points in our example); only afterwards it will sample the number of structures on which the spectrum analysis will be performed (here six structures corresponding to points 0, 4, 8, 12, 16, 20).

  #. Down-sampling issues might arise from the number of points that are actually printed during the MD calculations. Some programs, indeed, offer the possibility to print (in the output file) only one point out of ten (or more) calculated. In this case, applying a stride: 4 would in practice mean that you are sampling 1 point out of 40 points in the trajectory.

Setting up the calculation 
---------------------------

Once all settings of your yml input have been customized, you are ready to launch your single point calculation.

- First, activate the conda environment with QMFlows:

  ``conda activate qmflows``
  
- Then, load the module with your version of cp2k, for example:

  ``module load CP2K/7.1.0``
  
- Finally, use the command run_workflow.py to submit your calculation:

  ``run_workflow.py -i input_test_absorption_spectrum_basic.yml``
  
for the basic example.

Results 
-------

Once your calculation has finished successfully, you will find one (or more) *output_n_stda.txt* file(s) in your scratch directory (with *n* being the index of the geometry at which the spectrum analysis has been performed). The first two lines of the file *output_0_stda.txt* generated in our basic example are reported below.

::

    # state    energy       f      t_dip_x    t_dip_y    t_dip_y    weight   from   energy  to     energy     delta_E
        1      4.566    0.03832   -0.51792   -0.25870    0.08573    0.50158  20     -5.175  21     -1.261      3.914

For each excited state (line), the first six columns contain, from left to right:

- *# state*: Assigned index, in ascending order of energy. Here, the lowest excitation is reported and corresponds to # state 1.
- *energy*: Transition energy, in eV.
- *f*: Oscillator strength, dimensionless.
- *t_dip_x*, *t_dip_y*, *t_dip_z*: Transition dipole moment components along x, y and z.

The next six columns report some useful information about the dominant single orbital transition for the excited state under examination:

- *weight*: Weight in the overall transition. Always 1.0000 in the Single Orbital approximation.
- *from*: Index of the initial occupied orbital in the active space.
- *energy*: Energy of the initial occupied orbital.
- *to*: Index of the final virtual orbital in the active space.
- *energy*: Energy of the final virtual orbital.
- *delta_E*:Energy of the dominant single orbital transition. Corresponds to the excited state energy in the Single Orbital approximation.

Copy the output file(s) to your working directory and plot the absorption spectrum using the script convolution.py_:

  ``convolution.py -nm True``
  
In case of multiple output files, the returned absorption spectrum is an average over all sampled strucutures, unless you define the index of a specific sample using the -n option.

.. _convolution.py: https://github.com/SCM-NV/nano-qmflows/blob/master/scripts/qmflows/convolution.py

Reporting a bug or requesting a feature
---------------------------------------
To report an issue or request a new feature you can use the github issues_ tracker.

.. _HDF5: http://www.h5py.org/
.. _issues: https://github.com/SCM-NV/nano-qmflows/issues
.. _QMflows: https://github.com/SCM-NV/qmflows
.. _PYXAID: https://www.acsu.buffalo.edu/~alexeyak/pyxaid/overview.html
.. _YAML: https://pyyaml.org/wiki/PyYAML


