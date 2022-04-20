# Project Authors

MANDYOC is being developed and currently maintained by

- [Victor Sacek](https://github.com/victorsacek): sacek@usp.br -
IAG, Universidade de São Paulo, Brazil - (ORCID: 0000-0001-9598-5081)

The following people have made significantly contributions to the project
(in alphabetical order by last name) are considered "The MANDYOC Developers":

- [Jamison Assunção](https://github.com/jamisonassuncao): jamison.assuncao@usp.br - IAG, Universidade de São Paulo, Brazil - (ORCID: 0000-0003-2822-2417)
- [Agustina Pesce](https://github.com/aguspesce): pesce.agustina@gmail.com - Instituto Geofísico Sismológico Ing. Volponi, Universidad Nacional de San Juan, Argentina; CONICET, Argentina - (ORCID: 0000-0002-5538-8845)
- [Rafael Monteiro da Silva](https://github.com/rafaelmds): rafael.m.silva@alumni.usp.br - IAG, Universidade de São Paulo, Brazil - (ORCID: 0000-0001-8645-2443)
# Contributing to MANDYOC

Contributions to MANDYOC are welcome.
If you have an issue, a bug, a code contribution or a documentation
contribution, **thanks for helping to improve MANDYOC!**

**How to contribute to MANDYOC?**

- Submitting bug reports and feature requests
- Writing tutorials or examples
- Fixing typos and improving to the documentation
- Writing code for everyone to use

## Reporting a Bug

Find the Issues and click New Issue.
Remember to _choose bug label_ for this type of issues.

When creating a bug report, please be as specific as possible when describing
how to reproduce an issue, and include both the intended/expected result and
what you are actually getting.
Attach as much as possible of the following information to your issue:

- a minimal parameter file that reproduces the issue,
- the log.txt file that was created during the model run,
- the error message you saw on your screen,
- any information that helps us understand why you think this is a bug, and how to reproduce it.

## Making MANDYOC Better

**First off, thank you for considering contributing to our project!**

If you want to make some kind of contribution to MANDYOC, please note the
following general guidelines.

### General guidelines

We _follow the
[git pull request workflow](https://www.asmeurer.com/git-workflow/) to make
changes to our code base_.
Every change made goes through a pull request, even our own, so that our
continuous integration services have a change to check that the code is up to
standards and passes all our tests.
This way, the master/main branch is always stable.

General guidelines for pull requests (PRs):

- **Open an issue first describing what you want to do**.
  If there is already an issue that matches your PR, leave a comment there
  instead to let us know what you plan to do.

- **Create a fork** of the repository on your personal account.
  **Then, clone it** in your computer:

  ```bash
  git clone your-fork-url
  ```

- **Create a separate branch** (sometimes called a feature branch) on which
  you do your modifications:

  ```bash
  git checkout -b branch-name
  ```

- Once you have created your branch, **make your changes and commit them**.
  Remember to keep your commits atomic, that is, each commit should represent a
  single unit of change.
  Also, remember to write helpful commit messages, so that someone can
  understand what the commit does.

  ```bash
  git add filename

  git commit
  ```

  After that **push up your changes** to your fork:

  ```bash
  git push
  ```

- **Make a pull request**.
  Ensure the PR description clearly describes the problem and solution.
  Include the relevant issue number.

**Remember to test your changes**. We have a regression test to check that the Mandyoc output results are equal to the expected result.
To make this test you need to install `pytest` and run:

```bash
make test
```

_NOTE_: you might need to manually set the environment variables `PETSC_DIR` (and `PETSC_ARCH', if you are not using the prefix installation).

## Do you have questions about MANDYOC code?

Ask any question you have by **opening up an Issue**, and labeling it as
a _question_.
# Contributor Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age,
body size, disability, ethnicity, gender identity and expression, level of
experience, nationality, personal appearance, race, religion, or sexual
identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual attention or
  advances
- Trolling, insulting/derogatory comments, and personal or political attacks
- Public or private harassment
- Publishing others' private information, such as a physical or electronic
  address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in a
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
We see all of these actions as last resorts.
Our goal is to maintain a friendly and inclusive community where everyone feels
welcome to participate and express themselves, but conduct by individuals that
jeopardizes the harmony of the project will need to be addressed.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting [Victor Sacek](sacek@usp.br).

All complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. The project team
is obligated to maintain confidentiality with regard to the reporter of an
incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the
[Contributor Covenant](http://contributor-covenant.org),
version 1.4, available at
[http://contributor-covenant.org/version/1/4](http://contributor-covenant.org/version/1/4/)
---
title: "Mandyoc: A finite element code to simulate thermochemical convection in parallel"
tags:
  - PETSc
  - mantle convection
  - lithosphere geodynamics
  - finite element

authors:
  - name: Victor Sacek
    orcid: 0000-0001-9598-5081
    affiliation: 1
  - name: Jamison Assunção
    orcid: 0000-0003-2822-2417
    affiliation: 1
  - name: Agustina Pesce
    orcid: 0000-0002-5538-8845
    affiliation: "2,3"
  - name: Rafael Monteiro da Silva
    orcid: 0000-0001-8645-2443
    affiliation: 1
affiliations:
  - name: Instituto de Astronomia, Geofísica e Ciências Atmosféricas, Universidade de São Paulo, Brazil
    index: 1
  - name: Instituto Geofísico Sismológico Ing. Volponi, Universidad Nacional de San Juan, Argentina
    index: 2
  - name: CONICET, Argentina
    index: 3
date: 10 June 2021
bibliography: paper.bib
# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

`Mandyoc` is a 2-D finite element code written in C dedicated to simulating thermochemical convection in the interior of terrestrial planets.
Different linear and non-linear rheologies can be adopted, appropriately simulating the strain and stress pattern in the Earth's crust and mantle, both in extensional and collisional tectonics.
Additionally, the code allows variations of boundary condition for the velocity field in space and time, simulating different pulses of tectonism in the same numerical scenario.

# Statement of need

`Mandyoc`, the acronym for `MANtle DYnamics simulatOr Code`, is designed to simulate Stokes flow type thermochemical convection taking different compositional layers into account, and it is also appropriate to simulate Earth's lithospheric dynamics on geological timescales. Similar codes are available for geodynamic problems, like `ASPECT` [@KHB12;@aspect2017;@Bangerth2021;@wolfgang_bangerth_2021_5131909] and `Underworld` [@moresi2007computational;@beucher2019uwgeodynamics]. Therefore, `Mandyoc` is an alternative to preexistent softwares.

One advantage of `Mandyoc` is the possibility to create scenarios with velocity boundary conditions variable in space and time, allowing the user to simulate different tectonic pulses in the same model run.
Additionally, the current version incorporates surface processes, imposing rates of erosion and sedimentation on the top of the free surface [@silva2022influence]. Recently, other thermomechanical codes available for the scientific community incorporated the interaction with surface processes [e.g., @beucher2020morphotectonic; @neuharth2022evolution] taking into account the simulation of fluvial and hillslope processes. In the present version of `Mandyoc`, only predefined erosion/sedimentation rate (variable in space and time) is possible. 

Previous versions of the code were used to study the evolution of continental margins, showing the interaction of the continental lithosphere with the asthenospheric mantle [@sacek2017post;@salazar2021lateral].

# Mathematics

`Mandyoc` solves the equations for conservation of mass, momentum and energy using the Finite Element Method assuming the extended Boussinesq approximation, respectively:

\begin{align*}
u_{i,i} &= 0 \\
\sigma_{ij,j} + g_i \rho &= 0 \\
  \frac{\partial T}{\partial t} + u_i T_{,i} &=
  \kappa T_{,ii} + \frac{H}{c_p \rho} + \frac{u_i g_i \alpha T}{c_p}
\end{align*}

where

\begin{align*}
\sigma_{ij} &= -P \delta_{ij} + \eta \left( u_{i,j} + u_{j,i} \right), \\
\rho &= \rho_0 \left( 1 - \alpha (T - T_0) \right),
\end{align*}

$u_i$ is the component $i$ of the velocity field, $T$ is temperature, $t$ is time, $\kappa$ is the thermal diffusivity, $H$ is the volumetric heat production, $c_p$ is the specific heat capacity, $g$ is gravity, $\rho$ is the effective rock density dependent on temperature and composition, $\rho_0$ is the reference rock density at temperature $T_0$, $\alpha$ is the coefficient of thermal expansion, $P$ is the dynamic pressure, $\eta$ is the effective viscosity, and $\delta_{ij}$ is the Kronecker delta.

The code is fully parallelized using the Portable, Extensible Toolkit for Scientific Computation (PETSc) [@petsc-efficient;@petsc-user-ref;@petsc-web-page].
The present version of the code can simulate thermochemical convection using different rheological formulations: Newtonian flow, non-linear viscous flow or visco-plastic deformation.
For example, the lithosphere can be simulated as a combination of different visco-plastic layers in which the effective viscosity depends on a nonlinear power law viscous rheology and a plastic yield criterion, like the Drucker-Prager criterion.
Additionally, strain softening is implemented to facilitate the localization of strain in the plastic regime during, for example, lithospheric stretching.

The composition and strain history is tracked by particles present in the interior of the finite element.
The exchange of particles among the subdomains of the model is efficiently parallelized in PETSc using DMSwarm [@may2017dmswarm].

The free surface of the Earth can be simulated and is numerically stabilized using the Free Surface Stabilization Algorithm [@kaus2010stabilization].
Surface processes of erosion and sedimentation can also be incorporated in the thermo-mechanical model, allowing the coupling between the Earth's interior dynamics and the processes occurring at its surface.
Complex boundary conditions for the velocity field, variable both in space and time, can be adopted by the user to simulate different episodes of tectonism.
Different benchmarks are available in the repository and can be reproduced by the user, including thermochemical convection [@van1997comparison] and plume-lithosphere interaction [@crameri2012comparison].

As an example application of `Mandyoc`, \autoref{fig:rift} presents snapshots of one numerical scenario of lithospheric stretching imposing a divergent flow direction, resulting in rifting and break-up.
In this example, the upper crust, lower crust, lithospheric mantle and asthenosphere present different rheology and density, resulting in faulting mainly in the upper crust and part of the lithospheric mantle.
Additionally, deformations in the lower crust and at the base of the lithospheric mantle are accommodated by ductile creep flow. This example can be reproduced from the repository.

![`Mandyoc` example of application of the thermo-mechanical model to simulate the stretching of the lithosphere, assuming different rheologies. The scales of gray represent cumulative strain in the different materials. Details can be found in the repository.\label{fig:rift}](JOSS_figure.png)

# Acknowledgements

We acknowledge contributions from Dave May for the help with the implementation of the multigrid algorithm.
This project was sponsored by FAPESP (Processes 2017/24870-5 and 2019/23246-1) and Petrobras (Process 2017/00461-9).

# References
.. image:: https://joss.theoj.org/papers/10.21105/joss.04070/status.svg
   :target: https://doi.org/10.21105/joss.04070

.. raw:: html

    <div class="banner">
        <h1>Mandyoc</h1>
        <p> <strong>MANtle DYnamics simulatOr Code</strong></p>
    </div>


About
-----------
*Mandyoc* is a finite element code written on top of the `PETSc`_ library  to simulate thermo-chemical convection of the Earth's mantle.
Different linear and non-linear rheologies can be adopted, appropriately simulating the strain and stress pattern in the Earth's crust and mantle, both in extensional or collisional tectonics.

Documentation
----------------------------

On the `Documentation page <https://ggciag.github.io/mandyoc/>`__, you will be able to find relevant information on how the code  works, how it was implemented, how to install it, the description of the parameters and input files, some applications and some useful examples.

Need any help?
------------------------

* Most discussion happens `on Github <https://github.com/ggciag/mandyoc>`__.
  Feel free to `open an issue
  <https://github.com/ggciag/mandyoc/issues/new>`__ or comment  on any open issue or pull request.
*  You can also send an email to sacek@usp.br or jamison.assuncao@usp.br

Contributing
---------------------

Code of conduct
+++++++++++++++++

Please note that this project is released with a `Contributor Code of
Conduct <https://github.com/ggciag/mandyoc/blob/main/CODE_OF_CONDUCT.md>`__.
By participating in this project you agree to abide by its terms.

Contributing Guidelines
++++++++++++++++++++++++

**Contributions to Mandyoc are welcome**.
If you have an issue, a bug, a code contribution or a documentation contribution, **thanks for helping to improve Mandyoc!**
Check out our `Contributing guide
<https://github.com/ggciag/mandyoc/blob/main/CONTRIBUTING.md>`__.


License
-----------

This is free software, you can redistribute it and/or modify it under the terms
of the **BSD 3-clause License**. A copy of this license is provided in
`LICENSE <https://github.com/ggciag/mandyoc/blob/main/LICENSE>`__.

.. _PETSc: https://www.mcs.anl.gov/petsc/
.. mandyoc-docs documentation master file, created by
   sphinx-quickstart on Mon Nov  2 14:45:10 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. include:: ../README.rst

.. toctree::
    :hidden:
    :numbered:
    :maxdepth: 2
    :caption: Getting started

    files/theory
    files/implementation
    files/installation
    files/parameter-file
    files/input
    files/benchmarks
    files/references

.. toctree::
    :hidden:
    :numbered:
    :maxdepth: 2
    :caption: Getting help and contributing

    How to contribute <https://github.com/ggciag/mandyoc/blob/main/CONTRIBUTING.md>
    Code of Conduct <https://github.com/ggciag/mandyoc/blob/main/CODE_OF_CONDUCT.md>
    Source code <https://github.com/ggciag/mandyoc>
References
==========

.. bibliography::
.. _inputfiles:

Input files
===========

*Mandyoc* accepts ASCII files to use as initial conditions for the simulation, such as an initial temperature file ``input_temperature_0.txt``, a file containing the initial interfaces geometry ``interfaces.txt`` and a file containing the initial velocity field at the faces of the model ``input_velocity_0.txt``. The next subsections will help you creating each one of these files.


ASCII temperature file
----------------------

..
   For both 2-D and 3-D grids, the initial temperature configuration can be provided as an ASCII file called ``intial-temperature.txt``. Considering :math:`x` to be the longitudinal direction, :math:`y` the vertical direction, and :math:`z` the latitudinal direction, the next subsections will guide you through the understanding of the initial temperature file structure.

For a 2-D grid, the initial temperature configuration can be provided as an ASCII file called ``input_temperature_0.txt``. Considering :math:`x` to be the longitudinal direction and :math:`y` the vertical direction, the next subsections will guide you through the understanding of the initial temperature file structure.

2-D temperature grid
********************

Consider a 2-D grid where the number of nodes in the :math:`x` and :math:`y` directions are :math:`nx\geq 2` and :math:`ny\geq 2`, respectively, with :math:`(nx-1)` elements in the :math:`x` direction and :math:`(ny-1)` elements in the :math:`y` direction. Each node in the grid can be identified with a coordinates pair :math:`(x_n, y_m)`, where :math:`n` and :math:`m` are natural numbers, and :math:`x_n` and :math:`y_m` are the :math:`x` and :math:`y` positions (in meters) of any grid node. :numref:`coordinates2d` shows a scheme for a 2-D grid, where the red dots represent the grid nodes and every node is labeled with a :math:`(x_n, y_m)` pair. The dotted lines represent any number of intermediate nodes.

.. Consider a 2-D grid with :math:`(nx-1)>2` elements in the :math:`x` direction and :math:`(ny-1)>2` elements in the :math:`y` direction. The number of nodes in the :math:`x` and :math:`y` directions are :math:`nx` and :math:`ny`, respectively, and each node can be identified with a pair :math:`(x_n, y_m)`, where :math:`n \leq (nx-1) \in \mathbb{N}` and :math:`m \leq (ny-1) \in \mathbb{N}`, where :math:`x_n` and :math:`y_m` are the :math:`x` and :math:`y` positions (in meters) of any grid node. :numref:`coordinates2d` shows a scheme for a 2-D grid, where the red dots represent the grid nodes, the the dotted lines represent any number of element nodes, and every node is labeled with a :math:`(x_n, y_m)` pair.

.. _coordinates2d:

.. figure:: figs/coordinates-2d.png
   :width: 90%
   :align: center
   :alt: Coordinates

   2-D grid scheme. The red dots represent the grid nodes where the values for the temperature must be defined, the dotted lines represent any number of grid nodes, for simplification, and the labels indicate the :math:`x` and :math:`y` position of each node.

..
   .. note::
      Because of the way *Mandyoc* deals with the coordinates in 2-D, the node at :math:`(x_0,y_0)` is always at the origin :math:`(0,0)`. This is also true for the 3-D grid, where the :math:`(x_0,y_0,z_0)` is at :math:`(0,0,0)`.

.. note::
   Because of the way *Mandyoc* deals with the coordinates in 2-D, the node at :math:`(x_0,y_0)` is always at the origin :math:`(0,0)`.

The example below shows how the ``input_temperature_0.txt`` file must be written for the 2-D grid in :numref:`coordinates2d`: after writing four lines of comments, the file must include the temperature :math:`T(x_n, y_m)` at each node starting at the bottom left of the model and going up in the :math:`y` direction. Intuitively, the temperature is given in horizontal layers for every :math:`x_n`.

.. literalinclude:: src/initial-temperature-2d.txt
   :language: text
   :linenos:

.. tip::
   It is always helpful to use the four comment lines to write down any useful information for later. Because these lines are simply skipped, there is no rule to write them, the '#' symbol is used only by convention.

..
   3-D temperature grid
   ********************

   Constructing the ``input_temperature.txt`` file for a three dimensional grid (:numref:`coordinates3d`) is very similar to the 2-D case, except it is necessary to provide the temperature of each horizontal layer of nodes :math:`xz` for each :math:`y=y_0,y_1,...,y_{ny-1}`. For the 3-D grid in :numref:`coordinates3d`, the number of nodes in the :math:`z` direction is called :math:`nz` and the temperature at each node is :math:`T(x_n, y_m, z_p)`, where :math:`p \leq (nz-1) \in \mathbb{N}`.

   .. _coordinates3d:

   .. figure:: figs/coordinates-3d.png
      :width: 90%
      :align: center
      :alt: Coordinates

      3-D grid scheme. The red dots represent the grid nodes where the values for the temperature must be defined. The dotted lines represent any number of grid nodes, for simplification.

   The example below shows the ``input_temperature.txt`` file for the generic 3-D grid in :numref:`coordinates3d`. Following the four lines of comments, the temperature value at each node must be provided in each line.

   .. literalinclude:: src/initial-temperature-3d.txt
      :language: text
      :linenos:

ASCII interfaces file
---------------------

..
   The ``input_interfaces.txt`` file, for both 2-D and 3-D grids, starts with seven lines of variables that are used by the rheology models. The variables that are in the interfaces are those that might be used in :eq:`power-law` or :eq:`frank-kamenetskii`: :math:`C`, :math:`\rho_r`, :math:`H`, :math:`A`, :math:`n`, :math:`Q`, and :math:`V`. Each one of these variables possesses a number of values (columns) that are assigned to each lithological unit in the model. Every interface is a boundary between two lithological units, therefore the number of columns is 1 plus the number of interfaces set in the ``param.txt`` file (see the :doc:`parameter file section<parameter-file>`), such that if there are :math:`i` interfaces, :math:`i+1` values for each variable **must** be provided in the ``input_interfaces.txt`` file. More will be discussed about the order of these values shortly.

The ``interfaces.txt`` file, for a 2-D grid, starts with seven lines of variables that are used by the rheology models. The variables that are in the interfaces are those that might be used in :eq:`power-law` or :eq:`frank-kamenetskii`: :math:`C`, :math:`\rho_r`, :math:`H`, :math:`A`, :math:`n`, :math:`Q`, and :math:`V`. Each one of these variables possesses a number of values (columns) that are assigned to each lithological unit in the model. Every interface is a boundary between two lithological units, therefore the number of columns is 1 plus the number of interfaces set in the ``param.txt`` file (see the :doc:`parameter file section<parameter-file>`), such that if there are :math:`i` interfaces, :math:`i+1` values for each variable **must** be provided in the ``interfaces.txt`` file. More will be discussed about the order of these values shortly.


2-D Initial interfaces
**********************

Below, the example corresponds to a 2-D grid with two interfaces and, therefore, three lithological units. The first column contains the vertical positions **in meters** of every grid node :math:`y_m` that corresponds to the  **deepest** interface boundary, starting at :math:`x=x_0` on line 8, and ending at :math:`x=x_{nx-1}` on line :math:`nx+7`. The second column contains the vertical position of every :math:`y_m` that corresponds to the second interface boundary. When defining the interfaces, it is rather common for them to "touch". Because of that, all the interfaces must be provided in a "tetris" manner, where interfaces that are collinear in parts fit the interface below.

.. note::
   Because the interfaces are defined linearly between the nodes, it is important to define them properly, so every point inside the grid can be attributed to a lithological unit.

.. literalinclude:: src/initial-interfaces-2d.txt
   :language: text
   :linenos:

The values of the variables in the first seven lines are the values of the lithological units bound by the interfaces. For the 2-D grid with two interfaces, the first values of each variable refers to the lithological unit below the first interface, the second value of each variable refers to the lithological unit between the first and second interfaces, and the third value of each variable refers to the lithological unit above the second interface. If more interfaces are added, there will be more units bounded by an upper and a lower interface.

..
   3-D Initial interfaces
   **********************

   The 3-D grid ``input_interfaces.txt``  file is very similar to the 2-D grid one. The difference lies that after completing the interfaces for :math:`z=z_0`, another set of values is added subsequently to :math:`z=z_1` and so on, until :math:`z=z_{nz-1}`. The ``input_interfaces.txt`` file below shows and example of such case.

   .. literalinclude:: src/initial-interfaces-3d.txt
      :language: text
      :linenos:


ASCII velocities files
----------------------

..
   For both 2-D and 3-D grids, the ``input_velocity.txt`` file possesses four lines of headers followed by the velocity data for each node of the model.

An initial velocity file allows an initial velocity field to be established so the simulated fluid presents an initial momentum. Respecting the conservation of mass equation when designing the initial velocity field is fundamental so the simulation can run properly.

Additionally, the velocity boundary conditions defined in the ``param.txt`` (see the :doc:`parameter file section<parameter-file>`) will be used to define the normal and tangential velocities on the faces of the model during the simulation.

2-D initial velocity
********************

When ``velocity_from_ascii = True`` is set in the ``param.txt`` file, an initial velocity field must be provided in an ASCII file called ``input_velocity_0.txt``.

For a 2-D grid, the ``input_velocity_0.txt`` file possesses four lines of headers followed by the velocity data for each node of the model. The velocity must be defined for both :math:`x` and :math:`y` directions. Taking the example in the file below, the ``input_velocity_0.txt`` file must be written in such a way that the values for the velocity in the :math:`x` direction :math:`v_x` and the velocity in the :math:`y` direction :math:`v_y` for the node :math:`(x_0,y_0)` corresponds to lines 5 and 6, respectively, the velocity values for the :math:`(x_1,y_0)` corresponds to lines 7 and 8, and so on until :math:`(x_{nx-1},y_0)`. Once all the nodes in the :math:`y_0` are written, the values for all :math:`x_n` are inserted for :math:`y_1`, and so on until :math:`y_{ny-1}`.

.. literalinclude:: src/initial-velocity-2d.txt
   :language: text
   :linenos:

2-D multiple initial velocities
*******************************

The user can provide multiple velocity fields to be set in different instants of the simulation. By setting ``multi_velocity = True`` in the ``param.txt`` file, an ASCII file called ``multi_veloc.txt`` must be provided following the structure shown in the example file below. The first line of the file contains the amount of instants where a different velocity field will be used and the other lines contain the time, in millions of years, that each velocity field will be adopted. In the example below, a different velocity field will be used at three instants, the files must follow the structure presented in the last section and should be named in sequence: ``input_velocity_1.txt``, ``input_velocity_2.txt`` and ``input_velocity_3.txt``. The names of the files are labeled after the three different instants set in the first line of the ``multi_veloc.txt`` file.

.. literalinclude:: src/multi_veloc.txt
   :language: text
   :linenos:

.. note::
   The ``input_velocity_0.txt`` is used for the beginning of the simulation. The multiple initial velocities that can be set during the simulation start at ``input_velocity_1.txt`` and end at ``input_velocity_X.txt``, where :math:`X` is the number of different instants that the velocity field will be changed. 

2-D re-scalable boundary condition
**********************************

When ``variable_bcv = True`` is set in the ``param.txt``file, an ASCII file called ``scale_bcv.txt`` is used to define a sequence of re-scaling steps for the velocity boundary conditions. The example file below is an example of how such file can be written. In the example file, the first line corresponds to the number of instants during the simulation that the velocity field will be re-scaled. The first column contains the time, in millions of years, that the velocity field will be changed and the second column contains the scale value. For the example below, in the instant :math:`10.0` Myr the velocity field will be multiplied by :math:`0.5`, effectively halving it; at :math:`25.0` Myr the current velocity field will have its direction changed, hence multiplied by :math:`-1.0`; and finally at :math:`50` My the velocity field will have its direction changed again and doubled by multiplying it by :math:`-2.0`.

.. literalinclude:: src/scale_bcv.txt
   :language: text
   :linenos:


..
   3-D initial velocity
   ********************

   Similarly to the 2-D grid, the ``input_velocity_0.txt`` file for a 3-D grid must contain the velocity values for each :math:`x`, :math:`y` and :math:`z` components. Therefore, the resulting files differs from the 2-D case as shown in the example below, which specifies the velocity for each node in the :math:`xz` plane for each depth :math:`y_m`.

   .. literalinclude:: src/initial-velocity-3d.txt
      :language: text
      :linenos:


How to install
==============

*Mandyoc* installation is very simple and it consists of installing both `PETSc`_
and *Mandyoc*.

.. warning::
	The following installation steps work for both Linux and macOS machines
	**only** and no tests were made to install *Mandyoc* on Windows machines yet.

.. _Dependencies:

Dependencies
------------

To build *Mandyoc*, the following requirements are needed:

* PETSc_ (currently tested on version v3.15.5)
* gcc
* make
* git (recommended, but not strictly needed)

If you do not already have a PETSc installation, you may need a Fortran compiler.

Additionally, the following additional software is needed to run the examples
that come with *Mandyoc*:

* Python 3.5+

With required python packages:

* numpy
* matplotlib
* pandas
* jupyterlab

To run the tests, some additional python packages are required:

* pytest

To build the documentation, further python packages are necessary:

* sphinx
* sphinx_rtd_theme

PETSc Installation
------------------

*Mandyoc* requires the PETSc library to run.
Check out the `PETSc installation`_ for details.

.. note::

	The following steps its a example installation with minimum requirements
	to run *Mandyoc*.

The first step is to **download** PETSc (v3.15.5) release from `PETSc website`_
or **clone** the repository into your machine.
*Mandyoc* might work with latest release of PETSc, but this is not guaranteed
since new verions might introduce breaking changes.

Clone the repository to your desired location::

	git clone -b v3.15.5 https://gitlab.com/petsc/petsc.git

Second, **configure the PETSc build** and set up the installation directory.

.. code-block:: bash

	cd path/to/petsc
	./configure \
	  PETSC_DIR=/path/to/petsc \
	  PETSC_ARCH=arch-label-optimized \
	  --with-debugging=0 \
	  --with-cc=gcc \
	  --with-cxx=g++ \
	  --with-fc=gfortran \
	  --download-fblaslapack \
	  --download-mpich \
	  --download-mumps \
	  --download-scalapack \
	  --download-parmetis \
	  --download-metis

.. note::
	This example installation uses the ``gfortran`` compiler.
	You may use another Fortran compiler changing the ``--with-fc`` flag.
	It is also possible to configure PETSc without Fortran by using ``--with-fc=0``
	and changing ``--download-fblaslapack`` to ``--download-f2cblaslapack``.

.. note::

	By default, *Mandyoc* uses direct solvers (LU and Cholesky) provided by `MUMPS`_.
	This requires additional external packages. Refer to `PETSc documentation`_
	for further information.

.. note::

	If you want to build a development version of *Mandyoc*
	its recommended to build a **debug version** of PETSc
	by setting ``--with-debugging=1``.
	In this case, you may set ``PETSC_ARCH=arch-label-debug``.

.. note::

	If you prefer *openmpi*, you need to swith ``--download-mpich`` to ``--download-openmpi``.

**Check** the installation with:

.. code-block::

	make all check

Or follow the instructions that pop up on the terminal.

For further information about the PETSc library, check out the `PETSc website`_.

Finally, add a symlink of `mpirun` to `~/.local/bin`:

.. code-block::

	ln -s /path/to/pets/arch-label-optimized/bin/mpirun ~/.local/bin/mpirun


*Mandyoc* Installation
----------------------

To install the *Mandyoc* in your machine, you need to **clone or download the latest release** of the code from the `Mandyoc repository`_.
To clone the repository, navigate to the directory you wish to install *Mandyoc* and type:

.. code-block:: bash

   git clone https://github.com/ggciag/mandyoc

Before to install Mandyoc, you mast *set an environment variable* which indicates the path to PETSc installation folder:

.. code-block:: bash

	export PETSC_DIR=/path/to/petsc

*Build Mandyoc* by running:

.. code-block::

	make all

Next, *install Mandyoc* with:

.. code-block::

	make install

By default, it will be installed in ``~/.local/bin``.

.. note::

	Make sure the directory ``~/.local/bin`` exists, otherwise the above command will fail.
	You can change the installation location setting ``INSTALL_PATH`` variable by running:

	.. code-block::

		make INSTALL_PATH=/path/to/install/mandyoc install

.. note::

	To print *Mandyoc* runtime options, run mandyoc with ``-flags`` command line
	argument.

**Check** Mandyoc installation with:

.. code-block::

	make test

.. note::

	You need python and some python packages to run the last commmand succesfully.
	Check out requirements in `Dependencies`_ section.

Docker Container
----------------

We provide a `Docker container`_ image for *Mandyoc*.
Docker is an implementation of container virtualization.
Citing their documentation "it is a lightweight, standalone, executable package of software
that includes everything needed to run an application:
code, runtime, system tools, system libraries and settings".

Visit the `Dockerhub Mandyoc repository`_ to find out more on how to use the container to run *Mandyoc*.

.. note::

	To use the *Mandyoc* docker image, it is required to install the Docker Engine.
	Find out more on `Install Docker Engine`_ page.

Examples
--------

The benchmarks and other experiments are located in the `examples <https://github.com/ggciag/mandyoc/tree/main/examples>`_ folder of the Mandyoc repository.

Inside each example folder, you find a Jupyter notebook with detailed explanation and instructions on how to run the experiment.



.. _PETSc: https://petsc.org/release/
.. _PETSc installation: https://petsc.org/release/install/
.. _PETSc website: https://petsc.org/release/download/
.. _PETSc documentation: https://petsc.org/main/docs/manualpages/Mat/MATSOLVERMUMPS.html
.. _Mandyoc repository: https://github.com/ggciag/mandyoc
.. _MUMPS: http://mumps.enseeiht.fr/
.. _Docker container: https://www.docker.com/resources/what-container
.. _Dockerhub Mandyoc repository: https://hub.docker.com/r/ggciag/mandyoc
.. _Install Docker Engine: https://docs.docker.com/engine/install/
.. _parameterfile:

Parameter file 
==============

.. role:: raw-html(raw)
    :format: html

The parameter file ``param.txt`` contains the information that is necessary for the simulation to run. 

#. Geometry

    * nx
        default: :raw-html:`<br />` 
        type: integer :raw-html:`<br />` 
        unit: :raw-html:`<br />` 
        definition: number of nodes in the horizontal direction
    * nz
        default: :raw-html:`<br />` 
        type: integer :raw-html:`<br />` 
        unit: :raw-html:`<br />` 
        definition: number of nodes in the vertical direction
    * lx                                  
        default: :raw-html:`<br />` 
        type: real :raw-html:`<br />` 
        unit: m :raw-html:`<br />` 
        definition: extent in the horizontal direction
    * lz                                  
        default: :raw-html:`<br />`
        type: real :raw-html:`<br />`
        unit: m :raw-html:`<br />`
        definition: extent in the vertical direction

#. Simulation options

    * solver                              
        default: direct :raw-html:`<br />`
        type: direct/iterative :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the solver to be direct or iterative
    * denok                               
        default: 1.0e-4 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: tolerance criterion for the Uzawa's scheme
    * rtol                                
        default: 1.0e-5 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: the absolute size of the residual norm (relevant only for iterative methods)
    * RK4                                 
        default: Euler  :raw-html:`<br />`
        type: Euler/Runge-Kutta :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: particles advection method
    * Xi_min                              
        default: 1.0e-14        :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: tolerance criterion for the convergence of the non-linear flow
    * random_initial_strain               
        default: 0.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: non-dimensional value for the initial strain perturbation for the entire domain
    * pressure_const                      
        default: -1.0 (i.e. not used) :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa :raw-html:`<br />`
        definition: set constant pressure value for the domain (relevant when 2-D is plain view)
    * initial_dynamic_range               
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: method to smoothen convergence of the velocity field in scenarios with wide viscosity range, see Gerya (2019) :cite:`gerya2019`
    * periodic_boundary                   
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows simulation with periodic boundary in the horizontal direction
    * high_kappa_in_asthenosphere         
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: mimics high heat transport in the asthenosphere increasing its thermal diffusivity coefficient
    * basal_heat                          
        default: -1.0 (i.e. not used)  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: W/m^2 :raw-html:`<br />`
        definition: set basal heat flux value

#. Particles options

    * particles_per_element               
        default: 81           :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in each element
    * particles_per_element_x             
        default: 0 (automatic calculation) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in the horizontal direction
    * particles_per_element_z             
        default: 0 (automatic calculation) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: number of Lagrangian particles in the vertical direction
    * particles_perturb_factor            
        default: 0.5  :raw-html:`<br />`
        type: real number between 0 and 1 :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: indicates the amount of perturbation of the initial location of the particles relative to a regular grid distribution. 

#. Surface processes

    * sp_surface_tracking                 
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows free surface tracking across time and outputs it
    * sea_level                           
        default: 0.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m :raw-html:`<br />`
        definition: sea level used to limit the influence of the surface process
    * sp_surface_processes                
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows erosion and sedimentation simulation
    * sp_dt                               
        default: 0              :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: time step for surface processes simulation
    * a2l                                 
        default: True  :raw-html:`<br />`
        type:  True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows the conversion of air particles to land particles during sedimentation
    * sp_mode                             
        default: 1              :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: specify the surface processes method
    * free_surface_stab                   
        default: True  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if the free surface stabilization algorithm (FSSA) is used, see Kaus et al. (2010) :cite:`kaus2010`
    * theta_FSSA                          
        default: 0.5            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: weight of the influence of the FSSA method (only relevant when <free_surface_stab> is True)
    * sticky_blanket_air
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows the increase of viscosity for the first air layer of particles
    * precipitation_profile_from_ascii    
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if precipitation profile along the horizontal axis is read from an ASCII file
    * climate_change_from_ascii           
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: if True, re-scales through time the precipitation profile using an ASCII file

#. Time constrains
    
    * step_max                            
        default: :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: steps :raw-html:`<br />`
        definition: maximum time-step of the simulation
    * time_max                            
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time of the simulation
    * dt_max                              
        default: :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time between steps of the simulation 
    * step_print                          
        default:              :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: steps :raw-html:`<br />`
        definition: make output files every <step_print>
    * sub_division_time_step              
        default: 1.0            :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: re-scale value for the calculated time-step
    * initial_print_step                  
        default: 0 (i.e. not used) :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: <step_print> used until <initial_print_max_time>
    * initial_print_max_time              
        default: 1.0e6 :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: years :raw-html:`<br />`
        definition: maximum time to make output files every <initial_print_step>

#. Viscosity

    * viscosity_reference                 
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: reference mantle viscosity 
    * viscosity_max                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: maximum viscosity during simulation 
    * viscosity_min                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: Pa.s :raw-html:`<br />`
        definition: minimum viscosity during simulation 
    * viscosity_per_element               
        default: constant  :raw-html:`<br />`
        type: constant/variable :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: sets if viscosity is constant or linearly variable for every element
    * viscosity_mean_method               
        default: harmonic  :raw-html:`<br />`
        type: harmonic/arithmetic :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: defines method do calculate the viscosity for each element
    * viscosity_dependence                
        default: depth  :raw-html:`<br />`
        type: pressure/depth :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: defines if viscosity depends on pressure or depth

#. External ASCII inputs/outputs

    * interfaces_from_ascii               
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if interfaces between lithologies are read from an ASCII file (interfaces.txt)
    * n_interfaces                        
        default:        :raw-html:`<br />`
        type: :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the number of interfaces to be read from the interfaces ASCII file (interfaces.txt) 
    * temperature_from_ascii              
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if initial temperature is read from an ASCII file (input_temperature_0.txt)
    * velocity_from_ascii                 
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if initial velocity field is read from an ASCII file (input_velocity_0.txt)
    * variable_bcv                        
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: allows velocity field re-scaling through time according to an ASCII file (scale_bcv.txt)
    * multi_velocity                      
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if boundary velocities can change with time from ASCII file(s) (multi_veloc.txt and additional input_velocity_[X].txt files)
    * binary_output                       
        default: False  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if output is in binary format
    * print_step_files                    
        default: True  :raw-html:`<br />`
        type: True/False :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if the particles position are printed to an output file

#. Physical parameters

    * temperature_difference              
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: K :raw-html:`<br />`
        definition: temperature difference between the top and bottom of the model (relevant if <temperature_from_ascii> is False) 
    * thermal_expansion_coefficient       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: 1/K :raw-html:`<br />`
        definition: value for the coefficient of thermal expansion
    * thermal_diffusivity_coefficient     
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m^2/s :raw-html:`<br />`
        definition: value for the coefficient of thermal diffusivity 
    * gravity_acceleration                
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: m/s^2 :raw-html:`<br />`
        definition: value for the gravity acceleration 
    * density_mantle                      
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: kg/m^3 :raw-html:`<br />`
        definition: value for the mantle reference density 
    * heat_capacity                       
        default:  :raw-html:`<br />`
        type: real number :raw-html:`<br />`
        unit: J/K :raw-html:`<br />`
        definition: value for the heat capacity 
    * non_linear_method                   
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if non linear method is used for the momentum equation 
    * adiabatic_component                 
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if adiabatic heating/cooling is active 
    * radiogenic_component                
        default:  :raw-html:`<br />`
        type: on/off :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set if radiogenic heating is active 

#. Velocity boundary conditions

    * top_normal_velocity                 
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the top side of the model to be fixed or free 
    * top_tangential_velocity             
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the top side of the model to be fixed or free
    * bot_normal_velocity                 
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the bottom side of the model to be fixed or free
    * bot_tangential_velocity             
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the bot side of the model to be fixed or free
    * left_normal_velocity                
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the left side of the model to be fixed or free
    * left_tangential_velocity            
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the left side of the model to be fixed or free
    * right_normal_velocity               
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the normal velocity on the right side of the model to be fixed or free
    * right_tangential_velocity           
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set the tangential velocity on the right side of the model to be fixed or free

#. Temperature boundary conditions

    * top_temperature                     
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the top side of the model to be fixed or free
    * bot_temperature                     
        default:            :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the bottom side of the model to be fixed or free
    * left_temperature                    
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the left side of the model to be fixed or free
    * right_temperature                   
        default:             :raw-html:`<br />`
        type: fixed/free :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: set temperature on the right side of the model to be fixed or free
    * rheology_model                      
        default:  :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: flag number of a pre-defined rheology model to use during simulation
    * T_initial                           
        default:  :raw-html:`<br />`
        type: integer :raw-html:`<br />`
        unit: :raw-html:`<br />`
        definition: flag number of a pre-defined temperature model to use during simulation (relevant when <temperature_from_ascii> is False)

Below there is an example of a parameter file ``param.txt`` used to make a simulation with *Mandyoc*.

.. literalinclude:: src/param.txt
   :language: text
   :linenos:



.. _basictheory:

Basic theory
============

The *Mandyoc* code simulates thermochemical convection of the Earth's mantle. The following sections explain which equations are being solved by the code and the numerical approach that was used.

Basic equations
---------------

To simulate mantle thermochemical convection, we adopted the formulation for non-Newtonian fluids together with the Boussinesq approximation to solve the following equations of conservation of mass (:eq:`mass-conservation`), momentum (:eq:`momentum-conservation`) and energy (:eq:`energy-conservation`) :cite:`zhong2007`.

.. math:: 
	:label: mass-conservation

	u_{i,i}=0 

.. math:: 
	:label: momentum-conservation

	\sigma_{ij,j}+g_{i}\rho_{0}(1-\alpha(T-T_{0}))=0

.. math::
	:label: energy-conservation

	\frac{\partial T}{\partial t} + u_{i}T_{,i}=\kappa T_{,ii} + \frac{H}{c_p}-\frac{\alpha T g u_{e}}{c_{p}}
  
where :math:`u` is the velocity in the :math:`i` direction, :math:`g` is the gravity acceleration, :math:`\rho_{0}` is the reference rock density at temperature :math:`T_0`, :math:`\alpha` is the coefficient of volumetric expansion, :math:`T` is the temperature, :math:`\kappa` is the thermal diffusivity, :math:`H` is the rate of radiogenic heat per unit of mass, :math:`c_{p}` is the specific heat, :math:`\delta_{ij}` is the Kronecker delta, and :math:`\sigma_{ij}` is the stress tensor:

.. math::
	:label: stress-tensor

	\sigma_{ij}=-P\delta_{ij}+\eta (u_{i,j}+u_{j,i})

where :math:`P` is the dynamic pressure and :math:`\eta` is the effective viscosity of the rock.

.. note::
	For the adopted Einstein notation, repeated indexes in the same term represent summation and indexes after comma are partial derivative to a spatial coordinate.

Numerical approach
------------------

The equations of conservation of mass, momentum and energy are solved using the finite element method :cite:`zhong2007`. *Mandyoc* uses hexahedral elements for three dimensional grids and quadrilateral elements for two dimensional grids :cite:`hughes2012`. The :ref:`massmomentumimplementation` subsection presents the numerical methods used to solve the mass and momentum equations and the :ref:`energyimplementation` subsection shows how an implicit formulation was used to solve the energy equation :cite:`braun2003`.

To simulate any scenario, the user **must** provide the parameter file ``param.txt`` and, if necessary, the ASCII files with the initial temperature field, velocity field and/or the initial interfaces of the model. To see how these files can be created/modified, see the section :ref:`parameterfile` and :ref:`inputfiles`. The flowchart in :numref:`mandyocscheme` summarizes the steps *Mandyoc* takes to solve the conservation equations and perform a simulation.

.. _mandyocscheme:

.. figure:: figs/mandyoc-scheme.png
	:width: 100%
	:align: center
	:alt: Flowchart

	Flowchart showing the steps *Mandyoc* takes to solve the equations of conservation of mass, momentum and energy.
	
:numref:`mandyocscheme` shows that once the code starts running and the input files are read (``param.txt`` and the ASCII input files), *Mandyoc* uses the effective viscosity field :math:`\eta` (:eq:`effective-eta`) to calculate the velocity field :math:`u` and checks if the convergence condition satisfies the tolerance :math:`tol` as shown in :eq:`tol` :cite:`thieulot2014`.

.. math::
	:label: tol

	\chi_{f}=1-\frac{(\langle f^{i}- \langle f^{i} \rangle \rangle)\cdot (\langle f^{i+1}- \langle f^{i+1} \rangle \rangle)}{|\langle f^{i}- \langle f^{i} \rangle \rangle|\cdot |\langle f^{i+1}- \langle f^{i+1} \rangle \rangle|}\leq tol

where :math:`f` is a vector with the components of the velocity :math:`u` at all mesh nodes, :math:`i` is the iteration number, the :math:`\langle f \rangle` represents the mean value of :math:`f`, and :math:`tol` is the tolerance parameter.

While the minimum tolerance is not reached, *Mandyoc* utilizes the Uzawa's method to iteratively calculate new :math:`u` and :math:`P` fields. The updated fields modify the viscosity field :math:`\eta`, which in turn disturbs the velocity field again. These fields are updated until tolerance is reached. By default, the tolerance value for *Mandyoc* is :math:`10^{-6}`.

Additionally, the compositional factor :math:`C` is evaluated for an advection as in the equation below. Its solution is calculated placing randomly a number of particles within each finite element of the mesh, which are displaced based on the adjacent node velocity values :cite:`tackley2003`. The individual value for each particle is obtained by linear interpolation of the node values.

.. math::
	:label: advection

	\frac{\partial C}{\partial t} + u_{i}C_{,i} = 0 

Once the velocity field is solved, *Mandyoc* computes the temperature field as a function of the :math:`u`, :math:`\kappa` and :math:`H`. In the next steps, the surface processes are computed and if the maximum time step or the maximum simulation time was not reached, the code updates the time and goes back to compute a new velocity field.
.. _numericalimplementation:

Numerical implementation
========================

The following sections show how the numerical methods were implemented to solve the equations of conservation of mass, momentum and energy that were presented in the :ref:`basictheory` section.

.. _massmomentumimplementation:

Mass and momentum equations
---------------------------

For a solution domain :math:`\Omega`, finding the flow velocity :math:`u_i=g_i+v_i` and pressure :math:`P` states the Galerkin weak formulation for Stokes' flow, where :math:`g_i` is a specified boundary velocity, :math:`v_i` belongs to a set of functions :math:`\mathcal{V}` in which every function is equal to zero on the boundary where the :math:`ith` components of the velocity is :math:`g_i`, and :math:`P` belongs to a set of functions :math:`\mathcal{P}` such that the equations of conservation of mass and momentum can be written as follows :cite:`zhong2007`:

.. math::
    :label: weak-formulation-1

    \int_{\Omega}{w_{i,j}\sigma_{ij}d\Omega} -
    \int_{\Omega}{qu_{i,i}d\Omega} =
    \int_{\Omega}{w_i f_id\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{w_i h_i d\Gamma}}

where :math:`w_i \in \mathcal{V}` and :math:`q \in \mathcal{P}` are weighting functions and the boundary conditions are defined:

.. math::
    :label: boundary-conditions

    u_i = g_i \text{ on } \Gamma_{g_i}, \sigma_{ij}n_j = h_i \text{ on } \Gamma_{h_i}

where :math:`\Gamma_{h_i}` is the boundary where the :math:`ith` components of the forces are set to be :math:`h_i`, :math:`n_j` is the normal vector at the boundary :math:`\Gamma_{h_i}`, and :math:`f_i=g\rho_0 \alpha \delta_{i3}` :cite:`hughes2000`. The weak formulation can be re-written as below:

.. math::
    :label: weak-formulation-2

    \int_{\Omega}{w_{i,j}c_{ijkl}v_{k,l}d\Omega} -
    \int_{\Omega}{qv_{i,i}d\Omega} -
    \int_{\Omega}{w_{i,i}Pd\Omega} = \\
    \int_{\Omega}{w_{i}f_{i}d\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{w_{i}h_{i}d\Gamma}} - 
    \int_{\Omega}{w_{i,j}c_{ijkl}g_{k,l}d\Omega}
    
where :math:`c_{ijkl}=\eta(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})` is obtained from the stress tensor equation (see :eq:`stress-tensor` in the :doc:`basic theory section<theory>`).

The velocity field, the pressure field and the weighting functions shape functions that interpolate the grid points at every node are given:

.. math::
    :label: v-shape-function

    \mathbf{v} = 
    v_i \mathbf{e}_i = 
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A v_{iA}\mathbf{e}_i}

.. math::
    :label: w-shape-function

    \mathbf{w} =
    w_i \mathbf{e}_i =
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A w_{iA}\mathbf{e}_i}

.. math::
    :label: g-shape-function

    \mathbf{g} =
    \sum_{A\in \Omega^{v}-\Gamma^{v}_{g_i}}{N_A g_{iA}\mathbf{e}_i}

.. math::
    :label: p-shape-functionn

    P = \sum_{B\in \Omega^{p}}{M_B P_B}

.. math::
    :label: q-shape-function

    q = \sum_{B\in \Omega^{p}}{M_B q_B}

where :math:`N_A` is the shape function for the velocity at node A, :math:`M_B` is the shape function for the pressure at node B, :math:`\Omega^{v}` is the velocity nodes set, :math:`\Omega^{p}` is the pressure nodes set and :math:`\Gamma^{g}_{g_i}` is the velocity nodes set along the boundary :math:`\Gamma_{g_i}`. *Mandyoc* defines :math:`\Omega^{p}` at the center of each element, while :math:`\Omega^{v}` is defined at every element vertex. This avoids spurious flow solutions and numerical instabilities :cite:`zhong2007` and it keeps the velocity shape functions one order higher than the pressure shape functions, a common strategy used in finite element modeling of incompressible media :cite:`hughes2000`.

From the shape functions above and the Galerkin weak formulation (:eq:`weak-formulation-2`), the following expression can be obtained:

.. math::
    :label: implication-1

    \sum_{B\in\Omega^{v}-\Gamma^{v}_{g_j}}{\Big( \mathbf{e}^{T}_{i} \int_{\Omega}{B^T_A D B_B d\Omega \mathbf{e}_j v_{jB}} \Big)} - 
    \sum_{B\in \Omega^p}{\Big( \mathbf{e}_i \int_{\Omega}{N_{A,i} M_B d\Omega P_{B}} \Big)} = \\
    \int_{\Omega}{N_A \mathbf{e}_i f_i d\Omega} +
    \sum_{i=1}^{n_{sd}}{\int_{\Gamma_{h_i}}{N_A \mathbf{e}_i h_i d\Gamma}} -
    \sum_{B\in \Gamma^{v}_{g_j}}{\Big( \mathbf{e}^{T}_{i} \int_{\Omega}{B^{T}_{A} D B^{T}_{B} d\Omega \mathbf{e}_j g_{jB}} \Big)}
    
and also:

.. math::
    :label: implication-2

    \sum_{B\in \Omega^{v} - \Gamma^{v}_{g_j}}{\int_{\Omega}{M_A N_{Bj} d\Omega \mathbf{e}_j v_{jB}}} = 0

The matrix representation of these two equations above can be presented:

.. math::
    :label: matrix-representation

    \begin{bmatrix}
        K & G \\
        G^T & 0
    \end{bmatrix}
    \left\{
        \begin{array}{c}
            V \\
            P
        \end{array}
    \right\} = 
    \left\{
        \begin{array}{c}
            F \\
            0
        \end{array}
    \right\} 

where :math:`V` is the vector of velocity values at :math:`\Omega^v`, :math:`P` is the vector of pressure values at :math:`\Omega^p`, :math:`F` is the resulting vector of the rigth-hand side of equations :eq:`implication-1` or :eq:`implication-2`, :math:`K` is the stiffness matrix, :math:`G` is the discrete gradient operator, and :math:`G^T` is the discrete divergence operator :cite:`zhong2007`. :math:`K`, :math:`G` and :math:`G^T` are derived from the first and second terms of :eq:`implication-1` and :eq:`implication-2`.

The matrix operator :math:`B` from :eq:`implication-1` contains the spatial derivatives of the shapes function :math:`N` and, for a 2-D plane, can be written as:

.. math::
    :label: B_A

    B_A = 
    \begin{bmatrix}
        N_{A,1} & 0 \\
        0 & N_{A,2} \\
        N_{A,2} & N_{A,1}
    \end{bmatrix}

Again from :eq:`implication-1` and for 2-D plane strain problems, the effective viscosity matrix :math:`D` can be written:

.. math::
    :label: D

    D = 
    \begin{bmatrix}
        2\eta & 0 & 0 \\
        0 & 2\eta & 0 \\
        0 & 0 & \eta
    \end{bmatrix}

The stiffness matrix :math:`K` and the gradient operator :math:`G` can be written:

.. math::
    :label: K

    K_{lm} = \mathbf{e}^T_{i} \int_{\Omega}{B^T_A D B_B d\Omega \mathbf{e}_j}

.. math::
    :label: G

    G_{lm} = \mathbf{e}_i \int_{\Omega}{N_A M_B d\Omega \mathbf{e}_j}

where the subscripts :math:`A` and :math:`B` are the global velocity node numbers, :math:`i` and :math:`j` are the degree of freedom per grid node, ranging from :math:`1` to :math:`n_{sd}`, :math:`l` and :math:`m` are global equation numbers for the velocity ranging from :math:`1` to :math:`n_v n_{sd}`, where :math:`n_{v}` is the number of velocity nodes in the grid.

.. _energyimplementation:

Energy equation
---------------

It is possible to represent the energy conservation equation (:eq:`energy-conservation`) as a finite element problem for a solution domain :math:`\Omega_V` as follow:

.. math::
    :label: fe-energy-conservation

    \mathbf{M} \mathbf{\dot{a}}_T + (\mathbf{K_a} + \mathbf{K_c})\mathbf{a}_T = \mathbf{F}

where :math:`\mathbf{M}`, :math:`\mathbf{K}_a`, :math:`\mathbf{K}_c` and :math:`\mathbf{F}` are written below:

.. math::
    :label: M

    \mathbf{M} = \int_{\Omega_V}{\mathbf{N}^T_V \rho_0 c_p \mathbf{N}_V d\Omega_V}

.. math::
    :label: Ka

    \mathbf{K}_a = \int_{\Omega_V}{\mathbf{N}^T_V  \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_V d\Omega_v}

.. math::
    :label: Kc

    \mathbf{K}_c = \int_{\Omega_V}{\mathbf{B}^T_V \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_v d\Omega_v}

.. math::
    :label: F

    \mathbf{F} = \int_{\Omega_V}{\mathbf{N}^T_V \Big( \frac{H}{c_p} - \frac{\alpha T g u_3}{c_p} \Big) d\Omega_v}

where :math:`\mathbf{N}_V` is a row vector of shape functions, :math:`\mathbf{a}_T` is a column vector of the unknown temperature parameters, :math:`\mathbf{\dot{a}}_T` is its time derivative, and :math:`\mathbf{B}_V \equiv \nabla \mathbf{N}_V`. 

.. note::
    The superscript :math:`T` represents the transpose of the matrix, while the non-superscript :math:`T` represents the temperature.

:math:`\mathbf{M}` and :math:`\mathbf{K}_c` are symmetric, but :math:`\mathbf{K}_a` is not. This asymmetry decreases the accuracy of the solution when advection is more dominant than conduction (:cite:`zienkiewicz2000`, chapter 2). To increase numerical accuracy and stability, *Mandyoc* uses the streamline upwind Petrov-Galerkin process to modify :math:`\mathbf{K}_a` to :math:`\mathbf{K}_a^*` :cite:`zienkiewicz2000,hughes1979,hughes1982`:

.. math::
    :label: Ka-star

    \mathbf{K}_a^* = \int_{\Omega_V}{\mathbf{N}^{*T}_V \rho_0 c_p \mathbf{v} \cdot \mathbf{B}_V d\Omega_V}

where :math:`N^{*}_{Vi}` is:

.. math::
    :label: N-star

    N^{*}_{Vi} = N_{Vi} + \frac{\alpha_{opt} h^e \mathbf{v} \cdot \nabla N_{Vi}}{2|\mathbf{v}|} 

where :math:`h^e` is the characteristic element size in the advection velocity direction (:math:`\mathbf{v}`), and :math:`\alpha_{opt}` is:

.. math::
    :label: alpha-opt

    \alpha_{opt} = \coth{Pe} - \frac{1}{Pe} \text{, where } Pe = \frac{|\mathbf{v}|h^e}{2 \kappa \rho_0 c_p}

The time discretization is made by an implicit scheme :cite:`braun2003`:

.. math::
    :label: time-discretization

    \frac{\mathbf{a}_T(t+\Delta t) - \mathbf{a}_T(t)}{\Delta t} = 
    \theta \mathbf{\dot{a}}_T(t+\Delta t)+(1-\theta)\mathbf{a}_T(t)

where :math:`\theta=0.5` is a weighting parameter. Multiplying both sides of :eq:`time-discretization` by :math:`\mathbf{M}(t+\Delta t)` and assuming :math:`\mathbf{M}(t+\Delta t) \approx \mathbf{M}(t)`:

.. math::
    :label: time-discretization-2

    \mathbf{M}(t+\Delta t) \frac{\mathbf{a}_T(t+\Delta t) - \mathbf{a}_T(t)}{\Delta t} = 
    \theta [\mathbf{F}(t+\Delta t) - \mathbf{K}_T(t+\Delta t) \mathbf{a}_T(t+\Delta t)] + \\
    (1-\theta)[\mathbf{F}(t)-\mathbf{K}_T(t)\mathbf{a}_T(t)]

where :math:`\mathbf{K}_T=\mathbf{K}^*_a+\mathbf{K}_c`. 

Rearranging :eq:`time-discretization-2` allows to rewrite it in a numerical form to solve the energy equation.

.. math::
    :label: numerical-form-energy-equation

    [\mathbf{M}(t+\Delta t) + \Delta t \theta \mathbf{K}_T(t+\Delta t)]\mathbf{a}_T(t)(t+\Delta t) = \\
    [\mathbf{M}(t+\Delta t)- \Delta t(1-\theta)\mathbf{K}_T(t)]\mathbf{a}_T(t) + \\
    \Delta t[\theta \mathbf{F}(t+\Delta t) + (1-\theta)\mathbf{F}(t)]
    
Free surface
------------

*Mandyoc* uses the *Free Surface Stabilization Algorithm* :cite:`kaus2010` to modify the Stokes equation and avoid numerical instabilities that can occur on the surface of the model. The up-and-down oscillations around the steady state ("sloshing instability" or "drunken sailer effect") would require small time steps, which would increase running time by unfeasible amounts.

The modification is done on the stiffness matrix :math:`K_e`, such that :math:`\tilde{K}_e=K_e+L_e`, where the correction :math:`L_e` is evaluated at the boundary :math:`\Gamma_e` of each finite element. The correction is given by the :eq:`Le` below:

.. math::
    :label: Le

    L_e = \int_{\Gamma_e}{\mathbf{N} \Theta \Delta\rho \Delta t \mathbf{g} \mathbf{n} d\Gamma}

where :math:`\mathbf{N}` is the element shape function, :math:`0 \leq \Theta \leq 1` is a wight factor of the correction term, :math:`\Delta \rho` is the density contrast between the two mediums, :math:`\Delta t` is the numerical integration time step, :math:`\mathbf{g}` is the gravity acceleration vector, and :math:`\mathbf{n}` is the normal vector to the element.

.. _rheologysection:

Rheology
--------

Considering a visco-plastic model, the effective viscosity :math:`\eta` follows the formulation described by Moresi and Solomatov (1998) :cite:`moresi1998`, which combines plastic deformation and viscous deformation:

.. math::
    :label: effective-eta
    
    \eta 
    = \min{(\eta_{plas},\eta_{visc})}
    = \min{\bigg(\frac{\tau_{yield}}{2\dot{\varepsilon}_{II}},\eta_{visc}\bigg)}

where :math:`\tau_{yield}` is the rupture tension and :math:`\dot{\varepsilon}_{II}=(\dot{\varepsilon}_{ij}'\dot{\varepsilon}_{ij}'/2)^{1/2}` is the second invariant of the deviatoric strain rate tensor, and :math:`\eta_{plas}` and :math:`\eta_{visc}` are the plastic and ductile viscosities.

Plastic deformation
*******************

The plastic deformation can be calculated using the Byerlee Law :cite:`byerlee1968` to compute :math:`\tau_{yield}` and :math:`\eta_{plas}`. :eq:`byerlee-law` shows the relationship implemented in *Mandyoc*.

.. math::
    :label: byerlee-law

    \tau_{yield} = c_{0}+\mu \rho g z

where :math:`c_{0}` in the internal cohesion, :math:`\mu` is the friction coefficient, :math:`\rho` is the density and :math:`z` is the depth.

Alternatively, the user can choose to use the the Druker-Prager criterion :cite:`drucker-prager1952`, which is presented in :eq:`drucker-prager`.

.. math::
    :label: drucker-prager

    \tau_{yield} = c_0 \cos{\varphi} + P \sin{\varphi}

where :math:`\varphi` is the internal angle of friction.

Viscous deformation
*******************

*Mandyoc* contains several rheology models that the user can choose for viscous deformation. Among them, two will be discussed here. 

The ductile rheology can be simulated using the Frank-Kamenetskii approximation, following the formulation described by Solomatov and Moresi (2000) :cite:`solomatov2000`, where the viscosity is a function of the temperature :math:`T` as in the equation below. Such formulation was also used by Sacek (2017) :cite:`sacek2017` in an earlier *Mandyoc* version.

.. math::
    :label: frank-kamenetskii

    \eta_{visc}(T) = C \eta_r b^* \exp{(-\gamma T)}

where :math:`\eta_r` is the reference viscosity, :math:`C`is a compositional factor to scale the effective viscosity, and :math:`b^*` and :math:`\gamma = E_a / RT^2_b` are constants, which in turn, :math:`E_a` is the activation energy, :math:`R` is the gas constant, and :math:`T_b` is the basal temperature.

Additionally, the rheology can also be considered to follow a power law, as a function of the temperature :math:`T`, compositional factor :math:`C`, pressure :math:`P` and strain rate :math:`\varepsilon` as follows:

.. math::
    :label: power-law

    \eta_{visc} = C A^{\frac{-1}{n}} \dot{\varepsilon}^{\frac{1-n}{n}} \exp{\frac{Q+V P}{nRT}}

where :math:`A` is a pre-exponential scale factor, :math:`n` is the power law exponent, :math:`\dot{\varepsilon}` is the square root of the second invariant of the strain rate tensor, :math:`Q` is the activation energy, and :math:`V` is the activation volume. The values of :math:`A`, :math:`n`, :math:`Q`, and :math:`V` are measured under laboratory conditions :cite:`karato1993,gleason1995`.


Non-linear iterations
*********************

When the non-linear option is chosen by the user, the effective viscosity is dependent on the velocity field, which is iteratively updated following the algorithm described by Thieulot (2014) :cite:`thieulot2014`. In this algorithm, the velocity and effective viscosity field are iteratively updated until the following convergence criterion is satisfied:

.. math::
    :label: non-linear

    \chi_f = 1 - \frac{\langle (f^i-\langle f^i\rangle) \cdot (f^{i+1}-\langle f^{i+1}\rangle) \rangle}{|f^i-\langle f^i\rangle| \ |f^{i+1}-\langle f^{i+1}\rangle|} \le tol 

where :math:`f` represents an array with all the nodal values of the velocity components, :math:`tol` is a tolerance factor, and :math:`\langle f\rangle` is the mean value of :math:`f`. The superscript :math:`i` and :math:`i+1` indicate two consecutive iterations in the same time step.  In each iteration, the momentum and mass equations are calculated with an updated effective viscosity field.





.. _benchmarks:

Benchmarks
==========

In order to test the accuracy of the *Mandyoc* code, its results can be compared to benchmark studies. The following subsections will present the procedures and results for some well established modeling problems. Such cases should provide information about the applicability and performance of the code as well as some of its limitations.

The scripts to build and run these numerical experiments are located in the `examples <https://github.com/ggciag/mandyoc/tree/main/examples>`_ folder of the Mandyoc repository.
Inside each example folder, you find a Jupyter notebook with detailed explanation and instructions on how to run the experiment.

van Keken et al. (1997) :cite:`vankeken1997`
--------------------------------------------

.. note::

  Mandyoc version used for this benchmark: v0.1.4.

The set of simulations proposed by van Keken et al. (1997) :cite:`vankeken1997` compares several methods of studying two dimensional thermochemical convection, where the Boussinesq approximation and infinite Prandtl number are used.

For this benchmark, the first case is presented with the three distinct variations, as proposed by the article. The simulation consists of two layers, where a buoyant thin layer is under a denser thicker package. The problem can be interpreted as a salt layer under a sediment package, and the interface between the layers is defined by the :eq:`interfacerayleigh` below.

.. math::
    :label: interfacerayleigh

    y=-0.8 \lambda_y + 0.02 \cos{\frac{\pi x }{\lambda_x}}

where :math:`\lambda_x` and :math:`\lambda_y` are the horizontal and vertical lengths of the simulated 2-D box, respectively.

The simulations are carried out in a Cartesian box where the fluid is isothermal and Rayleigh-Taylor instability is expected for the proposed setup. The table below lists the parameters used to run this scenario together with the distinction between the three proposed simulation (*case 1a*, *1b* and *1c*).

.. list-table:: Parameters used for the Rayleigh-Taylor instability simulation.
    :header-rows: 1
    :widths: 30 20 20
    :align: center

    * - Parameter
      - .. centered:: Symbol
      - Value
    * - Horizontal length
      - .. centered:: :math:`\lambda_x`
      - 1.0000
    * - Vertical length
      - .. centered:: :math:`\lambda_y`
      - 0.9142
    * - Thermal diffusion coefficient
      - .. centered:: :math:`\kappa`
      - :math:`1.0\times 10^{-6}`
    * - Gravity acceleration
      - .. centered:: :math:`g`
      - :math:`10`
    * - Reference viscosity
      - .. centered:: :math:`\eta_r`
      - :math:`1.0\times 10^{21}`
    * - Buoyant layer viscosity
      - .. centered:: :math:`\eta_0`
      - | :math:`1.00\times\eta_r` (case 1a)
        | :math:`0.10\times\eta_r` (case 1b)
        | :math:`0.01\times\eta_r` (case 1c)

Results for *case 1a*
*********************

For the *case 1a* where :math:`\eta_0/\eta_r=1.00`, :numref:`vankekenCase1aEvolution` below compares the evolution of the isoviscous Rayleigh-Taylor instability between the van Keken et al. (1997) :cite:`vankeken1997` and the *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1aEvolution:

.. figure:: figs/vankeken-snaps-1a.png
  :align: center
  :width: 80%
  :alt: Results

  Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=1.00`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

.. note::
  Because of the different methods used by van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*, the *Mandyoc* results for the evolution of the isoviscous Rayleigh-Taylor instability presents its data colored instead of contoured.

:numref:`vankekenCase1aGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1aGraph:

.. figure:: figs/vrms-1a.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=1.00`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Results for *case 1b*
*********************

For the *case 1b* where :math:`\eta_0/\eta_r=0.10`, :numref:`vankekenCase1bEvolution` compares the evolution of the isoviscous Rayleigh-Taylor instability between van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1bEvolution:

.. figure:: figs/vankeken-snaps-1b.png
  :align: center
  :width: 80%
  :alt: Results

  Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.10`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

:numref:`vankekenCase1bGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1bGraph:

.. figure:: figs/vrms-1b.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.10`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Results for *case 1c*
*********************

For the *case 1c* where :math:`\eta_0/\eta_r=0.01`, :numref:`vankekenCase1cEvolution` compares the evolution of the isoviscous Rayleigh-Taylor instability between van Keken et al. (1997) :cite:`vankeken1997` and *Mandyoc*. The time steps shown for the *Mandyoc* code are the closest the simulation could provide, considering the chosen simulation parameters.

.. _vankekenCase1cEvolution:

.. figure:: figs/vankeken-snaps-1c.png
  :align: center
  :width: 80%
  :alt: Results

  Evolution of the isoviscous Rayleigh-Taylor instability for :math:`\eta_0/\eta_r=0.01`. The best result presented by van Keken et al. (1997) :cite:`vankeken1997` are on the left and the *Mandyoc* results are on the right.

:numref:`vankekenCase1cGraph` below compares the change of the :math:`v_{rms}` with time, showing the results from van Keken et al. (1997) :cite:`vankeken1997` in gray and *Mandyoc* in black.

.. _vankekenCase1cGraph:

.. figure:: figs/vrms-1c.png
  :align: center
  :width: 100%
  :alt: Results

  Evolution of the :math:`v_{rms}` for :math:`\eta_0/\eta_r=0.01`. The van Keken et al. (1997) :cite:`vankeken1997` result is shown in black and the *Mandyoc* code result is shown in gray.

Basic scaling
*************

Below, :numref:`vankekenCase1aScalingGraph` shows the results of a basic scaling test for Mandyoc (v.0.1.4) running on Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz.

.. _vankekenCase1aScalingGraph:

.. figure:: figs/scaling_vanKeken1997_case1a.png
  :align: center
  :width: 85%
  :alt: Scaling results for vanKeken1997_case1a

  Results of scaling test for vanKeken1997_case1a.

Crameri et al. (2012) :cite:`crameri2012`
-----------------------------------------

.. note::

  Mandyoc version used for this benchmark: v0.1.4.

The *Case 2* experiment presented by Crameri et al. (2012) :cite:`crameri2012` evaluates the *sticky air* method to obtain a numerical surface topography in geodynamic modelling.

The experiment analyses the change in topography due to the rising of a mantle plume.
The model setup (:numref:`crameri_setup`) consists of a :math:`2800 \, \mathrm{km}` by :math:`850 \, \mathrm{km}` box with a :math:`150 \, \mathrm{km}` sticky air layer on the top of the model.
The mantle thickness is :math:`600 \, \mathrm{km}` with a :math:`100 \, \mathrm{km}` thick lithosphere.
The lithosphere density is :math:`3300 \, \mathrm{kg/m}^3` with viscosity :math:`10^{23} \, \mathrm{Pa\,s}`,
the mantle density is :math:`3300 \, \mathrm{kg/m}^3` with viscosity :math:`10^{21} \, \mathrm{Pa\,s}`
and the mantle plume density is :math:`3200 \, \mathrm{kg/m}^3` with viscosity :math:`10^{20} \, \mathrm{Pa\,s}`.
Initially, the center of the plume is horizontally centered and :math:`300 \, \mathrm{km}` above the base of the model.
At the top, the sticky air layer has density :math:`0 \, \mathrm{kg/m}^3` with viscosity :math:`10^{19} \, \mathrm{Pa\,s}`.
A free slip boundary condition is applied to the upper boundary of the sticky air layer and the vertical sides of the model and the base is kept fixed.
There is no temperature difference, and the geodynamic evolution is guided solely by compositional density differences.

.. _crameri_setup:

.. figure:: figs/crameri-et-al-2012-case-2-setup.png
	:width: 90%
	:align: center
	:alt: Crameri case 2 model setup

	*Case 2* model setup to evaluate the sticky air method. Extracted from Crameri et al. (2012) :cite:`crameri2012`.

From the results of this experiment reproduced in MANDYOC we obtain the maximum topography with time, similar to Fig. 6a of Crameri et al. (2012) :cite:`crameri2012`, presented in :numref:`maximum_topography`.
The models used for comparison are: UNDERWORLD :cite:`moresi2003`, STAGYY :cite:`tackley1993` and I2VIS :cite:`gerya2003`.
The data used for comparison was extract from Crameri et al. (2012), which used UNDERWORLD v1.5.0.
No version information for STAGYY and I2VIS codes was reported.

.. _maximum_topography:

.. figure:: figs/crameri-et-al-2012-case-2-comparison.png
   :width: 100%
   :align: center
   :alt: Comparison of MANDYOC results

   Comparison of the maximum topography with time for the *Case 2* (:numref:`crameri_setup`) model setup from Crameri et al. (2012) :cite:`crameri2012`.

Basic scaling
*************

Below, :numref:`crameriScalingGraph` shows the results of a basic scaling test for Mandyoc (v.0.1.4) running on Intel(R) Core(TM) i7-10700 CPU @ 2.90GHz.

.. _crameriScalingGraph:

.. figure:: figs/scaling_Crameri2012_case2.png
  :align: center
  :width: 85%
  :alt: Scaling results for Crameri2012_case2

  Results of scaling test for Crameri2012_case2.
