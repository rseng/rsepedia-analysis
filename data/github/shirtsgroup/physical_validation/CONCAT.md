[![GitHub Actions Test Status](https://github.com/shirtsgroup/physical_validation/actions/workflows/continous_integration.yaml/badge.svg)](https://github.com/shirtsgroup/physical_validation/actions/workflows/continous_integration.yaml)
[![GitHub Actions Lint Status](https://github.com/shirtsgroup/physical_validation/actions/workflows/lint.yaml/badge.svg)](https://github.com/shirtsgroup/physical_validation/actions/workflows/lint.yaml)
[![Documentation Status](https://readthedocs.org/projects/physical-validation/badge/?version=latest)](https://physical-validation.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/shirtsgroup/physical_validation/branch/master/graph/badge.svg)](https://codecov.io/gh/shirtsgroup/physical_validation)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/shirtsgroup/physical_validation.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/shirtsgroup/physical_validation/context:python)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/shirtsgroup/physical_validation.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/shirtsgroup/physical_validation/alerts/)  
[![DOI](https://zenodo.org/badge/90674371.svg)](https://zenodo.org/badge/latestdoi/90674371) [![DOI](https://joss.theoj.org/papers/10.21105/joss.03981/status.svg)](https://doi.org/10.21105/joss.03981)  
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

physical_validation: A Python package to assess the physical validity of molecular simulation results
=====================================================================================================

`physical_validation` is a package testing results obtained by molecular simulations for their physical validity.

Please cite 
> Merz PT, Shirts MR (2018) Testing for physical validity in molecular simulations. PLoS ONE 13(9): e0202764. https://doi.org/10.1371/journal.pone.0202764  
> Merz et al., (2022). physical_validation: A Python package to assess the physical validity of molecular simulation results. Journal of Open Source Software, 7(69), 3981, https://doi.org/10.21105/joss.03981

`physical_validation` incorporates the functionality of [checkensemble](https://github.com/shirtsgroup/checkensemble).

This software is developed in the Shirts group at University of Colorado in Boulder.

Documentation
-------------
Please check
[https://physical-validation.readthedocs.io](http://physical-validation.readthedocs.io)
for the full reference.
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into this openff-evaluator. 

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of openff-evaluator on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the openff-evaluator devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)
## Description
Provide a brief description of the PR's purpose here.

## Todos
Notable points that this PR has either accomplished or will accomplish.
  - [ ] TODO 1

## Questions
- [ ] Question1

## Status
- [ ] Ready to go---
title: 'physical_validation: A Python package to assess the physical validity of molecular simulation results'
tags:
  - Python
  - molecular simulation
  - molecular dynamics
  - molecular mechanics
  - statistical mechanics
  - physical validation
authors:
  - name: Pascal T. Merz
    orcid: 0000-0002-7045-8725
    affiliation: 1
  - name: Wei-Tse Hsu
    affiliation: 1
  - name: Matt W. Thompson
    affiliation: 1
  - name: Simon Boothroyd
    affiliation: 2
  - name: Chris C. Walker
    affiliation: 1
  - name: Michael R. Shirts
    orcid: 0000-0003-3249-1097
    affiliation: 1
affiliations:
 - name: Department of Chemical and Biological Engineering, University of Colorado Boulder, Boulder, CO 80309, United States of America
   index: 1
 - name: Boothroyd Scientific Consulting Ltd., 71-75 Shelton Street, London, United Kingdom
   index: 2
date: 7 Jun 2021
bibliography: paper.bib
---

# Summary

Molecular simulations such as molecular dynamics (MD) and Monte Carlo (MC)
simulations are powerful tools allowing the prediction of experimental
observables in the study of systems such as proteins, membranes, and polymeric
materials.
The quality of predictions based on molecular simulations depend on the
validity of the underlying physical assumptions.
`physical_validation` allows users of molecular simulation programs to perform
simple yet powerful tests of physical validity on their systems and setups.
It can also be used by molecular simulation package developers to run
representative test systems during development, increasing code correctness.
The theoretical foundation of the physical validation tests were established
by @Merz:2018, in which the `physical_validation` package was first
mentioned.

# Statement of need

For most of the history of molecular simulation-based research in chemistry, biophysics, 
physics and engineering, most users of molecular simulation packages were experts that 
contributed to the code bases themselves or were very familiar with the methodology used.
Increased popularity of molecular simulation methods has led to a significantly
increased user base and to an explosion of available methods.
The simulation packages are faster and more powerful than ever, and even more than before require
expertise to avoid using combinations of methods and parameters that could violate
physical assumptions or affect reproducibility.
Unphysical artifacts have frequently been reported to significantly affect physical observables
such as the folding of proteins or DNA, the properties of lipid bilayers, the
dynamics of peptides and polymers, or properties of simple liquids (see @Merz:2018
for further references).

# Functionality

`physical_validation` tackles the problem of robustness in molecular
simulations at two levels.
The first level is the end-user level.
`physical_validation` allows users to test their simulation results for a number
of deviations from physical assumptions such as the distribution of the kinetic
energy, the equipartition of kinetic energy throughout the system, the sampling
of the correct ensemble in the configurational quantities, and the precision of
the integrator.
The combination of these tests allow to cover a wide range of potential
unphysical simulation conditions [@Merz:2018].
This increases the confidence that users can have in their simulation results
independently of and in addition to any code correctness tests provided by the developers of their
simulation package.
These validation tools explain their assumptions and conclusions using
comprehensive output and figures, making their use suitable also for users
new to the field of molecular simulations.
Since `physical_validation` also returns its conclusions in machine-readable
form, it can be included in pipelines allowing results to be tested for
physical validity without user interaction.
The second level of usage is by code developers. Unphysical behavior might not only result
from poor or incompatible parameters and models, it might also stem from
coding errors in the simulation programs.
`physical_validation` can therefore be used to regularly run representative simulations
as end-to-end tests in a continuous integration setup, ensuring that code
changes do not introduce bugs that lead to unphysical results.
GROMACS [@Abraham:2015], one of the leading MD packages, has been using `physical_validation`
since 2019 to test every release version with a comprehensive set of
simulations covering all major code paths.

# Relation to prior work

@Shirts:2013 and @Merz:2018 laid the theoretical foundation for
the `physical_validation` package. @Shirts:2013 introduced the
ensemble validation tests, and implemented them in a simple python
script which was made available accompanied by some examples on
github.com/shirtsgroup/checkensemble. This script was used as a base
for the ensemble validation tests in `physical_validation`.
@Merz:2018 built upon the previous work by showing that combining
the ensemble tests with kinetic energy distribution and equipartition
checks as well as integrator convergence tests could detect many types
of unphysical simulation conditions. @Merz:2018 first mentioned
`physical_validation` and its use in the validation of GROMACS
releases.

In the three years since the publication, the software has matured
into a stable release. The ensemble tests now also support $\mu VT$
ensembles, covering the full set of ensembles described by
@Shirts:2013. The user interface, the screen output and the plotting
functionality were polished based on user feedback. The API was
improved and is now considered stable, and the package can be
installed using `conda`, both of which were much requested features
from users looking to use the package in pipelines automating
simulation protocols. While the version published in 2018 had no test
coverage, the stable release is extensively covered by both unit and
regression tests, reaching a test coverage of above 90%. Finally, the
documentation was significantly improved based on user feedback.

# Acknowledgements

* Research reported in this publication was supported in part by the
  National Institute of General Medical Science of the National
  Institutes of Health under award number R01GM115790 (funding PTM)
  and R01RGM132386 (funding MTT), and also in part by the National
  Science Foundation under grant OAC-1835720 (funding WTH) and the
  U.S. Department of Energy, Office of Science, Office of Basic Energy
  Sciences, Materials Sciences and Engineering (MSE) Division, under
  Award Number DE-SC0018651 (funding CCW).
* The Molecular Sciences Software Institute (MolSSI) for a MolSSI
  Software Fellowship to Pascal Merz
* Can Pervane for helpful discussions in the early stages of the
  project
* Nate Abraham for careful reading of the documentation
* Lenny Fobe for help in the setup of the CI

# References
.. _doc_simulation_data:

Creation of :class:`.SimulationData` objects
============================================

The data of simulations to be validated need to be represented by objects
of the  :class:`.SimulationData` type.
The  :class:`.SimulationData` objects are consisting of information about the
simulation and the system. This information is collected in objects of different
classes, namely

* :obj:`.SimulationData.units` of type :class:`.UnitData`:
  Information on the units used by the simulation program.
* :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`:
  Information describing the sampled ensemble.
* :obj:`.SimulationData.system` of type :class:`.SystemData`:
  Information on the system (numbers of atoms, molecules, constraints, etc.).
* :obj:`.SimulationData.observables` of type :class:`.ObservableData`:
  Trajectories of observables along the simulation.
* :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`:
  Position / velocity / force trajectories along the simulation.
* :obj:`.SimulationData.dt` of type :code:`float`:
  The time step at which the simulation was performed.

The :class:`.SimulationData` objects can either be constructed
directly from arrays and numbers, or (partially) automatically via parsers.

Create SimulationData objects from python data
----------------------------------------------

Example usage, system of 900 water molecules in GROMACS units simulated in NVT:
::

   import numpy as np
   import physical_validation

   simulation_data = physical_validation.data.SimulationData()

   num_molecules = 900
   simulation_data.system = physical_validation.data.SystemData(
       # Each water molecule has three atoms
       natoms=num_molecules * 3,
       # Each molecule has three constraints
       nconstraints=num_molecules * 3,
       # In this simulation, translational center of mass motion was removed
       ndof_reduction_tra=3,
       # Rotational center of mass motion was not removed
       ndof_reduction_rot=0,
       # Repeat weight of one oxygen and two hydrogen atoms 900 times
       mass=np.tile([15.9994, 1.008, 1.008], num_molecules),
       # Denotes the first atom of each molecules: [0, 3, 6, ...]
       molecule_idx=np.linspace(0, num_molecules * 3, num_molecules, endpoint=False, dtype=int),
       # Each molecule has three constraints
       nconstraints_per_molecule=3 * np.ones(num_molecules),
   )

   # Set GROMACS units
   simulation_data.units = physical_validation.data.UnitData.units("GROMACS")

   # Simulation was performed under NVT conditions
   simulation_data.ensemble = physical_validation.data.EnsembleData(
       ensemble='NVT',
       natoms=num_molecules * 3,
       volume=3.01125 ** 3,
       temperature=298.15,
   )

   # This snippet is assuming that `kin_ene`, `pot_ene` and `tot_ene` are lists
   # or numpy arrays filled with the time series of kinetic, potential and total energy
   # of a simulation run. These might be obtained, e.g., from the python
   # API of a simulation code, or from other python-based analysis tools.
   simulation_data.observables = physical_validation.data.ObservableData(
       kinetic_energy=kin_ene,
       potential_energy=pot_ene,
       total_energy=tot_ene,
   )

   # We are further assuming that `positions` and `velocities` are arrays
   # of shape (number of frames) x (number of atoms) x 3, where the last
   # number stands for the 3 spatial dimensions. Again, these arrays would
   # most likely have been obtained from a python interface of the simulation
   # package or from other python-based analysis tools
   simulation_data.trajectory = physical_validation.data.TrajectoryData(
       position=positions,
       velocity=velocities,
   )

Package-specific instructions
-----------------------------

GROMACS
~~~~~~~
GROMACS does not offer a well-established Python interface to read out
energies or trajectories. :code:`physical_validation` therefore offers a parser,
which will return a fully populated :class:`.SimulationData` object by
reading in GROMACS input and output files.

The :class:`.GromacsParser` takes the GROMACS input files :code:`mdp` (run options)
and :code:`top` (topology file) to read the details about the system, the ensemble
and the time step. The observable trajectory is extracted from an :code:`edr`
(binary energy trajectory), while the position and velocity trajectory can
be read either from a :code:`trr` (binary trajectory) or a :code:`gro` (ASCII trajectory)
file. The constructor optionally takes the path to a gromacs binary as well
as the path to the topology library as inputs. The first is necessary to
extract information from binary files (using :code:`gmx energy` and :code:`gmx dump`),
while the second becomes necessary if the :code:`top` file contains :code:`#include` statements
which usually rely on GROMACS environment variables. The parser is able to
find GROMACS installations which are in the path (e.g. after sourcing the
:code:`GMXRC` file) and the corresponding topology library automatically.

Example usage:
::

   import physical_validation

   parser = physical_validation.data.GromacsParser()

   res = parser.get_simulation_data(
        mdp='mdout.mdp',
        top='system.top',
        gro='system.gro',
        edr='system.edr'
   )

.. note:: Always double-check the results received from the automatic parser.
   Since this is not an official GROMACS tool, it is very likely that some
   special cases or changes in recent versions might not be interpreted
   correctly.

LAMMPS
~~~~~~
To analyze simulations performed with LAMMPS, we strongly suggest using its
Python interface `Pizza.py <https://pizza.sandia.gov/index.html>`_ to create
a SimulationData object as explained in `Create SimulationData objects from python data`_.
Note that :class:`.UnitData` offers access to a UnitData
object representing the LAMMPS :code:`real` units by using
:code:`.data.UnitData.units("LAMMPS real")`.

As an alternative, :code:`physical_validation` ships with a LAMMPS parser, which tries
to read part of the system information, the observable and position / velocity
trajectories from LAMMPS output files.

Example usage:
::

   import physical_validation

   parser = physical_validation.data.LammpsParser()

   res = parser.get_simulation_data(
       # The LAMMPS parser cannot infer the ensemble from the LAMMPS files,
       # so we pass an EnsembleData object with the information matching the simulation
       ensemble=physical_validation.data.EnsembleData(
           ensemble="NVT",
           natoms=900,
           volume=20**3,
           temperature=300
       ),
       in_file=dir_1 + '/water.in',
       log_file=dir_1 + '/log.lammps',
       data_file=dir_1 + '/water.lmp',
       dump_file=dir_1 + '/dump.atom'
   )

.. warning:: The LAMMPS parser is in an early development stage. It 
   is part of the :code:`physical_validation`
   package in the hope that it is helpful to someone, but it is very
   likely to go wrong in a number of cases. Please check any object data
   create by the LAMMPS parser carefully.

Flatfile parser
---------------

For MD packages not supported by the package-specific parsers, it is possible
to create the :class:`.SimulationData` objects via the
:class:`.FlatfileParser`. This parser fills the
:obj:`.SimulationData.trajectory` object via 3-dimensional ASCII files
containing the position and velocity trajectories, and the
:obj:`.SimulationData.observables` via 1-dimensional ASCII files containing
the trajectories for the observables of interest. As the details on the
units, the simulated system and the sampled ensemble can not easily be read
from such files, this information has to be provided by the user by passing
objects of the respective data structures. See
:attr:`.FlatfileParser.get_simulation_data` for more details on the
:class:`.SimulationData` creation via the flat file parser, and
:ref:`simulationdata_details` for details on which test requires which
information.

Example usage, system of 900 water molecules in GROMACS units simulated in
NVT (note that this example leaves some fields in :class:`.SystemData`
empty, as well as the trajectory of some observables and the position and
velocities):
::

   import physical_validation as pv

   parser = pv.data.FlatfileParser()

   system = pv.data.SystemData(
       natoms=900*3,
       nconstraints=900*3,
       ndof_reduction_tra=3,
       ndof_reduction_rot=0
   )

   # We need to specify the units in which the simulation was performed,
   # specifically the value of k_B in the used energy units, the conversion
   # factor of the simulation units to the physical validation units
   # (*_conversion keywords), and a string representation of the simulation
   # units (*_str keywords, used for output only).
   # See documentation below about UnitData object for more details.
   units = pv.data.UnitData(
       kb=8.314462435405199e-3,
       energy_str='kJ/mol',
       energy_conversion=1.0,
       length_str='nm',
       length_conversion=1.0,
       volume_str='nm^3',
       volume_conversion=1.0,
       temperature_str='K',
       temperature_conversion=1.0,
       pressure_str='bar',
       pressure_conversion=1.0,
       time_str='ps',
       time_conversion=1.0
   )

   ensemble = pv.data.EnsembleData(
       ensemble='NVT',
       natoms=900*3,
       volume=3.01125**3,
       temperature=298.15
   )

   res = parser.get_simulation_data(
       units=units, ensemble=ensemble, system=system,
       kinetic_ene_file='kinetic.dat',
       potential_ene_file='potential.dat',
       total_ene_file='total.dat'
   )

Additional examples
-------------------

Use :code:`MDAnalysis` to create mass vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using :code:`MDAnalysis`, creating a mass vector which can be fed to
:attr:`.SystemData.mass` is straightforward. See the following snippet
for an example using a GROMACS topology:
::

   import MDAnalysis as mda
   import numpy as np

   u = mda.Universe('system.gro')
   mass=np.array([u.atoms[i].mass for i in range(len(u.atoms))])

Use :code:`MDAnalysis` to define molecule groups for equipartition testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`physical_validation.kinetic_energy.equipartition` allows to specify
molecule groups which can be tested for equipartition. The segments used in
:code:`MDAnalysis` can easily be used to define molecule groups as input to
the equipartition check:
::

   import MDAnalysis as mda
   import numpy as np

   u = mda.Universe('system.tpr', 'system.gro')
   molec_groups = []
   for i in range(len(u.segments)):
       seg = u.segments[i]
       molec_groups.append(np.array([seg.atoms[j].index for j in range(len(seg.atoms))]))

Use :code:`MDAnalysis` to read position and velocity trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:code:`MDAnalysis` also makes it easy to create :class:`.TrajectoryData` objects
which require position and velocity trajectories as inputs. Given a
:code:`Universe` object which contains a trajectory, we can simply use a list
comprehension to create a full trajectory in memory:
::

   import MDAnalysis as mda
   import numpy as np
   import physical_validation

   u = mda.Universe('system.tpr', 'system.trr')
   trajectory = physical_validation.data.TrajectoryData(
       position=[frame.positions for frame in u.trajectory],
       velocity=[frame.velocities for frame in u.trajectory])

We can also use the atom selector to only feed part of the trajectory
to the :code:`physical_validation` tests.
This is useful if we want to analyze the equipartition of parts of the
system only (e.g. the solute) which can massively speed up the
validation check. Note that we have to adapt the :class:`.SystemData` object
accordingly to inform :code:`physical_validation` that we are only analyzing
part of the system.
::

   import MDAnalysis as mda
   import numpy as np
   import physical_validation

   u = mda.Universe('system.tpr', 'system.trr')
   protein = u.select_atoms('protein')
   trajectory = physical_validation.data.TrajectoryData(
       position=[protein.positions for _ in u.trajectory],
       velocity=[protein.velocities for _ in u.trajectory])

.. note:: :code:`MDAnalysis` uses Å (ångström) as a length unit. Don't forget to
          choose the :class:`.UnitData` accordingly!


.. _simulationdata_details:

Data contained in :class:`.SimulationData` objects
==================================================

Units: :obj:`.SimulationData.units` of type :class:`.UnitData`
--------------------------------------------------------------
Attributes:

* :attr:`.UnitData.kb`, :code:`float`
* :attr:`.UnitData.energy_conversion`, :code:`float`
* :attr:`.UnitData.length_conversion`, :code:`float`
* :attr:`.UnitData.volume_conversion`, :code:`float`
* :attr:`.UnitData.temperature_conversion`, :code:`float`
* :attr:`.UnitData.pressure_conversion`, :code:`float`
* :attr:`.UnitData.time_conversion`, :code:`float`
* :attr:`.UnitData.energy_str`, :code:`str`
* :attr:`.UnitData.length_str`, :code:`str`
* :attr:`.UnitData.volume_str`, :code:`str`
* :attr:`.UnitData.temperature_str`, :code:`str`
* :attr:`.UnitData.pressure_str`, :code:`str`
* :attr:`.UnitData.time_str`, :code:`str`

The information about units consists of different parts:

* The value of :math:`k_B` in the used energy units,
* the conversion factor to :code:`physical_validation` units (kJ/mol, nm, nm^3, K, bar, ps, the same as the GROMACS default units), and
* the name of the units (:code:`energy_str`, :code:`length_str`, :code:`volume_str`, :code:`temperature_str`, :code:`pressure_str`, :code:`time_str`).

The names are only used for output (console printing and plotting), and are optional.
The conversion factors and kB are, on the other hand, used in computations and need
to be given. To avoid silent errors, these keywords to not have defaults and must be specified. 

Needed by

  * :func:`physical_validation.ensemble.check`
  * :func:`physical_validation.ensemble.estimate_interval`
  * :func:`physical_validation.kinetic_energy.distribution`, only

    - :attr:`.UnitData.kb`

Ensemble: :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`
------------------------------------------------------------------------
Attributes:

* :attr:`.EnsembleData.ensemble`, :code:`str`
* :attr:`.EnsembleData.natoms`, :code:`int`
* :attr:`.EnsembleData.mu`, :code:`float`
* :attr:`.EnsembleData.volume`, :code:`float`
* :attr:`.EnsembleData.pressure`, :code:`float`
* :attr:`.EnsembleData.energy`, :code:`float`
* :attr:`.EnsembleData.temperature`, :code:`float`

The ensemble is a string indicating the thermodynamical ensemble a simulation was
performed in, and is any of :code:'NVE', :code:'NVT', :code:'NPT', :code:'muVT'.

Depending on the ensemble, :class:`.EnsembleData` then holds additional information defining
the ensemble, such as the number of particles N, the chemical potential mu, the
volume V, the pressure P, the constant energy E or the temperature T. While any
of these additional information are technically optional, most of them are needed by certain
tests, such that not fully defining the ensemble results in warnings. The notable
exception to this rule is the constant energy E for NVE, which is not needed
by any test and can hence be omitted without raising a warning.

Needed by
  * :func:`physical_validation.kinetic_energy.distribution`
  * :func:`physical_validation.ensemble.check`

System: :obj:`.SimulationData.system` of type :class:`.SystemData`
------------------------------------------------------------------
Attributes:

    * :attr:`.SystemData.natoms`, the total number of atoms in the system;
      e.g. for a system containing 100 water molecules: :code:`system_data.natoms = 300`
    * :attr:`.SystemData.nconstraints`, the total number of constraints in the
      system, not including the global translational and rotational constraints
      (see next two attributes); e.g. for a system containing 100 *rigid* water molecules:
      :code:`system_data.nconstraints = 300`
    * :attr:`.SystemData.ndof_reduction_tra`, global reduction of translational
      degrees of freedom (e.g. due to constraining the center of mass of the system)
    * :attr:`.SystemData.ndof_reduction_rot`, global reduction of rotational
      degrees of freedom (e.g. due to constraining the center of mass of the system)
    * :attr:`.SystemData.mass`, a list of the mass of every atom in the system;
      e.g. for a single water molecule: :code:`system_data.mass = [15.9994, 1.008, 1.008]`
    * :attr:`.SystemData.molecule_idx`, a list of the index of the first atom of every
      molecule (this assumes that the atoms are sorted by molecule); e.g. for a system
      containing 3 water molecules: :code:`system_data.molecule_idx = [0, 3, 6]`
    * :attr:`.SystemData.nconstraints_per_molecule`, a list of the number of
      constraints in every molecule; e.g. for a system containing 3 *rigid* water
      molecules: :code:`system_data.nconstraints_per_molecule = [3, 3, 3]`
    * :attr:`.SystemData.bonds`, a list containing all bonds in the system;
      e.g. for a system containing 3 water molecules:
      :code:`system_data.bonds = [[0, 1], [0, 2], [3, 4], [3, 5], [6, 7], [6, 8]]`
    * :attr:`.SystemData.constrained_bonds`, a list containing only the constrained
      bonds in the system, must be a subset of :attr:`.SystemData.bonds` (and equal, if
      all bonds are constrained).

.. todo:: Currently, there is some redundancy in the attributes listed above. The
   :attr:`.SystemData.bonds` and :attr:`.SystemData.constrained_bonds` are
   reserved for future use - included already in the information about the system,
   but not yet used by any tests included in the currently published package. In a
   future version, the :class:`.SystemData` should be streamlined to make the object
   initialization easier.

Needed by

  * :func:`physical_validation.kinetic_energy.distribution`, partially:

    - :attr:`.SystemData.natoms`,
    - :attr:`.SystemData.nconstraints`,
    - :attr:`.SystemData.ndof_reduction_tra`,
    - :attr:`.SystemData.ndof_reduction_rot`

  * :func:`physical_validation.kinetic_energy.equipartition`, all attributes except
    :attr:`.SystemData.bonds` and :attr:`.SystemData.constrained_bonds`.

Observables: :obj:`.SimulationData.observables` of type :class:`.ObservableData`
--------------------------------------------------------------------------------
Attributes:

  * :attr:`.ObservableData.kinetic_energy`, the kinetic energy trajectory (nframes x 1),
    also accessible via :code:`observable_data['kinetic_energy']`
  * :attr:`.ObservableData.potential_energy`, the potential energy trajectory (nframes x 1),
    also accessible via :code:`observable_data['potential_energy']`
  * :attr:`.ObservableData.total_energy`, the total energy trajectory (nframes x 1),
    also accessible via :code:`observable_data['total_energy']`
  * :attr:`.ObservableData.volume`, the volume trajectory (nframes x 1),
    also accessible via :code:`observable_data['volume']`
  * :attr:`.ObservableData.pressure` the pressure trajectory (nframes x 1),
    also accessible via :code:`observable_data['pressure']`
  * :attr:`.ObservableData.temperature` the temperature trajectory (nframes x 1),
    also accessible via :code:`observable_data['temperature']`
  * :attr:`.ObservableData.constant_of_motion` the constant of motion trajectory (nframes x 1),
    also accessible via :code:`observable_data['constant_of_motion']`
  * :attr:`.ObservableData.number_of_species` the trajectory of the number of molecules of a species,
    used for muVT, (nframes x num_species),
    also accessible via :code:`observable_data['number_of_species']`

Needed by

  * :func:`physical_validation.kinetic_energy.distribution`

    - :attr:`.ObservableData.kinetic_energy`

  * :func:`physical_validation.ensemble.check`

    - :attr:`.ObservableData.total_energy`, or
    - :attr:`.ObservableData.potential_energy`,
    - :attr:`.ObservableData.volume` (for NPT),
    - :attr:`.ObservableData.number_of_species` (for muVT)

  * :func:`physical_validation.integrator.convergence`

    - :attr:`.ObservableData.constant_of_motion`

Atom trajectories: :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`
-------------------------------------------------------------------------------------
Attributes:

  * :attr:`.TrajectoryData.position`, the position trajectory (nframes x natoms x 3),
    also accessible via :code:`trajectory_data['position']`
  * :attr:`.TrajectoryData.velocity`, the velocity trajectory (nframes x natoms x 3),
    also accessible via :code:`trajectory_data['velocity']`

Needed by

  * :func:`physical_validation.kinetic_energy.equipartition`


Time step: :obj:`.SimulationData.dt` of type :code:`float`
----------------------------------------------------------
The timestep used during the simulation run, a single :code:`float` value.

Needed by

  * :func:`physical_validation.integrator.convergence`
physical_validation\.util subpackage
====================================

.. warning:: This subpackage is intended for internal use. Documentation
    and input validation is reduced.

physical_validation\.util\.kinetic_energy module
------------------------------------------------
.. automodule:: physical_validation.util.kinetic_energy
    :members:
    :undoc-members:

physical_validation\.util\.ensemble module
------------------------------------------
.. automodule:: physical_validation.util.ensemble
    :members:
    :undoc-members:

physical_validation\.util\.integrator module
--------------------------------------------
.. automodule:: physical_validation.util.integrator
    :members:
    :undoc-members:

physical_validation\.util\.trajectory module
--------------------------------------------
.. automodule:: physical_validation.util.trajectory
    :members:
    :undoc-members:

physical_validation\.util\.plot module
--------------------------------------
.. automodule:: physical_validation.util.plot
    :members:
    :undoc-members:

physical_validation\.util\.error module
---------------------------------------
.. automodule:: physical_validation.util.error
    :members:
    :undoc-members:

physical_validation\.util\.gromacs_interface module
---------------------------------------------------
.. automodule:: physical_validation.util.gromacs_interface
    :members:
    :undoc-members:
Introduction
============

Advances in recent years have made molecular dynamics (MD) simulations a
powerful tool in molecular-level research, allowing the prediction of
experimental observables in the study of systems such as proteins, drug
targets or membranes. The quality of any prediction based on MD results
will, however, strongly depend on the validity of underlying physical
assumptions.

This package is intended to help detect (sometimes hard-to-spot)
unphysical behavior of simulations, which may have statistically important
influence on their results. It is part of a two-fold approach to
increase the robustness of molecular simulations.

First, it empowers users of MD programs to test the physical validity on
their respective systems and setups. The tests range from simple
post-processing analysis to more involved tests requiring additional
simulations. These tests can significantly increase the
reliability of MD simulations by catching a number of common simulation
errors violating physical assumptions, such as non-conservative
integrators, deviations from the specified Boltzmann ensemble, or lack of ergodicity
between degrees of freedom. To make usage as easy as possible,
parsers for the outputs of several popular MD programs are provided

Second, it can be integrated in MD code testing environments. While
unphysical behavior can be due to poor or incompatible choices of
parameters by the user, it can also originate in coding errors
within the program. Physical validation tests can be integrated in the
code-checking mechanism of MD software packages to facilitate the
detection of such bugs. The :code:`physical_validation` package is currently
used in the automated code-testing facility of the GROMACS software
package, ensuring that every major releases passes a number of physical
sanity checks performed on selected representative systems before
shipping.

.. note:: The physical validation tests have been described in [Merz2018]_.

.. note:: We are always looking to enlarge our set of tests. If you are a
   MD user or developer and have suggestions for physical validity tests
   missing in this package, we would love to hear from you! Please
   consider getting in touch with us via our `github repository`_.


Installation
============

pip
---
The most recent release of `physical_validation` can be installed from `PyPI`_
via :code:`pip`
::

   pip install physical_validation

conda
-----
The most recent release of `physical_validation` can also be installed using
:code:`conda`
::

   conda install -c conda-forge physical_validation

Development version
-------------------

The latest version is available on our `github repository`_. You can install
it via :code:`pip`
::

   pip install git+https://github.com/shirtsgroup/physical_validation.git


Simulation data
===============

The data of simulations to be validated are represented by objects
of the :class:`.SimulationData` type.

The :class:`.SimulationData` objects contain information about the simulation
and the system. This information is collected in objects of different
classes, namely:

* :obj:`.SimulationData.units` of type :class:`.UnitData`:
  Information on the units used by the simulation program.
* :obj:`.SimulationData.ensemble` of type :class:`.EnsembleData`:
  Information on the sampled ensemble. This includes the temperature, pressure, and chemical potential,
  with specific requirements depending on the ensemble specified.
* :obj:`.SimulationData.system` of type :class:`.SystemData`:
  Information on the system (number of atoms, molecules, constraints, etc.).
* :obj:`.SimulationData.observables` of type :class:`.ObservableData`:
  Trajectories of observables along the simulation, such as energy or volume. 
* :obj:`.SimulationData.trajectory` of type :class:`.TrajectoryData`:
  Position / velocity / force trajectories along the simulation.
* :obj:`.SimulationData.dt` of type :code:`float`:
  The time step at which the simulation was performed.

The different physical validation tests do not all require all data to be
able to run. Each :code:`physical_validation` function checks whether the required
information was provided, and raises an error if the information is
insufficient. :ref:`simulationdata_details` lists by which tests the single
members of :class:`.SimulationData` are required.

The :class:`.SimulationData` objects can either be constructed directly
from arrays and numbers, or (partially) automatically via parsers.
The preferred way to populate :class:`.SimulationData` objects is by
assigning its sub-objects explicitly with data obtained from the simulation
package. Many simulation packages have a well-defined Python interface which
allows to read observable, position and velocity trajectories into Python data
structures. The remaining information, such as details on the simulated
ensemble or the molecular system, is usually rather easy to fill in by hand.
The examples in this documentation follow this model.

Please see :ref:`doc_simulation_data` for more details on the
:class:`.SimulationData` type and on how to create objects from results
obtained from different simulation packages.


Kinetic energy validation
=========================
Kinetic energy validation includes testing the likelihood of a trajectory
to originate from the theoretically expected gamma distribution and
validating the temperature equipartition between groups of degrees
of freedom. For details on the employed algorithms, please check the
respective function documentations.

For both the full distribution test and the equipartition test, a strict
and a non-strict version are available. They are triggered using the
:code:`strict=[True|False]` keyword. The strict version does a full distribution
similarity analysis using the Kolmogorov-Smirnov (K-S) test. The K-S test
returns a p-value indicating the likelihood that the sample originates from
the expected distribution. Its sensitivity
increases with increasing sample size, and can flag even the smallest deviations
from the expected distribution at large sample sizes. When developing or
implementing new temperature control algorithms in a controlled testing
environment which keeps errors from other sources negligible, such a high
sensibility is desirable. In other
applications, however, a deviation insignificant in comparison with
other sources of inaccuracies might be enough to flag long simulation
trajectories of large systems as not having a gamma distribution. For
example, deviations from the desired kinetic energy distribution that
are smaller in magnitude than other well-controlled approximations, such as
the interaction cutoff or the treatment of bond constraints, might be enough
to flag large samples as not being properly distributed.

As an alternative to the strict test, the :code:`physical_validation` suite offers
the non-strict version. In this case, the mean and the standard deviation of
the sample are calculated and compared to the expected values. To make the
test easily interpretable, two distinct temperatures :math:`T_\mu` and
:math:`T_\sigma` are estimated from the kinetic energy distribution. They represent the
temperature at which the sample mean and standard would be physically expected.
An error estimate computed via bootstrapping of the provided kinetic energy samples is given for each of the
temperatures, giving information on the statistical significance of the results.

For more details about the difference between the strict test and non-strict test, please
see :func:`physical_validation.kinetic_energy.distribution`.

Full system distribution validation
-----------------------------------
Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.kinetic_energy.distribution`

Example
~~~~~~~
`Kinetic energy distribution example`_

.. _`Kinetic energy distribution example`: examples/kinetic_energy_distribution.ipynb

Equipartition validation
------------------------
Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.kinetic_energy.equipartition`

Example
~~~~~~~
`Kinetic energy equipartition example`_

.. _`Kinetic energy equipartition example`: examples/kinetic_energy_equipartition.ipynb


Ensemble validation
===================
As the distribution of configurational quantities like the potential
energy :math:`U`, the volume :math:`V` or (for the grand and semigrand canonical ensembles) 
the number of each species :math:`N_i` are in general not known analytically, testing the likelihood
of a trajectory sampling a given ensemble is less straightforward than
for the kinetic energy. However, generally, the *ratio* of the probability
distribution between samplings of the same system generated at different state
points (e.g. simulations run at at different temperatures or different pressures) is exactly known for each ensemble
[Merz2018]_, [Shirts2013]_.
Providing two simulations at different state points therefore allows a
validation of the sampled ensemble.

Note that the ensemble validation function is automatically inferring the
correct test based on the simulation input data (such as temperature and pressure) that are given as input.

Choice of the state points
--------------------------
As the above ensemble tests require two simulations at distinct
state points, the choice of interval between the two points is an
important question. Choosing two state points too far apart will result
in poor or zero overlap between the distributions, leading to very noisy
results (due to sample errors in the tails) or a breakdown of the method,
respectively. Choosing two state points very close to each others, on the
other hand, makes it difficult to distinguish the slope from statistical
error in the samples.

A rule of thumb states [Shirts2013]_ that the maximal efficiency of the
method is reached when the distance between the peaks of the distributions
are roughly equal to the sum of their standard deviations. For most systems
with the exception of extremely small or very cold systems, it is reasonable
to assume that the difference in standard deviations between the state points
will be negligable. This leads to two ways of calculating the intervals:

*Using calculated standard deviations*: Given a simulation at one state point,
the standard deviation of the distributions can be calculated numerically. The
suggested intervals are then given by

* :math:`\Delta T = 2 k_B T^2 / \sigma_E`, where :math:`\sigma_E` is the standard
  deviation of the energy distribution used in the test (potential energy, enthalpy,
  or total energy).
* :math:`\Delta P = 2 k_B T / \sigma_V`, where :math:`\sigma_V` is the standard
  deviation of the volume distribution.

*Using physical observables*: The standard deviations can also be estimated using
physical observables such as the heat capacity and the compressibility. The
suggested intervals are then given by:

* :math:`\Delta T = T (2 k_B / C_V)^{1/2}` (NVT), or
  :math:`\Delta T = T (2 k_B / C_P)^{1/2}` (NPT), where :math:`C_V` and :math:`C_P`
  denote the isochoric and the isobaric heat capacities, respectively.
* :math:`\Delta P = (2 k_B T / V \kappa_T)`, where :math:`\kappa_T` denotes the
  isothermal compressibility.

When setting :code:`verbosity >= 1` in :func:`physical_validation.ensemble.check`, the
routine is printing an estimate for the optimal spacing based on the distributions
provided. Additionally, :func:`physical_validation.ensemble.estimate_interval`
calculates the estimate given a single simulation result. This can be used to determine
at which state point a simulation should be repeated in order to efficiently check
its sampled ensemble.

Function reference
~~~~~~~~~~~~~~~~~~
:func:`physical_validation.ensemble.check`

:func:`physical_validation.ensemble.estimate_interval`

Example
~~~~~~~
`Ensemble validation example`_

.. _`Ensemble validation example`: examples/ensemble_check.ipynb


Integrator Validation
=====================
A symplectic integrator can be shown to conserve a constant of motion
(such as the energy in a microcanonical simulation) up to a fluctuation
that is quadratic in time step chosen. Comparing two or more
constant-of-motion trajectories realized using different time steps (but
otherwise unchanged simulation parameters) allows a check of the
symplecticity of the integration. Note that lack of symplecticity does not
necessarily imply an error in the integration algorithm, it can also hint
at physical violations in other parts of the model, such as non-continuous
potential functions, imprecise handling of constraints, etc.

Functions
---------
:func:`physical_validation.integrator.convergence`

Example
-------
`Integrator convergence example`_

.. _`Integrator convergence example`: examples/integrator_validation.ipynb


.. _`PyPI`: https://pypi.org/project/physical_validation

.. _`github repository`: https://github.com/shirtsgroup/physical_validation

.. [Merz2018] Merz PT, Shirts MR (2018)
   "Testing for physical validity in molecular simulations",
   PLOS ONE 13(9): e0202764.
   https://doi.org/10.1371/journal.pone.0202764

.. [Shirts2013] Shirts, M.R.
   "Simple Quantitative Tests to Validate Sampling from Thermodynamic Ensembles",
   J. Chem. Theory Comput., 2013, 9 (2), pp 909–926,
   http://dx.doi.org/10.1021/ct300688p
*****************************
Physical validation reference
*****************************

:code:`physical_validation` is a package aimed at testing results obtained
by molecular dynamics simulations for their physical validity.

.. note:: The physical validation methodology has been described in
   Merz PT, Shirts MR (2018), Testing for physical validity in molecular simulations.
   PLoS ONE 13(9): e0202764. https://doi.org/10.1371/journal.pone.0202764

.. note:: We are always looking to enlarge our set of tests. If you are a
   MD user or developer and have suggestions for physical validity tests
   missing in this package, we would love to hear from you! Please
   consider getting in touch with us via our `github repository`_.

.. toctree::
   userguide
   :maxdepth: 2
   :caption: User guide:

.. toctree::
   examples/kinetic_energy_distribution
   examples/ensemble_check
   examples/kinetic_energy_equipartition
   examples/integrator_validation
   examples/openmm_replica_exchange
   :maxdepth: 2
   :caption: Examples:

.. toctree::
   simulation_data
   :maxdepth: 2
   :caption: Data format and parsers:

.. toctree::
   physical_validation
   physical_validation.data
   physical_validation.util
   :maxdepth: 2
   :caption: Package reference:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


.. _`github repository`: https://github.com/shirtsgroup/physical-validation
physical_validation\.data subpackage
====================================

physical_validation\.data\.simulation_data module
-------------------------------------------------
.. automodule:: physical_validation.data.simulation_data
    :members:
    :undoc-members:

physical_validation\.data\.unit_data module
-------------------------------------------
.. automodule:: physical_validation.data.unit_data
    :members:
    :undoc-members:

physical_validation\.data\.ensemble_data module
-----------------------------------------------
.. automodule:: physical_validation.data.ensemble_data
    :members:
    :undoc-members:

physical_validation\.data\.trajectory_data module
-------------------------------------------------
.. automodule:: physical_validation.data.trajectory_data
    :members:
    :undoc-members:

physical_validation\.data\.observable_data module
-------------------------------------------------
.. automodule:: physical_validation.data.observable_data
    :members:
    :undoc-members:

physical_validation\.data\.system_data module
---------------------------------------------
.. automodule:: physical_validation.data.system_data
    :members:
    :undoc-members:

physical_validation\.data\.parser module
----------------------------------------
.. automodule:: physical_validation.data.parser

.. autoclass:: physical_validation.data.parser.Parser
    :members:
    :undoc-members:

physical_validation\.data\.gromacs_parser module
------------------------------------------------
.. automodule:: physical_validation.data.gromacs_parser

.. autoclass:: physical_validation.data.gromacs_parser.GromacsParser
    :members:
    :undoc-members:

physical_validation\.data\.flatfile_parser module
-------------------------------------------------
.. automodule:: physical_validation.data.flatfile_parser

.. autoclass:: physical_validation.data.flatfile_parser.FlatfileParser
    :members:
    :undoc-members:
physical_validation package
===========================

.. automodule:: physical_validation
    :members:

physical_validation\.kinetic_energy module
------------------------------------------

.. automodule:: physical_validation.kinetic_energy
    :members:

physical_validation\.ensemble module
------------------------------------

.. automodule:: physical_validation.ensemble
    :members:

physical_validation\.integrator module
--------------------------------------

.. automodule:: physical_validation.integrator
    :members:
