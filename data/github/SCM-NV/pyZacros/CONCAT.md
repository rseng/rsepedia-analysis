![Python application](https://github.com/SCM-NV/pyZacros/workflows/Build/badge.svg?branch=master)

What is pyZacros
----------------

pyZacros (Python Library for Automating Zacros Simulation) is a collection of tools that aims to provide a powerful, flexible, and easily extendable Python interface to Zacros. It is designed as an extension of the python library [PLAMS](https://github.com/SCM-NV/PLAMS). Thereby, pyZacros inherits from PLAMS the robust way of managing the inputs file preparation, job execution, file management, and output file processing. Above and above that, it also offers the possibility of postprocessing the results and building very advanced data workflows.

The normal Zacros workflow has the following steps:

1. First, the subject of the problem (description of the system, and set of the desired simulation parameters) has to be written
   to input text files (i.e. ``energetics_input.dat``, ``mechanism_input.dat``, ``lattice_input.dat``, ``simulation_input.dat``).
2. The ``zacros.x`` program is executed and produces output text files (i.e. ``general_output.txt``, ``procstat_output.txt``,
   ``history_output.txt``, ``specnum_output.txt``, ``lattice_output.txt``).
3. Those output files may contain already the required information or at least contain enough information to get it after
   a postprocessing step.
4. This resultant information may be used to define parameters for further calculations.

pyZacros helps with the automation of all these steps described above directly from a python script and takes responsibility for tiresome and monotonous technical details allowing the user to focus on real science and your problem.

What can be done with pyZacros
------------------------------

As an extension of PLAMS, pyZacros is also designed under the same key design principle ... *flexibility*.
If something (and by something we mean: adjusting an input parameter, executing some program with particular options, extracting a value from output etc.) can be done by hand, it can be done with pyZacros.
The internal structure of the library was designed in a highly modular, especially an object-oriented manner. In particular, there are classes to represent species, clusters, elementary Reactions, among others that are easy to set up and use.

The most important features of pyZacros:

* Preparing, running and examining results of a Zacros jobs from within a single Python script
* Convenient automatic file and folder management
* Running jobs in parallel without a need to prepare a special parallel script
* Integration with popular job schedulers (OGE, SLURM, TORQUE)
* Prevention of multiple runs of the same job
* Easy data transfer between separate runs
* An efficient and simple way for restarting calculation in case of a crash or just to extend the simulation.
* Almost full coverage of all input options and output data in Zacros.
* Visualization and interactive building of the lattice of the system.
* Default plot functions to visualize results like adlayer configurations, process statistics, and species numbers.
* Reconstruction of the pyZacros objects from Zacros calculation which were not managed by pyZacros.

Simple example
--------------

Here we show a simple pyZacros script that reproduces the Zacros tutorial
[Ziff-Gulari-Barshad (ZGB) Model in Zacros](https://zacros.org/tutorials/4-tutorial-1-ziff-gulari-barshad-model-in-zacros). 

The ZGB model includes (see the script below):

1. Three gas species: CO, O2, and CO2. (Lines 5-7)
2. Three surface species: \*, CO\*, O\*. (Lines 10-12)
3. A rectangular lattice with a single site type. (Lines 15-16)
4. Two clusters are included in the cluster-expansion Hamiltonian approach for the energetics. The CO* and O* individual
   adsorbates (without lateral interactions) with 1.3 eV and 2.3 eV, binding energies, respectively. (Lines 20-21)
5. Three irreversible events: non-dissociative adsorption of CO, dissociative adsorption of O2, and fast reaction between
   an O adatom and a CO adsorbate. (Lines 24-31)

```python {.line-numbers}
import scm.plams
import scm.pyzacros as pz

# Gas species:
CO_g = pz.Species("CO")
O2_g = pz.Species("O2")
CO2_g = pz.Species("CO2", gas_energy=-2.337)

# Surface species:
s0 = pz.Species("*", 1)      # Empty adsorption site
CO_s = pz.Species("CO*", 1)
O_s = pz.Species("O*", 1)

# Lattice setup:
lattice = pz.Lattice( lattice_type=pz.Lattice.RECTANGULAR,
                        lattice_constant=1.0, repeat_cell=[10,10] )
lattice.plot()

# Clusters:
CO_p = pz.Cluster( species=[CO_s], cluster_energy=-1.3 )
O_p = pz.Cluster( species=[O_s], cluster_energy=-2.3 )

# Elementary Reactions
CO_ads = pz.ElementaryReaction( initial=[s0, CO_g], final=[CO_s],
                                reversible=False, pre_expon=10.0, activation_energy=0.0 )

O2_ads = pz.ElementaryReaction( initial=[s0, s0, O2_g], final=[O_s, O_s], neighboring=[(0, 1)],
                                reversible=False, pre_expon=2.5, activation_energy=0.0 )

CO_oxi = pz.ElementaryReaction( initial=[CO_s, O_s], final=[s0, s0, CO2_g], neighboring=[(0, 1)],
                                reversible=False, pre_expon=1.0e+20, activation_energy=0.0)

scm.plams.init()

# Settings:
sett = pz.Settings()
sett.temperature = 500.0
sett.pressure = 1.0
sett.snapshots = ('time', 5.e-1)
sett.process_statistics = ('time', 1.e-2)
sett.species_numbers = ('time', 1.e-2)
sett.max_time = 25.0

sett.molar_fraction.CO = 0.45
sett.molar_fraction.O2 = 0.55

myJob = pz.ZacrosJob( settings=sett, lattice=lattice,
                        mechanism=[CO_ads, O2_ads, CO_oxi],
                        cluster_expansion=[CO_p, O_p] )

results = myJob.run()

print( "nCO2 = ", results.provided_quantities()["CO2"][-10:] )
results.plot_molecule_numbers( results.gas_species_names() )
results.plot_molecule_numbers( results.surface_species_names() )

scm.plams.finish()
```

Don't worry if something in the above code is incomprehensible or confusing.
Everything you need to know to understand how pyZacros works and how to write your own scripts is explained
in next chapters of this documentation.

By executing the above script, you are going to see a visual representation of the lattice (see script's line 17) that should
be similar to the image below:

<p align="center">
    <img src="doc/images/ZGB-lattice.png" style="width:40%">
</p>

Then, you should see the plot of the number of molecules of each kind as a function of time during the simulation. We have split this information into two Figures for clarity, one for gas-phase species and the other one for surface species, as follows (see script's line 54-55):

<p align="center">
    <img src="doc/images/ZGB-mol_gas_nums.png" style="width:40%">
    <img src="doc/images/ZGB-mol_surf_nums.png" style="width:40%">
</p>
   
During the execution the following information is written to the standard output:

```
[02.11|12:07:12] PLAMS working folder: /home/user/plams_workdir
[02.11|12:07:12] JOB plamsjob STARTED
[02.11|12:07:12] JOB plamsjob RUNNING
[02.11|12:07:12] JOB plamsjob FINISHED
[02.11|12:07:12] JOB plamsjob SUCCESSFUL
nCO2 = [2825, 2827, 2828, 2829, 2829, 2830, 2830, 2832, 2832, 2834]
[02.11|12:07:40] PLAMS run finished. Goodbye
```
   
It indicates that pyZacros created a uniquely named working folder (``plams_workdir``) and then ran the Zacros calculation in a separate
subfolder of the working folder (``plamsjob``). All the files created by each Zacros run are saved in the corresponding subfolder for future reference. However, notice that you can access the results directly from the python script. To illustrate this, see line 54 of the script that produces line 6 in the standard output, which prints the number of CO2 molecules produced in the last ten-time steps of the simulation.

Further reading
--------------------

You can find full [pyZacros documentation](https://www.scm.com/doc/pyzacros/index.html) hosted on our website, together with some [tutorials](https://www.scm.com/doc/Tutorials/Scripting/Scripting.html).

Thank you for your interest in contributing to the pyZacros project!

# How to contribute

We want to keep it as easy as possible to contribute changes that
get things working in your environment. There are a few guidelines that we
need contributors to follow so that we can have a chance of keeping on
top of things.

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Fork the repository on GitHub

## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the master branch.
  * Only target release branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master; `git checkout -b
    fix/master/my_contribution master`. Please avoid working directly on the
    `master` branch.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before committing.
* Make sure your commit messages are informative.
* Make sure you have added the necessary tests for your changes.
* Run _all_ the tests to assure nothing else was accidentally broken.

## Submitting Changes

* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the repository in the SCM-NV organization.
* The core team looks at Pull Requests on a regular basis.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [GitHub pull request documentation](https://help.github.com/send-pull-requests/)
Citations
*********

When you publish results in the scientific literature that were obtained using the PLAMS library, please include a reference:

| **PLAMS, SCM, Theoretical Chemistry, Vrije Universiteit, Amsterdam, The Netherlands** `<https://www.scm.com>`__, `<https://github.com/SCM-NV/PLAMS>`__

Additionally, you may add a list of authors and contributors:

| **Main author:** Michał Handzlik.
| **Main contributors:** Bas van Beek, Patrick Melix, Robert Rüger, Tomáš Trnka, Lars Ridder, Felipe Zapata.
.. _index:

================================================
Python Library for Automating Zacros Simulations
================================================

.. toctree::
    :caption: Table of contents
    :maxdepth: 2

    intro
    components/components
    examples/examples
    citations
.. _intro:

Introduction
============


What is pyZacros
----------------

pyZacros (Python Library for Automating Zacros Simulation) is a collection of tools that aims to provide a powerful, flexible, and easily extendable Python interface to Zacros. It is designed as an extension of the python library `PLAMS <https://github.com/SCM-NV/PLAMS>`_. Thereby, pyZacros inherits from PLAMS the robust way of managing the inputs file preparation, job execution, file management, and output file processing. Above and above that, it also offers the possibility of postprocessing the results and building very advanced data workflows.

The normal Zacros workflow has the following steps:

1. First, the subject of the problem (description of the system, and set of the desired simulation parameters) has to be written
   to input text files (i.e. ``energetics_input.dat``, ``mechanism_input.dat``, ``lattice_input.dat``, ``simulation_input.dat``).
2. The ``zacros.x`` program is executed and produces output text files (i.e. ``general_output.txt``, ``procstat_output.txt``,
   ``history_output.txt``, ``specnum_output.txt``, ``lattice_output.txt``).
3. Those output files may contain already the required information or at least contain enough information to get it after
   a postprocessing step.
4. This resultant information may be used to define parameters for further calculations.

pyZacros helps with the automation of all these steps described above directly from a python script and takes responsibility for tiresome and monotonous technical details allowing the user to focus on real science and your problem.

What can be done with pyZacros
------------------------------

As an extension of PLAMS, pyZacros is also designed under the same key design principle ... *flexibility*.
If something (and by something we mean: adjusting an input parameter, executing some program with particular options, extracting a value from output etc.) can be done by hand, it can be done with pyZacros.
The internal structure of the library was designed in a highly modular, especially an object-oriented manner. In particular, there are classes to represent species, clusters, elementary Reactions, among others that are easy to set up and use.

The most important features of pyZacros:

* Preparing, running and examining results of a Zacros jobs from within a single Python script
* Convenient automatic file and folder management
* Running jobs in parallel without a need to prepare a special parallel script
* Integration with popular job schedulers (OGE, SLURM, TORQUE)
* Prevention of multiple runs of the same job
* Easy data transfer between separate runs
* An efficient and simple way for restarting calculation in case of a crash or just to extend the simulation.
* Almost full coverage of all input options and output data in Zacros.
* Visualization and interactive building of the lattice of the system.
* Default plot functions to visualize results like adlayer configurations, process statistics, and species numbers.
* Reconstruction of the pyZacros objects from Zacros calculation which were not managed by pyZacros.

.. _simple_example:

Simple example
--------------

Here we show a simple pyZacros script that reproduces the Zacros tutorial
`Ziff-Gulari-Barshad (ZGB) Model in Zacros <https://zacros.org/tutorials/4-tutorial-1-ziff-gulari-barshad-model-in-zacros>`_.

The ZGB model includes (see the script below):

1. Three gas species: CO, O\ :sub:`2`, and CO\ :sub:`2`. (Lines 5-7)
2. Three surface species: \*, CO\*, O\*. (Lines 10-12)
3. A rectangular lattice with a single site type. (Lines 15-16)
4. Two clusters are included in the cluster-expansion Hamiltonian approach for the energetics. The CO* and O* individual
   adsorbates (without lateral interactions) with 1.3 eV and 2.3 eV, binding energies, respectively. (Lines 20-21)
5. Three irreversible events: non-dissociative adsorption of CO, dissociative adsorption of O2, and fast reaction between
   an O adatom and a CO adsorbate. (Lines 24-31)

.. code-block:: python
   :linenos:

   import scm.plams
   import scm.pyzacros as pz

   # Gas species:
   CO_g = pz.Species("CO")
   O2_g = pz.Species("O2")
   CO2_g = pz.Species("CO2", gas_energy=-2.337)

   # Surface species:
   s0 = pz.Species("*", 1)      # Empty adsorption site
   CO_s = pz.Species("CO*", 1)
   O_s = pz.Species("O*", 1)

   # Lattice setup:
   lattice = pz.Lattice( lattice_type=pz.Lattice.RECTANGULAR,
                         lattice_constant=1.0, repeat_cell=[10,10] )
   lattice.plot()

   # Clusters:
   CO_p = pz.Cluster( species=[CO_s], cluster_energy=-1.3 )
   O_p = pz.Cluster( species=[O_s], cluster_energy=-2.3 )

   # Elementary Reactions
   CO_ads = pz.ElementaryReaction( initial=[s0, CO_g], final=[CO_s],
                                   reversible=False, pre_expon=10.0, activation_energy=0.0 )

   O2_ads = pz.ElementaryReaction( initial=[s0, s0, O2_g], final=[O_s, O_s], neighboring=[(0, 1)],
                                   reversible=False, pre_expon=2.5, activation_energy=0.0 )

   CO_oxi = pz.ElementaryReaction( initial=[CO_s, O_s], final=[s0, s0, CO2_g], neighboring=[(0, 1)],
                                   reversible=False, pre_expon=1.0e+20, activation_energy=0.0)

   scm.plams.init()

   # Settings:
   sett = pz.Settings()
   sett.temperature = 500.0
   sett.pressure = 1.0
   sett.snapshots = ('time', 5.e-1)
   sett.process_statistics = ('time', 1.e-2)
   sett.species_numbers = ('time', 1.e-2)
   sett.max_time = 25.0

   sett.molar_fraction.CO = 0.45
   sett.molar_fraction.O2 = 0.55

   myJob = pz.ZacrosJob( settings=sett, lattice=lattice,
                         mechanism=[CO_ads, O2_ads, CO_oxi],
                         cluster_expansion=[CO_p, O_p] )

   results = myJob.run()

   print( "nCO2 = ", results.provided_quantities()["CO2"][-10:] )
   results.plot_molecule_numbers( results.gas_species_names() )
   results.plot_molecule_numbers( results.surface_species_names() )

   scm.plams.finish()


Don't worry if something in the above code is incomprehensible or confusing.
Everything you need to know to understand how pyZacros works and how to write your own scripts is explained
in next chapters of this documentation.

By executing the above script, you are going to see a visual representation of the lattice (see script's line 17) that should
be similar to the image below:

.. image:: ../images/ZGB-lattice.png
   :scale: 60 %
   :align: center

Then, you should see the plot of the number of molecules of each kind as a function of time during the simulation. We have split this information into two Figures for clarity, one for gas-phase species and the other one for surface species, as follows (see script's line 54-55):

.. image:: ../images/ZGB-mol_gas_nums.png
   :scale: 55 %

.. image:: ../images/ZGB-mol_surf_nums.png
   :scale: 55 %

During the execution the following information is written to the standard output:

.. code-block:: none
   :linenos:

   [02.11|12:07:12] PLAMS working folder: /home/user/plams_workdir
   [02.11|12:07:12] JOB plamsjob STARTED
   [02.11|12:07:12] JOB plamsjob RUNNING
   [02.11|12:07:12] JOB plamsjob FINISHED
   [02.11|12:07:12] JOB plamsjob SUCCESSFUL
   nCO2 = [2825, 2827, 2828, 2829, 2829, 2830, 2830, 2832, 2832, 2834]
   [02.11|12:07:40] PLAMS run finished. Goodbye

It indicates that pyZacros created a uniquely named working folder (``plams_workdir``) and then ran the Zacros calculation in a separate
subfolder of the working folder (``plamsjob``). All the files created by each Zacros run are saved in the corresponding subfolder for future reference. However, notice that you can access the results directly from the python script. To illustrate this, see line 54 of the script that produces line 6 in the standard output, which prints the number of CO2 molecules produced in the last ten-time steps of the simulation.
.. _species_and_specieslist:

Species / Species List
----------------------

In a Zacros simulation, species are necessary to describe the chemistry involved in the system.
If one atom is identical to another, we can say they are the same chemical species. In the same way,
if one molecule is identical to another, we can say they are the same chemical species. Thus, we can
highlight two essential kinds of species: 1) gas species and 2) surface or adsorbed species. A
molecule in the gas phase is a species itself, but when it interacts with a catalytic surface,
its properties generally change enough to be different from its gas counterpart. However, once it
is adsorbed on the catalytic surface, its properties may not change enough to change its identity
as it moves on the surface. In that sense, for example, for a CO molecule, we can identify two kinds
of species: 1) CO in the gas phase and 2) CO adsorbed. Here, the CO is the same species independently
if it is adsorbed, e.g., on an fcc or an hcp site.

For any kind of species, the only required parameter is the symbol (e.g., ``"CO"``), and it can be created with the sentence
``pz.Species("CO")``. The ``Species`` constructor allows specifying different parameters like denticity, gas energy, kind, and mass.
By default, pyZacros parses the species symbol to get these parameters. If the symbol contains the character ``*``,
it assumes that the species is a surface species (``kind=pz.Species.SURFACE``) with a denticity given by the number
of times that ``*`` is found in the symbol; i.e. ``pz.Species("O2**")`` is equivalent to
``pz.Species("O2**",denticity=2,kind=pz.Species.SURFACE)``. On the other hand, if the symbol doesn't contain any ``*``,
it assumes that the species is a gas species (``kind=pz.Species.GAS``), with the mass calculated from the parsing of
the symbol as a chemical formula string; i.e. ``pz.Species("O2")`` is equivalent
to ``pz.Species("O2",kind=pz.Species.GAS,mass=31.9898)``.

For our example (see :ref:`use case system <use_case_model_zgb>`), we need to create
three gas species (CO, O\ :sub:`2`, and CO\ :sub:`2`), and three surface species (\*, CO\*, O\*).
This can be achieved by using the lines 1-10 of following code:

.. code-block:: python
  :linenos:

  # Gas species
  CO_g = pz.Species("CO")
  O2_g = pz.Species("O2")
  CO2_g = pz.Species("CO2", gas_energy=-2.337)

  # Surface species
  s0 = pz.Species("*")   # Empty adsorption site
  CO_s = pz.Species("CO*")
  O_s = pz.Species("O*", denticity=1)

  # Species List
  spl = pz.SpeciesList([CO_g,O2_g,CO2_g,s0,CO_s])
  spl.append( O_s )
  print(spl)

Notice that the species symbol ``*`` is reserved for the empty site species (or pseudo-species).

In this script, we have also introduced the class ``SpeciesList``, which represents nothing more
than a list of species (see lines 11-13). But, in particular, it is helpful to look at the Zacros
code that will be generated based on it by using the ``print()`` function (see line 14). The execution
of the previous script displays the following on the standard output:

.. code-block:: none

  n_gas_species    3
  gas_specs_names              CO           O2          CO2
  gas_energies        0.00000e+00  0.00000e+00 -2.33700e+00
  gas_molec_weights   2.79949e+01  3.19898e+01  4.39898e+01
  n_surf_species    2
  surf_specs_names         CO*        O*
  surf_specs_dent            1         1

Please consult Zacros' user guide for more details about the specific meaning of the keywords used in the previous lines.

API
~~~

.. currentmodule:: scm.pyzacros.core.Species
.. autoclass:: Species
  :exclude-members: __eq__, __init__, __hash__, __str__, __weakref__


.. currentmodule:: scm.pyzacros.core.SpeciesList
.. autoclass:: SpeciesList
  :exclude-members: __init__, __hash__, __str__, _SpeciesList__updateLabel
.. _zacrosresults:

ZacrosResults
-------------

For an explanation purpose let us assume that ``/home/user/xyz`` contains three files: ``ammonia.xyz``, ``ethanol.xyz``, ``water.xyz``.
When you run this script the standard output will look something like:

.. code-block:: python
  :linenos:

   results = job.run()

   if( job.ok() ):
      provided_quantities = results.provided_quantities()
      print("nCO2 =", provided_quantities['CO2'])

      results.plot_molecule_numbers( results.gas_species_names() )
      results.plot_lattice_states( results.lattice_states() )

      pstat = results.get_process_statistics()
      results.plot_process_statistics( pstat, key="number_of_events" )

For an explanation purpose let us assume that ``/home/user/xyz`` contains three files: ``ammonia.xyz``, ``ethanol.xyz``, ``water.xyz``.
When you run this script the standard output will look something like:

.. code-block:: none

   [05.11|10:22:27] JOB plamsjob STARTED
   [05.11|10:22:27] JOB plamsjob RUNNING
   [05.11|10:22:27] JOB plamsjob FINISHED
   [05.11|10:22:27] JOB plamsjob SUCCESSFUL
   nCO2 = [0, 28, 57, 85, 118, 139, 161, 184, 212, 232, 264]
   [05.11|10:22:27] PLAMS run finished. Goodbye

.. code-block:: python

   results.plot_molecule_numbers( results.gas_species_names() )

.. image:: ../../images/mol_gas_nums.png
   :scale: 100 %
   :align: center

.. code-block:: python

   results.plot_lattice_states( results.lattice_states() )

.. image:: ../../images/lattice_state.gif
   :scale: 100 %
   :align: center

.. code-block:: python

   pstat = results.get_process_statistics()
   results.plot_process_statistics( pstat, key="number_of_events" )

.. image:: ../../images/number_of_events.gif
   :scale: 80 %
   :align: center

API
~~~

.. currentmodule:: scm.pyzacros.core.ZacrosResults
.. autoclass:: ZacrosResults
    :exclude-members: _ZacrosResults__plot_process_statistics

.. _lattice:

Lattice
-------

The Lattice class defines the lattice structure on which species can bind, diffuse and react. There
are several ways to specify the lattice structure. They are defined in a correspondence one-to-one
with the conventions used in the Zacros' input files. See the API section below for a detailed description
of these three ways: 1) Default Lattices, 2) Unit-Cell-Defined Periodic Lattices, and 3) Explicitly Defined Custom Lattices.

Following our example (see :ref:`use case system <use_case_model_zgb>`), we just need a single-site lattice with
a coordination number of 3, a lattice constant equal to ``1.0``, and a modest number of copies of the unit cell ``10x3``:

.. code-block:: python
  :linenos:

   # Lattice setup
   lat = pz.Lattice( lattice_type=pz.Lattice.TRIANGULAR,
                     lattice_constant=1.0, repeat_cell=[10,3] )

   print(lat)

   lattice.plot()

The previous lines produce the following output:

.. code-block:: none

   lattice default_choice
   triangular_periodic 1.0 10 3
   end_lattice

In addition to the capabilities of building lattices, pyZacros also offers a way to visualize them by calling
the function ``plot()``. e.g., see line 7 of the script above. This line produces the following figure:

.. image:: ../../images/lattice.png
   :scale: 100 %
   :align: center

API
~~~

.. currentmodule:: scm.pyzacros.core.Lattice
.. autoclass:: Lattice
  :exclude-members: __init__, __str__, __weakref__, _Lattice__fromDefaultLattices, _Lattice__fromExplicitlyDefined, _Lattice__fromUnitCellDefined
.. _zacrosjob:

ZacrosJob
---------

The Settings class provides a general-purpose data container for various kinds of information that need to be stored and processed by the pyZacros and PLAMS environment. It is formally identical that the Settings class in PLAMS. Please, see all details here `PLAMS.Job <../../plams/components/jobs.html>`_.


.. code-block:: python
  :linenos:

   job = pz.ZacrosJob( settings=sett, lattice=lat,
                       mechanism=[CO_ads, O2_ads, CO_oxi],
                       cluster_expansion=[CO_p, O_p] )

   print(job)

For an explanation purpose let us assume that ``/home/user/xyz`` contains three files: ``ammonia.xyz``, ``ethanol.xyz``, ``water.xyz``.
When you run this script the standard output will look something like:

.. code-block:: none

   ---------------------------------------------------------------------
   simulation_input.dat
   ---------------------------------------------------------------------
   random_seed         953129
   temperature          500.0
   pressure               1.0

   snapshots                 on time       0.1
   process_statistics        on time       0.1
   species_numbers           on time       0.1
   event_report      off
   max_steps         infinity
   max_time          1.0

   n_gas_species    3
   gas_specs_names              CO           O2          CO2
   gas_energies        0.00000e+00  0.00000e+00 -2.33700e+00
   gas_molec_weights   2.79949e+01  3.19898e+01  4.39898e+01
   gas_molar_fracs     4.50000e-01  5.50000e-01  0.00000e+00

   n_surf_species    2
   surf_specs_names         CO*        O*
   surf_specs_dent            1         1

   finish
   ---------------------------------------------------------------------
   lattice_input.dat
   ---------------------------------------------------------------------
   lattice default_choice
     triangular_periodic 1.0 10 3
   end_lattice
   ---------------------------------------------------------------------
   energetics_input.dat
   ---------------------------------------------------------------------
   energetics

   cluster CO*_0-0
     sites 1
     lattice_state
       1 CO* 1
     site_types 1
     graph_multiplicity 1
     cluster_eng -1.30000e+00
   end_cluster

   cluster O*_0-0
     sites 1
     lattice_state
       1 O* 1
     site_types 1
     graph_multiplicity 1
     cluster_eng -2.30000e+00
   end_cluster

   end_energetics
   ---------------------------------------------------------------------
   mechanism_input.dat
   ---------------------------------------------------------------------
   mechanism

   step CO_adsorption
     gas_reacs_prods CO -1
     sites 1
     initial
       1 * 1
     final
       1 CO* 1
     site_types 1
     pre_expon  1.00000e+01
     activ_eng  0.00000e+00
   end_step

   step O2_adsorption
     gas_reacs_prods O2 -1
     sites 2
     neighboring 1-2
     initial
       1 * 1
       2 * 1
     final
       1 O* 1
       2 O* 1
     site_types 1 1
     pre_expon  2.50000e+00
     activ_eng  0.00000e+00
   end_step

   step CO_oxidation
     gas_reacs_prods CO2 1
     sites 2
     neighboring 1-2
     initial
       1 CO* 1
       2 O* 1
     final
       1 * 1
       2 * 1
     site_types 1 1
     pre_expon  1.00000e+20
     activ_eng  0.00000e+00
   end_step

   end_mechanism

API
~~~

.. currentmodule:: scm.pyzacros.core.ZacrosJob
.. autoclass:: ZacrosJob
    :exclude-members: _result_type, __init__, _get_ready, __str__
.. _latticestate:

LatticeState
------------

A significant result from a KMC simulation is how the different sites in the lattice are populated as a function of time.
During a Zacros simulation, Zacros takes snapshots of the lattice state and writes them in the file ``history_output.txt``.
Parallelly pyZacros stores each snapshot as a ``LatticeState`` object (see ZacrosResults) for further analysis and/or use
on the python side. Another important application of LatticeState objects is that they can be used as initial states, either
as objects from a previous simulation or built from scratch.

For our example (see :ref:`use case system <use_case_model_zgb>`), we will use the ``LatticeState`` class to build a
lattice state from scratch, and use it as the initial state.

.. note::

    pyZacros (as Zacros does), by default, will start with an empty lattice if not stated otherwise.

We are going to make the initial state as a randomly populated lattice by ``CO*`` and ``O*`` with a given coverage:

.. code-block:: python
  :linenos:

   # LatticeState setup (initial state)
   ist = pz.LatticeState(lat, surface_species=spl.surface_species())
   ist.fill_sites_random(site_name='StTp1', species='CO*', coverage=0.1)
   ist.fill_sites_random(site_name='StTp1', species='O*', coverage=0.1)

   print(ist)

   ist.plot()

Similarly to the other classes, the function ``print()`` (see line 6) allows taking a look at the Zacros code that
will be internally generated, which for this example is the following:

.. code-block:: none

      initial_state
        # species * CO* O*
        # species_numbers
        #   - CO*  12
        #   - O*  11
        seed_on_sites CO* 1
        seed_on_sites CO* 4
        seed_on_sites O* 6
        seed_on_sites O* 10
        seed_on_sites O* 20
        seed_on_sites CO* 30
        seed_on_sites CO* 43
        seed_on_sites O* 48
        seed_on_sites O* 52
        seed_on_sites CO* 55
        seed_on_sites O* 58
        seed_on_sites CO* 62
        seed_on_sites CO* 69
        seed_on_sites CO* 70
        seed_on_sites O* 72
        seed_on_sites CO* 73
        seed_on_sites CO* 78
        seed_on_sites CO* 93
        seed_on_sites O* 99
        seed_on_sites O* 106
        seed_on_sites O* 109
        seed_on_sites O* 110
        seed_on_sites CO* 115
      end_initial_state

Please consult Zacros' user guide for more details about the specific meaning of the keywords used in the previous lines.

Finally, to visualize the lattice you can make use of the function ``plot()`` (see line 8). The result is as follows:

.. image:: ../../images/lattice_initial_state.png
   :scale: 100 %
   :align: center

.. note::

    To visualize the previous figure, be sure you have `matplotlib <https://matplotlib.org/>`_ installed.

API
~~~

.. currentmodule:: scm.pyzacros.core.LatticeState
.. autoclass:: LatticeState
   :exclude-members: __init__, __str__, __weakref__, _updateSpeciesNumbers
.. _settings:

Settings
--------

The Settings class provides a general-purpose data container for various kinds of information that need to be stored and processed by the pyZacros and PLAMS environment. It is formally identical that the Settings class in PLAMS. Please, see all details here `PLAMS.Settings <../../plams/components/settings.html>`_.

For our example (see :ref:`use case system <use_case_model_zgb>`), the following lines create the conditions we need for our calculation like e.g., the temperature (500 K) and pressure (1.0 atm.):

.. code-block:: python
  :linenos:

  # Settings:
  sett = pz.Settings()
  sett.random_seed = 953129
  sett.temperature = 500.0
  sett.pressure = 1.0
  sett.snapshots = ('time', 0.1)
  sett.process_statistics = ('time', 0.1)
  sett.species_numbers = ('time', 0.1)
  sett.event_report = 'off'
  sett.max_steps = 'infinity'
  sett.max_time = 1.0

  sett.molar_fraction.CO = 0.45
  sett.molar_fraction.O2 = 0.55

  print(sett)

As in the other pyZacros objects, the function ``print()`` (see line 16) shows the Settings object as it is going to be used
in the Zacros input files. The previous lines produce the following output:

.. code-block:: none

  random_seed         953129
  temperature          500.0
  pressure               1.0

  snapshots                 on time       0.1
  process_statistics        on time       0.1
  species_numbers           on time       0.1
  event_report      off
  max_steps         infinity
  max_time          1.0

Notice that the CO and O\ :sub:`2` molar fractions (keywords ``sett.molar_fraction.CO`` and ``sett.molar_fraction.O2``) are not printed out in the previous configuration block. This is because information about species involved in the mechanism and clusters is needed to generate the corresponding block in the ``simulation_input.dat`` zacros file. This information is going to print out from the :ref:`ZacrosJob <zacrosjob>` object.
Please consult Zacros' user guide for more details about the specific meaning of the keywords used above.

API
~~~

.. currentmodule:: scm.pyzacros.core.Settings
.. autoclass:: Settings
    :exclude-members: __init__, __str__

.. _components_overview:

Components overview
===================

.. note::

    In this documentation, we assume you are familiarized with both projects PLAMS and Zacros. If not, first, please
    take a look at our comprehensive documentation about PLAMS on this
    link: `PLAMS Documentation <../../plams/index.html>`_,
    and the documentation about Zacros on its official web page:
    `Zacros Website <https://zacros.org>`_

This chapter contains a description of all components (classes) that can be used within pyZacros scripts.
The image below shows all classes available in pyZacros.

.. image:: ../../images/architecture.png
   :scale: 60 %
   :align: center

The classes represented in gray boxes are extensions of PLAMS. Settings, ZacrosJob, and ZacrosResults are subclasses
of the PLAMS classes `Settings <../../plams/components/settings.html>`_,
`SingleJob <../../plams/components/jobs.html#scm.plams.core.basejob.SingleJob>`_,
and `Results <../../plams/components/results.html>`_ respectively. Thus, these classes
inherit from PLAMS the robust way of managing the inputs file preparation, job execution, file management,
and output file processing. In a nutshell, the class :ref:`Settings <settings>` is used for establishing the parameters of the
calculation. The :ref:`ZacrosJob <zacrosjob>` class is the primary piece of computational work, and it takes care of running jobs.
Finally, the :ref:`ZacrosResults <zacrosresults>` class takes care of the job results after the execution is finished; it gathers
information about output files, helps to manage them, and extracts data of interest from them.

On the other side, the rest of the classes are specifically designed to define a system in Zacros. The Zacros
package implements a Graph-Theoretical KMC on-lattice methodology coupled with cluster expansion Hamiltonians
for the adlayer energetics and Brønsted-Evans-Polanyi relations for the activation energies of elementary events.
Thus, every system in Zacros needs at least the definition of a set of clusters to evaluate the system's energy
(`ClusterExpansion <clusters.html>`_), a set of elementary events
(`Mechanism <mechanism.html>`_), a lattice representing the catalytic surface
(`Lattice <lattice.html>`_), and possibly an initial state configuration
(`LatticeState <latticestate.html>`_).

.. _use_case_model_zgb:

In the following sections, you can find the API specifications of each particular component, an explanation of its
role in the whole environment, and one example of usage. In particular, we will explain the different components
using as a use case the following example, which is known as the
`Ziff-Gulari-Barshad (ZGB) Model <https://zacros.org/tutorials/4-tutorial-1-ziff-gulari-barshad-model-in-zacros>`_:

.. image:: ../../images/reaction_example.png
   :scale: 60 %
   :align: center

The ZGB model includes:

* Three gas species (CO, O\ :sub:`2`, and CO\ :sub:`2`), and three surface species (\*, CO\*, O\*).
* A lattice with a single site type.
* Two clusters: The CO* and O* individual adsorbates.
* Three irreversible events:

  1) Non-dissociative adsorption of CO
  2) Dissociative adsorption of O\ :sub:`2`
  3) Fast reaction between an O adatom and a CO adsorbate to produce CO\ :sub:`2`

.. toctree::
    :maxdepth: 1

    species
    lattice
    cluster
    mechanism
    latticestate
    settings
    zacrosjob
    zacrosresults
.. _reactions_and_mechanism:

Elementary Reaction / Mechanism
-------------------------------

``ElementaryReaction`` class allows to define reversible and/or irreversible elementary steps for adsorption
of molecules on surface sites, desorption from there, diffusion from one site to a neighboring site, or reactions
between adsorbed particles and gas species. ``ElementaryReaction`` class has a one-to-one relation with the
sections ``steps`` and ``reversible_step`` in Zacros' input files.

For our example (see :ref:`use case system <use_case_model_zgb>`), we need to create
three irreversible events:

1. Non-dissociative adsorption of CO: ``CO(g) + * 🠒 CO*``
2. Dissociative adsorption of O\ :sub:`2`: ``O2(g) + * + * 🠒 O* + O*``
3. Fast reaction between an O adatom and a CO adsorbate to produce CO\ :sub:`2`: ``O* + CO* 🠒 * + * + CO2(g)``

This can be achieved by using the lines 1-10 of following code:

.. code-block:: python
  :linenos:

  # Elementary Reactions
  CO_ads = pz.ElementaryReaction( initial=[s0, CO_g], final=[CO_s],
                                  reversible=False, pre_expon=10.0,
                                  label="CO_adsorption" )

  O2_ads = pz.ElementaryReaction( initial=[s0, s0, O2_g], final=[O_s, O_s],
                                  neighboring=[(0,1)],
                                  reversible=False, pre_expon=2.5,
                                  label="O2_adsorption" )

  CO_oxi = pz.ElementaryReaction( initial=[CO_s, O_s], final=[s0, s0, CO2_g],
                                  neighboring=[(0,1)],
                                  reversible=False, pre_expon=1.0e+20,
                                  label="CO_oxidation")

  mech = pz.Mechanism([O2_ads, CO_ads, CO_oxi])

  print(mech)

Here, we have also introduced the class ``Mechanism``, which is formally a list of ``ElementaryReaction`` objects
(see line 16). But, it allows looking at the Zacros code that will be generated by using the ``print()`` function (see line 18).
The execution of the previous lines show the following on the standard output:

.. code-block:: none

  mechanism

    step O2_adsorption
      gas_reacs_prods O2 -1
      sites 2
      neighboring 1-2
      initial
        1 * 1
        2 * 1
      final
        1 O* 1
        2 O* 1
      site_types 1 1
      pre_expon  2.50000e+00
      activ_eng  0.00000e+00
    end_step

    step CO_adsorption
      gas_reacs_prods CO -1
      sites 1
      initial
        1 * 1
      final
        1 CO* 1
      site_types 1
      pre_expon  1.00000e+01
      activ_eng  0.00000e+00
    end_step

    step CO_oxidation
      gas_reacs_prods CO2 1
      sites 2
      neighboring 1-2
      initial
        1 CO* 1
        2 O* 1
      final
        1 * 1
        2 * 1
      site_types 1 1
      pre_expon  1.00000e+20
      activ_eng  0.00000e+00
    end_step

  end_mechanism

Please consult Zacros' user guide for more details about the specific meaning of the keywords used in the previous lines.

API
~~~

.. currentmodule:: scm.pyzacros.core.ElementaryReaction
.. autoclass:: ElementaryReaction
  :exclude-members: __init__, __eq__, __hash__, _ElementaryReaction__updateLabel, __weakref__, __str__

.. currentmodule:: scm.pyzacros.core.Mechanism
.. autoclass:: Mechanism
  :exclude-members: __init__, __str__
.. _cluster_and_clusterexpansion:

Cluster / Cluster Expansion
---------------------------

In a Zacros simulation, clusters  (also
referred to as patterns) are used as the base of a cluster expansion Hamiltonian for
calculating the energy of a given lattice configuration. The energy is calculated based
on the number of times that each cluster/pattern is found on the catalytic surface during
the simulation. A cluster consists of a collection of binding sites, their connectivity,
the surface species bonded to those sites, and the energetic contribution thereof.

For our example (see :ref:`use case system <use_case_model_zgb>`), the following lines create the
two needed clusters; the CO* and O* individual adsorbates:

.. code-block:: python
  :linenos:

   # Clusters
   CO_p = pz.Cluster( species=[CO_s], cluster_energy=-1.3 )
   O_p = pz.Cluster( species=[O_s], cluster_energy=-2.3 )
   print(CO_p)

which produce the following output:

.. code-block:: none

   cluster CO*_0-0
     sites 1
     lattice_state
       1 CO* 1
     site_types 1
     graph_multiplicity 1
     cluster_eng -1.30000e+00
   end_cluster

Please consult Zacros' user guide for more details about the specific meaning of the keywords used in the previous lines.
Notice that the function ``print()`` in line 4 shows the cluster ``CO_p`` as it is going to be used in the Zacros input files.
The label ``CO*_0-0`` is automatically generated except if the user specifies it by the parameter ``label`` in the constructor.
This label is used as a unique identifier to avoid duplicates.

.. note::

   pyZacros lists are always numbered from 0 to be consistent with the Python language. However, notice that Zacros
   input files require all elements should be numbered from 1; pyZacros takes care internally of this transformation.

The ``ClusterExpansion`` object is formally a list of clusters and as such inherits all properties of Python lists.
The following lines illustrate an example:

.. code-block:: python
  :linenos:

   # Cluster Expansion
   ce = pz.ClusterExpansion([CO_p, O_p])
   print(ce)

which produce the following output:

.. code-block:: none

   energetics

     cluster CO*_0-0
     sites 1
     lattice_state
        1 CO* 1
     site_types 1
     graph_multiplicity 1
     cluster_eng -1.30000e+00
     end_cluster

     cluster O*_0-0
     sites 1
     lattice_state
        1 O* 1
     site_types 1
     graph_multiplicity 1
     cluster_eng -2.30000e+00
     end_cluster

   end_energetics

These lines will be used to create the Zacros input file ``energetics_input.dat``.

API
~~~

.. currentmodule:: scm.pyzacros.core.Cluster
.. autoclass:: Cluster
   :exclude-members: __init__, __len__, __eq__, __hash__, __str__, _Cluster__updateLabel, __weakref__


.. currentmodule:: scm.pyzacros.core.ClusterExpansion
.. autoclass:: ClusterExpansion
   :exclude-members: __init__, __str__
Phase Transitions in the ZGB model & Parameter Continuation.
------------------------------------------------------------
The Ziff-Gulari-Barshad (ZGB) model.
------------------------------------
.. |br| raw:: html

      <br>

Examples
========

In this chapter we present example PLAMS scripts covering various applications, from very simple tasks (like running the same calculation for multiple molecules) to more advanced dynamic workflows.

Most of these example scripts use computational engines from the Amsterdam Modeling Suite, and you will need a license to run them. Contact license@scm.com for further questions.

In order to run the examples, the ``AMSBIN`` environment variable should be properly set. You can test this by typing ``$AMSBIN/plams -h`` in a terminal: this should print PLAMS' help message. If this is not the case (e.g. you get 'No such file or directory'), you need to set up the environmental variable ``$AMSBIN`` (see the `Linux Quickstart guide <../../Installation/Linux_Quickstart_Guide.html>`__ for details).

Simple examples
---------------

.. toctree::
   :hidden:

   o_pt111
   co_tutorial
   zgb

.. |example1| image:: ../../images/example_O+Pt111.png
   :scale: 35 %
   :target: o_pt111.html

.. |example2| image:: ../../images/example_CO-tutorial.png
   :scale: 35 %
   :target: co_tutorial.html

.. |example3| image:: ../../images/example_ZGB.gif
   :scale: 35 %
   :target: zgb.html

.. csv-table::
   :header: |example1|, |example2|, |example3|

   "KMC lattice from first |br| principles. O+Pt(111)", "Water-gas shift reaction on Pt(111)", "Ziff-Gulari-Barshad (ZGB) model"

Advanced examples
-----------------

.. toctree::
   :hidden:

   zgb_pts
   zgb_ss
   zgb_ss_pc

.. |example4| image:: ../../images/example_ZGB-PhaseTransitions.png
   :scale: 32 %
   :target: zgb_pts.html

.. |example5| image:: ../../images/example_ZGB-PhaseTransitions-SS.png
   :scale: 32 %
   :target: zgb_ss.html

.. |example6| image:: ../../images/example_ZGB-PhaseTransitions.png
   :scale: 32 %
   :target: zgb_ss_pc.html

.. csv-table::
   :header: |example4|, |example5|, |example5|

   "ZGB model |br| + Phase Transitions. |br| |br| |br|", "ZGB model |br| + Phase Transitions |br| + Stationary Simulation.  |br| |br|", "ZGB model |br| + Phase Transitions |br| + Stationary Simulation |br| + Parameter Continuation."

.. |br| raw:: html

      <br>

Phase Transitions in the ZGB model: Stationary Simulation.
----------------------------------------------------------

.. _code_plot_coverage_zbg_ss:
.. code-block:: python
  :caption: **Code: Visualizing coverage results**
  :linenos:

  # Lattice states for x_CO=0.54 and x_CO=0.55
  results[33].last_lattice_state().plot()
  results[34].last_lattice_state().plot()


.. |latticeState1| image:: ../../images/example_ZGB-PhaseTransitions-SS-ls1.png
   :scale: 60 %

.. |latticeState2| image:: ../../images/example_ZGB-PhaseTransitions-SS-ls2.png
   :scale: 60 %

.. csv-table:: **Views of the Catalyst Surface**
   :header: |latticeState1|, |latticeState2|

   "A view of the catalyst surface at |br| partial pressure of CO = 0.54. Steady-state.", "A view of the catalyst surface at |br| partial pressure of CO = 0.55. Non-steady-state."


.. _code_plot_mol_num_zgb_ss:
.. code-block:: python
  :caption: **Code: Visualizing Molecule Numbers and Its First Derivative**
  :linenos:

  # Molecule numbers for x_CO=0.54 and x_CO=0.55
  results[33].plot_molecule_numbers( ["CO2"], normalize_per_site=True )
  results[34].plot_molecule_numbers( ["CO2"], normalize_per_site=True )

  # First Derivative. Molecule numbers for x_CO=0.54 and CO=0.55
  results[33].plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )
  results[34].plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )


.. |molnum1| image:: ../../images/example_ZGB-PhaseTransitions-SS-mn1.png
   :scale: 60 %

.. |molnum2| image:: ../../images/example_ZGB-PhaseTransitions-SS-mn2.png
   :scale: 60 %

.. |dmolnum1| image:: ../../images/example_ZGB-PhaseTransitions-SS-dmn1.png
   :scale: 60 %

.. |dmolnum2| image:: ../../images/example_ZGB-PhaseTransitions-SS-dmn2.png
   :scale: 60 %


.. _figure_mol_numbers_zgb_ss:
.. csv-table:: **Molecule Numbers and Its First Derivative**
   :header: |molnum1| |br| |dmolnum1|, |molnum2| |br| |dmolnum2|

   "CO\ :sub:`2` production for CO = 0.54. Steady-state", "CO\ :sub:`2` for CO = 0.55. Non-steady-state"
KMC lattice from first principles. O+Pt(111)
--------------------------------------------
Water-gas shift reaction on Pt(111).
------------------------------------
.. |br| raw:: html

      <br>

Phase Transitions in the ZGB model.
-----------------------------------

This example is inspired in the seminal paper: Kinetic Phase Transitions in an Irreversible Surface-Reaction Model by
Robert M. Ziff, Erdagon Gulari, and Yoav Barshad in 1986 (`Phys. Rev. Lett. 56, 25 <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.56.2553>`_).
This paper is the origin of the ZGB model we have been discussing above. While the model leaves out many important steps of the real system, it exhibits interesting steady-state off-equilibrium behavior and two types of phase transitions, which actually occur in real systems. Please refer to the original paper for more details. In this example, we will analyze the effect of changing the composition of the gas phase, namely partial pressures for O\ :sub:`2` and CO, in the CO\ :sub:`2` Turnover frequency (TOF).

You can download the full example script following this link :download:`ZiffGulariBarshad-PhaseTransitions.py <ZiffGulariBarshad-PhaseTransitions.py>`.

The first part of the script consists of the definition of the system (see :ref:`code <code_def_sys>`). Species, lattice,
cluster expansion, and mechanisms are defined. This is identical to the :ref:`use case system <use_case_model_zgb>`
described along the :ref:`Components overview section<components_overview>`. Thus, please refer to this section for details.

.. _code_def_sys:
.. code-block:: python
  :caption: **Code: Definition of the System**
  :linenos:

  # Gas-species:
  CO_gas = pz.Species("CO")
  O2_gas = pz.Species("O2")
  CO2_gas = pz.Species("CO2", gas_energy=-2.337)

  # Surface species:
  s0 = pz.Species("*", 1)
  CO_ads = pz.Species("CO*", 1)
  O_ads = pz.Species("O*", 1)

  # Lattice setup:
  lattice = pz.Lattice( lattice_type=pz.Lattice.RECTANGULAR,
                            lattice_constant=1.0, repeat_cell=[50,50] )

  # Cluster Expansion:
  CO_point = pz.Cluster(species=[CO_ads], cluster_energy=-1.3)
  O_point = pz.Cluster(species=[O_ads], cluster_energy=-2.3)

  cluster_expansion = [CO_point, O_point]

  # Mechanism:
  CO_adsorption = pz.ElementaryReaction(initial=[s0,CO_gas],
                                        final=[CO_ads],
                                        reversible=False,
                                        pre_expon=10.0,
                                        activation_energy=0.0)

  O2_adsorption = pz.ElementaryReaction(initial=[s0,s0,O2_gas],
                                        final=[O_ads,O_ads],
                                        neighboring=[(0, 1)],
                                        reversible=False,
                                        pre_expon=2.5,
                                        activation_energy=0.0)

  CO_oxidation = pz.ElementaryReaction(initial=[CO_ads, O_ads],
                                       final=[s0, s0, CO2_gas],
                                       neighboring=[(0, 1)],
                                       reversible=False,
                                       pre_expon=1.0e+20,
                                       activation_energy=0.0)

  mechanism = [CO_adsorption, O2_adsorption, CO_oxidation]

The second part corresponds to the calculations settings (see :ref:`code <code_settings>`). It starts with the line ``scm.pyzacros.init()``, which initializes the pyZacros and PLAMS environment. Then, in lines 3 to 7, we configure the parallel execution of the calculations. These lines mean running as many zacros jobs simultaneously as many CPUs are on the system. In particular, line 6 establishes that only one processor will be used for each zacros instance. Lines 9 to 18 are fundamentally the same used in the :ref:`use case system <use_case_model_zgb>` described along the :ref:`Components overview section<components_overview>`. Thus, please refer to this section for details.

.. _code_settings:
.. code-block:: python
  :caption: **Code: Calculation Settings**
  :linenos:

  scm.pyzacros.init()

  # Parallel Settings: Run as many job simultaneously as there are cpu on the system
  maxjobs = multiprocessing.cpu_count()
  scm.plams.config.default_jobrunner = scm.plams.JobRunner(parallel=True, maxjobs=maxjobs)
  scm.plams.config.job.runscript.nproc = 1
  print('Running up to {} jobs in parallel simultaneously'.format(maxjobs))

  # Calculation Settings:
  sett = pz.Settings()
  sett.molar_fraction.CO = 0.45
  sett.molar_fraction.O2 = 0.55
  sett.random_seed = 953129
  sett.temperature = 500.0
  sett.pressure = 1.0
  sett.snapshots = ('time', 0.5)
  sett.species_numbers = ('time', 0.01)
  sett.max_time = 10.0

The next block executes the zacros calculations  (see :ref:`code <code_run>`). Lines 1-2 define the grid of partial pressures of CO to study. In this case, from 0.2 up to 0.8. Line 4 defines the results list, initially empty. From lines 5 up to 14, we have the loop that submits one zacros calculation for each value of CO partial pressure. First, we establish the composition in the settings object by selecting the partial pressure of CO and O\ :sub:`2` (``sett.molar_fraction.CO`` and ``sett.molar_fraction.O2`` respectively. Notice that we assumed that the gas phase is composed only of CO and O\ :sub:`2`. Thus, x\ :sub:`CO` +x\ :sub:`O_2` =1). Lines 9 to 12 initialize the ZacrosJob, and line 14 collects the corresponding results into the ``results`` list. ``job.run()`` will return a ``ZacrosResults`` object. The full loop will execute all jobs in groups of ``maxjobs`` jobs.

.. _code_run:
.. code-block:: python
  :caption: **Code: Running the Calculations**
  :linenos:

  dx = 0.01
  x_CO = numpy.arange(0.2,0.8,dx)

  results = []
  for x in x_CO:
     sett.molar_fraction.CO = x
     sett.molar_fraction.O2 = 1.0-x

     job = pz.ZacrosJob( settings=sett,
                          lattice=lattice,
                          mechanism=mechanism,
                          cluster_expansion=cluster_expansion )

     results.append( job.run() )

Now we move to analyze the results (see :ref:`code <code_results>`). Lines 1 to 3 defines vectors to store important results. Specifically, the average coverage of O and CO species on the surface (``cf_O`` and ``cf_CO`` respectively) and the turnover frequency (TOF) of CO2 (``TOF_CO2``). The loop starting at line 5 fill these vectors by iterating through each element of ``x_CO``. Line 6 is crucial because it forces to wait for the job to finish and checks if the status is successful. Only if both conditions are successful can it proceed to access the results.

Lines 5 to 11, calculate the coverage fractions using the last five lattice states, and line 13 calculates the TOFs by utilizing the function ``get_TOFs()``. Roughly, the TOF is the slope of the regression line for the number of molecules produced as a function of time (we will go in deeper about this concept in the next example). Lines 15 to 17 just save the calculated values into the results vectors, and line 19 waits for all threads to finish and clean the pyZacros and PLAMS environment. Finally, lines 21 to 26 print the results nicely to standard output.

.. _code_results:
.. code-block:: python
  :caption: **Code: Getting the Results**
  :linenos:

  cf_O = []
  cf_CO = []
  TOF_CO2 = []

  for i,x in enumerate(x_CO):
     if( results[i].ok() ):
        acf = { "O*":0.0, "CO*":0.0 }
        for lattice_state in results[i].lattice_states(last=5):
           fractions = lattice_state.coverage_fractions()
           acf["O*"] += fractions["O*"]/5
           acf["CO*"] += fractions["CO*"]/5

        TOFs,_,_ = results[i].get_TOFs()

        cf_O.append( acf["O*"] )
        cf_CO.append( acf["CO*"] )
        TOF_CO2.append( TOFs["CO2"] )

  scm.pyzacros.finish()

  print("----------------------------------------------")
  print("%4s"%"cond", "%8s"%"x_CO", "%10s"%"acf_O", "%10s"%"acf_CO", "%10s"%"TOF_CO2")
  print("----------------------------------------------")

  for i,x in enumerate(x_CO):
     print("%4d"%i,"%8.2f"%x_CO[i],"%10.6f"%cf_O[i],"%10.6f"%cf_CO[i],"%10.6f"%TOF_CO2[i])


If the script work successfully, you would see the following output:

.. _code_output:
.. code-block:: none
  :caption: **Execution: Output**
  :linenos:

  $ amspython ZiffGulariBarshad-PhaseTransitions.py
  [26.11|12:15:51] PLAMS working folder: /home/user/pyzacros/examples/plams_workdir
  Running up to 8 jobs in parallel simultaneously
  [26.11|12:15:51] JOB plamsjob STARTED
  [26.11|12:15:51] JOB plamsjob STARTED
  [26.11|12:15:51] Renaming job plamsjob to plamsjob.002
  [26.11|12:15:51] JOB plamsjob STARTED
  [26.11|12:15:51] Renaming job plamsjob to plamsjob.003
  [26.11|12:15:51] JOB plamsjob STARTED
  [26.11|12:15:51] JOB plamsjob RUNNING
  [26.11|12:15:51] Renaming job plamsjob to plamsjob.004
  [26.11|12:15:51] JOB plamsjob STARTED
  [26.11|12:15:51] JOB plamsjob.002 RUNNING
  ...
  [26.11|12:16:08] JOB plamsjob.057 SUCCESSFUL
  [26.11|12:16:08] JOB plamsjob.056 SUCCESSFUL
  [26.11|12:16:08] JOB plamsjob.058 SUCCESSFUL
  [26.11|12:16:08] JOB plamsjob.059 SUCCESSFUL
  [26.11|12:16:09] JOB plamsjob.060 SUCCESSFUL
  [26.11|12:16:09] JOB plamsjob.061 SUCCESSFUL
  [26.11|12:16:09] JOB plamsjob.062 SUCCESSFUL
  [26.11|12:39:42] PLAMS run finished. Goodbye
  -----------------------------------------
      x_CO      acf_O     acf_CO    TOF_CO2
  -----------------------------------------
      0.20   0.998000   0.000000   0.040744
      0.21   0.999520   0.000000   0.036692
      0.22   1.000000   0.000000   0.042709
      0.23   0.998400   0.000000   0.041491
      0.24   0.997360   0.000000   0.051405
      0.25   0.993360   0.000000   0.074524
      0.26   0.998400   0.000000   0.059448
      0.27   0.997280   0.000000   0.075712
      0.28   0.997440   0.000000   0.085320
      0.29   0.993440   0.000080   0.102368
      0.30   0.993120   0.000000   0.114813
      0.31   0.995040   0.000000   0.120510
      0.32   0.991120   0.000000   0.150822
      0.33   0.988640   0.000000   0.150280
      0.34   0.986640   0.000000   0.210519
      0.35   0.974160   0.000080   0.266200
      0.36   0.961360   0.000240   0.301795
      0.37   0.956160   0.000320   0.341076
      0.38   0.933280   0.000400   0.383555
      0.39   0.925680   0.000320   0.511855
      0.40   0.897760   0.000880   0.551203
      0.41   0.862640   0.002160   0.619353
      0.42   0.867040   0.001280   0.737965
      0.43   0.820560   0.001680   0.881659
      0.44   0.815760   0.002160   0.979467
      0.45   0.743920   0.003680   1.266927
      0.46   0.719840   0.006320   1.311960
      0.47   0.653200   0.011520   1.495406
      0.48   0.648240   0.009360   1.712626
      0.49   0.602320   0.016240   1.847959
      0.50   0.561440   0.020480   2.107661
      0.51   0.540320   0.025440   2.248969
      0.52   0.450880   0.057120   2.500418
      0.53   0.396160   0.078080   2.759625
      0.54   0.073440   0.708800   2.168947
      0.55   0.019040   0.896560   1.873619
      0.56   0.000000   0.998720   0.879270
      0.57   0.000000   1.000000   0.358375
      0.58   0.000000   1.000000   0.225387
      0.59   0.000000   1.000000   0.148030
      0.60   0.000000   1.000000   0.132571
      0.61   0.000000   1.000000   0.085284
      0.62   0.000000   1.000000   0.064224
      0.63   0.000000   1.000000   0.040768
      0.64   0.000000   1.000000   0.036527
      0.65   0.000000   1.000000   0.029231
      0.66   0.000000   1.000000   0.028916
      0.67   0.000000   1.000000   0.022165
      0.68   0.000000   1.000000   0.015293
      0.69   0.000000   1.000000   0.012087
      0.70   0.000000   1.000000   0.011946
      0.71   0.000000   1.000000   0.010444
      0.72   0.000000   1.000000   0.007646
      0.73   0.000000   1.000000   0.006830
      0.74   0.000000   1.000000   0.006555
      0.75   0.000000   1.000000   0.004735
      0.76   0.000000   1.000000   0.004933
      0.77   0.000000   1.000000   0.003422
      0.78   0.000000   1.000000   0.002669
      0.79   0.000000   1.000000   0.003086
      0.80   0.000000   1.000000   0.002969
      0.81   0.000000   1.000000   0.002624

The above results are the final aim of the calculation. However, one can take advantage of python libraries to visualize them. Here, we use matplotlib. Please check the matplotlib documentation for more details at `https://matplotlib.org <https://matplotlib.org>`_. The following lines of code allow visualizing the effect of changing the CO partial pressure on the average coverage of O and CO and the production rate of CO\ :sub:`2`.

.. _code_plot_cov_tof_results:
.. code-block:: python
  :caption: **Code: Visualizing the Coverage and TOF Results**
  :linenos:

  # Coverage and TOF plot
  fig = plt.figure()

  ax = plt.axes()
  ax.set_xlabel('Partial Pressure CO', fontsize=14)
  ax.set_ylabel("Coverage Fraction (%)", color="blue", fontsize=14)
  ax.plot(x_CO, cf_O, color="blue", linestyle="-.", lw=2, zorder=1)
  ax.plot(x_CO, cf_CO, color="blue", linestyle="-", lw=2, zorder=2)
  plt.text(0.3, 0.9, 'O', fontsize=18, color="blue")
  plt.text(0.7, 0.9, 'CO', fontsize=18, color="blue")

  ax2 = ax.twinx()
  ax2.set_ylabel("TOF (mol/s/site)",color="red", fontsize=14)
  ax2.plot(x_CO, TOF_CO2, color="red", lw=2, zorder=5)
  plt.text(0.37, 1.5, 'CO$_2$', fontsize=18, color="red")

  plt.show()


.. _figure_cov_tof_results:
.. image:: ../../images/example_ZGB-PhaseTransitions.png
   :scale: 60 %
   :align: center

This model assumes that when gas-phase molecules of CO and O\ :sub:`2` are adsorbed immediately on empty sites,
and when the 0 and CO occupy adjacent sites, they react immediately. This model is intrinsically irreversible
because the molecules are sticky to their original sites and remain stationary until they are removed by a reaction.
The last figure shows three regions:

1. Oxygen poisoned state, x\ :sub:`CO` <0.32.
2. Reactive state 0.32<x\ :sub:`CO` <0.55.
3. CO poisoned state x\ :sub:`CO` >0.55.

The first transition at x\ :sub:`CO` =0.32 is continuous, and therefore it is of the second order. The second transition at x\ :sub:`CO` =0.55 occurs abruptly, implying that this is of a first-order transition. As you increase the KMC simulation time, the transition becomes more abrupt. We will discuss this effect in the next example.

pyZacros also offers some predefined plot functions that use matplotlib as well. For example, it is possible to see a typical reactive state configuration (x\ :sub:`CO` =0.54) and one in the process of being poisoned by CO (x\ :sub:`CO` =0.55). Just get the last lattice state with the ``last_lattice_state()`` function and visualize it with ``plot()``. See the code and figures below. The state at x\ :sub:`CO` =0.54 is a prototypical steady-state, contrary to the one at x\ :sub:`CO` =0.55, which is otherwise a good example where we can see the two phases coexisting.

.. _code_plot_coverage_zbg_pts:
.. code-block:: python
  :caption: **Code: Visualizing coverage results**
  :linenos:

  # Lattice states for x_CO=0.54 and x_CO=0.55
  results[33].last_lattice_state().plot()
  results[34].last_lattice_state().plot()


.. |latticeState1| image:: ../../images/example_ZGB-PhaseTransitions-ls1.png
   :scale: 60 %

.. |latticeState2| image:: ../../images/example_ZGB-PhaseTransitions-ls2.png
   :scale: 60 %

.. csv-table:: **Views of the Catalyst Surface**
   :header: |latticeState1|, |latticeState2|

   "A view of the catalyst surface at |br| partial pressure of CO = 0.54. Steady-state.", "A view of the catalyst surface at |br| partial pressure of CO = 0.55. Non-steady-state."

In the previous paragraph, we introduced the concept of steady-state. However, let's define it slightly more formally. For our study system, the steady-state for a given composition is characterized when the derivative of the CO2 production (TOF) with respect to time is zero and remains so:

.. math::

  \frac{d}{dt}TOF_{\text{CO}_2} = 0, \,\,\text{for all present and future}\,\, t

pyZacros also offers the function ``plot_molecule_numbers()`` to visualize the molecule numbers and its first derivative as a function of time. See code and figures below:

.. _code_plot_mol_num_zgb_pts:
.. code-block:: python
  :caption: **Code: Visualizing Molecule Numbers and Its First Derivative**
  :linenos:

  # Molecule numbers for x_CO=0.54 and x_CO=0.55
  results[33].plot_molecule_numbers( ["CO2"], normalize_per_site=True )
  results[34].plot_molecule_numbers( ["CO2"], normalize_per_site=True )

  # First Derivative. Molecule numbers for x_CO=0.54 and CO=0.55
  results[33].plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )
  results[34].plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )


.. |molnum1| image:: ../../images/example_ZGB-PhaseTransitions-mn1.png
   :scale: 60 %

.. |molnum2| image:: ../../images/example_ZGB-PhaseTransitions-mn2.png
   :scale: 60 %

.. |dmolnum1| image:: ../../images/example_ZGB-PhaseTransitions-dmn1.png
   :scale: 60 %

.. |dmolnum2| image:: ../../images/example_ZGB-PhaseTransitions-dmn2.png
   :scale: 60 %


.. _figure_mol_numbers_zgb_pts:
.. csv-table:: **Molecule Numbers and Its First Derivative**
   :header: |molnum1| |br| |dmolnum1|, |molnum2| |br| |dmolnum2|

   "CO\ :sub:`2` production for CO = 0.54. Steady-state", "CO\ :sub:`2` for CO = 0.55. Non-steady-state"

From the figures above, it is clear that we have reached a steady-state for x\ :sub:`CO` =0.54. Notice that the first derivative is approximately constant at 2.7 mol/s/site within a tolerance of 5 mol/s/site. Contrary, this is not the case of x\ :sub:`CO` =0.55, where the first derivative continuously decreases.

In the next example, we will modify the script presented here to reach a steady-state configuration for every composition.

As a final note, you can use the following script to visualize the results without running the full calculation:

.. code-block:: python
  :caption: **Code: Visualizing the Results**
  :linenos:

  import scm.pyzacros as pz

  # xCO=0.54
  job = pz.ZacrosJob.load_external( path="plams_workdir/plamsjob.034" )
  job.results.last_lattice_state().plot()
  job.results.plot_molecule_numbers( ["CO2"], normalize_per_site=True )
  job.results.plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )

  # xCO=0.55
  job = pz.ZacrosJob.load_external( path="plams_workdir/plamsjob.035" )
  job.results.last_lattice_state().plot()
  job.results.plot_molecule_numbers( ["CO2"], normalize_per_site=True )
  job.results.plot_molecule_numbers( ["CO2"], normalize_per_site=True, derivative=True )
