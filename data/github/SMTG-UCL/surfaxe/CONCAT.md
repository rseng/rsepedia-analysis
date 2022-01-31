---
title: 'Surfaxe: Systematic surface calculations'
tags:
  - Python
  - materials design
  - high-throughput screening
  - surfaces
authors:
  - name: Katarina Brlec
    orcid: 0000-0003-1485-1888
    affiliation: 1, 2
  - name: Daniel W. Davies
    orcid: 0000-0003-4094-5992
    affiliation: 1, 2 
  - name: David O. Scanlon
    orcid: 0000-0001-9174-8601
    affiliation: 1, 2, 3
affiliations:
 - name: Department of Chemistry, University College London, 20 Gordon Street, London WC1H 0AJ, United Kingdom
   index: 1
 - name: Thomas Young Centre, University College London, Gower Street, London WC1E 6BT, United Kingdom
   index: 2
 - name: Diamond Light Source Ltd., Diamond House, Harwell Science and Innovation Campus, Didcot, Oxfordshire OX11 0DE, UK
   index: 3
date: 01 January 2021
bibliography: paper.bib
---

# Summary 
Surface science is key to understanding the properties of a wide range of materials for energy applications, from catalysts to solar cells to battery components. Computational modelling based on quantum mechanics is often used to calculate surface properties of materials, which in turn determine their stability and performance. The maturity of these "first-principles" methods, coupled with the huge amount of computational power accessible today, means they can now be used predictively in high-throughput screening workflows to suggest new materials for specific applications before they are synthesised. The `surfaxe` package provides a framework for such screening workflows, automating each stage of the process. 

# Statement of need 
Accurate descriptions of electronic structure are needed to calculate surface properties including surface formation energies, adsorption energies and absolute electron energies (ionisation potential, electron affinity, and work function). However, surface calculations diverge significantly from the typical setup in periodic codes where the bulk crystal is described by repeat boundary conditions in all three dimensions. To reveal a surface, the bulk must be cleaved into a slab and periodicity reduced to just two dimensions. Additional parameters must then be considered such as variance in slab thickness, vacuum size, surface termination (dangling bonds), and net electrostatic dipole moments, all of which complicate the calculation workflow and make reliable determination of properties more difficult.

# Surfaxe
The aims of `surfaxe` are:

- To act as a framework for the automation of surface calculations, with particular emphasis on ensuring that properties are converged with respect to the additional parameters that are introduced compared to bulk calculations
- To increase the efficiency and reproducibility of surface calculations by automating the generation of input files and processing of output files for density functional theory (DFT) codes
- To provide a toolbox of intuitive analytical tools to calculate performance-critical materials properties and directly generate publication-quality plots

The code makes extensive use of existing Python Materials Genomics (pymatgen) [@pymatgen] surface modules with full functionality retained in `surfaxe`. As well as a fully flexible Python API, `surfaxe` has a lightweight command line interface.
The modularity of `surfaxe` closely follows a best-practice workflow for the calculation of surface properties, with key features including:

- Automatic cleaving of slabs from the bulk crystal and organising them into a directory structure with all necessary calculation input files -- generation module
- Analyses of atomic displacements and coordination environments, bond lengths and electrostatic potential through the slab (Figure 1a and b) -- analysis module
- Processing of raw DFT outputs to determine surface energy variation with slab and vacuum thickness (Figure 1c) -- convergence module
- Automatic extraction of surface energy, vacuum and core energy levels, along with the important calculation parameters -- data module

In addition to pymatgen, existing packages related to surface calculations include the Atomic Simulation Environment (ase) [@ase], which is a large materials informatics library, and smaller packages to aid with specific post-processing tasks: MacroDensity for plotting of potentials [@macrodensity], WullfPack for plotting of Wulff shapes [@wulffpack], and bapt for plotting band alignments [@bapt]. While these toolkits are extremely useful, `surfaxe` is distinct with its focus  on the rigorous convergence of properties, the enabling of reproducible workflows, and the production of processed datasets and plots at the command line. Lastly, `surfaxe` is built on the pymatgen ecosystem, so full integration with the workflow packages FireWorks [@fireworks] and AiiDA [@aiida] is possible for managing calculations on high-performance computing clusters. 

![Example analysis: a) average bond length, b) electrostatic potential as a function of lattice parameter perpendicular to the surface, and c) a typical surface energy convergence plot with respect to slab and vacuum thickness. \label{fig1}](figures/joss_fig1.png)

# Acknowledgements
The development of this code has benefited from useful discussions with Seán Kavanagh, Graeme W. Watson, Luisa Herring-Rodriguez, Christopher N. Savory, Bonan Zhu, and Maud Einhorn. KB, DWD, and DOS acknowledge support from the European Research Council, ERC, (Grant 758345).



# References
[![Build status](https://github.com/smtg-ucl/surfaxe/actions/workflows/tests.yml/badge.svg)](https://github.com/SMTG-UCL/surfaxe/actions)
[![Documentation Status](https://readthedocs.org/projects/surfaxe/badge/?version=latest)](https://surfaxe.readthedocs.io/en/latest/?badge=stable) 
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03171/status.svg)](https://doi.org/10.21105/joss.03171)

<img src='example_data/figures/surfaxe_header_v2.png' alt='surfaxe logo header' width='600'/>

Calculating the surface properties of crystals from first principles typically introduces several extra parameters including slab thickness, vacuum size, Miller index, surface termination and more.
These factors all influence key properties of interest, making it a challenge to carry out simulations repeatably and draw reliable conclusions.
Surfaxe is a Python package for automating and simplifying density functional theory (DFT) calculations of surface properties, as well as providing analytical tools for bulk and surface calculations.

The code is organised according to the best-practice workflow below:

<img src="example_data/figures/surfaxe_workflow_v4.png" alt="drawing" width="600"/>

The main features include:

1. **Slab generation:** Automated generation of surface slabs from command line.

  * All unique zero-dipole symmetric terminations of slabs are cleaved from a bulk structure.
  * Slabs can be organised into separate folders, optionally with all the required input files to run each calculation.

2. **Raw data processing:** Extracting data from convergence tests.

  * Parsing the convergence testing folders created using the slab generation scripts.
  * Producing dataframes and csv files of summary data.
  * Plotting scripts visualising convergence with respect to slab and vacuum thickness.

3. **Analysis:** Various scripts for surface and bulk calculations.

  * Calculation of planar and macroscopic average of the electrostatic potential through the slab to determine absolute electron energies (ionisation potential, electron affinity).
  * Nearest neighbour atom determination and bond distance analysis (useful for geometry relaxation convergence checks).

Surfaxe primarily supports the [VASP](https://www.vasp.at/) DFT code, however most of the `generation` module is code-agnostic. In the future we would like to add support for more periodic codes in the other modules.

## Example outputs

**Analysis of average bond lengths as a function of slab thickness**
![Bond analysis example](example_data/figures/bond_analysis_plot.png)

**Surface energy convergence checks with respect to vacuum and slab thickness**
![Surface energy convergence example](example_data/figures/010_surface_energy.png)

See the [tutorials directory](tutorials/) for more examples.

## Installation

Surfaxe is a Python 3 package and requires pymatgen and other standard scientific Python packages.

Recommended installation is to git clone and install with `pip` in a stable virtual environment:

```sh
git clone https://github.com/SMTG-UCL/surfaxe.git
cd surfaxe
pip install -e .
```

The `-e` option creates links to the source folder so any changes to the code are reflected on the path.

For the code to generate VASP input files along with the surface slabs, POTCARs need to be [set up with pymatgen](https://pymatgen.org/installation.html#potcar-setup).

## Usage

### Quick start

Surfaxe can be used via the command line and via Python API. [The docs](https://surfaxe.readthedocs.io/en/latest/) include information on both, and the built-in `-h` option is available in the command line interface for each of the scripts.

We recommend starting off by looking at the dedicated [tutorials](https://github.com/SMTG-UCL/surfaxe/tree/master/tutorials). These Jupyter notebooks will guide you through most of the functionality of the package.

The tutorials can also be run interactively on Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SMTG-UCL/surfaxe/HEAD?filepath=tutorials)

### Command line interface

The scripts can be separated into four modules that follow a typical surfaces workflow; these are generation, convergence, analysis and data, with added plotting functionality. The vast majority of `surfaxe` functionality is available via the command line interface, but the python API allows for more flexibility.

Generation:

* `surfaxe-generate`: Generates all unique symmetric zero-dipole surface slabs for one or more specified Miller indices or up to a maximum Miller index specified, in any format supported by pymatgen. Optionally provides all VASP input files.

Convergence:

* `surfaxe-parse-energies`: Extracts the relevant data from the convergence folders set up with `surfaxe-generate` where calculations were run with VASP. Plots convergence graphs of the variation of surface energy with respect to slab and vacuum thickness. Can optionally parse core atom and vacuum energies as well.
* `surfaxe-parse-structures`: Collects the structures' metadata into a json file and optionally performs bond analysis as in `surfaxe-bonds`.

Analysis:

* `surfaxe-bonds`: Parses the structure, looking for bonds between specified atoms.
* `surfaxe-simplenn`: Predicts the coordination environment of atoms for simple structures.
* `surfaxe-complexnn`: Predcts the coordination environment of atoms in more complex structures where the default prediction algorithm fails.
* `surfaxe-potential`: Calculates the planar potential of the slab along c-axis, the gradient of the planar potential and optionally macroscopic potential.
* `surfaxe-cartdisp`: Calculates the Cartesian displacements of atoms during relaxation from intial and final structures.

Plotting:

* `surfaxe-plot-surfen` and `surfaxe-plot-enatom`: Plot the surface energy and energy per atom based on data from `surfaxe-parsefols` with individual customisability
* `surfaxe-plot-bonds`: Plots the bond distance with respect to fractional coordinate, based on `surfaxe-bonds`
* `surfaxe-plot-potential`: Plots the planar (and macroscopic) potential based on data already analysed with `surfaxe-potential`

Data:

* `surfaxe-core`: Collects the core energy level from the middle of a surface slab, based on supplied bulk core atom and the list of its nearest neighbours.
* `surfaxe-vacuum`: Collects the vacuum potential level from a VASP LOCPOT file.

Pymatgen issues warnings whenever the hash in a VASP POTCAR present does not match the one in their library. Only one warning of the same type will be issued by default. All warnings can be suppressed completely by adding the following to your script:

```python
import warnings
warnings.filterwarnings('ignore')
```

## Development notes

### Bugs, features and questions

Please use the Issue Tracker to report bugs or request features in the first instance.

Contributions to interface with this package are most welcome. Please use the ["Fork and Pull"](https://guides.github.com/activities/forking/) workflow to make contributions and stick as closely as possible to the following:

* Code style should comply with [PEP8](http://www.python.org/dev/peps/pep-0008) where possible. [Google's house style](https://google.github.io/styleguide/pyguide.html)
is also helpful, including a good model for docstrings.
* Please use comments liberally!
* Add tests wherever possible, and use the test suite to check if everything still works.

### Tests

Unit tests are in the `tests` directory and can be run from the top directory using `pytest`. Please run these tests whenever submitting bug fix pull requests and include new tests for new features as appropriate.

We also use CI build and testing using [GitHub Actions](https://github.com/SMTG-UCL/surfaxe/actions).

## License and how to cite

Surfaxe is free to use under the MIT License. If you use it in your research, please cite
> K. Brlec, D. W. Davies and D. O. Scanlon, *Surfaxe: Systematic surface calculations.* Journal of Open Source Software, 6(61), 3171, (2021) [DOI: 10.21105/joss.03171](https://joss.theoj.org/papers/10.21105/joss.03171)

## Dependencies

Surfaxe relies primarily on [Pymatgen](pymatgen.org) for manipulating crystal structures and interfacing with the DFT codes.

## Detailed requirements

Surfaxe is compatible with python 3.7+ and requires the following packages:

* [Pymatgen](https://pymatgen.org/)
* [Pandas](https://pandas.pydata.org/)
* [Matplotlib](https://matplotlib.org/)
* [Numpy](https://numpy.org/)
* [Scikit-learn](https://scikit-learn.org/stable/)

## Contributors

List of contributors:

* Katarina Brlec (@brlec)
* Daniel Davies (@dandavies99)
* David Scanlon (@scanlond)

Acknowledgements:

* Surfaxe has benefited from useful discussions with Adam Jackson (@ajjackson), Seán Kavanagh (@kavanase), Graeme Watson (@wantsong), Luisa Herring-Rodriguez (@zccalgh), Christopher Savory (@cnsavory), Bonan Zhu (@zhubonan) and Maud Einhorn (@maudeinhorn).
* Thanks to Keith Butler (@keeeto) for providing a starting point and validation examples for calculating the planar average electrostatic potential via the [macrodensity](https://github.com/WMD-group/MacroDensity) package.
* We thank @eihernan, @pzarabadip and @danielskatz for taking the time to review the code and offering valuable suggestions for improvements.
surfaxe.convergence module
==========================

.. automodule:: surfaxe.convergence
    :members:
    :undoc-members:
    :show-inheritance:

Command Line Interface (CLI)
============================

While :mod:`surfaxe` has a full python API (see our tutorials page for example usage) it also has an
intuitive command line interface (CLI). Below are some simple examples of what you can do with
:mod:`surfaxe` at the command line. 

You can get a full list of accepted flags and what they do for each command using 
the :mod:`--help` or :mod:`-h` flag, e.g.:

.. code:: 

    $ surfaxe-bonds -h

    > usage: surfaxe-bonds [-h] [-s STRUCTURE] [-b BOND [BOND ...]]
                     [--oxi-list OX_STATES_LIST [OX_STATES_LIST ...]]
                     [--oxi-dict OX_STATES_DICT] [--no-csv]
                     [--csv-fname CSV_FNAME] [--no-plot]
                     [--plt-fname PLT_FNAME] [-c COLOR] [--width WIDTH]
                     [--height HEIGHT] [--dpi DPI] [--yaml]

    Parses the structure looking for bonds between atoms. Check the validity of
    the nearest neighbour method on the bulk structure before using it on slabs.

    optional arguments:
    -h, --help            show this help message and exit
    -s STRUCTURE, --structure STRUCTURE
                          Filename of structure file in any format supported by
                          pymatgen (default: POSCAR
    -b BOND [BOND ...], --bond BOND [BOND ...]
                          List of elements e.g. Ti O for a Ti-O bond
    --oxi-list OX_STATES_LIST [OX_STATES_LIST ...]
                          Add oxidation states to the structure as a list 
                          e.g. 3 3 -2 -2 -2
    --oxi-dict OX_STATES_DICT
                          Add oxidation states to the structure as a dictionary
                          e.g. Fe:3,O:-2
    --no-csv              Prints data to terminal
    --csv-fname CSV_FNAME
                          Filename of the csv file (default: bond_analysis.csv)
    --no-plot             Turns off plotting 
    --plt-fname PLT_FNAME
                          Filename of the plot (default: bond_analysis.png)
    -c COLOR, --color COLOR
                          Color of the marker in any format supported by mpl
                          e.g. "#eeefff" hex colours starting with # need to be
                          surrounded with quotation marks
    --width WIDTH         Width of the figure in inches (default: 6)
    --height HEIGHT       Height of the figure in inches (default: 5)
    --dpi DPI             Dots per inch (default: 300)
    --yaml YAML           Read all args from a yaml config file. Completely 
                          overrides any other flags set

The behaviour of default parameters of the functions is extensively documented in 
the :mod:`surfaxe` python package section of the docs. 

=======================
Pre-processing commands
=======================

**surfaxe-generate-slabs**: Generates all unique slabs with specific Miller indices or 
up to a maximum Miller index for a set of slab and vacuum thicknesses. 

Example: :mod:`surfaxe-generate -s bulk_structure.cif --hkl 1,1,0 -t 20 40 -v 20 40 -f` generates
all slabs for the (1,1,0) direction for minimum slab and vacuum thicknesses of 20 Å and 40 Å. 
The :mod:`-f` option organises these into subdirectories with all required VASP input 
files required to run singleshot calculations uisng default settings. It includes all combinations 
for zero-dipole terminations with inversion symmetry. 
The directory structure produced is:

.. code::

    100/              <-- Miller index
      ├── 20_20_0/    <-- slab-thickness_vacuum-thickness_termination-number
      ├── 20_40_0/   
      ├── 40_20_0/
      └── 40_40_0/
        ├── POSCAR    <-- VASP files 
        ├── INCAR
        ├── POTCAR
        └── KPOINTS

*Note: The hkl flag must be comma-separated with no spaces and the list of thicknesses and 
vacuums must be space-separated.*

*Note: To use the :mod:`-f` option you must first set up the 
`pymatgen POTCAR environment <https://pymatgen.org/installation.html#potcar-setup>`_.* 

Similarly, to above the script can be modified to consider multiple Miller indices. 

Example: :mod:`surfaxe-generate -s bulk_structure.cif --hkl 1,1,0 1,1,1 -t 20 40 -v 20 40 -f` 
generates all (1,1,0) and (1,1,1) slabs with minimum slab and vacuum thicknesses of 20 Å and 40 Å. 

*Note: h,k,l are comma-separated with no spaces, while the two (or more) Miller indices are space-separated.*

Lastly, a maximum hkl value can be supplied as an integer so that the script finds all 
zero-dipole slabs up to that maximum Miller index. 

Example: :mod:`surfaxe-generate -s SnO2.cif --hkl 2 -t 20 40 -v 30` generates all slabs with Miller 
indices up to a maximum value of 2, with minimum slab thicknesses of 20 Å and of 40 Å, and 
minimum vacuum thickness of 30 Å. 

========================
Post-processing commands
========================

**surfaxe-parse-energies**: Parses data produced by electronic structure codes once calculations
have been run in then directory structures produced by the pre-processing commands. 
Can optionally collect vacuum and core energies. 

Example: :mod:`surfaxe-parse-energies --hkl 0,0,1 -b 8.83099` saves a csv file of surface energies
and energies per atom for each slab-vacuum combination. See the Tutorials directory for examples. 

**surfaxe-plot-surfen** and **surfaxe-plot-enatom** can be used to customise the surface 
energy and energy per atom plots independetnly based on the data already collated 
with **surfaxe-parse-energies**. 

**surfaxe-parse-structures**: Parses the (relaxed) structures from convergence calculations
and collates them into the same json format as is created when surface slabs are generated. Can 
optionally perform bond analysis for multiple specified bonds. Useful for comparison of relaxed 
and unrelaxed surfaces slabs and determination of convergence.

=================
Analysis commands
=================

**surfaxe-potential**: Reads the local electrostatic potential file and plots the planar 
and macroscopic averages normal to the surface. Currently only the VASP LOCPOT 
file is supported as input. 

Example: :mod:`surfaxe-potential -l LOCPOT -v 11.5` produces a plot assuming a lattice vector of 
11.5 Angstroms and saves the plot data to a csv file. 

**surfaxe-bonds**: Analyse bonding in the structure using Pymatgen's local_env module.
Average bond lengths for each pair of species of interest can be plotted as a function 
of c lattice vector (normal to the slab surface). This can be useful for checking whether
the center of the slab has converged, where bond distances should be bulk-like. 

Example: :mod:`surfaxe-bonds -s CONTCAR -b Sn O` plots the average Sn-O bond length from the 
VASP output structure file. A csv file of the data plotted is also produced. 

**surfaxe-plot-potential** and **surfaxe-plot-bonds** can be used to generate the  
plots based on the data collated with **surfaxe-potential** and **surfaxe-bonds**, 
allowing customisation of plots without having to re-analyse the data. All plotting 
functionality is accessible through the main functions as well. 

**surfaxe-simplenn** and **surfaxe-complexnn**: Analyse the bonding in the slab, again using Pymatgen 
functions. *simplenn* is faster, but less reliable for systems with more complex bonding.
*complexnn* is more robust but requires a dictionary of cutoff bond lengths to be supplied
for each pair of species. See the analysis tutorial for further explanation. 

Example: :mod:`surfaxe-complexnn -s CONTCAR_bivo4 -b Bi3+ O2- 2.46 V5+ O2- 1.73` will 
analyse the coordination of atoms in this BiVO4 slab and save them to a csv file. 

=============
Data commands
=============

There are some simple convenience commands that can also be used to extract key values from
raw data files produced by solid state codes. Currently only commands relating to VASP output
files are included, which rely on the surfaxe :mod:`vasp_data` module. We hope to expand this
in the future. 

**surfaxe-vacuum** and **surfaxe-core** can be used to extract vacuum and core energies, respectively, 
that are needed to calculate absolute electron energies (ionisation potential and electron affinity). 
See the `Macrodensity <https://www.github.com/WMD-group/macrodensity>`_ tutorials for more information
on the steps needed to do this. 

================
YAML input files
================

Most CLI commands allow use of YAML input files containing all the arguments which cannot be 
used in conjunction with other command line argument flags. This is done by specifying 
the :mod:`--yaml` flag which overrides any other flags set in command line by loading the 
:mod:`surfaxe_config.yaml` file.

Sample YAML input files for each of the functions, with defaults and comments are in 
the :mod:`surfaxe/cli/templates` folder. 
All :mod:`**kwargs` of the main function can be passed in the YAML file.  

Example: Generation of (1,0,1) CdTe slabs could easily customised so that all VASP 
input files are created with specific INCAR tags using the following config.yaml file: 

.. code-block:: yaml 

    structure: CdTe.cif 
    hkl: (1,0,1) 
    thicknesses: [20, 40] 
    vacuums: [20, 40] 
    make_fols: True 
    make_files: True 
    max_size: 500 
    center_slab: True 
    ox_states: 
      Cd: 2
      Te: -2
    fmt: poscar 
    name: POSCAR 
    config_dict: PBE_config.json 
    user_incar_settings: 
      ENCUT: 460
      KPAR: 3
    user_kpoints_settings: 
      reciprocal_density: 35

The slabs would then be generated using :mod:`surfaxe-gethkl --yaml config.yaml`Tutorials
=========

We recommend starting off by looking at the `dedicated tutorials. <https://github.com/SMTG-UCL/surfaxe/tree/master/tutorials>`_ 
These Jupyter notebooks will guide you through most of the functionality of the package. 

The tutorials can also be run interactively `on Binder. <https://mybinder.org/v2/gh/SMTG-UCL/surfaxe/HEAD?filepath=tutorials>`_


================================
Using configuration dictionaries
================================

One of the most powerful parts of surfaxe is its ability to make all VASP input 
files needed for convergence testing. To do so surfaxe makes use of configuration
dictionaries (config dicts for short). These are python dictionaries that contain 
information used to set up INCAR, KPOINTS and POTCAR files. 

For example, if we were interested in setting up a single shot PBEsol calculation 
on SnO2 slabs, we could set up the config dict as follows:

.. code-block:: python

    config_dict = {
    "INCAR": {
        "ALGO": "Normal",
        "EDIFF": 1e-06,
        "EDIFFG": -0.01,
        "ENCUT": 500,
        "GGA": "PS",
        "ISMEAR": 0,
        "ISYM": 2,
        "IWAVPR": 1,
        "LASPH": true,
        "LORBIT": 11,
        "LREAL": "auto",
        "NELM": 200,
        "NSW": 0,
        "PREC": "Accurate",
        "SIGMA": 0.02
    },
    "KPOINTS": {
        "reciprocal_density": 55
    },
    "POTCAR": {
        "Sn": "Sn_d", 
        "O" : "O"
    }
    }

Alternatively, one of the ready-made :mod:`surfaxe` config dicts (:mod:`PBEsol.json`, 
:mod:`PBEsol_relax.json`, :mod:`PBE.json`, :mod:`PBE_relax.json` or :mod:`HSE06.json`) 
can be used and further modified using :mod:`user_incar_settings`, 
:mod:`user_kpoints_settings` and :mod:`user_potcar_settings`. The :mod:`relax` config dicts 
contain additional parameters necessary for geometric relaxations of slabs. 
The POTCAR functional (i.e. PBE, PBE_54) can be chosen with :mod:`user_potcar_functional`.  

`Pymatgen documentation <https://pymatgen.org/pymatgen.io.vasp.sets.html#pymatgen.io.vasp.sets.DictSet>`_ 
covers exact behaviour of the :mod:`user_incar_settings`, :mod:`user_kpoints_settings` and :mod:`user_potcar_settings` 
and all additional keyword arguments that can be supplied to slab generation scripts.  
surfaxe.io module
==========================

.. automodule:: surfaxe.io
    :members:
    :undoc-members:
    :show-inheritance:
.. _surfaxe_module:

surfaxe Python package
======================

Please peruse the submodules at your leisure

 .. automodule:: surfaxe
     :members:
     :undoc-members:
     :show-inheritance:

Submodules
----------

.. toctree::

  surfaxe.generation  
  surfaxe.convergence
  surfaxe.analysis
  surfaxe.io
  surfaxe.vasp_data
surfaxe.generation module
==========================

.. automodule:: surfaxe.generation
    :members:
    :undoc-members:
    :show-inheritance:
surfaxe.analysis module
==========================

.. automodule:: surfaxe.analysis
    :members:
    :undoc-members:
    :show-inheritance:
surfaxe.vasp_data module
==========================

.. automodule:: surfaxe.vasp_data
    :members:
    :undoc-members:
    :show-inheritance:
Welcome to surfaxe!
===================

View the code on Github `here <http://github.com/SMTG-UCL/surfaxe>`_.

Contents:

.. toctree::
   :maxdepth: 4
   
   introduction
   command_line_examples
   surfaxe
   tutorials



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Getting started 
===============

============
Introduction
============

The purpose of :mod:`surfaxe` is to simplify the calculation 
of surface properties of inorganic compounds using first-principles codes. 
The package includes pre-processing tools to cleave surfaces from the bulk, automatically 
set up calculation directories and facilitate convergence testing, as well as
post-processing tools to calcualte surface energies, electrostatic potentials, and produce
publication-qaulity plots. 

============
Installation
============

Installation will be via :mod:`pip` in the near future. For now, please clone the repository 
and install the latest stable version:

.. code:: bash

    git clone https://github.com/SMTG-UCL/surfaxe.git
    cd surfaxe
    pip install -e .

The :mod:`-e` is optional and will install the project in developer (editable) mode.

For the code to generate VASP input files along with the surface slabs, 
`pymatgen POTCAR environment <https://pymatgen.org/installation.html#potcar-setup>`_
must be set up correctly. 

=====
Usage
=====

In general, there are two ways to use :mod:`surfaxe`: 
(i) at the command line or (ii) using python in scripts or notebooks. 

Take a look at our `command line interface (CLI) examples <command_line_examples.html>`_ for an overview
of the CLI tools available. 

There is full documentation for all modules (have a browse of the side bar on the left)  
if you would rather use the python API directly. 
