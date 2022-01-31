<p align="center"> 
<img src="https://github.com/symmy596/Polypy/blob/master/docs/source/Figures/polypy_1.png?raw=true"/>
</p>

[![status](https://joss.theoj.org/papers/e17ff370f6ef5fa95bea0fea24cb856c/status.svg)](https://joss.theoj.org/papers/e17ff370f6ef5fa95bea0fea24cb856c)
[![PyPI version](https://badge.fury.io/py/polypy.svg)](https://badge.fury.io/py/polypy)
[![Build Status](https://travis-ci.com/symmy596/PolyPy.svg?branch=master)](https://travis-ci.com/symmy596/PolyPy)
[![Build status](https://ci.appveyor.com/api/projects/status/eo426m99lmkbh5rx?svg=true)](https://ci.appveyor.com/project/symmy596/polypy)
[![Documentation Status](https://readthedocs.org/projects/polypy/badge/?version=latest)](https://polypy.readthedocs.io/en/latest/?badge=latest)
<a href='https://coveralls.io/github/symmy596/PolyPy?branch=master'><img src='https://coveralls.io/repos/github/symmy596/PolyPy/badge.svg?branch=master' alt='Coverage Status' /></a> 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4568493.svg)](https://doi.org/10.5281/zenodo.4568493)

This is the documentation for the open-source Python project - `polypy`,
A library designed to facilitate the analysis of [DL_POLY](https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx) and [DL_MONTE](https://www.ccp5.ac.uk/DL_MONTE) calculations.
polypy is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with.
It is hoped that this tool will bring some benfits to the solid state community and facilitate data analysis and the generation of publication ready plots (powered by Matplotlib.)

The main features include:

1. **Method to analyse the number denisty of a given species in one and two dimensions.**  

   - Generate a plot of the total number of species in bins perpendicular to a specified direction.  
   - Generate a plot of the total number of species in cuboids parallel to a specified direction.  

2. **Method calculate the charge density from the number density.**  

   - Convert number densities of all species in bins perpendicular to a specified direction into the charge density.  

3. **Calculate the electric field and electrostatic potential from the charge density.**  

   - Solves the Poisson Boltzmann equation to convert the charge density into the electric field and the electrostatic potential.

4. **Calculate the diffusion coefficient for a given species from a mean squared displacement.**

   - Carries out a mean squared displacement on an MD trajectory.
   - Calculates the diffusion coefficient.
   - Uses the density analysis and the diffusion coefficient to calculate the ionic conductivity. 
   

<p align="center"> 
<img src="https://github.com/symmy596/PolyPy/blob/master/docs/source/Figures/Show_off.png?raw=true" width="100%"/>
<figcaption>(a) Mean squared displacement for calcuim fluoride. (b) System volume of calcium fluoride during a molecular dynamics simulation. (c) Cerium and oxygen number density at a grain boundary. (d) Electrostatic potential across a grain boundary of cerium oxide. (e) Center of mass of cerium (blue) and oxygen (orange) in a cerium oxide grain boundary in two dimensions.</figcaption>
</p>


The code has been developed to analyse DL_POLY and DL_MONTE calculations however other codes can be incorporated if there is user demand. Other formats, such as pdb or xyz can be converted to `DL_POLY` format with codes such as [atomsk](https://atomsk.univ-lille.fr/) and then analysed with `polypy`. Users are welcome to increase the file coverage by adding a reading function for a different format. This can be accomplished by adding to the `read` module which has a class for each unique file type that converts it to a `polypy.read.trajectory` object. 

`polypy` was developed during a PhD project and as such the functionality focuses on the research questions encountered during that project, which we should clarify are wide ranging. Code contributions aimed at expanding the code to new of problems are encouraged.

`polypy` is free to use.

## Usage

A full list of examples can be found in the examples folder of the git repository, these include both the Python scripts and jupyter notebook tutorials which combine the full theory with code examples. It should be noted however that DL_POLY HISTORY files and DL_MONTE ARCHIVE files are sizable (1-5GB) and as such only short example trajectories are included in this repository. Notebooks are provided here to illustrate the theory but are not practicle.

## Installation

`polypy` is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

    pip install -r requirements.txt

    python setup.py build

    python setup.py install

Or alternatively install with pip

    pip install polypy

Using conda, 

    conda skeleton pypi polypy

    conda build polypy
    
    conda install --use-local polypy


### Tests

Tests can be run by typing:

    python setup.py test


in the root directory. 


### Documentation

To build the documentation from scratch
  
    cd docs
    make html

An online version of the documentation can be found [here](https://polypy.readthedocs.io/en/latest/index.html). The documentation contains an extensive explanation of the underlying theory, function documentation and tutorials. 


### License

`polypy` is made available under the MIT License.

### Detailed requirements

`polypy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- matplotlib
- numpy
- scipy
- coveralls
- coverage
- seaborn
- pandas
- jupyter
- nbsphinx
- jupyter-sphinx==0.2.4
- sphinx_rtd_theme

## Contributing

### Contact

If you have questions regarding any aspect of the software then please get in touch with the developer Adam Symington via email - ars44@bath.ac.uk.
Alternatively you can create an issue on the [Issue Tracker](https://github.com/symmy596/PolyPy/issues).

### Bugs

There may be bugs. If you think you've caught one, please report it on the [Issue Tracker](https://github.com/symmy596/PolyPy/issues).
This is also the place to propose new ideas for features or ask questions about the design of `polypy`. Poor documentation is considered a bug
so feel free to request improvements.

### Code contributions

We welcome help in improving and extending the package. This is managed through Github pull requests; for external contributions we prefer the
["fork and pull"](https://guides.github.com/activities/forking/)
workflow while core developers use branches in the main repository:

   1. First open an Issue to discuss the proposed contribution. This
      discussion might include how the changes fit polypy's scope and a
      general technical approach.
   2. Make your own project fork and implement the changes
      there. Please keep your code style compliant with PEP8.
   3. Open a pull request to merge the changes into the main
      project. A more detailed discussion can take place there before
      the changes are accepted.

For further information please contact Adam Symington, ars44@bath.ac.uk

## Future

Listed below are a series of useful additions that we would like to make to the codebase. Users are encouraged to fork the repository and work on any of these problems. Indeed, if functionality is not listed below you are more than welcome to add it. 

- RDF
- Diagonal slices
- Regional MSDs in a cube


## Author

* Adam R.Symington
  
## Research 

- [Defect segregation facilitates oxygen transport at fluorite UO2 grain boundaries](https://royalsocietypublishing.org/doi/full/10.1098/rsta.2019.0026)
- [The role of dopant segregation on the oxygen vacancy distribution and oxygen diffusion in CeO2 grain boundaries](https://iopscience.iop.org/article/10.1088/2515-7655/ab28b5/meta)
- [Quantifying the impact of disorder on Li-ion and Na-ion transport in perovskite titanate solid electrolytes for solid-state batteries](https://pubs.rsc.org/en/content/articlehtml/2020/ta/d0ta05343k)
- [Elucidating the nature of grain boundary resistance in lithium lanthanum titanate](https://pubs.rsc.org/en/content/articlelanding/2021/TA/D0TA11539H#!divAbstract)

## Acknowledgements
  
This package was written during a PhD project that was funded by AWE and EPSRC (EP/R010366/1). The `polypy` software package was developed to analyse data generated using the Balena HPC facility at the University of Bath and the ARCHER UK National Supercomputing Service (http://www.archer.ac.uk) via our membership of the UK's HEC Ma-terials Chemistry Consortium funded by EPSRC (EP/L000202).The author would like to thank Andrew R. McCluskey, Benjamin Morgan, Marco Molinari, James Grant and Stephen C. Parker for their help and guidance during this PhD project.
---
title: 'polypy - Analysis Tools for Solid State Molecular Dynamics and Monte Carlo Trajectories'
tags:
- Chemistry
- Physics
- Materials Science
- Solid State Chemistry
- Simulation
- Molecular Dynamics
- Monte Carlo
authors:
- name: Adam R. Symington
  orcid: 0000-0001-6059-497X
  affiliation: "1"
affiliations:
- name: Department of Chemistry, University of Bath
  index: 1
date: 29 September 2020
bibliography: paper.bib
---

# Summary

A large number of research questions in solid state chemistry can be addressed using molecular dynamics and Monte Carlo simulations. These simulations allow many material properties to be calculated for direct comparison with experiment. These include the diffusion coefficients, ionic conductivities, charge density, electric field, and electrostatic potential. The diffusion coefficient and ionic conductivity are of particular importance for the study of battery materials (e.g., Li-ion / Na-ion diffusion [@LLTO; @LLTO_2]), solid oxide fuel cell materials (e.g. O-ion diffusion [@CeO2]) and many other applications in solid state chemistry. The charge density, electric field, and electrostatic potential are of interest to problems relating to interfaces in solid state chemistry, e.g., space charge theory .[@Maier_BerBunsenges1984; @Maier_JPhysChemSol1985; @Maier_SolStatIonics2003; @ChiangEtAl_ApplPhysLett1996; @KimAndMaier_JElectrochemSoc2002] Finally, calculating the distribution of defects in a material is useful for the study of segregation behaviour [@UO2; @CeO2] or adsorption behaviour.[@Nora]

In a molecular dynamics simulation the positions of atoms throughout time are being simulated. A molecular dynamics trajectory is a snapshot of the positions occupied by each atom in the simulation as a function time. For example, the trajectory of a single atom would show, sequentially, all of the positions occupied by that atom throughout the simulation. In a Monte Carlo simulation the positions of atoms are updated randomly to provide a statistical ensemble describing the material. A Monte Carlo trajectory is similar although the simulation is not time resolved and the atom positions are simply a function of simulation step, not simulation timestep. The positions of the atoms allow the particle density of each atom to be determined and from these, the electrostatic potential, electric field and charge density can be calculated. A mean squared displacement can be performed on the molecular dynamics trajectories and from these, the diffusion coefficients and ionic conductivities can be calculated. Diffusion coefficients and ionic conductivities can then be used to estimate the activation energy for diffusion using the Arrhenius relationship. 

The `polypy` code is designed to solve the following problems.

- Read DL_POLY [@Smith] and DL_MONTE [@Purton] trajectories.
- Calculate the particle density of all species in a trajectory in one and two dimensions.
- Calculate the charge density in one and two dimensions.
- Calculate the electric field and electrostatic potential in one dimension.
- Calculate the mean squared displacement for a given atom and use this to calculate the diffusion coefficient and ionic conductivity.
- Calculate the volume as a function of simulation timestep.
- Generate publication ready figures.

`polypy` has been used to study Li-ion diffusion in lithium lanthanum titanate, [@LLTO; @LLTO_2] oxygen diffusion and cation migration in both uranium oxide and cerium oxide, [@UO2; @CeO2] thus there is a clear research application. 

![Figure 1 - (a) Mean squared displacement for fluorine diffusion in calcium fluoride. The msd has been plotted in one, two, and three dimensions. (b) The evolution of the system volume during a molecular dynamics simulation. (c) The particle density of cerium (blue) and oxygen (orange) atoms at a grain boundary in cerium oxide. (d) The electrostatic potential at a cerium oxide grain boundary. (e) The center of mass of cerium (blue) and oxygen (orange) atoms in two dimensions, at a grain boundary in cerium oxide.](fig_1.png)

# `polypy`

`polypy` is a Python module for analysing molecular dynamics and Monte Carlo trajectories generated from the DL_POLY [@Smith] and DL_MONTE [@Purton] codes. The code reads DL_POLY HISTORY and CONFIG files, DL_MONTE ARCHIVE files, and stores the data in a `polypy.read.Trajectory` object that is then used by the various data analysis modules.

The `polypy.density.Density` module generates a three dimensional grid and counts the number of times a given atom spends at each grid point during the simulation. This is then used to generate the particle density of a given atom in one and two dimensions. From here, the charge density in one and two dimensions, the electric field in one dimension, and electrostatic potential in one dimension can be calculated using the `polypy.analysis` module. 

The `polypy.msd` module performs a mean squared displacement calculation. From the mean squared displacement, the three, two and one dimensional diffusion coefficient, and ionic coefficient can be calculated. 

A module allowing easy generation of publication plots from the calculated data is available. The outputs are returned in a sensible form, allowing further manipulation and plotting.
The repository contains examples of the core functionality as well as tutorials, implemented in Jupyter notebooks to explain the full theory. Furthermore, a detailed description of theory is also available within the documentation. `polypy` is aimed towards theoretical solid state physicists and chemists who have a basic familiarity with Python, although the examples contained in the repository are designed to help less experienced users make use of the code.

# Statement of Need

There are a number of codes designed to analyse molecular dynamics trajectories that currently exist [@mdanalysis; @chemfiles; @dlanalyser]. DL_ANALYSER [@dlanalyser] is available under license for the analysis of DL_POLY simulations and chemfilles [@chemfiles] is available for the analysis of a wide range of file types. `MDAnalysis` [@mdanalysis] is the most widely used molecular dynamics analysis code and some of the functionality in `polypy` is already present in `MDAnalysis`. The `MDAnalysis.analysis.lineardensity` module calculates the charge density in different dimensions, although according to the documentation, is limited to orthorombic, fixed volume cells. `polypy` is designed to work with several simulation ensembles including, NPT, NVT, semi grand, and grand canonical. Furthermore, the calculation of the electric field and electrostatic potential is unique to `polypy`. `MDAnalysis` and `polypy` are both capable of calculating mean squared displacements. `polypy` goes a step further by allowing the calculation of diffusion coefficients and conductivities within localised regions of a structure, e.g., a grain boundary [@UO2: @CeO2] or local structural environments [@LLTO]. `polypy` is also unique in the sense that it is designed for the analysis of both molecular dynamics and Monte Carlo trajectories.

In summary the features that are unique to `polypy` are as follows

- The analysis of both molecular dynamics and Monte Carlo trajectories. 
- The calculation of the electric field and electrostatic potential.
- Regional mean squared displacements.

# Acknowledgements
  
This package was written during a PhD project that was funded by AWE and EPSRC (EP/R010366/1). The `polypy` software package was developed to analyse data generated using the Balena HPC facility at the University of Bath and the ARCHER UK National Supercomputing Service (http://www.archer.ac.uk) via our membership of the UK's HEC Materials Chemistry Consortium funded by EPSRC (EP/L000202). The author would like to thank Andrew R. McCluskey, Benjamin Morgan, Marco Molinari, James Grant, and Stephen C. Parker for their help and guidance during this PhD project.

# References
Atomic Density
==============

Density Analysis
----------------

Understanding the positions of atoms in a material is incredibly useful
when studying things like atomic structure and defect segregation.
Consider a system with an interface, it may be interesting to know how
the distributions of the materials atoms change at that interface, e.g
is there an increase or decrease in the amount of a certain species at
the interface and does this inform you about any segregation behaviour?

This module of polypy allows the positions of atoms in a simulation to
be evaluated in one and two dimensions, this can then be converted into
a charge density and (in one dimension) the electric field and
electrostatic potential.

.. code:: ipython3

    from polypy.read import History
    from polypy.read import Archive
    from polypy.density import Density
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    
    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    warnings.filterwarnings('ignore')

In this tutorial, we will use polypy to analyse a molecular dynamics
simulation of a grain boundary in fluorite cerium oxide and a Monte
Carlo simulation of Al, Li and Li vacancy swaps in lithium lanthanum
titanate.

Example 1 - Cerium Oxide Grain Boundary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we will use ``polypy`` to analyse a molecular dynamics
simulation of a grain boundary in cerium oxide.

The first step is to read the data. We want the data for both species so
need to provide a list of the species.

::

   ["CE", "O"]

Note. In all examples, an ``xlim`` has been specified to highlight the
grain boundary. Feel free to remove the ``ax.set_xlim(42, 82)`` to see
the whole plot.

.. code:: ipython3

    history = History("../example_data/HISTORY_GB", ["CE", "O"])
    
    print(np.amin(history.trajectory.cartesian_trajectory))
    print(np.amax(history.trajectory.fractional_trajectory))


.. parsed-literal::

    -63.929
    0.9999993486383602


The next step is to create the density object for both species. In this
example we create a seperate object for the cerium and oxygen atoms and
we will be analysing the positions to a resolution of 0.1 angstroms.

.. code:: ipython3

    ce_density = Density(history.trajectory, atom="CE", histogram_size=0.1)
    o_density = Density(history.trajectory, atom="O", histogram_size=0.1)

All subsequent analysis is performed on these two objects.

One Dimension
-------------

Particle Density
~~~~~~~~~~~~~~~~

The ``one_dimensional_density`` function will take a direction which
corresponds to a dimension of the simulation cell. For example, ‘x’
corresponds to the first lattice vector. The code will calculate the
total number of a species in 0.1 angstrom histograms along the first
cell dimension.

The function will return the positions of the histograms and the total
number of species. These can then be plotted with the
``one_dimensional_density_plot`` function which takes a list of
histogram values, a list of particle densities and a list of labels.

.. code:: ipython3

    cx, cy, c_volume = ce_density.one_dimensional_density(direction="z")
    ox, oy, o_volume = o_density.one_dimensional_density(direction="z")
    
    ax = plotting.one_dimensional_density_plot([cx, ox], [cy, oy], ["Ce", "O"])
    ax.set_xlim(42, 82)
    plt.show()




.. image:: Figures/output_9_0.png


Charge Density
~~~~~~~~~~~~~~

The particle densities can be combined with the atom charges to generate
the one dimensional charge density according to

.. math::
    \rho_q(z) = \sum_{i} q_i \rho_i(z)


where :math:`\rho_{i}` is the particle density of atom i and
:math:`q_{i}` is its charge.

The ``OneDimensionalChargeDensity`` class is used for the charge
density, electric field and electrostatic potential. It requires a list
of particle densities, list of charges, the histogram volume and the
total number of timesteps.

.. code:: ipython3

    charge = analysis.OneDimensionalChargeDensity(ox, [oy, cy], [-2.0, 4.0], c_volume, history.trajectory.timesteps)
    
    dx, charge_density = charge.calculate_charge_density()
    
    ax = plotting.one_dimensional_charge_density_plot(dx, charge_density)
    ax.set_xlim(42, 82)
    
    plt.show()



.. image:: Figures/output_11_0.png


Electric Field and Electrostatic Potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The charge density can be converted into the electric field and the
electrostatic potential.

.. math::
    E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz'

.. math::
    \Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz'


where :math:`\rho_{i}` is the charge density and :math:`\epsilon_{0}` is
the permittivity of free space The ``calculate_electric_field`` and
``calculate_electrostatic_potential`` functions will return the electric
field and the electrostatic potential.

.. code:: ipython3

    dx, electric_field = charge.calculate_electric_field()
    
    ax = plotting.electric_field_plot(dx, electric_field)
    ax.set_xlim(42, 82)
    plt.show()



.. image:: Figures/output_13_0.png


.. code:: ipython3

    dx, electrostatic_potential = charge.calculate_electrostatic_potential()
    
    ax = plotting.electrostatic_potential_plot(dx, electrostatic_potential)
    ax.set_xlim(42, 82)
    
    plt.show()



.. image:: Figures/output_14_0.png


Two Dimensions
--------------

Particle Density
~~~~~~~~~~~~~~~~

The particle density can be evaluated in two dimensions. The
``two_dimensional_density`` function will calculate the total number of
species in histograms. The coordinates in x and y of the box are
returned and a grid of species counts are returned.

In this example, the colorbar has been turned off, we are using a grey
palette and the data is being plotted on a log scale.

.. code:: ipython3

    cx_2d, cy_2d, cz_2d, c_volume = ce_density.two_dimensional_density(direction="x")
    ox_2d, oy_2d, oz_2d, o_volume = o_density.two_dimensional_density(direction="x")
    
    fig, ax = plotting.two_dimensional_density_plot(cx_2d, cy_2d, cz_2d, colorbar=False, palette="Greys", log=True)
    ax.set_xlim(42, 82)
    ax.axis('off')
    plt.show()
    
    fig, ax = plotting.two_dimensional_density_plot(ox_2d, oy_2d, oz_2d, colorbar=False, palette="Greys", log=True)
    ax.set_xlim(42, 82)
    ax.axis('off')
    plt.show()




.. image:: Figures/output_16_0.png



.. image:: Figures/output_16_1.png


Charge Density
~~~~~~~~~~~~~~

In the same fashion as the one dimensional case, the charge density can
be evaluated in two dimensions using the
``two_dimensional_charge_density`` function. This function requires the
two dimensional array of atom positions, the atom charges, the volume at
each grid point and the total number of timesteps in the simulation.

.. code:: ipython3

    charge_density = analysis.two_dimensional_charge_density([oz_2d, cz_2d], [-2.0, 4.0], o_volume, history.trajectory.timesteps)
    
    fig, ax = plotting.two_dimensional_charge_density_plot(ox_2d, oy_2d, charge_density, palette='bwr')
    ax.set_xlim(42, 82)
    plt.show()



.. image:: Figures/output_18_0.png


One and Two Dimensions
----------------------

The contour plots can give a good understanding of the average positions
of the atoms (or the location of the lattice sites) however it does not
give a good representation of how many species are actually there. The
``combined_density_plot`` function will evaluate the particle density in
one and two dimensions and then overlay the two on to a single plot,
allowing both the lattice sites, and total density to be viewed.

In this example we are using an orange palette and orange line color for
the cerium atoms, a blue palette and blue line for the oxygen positions
and the data is plotted on a log scale.

.. code:: ipython3

    fig, ax = plotting.combined_density_plot(cx_2d, cy_2d, cz_2d, palette="Oranges", linecolor="orange", log=True)
    for axes in ax:
        axes.set_xlim(42, 82)
    plt.show()
    
    fig, ax = plotting.combined_density_plot(ox_2d, oy_2d, oz_2d, palette="Blues", linecolor="blue", log=True)
    for axes in ax:
        axes.set_xlim(42, 82)
    plt.show()



.. image:: Figures/output_20_0.png



.. image:: Figures/output_20_1.png


Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~

Finally, ``polypy.plotting`` has some functions that will generate a
single contour plot for all species. This function requires the a list
of x axes, a list of y axes, a list of two dimensional arrays
corresponding to the x and y axes and a list of color palettes.

.. code:: ipython3

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], 
                                                                     [cz_2d, oz_2d], ["Blues", "Oranges"], 
                                                                     log=True)
    ax.set_xlim(42, 82)
    plt.show()



.. image:: Figures/output_22_0.png


When analysing things like the electrostatic potential, it is useful to
be able to view how the electrostatic potential changes with structure,
it is very easy to use the ``polypy.plotting`` functions in conjunction
with matplotlib to visualise the relationships.

.. code:: ipython3

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], 
                                                                     [cz_2d, oz_2d], ["Blues", "Oranges"], 
                                                                     log=True)
    ax.set_xlim(42, 82)
    ax2 = ax.twinx()
    ax2.plot(dx, electrostatic_potential, color="green")
    ax2.set_ylabel("Electrostatic Potential (V)")
    plt.show()



.. image:: Figures/output_24_0.png


Finally, ``polypy.plotting`` can generate a contour plot showing the
number density in one and two dimensions in a single plot. This function
requires the a list of x axes, a list of y axes, a list of two
dimensional arrays corresponding to the x and y axes, a list of color
palettes, a list of labels and a list of line colors.

.. code:: ipython3

    fig, ax = plotting.combined_density_plot_multiple_species(x_list=[cx_2d, ox_2d], 
                                                              y_list=[cy_2d, oy_2d], 
                                                              z_list=[cz_2d, oz_2d], 
                                                              palette_list=["Blues", "Oranges"], 
                                                              label_list=['Ce', 'O'], 
                                                              color_list=["blue", "orange"],
                                                              log=True)
    for axes in ax:
        axes.set_xlim(42, 82)
    plt.show()



.. image:: Figures/output_26_0.png


Example 2 - Li, Al and Li vacancy swaps
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this example we will analyse a Monte Carlo simulation of Al doped
lithium lanthanum titanate. It is possible to use molecular dynamics
simulations to study defect segregation if the defects have a relatively
high diffusion coefficient. One could randomly dope a configuration, run
a long molecular dynamics simulation and then analyse the evolution of
the defect locations. When the diffusion coefficient of your defect is
very low, it is not possible to use molecular dynamics simulations to
study defect segregation because you would need a huge MD simulation, in
order to record enough statistics. Monte Carlo simulations allow you to
perform unphysical moves and with a comparitively small Monte Carlo
simulation, you can generate enough statistics to reliably study things
like defect segregation.

In this example, we are analysing a MC simulation of Al in LLZO.
:math:`Al^{3+}` has been doped on the :math:`Li^{+}` sites and charge
compensating Li vacancies have been added. Ultimately, we want to
calculate how the Al doping effects the Li conductivity, however without
a representative distribution of Al/Li/Li vacancies we can’t calculate a
representative conductivity. After 10 ns of MD, the distribution of Al
was unchanged, so Monte Carlo simulations with swap moves are needed to
shake up the distribution. The following swap moves were used;

-  Al <-> Li
-  Al <-> :math:`V_{Li}`
-  Li <-> :math:`V_{Li}`

ARCHIVE_LLZO is a short MC trajectory that we will analyse.

First we will extract and plot the configuration at the first timestep
and then we will plot the positions across the whole simulation to see
how the distributions have changed.

.. code:: ipython3

    archive = Archive("../example_data/ARCHIVE_LLZO", ["LI", "AL", "LV"])
    config_1 = archive.trajectory.get_config(1)

Timestep 1
''''''''''

.. code:: ipython3

    li_density = Density(config_1, atom="LI", histogram_size=0.1)
    al_density = Density(config_1, atom="AL", histogram_size=0.1)
    lv_density = Density(config_1, atom="LV", histogram_size=0.1)

.. code:: ipython3

    lix, liy, li_volume = li_density.one_dimensional_density(direction="y")
    alx, aly, al_volume = al_density.one_dimensional_density(direction="y")
    lvx, lvy, lv_volume = lv_density.one_dimensional_density(direction="y")
    
    ax = plotting.one_dimensional_density_plot([lix, lvx, alx], [liy, lvy, aly], ["Li", "$V_{Li}$", "Al"])
    plt.show()



.. image:: Figures/output_31_0.png


Full Simulation
'''''''''''''''

Disclaimer. This is a short snapshot of a simulation and is not fully
equilibriated, however it provides an example of the ``polypy``
functionailty.

Interestingly, what we find is that the Al, Li and :math:`V_{Li}` tend
to distribute in an even pattern within the structure. This is in sharp
contrast to the distribution at the start of the simulation.

.. code:: ipython3

    li_density = Density(archive.trajectory, atom="LI", histogram_size=0.1)
    al_density = Density(archive.trajectory, atom="AL", histogram_size=0.1)
    lv_density = Density(archive.trajectory, atom="LV", histogram_size=0.1)

.. code:: ipython3

    lix, liy, li_volume = li_density.one_dimensional_density(direction="y")
    alx, aly, al_volume = al_density.one_dimensional_density(direction="y")
    lvx, lvy, lv_volume = lv_density.one_dimensional_density(direction="y")
    
    ax = plotting.one_dimensional_density_plot([lix, lvx, alx], [liy, lvy, aly], ["Li", "$V_{Li}$", "Al"])
    plt.show()



.. image:: Figures/output_34_0.png


.. code:: ipython3

    lix_2d, liy_2d, liz_2d, li_volume = li_density.two_dimensional_density(direction="z")
    alx_2d, aly_2d, alz_2d, al_volume = al_density.two_dimensional_density(direction="z")
    lvx_2d, lvy_2d, lvz_2d, lv_volume = lv_density.two_dimensional_density(direction="z")


.. code:: ipython3

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([alx_2d, lvx_2d], [aly_2d, lvy_2d], 
                                                                     [alz_2d, lvz_2d], ["Blues", "Oranges"], 
                                                                     log=True, figsize=(6, 6))
    plt.show()



.. image:: Figures/output_36_0.png


.. code:: ipython3

    fig, ax = plotting.combined_density_plot_multiple_species(x_list=[lix_2d, alx_2d, lvx_2d],
                                                              y_list=[liy_2d, aly_2d, lvy_2d], 
                                                              z_list=[liz_2d, alz_2d, lvz_2d],
                                                              palette_list=["Greens", "Blues", "Oranges"], 
                                                              label_list=["Li", 'Al', '$V_{Li}$'], 
                                                              color_list=["green", "blue", "orange"],
                                                              log=True, figsize=(6, 6))
    plt.show()



.. image:: Figures/output_37_0.png


polypy\.msd
===========

.. automodule:: polypy.msd
    :members:
    :undoc-members:
    :show-inheritance:Tutorial 1 - Reading data
-------------------------

The HISTORY, ARCHIVE and CONFIG classes expects two things, the filename
of the history file and a list of atoms to read. They will return a
``polypy.read.Trajectory`` object, which stores the the atom labels
(``Trajectory.atom_list``), datatype (``Trajectory.datatype``),
cartesian coordinates (``Trajectory.cartesian_coordinates``), fractiona
coordinates (``Trajectory.fractional_coordinates``), reciprocal lattice
vectors (``Trajectory.reciprocal_lv``), lattice vectors
(``Trajectory.lv``) cell lengths (``Trajectory.cell_lengths``), total
atoms in the file (``Trajectory.atoms_in_history``), timesteps
(``Trajectory.timesteps``), total atoms per timestep
(``Trajectory.total_atoms``).

HISTORY Files
~~~~~~~~~~~~~

.. code:: ipython3

    from polypy import read as rd

.. code:: ipython3

    history = rd.History("../example_data/HISTORY_CaF2", ["CA", "F"])

.. code:: ipython3

    print(history.trajectory.fractional_trajectory)


.. parsed-literal::

    [[0.5170937  0.51658126 0.51643485]
     [0.51658126 0.61669107 0.61654466]
     [0.61669107 0.51658126 0.61691069]
     ...
     [0.46866197 0.25395423 0.58485915]
     [0.37035211 0.58795775 0.45221831]
     [0.36552817 0.48637324 0.17484859]]


.. code:: ipython3

    print(history.trajectory.timesteps)


.. parsed-literal::

    500


.. code:: ipython3

    print(history.trajectory.atoms_in_history)


.. parsed-literal::

    750000


.. code:: ipython3

    print(history.trajectory.total_atoms)


.. parsed-literal::

    1500


It is often necessary to remove the equilibriation timesteps from the
simulation. This can be accomlished with the remove_initial_timesteps
method to remove timesteps at the start of the simulation and the
remove_final_timesteps, to remove timesteps at the end of the
simulation.

.. code:: ipython3

    new_history = history.trajectory.remove_initial_timesteps(10)
    print(new_history.timesteps)


.. parsed-literal::

    490


.. code:: ipython3

    new_history = new_history.remove_final_timesteps(10)
    print(new_history.timesteps)


.. parsed-literal::

    480


It is possible to return the trajectory for a single timestep within the
history file or to return the trajectory for a single atom.

.. code:: ipython3

    config_ca = history.trajectory.get_atom("CA")
    
    print(config_ca.fractional_trajectory)


.. parsed-literal::

    [[0.5170937  0.51658126 0.51643485]
     [0.51658126 0.61669107 0.61654466]
     [0.61669107 0.51658126 0.61691069]
     ...
     [0.31458099 0.41869718 0.41764085]
     [0.42742958 0.32461268 0.42507042]
     [0.42485915 0.42183099 0.31564789]]


.. code:: ipython3

    config_1 = history.trajectory.get_config(1)
    
    print(config_1.fractional_trajectory)


.. parsed-literal::

    [[0.53227339 0.51016082 0.50950292]
     [0.52116228 0.62894737 0.61761696]
     [0.62240497 0.50526316 0.6056652 ]
     ...
     [0.39444444 0.44974415 0.45102339]
     [0.45599415 0.37865497 0.39890351]
     [0.36343202 0.49309211 0.3690424 ]]


CONFIG Files
~~~~~~~~~~~~

.. code:: ipython3

    config = rd.Config("../example_data/CONFIG", ["CA", "F"])

.. code:: ipython3

    print(config.trajectory.fractional_trajectory)


.. parsed-literal::

    [[0.51666667 0.51666667 0.51666667]
     [0.51666667 0.61666667 0.61666667]
     [0.61666667 0.51666667 0.61666667]
     ...
     [0.36666667 0.46666667 0.46666667]
     [0.46666667 0.36666667 0.36666667]
     [0.36666667 0.46666667 0.36666667]]


DLMONTE
~~~~~~~

.. code:: ipython3

    archive = rd.Archive("../example_data/ARCHIVE_Short", ["AL"])

.. code:: ipython3

    print(archive.trajectory.timesteps)


.. parsed-literal::

    1000


Theory
======

`polypy` is a Python module to analyse DLPOLY and DLMONTE trajectory files. Before using this code you will need to generate the relevant data. `polypy` is aimed at the solid state community and there are a wide range of applications. 

Charge Density
--------------

Using the density module you can calculate the number density of atoms

.. math::
    \rho_{q}(z) = \sum_{i} q_{i} \rho_{i}(z)

where :math:`\rho_{i}` is the density of atom i and :math:`q_{i}` is its charge.    

Electric Field and Electrostatic Potential
------------------------------------------

The charge density can be converted into the electric field and the electrostatic potential.

The electric field is calculated according to 

.. math::
    E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz'

where :math:`\epsilon_{0}` is the permittivity of the vacuum and :math:`\rho_{q}` is the charge density.  
The electrostatic potential is calculated according to

.. math::
    \Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz'

Mean Squared Displacement
-------------------------

Molecules in liquds, gases and solids do not stay in the same place and move constantly. Think about a drop of dye in a glass of water, as time passes the dye distributes throughout the water. This process is called diffusion and is common throughout nature and an incredibly relevant property for materials scientists who work on things like batteries.  

Using the dye as an example, the motion of a dye molecule is not simple. As it moves it is jostled by collisions with other molecules, preventing it from moving in a straight path. If the path is examined in close detail, it will be seen to be a good approximation to a random walk. In mathmatics a random walk is a series of steps, each taken in a random direction. This was analysed by Albert Einstein in a study of Brownian motion and he showed that the mean square of the distance travelled by a particle following a random walk is proportional to the time elapsed. 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>,


where :math:`\Big \langle r^2 \big \rangle` is the mean squared distance, t is time, :math:`D_t` is the diffusion rate and C is a constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i. 

What is the mean squared displacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Going back to the example of the dye in water, lets assume for the sake of simplicity that we are in one dimension. Each step can either be forwards or backwards and we cannot predict which. From a given starting position, what distance is our dye molecule likely to travel after 1000 steps? This can be determined simply by adding together the steps, taking into account the fact that steps backwards subtract from the total, while steps forward add to the total. Since both forward and backward steps are equally probable, we come to the surprising conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive numbers to our total which now increases linearly with time. Based upon equation 1 it should now be clear that a plot of :math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line, the gradient of which is equal to 6D. Giving us direct access to the diffusion coefficient of the system. The state of the matter effects the shape of the MSD plot, solids, where little to no diffusion is occuring, has a flat MSD profile. In a liquid however, the particles diffusion randomly and the gradient of the curve is proportional to the diffusion coefficient. 

.. image:: Figures/States_of_Matter.png
    :height: 300px
    :align: center

The following example is for fluorine diffusion in :math:`CaF_2`.

.. image:: Figures/MSD_3.png
    :align: center


Ionic Conductivity
------------------

Usefully, as we have the diffusion coefficient, the number of particles (charge carriers) and the ability to calculate the volume, we can convert this data into the ionic conductivity and then the resistance. 

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} Hr

where :math:`\sigma` is the ionic conductivity, D is the diffusion coefficient,:math:`C_F` is the concentration of charge carriers, which in this case if F ions, :math:`e^2` is the charge of the diffusing species, :math:`k_B` is the Boltzmann constant, T is the temperature and Hr is the Haven ratio.

The resistance can then be calculated according to 

.. math::
    \Omega = \frac{1}{\sigma} 


Arrhenius
---------

It is possible to calculate the diffusion coefficients over a large temperature range and then use the Arrhenius equation to calculate the activation energy for diffusion. Common sense and chemical intuition suggest that the higher the temperature, the faster a given chemical reaction will proceed. Quantitatively this relationship between the rate a reaction proceeds and its temperature is determined by the Arrhenius Equation. At higher temperatures, the probability that two molecules will collide is higher. This higher collision rate results in a higher kinetic energy, which has an effect on the activation energy of the reaction. The activation energy is the amount of energy required to ensure that a reaction happens.  
  
.. math::
    k = A e^{(-Ea / RT)}
  
where k is the rate coefficient, A is a constant, Ea is the activation energy, R is the universal gas constant, and T is the temperature (in kelvin).polypy\.utils
=============

.. automodule:: polypy.utils
    :members:
    :undoc-members:
    :show-inheritance:
    
polypy\.density
===============

.. automodule:: polypy.density
    :members:
    :undoc-members:
    :show-inheritance:Mean Squared Displacement MSD
=============================

Mean Squared Displacement (MSD)
-------------------------------

Molecules in liquds, gases and solids do not stay in the same place and
move constantly. Think about a drop of dye in a glass of water, as time
passes the dye distributes throughout the water. This process is called
diffusion and is common throughout nature and an incredibly relevant
property for materials scientists who work on things like batteries.

Using the dye as an example, the motion of a dye molecule is not simple.
As it moves it is jostled by collisions with other molecules, preventing
it from moving in a straight path. If the path is examined in close
detail, it will be seen to be a good approximation to a random walk. In
mathmatics a random walk is a series of steps, each taken in a random
direction. This was analysed by Albert Einstein in a study of Brownian
motion and he showed that the mean square of the distance travelled by a
particle following a random walk is proportional to the time elapsed.

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>

where :math:`\Big \langle r^2 \big \rangle` is the mean squared
distance, t is time, :math:`D_t` is the diffusion rate and C is a
constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a
function of time, the gradient of the curve obtained is equal to 6 times
the self-diffusion coefficient of particle i. The state of the matter
effects the shape of the MSD plot, solids, where little to no diffusion
is occuring, has a flat MSD profile. In a liquid however, the particles
diffusion randomly and the gradient of the curve is proportional to the
diffusion coefficient.

What is the mean squared displacement
-------------------------------------

Going back to the example of the dye in water, lets assume for the sake
of simplicity that we are in one dimension. Each step can either be
forwards or backwards and we cannot predict which. From a given starting
position, what distance is our dye molecule likely to travel after 1000
steps? This can be determined simply by adding together the steps,
taking into account the fact that steps backwards subtract from the
total, while steps forward add to the total. Since both forward and
backward steps are equally probable, we come to the surprising
conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive
numbers to our total which now increases linearly with time. Based upon
equation 1 it should now be clear that a plot of
:math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line,
the gradient of which is equal to 6D. Giving us direct access to the
diffusion coefficient of the system.

.. code:: ipython3

    from polypy import read as rd
    from polypy.msd import MSD 
    from polypy.msd import RegionalMSD 
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    import numpy as np
    import matplotlib.pyplot as plt

This example will use a short (50,000 steps), pre-prepared trajectory of
bulk :math:`CaF_2`. In reality we probably want a considerably longer
simulation (~10,000,000 steps). Such simulations generate huge files
(5GB) and the analysis would take too long for this tutorial.

The first step is to read the history file to generate the data. The
``HISTORY`` class expects two things, the filename of the history file
and a list of atoms to read. It will return a ``polypy.read.Trajectory``
object, which stores the the atom labels (``Trajectory.atom_labels``),
datatype (``Trajectory.data_type``), cartesian coordinates
(``Trajectory.cartesian_coordinates``), fractiona coordinates
(``Trajectory.fractional_coordinates``), reciprocal lattice vectors
(``Trajectory.reciprocal_lv``), lattice vectors (``Trajectory.lv``) cell
lengths (``Trajectory.cell_lengths``), total atoms in the file
(``Trajectory.atoms_in_history``), timesteps (``Trajectory.timesteps``),
total atoms per timestep (``Trajectory.total_atoms``).

.. code:: ipython3

    history_caf2 = rd.History("../example_data/HISTORY_CaF2", ["F"])

Once the data has been read into the code the MSD calculation can be
performed using the ``MSD`` class. The code will return a
``polypy.MSD.MSDContainer`` object, which contains the MSD information.

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=2)
    
    output = f_msd.msd()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    
    plt.show()



.. image:: Figures/output_6_0.png

MSD calculations require a large number of statistics to be considered representative. A full msd will use every single frame of the trajectory as a starting point and effectively do a seperate msd from each starting point, these are then averaged to give the final result.  An MSD is technically an ensemble average over all sweeps and number of particles. 
The sweeps paramter is used to control the number of frames that are used as starting points in the calculation. For simulations with lots of diffusion events, a smaller number will be sufficient whereas simulations with a small number of diffusion events will require a larger number. 

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=10)
    
    output = f_msd.msd()
    
    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_7_0.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.6078332646337548
    One Dimensional Diffusion Coefficient in X 1.6045620180115865
    One Dimensional Diffusion Coefficient in Y 1.6856414148385679
    One Dimensional Diffusion Coefficient in Z 1.5332963610511103


Note:
An MSD is supposed to be linear only after a ballistic regime and it usually lacks statistics for longer times. Thus the linear fit to extract the slope and thus the diffusion coefficient should be done on a portion of the MSD only.
This can be accomplished using the `exclude_initial` and `exclude_final` parameters

.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient(exclude_initial=50, 
                                                                                    exclude_final=50))


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.5912662736409342
    One Dimensional Diffusion Coefficient in X 1.5862517497696607
    One Dimensional Diffusion Coefficient in Y 1.6753802400942055
    One Dimensional Diffusion Coefficient in Z 1.5121668310589353


Arrhenius
---------

It is then possible to take diffusion coefficients, calculated over a
large temperature range and, using the Arrhenius equation calculate the
activation energy for diffusion. Common sense and chemical intuition
suggest that the higher the temperature, the faster a given chemical
reaction will proceed. Quantitatively this relationship between the rate
a reaction proceeds and its temperature is determined by the Arrhenius
Equation. At higher temperatures, the probability that two molecules
will collide is higher. This higher collision rate results in a higher
kinetic energy, which has an effect on the activation energy of the
reaction. The activation energy is the amount of energy required to
ensure that a reaction happens.

.. math::
    k = A * e^{(-Ea / RT)}

where k is the rate coefficient, A is a constant, Ea is the activation
energy, R is the universal gas constant, and T is the temperature (in
kelvin).

Ionic Conductivity
------------------

Usefully, as we have the diffusion coefficient, the number of particles
(charge carriers) and the ability to calculate the volume, we can
convert this data into the ionic conductivity and then the resistance.

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} 


where :math:`\sigma` is the ionic conductivity, D is the diffusion
coefficient, :math:`C_F` is the concentration of charge carriers, which
in this case if F ions, :math:`e^2` is the charge of the diffusing
species, :math:`k_B` is the Boltzmann constant and T is the temperature.

The resitance can then be calculated according to

.. math::
    \Omega = \frac{1}{\sigma} 

So the first step is to calculate the volume, the system voume module
will do this from the given data.

.. code:: ipython3

    volume, step = analysis.system_volume(history_caf2.trajectory)
    average_volume = np.mean(volume[:50])

The number of charge carriers is just the total number of atoms.

.. code:: ipython3

    sigma = analysis.conductivity(history_caf2.trajectory.total_atoms, 
                            average_volume, 
                            output.xyz_diffusion_coefficient(), 
                            1500, 1)

.. code:: ipython3

    print("Ionic Conductivity :", sigma)


.. parsed-literal::

    Ionic Conductivity : 0.0008752727736501591


.. code:: ipython3

    print("Resistivity :", (1 / sigma)) 


.. parsed-literal::

    Resistivity : 1142.5009781004494


Simulation Length
-----------------

It is important to consider the lenght of your simulation (Number of
steps). The above examples use a short trajectory but it is at a
sufficient temperature that there are enough diffusion events to get a
good MSD plot. The following example is of a very short simulation, you
will hopefully note that the MSD plot is clearly not converged.

.. code:: ipython3

    history_short = rd.History("../example_data/HISTORY_short", atom_list=["F"])

.. code:: ipython3

    f_msd_short = MSD(history_short.trajectory, sweeps=2)
    
    output = f_msd_short.msd()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_19_0.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.58656319093229
    One Dimensional Diffusion Coefficient in X 1.5739020833099904
    One Dimensional Diffusion Coefficient in Y 1.630216356788139
    One Dimensional Diffusion Coefficient in Z 1.5555711326987387


Amusingly, this actually does not seem to have a huge effect on the
diffusion coefficient compared to the longer simulation. However these
trajectories are from a CaF\ :math:`_2` simulation at 1500 K and there
are thus a large number of diffusion events in the short time frame.

State of Matter
---------------

It is possible to identify the phase of matter from the MSD plot.

.. raw:: html

   <center>

 Figure 1. The anticipated MSD form for each state of matter.

.. raw:: html

   </center>

The fluorine diffusion discussed already clearly shows that the fluorine
sub lattice has melted and the diffusion is liquid like. Whereas,
carrying out the same analysis on the calcium sub lattice shows that
while the fluorine sub lattice has melted, the Calcium sub lattice is
still behaving like a solid.

.. code:: ipython3

    f_msd = MSD(history_caf2.trajectory, sweeps=2)
    
    output = f_msd.msd()
    
    ax = plotting.msd_plot(output)
    plt.show()



.. image:: Figures/output_23_0.png


Regional MSD Calculations
-------------------------

Often in solid state chemistry simulations involve defects, both
structural e.g. grain boundaries, dislocations and surface, and chemical
e.g. point defects. It is important to try and isolate the contributions
of these defects to the overall properties. Regarding diffusion, it
could be imagined that a certain region within a structure will have
different properties compared with the stoichiometric bulk, e.g. a grain
boundary vs the grains, or the surface vs the bulk. ``polypy`` has the
capability to isolate trajectories that pass within certain regions of a
structure and thus calculate a diffusion coefficient for those regions.

In this example we will calculate the diffusion coefficient in a box
between -5.0 and 5.0 in the dimension of the first lattice vector.

.. code:: ipython3

    f_msd = RegionalMSD(history_caf2.trajectory, -5, 5, dimension="x")
    output = f_msd.analyse_trajectory()

.. code:: ipython3

    ax = plotting.msd_plot(output)
    
    plt.show()



.. image:: Figures/output_26.png


.. code:: ipython3

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())


.. parsed-literal::

    Three Dimensional Diffusion Coefficient 1.597047044241002
    One Dimensional Diffusion Coefficient in X 1.6120172452124801
    One Dimensional Diffusion Coefficient in Y 1.671268658071343
    One Dimensional Diffusion Coefficient in Z 1.5078552294391845


DLMONTE
^^^^^^^

.. code:: ipython3

    archive = rd.Archive("../example_data/ARCHIVE_LLZO", atom_list=["O"])

.. code:: ipython3

    f_msd = MSD(archive.trajectory, sweeps=2)


::


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    <ipython-input-20-2e636209fda5> in <module>
    ----> 1 f_msd = MSD(archive.trajectory, sweeps=2)
    

    /opt/anaconda3/lib/python3.7/site-packages/polypy-0.7-py3.7.egg/polypy/msd.py in __init__(self, data, sweeps)
        153             raise ValueError("ERROR: MSD can only handle one atom type. Exiting")
        154         if data.data_type == "DL_MONTE ARCHIVE":
    --> 155             raise ValueError("DLMONTE simulations are not time resolved")
        156         self.distances = []
        157         self.msd_information = MSDContainer()


    ValueError: DLMONTE simulations are not time resolved

Tutorials
=========

These tutorials are replicated in jupyter notebook form and contained within examples. All of these examples can be found in `examples/notebooks <https://github.com/symmy596/PolyPy/tree/master/examples/notebooks>`_.

All tutorials use fluorite :math:`CeO_2` grain boundary as an example. Due to the large size of DL_POLY and DL_MONTE trajectory files, the tutorial notebooks contained within the git repository use a very short :math:`CaF_2` trajectory. 

.. toctree::
    :maxdepth: 1

    reading_data.rst
    volume.rst
    density_tutorial.rst
    msd_tutorial.rst
    Mean Squared Displacement MSD
=============================

Molecules in liquds, gases and solids do not stay in the same place and move constantly. Think about a drop of dye in a glass of water, as time passes the dye distributes throughout the water. This process is called diffusion and is common throughout nature and an incredibly relevant property for materials scientists who work on things like batteries.  

Using the dye as an example, the motion of a dye molecule is not simple. As it moves it is jostled by collisions with other molecules, preventing it from moving in a straight path. If the path is examined in close detail, it will be seen to be a good approximation to a random walk. In mathmatics a random walk is a series of steps, each taken in a random direction. This was analysed by Albert Einstein in a study of Brownian motion and he showed that the mean square of the distance travelled by a particle following a random walk is proportional to the time elapsed. 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = 6 D_t + C 

where 

.. math::
    \Big \langle r_{i}^{2} \big \rangle = \frac{1}{3} \Big< | r_{i}(t) - r_{i}(0) |^2 \Big>

where :math:`\Big \langle r^2 \big \rangle` is the mean squared distance, t is time, :math:`D_t` is the diffusion rate and C is a constant. If :math:`\Big \langle r_{i}^{2} \big \rangle` is plotted as a function of time, the gradient of the curve obtained is equal to 6 times the self-diffusion coefficient of particle i. 
The state of the matter effects the shape of the MSD plot, solids, where little to no diffusion is occuring, has a flat MSD profile. In a liquid however, the particles diffusion randomly and the gradient of the curve is proportional to the diffusion coefficient. 

What is the mean squared displacement?
--------------------------------------

Going back to the example of the dye in water, lets assume for the sake of simplicity that we are in one dimension. Each step can either be forwards or backwards and we cannot predict which. From a given starting position, what distance is our dye molecule likely to travel after 1000 steps? This can be determined simply by adding together the steps, taking into account the fact that steps backwards subtract from the total, while steps forward add to the total. Since both forward and backward steps are equally probable, we come to the surprising conclusion that the probable distance travelled sums up to zero.

By adding the square of the distance we will always be adding positive numbers to our total which now increases linearly with time. Based upon equation 1 it should now be clear that a plot of :math:`\Big \langle r_{i}^{2} \big \rangle` vs time with produce a line, the gradient of which is equal to 6D. Giving us direct access to the diffusion coefficient of the system. 

Usage
~~~~~

.. code-block:: python

    from polypy import read as rd
    from polypy.msd import MSD 
    from polypy.msd import RegionalMSD 
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    import numpy as np
    import matplotlib.pyplot as plt

This example will use a short (50,000 steps), pre-prepared trajectory of bulk $CaF_2$. In reality we probably want a considerably longer simulation (~10,000,000 steps). Such simulations generate huge files (5GB) and the analysis would take too long for this tutorial. 

The first step is to read the history file to generate the data. The `HISTORY` class expects two things, the filename of the history file and a list of atoms to read. It will return a `polypy.read.Trajectory` object, which stores the the atom labels (`Trajectory.atom_labels`), datatype (`Trajectory.data_type`), cartesian coordinates (`Trajectory.cartesian_coordinates`), fractiona coordinates (`Trajectory.fractional_coordinates`), reciprocal lattice vectors (`Trajectory.reciprocal_lv`), lattice vectors (`Trajectory.lv`) cell lengths (`Trajectory.cell_lengths`), total atoms in the file (`Trajectory.atoms_in_history`), timesteps (`Trajectory.timesteps`), total atoms per timestep (`Trajectory.total_atoms`).

.. code-block:: python

    history = rd.History("../example_data/HISTORY", ["F"])

Once the data has been read into the code the MSD calculation can be performed.

.. code-block:: python

    f_msd = MSD(history.trajectory, sweeps=2)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_1.png
    :align: center

.. code-block:: python

    f_msd = MSD(history.trajectory, sweeps=10)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)

    plt.show()

.. image:: Figures/MSD_2.png
    :align: center

Using the data the diffusion coefficient can then be calculated from the slopes. 

.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.6078332646337548
| One Dimensional Diffusion Coefficient in X 1.6045620180115865
| One Dimensional Diffusion Coefficient in Y 1.6856414148385679
| One Dimensional Diffusion Coefficient in Z 1.5332963610511103


Arrhenius
~~~~~~~~~

It is then possible to take diffusion coefficients, calculated over a large temperature range and, using the Arrhenius equation calculate the activation energy for diffusion. Common sense and chemical intuition suggest that the higher the temperature, the faster a given chemical reaction will proceed. Quantitatively this relationship between the rate a reaction proceeds and its temperature is determined by the Arrhenius Equation. At higher temperatures, the probability that two molecules will collide is higher. This higher collision rate results in a higher kinetic energy, which has an effect on the activation energy of the reaction. The activation energy is the amount of energy required to ensure that a reaction happens.  
  
.. math::
    k = A * e^{(-Ea / RT)}
  
where k is the rate coefficient, A is a constant, Ea is the activation energy, R is the universal gas constant, and T is the temperature (in kelvin).


Ionic Conductivity
~~~~~~~~~~~~~~~~~~

Usefully, as we have the diffusion coefficient, the number of particles (charge carriers) and the ability to calculate the volume, we can convert this data into the ionic conductivity and then the resistance. 

.. math::
    \sigma = \frac{D C_F e^2}{k_B T} 

where :math:`\sigma` is the ionic conductivity, D is the diffusion coefficient, :math:`C_F` is the concentration of charge carriers, which in this case if F ions, :math:`e^2` is the charge of the diffusing species, :math:`k_B` is the Boltzmann constant and T is the temperature. 

The resitance can then be calculated according to 

.. math::
    \Omega = \frac{1}{\sigma} 

So the first step is to calculate the volume, the system volume module will do this from the given data. 

.. code-block:: python

    volume, step = analysis.system_volume(history.trajectory)
    average_volume = np.mean(volume[:50])

The number of charge carriers is just the total number of atoms.

.. code-block:: python

    sigma = analysis.conductivity(history.trajectory.total_atoms, 
                                  average_volume, 
                                  output.xyz_diffusion_coefficient(), 
                                  1500)
    print("Ionic Conductivity :", sigma)

Ionic Conductivity : 0.0008752727736501591

.. code-block:: python

    print("Resistivity :", (1 / sigma))     

Resistivity : 1142.5009781004494

Simulation Length
~~~~~~~~~~~~~~~~~

It is important to consider the lenght of your simulation (Number of steps). The above examples use a short trajectory but it is at a sufficient temperature that there are enough diffusion events to get a good MSD plot. The following example is of a very short simulation, you will hopefully note that the MSD plot is clearly not converged.

.. code-block:: python

    data_short = rd.History("../example_data/HISTORY_short", atom_list=["F"])
    f_msd_short = MSD(data_short.trajectory, sweeps=2)

    output = f_msd_short.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_3.png
    :align: center

.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.58656319093229
| One Dimensional Diffusion Coefficient in X 1.5739020833099904
| One Dimensional Diffusion Coefficient in Y 1.630216356788139
| One Dimensional Diffusion Coefficient in Z 1.5555711326987387


State of Matter
~~~~~~~~~~~~~~~

It is possible to identify the phase of matter from the MSD plot.

.. image:: Figures/States_of_Matter.png
    :height: 300px
    :align: center

The Fluorine diffusion discussed already clearly shows that the fluorine sub lattice has melted and the diffusion is liquid like. Whereas, carrying out the same analysis on the Calcium sub lattice shows that while the fluorine sub lattice has melted, the Calcium sub lattice is still behaving like a solid. 

.. code-block:: python

    history = rd.History("../example_data/HISTORY", ["CA"])

    f_msd = MSD(history.trajectory, sweeps=2)

    output = f_msd.msd()

    ax = plotting.msd_plot(output)
    plt.show()

.. image:: Figures/MSD_4.png
    :align: center

Regional MSD
~~~~~~~~~~~~

Often in solid state chemistry simulations involve defects, both structural e.g. grain boundaries, dislocations and surface, and chemical e.g. point defects. It is important to try and isolate the contributions of these defects to the overall properties. Regarding diffusion, it could be imagined that a certain region within a structure will have different properties compared with the stoichiometric bulk, e.g. a grain boundary vs the grains, or the surface vs the bulk. `polypy` has the capability to isolate trajectories that pass within certain regions of a structure and thus calculate a diffusion coefficient for those regions. 

.. code-block:: python

    history = rd.History("../example_data/HISTORY", atom_list=["F"])

    f_msd = RegionalMSD(history.trajectory, -5, 5)
    output = f_msd.analyse_trajectory()

    ax = plotting.msd_plot(output)

    plt.show()

.. image:: Figures/MSD_5.png
    :align: center


.. code-block:: python

    print("Three Dimensional Diffusion Coefficient", output.xyz_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in X", output.x_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Y", output.y_diffusion_coefficient())
    print("One Dimensional Diffusion Coefficient in Z", output.z_diffusion_coefficient())

| Three Dimensional Diffusion Coefficient 1.597047044241002
| One Dimensional Diffusion Coefficient in X 1.6120172452124801
| One Dimensional Diffusion Coefficient in Y 1.671268658071343
| One Dimensional Diffusion Coefficient in Z 1.5078552294391845
Atomic Density
==============

Understanding the positions of atoms in a material is incredibly useful when studying things like atomic structure and defect segregation. Consider a system with an interface, it may be interesting to know how the distributions of the materials atoms change at that interface, e.g is there an increase or decrease in the amount of a certain species at the interface and does this inform you about any segregation behaviour? 

This module of polypy allows the positions of atoms in a simulation to be evaluated in one and two dimensions, this can then be converted into a charge density and (in one dimension) the electric field and electrostatic potential.

.. code-block:: python

    from polypy.read import History
    from polypy.read import Archive
    from polypy.density import Density
    from polypy import analysis
    from polypy import utils as ut
    from polypy import plotting
    import numpy as np
    import matplotlib.pyplot as plt

In this tutorial, we will use polypy to analyse a molecular dynamics simulation of a grain boundary in fluorite cerium oxide and a Monte Carlo simulation of Al, Li and Li vacancy swaps in lithium lanthanum titanate. 


Example 1 - Cerium Oxide Grain Boundary
---------------------------------------

In this example we will use :py:attr:`polypy` to analyse a molecular dynamics simulation of a grain boundary in cerium oxide.

The first step is to read the data. We want the data for both species so need to provide a list of the species.  ["CE", "O"]

Note. In all examples, an `xlim` has been specified to highlight the grain boundary. Feel free to remove the `ax.set_xlim(42, 82)` to see the whole plot.

.. code-block:: python

    history = History("HISTORY_GB", ["CE", "O"])

The next step is to create the density object for both species. In this example we create a seperate object for the cerium and oxygen atoms and we will be analysing the positions to a resolution of 0.1 angstroms.

.. code-block:: python

    ca_density = Density(history.trajectory, atom="CE")
    f_density = Density(history.trajectory, atom="O")

The :py:attr:`OneDimensionalChargeDensity` class will take a direction which corresponds to a dimension of the simulation cell. For example, 'x' corresponds to the first lattice vector. The code will calculate the total number of a species in 0.1 A histograms along the first cell dimension.

The function will return the positions of the histograms and the total number of species. To be clear These can then be plotted with the :py:attr:`one_dimensional_plot` function which takes a list of histograms values, a list of y values and a list of labels. 

.. code-block:: python

    cx, cy, c_volume = ce_density.one_dimensional_density(direction="z")
    ox, oy, o_volume = o_density.one_dimensional_density(direction="z")

    ax = plotting.one_dimensional_density_plot([ox, cx], [oy, cy], ["O", "Ce"])
    ax.set_xlim(42, 82)
    plt.show()

.. image:: Figures/Density_GB_1.png
    :align: center

Charge Density
~~~~~~~~~~~~~~

The particle densities can be combined with the atom charges to generate the one dimensional charge density according to 

.. math::
    \rho_q(z) = \sum_{i} q_i \rho_i(z)

where :math:`\rho_{i}` is the density of atom i and :math:`q_{i}` is its charge.  

The :py:attr:`OneDimensionalChargeDensity` class is used for the charge density, electric field and electrostatic potential. It requires a list of particle densities, list of charges, the histogram volume and the total number of timesteps.

.. code-block:: python

    charge = analysis.OneDimensionalChargeDensity(ox, [oy, cy], [-2.0, 4.0], c_volume, history.trajectory.timesteps)
    dx, charge_density = charge.calculate_charge_density()

    ax = plotting.one_dimensional_charge_density_plot(dx, charge_density)
    ax.set_xlim(42, 82)
    plt.show()

.. image:: Figures/Density_GB_2.png
    :align: center

Electric Field and Electrostatic Potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The charge density can be converted into the electric field :math:`E(z)` and the electrostatic potential :math:`\Delta_{\psi}(z)`.

.. math::
    E(z) = \frac{1}{- \epsilon_{0}} \int_{z_{0}}^{z} \rho_{q}(z')dz'

.. math::
    \Delta_{\psi}(z) = \int_{z_{0}}^{z} E(z')dz'

where :math:`\rho_{i}` is the charge density and :math:`\epsilon_{0}` is the permittivity of free space
The :py:attr:`calculate_electric_field` and :py:attr:`calculate_electrostatic_potential` functions will take the bin positions, and the charge density and return the electric field and the electrostatic potential. 

.. code-block:: python

    dx, electric_field = charge.calculate_electric_field()

    ax = plotting.electric_field_plot(dx, electric_field)
    ax.set_xlim(42, 82)
    plt.show()

.. image:: Figures/Density_GB_3.png
    :align: center

.. code-block:: python

    dx, electrostatic_potential = charge.calculate_electrostatic_potential()

    ax = plotting.electrostatic_potential_plot(dx, electrostatic_potential)
    ax.set_xlim(42, 82)
    plt.show()

.. image:: Figures/Density_GB_4.png
    :align: center


Two Dimensions
~~~~~~~~~~~~~~

The particle density can be evaluated in two dimensions. The :py:attr:`two_dimensional_density` function will calculate the total number of species in histograms. The coordinates in x and y of the box are returned and a grid of species counts are returned.

In this example, the colorbar has been turned off, we are using a grey palette and the data is being plotted on a log scale. 

.. code-block:: python

    cx_2d, cy_2d, cz_2d, c_volume = ce_density.two_dimensional_density(direction="x")
    ox_2d, oy_2d, oz_2d, o_volume = o_density.two_dimensional_density(direction="x")

    fig, ax = plotting.two_dimensional_density_plot(cx_2d, cy_2d, cz_2d)
    plt.show()

.. image:: Figures/Density_GB_5.png
    :align: center

.. code-block:: python

    fig, ax = plotting.two_dimensional_density_plot(ox_2d, oy_2d, oz_2d)
    plt.show()

.. image:: Figures/Density_GB_6.png
    :align: center


Charge Density
~~~~~~~~~~~~~~

In the same fashion as the one dimensional case, the charge density can be evaluated in two dimensions using the :py:attr:`two_dimensional_charge_density` function. This function requires the two dimensional array of atom positions, the atom charges, the volume at each grid point and the total number of timesteps in the simulation. 

.. code-block:: python

    charge_density = analysis.two_dimensional_charge_density([oz_2d, cz_2d], [-2.0, 4.0], o_volume)

    fig, ax = plotting.two_dimensional_charge_density_plot(ox_2d, oy_2d, charge_density)
    plt.show()

.. image:: Figures/Density_GB_7.png
    :align: center

One and Two Dimensions
----------------------

The contour plots can give a good understanding of the average positions of the atoms (or the location of the lattice sites) however it does not give a good representation of how many species are actually there. The :py:attr:`combined_density_plot` function will evaluate the particle density in one and two dimensions and then overlay the two on to a single plot, allowing both the lattice sites, and total density to be viewed.

In this example we are using an orange palette and orange line color for the cerium atoms, a blue palette and blue line for the oxygen positions and the data is plotted on a log scale. 

.. code-block:: python

    fig, ax = plotting.combined_density_plot(cx_2d, cy_2d, cz_2d)
    plt.show()

.. image:: Figures/Density_GB_8.png
    :align: center

.. code-block:: python

    fig, ax = plotting.combined_density_plot(fx_2d, fy_2d, fz_2d)
    plt.show()

.. image:: Figures/Density_GB_9.png
    :align: center

All Together
------------

Finally, :py:attr:`polypy.plotting` has some functions that will generate a single contour plot for all species. This function requires the a list of x axes, a list of y axes, a list of two dimensional arrays corresponding to the x and y axes and a list of color palettes. 

.. code-block:: python

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], [cz_2d, oz_2d], ["Blues", "Oranges"], log=True)
    plt.show()

.. image:: Figures/Density_GB_10.png
    :align: center

When analysing things like the electrostatic potential, it is useful to be able to view how the electrostatic potential changes with structure, it is very easy to use the :py:attr:`polypy.plotting` functions
in conjunction with matplotlib to visualise the relationships.

.. code-block:: python

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], [cz_2d, oz_2d], ["Blues", "Oranges"], log=True)
    ax.set_xlim(42, 82)
    ax2 = ax.twinx()
    ax2.plot(dx, electrostatic_potential, color="green")
    ax2.set_ylabel("Electrostatic Potential (V)")
    plt.show()

.. image:: Figures/Density_GB_EP.png
    :align: center

Finally, :py:attr:`polypy.plotting` can generate a contour plot showing the number density in one and two dimensions in a single plot.

.. code-block:: python

    fig, ax = plotting.combined_density_plot_multiple_species([cx_2d, ox_2d], [cy_2d, oy_2d], [cz_2d, oz_2d], ["Blues", "Oranges"], log=True)
    plt.show()

.. image:: Figures/Density_GB_11.png
    :align: center


Example 2 - Li, Al and Li vacancy swaps
---------------------------------------

In this example we will analyse a Monte Carlo simulation of Al doped lithium lanthanum titanate. It is possible to use molecular dynamics simulations to study defect segregation if the defects have a relatively high diffusion coefficient. One could randomly dope a configuration, run a long molecular dynamics simulation and then analyse the evolution of the defect locations. When the diffusion coefficient of your defect is very low, it is not possible to use molecular dynamics simulations to study defect segregation because you would need a huge MD simulation, in order to record enough statistics. Monte Carlo simulations allow you to perform unphysical moves and with a comparitively small Monte Carlo simulation, you can generate enough statistics to reliably study things like defect segregation. 

In this example, we are analysing a MC simulation of Al in LLZO. $Al^{3+}$ has been doped on the $Li^{+}$ sites and charge compensating Li vacancies have been added. Ultimately, we want to calculate how the Al doping effects the Li conductivity, however without a representative distribution of Al/Li/Li vacancies we can't calculate a representative conductivity. After 10 ns of MD, the distribution of Al was unchanged, so Monte Carlo simulations with swap moves are needed to shake up the distribution. The following swap moves were used;

- Al <-> Li
- Al <-> $V_{Li}$
- Li <-> $V_{Li}$

ARCHIVE_LLZO is a short MC trajectory that we will analyse. 

First we will extract and plot the configuration at the first timestep and then we will plot the positions across the whole simulation to see how the distributions have changed. 

.. code-block:: python

    archive = Archive("../example_data/ARCHIVE_LLZO", ["LI", "AL", "LV"])
    config_1 = archive.trajectory.get_config(1)

Timestep 1
~~~~~~~~~~

.. code-block:: python

    li_density = Density(config_1, atom="LI", histogram_size=0.1)
    al_density = Density(config_1, atom="AL", histogram_size=0.1)
    lv_density = Density(config_1, atom="LV", histogram_size=0.1)

    lix, liy, li_volume = li_density.one_dimensional_density(direction="y")
    alx, aly, al_volume = al_density.one_dimensional_density(direction="y")
    lvx, lvy, lv_volume = lv_density.one_dimensional_density(direction="y")

    ax = plotting.one_dimensional_density_plot([lix, lvx, alx], [liy, lvy, aly], ["Li", "$V_{Li}$", "Al"])
    plt.show()

.. image:: Figures/Density_12.png
    :align: center

Full Simulation
~~~~~~~~~~~~~~~

Disclaimer. This is a short snapshot of a simulation and is not fully equilibriated, however it provides an example of the `polypy` functionailty. 

Interestingly, what we find is that the Al, Li and $V_{Li}$ tend to distribute in an even pattern within the structure. This is in sharp contrast to the distribution at the start of the simulation.


.. code-block:: python

    li_density = Density(archive.trajectory, atom="LI", histogram_size=0.1)
    al_density = Density(archive.trajectory, atom="AL", histogram_size=0.1)
    lv_density = Density(archive.trajectory, atom="LV", histogram_size=0.1)

    lix, liy, li_volume = li_density.one_dimensional_density(direction="y")
    alx, aly, al_volume = al_density.one_dimensional_density(direction="y")
    lvx, lvy, lv_volume = lv_density.one_dimensional_density(direction="y")

    ax = plotting.one_dimensional_density_plot([lix, lvx, alx], [liy, lvy, aly], ["Li", "$V_{Li}$", "Al"])
    plt.show()

.. image:: Figures/Density_13.png
    :align: center

.. code-block:: python

    lix_2d, liy_2d, liz_2d, li_volume = li_density.two_dimensional_density(direction="z")
    alx_2d, aly_2d, alz_2d, al_volume = al_density.two_dimensional_density(direction="z")
    lvx_2d, lvy_2d, lvz_2d, lv_volume = lv_density.two_dimensional_density(direction="z")

    fig, ax = plotting.two_dimensional_density_plot_multiple_species([alx_2d, lvx_2d], [aly_2d, lvy_2d], 
                                                                    [alz_2d, lvz_2d], ["Blues", "Oranges"], 
                                                                    log=True, figsize=(6, 6))
    plt.show()

.. image:: Figures/Density_14.png
    :align: center


.. code-block:: python

    fig, ax = plotting.combined_density_plot_multiple_species(x_list=[lix_2d, alx_2d, lvx_2d],
                                                            y_list=[liy_2d, aly_2d, lvy_2d], 
                                                            z_list=[liz_2d, alz_2d, lvz_2d],
                                                            palette_list=["Greens", "Blues", "Oranges"], 
                                                            label_list=["Li", 'Al', '$V_{Li}$'], 
                                                            color_list=["green", "blue", "orange"],
                                                            log=True, figsize=(6, 6))
    plt.show()

.. image:: Figures/Density_15.png
    :align: centerpolypy\.analysis
================

.. automodule:: polypy.analysis
    :members:
    :undoc-members:
    :show-inheritance:
    
Installation
============

:py:mod:`polypy` can be installed from the PyPI package manager with :py:mod:`pip`

.. code-block:: bash 

   pip install polypy

Alternatively, if you would like to download and install the latest development build, it can be found `Github`_ along with specific installation instructions. 

.. _Github: https://github.com/symmy596/polypyVolume
======

.. code-block:: python

    from polypy import analysis
    from polypy import plotting
    from polypy import read as rd
    import matplotlib.pyplot as plt


    history = rd.History("../example_data/HISTORY_CaF2", ["CA"])

    volume, step = analysis.system_volume(history.trajectory)

    ax = plotting.volume_plot(step, volume)
    plt.show()


.. image:: Figures/volume.png
    :align: centerAPI
===

.. toctree::
   :maxdepth: 4

   analysis.rst
   density.rst
   msd.rst
   plotting.rst
   read.rst
   utils.rst
polypy\.read
============

.. automodule:: polypy.read
    :members:
    :undoc-members:
    :show-inheritance:
Using polypy
============

The are number of ways to get and use polypy

- Fork the code: please feel free to fork the code on `Github <https://github.com/symmy596/PolyPy>`_
  and add functionality that interests you.
- Run it locally: polypy is available through the pip package manager.
- Get in touch: Adam R.Symington (ars44@bath.ac.uk) is always keen to chat to potential users.

.. include:: README.rst

.. toctree::
   :maxdepth: 1

   installation.rst
   theory.rst
   using_polypy.rst
   tutorials.rst
   modules



indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
polypy\.plotting
================

.. automodule:: polypy.plotting
    :members:
    :undoc-members:
    :show-inheritance:
    

.. image:: Figures/polypy_1.png
    :align: center

This is the documentation for the open-source Python project - `polypy`.
A library designed to facilitate the analysis of `DL_POLY <https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx>`_ and `DL_MONTE <https://www.ccp5.ac.uk/DL_MONTE>`_ calculations.
`polypy` is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with.
It is hoped that this tool will bring some benfits to the solid state community and facilitate data analysis and the generation of publication ready plots (powered by Matplotlib.)

The main features include:

1. **Method to analyse the number denisty of a given species in one and two dimensions.**  

   - Generate a plot of the total number of species in bins perpendicular to a specified direction.  
   - Generate a plot of the total number of species in cuboids parallel to a specified direction.  

2. **Method calculate the charge density from the number density.**  

   - Convert number densities of all species in bins perpendicular to a specified direction into the charge density.  

3. **Calculate the electric field and electrostatic potential from the charge density.**  

   - Solves the Poisson Boltzmann equation to convert the charge density into the electric field and the electrostatic potential.

4. **Calculate the diffusion coefficient for a given species from a mean squared displacement.**

   - Carries out a mean squared displacement on an MD trajectory.
   - Calculates the diffusion coefficient.
   - Uses the density analysis and the diffusion coefficient to calculate the ionic conductivity. 


.. image:: Figures/Show_off.png
    :align: center

The code has been developed to analyse DL_POLY and DL_MONTE calculations however other codes can be incorporated if there is user demand. 
`polypy` was developed during a PhD project and as such the functionality focuses on the research questions encountered during that project, which we should clarify
are wide ranging. Code contributions aimed at expanding the code to new of problems are encouraged. The code has been developed to analyse DL_POLY and DL_MONTE calculations however other codes can be incorporated if there is user demand. Other formats, such as pdb or xyz can be converted to `DL_POLY` format with codes such as `<atomsk <https://atomsk.univ-lille.fr/>`_. and then analysed with `polypy`. Users are welcome to increase the file coverage by adding a reading function for a different format. This can be accomplished by adding to the `read` module which has a class for each unique file type that converts it to a `polypy.read.trajectory` object. 

`polypy` is free to use.

Usage
-----

A full list of examples can be found in the examples folder of the git repository, these include both the Python scripts and jupyter notebook tutorials which combine the full theory with code examples. It should be noted however that DL_POLY HISTORY files and DL_MONTE ARCHIVE files are sizable (1-5GB) and as such only short example trajectories are included in this repository. Notebooks are provided here to illustrate the theory but are not practicle.

Installation
------------

`polypy` is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

.. code-block:: bash 

    pip install -r requirements.txt

    python setup.py build

    python setup.py install


Or alternatively install with pip

.. code-block:: bash

    pip install polypy

Using conda, 

.. code-block:: bash

    conda skeleton pypi polypy

    conda build polypy
    
    conda install --use-local polypy

Tests
-----

Tests can be run by typing:

    python setup.py test

in the root directory. 

Documentation
-------------

To build the documentation from scratch
  
.. code-block:: bash

    cd docs

    make html

License
-------

`polypy` is made available under the MIT License.

Detailed requirements
---------------------

`polypy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- Numpy
- Scipy
- Matplotlib

Contributing
------------

Contact
~~~~~~~

If you have questions regarding any aspect of the software then please get in touch with the developer Adam Symington via email - ars44@bath.ac.uk.
Alternatively you can create an issue on the `Issue Tracker <https://github.com/symmy596/PolyPy/issues>`_.

Bugs
~~~~

There may be bugs. If you think you've caught one, please report it on the `<Issue Tracker <https://github.com/symmy596/PolyPy/issues>`_.
This is also the place to propose new ideas for features or ask questions about the design of `polypy`. Poor documentation is considered a bug
so feel free to request improvements.

Code contributions
~~~~~~~~~~~~~~~~~~

We welcome help in improving and extending the package. This is managed through Github pull requests; for external contributions we prefer the
`"fork and pull" <https://guides.github.com/activities/forking/>`__
workflow while core developers use branches in the main repository:

   1. First open an Issue to discuss the proposed contribution. This
      discussion might include how the changes fit surfinpy's scope and a
      general technical approach.
   2. Make your own project fork and implement the changes
      there. Please keep your code style compliant with PEP8.
   3. Open a pull request to merge the changes into the main
      project. A more detailed discussion can take place there before
      the changes are accepted.

For further information please contact Adam Symington, ars44@bath.ac.uk

Future
~~~~~~

Listed below are a series of useful additions that we would like to make to the codebase. Users are encouraged to fork the repository and work on any of these problems. Indeed, if functionality is not listed below you are more than welcome to add it. 

- RDF
- Diagonal slices
- Regional MSDs in a cube

Acknowledgements
~~~~~~~~~~~~~~~~
 
This package was written during a PhD project that was funded by AWE and EPSRC (EP/R010366/1). The `polypy` software package was developed to analyse data generated using the Balena HPC facility at the University of Bath and the ARCHER UK National Supercomputing Service (http://www.archer.ac.uk) via our membership of the UK's HEC Ma-terials Chemistry Consortium funded by EPSRC (EP/L000202).The author would like to thank Andrew R. McCluskey, Benjamin Morgan, Marco Molinari, James Grant and Stephen C. Parker for their help and guidance during this PhD project.


API
~~~
