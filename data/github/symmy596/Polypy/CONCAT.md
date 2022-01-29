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
