We aspire to deal with all pull requests in a timely manner. Please be patient.

Ideally pull requests will respect the underlying intrastructure of surfinpy and agree with PEP8, etc. However, please do not let this put novice programmers off contributing as we will do everything we can to help. # SurfinPy


<p align="center"> 
<img src="https://github.com/symmy596/SurfinPy/blob/master/logo/Logo.png?raw=true" width="40%"/>
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2573646.svg)](https://doi.org/10.5281/zenodo.2573646)
[![status](http://joss.theoj.org/papers/368e55451d3fd6ae4b939e1b8e7843ba/status.svg)](http://joss.theoj.org/papers/368e55451d3fd6ae4b939e1b8e7843ba)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04014/status.svg)](https://doi.org/10.21105/joss.04014)
[![Documentation Status](https://readthedocs.org/projects/surfinpy/badge/?version=latest)](https://surfinpy.readthedocs.io/en/latest/)
[![Build Status](https://travis-ci.com/symmy596/SurfinPy.svg?branch=master)](https://travis-ci.com/symmy596/SurfinPy)
[![Build status](https://ci.appveyor.com/api/projects/status/nww04m6tp3335jjr?svg=true)](https://ci.appveyor.com/project/symmy596/surfinpy)
[![PyPI version](https://badge.fury.io/py/surfinpy.svg)](https://badge.fury.io/py/surfinpy)
<a href='https://coveralls.io/github/symmy596/SurfinPy?branch=master'><img src='https://coveralls.io/repos/github/symmy596/SurfinPy/badge.svg?branch=master' alt='Coverage Status' /></a>
<a href='https://gitter.im/Surfinpy/Lobby'>
<img src='https://badges.gitter.im/gitterHQ/gitter.png' alt='Gitter chat' /></a>
 
This is the documentation for the open-source Python project - `surfinpy`.
A library designed to facilitate the generation of publication ready phase diagrams from ab initio calculations for both surface and bulk materials.
surfinpy is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with. 
It is hoped that this tool will bring some benfits to the solid state community and facilitate the generation of publication ready phase diagrams (powered by Matplotlib.)
 
The main features include:

1. **Method to generate surface phase diagrams as a function of chemical potential.**  

   - Generate a diagram as a function of the chemical potential of two adsorbing species e.g. water and carbon dioxide.  
   - Generate a diagram as a function of the chemical potential of one adsorbing species and a surface species e.g. water and oxygen vacancies.  
   - Use experimental data combined with ab initio data to generate a temperature dependent phase diagram.  

2. **Method to generate surface phase diagrams as a function of temperature and pressure.**  

   - Use experimental data combined with ab initio data to generate a pressure vs temperature plot showing the state of a surface as a function of temperature and pressure of one species.

3. **Use calculated surface energies to built crystal morphologies.**  

   - Use the surface energies produced by `surfinpy` alongside Pymatgen to built particle morphologies.  
   - Evaulate how a particles shape changes with temperature and pressure.

4. **Method to generate bulk phase diagrams as a function of chemical potential.**  

   - Generate a diagram as a function of the chemical potential of two species e.g. water and carbon dioxide.  
   - Use experimental data combined with ab initio data to generate a temperature dependent phase diagram.  

5. **Method to generate bulk phase diagrams as a function of temperature and pressure.**  

   - Use experimental data combined with ab initio data to generate a pressure vs temperature plot showing the phase space as a function of temperature and pressure.  

6. **Method to include vibrational properties in a phase diagram**
   
   - Module to calculate the zero point energy and vibrational entropy
   - Encorporate the zero point energy and/or the vibrational entropy into a phase diagram.

The code has been developed to analyse VASP calculations but is compatible with other ab initio codes. 
`surfinpy` was developed across several PhD projects and as such the functionality focuses on the research questions encountered during those projects, which we should clarify 
are wide ranging. Code contributions aimed at expanding the code to new problems are encouraged.

`surfinpy` is free to use.

## Usage

A full list of examples can be found in the examples folder of the git repository, these include jupyter notebook tutorials which combine the full theory with code examples.

## Installation

surfinpy is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

    pip install -r requirements.txt -e .

Or for jupyter compatable use

    pip install -r requirements.txt -e .[Tutorials]

    python setup.py test


Or alternatively install with pip

    pip install surfinpy


### Documentation

To build the documentation from scratch
  
    cd docs
    make html

Alternativly, documentation can be found [here](https://surfinpy.readthedocs.io/en/latest/).

### License

`surfinpy` is made available under the MIT License.

### Detailed requirements

`surfinpy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- Numpy
- Scipy
- Matplotlib
- Pymatgen
- Seaborn
- Pyyaml
- Jupyter (Examples using Jupyter Notebooks, use Tutorials install)


## Contributing

### Contact

If you have questions regarding any aspect of the software then please get in touch with the development team via email Adam Symington (symmy596@gmail.com), Joshua Tse (joshua.s.tse@gmail.com). 
Alternatively you can create an issue on the [Issue Tracker](https://github.com/symmy596/SurfinPy/issues) or you can discuss your questions on our [gitter channel](https://gitter.im/Surfinpy/Lobby).

### Bugs

There may be bugs. If you think you've caught one, please report it on the [Issue Tracker](https://github.com/symmy596/SurfinPy/issues).
This is also the place to propose new ideas for features or ask questions about the design of `surfinpy`. Poor documentation is considered a bug
so feel free to request improvements.

### Code contributions

We welcome help in improving and extending the package. This is managed through Github pull requests; for external contributions we prefer the
["fork and pull"](https://guides.github.com/activities/forking/)__
workflow while core developers use branches in the main repository:

   1. First open an Issue to discuss the proposed contribution. This
      discussion might include how the changes fit surfinpy's scope and a
      general technical approach.
   2. Make your own project fork and implement the changes
      there. Please keep your code style compliant with PEP8.
   3. Open a pull request to merge the changes into the main
      project. A more detailed discussion can take place there before
      the changes are accepted.

For further information please contact Adam Symington - symmy596@gmail.com, Joshua Tse - joshua.s.tse@gmail.com

## Research

- [Strongly Bound Surface Water Affects the Shape Evolution of Cerium Oxide Nanoparticles](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.9b09046)
- [The energetics of carbonated PuO2 surfaces affects nanoparticle morphology: a DFT+U study](https://pubs.rsc.org/lv/content/articlelanding/2020/cp/d0cp00021c/unauth#!divAbstract)
- [Exploiting cationic vacancies for increased energy densities in dual-ion batteries](https://www.sciencedirect.com/science/article/abs/pii/S2405829719310153)
- [Thermodynamic Evolution of Cerium Oxide Nanoparticle Morphology Using Carbon Dioxide](https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.0c07437)

## Author

* Adam R.Symington
* Joshua Tse (Uniersity of Huddersfield)  

## Acknowledgements
  
* [Prof Stephen C.Parker](http://people.bath.ac.uk/chsscp/) - (Bath University)

All notable changes to this project will be documented in the file.

30 September 2020

- Data is no longer stored in dictionaries. A new `surfinpy.data` module has been developed to store the information from ab initio calculations. 
- `surfinpy.bulk_mu_vs_mu` and `surfinpy.bulk_mu_vs_t` have been added to allow phase diagrams of bulk phases to be generated.
- `surfinpy.vibrational_data` has been added to allow the vibrational entropy and zero point energy to be calculated and added to the phase diagrams. 
- `surfinpy.plotting` has been added and `surfinpy.chemical_potential_plot` and `surfinpy.pvt_plot` have been moved here. 
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

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by contacting the project team at ars44@bath.ac.uk. The project team will review and investigate all complaints, and will respond in a way that it deems appropriate to the circumstances. The project team is obligated to maintain confidentiality with regard to the reporter of an incident. Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good faith may face temporary or permanent repercussions as determined by other members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4, available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/## Contributing

If you are interested in contributing to surfinpy. Please feel free; fork the code and go wild to your hearts content. 

If you want to find a good issue to get you in the door as a contributor, check out the issues marked as [good first issue](https://github.com/symmy596/SurfinPy/issues) in the GitHub issue tracker.
---
title: 'SurfinPy 2.0:  A Phase Diagram Generator for Surfaces and Bulk Phases'
tags:
- Chemistry
- Physicsmar
- Density Functional Theory
- Solid State Chemistry
- Simulation
- materials
authors:
- name: Joshua S. Tse
  orcid: 0000-0002-1320-557X
  affiliation: 1
- name: Marco Molinari
  orcid: 0000-0001-7144-6075
  affiliation: 1
- name: Stephen C. Parker
  orcid: 0000-0003-3804-0975
  affiliation: 2
- name: Adam R. Symington
  affiliation: 2
  orcid: 0000-0001-6059-497X
affiliations:
- name: Department of Chemistry, University of Huddersfield
  index: 1
- name: Department of Chemistry, University of Bath
  index: 2
date: 08 December 2020
bibliography: paper.bib
---

# Summary

`SurfinPy` is a Python module for generating phase diagrams from data derived from *ab initio* and/or classical methodologies.

The previous code release, reported by @Symington2019, calculated the surface free energy values under different external conditions, and used these values to generate phase diagrams. 
Surface phase diagrams have been used to provide an understanding of surface composition under various environmental conditions, thus giving crucial information for a range of surface science problems [@Symington2020c;@Symington2020a;@Moxon2020]. 

In this second `SurfinPy` release, the capability of the code has been expanded to generate phase diagrams for bulk phases, as well as surface phases. The code now has the ability to calculate free energy values of bulk phases under specific values of pressure and temperature, and use these to plot phase diagrams for bulk phases as a function  of the chemical potential (and/or pressure) (\autoref{fig:example}) of two species, and as a function  of the chemical potential and temperature. 
Instead of an absolute value, temperature ranges can now be provided, enabling the ability to plot pressure as a function of temperature (or vice-versa), giving results that are comparable to experimental data where available.
Another notable addition to this release is the ability to calculate vibrational properties for bulk phases. 
The vibrational modes for bulk phases are used to calculate the zero point energy and thus the vibrational entropy. 
The code allows for the inclusion of these values into the generation of phase diagrams, which removes the approximation that entropy of bulk phases has little contribution to the free energy, and may improve the accuracy of the methodology.

A significant update to the original code has also been made to improve performance in terms of speedup, streamline workflow and enhanced plotting options.
Finally, eleven tutorials have been developed to highlight the full functionality of this new SurfinPy release. These are all available in Jupyter notebooks in the repository.

![An example of a phase diagram as a function of chemical potential.\label{fig:example}](surfinpy.pdf)



# Statement of Need

`SurfinPy` is a Python module for generating phase diagrams from *ab initio* and/or classical data.

With this release SurfinPy is no longer limited to the surface chemistry and particle morphology, but expands on the chemistry of bulk phases, which makes it ideal for applications to the broad spectrum of research questions in materials science, and solid-state chemistry and physics. 
This allows for the exploration of the phase stability of solid-state systems (bulk and surface phases) of different compositions as a function of external conditions. 

Other codes capable to produce phase diagrams, e.g., pymatgen [@ONG2013314] and ASE [@Hjorth_Larsen_2017] are available.  However, our code is self-contained allowing for the generation of both bulk and surface phase diagrams, while offering easier and enhanced plotting capability to compare phases as a function of chemical potential of different species and temperatures.  Additionally, unlike other codes detailed tutorials are available, offering a more tailored and focused experience compared to other codes.

# Acknowledgements
  
The authors acknowledge support from the EPSRC (EP/K025597/1 and EP/R010366/1), and the Royal Society (Newton Advanced Fellowship NA150190).

# References
# Tutorials and Examples
  
#### Author - Adam Symington

##### Contained here are 11 tutorials for using surfinpy in your research. 

It should be noted that these are generic examples to demonstrate the functionality of the code. 

- Notebooks/Surfaces/Tutorial 1 - Guide to generating a phase diagram for a surface as a function of the chemical potetnial of two variable species.
- Notebooks/Surfaces/Tutorial 2 - Guide to convert the chemcial potential values used in tutorial 1 into meaningful values.
- Notebooks/Surfaces/Tutorial 3 - Guide to convert chemical potential to pressure.
- Notebooks/Surfaces/Tutroial 4 - Guide for creating a phase diagram for a surface as a function of temperature and pressure. 
- Notebooks/Surfaces/Tutorial 5 - Guide to generate particle morphologies - requires an installed version of pymatgen. 

- Notebooks/Bulk/Tutorial 1 - Guide to generating a phase diagram for a bulk material as a function of the chemical potential of two variable species e.g. Bulk + $CO_2$ + $H_2O$.
- Notebooks/Bulk/Tutorial 2 - Guide to including temperature in a bulk phase diagram.
- Notebooks/Bulk/Tutorial 3 - Guide to convert chemical potential to pressure.
- Notebooks/Bulk/Tutroial 4 - Guide for creating a phase diagram for a bulk material as a function of temperature and pressure. 
- Notebooks/Bulk/Tutorial 5 - Guide to calculating and including vibrational entropy in a phase diagram.
- Notebooks/Bulk/Tutorial 6 - Guide to calculating and including vibrational entropy in a temperature dependent phase diagram.

In the Scripts/ you will find the tutorials but in python script format. ---
title: 'surfinpy: A Surface Phase Diagram Generator'
tags:
- Chemistry
- Physics
- Density Functional Theory
- Solid State Chemistry
- Simulation
- materials
authors:
- name: Adam R. Symington
  orcid: 0000-0001-6059-497X
  affiliation: "1"
- name: Joshua Tse
  orcid: 0000-0002-1320-557X
  affiliation: 2
- name: Marco Molinari
  orcid: 0000-0001-7144-6075
  affiliation: 2
- name: Arnaud Marmier
  orcid: 0000-0003-3836-0004
  affiliation: 3
- name: Stephen C. Parker
  orcid: 0000-0003-3804-0975
  affiliation: 1
affiliations:
- name: Department of Chemistry, University of Bath
  index: 1
- name: Department of Chemistry, University of Huddersfield
  index: 2
- name: FET - Engineering, Design and Mathematics, University of the West of England
  index: 3
date: 25 January 2019
bibliography: paper.bib
---

# Summary

A surface phase diagram is a graphical representation of the different physical states of a surface under different conditions.
The surface represents the first point of contact between the material and the environment.
Thus, understanding the state of surface is crucial for a wide range of problems in materials science concerning the relationship between
the state of the surface and the surrounding environmental condtions.
Examples include particle morphologies in solid state batteries [@Canepa2018];
determining the concentration of adsorbed water at a surface depending on synthesis conditions [@Molinari2012] [@Tegner2017];
catalytic reactions [@Reuter2003]; or determing the effect of dopants and impurities on the surface stability.  

Computational modelling can be used to generate surface phase diagrams from energy minimisation data.
One common research question is how water adsorption affects the surface and material properties.
The conventional starting point is to perform a series of energy minimisation calculations with varying concentrations of water on several different slabs.
From the energies, the surface free energy of each calculation (phase) as a function of temperature and pressure can be calculated using a well-established approach [@Molinari2012].
Once the free energy is known under different constants, the phase which is most stable at a specific temperature and pressure, and thus a phase diagram, can be generated.

A further degree of complexity can be introduced by considering surface defects, e.g., vacancies or interstitials, or other adsorbants, e.g., carbon dioxide.
Using surface defects as an example, it is important to consider the relationship between the defective surface, the stoichiometric surface and the adsorbing water molecules.
A surface phase diagram can be constructed as a function of the chemcial potential of the adsorbing species (water) and the surface defect
(e.g., oxygen, if oxygen vacancies are being considered). This is done using the method of Marmier & Parker[-@Marmier2004].

![An example phase diagram as a function of chemical potential (a), and as a function of temperature and pressure (b).\label{fig:example}](Figure_1.png)

# `surfinpy`

`surfinpy` is a Python module for generating surface phase diagrams from DFT data.
It contains two core modules for generating surface phase diagrams using both the methods employed in @Molinari2012 and @Marmier2004.
These allow fast generation of temperature vs. pressure phase diagrams and phase diagrams as a function of chemcial potential of species A and B.
The plotting classes take the outputs of the calculation modules and generate phase diagrams using `matplotlib`.
`surfinpy` is aimed towards theoretical solid state physicist who have a basic familiarity with Python.
The repository contains examples of the core functionality as well as tutorials, implemented in Jupyter notebooks to explain the full theory.
Furthermore, a detailed description of theory is also available within the documentation.

# Acknowledgements
  
ARS would like to thank Andrew R. McCluskey for his guidance through this project. This package was written during a PhD funded by AWE and EPSRC (EP/R010366/1). The input
data for the development and testing of this project was generated using ARCHER UK National Supercomputing Service (http://www.archer.ac.uk) via our membership of
the UK's HEC Ma-terials Chemistry Consortium funded by EPSRC (EP/L000202).

# ReferencesPressure vs Temperature
=======================

`Surfinpy` has the functionality to generate phase diagrams as a function of pressure vs temperature based upon the methodology used in `Molinari et al
<https://pubs.acs.org/doi/abs/10.1021/jp300576b>`_ according to

.. math::
    \gamma_{adsorbed, T, P} = \gamma_{bare} + ( C ( E_{ads, T} - RTln(\frac{p}{p^o})

where :math:`\gamma_{adsorbed, T, p}` is the surface energy of the surface with adsorbed species at temperature (T) and pressure (P),
:math:`\gamma_{bare}` is the surface energy of the bare surface, C is the coverage of adsorbed species, :math:`E_{ads}` is the adsorption energy,

.. math::
    E_{ads, T} =  E_{slab, adsorbant} - (E_{slab, bare} + n_{H_2O} E_{H_2O, T}) / n_{H_2O}

where :math:`E_{slab, adsorbant}` is the energy of the surface and the adsorbed species, :math:`n_{H_2O}` is he number of adsorbed species,

.. math::
    E_{H_2O, (T)} = E_{H_2O, (g)} - TS_{(T)}

where :math:`S_{(T)}` is the experimental entropy of gaseous water in the standard state.

Usage
~~~~~

.. code-block:: python

    from surfinpy import utils as ut
    from surfinpy import p_vs_t

    adsorbant = -14.00
    SE = 1.40

    stoich = {'Cation': 24, 'X': 48, 'Y': 0, 'Area': 60.22,
              'Energy': -575.00, 'Label': 'Bare'}
    H2O =    {'Cation': 24, 'X': 48, 'Y': 2, 'Area': 60.22,
              'Energy': -605.00, 'Label': '1 Water'}
    H2O_2 =  {'Cation': 24, 'X': 48, 'Y': 8, 'Area': 60.22,
              'Energy': -695.00, 'Label': '2 Water'}
    data = [H2O, H2O_2]

    coverage = ut.calculate_coverage(data)

    thermochem = ut.read_nist("H2O.txt")

    system = p_vs_t.calculate(stoich, data, SE,
                              adsorbant,
                              thermochem,
                              coverage)
    system.plot()

.. image:: Figures/Tutorial_2/First.png
    :height: 300px
    :align: center


Alternatively you can also tweak the style

.. code-block:: python

    system.plot(output="dark_pvt.png",
                set_style="dark_background",
                colourmap="PiYG")

.. image:: Figures/Tutorial_2/Second.png
    :height: 300px
    :align: center
Gallery
=======

The gallery is a preview of some of the plots available in ``surfinpy``. Clicking on a plot will provide a link to a tutorial 
for generating the plot. 

Surfaces
--------

Chemical Potential
~~~~~~~~~~~~~~~~~~

The following are examples of a phase diagram as a function of chemical potential. The first is the default output 
and the rest are generated by playing with the style and colourmap.

.. image:: Figures/Surfaces_1.png
    :height: 300px
    :align: center
    :target: tutorial_1.html

.. image:: Figures/Surfaces_2.png
    :height: 300px
    :align: center
    :target: tutorial_1.html


Pressure
~~~~~~~~

Chemical potential can be converted to pressure and a diagram with pressure of species A/B displayed.

.. image:: Figures/Surfaces_5.png
    :height: 300px
    :align: center
    :target: tutorial_1.html#Pressure



Chemical potential and pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`surfinpy` can produce a plot with the chemical potential of A/B on axes X/Y and the pressure of 
A/B on axes X2/Y2. 

.. image:: Figures/Surfaces_4.png
    :height: 300px
    :align: center
    :target: tutorial_1.html#Pressure.. 
    
.. image:: Figures/Surfaces_6.png
    :height: 300px
    :align: center
    :target: tutorial_1.html#Pressure


Temperature vs Pressure
~~~~~~~~~~~~~~~~~~~~~~~

`surfinpy` can produce simple pvt plots showing the relationship between a single species "A" at the surface e.g. water.

.. image:: Figures/Surfaces_7.png
    :height: 300px
    :align: center
    :target: tutorial_2.html


Particle Morphology
~~~~~~~~~~~~~~~~~~~

`surfinpy` provides examples of how to use the surface energy calculation alongside pymatgen to generate particle morphologies 
at different temperature and pressure values.

.. image:: Figures/Tutorial_3/Wulff.png
    :height: 300px
    :align: center
    :target: tutorial_3.html


Bulk
----

Chemical Potential
~~~~~~~~~~~~~~~~~~

The following are examples of a phase diagram as a function of chemical potential.

.. image:: Figures/Bulk_1.png
    :height: 300px
    :align: center
    :target: tutorial_4.html

Temperature
~~~~~~~~~~~

The following are examples of a phase diagram as a function of chemical potential with a temperature contribution introduced.

.. image:: Figures/Bulk_2.png
    :height: 300px
    :align: center
    :target: tutorial_4.html

Pressure
~~~~~~~~

The following are examples of a phase diagram as a function of pressure.


.. image:: Figures/Bulk_4.png
    :height: 300px
    :align: center
    :target: tutorial_4.html


Pressure vs Temperature
~~~~~~~~~~~~~~~~~~~~~~~

The following are examples of a phase diagram as a function of chemical potential, chemical potential and pressure, and temperature.

.. image:: Figures/Bulk_6.png
    :height: 300px
    :align: center
    :target: tutorial_5.html

.. image:: Figures/Bulk_7.png
    :height: 300px
    :align: center
    :target: tutorial_5.html


Vibrational Entropy
~~~~~~~~~~~~~~~~~~~

The following are examples of how to include the effects of vibrational entropy to the phase diagrams.

.. image:: Figures/Bulk_8.png
    :height: 300px
    :align: center
    :target: tutorial_6.html

.. image:: Figures/Bulk_9.png
    :height: 300px
    :align: center
    :target: tutorial_6.html


.. image:: Figures/Bulk_10.png
    :height: 300px
    :align: center
    :target: tutorial_6.html
surfinpy\.utils
===============

The utils module contains functions that are common and find various uses throughout the code. 

.. automodule:: surfinpy.utils
    :members:
    :undoc-members:
    :show-inheritance:
    
Chemical Potential
==================

The physical quantity that is used to define the stability of a surface with a given composition is its surface energy :math:`\gamma` (J :math:`m^{-2}`).
At its core, `surfinpy` is a code that calculates the surface energy of different slabs and uses these surface energies to build a phase diagram.
In this explanation of theory we will use the example of water adsorbing onto a surface of :math:`TiO_2` containing oxygen vacancies.
In such an example there are two variables, water concentration and oxygen vacancy concentration. We are able to calculate the surface energy according to

.. math::
    \gamma_{Surf} = \frac{1}{2A} \Bigg( E_{TiO_2}^{slab} - \frac{nTi_{slab}}{nTi_{Bulk}} E_{TiO_2}^{Bulk} \Bigg) - \Gamma_O \mu_O - \Gamma_{H_2O} \mu_{H_2O} ,

where A is the surface area, :math:`E_{TiO_2}^{slab}` is the DFT energy of the slab, :math:`nTi_{Slab}` is the number of cations in the slab,
:math:`nTi_{Bulk}` is the number of cations in the bulk unit cell, :math:`E_{TiO_2}^{Bulk}` is the DFT energy of the bulk unit cell and

.. math::
    \Gamma_O = \frac{1}{2A} \Bigg( nO_{Slab} - \frac{nO_{Bulk}}{nTi_{Bulk}}nTi_{Slab}  \Bigg) ,

.. math::
    \Gamma_{H_2O} = \frac{nH_2O}{2A} ,

where :math:`nO_{Slab}` is the number of anions in the slab, :math:`nO_{Bulk}` is the number of anions in the bulk and :math:`nH_2O` is the number of adsorbing water molecules.
:math:`\Gamma_O` / :math:`\Gamma_{H_2O}` is the excess oxygen / water at the surface and :math:`\mu_O` / :math:`\mu_{H_2O}` is the oxygen / water chemical potential.
Clearly :math:`\Gamma` and :math:`mu` will only matter when the surface is non stoichiometric.

Usage
~~~~~

The first thing to do is input the data that we have generated from our DFT calculations.
The input data needs to be contained within a dictionary.
First we have created the dictionary for the bulk data, where ``Cation`` is the number of cations, ``Anion`` is the number of anions,
``Energy`` is the DFT energy and ``F-Units`` is the number of formula units.

.. code-block:: python

    bulk = {'Cation' : Cations in Bulk Unit Cell,
            'Anion' : Anions in Bulk Unit Cell,
            'Energy' :  Energy of Bulk Calculation,
            'F-Units' : Formula units in Bulk Calculation}

Next we create the slab dictionaries - one for each slab calculation or "phase". ``Cation`` is the number of cations,
``X`` is in this case the number of oxygen species (corresponding to the X axis of the phase diagram),
``Y`` is the number of in this case water molecules (corresponding to the Y axis of our phase diagram),
``Area`` is the surface area, ``Energy`` is the DFT energy, ``Label`` is the label for the surface (appears on the phase diagram) and
finally ``nSpecies`` is the number of adsorbing species (In this case we have a surface with oxygen vacancies and adsorbing water molecules -
so nSpecies is 1 as oxygen vacancies are not an adsorbing species, they are a constituent part of the surface).

.. code-block:: python

    surface = {'Cation': Cations in Slab,
               'X': Number of Species X in Slab,
               'Y': Number of Species Y in Slab,
               'Area': Surface area in the slab,
               'Energy': Energy of Slab,
               'Label': Label for phase,
               'nSpecies': How many species are non stoichiometric}


This data needs to be contained within a list. Don't worry about the order, `surfinpy` will sort that out for you.

We also need to declare the range in chemical potential that we want to consider.
Again these exist in a dictionary. ``Range`` corresponds to the range of chemical potential values to be considered and ``Label`` is the axis label.

.. code-block:: python

    deltaX = {'Range': Range of Chemical Potential,
              'Label': Species Label}


.. code-block:: python

    from surfinpy import mu_vs_mu

    bulk = {'Cation' : 1, 'Anion' : 2, 'Energy' : -780.0, 'F-Units' : 4}

    pure =     {'Cation': 24, 'X': 48, 'Y': 0, 'Area': 60.0,
                'Energy': -575.0,   'Label': 'Stoich',  'nSpecies': 1}
    H2O =      {'Cation': 24, 'X': 48, 'Y': 2, 'Area': 60.0,
                'Energy': -612.0,   'Label': '1 Water', 'nSpecies': 1}
    H2O_2 =    {'Cation': 24, 'X': 48, 'Y': 4, 'Area': 60.0,
                'Energy': -640.0,   'Label': '2 Water', 'nSpecies': 1}
    H2O_3 =    {'Cation': 24, 'X': 48, 'Y': 8, 'Area': 60.0,
                'Energy': -676.0,   'Label': '3 Water', 'nSpecies': 1}
    Vo =       {'Cation': 24, 'X': 46, 'Y': 0, 'Area': 60.0,
                'Energy': -558.0,   'Label': 'Vo', 'nSpecies': 1}
    H2O_Vo =   {'Cation': 24, 'X': 46, 'Y': 2, 'Area': 60.0,
                'Energy': -594.0,  'Label': 'Vo + 1 Water', 'nSpecies': 1}
    H2O_Vo_2 = {'Cation': 24, 'X': 46, 'Y': 4, 'Area': 60.0,
                'Energy': -624.0,  'Label': 'Vo + 2 Water', 'nSpecies': 1}
    H2O_Vo_3 = {'Cation': 24, 'X': 46, 'Y': 6, 'Area': 60.0,
                'Energy': -640.0, 'Label': 'Vo + 3 Water', 'nSpecies': 1}
    H2O_Vo_4 = {'Cation': 24, 'X': 46, 'Y': 8, 'Area': 60.0,
                'Energy': -670.0, 'Label': 'Vo + 4 Water', 'nSpecies': 1}

    data = [pure, H2O_2, H2O_Vo, H2O,  H2O_Vo_2, H2O_3, H2O_Vo_3,  H2O_Vo_4, Vo]

    deltaX = {'Range': [ -12, -6],  'Label': 'O'}
    deltaY = {'Range': [ -19, -12], 'Label': 'H_2O'}

This data will be used in all subsequent examples and will not be declared again. Once the data has been declared it is a simple
two line process to generate the diagram.

.. code-block:: python

    system = mu_vs_mu.calculate(data, bulk, deltaX, deltaY)
    system.plot_phase()

.. image:: Figures/Tutorial_1/First.png
    :height: 300px
    :align: center

Temperature
~~~~~~~~~~~

The previous phase diagram is at 0K. It is possible to use experimental data from the NIST_JANAF database to make the chemical potential a temperature dependent
term and thus generate a phase diagram at a temperature (T). Using oxygen as an example, this is done according to

.. math::
    \gamma_{Surf} = \frac{1}{2A} \Bigg( E_{TiO_2}^{slab} - \frac{nTi_{Slab}}{nTi_{Bulk}} E_{TiO_2}^{Bulk} \Bigg) - \Gamma_O \mu_O - \Gamma_{H_2O} \mu_{H_2O} - n_O \mu_O (T) - n_{H_2O} \mu_{H_2O} (T)

where

.. math::
    \mu_O (T)  = \frac{1}{2} \mu_O (T) (0 K , DFT) +  \frac{1}{2} \mu_O (T) (0 K , EXP) +  \frac{1}{2} \Delta G_{O_2} ( \Delta T, Exp),

:math:`\mu_O` (T) (0 K , DFT) is the 0K free energy of an isolated oxygen molecule evaluated with DFT, :math:`\mu_O` (T) (0 K , EXP) is the 0 K experimental
Gibbs energy for oxygen gas and $\Delta$ :math:`G_{O_2}` ( :math:`\Delta` T, Exp) is the Gibbs energy defined at temperature T as

.. math::
    \Delta G_{O_2} ( \Delta T, Exp)  = \frac{1}{2} [H(T, {O_2}) -  H(0 K, {O_2})] -  \frac{1}{2} T[S(T, {O_2}])

`surfinpy` has a built in function to read a NIST_JANAF table and calculate this temperature_correction for you. In the following example you will also
see an example of how you can tweak the style and colourmap of the plot.

.. code-block:: python

    from surfinpy import mu_vs_mu

    Oxygen_exp = mu_vs_mu.temperature_correction("O2.txt", 298)
    Water_exp = mu_vs_mu.temperature_correction("H2O.txt", 298)

    Oxygen_corrected = (-9.08 + -0.86 + Oxygen_exp)
    Water_corrected = -14.84 + 0.55 + Water_exp

    system =  mu_vs_mu.calculate(data, bulk, deltaX, deltaY,
                                 x_energy=Oxygen_corrected,
                                 y_energy=Water_corrected)
    system.plot_phase(temperature=298, set_style="fast",
                      colourmap="RdBu")

.. image:: Figures/Tutorial_1/Second.png
    :height: 300px
    :align: center


Pressure
~~~~~~~~

The chemical potential can be converted to pressure values according to

.. math::
    P = \frac{\mu_O}{k_B T}

where P is the pressure, :math:`\mu` is the chemical potential of oxygen, :math:`k_B` is the Boltzmann constant and T is the temperature.


.. code-block:: python

    from surfinpy import mu_vs_mu

    Oxygen_exp = mu_vs_mu.temperature_correction("O2.txt", 298)
    Water_exp = mu_vs_mu.temperature_correction("H2O.txt", 298)

    Oxygen_corrected = (-9.08 + -0.86 + Oxygen_exp)
    Water_corrected = -14.84 + 0.55 + Water_exp

    system =  mu_vs_mu.calculate(data, bulk, deltaX, deltaY,
                                 x_energy=Oxygen_corrected,
                                 y_energy=Water_corrected)
    system.plot_mu_p(output="Example_ggrd", colourmap="RdYlGn",
                     temperature=298)

.. image:: Figures/Tutorial_1/Third.png
    :height: 300px
    :align: center

.. code-block:: python

    system.plot_mu_p(output="Example_ggrd",
                     set_style="dark_background",
                     colourmap="RdYlGn",
                     temperature=298)

.. image:: Figures/Tutorial_1/Fourth.png
    :height: 300px
    :align: center

.. code-block:: python

    system.plot_pressure(output="Example_dark_rdgn",
                         set_style="dark_background",
                         colourmap="PuBu",
                         temperature=298)

.. image:: Figures/Tutorial_1/Filth.png
    :height: 300px
    :align: center
API
===

.. toctree::
   :maxdepth: 4

   data
   mu_vs_mu
   bulk_mu_vs_mu
   bulk_mu_vs_t
   vibrational_data
   p_vs_t
   plotting
   wulff
   utils
surfinpy\.bulk_mu_vs_mu
=======================

.. automodule:: surfinpy.bulk_mu_vs_mu
    :members:
    :undoc-members:
    :show-inheritance:surfinpy\.plotting
==================

.. automodule:: surfinpy.plotting
    :members:
    :undoc-members:
    :show-inheritance:Particle Morphology
===================

It is sometimes useful to use surface energies in order to generate particle morphologies.
This tutorial demonstrates how to obtain surface energies for surfaces containing adsorbed species using `surfinpy`.
With these you can then generate a wulff construction using `pymatgen <https://www.sciencedirect.com/science/article/pii/S0927025612006295?via%3Dihub>`_.
A Wulff construction is a method to determine the equilibrium shape of a crystal.
So by calculating the surface energies of multiple different surfaces, at different temperature and pressure values we can generate a particle morphology for the material,
in the presence of an adsorbing species, at a specific temperature and pressure.

`surfinpy` has a module called wulff that will return a surface energy at a given temperature and pressure value.
These can then be used in conjunction with Pymatgen for a wulff construction.
So first we need to declare the data for each surface and calculate the surface energies.
As an aside, it is possible to provide multiple coverages, the return will be an array of surface energies,
corresponding to each surface coverage, you would then select the minimum value with `np.amin()`

.. code-block:: python

    import numpy as np
    from surfinpy import p_vs_t as pt
    from surfinpy import wulff
    from surfinpy import utils as ut
    from pymatgen.core.surface import SlabGenerator,
                                      generate_all_slabs,
                                      Structure, Lattice
    from pymatgen.analysis.wulff import WulffShape

    adsorbant = -14.22
    thermochem = ut.read_nist('H2O.txt')

The first thing to do is calculate the surface energy at a temperature and pressure value for each surface.

.. code-block:: python

    SE = 1.44
    stoich =      {'M': 24, 'X': 48, 'Y': 0, 'Area': 60.22,
                  'Energy': -575.66, 'Label': 'Stoich'}
    Adsorbant_1 = {'M': 24, 'X': 48, 'Y': 2, 'Area': 60.22,
                   'Energy': -609.23, 'Label': '1 Species'}
    data = [Adsorbant_1]
    Surface_100_1 = wulff.calculate_surface_energy(stoich,
                                                   data,
                                                   SE,
                                                   adsorbant,
                                                   thermochem,
                                                   298,
                                                   0)

    SE = 1.06
    stoich =      {'M': 24, 'X': 48, 'Y': 0, 'Area': 85.12,
                   'Energy': -672.95, 'Label': 'Stoich'}
    Adsorbant_1 = {'M': 24, 'X': 48, 'Y': 2, 'Area': 85.12,
                   'Energy': -705.0, 'Label': '1 Species'}
    data = [Adsorbant_1]
    Surface_110_1 = wulff.calculate_surface_energy(stoich,
                                                   data,
                                                   SE,
                                                   adsorbant,
                                                   thermochem,
                                                   298,
                                                   0)

    SE = 0.76
    stoich =      {'M': 24, 'X': 48, 'Y': 0, 'Area': 77.14,
                   'Energy': -579.61, 'Label': 'Stoich'}
    Adsorbant_1 = {'M': 24, 'X': 48, 'Y': 2, 'Area': 77.14,
                   'Energy': -609.24, 'Label': '1 Species'}
    data = [Adsorbant_1]
    Surface_111_1 = wulff.calculate_surface_energy(stoich,
                                                   data,
                                                   SE,
                                                   adsorbant,
                                                   thermochem,
                                                   298,
                                                   0)

The with these surface energies we can build a particle morphology using pymatgen

.. code-block:: python

    lattice = Lattice.cubic(5.411)
    ceo = Structure(lattice,["Ce", "O"],
                   [[0,0,0], [0.25,0.25,0.25]])
    surface_energies_ceo = {(1,1,1): np.amin(Surface_111_1),
                            (1,1,0): np.amin(Surface_110_1),
                            (1,0,0): np.amin(Surface_100_1)}

    miller_list = surface_energies_ceo.keys()
    e_surf_list = surface_energies_ceo.values()

    wulffshape = WulffShape(ceo.lattice, miller_list, e_surf_list)
    wulffshape.show(color_set="RdBu", direction=(1.00, 0.25, 0.25))


.. image:: Figures/Tutorial_3/Wulff_2.png
    :height: 300px
    :align: center

Bulk Theory
===========

Bulk phase diagrams enable the comparison of the thermodynamic stability of various different bulk phases under different chemical potentials giving valuable insight in to the synthesis of solid phases.
This theory example will consider a series of bulk phases which can be defined through a reaction scheme across all phases,
thus for this example including MgO, :math:`H_2O` and :math:`CO_2` as reactions and A as a generic product.

.. math::
    x\text{MgO} + y\text{H}_2\text{O} + z\text{CO}_2 \rightarrow \text{A}

The system is in equilibrium when the chemical potentials of the reactants and product are equal; i.e. the change in Gibbs free energy is :math:`$\delta G_{T,p} = 0$`.

.. math::
	\delta G_{T,p} = \mu_A - x\mu_{\text{MgO}} - y\mu_{\text{H}_2\text{O}} - z\mu_{\text{CO}_2} = 0

Assuming that :math:`H_2O` and :math:`CO_2` are gaseous species, :math:`$\mu_{CO_2}$` and :math:`$\mu_{H_2O}$` can be written as

.. math::
	\mu_{\text{H}_2\text{O}} = \mu^0_{\text{H}_2\text{O}} + \Delta\mu_{\text{H}_2\text{O}}

and

.. math::
	\mu_{\text{CO}_2} = \mu^0_{\text{CO}_2} + \Delta\mu_{\text{CO}_2}

The chemical potential :math:`$\mu^0_x$` is the partial molar free energy of any reactants or products (x) in their standard states,
in this example we assume all solid components can be expressed as

.. math::
    \mu_{\text{component}} = \mu^0_{\text{component}}

Hence, we can now rearrange the equations to produce;

.. math::
	\mu^0_A - x\mu^0_{\text{MgO}} - y\mu^0_{\text{H}_2\text{O}} - z\mu^0_{CO_2} = y\Delta\mu_{H_2O} + z\Delta\mu_{CO_2}

As :math:`$\mu^0_A$` corresponds to the partial molar free energy of product A, we can replace the left side with the Gibbs free energy (:math:`$\Delta G_{\text{f}}^0$`).

.. math::
	\delta G_{T,p} = \Delta G_{\text{f}}^0 - y\Delta\mu_{\text{H}_2\text{O}} - z\Delta\mu_{\text{CO}_2}

At equilibrium :math:`$\delta G_{T,p} = 0$`, and hence

.. math::
	\Delta G_{\text{f}}^0 = y\Delta\mu_{\text{H}_2\text{O}} + z\Delta\mu_{\text{CO}_2}

Thus, we can find the values of :math:`$\Delta\mu_{H_2O}$` and :math:`$\Delta\mu_{CO_2}$` (or :math:`$(p_{H_2O})^y$` and :math:`$p_{CO_2}^z$` when Mg-rich phases are in thermodynamic equilibrium; i.e.
they are more or less stable than MgO.
This procedure can then be applied to all phases to identify which is the most stable, provided that the free energy :math:`$\Delta G_f^0$` is known for each Mg-rich phase.

The free energy can be calculated using

.. math::
    \Delta G^{0}_{f} = \sum\Delta G_{f}^{0,\text{products}} - \sum\Delta G_{f}^{0,\text{reactants}}

Where for this example the free energy (G) is equal to the calculated DFT energy (:math:`U_0`).

Temperature
~~~~~~~~~~~

The previous method will generate a phase diagram at 0 K. This is not representative of normal conditions.
Temperature is an important consideration for materials chemistry and we may wish to evaluate the phase thermodynamic stability at various synthesis conditions.

As before the free energy can be calculated using;

.. math::
    \Delta G^{0}_{f} = \sum\Delta G_{f}^{0,\text{products}} - \sum\Delta G_{f}^{0,\text{reactants}}

Where for this example the free energy (G) for solid phases is equal to is equal to the calculated DFT energy :math:`(U_0)`.
For gaseous species, the standard free energy varies significantly with temperature, and as DFT simulations are designed for condensed phase systems,
we use experimental data to determine the temperature dependent free energy term for gaseous species, where :math:`$S_{expt}(T)$` is specific entropy value for a given T and  :math:`$H-H^0(T)$` is the,
both can be obtained from the NIST database and can be calculated as;

.. math::
    G =  U_0 + (H-H^0(T) - T S_{\text{expt}}(T))

Pressure
--------

In the previous tutorials we went through the process of generating a simple phase diagram for bulk phases and introducing temperature dependence for gaseous species.
This useful however, sometimes it can be more beneficial to convert the chemical potentials (eVs) to partial pressure (bar).

Chemical potential can be converted to pressure values using

.. math::
    P & = \frac{\mu_O}{k_B T} ,

where P is the pressure, :math:`$\mu$` is the chemical potential of oxygen, $k_B$ is the Boltzmann constant and T is the temperature.
surfinpy\.data
==============

.. automodule:: surfinpy.data
    :members:
    :undoc-members:
    :show-inheritance:surfinpy
========

.. include:: README.rst

.. toctree::
   :maxdepth: 1

   theory.rst
   gallery.rst
   tutorials.rst
   using_surfinpy.rst
   modules.rst



indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Vibrational Entropy
===================

In this example we will expand this methodology to calculate the vibrational properties for solid phases (i.e. zero point energy, vibrational entropy)
and include these values in the generation of the phase diagrams.  This allows for a more accurate calculation of phase diagrams without the need to include experimental corrections for solid phases.

As with previous examples, the standard free energy varies significantly with temperature, and as DFT simulations are designed for condensed phase systems,
we use experimental data to determine the temperature dependent free energy term for gaseous species obtained from the NIST database.
In addition we also calculate the vibrational properties for the solid phases modifying the free energy (G) for solid phases to be;

.. math::
    \Delta G_f = U_0 + U_{\text{ZPE}} + A_{\text{vib}}

:math:`$U_0$` is the calculated internal energy from a DFT calculation, :math:`$U_{ZPE}$` is the zero point energy and :math:`$S_{vib}$` is the vibrational entropy.

.. math::
	U_{\text{ZPE}} = \sum_i^{3n} \frac{R \theta_i}{2}

where :math:`$A_{vib}$` is the vibrational Helmholtz free energy and defined as;

.. math::
	A_{\text{vib}} = \sum_i^{3n} RT \ln{(1-e^{-\theta_i/T})}


3n is the total number of vibrational modes, n is the number of species and :math:`$\theta_i$` is the characteristic vibrational temperature (frequency of the vibrational mode in Kelvin).

.. math::
	\theta_i = \frac{h\nu_i}{k_B}

.. code-block:: python

    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import data

    temperature_range = [298, 299]

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10, file = 'bulk_vib.yaml', entropy=True, zpe=True, temp_range=temperature_range)


In addition to entropy and zpe keyword you must provide the a file containing the vibrational modes and number of formula units used in that calculations.  You must create the yaml file using the following format

.. code-block:: python

    F-Units : number
    Frequencies :
    - mode1
    - mode2

Vibrational modes can be calculated via a density functional perturbation calculation or via the phonopy code.

.. code-block:: python

    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92,
                      label = "Bulk", entropy = True, zpe=True, file = 'ref_files/bulk_vib.yaml',
                      funits = 10, temp_range=temperature_range)

    A = data.DataSet(cation = 10, x = 5, y = 20, energy = -468,
                      label = "A", entropy = True, zpe=True, file = 'ref_files/A_vib.yaml',
                      funits = 5, temp_range=temperature_range)

    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -228,
                      label = "B", entropy = True, zpe=True, file = 'ref_files/B_vib.yaml',
                      funits =  10, temp_range=temperature_range)

    C = data.DataSet(cation = 10, x = 10, y = 30, energy = -706,
                      label = "C", entropy = True, zpe=True, file = 'ref_files/C_vib.yaml',
                      funits = 10, temp_range=temperature_range)

    D = data.DataSet(cation = 10, x = 10, y = 0, energy = -310,
                      label = "D", entropy = True, zpe=True,  file = 'ref_files/D_vib.yaml',
                      funits =  10, temp_range=temperature_range)

    E = data.DataSet(cation = 10, x = 10, y = 50, energy = -972,
                      label = "E", entropy = True, zpe=True, file = 'ref_files/E_vib.yaml',
                      funits =  10, temp_range=temperature_range)

    F = data.DataSet(cation = 10, x = 8, y = 10, energy = -398,
                      label = "F", entropy = True, zpe=True, file = 'ref_files/F_vib.yaml',
                      funits =  2, temp_range=temperature_range)


    data = [Bulk, A, B, C, D, E, F]

    x_energy=-20.53412969
    y_energy=-12.83725889


    CO2_exp = ut.fit_nist("CO2.txt")[298]
    Water_exp = ut.fit_nist("H2O.txt")[298]

    CO2_corrected = x_energy + CO2_exp
    Water_corrected = y_energy + Water_exp

    deltaX = {'Range': [ -1, 0.6],  'Label': 'CO_2'}
    deltaY = {'Range': [ -1, 0.6], 'Label': 'H_2O'}

    temp_298 = bmvm.calculate(data, bulk, deltaX, deltaY, CO2_corrected, Water_corrected)
    ax = temp_298.plot_mu_p(temperature=298, set_style="fast", colourmap="RdBu")

.. image:: Figures/Bulk_8.png
    :height: 300px
    :align: center


Temperature
-----------

In tutorial 5 we showed how SurfinPy can be used to calculate the vibrational entropy and zero point energy for solid phases and in tutorial 4 we showed how a temperature range can be used to calculate the phase diagram of temperature as a function of pressure.  In this example we will use both lesson from these tutorials to produce a phase diagram of temperature as a function of pressure including the vibrational properties for solid phases.  Again this produces results which are easily compared to experimental values in addition to increasing the level of theory used.

.. code-block:: python

    import matplotlib.pyplot as plt
    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import bulk_mu_vs_t as bmvt
    from surfinpy import data
    import numpy as np

    colors = ['#5B9BD5', '#4472C4', '#A5A5A5', '#772C24', '#ED7D31', '#FFC000', '#70AD47']

    temperature_range = [273, 373]

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10, file = 'bulk_vib.yaml', entropy=True, zpe=True, temp_range=temperature_range)


    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92., color=colors[0],
                    label = "Bulk", entropy = True, zpe=True, file = 'ref_files/bulk_vib.yaml',
                    funits = 10, temp_range=temperature_range)

    D = data.DataSet(cation = 10, x = 10, y = 0, energy = -310.,  color=colors[1],
                    label = "D", entropy = True, zpe=True,  file = 'ref_files/D_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -227.,  color=colors[2],
                    label = "B", entropy = True, zpe=True, file = 'ref_files/B_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    F = data.DataSet(cation = 10, x = 8, y = 10, energy = -398.,  color=colors[3],
                    label = "F", entropy = True, zpe=True, file = 'ref_files/F_vib.yaml',
                    funits =  2, temp_range=temperature_range)

    A = data.DataSet(cation = 10, x = 5, y = 20, energy = -467.,  color=colors[4],
                    label = "A", entropy = True, zpe=True, file = 'ref_files/A_vib.yaml',
                    funits = 5, temp_range=temperature_range)


    C = data.DataSet(cation = 10, x = 10, y = 30, energy = -705.,  color=colors[5],
                    label = "C", entropy = True, zpe=True, file = 'ref_files/C_vib.yaml',
                    funits = 10, temp_range=temperature_range)

    E = data.DataSet(cation = 10, x = 10, y = 50, energy = -971.,  color=colors[6],
                    label = "E", entropy = True, zpe=True, file = 'ref_files/E_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    data = [Bulk, A, B, C,  D, E, F]

    deltaX = {'Range': [ -1, 0.6],  'Label': 'CO_2'}
    deltaZ = {'Range': [ 273, 373], 'Label': 'Temperature'}
    x_energy=-20.53412969
    y_energy=-12.83725889
    mu_y = 0


    exp_x = ut.temperature_correction_range("CO2.txt", deltaZ)
    exp_y = ut.temperature_correction_range("H2O.txt", deltaZ)

    system = bmvt.calculate(data, bulk, deltaX, deltaZ, x_energy, y_energy, mu_y, exp_x, exp_y)
    ax = system.plot_mu_vs_t_vs_p(temperature=273)

.. image:: Figures/Bulk_9.png
    :height: 300px
    :align: center

When investigating the phase diagram for certain systems it could be beneficial to remove a kinetically inhibited but thermodynamically stable phase to investigate the metastable phase diagram.  Within SurfinPy this can be achieved via recreating the data list without the phase in question then recalculating the phase diagram, as below.

.. code-block:: python

    data = [Bulk, A, B, C, E, F]

    system = bmvt.calculate(data, bulk, deltaX, deltaZ, x_energy, y_energy, mu_y, exp_x, exp_y)
    ax = system.plot_mu_vs_t_vs_p(temperature=273)

.. image:: Figures/Bulk_10.png
    :height: 300px
    :align: center
Vibrational Theory
==================

The vibrational entropy allows for a more accurate calculation of phase diagrams without the need to include experimental corrections for solid phases.

The standard free energy varies significantly with temperature, and as DFT simulations are designed for condensed phase systems, 
we use experimental data to determine the temperature dependent free energy term for gaseous species obtained from the NIST database.  
In addition we also calculate the vibrational properties for the solid phases modifying the free energy (G) for solid phases to be;

.. math::
    \Delta G_f = U_0 + U_{\text{ZPE}} + A_{\text{vib}}

:math:`$U_0$` is the calculated internal energy from a DFT calculation, :math:`$U_{ZPE}$` is the zero point energy and :math:`$S_{vib}$` is the vibrational entropy.

.. math::
	U_{\text{ZPE}} = \sum_i^{3n} \frac{R \theta_i}{2}

where :math:`$A_{vib}$` is the vibrational Helmholtz free energy and defined as;

.. math::
	A_{\text{vib}} = \sum_i^{3n} RT \ln{(1-e^{-\theta_i/T})}

3n is the total number of vibrational modes, n is the number of species and :math:`$\theta_i$` is the characteristic vibrational temperature (frequency of the vibrational mode in Kelvin).

.. math::
	\theta_i = \frac{h\nu_i}{k_B}
surfinpy\.wulff
===============

The module required for the generation of wulff plots. An explanation of theory can be found `here <theory.html>`

.. automodule:: surfinpy.wulff
    :members:
    :undoc-members:
    :show-inheritance:
Theory
======

There is a significant amount of theory behind the methods in `surfinpy`. The following three pages provide an
explanation for the methods employed in the code.

`surfinpy` is a Python module to generate phase diagrams from energy minimisation data. 
Before using this code you will need to generate the relevant data. 


.. toctree::
   :maxdepth: 1

   surface_theory.rst
   bulk_theory.rst
   vibrational_theory.rst



.. image:: Figures/Logo.png
    :align: center
    :alt: Project Logo

.. image::  https://readthedocs.org/projects/surfinpy/badge/?version=latest
    :target: https://surfinpy.readthedocs.io/en/latest/
    :alt: Documentation Status

.. image:: https://travis-ci.com/symmy596/SurfinPy.svg?branch=master
    :target: https://travis-ci.com/symmy596/SurfinPy
    :alt: Build Status

This is the documentation for the open-source Python project, `SurfinPy`.
A library designed to facilitate the generation of publication ready phase diagrams from ab initio calculations.
`SurfinPy` is built on existing Python packages that those in the solid state physics/chemistry community should already be familiar with. 
This tool will bring some benfits to the solid state community and facilitate the generation of publication ready phase diagrams (powered by Matplotlib).

The main features include:

1. **Method to generate surface phase diagrams as a function of chemical potential.**  

   - Generate a diagram as a function of the chemical potential of two adsorbing species e.g. water and carbon dioxide.  
   - Generate a diagram as a function of the chemical potential of one adsorbing species and a surface species e.g. water and oxygen vacancies.  
   - Use experimental data combined with ab initio data to generate a temperature dependent phase diagram.  

2. **Method to generate surface phase diagrams as a function of temperature and pressure.**  

   - Use experimental data combined with ab initio data to generate a pressure vs temperature plot showing the state (e.g. composition) of a surface as a function of temperature and pressure of one species.

3. **Use calculated surface energies to built crystal morphologies.**  

   - Use the surface energies produced by `SurfinPy` alongside Pymatgen to built particle morphologies.  
   - Evaulate how a particles shape changes with temperature and pressure.

4. **Method to generate bulk phase diagrams as a function of chemical potential.**  

   - Generate a diagram as a function of the chemical potential of two species e.g. water and carbon dioxide.  
   - Use experimental data combined with ab initio data to generate a temperature dependent phase diagram.  

5. **Method to generate bulk phase diagrams as a function of temperature and pressure.**  

   - Use experimental data combined with ab initio data to generate a pressure vs temperature plot showing the phase space as a function of temperature and pressure.  

6. **Method to include vibrational properties in a phase diagram**
   
   - Module to calculate the zero point energy and vibrational entropy
   - Encorporate the zero point energy and/or the vibrational entropy into a phase diagram.

The code has been developed to be used with any ab initio code as it only requires a list of energies and vibrational frequencies.  
`surfinpy` was developed across several PhD projects and as such the functionality focuses on the research questions encountered during those projects, which we should clarify 
are wide ranging. Code contributions aimed at expanding the code to new problems are encouraged.

`surfinpy` is free to use.

Usage
-----

A full list of examples can be found in the examples folder of the git repository, these include jupyter notebook tutorials which combine the full theory with code examples.

Installation
------------

surfinpy is a Python 3 package and requires a typical scientific Python stack. Use of the tutorials requires Anaconda/Jupyter to be installed.

To build from source:

.. code-block:: bash

    pip install -r requirements.txt

    python setup.py build

    python setup.py install

    python setup.py test

Or alternatively install with pip

.. code-block:: bash

    pip install surfinpy

Documentation
-------------

To build the documentation from scratch 

.. code-block:: bash

    cd docs
    
    make html

License
-------

`surfinpy` is made available under the MIT License.


Detailed requirements
---------------------

`surfinpy` is compatible with Python 3.5+ and relies on a number of open source Python packages, specifically:

- Numpy
- Scipy
- Matplotlib
- Pymatgen

Contributing
------------

Contact
~~~~~~~

If you have questions regarding any aspect of the software then please get in touch with the development team via email Adam Symington (symmy596@gmail.com), Joshua Tse (joshua.s.tse@gmail.com). 
Alternatively you can create an issue on the `Issue Tracker <https://github.com/symmy596/SurfinPy/issues>`_ or you can discuss your questions on our `gitter channel <https://gitter.im/Surfinpy/Lobby>`_.

Bugs 
~~~~

There may be bugs. If you think you have caught one, please report it on the `Issue Tracker <https://github.com/symmy596/SurfinPy/issues>`_.
This is also the place to propose new ideas for features or ask questions about the design of `surfinpy`. Poor documentation is considered a bug 
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



For further information please contact Adam Symington - symmy596@gmail.com, Joshua Tse - joshua.s.tse@gmail.com


Research
--------

- `Strongly Bound Surface Water Affects the Shape Evolution of Cerium Oxide Nanoparticles <https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.9b09046>`__
- `The energetics of carbonated PuO2  surfaces affects nanoparticle morphology: a DFT+U study <https://pubs.rsc.org/lv/content/articlelanding/2020/cp/d0cp00021c/unauth#!divAbstract>`__
- `Exploiting cationic vacancies for increased energy densities in dual-ion batteries <https://www.sciencedirect.com/science/article/abs/pii/S2405829719310153>`__
- `Thermodynamic Evolution of Cerium Oxide Nanoparticle Morphology Using Carbon Dioxide <https://pubs.acs.org/doi/abs/10.1021/acs.jpcc.0c07437>`__ Using surfinpy
==============

The are number of ways to get and use surfinpy

- Fork the code: please feel free to fork the code on [GitHub](https://github.com/symmy596/SurfinPy)
  and add functionality that interests you. There are already plans for version 2, these will all be added
  to the issues section in due course. 
- Run it locally: surfinpy is available through the pip package manager.
- Get in touch: Adam R.Symington (ars44@bath.ac.uk) is always keen to chat to potential users.
Surface Theory
==============

`surfinpy` has the capability to generate phase diagrams as a function of chemical potential of two varying species e.g. water and carbon dioxide. In such an example
the user would require calculations with varying concentrations of water, carbon dioxide and water/carbon dioxide on a surface. Assuming that you have generated enough,
reliable data then you are ready to use `surfinpy`.

Surface Energy
~~~~~~~~~~~~~~

The physical quantity that is used to define the stability of a surface with a given composition is its surface energy :math:`\gamma` (J :math:`m^{-2}`).
At its core, `surfinpy` is a code that calculates the surface energy of different slabs at varying chemical potential and
uses these surface energies to construct a phase diagram.
In this explanation of theory we will use the example of water adsorbing onto a surface of :math:`TiO_2` containing oxygen vacancies.
In such an example there are two variables, water concentration and oxygen concentration. We are able to calculate the surface energy according to

.. math::
    \gamma_{Surf} = \frac{1}{2A} \Bigg( E_{TiO_2}^{slab} - \frac{nTi_{slab}}{nTi_{Bulk}} E_{TiO_2}^{Bulk} \Bigg) - \Gamma_O \mu_O - \Gamma_{H_2O} \mu_{H_2O} ,

where A is the surface area, :math:`E_{TiO_2}^{slab}` is the DFT energy of the slab, :math:`nTi_{Slab}` is the number of cations in the slab,
:math:`nTi_{Bulk}` is the number of cations in the bulk unit cell, :math:`E_{TiO_2}^{Bulk}` is the DFT energy of the bulk unit cell,

.. math::
    \Gamma_O = \frac{1}{2A} \Bigg( nO_{Slab} - \frac{nO_{Bulk}}{nTi_{Bulk}}nTi_{Slab}  \Bigg) ,

.. math::
    \Gamma_{H_2O} = \frac{nH_2O}{2A},

where :math:`nO_{Slab}` is the number of anions in the slab, :math:`nO_{Bulk}` is the number of anions in the bulk and :math:`nH_2O` is the number of adsorbing water molecules.
:math:`\Gamma_O` / :math:`\Gamma_{H_2O}` is the excess oxygen / water at the surface and :math:`\mu_O` / :math:`\mu_{H_2O}` is the oxygen / water chemical potential.
Clearly :math:`\Gamma` and :math:`mu` will only matter when the surface is non stoichiometric.

Temperature
~~~~~~~~~~~

The above phase diagram is at 0K. It is possible to use experimental data from the NIST_JANAF database to make the chemical potential a temperature dependent
term and thus generate a phase diagram at a temperature (T). This is done according to

.. math::
    \gamma_{Surf} = \frac{1}{2A} \Bigg( E_{TiO_2}^{slab} - \frac{nTi_{Slab}}{nTi_{Bulk}} E_{TiO_2}^{Bulk} \Bigg) - \Gamma_O \mu_O - \Gamma_{H_2O} \mu_{H_2O} - n_O \mu_O (T) - n_{H_2O} \mu_{H_2O} (T)

where

.. math::
    \mu_O (T)  = \frac{1}{2} \mu_O (T) (0 K , DFT) +  \frac{1}{2} \mu_O (T) (0 K , EXP) +  \frac{1}{2} \Delta G_{O_2} ( \Delta T, Exp),

:math:`\mu_O` (T) (0 K , DFT) is the 0K free energy of an isolated oxygen molecule evaluated with DFT, :math:`\mu_O` (T) (0 K , EXP) is the 0 K experimental
Gibbs energy for oxygen gas and $\Delta$ :math:`G_{O_2}` ( :math:`\Delta` T, Exp) is the Gibbs energy defined at temperature T as

.. math::
    \Delta G_{O_2} ( \Delta T, Exp)  = \frac{1}{2} [H(T, {O_2}) -  H(0 K, {O_2})] -  \frac{1}{2} T[S(T, {O_2}]).

This will generate a phase diagram at temperature (T)

Pressure
~~~~~~~~

Chemical potential can be converted to pressure values according to

.. math::
    P = \frac{\mu_X}{k_B T}

where P is the pressure, :math:`\mu` is the chemical potential of species X, :math:`k_B` is the Boltzmann constant and T is the temperature.

Pressure vs temperature
~~~~~~~~~~~~~~~~~~~~~~~

`Surfinpy` has the functionality to generate phase diagrams as a function of pressure vs temperature based upon the methodology used in `Molinari et al
<https://pubs.acs.org/doi/abs/10.1021/jp300576b>`_ according to

.. math::
    \gamma_{adsorbed, T, P} = \gamma_{bare} + ( C ( E_{ads, T} - RTln(\frac{p}{p^o})

where :math:`\gamma_{adsorbed, T, p}` is the surface energy of the surface with adsorbed species at temperature (T) and pressure (P),
:math:`\gamma_{bare}` is the surface energy of the bare surface, C is the coverage of adsorbed species, :math:`E_{ads}` is the adsorption energy,

.. math::
    E_{ads, T} =  E_{slab, adsorbant} - (E_{slab, bare} + n_{H_2O} E_{H_2O, T}) / n_{H_2O}

where :math:`E_{slab, adsorbant}` is the energy of the surface and the adsorbed species, :math:`n_{H_2O}` is he number of adsorbed species,

.. math::
    E_{H_2O, (T)} = E_{H_2O, (g)} - TS_{(T)}

where :math:`S_{(T)}` is the experimental entropy of gaseous water in the standard state.
surfinpy\.vibrational_data
==========================

.. automodule:: surfinpy.vibrational_data
    :members:
    :undoc-members:
    :show-inheritance:Tutorials
=========

These tutorials are replicated in jupyter notebook form and contained within `examples <https://github.com/symmy596/SurfinPy/tree/master/examples/>`_.
The accompanying python scripts are also included here. All of code examples within these tutorials can be found in
`examples/Scripts <https://github.com/symmy596/SurfinPy/tree/master/examples/Scripts>`_.


Surfaces
--------

.. toctree::
   :maxdepth: 1

   tutorial_1.rst
   tutorial_2.rst
   tutorial_3.rst

Bulk
----

.. toctree::
   :maxdepth: 1

   tutorial_4.rst
   tutorial_5.rst
   tutorial_6.rst
surfinpy\.bulk_mu_vs_t
======================

.. automodule:: surfinpy.bulk_mu_vs_t
    :members:
    :undoc-members:
    :show-inheritance:Chemical Potential
==================

In this tutorial we learn how to generate a basic bulk phase diagram from DFT energies.  This enables the comparison of the thermodynamic stability of various different bulk phases under different chemical potentials giving valuable insight in to the synthesis of solid phases.  This example will consider a series of bulk phases which can be defined through a reaction scheme across all phases, thus for this example including Bulk, :math:`H_2O` and :math:`CO_2` as reactions and A as a generic product.

.. math::
    x\text{Bulk} + y\text{H}_2\text{O} + z\text{CO}_2 \rightarrow \text{A}


The system is in equilibrium when the chemical potentials of the reactants and product are equal; i.e. the change in Gibbs free energy is :math:`$\delta G_{T,p} = 0$`.

.. math::
	\delta G_{T,p} = \mu_A - x\mu_{\text{Bulk}} - y\mu_{\text{H}_2\text{O}} - z\mu_{\text{CO}_2} = 0

Assuming that :math:`H_2O` and :math:`CO_2` are gaseous species, :math:`$\mu_{CO_2}$` and :math:`$\mu_{H_2O}$` can be written as

.. math::
	\mu_{\text{H}_2\text{O}} = \mu^0_{\text{H}_2\text{O}} + \Delta\mu_{\text{H}_2\text{O}}

and

.. math::
	\mu_{\text{CO}_2} = \mu^0_{\text{CO}_2} + \Delta\mu_{\text{CO}_2}

The chemical potential :math:`$\mu^0_x$` is the partial molar free energy of any reactants or products (x) in their standard states, in this example we assume all solid components can be expressed as

.. math::
    \mu_{\text{component}} = \mu^0_{\text{component}}

Hence, we can now rearrange the equations to produce;

.. math::
	\mu^0_A - x\mu^0_{\text{Bulk}} - y\mu^0_{\text{H}_2\text{O}} - z\mu^0_{CO_2} = y\Delta\mu_{H_2O} + z\Delta\mu_{CO_2}

As :math:`$\mu^0_A$` corresponds to the partial molar free energy of product A, we can replace the left side with the Gibbs free energy (:math:`$\Delta G_{\text{f}}^0$`).

.. math::
	\delta G_{T,p} = \Delta G_{\text{f}}^0 - y\Delta\mu_{\text{H}_2\text{O}} - z\Delta\mu_{\text{CO}_2}

At equilibrium :math:`$\delta G_{T,p} = 0$`, and hence

.. math::
	\Delta G_{\text{f}}^0 = y\Delta\mu_{\text{H}_2\text{O}} + z\Delta\mu_{\text{CO}_2}

Thus, we can find the values of :math:`$\Delta\mu_{H_2O}$` and :math:`$\Delta\mu_{CO_2}$` (or :math:`$(p_{H_2O})^y$` and :math:`$p_{CO_2}^z$` when solid phases are in thermodynamic equilibrium; i.e. they are more or less stable than Bulk.  This procedure can then be applied to all phases to identify which is the most stable, provided that the free energy :math:`$\Delta G_f^0$` is known for each solid phase.

The free energy can be calculated using

.. math::
    \Delta G^{0}_{f} = \sum\Delta G_{f}^{0,\text{products}} - \sum\Delta G_{f}^{0,\text{reactants}}

Where for this tutorial the free energy (G) is equal to the calculated DFT energy (:math:`U_0`).

.. code-block:: python

    import matplotlib.pyplot as plt
    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import data


The first thing to do is input the data that we have generated from our DFT simulations. The input data needs to be contained within a class. First we have created the class for the bulk data, where 'Cation' is the number of cations, 'Anion' is the number of anions, 'Energy' is the DFT energy and 'F-Units' is the number of formula units.

.. code-block:: python

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10)


Next we create the bulk phases classes - one for each phase (A-F). 'Cation' is the number of cations, 'x' is in this case the number of water species (corresponding to the X axis of the phase diagram), 'y' is the number of in this case :math:`CO_2` molecules (corresponding to the Y axis of our phase diagram), 'Energy' is the DFT energy and finally 'Label' is the label for the phase (appears on the phase diagram).


.. code-block:: python

    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92.0, label = "Bulk")
    A = data.DataSet(cation = 10, x = 5, y = 20, energy = -468.0, label = "A")
    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -228.0, label = "B")
    C = data.DataSet(cation = 10, x = 10, y = 30, energy = -706.0, label = "C")
    D = data.DataSet(cation = 10, x = 10, y = 0, energy = -310.0, label = "D")
    E = data.DataSet(cation = 10, x = 10, y = 50, energy = -972.0, label = "E")
    F = data.DataSet(cation = 10, x = 8, y = 10, energy = -398.0, label = "F")


Next we need to create a list of our data. Don't worry about the order, `surfinpy` will sort that out for you.

.. code-block:: python

    data = [Bulk, A, B, C,  D, E, F]


We now need to generate our X and Y axis, or more appropriately, our chemical potential values. These exist in a dictionary. 'Range' corresponds to the range of chemical potential values to be considered and 'Label' is the axis label.  Additionally, the x and y energy need to be specified.

.. code-block:: python

    deltaX = {'Range': Range of Chemical Potential,
              'Label': Species Label}


.. code-block:: python

    deltaX = {'Range': [ -3, 2],  'Label': 'CO_2'}
    deltaY = {'Range': [ -3, 2], 'Label': 'H_2O'}
    x_energy=-20.53412969
    y_energy=-12.83725889


And finally we can generate our plot using these 6 variables of data.

.. code-block:: python

    system = bmvm.calculate(data, bulk, deltaX, deltaY, x_energy, y_energy)

    ax = system.plot_phase()
    plt.show()

.. image:: Figures/Bulk_1.png
    :height: 300px
    :align: center


Temperature
-----------

In the previous example we generated a phase diagram at 0 K.  However, this is not representative of normal conditions.
Temperature is an important consideration for materials chemistry and we may wish to evaluate the phase thermodynamic stability at various synthesis conditions.
This example will again be using the :math:`Bulk-CO_2-H_2O` system.

As before the free energy can be calculated using;

.. math::
    \Delta G^{0}_{f} = \sum\Delta G_{f}^{0,\text{products}} - \sum\Delta G_{f}^{0,\text{reactants}}

Where for this tutorial the free energy (G) for solid phases  is equal to is equal to the calculated DFT energy :math:`(U_0)`.
For gaseous species, the standard free energy varies significantly with temperature, and as DFT simulations are designed for condensed phase systems,
we use experimental data to determine the temperature dependent free energy term for gaseous species,
where :math:`$S_{expt}(T)$` is specific entropy value for a given T and  :math:`$H-H^0(T)$` is the , both can be obtained from the NIST database and can be calculated as;

.. math::
    G =  U_0 + (H-H^0(T) - T S_{\text{expt}}(T))

.. code-block:: python

    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import data

.. code-block:: python

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10)

    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92.0, label = "Bulk")
    A = data.DataSet(cation = 10, x = 5, y = 20, energy = -468.0, label = "A")
    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -228.0, label = "B")
    C = data.DataSet(cation = 10, x = 10, y = 30, energy = -706.0, label = "C")
    D = data.DataSet(cation = 10, x = 10, y = 0, energy = -310.0, label = "D")
    E = data.DataSet(cation = 10, x = 10, y = 50, energy = -972.0, label = "E")
    F = data.DataSet(cation = 10, x = 8, y = 10, energy = -398.0, label = "F")
    data = [Bulk, A, B, C,  D, E, F]

    x_energy=-20.53412969
    y_energy=-12.83725889

In order to calculate :math:`$S_{expt}(T)$` for :math:`H_2O` and :math:`CO_2` we need to use experimental data from the NSIT JANAF database.
As a user you will need to download the tables for the species you are interested in (in this example water and carbon dioxide).
`surfinpy` has a function that can read this data, assuming it is in the correct format and calculate the temperature correction for you.
Provide the path to the file and the temperature you want.

.. code-block:: python

    CO2_exp = ut.fit_nist("CO2.txt")[298]
    Water_exp = ut.fit_nist("H2O.txt")[298]

    CO2_corrected = x_energy + CO2_exp
    Water_corrected = y_energy + Water_exp

    deltaX = {'Range': [ -3, 2],  'Label': 'CO_2'}
    deltaY = {'Range': [ -3, 2], 'Label': 'H_2O'}

CO2_corrected and H2O_corrected are now temperature dependent terms corresponding to a temperature of 298 K. The resulting phase diagram will now be at a temperature of 298 K.

.. code-block:: python

    system = bmvm.calculate(data, bulk, deltaX, deltaY, x_energy=CO2_corrected, y_energy=Water_corrected)

    system.plot_phase(temperature=298)

.. image:: Figures/Bulk_2.png
    :height: 300px
    :align: center

Pressure
--------

In the previous example we went through the process of generating a simple phase diagram for bulk phases and introducing temperature dependence for gaseous species.
This useful however, sometimes it can be more beneficial to convert the chemical potentials (eVs) to partial pressure (bar).

Chemical potential can be converted to pressure values using

.. math::
    P & = \frac{\mu_O}{k_B T} ,

where P is the pressure, :math:`$\mu$` is the chemical potential of oxygen, $k_B$ is the Boltzmann constant and T is the temperature.

.. code-block:: python

    import matplotlib.pyplot as plt
    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import data

    colors = ['#5B9BD5', '#4472C4', '#A5A5A5', '#772C24', '#ED7D31', '#FFC000', '#70AD47']


Additionally, `surfinpy` has the functionality to allow you to choose which colours are used for each phase.  Specify within the DataSet class color.

.. code-block:: python

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10)

    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92.0, color=colors[0], label = "Bulk")
    A = data.DataSet(cation = 10, x = 10, y = 0, energy = -310.0, color=colors[1], label = "A")
    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -228.0, color=colors[2], label = "B")
    C = data.DataSet(cation = 10, x = 8, y = 10, energy = -398.0, color=colors[3], label = "C")
    D = data.DataSet(cation = 10, x = 5, y = 20, energy = -468.0, color=colors[4], label = "D")
    E = data.DataSet(cation = 10, x = 10, y = 30, energy = -706.0, color=colors[5], label = "E")
    F = data.DataSet(cation = 10, x = 10, y = 50, energy = -972.0, color=colors[6], label = "F")

    data = [Bulk, A, B, C,  D, E, F]

    x_energy=-20.53412969
    y_energy=-12.83725889

    CO2_exp = ut.fit_nist("CO2.txt")[298]
    Water_exp = ut.fit_nist("H2O.txt")[298]

    CO2_corrected = x_energy + CO2_exp
    Water_corrected = y_energy + Water_exp

    deltaX = {'Range': [ -1, 0.6],  'Label': 'CO_2'}
    deltaY = {'Range': [ -1, 0.6], 'Label': 'H_2O'}

    system = bmvm.calculate(data, bulk, deltaX, deltaY, x_energy=CO2_corrected, y_energy=Water_corrected)

    system.plot_phase()

.. image:: Figures/Bulk_3.png
    :height: 300px
    :align: center

To convert chemical potential to pressure use the plot_pressure command and the temperature at which the pressure is calculated.  For this example we have used 298 K.

.. code-block:: python

    system.plot_pressure(temperature=298)

.. image:: Figures/Bulk_4.png
    :height: 300px
    :align: center
Pressure vs Temperature
=======================

In the previous example, we showed how experimental data could be used to determine the temperature dependent free energy term for gaseous species and then plot a phase diagram that represents 298 K.
This same method can be used in conjunction with a temperature range to produce a phase diagram of temperature as a function of pressure (or chemical potential).
This is an important step to producing relatable phase diagrams that can be compared to experimental findings.

To reiterate, the free energy can be calculated using;

.. math::
    \Delta G^{0}_{f} = \sum\Delta G_{f}^{0,\text{products}} - \sum\Delta G_{f}^{0,\text{reactants}}

Where for this tutorial the free energy (G) for solid phases  is equal to is equal to the calculated DFT energy (:math:`U_0`). For gaseous species,
the standard free energy varies significantly with temperature, and as DFT simulations are designed for condensed phase systems,
we use experimental data to determine the temperature dependent free energy term for gaseous species, where $S_{\text{expt}}(T)$ is specific entropy value for a given T and  $H-H^0(T)$ is the ,
both can be obtained from the NIST database and can be calculated as;

.. math::
    G =  U_0 + (H-H^0(T) - T S_{\text{expt}}(T))

.. code-block:: python

    from surfinpy import bulk_mu_vs_mu as bmvm
    from surfinpy import utils as ut
    from surfinpy import bulk_mu_vs_t as bmvt
    from surfinpy import data
    import numpy as np
    from seaborn import palettes

    import matplotlib.pyplot as plt
    colors = ['#5B9BD5', '#4472C4', '#A5A5A5', '#772C24', '#ED7D31', '#FFC000', '#70AD47']


temperature_range sets the temperature range which is calculated for the phase diagram and needs to be specified within the data ReferenceDataSet and DataSet.

.. code-block:: python

    temperature_range = [200, 400]

    bulk = data.ReferenceDataSet(cation = 1, anion = 1, energy = -92.0, funits = 10, file = 'bulk_vib.yaml', temp_range=temperature_range)


    Bulk = data.DataSet(cation = 10, x = 0, y = 0, energy = -92., color=colors[0],
                    label = "Bulk", file = 'ref_files/bulk_vib.yaml'',
                    funits = 10, temp_range=temperature_range)

    D = data.DataSet(cation = 10, x = 10, y = 0, energy = -310.,  color=colors[1],
                    label = "D", file = 'ref_files/D_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    B = data.DataSet(cation = 10, x = 0, y = 10, energy = -227.,  color=colors[2],
                    label = "B", file = 'ref_files/B_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    F = data.DataSet(cation = 10, x = 8, y = 10, energy = -398.,  color=colors[3],
                    label = "F", file = 'ref_files/F_vib.yaml',
                    funits =  2, temp_range=temperature_range)

    A = data.DataSet(cation = 10, x = 5, y = 20, energy = -467.,  color=colors[4],
                    label = "A", file = 'ref_files/A_vib.yaml',
                    funits = 5, temp_range=temperature_range)

    C = data.DataSet(cation = 10, x = 10, y = 30, energy = -705.,  color=colors[5],
                    label = "C", file = 'ref_files/C_vib.yaml',
                    funits = 10, temp_range=temperature_range)

    E = data.DataSet(cation = 10, x = 10, y = 50, energy = -971.,  color=colors[6],
                    label = "E", file = 'ref_files/E_vib.yaml',
                    funits =  10, temp_range=temperature_range)

    data = [Bulk, A, B, C,  D, E, F]

deltaZ specifies the temperature range which is plotted (Note that this must be the same as temperature_range).
mu_y is the chemical potential (eV) of third component, in this example we use a chemical potential of water = 0 eV which is equivalent to 1 bar pressure.

.. code-block:: python

    deltaX = {'Range': [ -1, 0.6],  'Label': 'CO_2'}
    deltaZ = {'Range': [ 200, 400], 'Label': 'Temperature'}
    x_energy=-20.53412969
    y_energy=-12.83725889
    mu_y = 0

    exp_x = ut.temperature_correction_range("CO2.txt", deltaZ)
    exp_y = ut.temperature_correction_range("H2O.txt", deltaZ)

    system = bmvt.calculate(data, bulk, deltaX, deltaZ, x_energy, y_energy, mu_y, exp_x, exp_y)
    ax = system.plot_mu_vs_t()

.. image:: Figures/Bulk_5.png
    :height: 300px
    :align: center

.. code-block:: python

    system.plot_p_vs_t(temperature=273, set_style="seaborn-dark-palette", colourmap="RdYlBu")

.. image:: Figures/Bulk_6.png
    :height: 300px
    :align: center


.. code-block:: python

    system.plot_mu_vs_t_vs_p(temperature=273, set_style="seaborn-dark-palette", colourmap="RdYlBu")

.. image:: Figures/Bulk_7.png
    :height: 300px
    :align: center
surfinpy\.mu_vs_mu
==================

Functions related to the generation of surface phase diagrams as a function of chemical potential.
An explanation of theory can be found `here <theory.html>`_

.. automodule:: surfinpy.mu_vs_mu
    :members:
    :undoc-members:
    :show-inheritance:surfinpy\.p_vs_t
================

Functions related to the generation of surface phase diagrams as a function of pressure and temperature.
An explanation of theory can be found `here <theory.html>`_

.. automodule:: surfinpy.p_vs_t
    :members:
    :undoc-members:
    :show-inheritance: