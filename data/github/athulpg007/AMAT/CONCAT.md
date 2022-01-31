# AMAT

Aerocapture Mission Analysis Tool (AMAT) is designed to provide rapid mission analysis capability for aerocapture and Entry, Descent, and Landing (EDL) mission concepts to the planetary science community. AMAT was developed at the [Advanced Astrodynamics Concepts (AAC)](https://engineering.purdue.edu/AAC/) research group at Purdue University.


See [AMAT documentation](https://amat.readthedocs.io) for more details. 

[![Documentation Status](https://readthedocs.org/projects/amat/badge/?version=master)](https://amat.readthedocs.io/en/master/?badge=master) [![DOI](https://joss.theoj.org/papers/10.21105/joss.03710/status.svg)](https://doi.org/10.21105/joss.03710) [![PyPI version](https://badge.fury.io/py/AMAT.svg)](https://badge.fury.io/py/AMAT)

If you find AMAT useful in your work, please consider citing us: Girija et al., (2021). AMAT: A Python package for rapid conceptual design of aerocapture and atmospheric Entry, Descent, and Landing (EDL) missions in a Jupyter environment. *Journal of Open Source Software*, 6(67), 3710, [DOI: 10.21105/joss.03710](https://doi.org/10.21105/joss.03710)


Please note the public release version has a minimal working interplanetary dataset to run some of the example Jupyter notebooks because of data sharing policies of external collaborators.

## Capabilities

AMAT allows the user to perform low-fidelity broad sweep parametric studies; as well as high fidelity Monte Carlo simulations to quantify aerocapture performance. AMAT supports analysis for all atmosphere-bearing destinations in the Solar System: Venus, Earth, Mars, Jupiter, Saturn, Titan, Uranus, and Neptune.

### Venus Aerocapture Trajectories
![Venus Aerocapture Trajectory and Heating](https://raw.githubusercontent.com/athulpg007/AMAT/master/plots/craig-lyne-altitude-higher-res-lower-size.png)
![Venus Aerocapture Trajectory and Heating](https://raw.githubusercontent.com/athulpg007/AMAT/master/plots/craig-lyne-heating-higher-res-lower-size.png)
### Neptune Aerocapture Feasibility Chart
![Neptune Aerocapture Feasibility](https://raw.githubusercontent.com/athulpg007/AMAT/master/plots/girijaSaikia2019b-higher-res-lower-size.png)
### Monte Carlo Simulations
![Monte Carlo Simulations](https://raw.githubusercontent.com/athulpg007/AMAT/master/plots/girijaSaikia2020b-fig-13-N5000.png)
![Monte Carlo Simulations](https://raw.githubusercontent.com/athulpg007/AMAT/master/plots/girijaSaikia2020b-apo-histogram-N5000.png)


## What kind of problems can AMAT solve?

AMAT can be used to quickly assess the feasibility of an aerocapture mission concept using aerocapture feasibiility charts, and perform trade studies involving vehicle type, control authority, thermal protection system materials, and useful delivered mass to orbit. AMAT can also be used to set up and run high-fidelity Monte Carlo simulations of aerocapture trajectories considering delivery errors, atmospheric uncertainties, and aerodynamic uncertainties to evaluate system performance under uncertainty.

### Examples

AMAT comes with a number of example [Jupyter notebooks](https://amat.readthedocs.io/en/master/jupyter_link.html) to help users get started. The examples inculde:

1. Venus aerocapture assesssment study
2. Titan aerocapture assessment study
3. Neptune aerocapture using blunt-body aeroshells
4. Drag modulation aerocapture asssessment
5. Planetary probe entry at various Solar System destinations

## Features

* Easy installation, only standard dependencies
* Many examples and recreation of results from published literature
* Extensive documented source code
* Can be used for probe and lander Entry, Descent, Landing (EDL) calculations
* Comes with a standard set of nominal atmospheric models

## Installation 

Note: AMAT is designed to work with Python 3.0 or greater. You must have a Python 3 installation in your system.

There are three ways to install AMAT. 

### Option 1 : Install from pip (recommended)

Note: Python Package Index limits the amount of additional data that can be packaged in the distribution, hence all data cannot be included in the built version. You will need to clone the GitHub repository to get the required data files, examples, and start using AMAT.

#### For Linux/MacOS machines:
  * ```$ pip install AMAT```
  * ```$ git clone https://github.com/athulpg007/AMAT.git```

If you are unable to clone the repository, you can download the repository as a .zip file from GitHub and extract it.

Once AMAT is installed, run an example Jupyter notebook to check everything works correctly.
  * ```$ cd AMAT/examples```
  * ```$ jupyter-notebook```

Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ```pip install jupyterlab pandas``` to install Jupyter and pandas if it is not already installed on your system. 


This will display the full list of example Jupyter notebooks included with AMAT.  Open and run the ```example-01-hello-world``` notebook to get started with AMAT.

#### For Windows machines:

You must have Anaconda installed. Open the Anaconda Prompt terminal:
  * ```$ pip install AMAT```

Open a Windows Powershell terminal and clone the GitHub reporistory. You must have Git installed.
  * ```$ git clone https://github.com/athulpg007/AMAT.git```

Run an example Jupyter notebook. From the Anaconda Prompt terminal:
  * ```$ cd AMAT/examples```
  * ```$ jupyter-notebook```

Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ```pip install jupyterlab pandas``` to install Jupyter and pandas if it is not already installed on your system. 


This will display the full list of example Jupyter notebooks included with AMAT. Open and run the ```example-01-hello-world``` notebook to get started with AMAT.


### Option 2 : Install from source

This will clone the repository from GitHub and install AMAT from the source code.

#### For Linux/MacOS machines:

Make sure you have numpy, scipy, matplotlib and pandas installed. Most likely you already have these installed. If not, use the following commands to install these dependenies first. Open a terminal window (on Linux machines) and type the following commands. You must have pip installed.
  * ```$ pip install numpy scipy matplotlib pandas jupyterlab```

Clone the GitHub repository and install AMAT.
  * ```$ git clone https://github.com/athulpg007/AMAT.git```
  * ```$ cd AMAT```
  * ```$ python setup.py install```
  * ```$ cd examples```
  * ```$ jupyter-notebook```

#### For Windows machines:

Open the Anaconda Prompt terminal to install the prerequisite packages.
  * ```$ pip install numpy scipy matplotlib```

Open a Windows Powershell terminal, clone the GitHub repository and install AMAT.
  * ```$ git clone https://github.com/athulpg007/AMAT.git```
  * ```$ cd AMAT```
  * ```$ python setup.py install```
  * ```$ cd examples```
  * ```$ jupyter-notebook```

Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ```pip install jupyterlab pandas``` to install Jupyter and pandas if it is not already installed on your system. 


To uninstall AMAT:

1. If you installed AMAT using pip:
  * ```$ pip uninstall AMAT```

2. If you installed AMAT from source, from the main AMAT directory:
  * ```$ python setup.py develop -u```

This will remove the AMAT installation from Python. You may simply delete the root folder where AMAT was installed to completely remove the files.


### Option 3 : Install in a virutalenv

If you plan to modifty the source code or add features, the recommended option is to to install it in a virtual environment. 

1. Change directory to where you want the virtual environment to be created.
  * ```$ cd home/path```

2. Create a virutal environment and activate it.

On Linux/MacOS machines:
  * ```$ python3 -m venv env1```
  * ```$ source env1/bin/activate```

On Windows machines (from Anaconda Prompt):
  * ```$ conda create --name env1```
  * ```$ conda activate env1```
  * ```$ conda install pip```

4. Follow the steps outlined in Option #2 (build from source) to clone the repository and install AMAT. If you make changes to the source code, remove the existing installation, update the setup file with a new version number, and re-install:
  * ```$ python setup.py develop -u```
  * ```$ python setup.py install```

If you want to create a new distrubution package:
  * ```$ python3 setup.py sdist bdist_wheel```

To re-make docs if you made changes to the source code (you must have Sphinx installed):
  * ```$ cd ~root/docs```
  * ```$ sphinx-apidoc -f -o source/ ../```
  * ```$ make html```

If you added a new AMAT module, appropriate changes must be made to docs/source/AMAT.rst.

## Usage

  * ```from AMAT.planet import Planet```
  * ```from AMAT.vehicle import Vehicle```
  * ```from AMAT.launcher import Launcher```

## License

AMAT is an open source project licensed under the GNU General Public License Version 3.

## Credits

Parts of the AMAT source code were originally developed in support of contracts between AAC and the Jet Propulsion Laboratory for various aerocapture mission studies between 2016 and 2020. Samples of atmospheric data from Global Reference Atmospheric Model (GRAM) software is used for illustration purpose only, and was developed by NASA Marshall Space Flight Center. The use of these GRAM models does not imply endorsement by NASA in any way whatsoever. A minimal working set of atmospheric profiles is included with AMAT to run the example notebooks. A minimal working interplanetary trajctory dataset is included with AMAT. The dataset was generated at Purdue University using the STOUR software package by Alec Mudek, and is also derived from trajectories published in the NASA Ice Giants Pre-Decadal Mission Study. The author plans to augment the interplanetary dataset with more publicly available information as it becomes available.

## Support and Contribution

If you wish to contribute or report an issue, feel free to [contact me](mailto:athulpg007@gmail.com) or to use the [issue tracker](https://github.com/athulpg007/AMAT/issues) and [pull requests](https://github.com/athulpg007/AMAT/pulls) from the [code repository](https://github.com/athulpg007/AMAT).

If you wish to make a contribution, you can do as follows:

 * fork the GitHub repository 
 * create a feature branch from *master* 
 * add your feature and document it
 * run the tests 
 * open a pull request

## Extras

The AMAT repository includes representative atmospheric profiles for Solar System bodies, an Excel sheet with a comprehensive literature review of aerocapture, sample feasibility charts for aerocapture at all destinations, reference journal articles (by the author), a PDF version of the documentation, and the author's Ph.D. dissertation which provides more details on the methods and algorithms implemented in AMAT.

## Reference Articles

Results from these articles are used as benchmark examples.

1. Craig, Scott, and James Evans Lyne. "Parametric Study of Aerocapture for Missions to Venus." Journal of Spacecraft and Rockets Vol. 42, No. 6, pp. 1035-1038. [DOI: 10.2514/1.2589](https://arc.aiaa.org/doi/10.2514/1.2589)

2. Putnam and Braun, "Drag-Modulation Flight-Control System Options for Planetary Aerocapture", Journal of Spacecraft and Rockets, Vol. 51, No. 1, 2014. [DOI:10.2514/1.A32589](https://arc.aiaa.org/doi/10.2514/1.A32589)

2. Lu, Ye, and Sarag J. Saikia. "Feasibility Assessment of Aerocapture for Future Titan Orbiter Missions." Journal of Spacecraft and Rockets Vol. 55, No. 5, pp. 1125-1135. [DOI: 10.2514/1.A34121](https://arc.aiaa.org/doi/10.2514/1.A34121)

3. Girija, A. P., Lu, Y., & Saikia, S. J. "Feasibility and Mass-Benefit Analysis of Aerocapture for Missions to Venus". Journal of Spacecraft and Rockets, Vol. 57, No. 1, pp. 58-73. [DOI: 10.2514/1.A34529](https://arc.aiaa.org/doi/10.2514/1.A34529)

4. Girija, A. P. et al. "Feasibility and Performance Analysis of Neptune
Aerocapture Using Heritage Blunt-Body Aeroshells", Journal of Spacecraft and Rockets, Vol. 57, No. 6, pp. 1186-1203. [DOI: 10.2514/1.A34719](https://arc.aiaa.org/doi/full/10.2514/1.A34719)---
title: "AMAT: A Python package for rapid conceptual design of aerocapture and atmospheric Entry, Descent, and Landing (EDL) missions in a Jupyter environment"
tags:
    - Python
    - aerocapture
    - atmospheric entry
    - mission design
    - planetary probe
authors:
    - name: Athul P. Girija
      orcid: 0000-0002-9326-3293
      affiliation: "1"
    - name: Sarag J. Saikia
      affiliation: "1"
    - name: James M. Longuski
      affiliation: "1"
    - name: James A. Cutts
      affiliation: "2"
affiliations:
    - name: School of Aeronautics and Astronautics, Purdue University, West Lafayette, IN 47907, United States
      index: 1
    - name: Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA 91109, United States
      index: 2
date: 31 August 2021
bibliography: paper.bib
---

# Summary
The Aerocapture Mission Analysis Tool (**AMAT**) is a free and open-source Python package for rapid conceptual design of aerocapture and atmospheric Entry, Descent, and Landing (EDL) missions in a Jupyter environment. Compared to conventional propulsive insertion, aerocapture allows interplanetary spacecraft to accomplish a near-propellantless method of orbit insertion at any planetary destination with a significant atmosphere. While aerocapture has been studied for many decades, to our knowledge, there are no publicly available tools for rapid conceptual design of aerocapture missions. AMAT aims to fill this gap by providing scientists and mission planners with an interactive design tool to quickly evaluate aerocapture mission concepts, perform parametric trade space studies, and select baseline mission concepts for higher-fidelity analysis. Atmospheric Entry, Descent, and Landing includes the hypersonic flight regime where the entry vehicle undergoes intense heating and rapid deceleration, followed by supersonic and subsonic flight, terminal descent phase, and final touchdown. AMAT provides capabilites for rapid end-to-end conceptual design of aerocapture and EDL mission concepts at any atmosphere-bearing Solar System destination.

# Statement of Need

Software tools for identification of trajectories and techniques that enhance or enable planetary missions or substantially reduce their cost are of great importance in realizing successful exploration missions [@Squyres2011]. While there have been some aerocapture mission analysis tools developed in the past such as ACAPS [@Leszczynski1998] and HyperPASS [@mcronald2003analysis], these tools are proprietary, outdated, and not publicly available. The NASA Ice Giants Pre-Decadal Mission Study in 2016, which studied mission concepts to Uranus and Neptune in the next decade highlighted the need for a design tool to evaluate aerocapture mission concepts [@Hofstadter2017; @Saikia2021]. A NASA-led study in 2016 to assess the readiness of aerocapture also underlined the need for an integrated aerocapture mission design tool which incorporates both interplanetary trajectory and vehicle design aspects [@Spilker2018]. While several software tools for conceptual design of EDL missions such as PESST [@otero2010planetary], TRAJ [@allen2005trajectory], and SAPE [@samareh2009multidisciplinary] have been developed in the past, these programs are not publicly available. AMAT aims to fill this gap for scientists, engineers, and academic researchers who want to quickly perform aerocapture mission analysis and atmospheric EDL as part of mission concept studies. 


# Rapid Design Capability

AMAT can be used to compute key aerocapture design parameters such as theoretical corridor width (TCW) and stagnation-point heat rate using just a few lines of code. Figure 1 shows an example of undershoot and overshoot trajectories for lift modulation aerocapture at Venus, and the corresponding stagnation-point heat rate computed by using AMAT.

![](https://i.imgur.com/3XPh6JY.png)
**Figure 1.** An example of a Venus aerocapture vehicle trajectory and the stagnation-point heating profile computed by using AMAT. Figure based on [@Craig2005].

AMAT can be used to create **aerocapture feasibility charts**, a graphical method to visualize the mission design trade space considering both interplanetary and vehicle design aspects [@Lu2018Titan]. Aerocapture feasibility charts help the mission designer assess the feasible range of vehicle lift-to-drag ratio ($L/D$) or ballistic-coefficient ratio ($\beta_2/\beta_1$) for lift and drag modulation respectively, as well as the feasible range of interplanetary arrival $V_{\infty}$ considering corridor width, deceleration, peak-heat rate, and total heat load. Figure 2 shows an example feasibility chart for lift modulation aerocapture at Neptune created by using AMAT.

![](https://i.imgur.com/BNINxh4.png)
**Figure 2.** Interplanetary trajectory trade space (left, generated using STOUR code) and aerocapture feasibility chart for Neptune (right, generated using AMAT). Figure based on the mission design presented in [@Saikia2021; @Girija2020b].  

AMAT can also be used to quickly set up and compute single-event jettison drag modulation aerocapture corridor bounds, and propagate guided trajectories. Figure 3 shows an example of trajectories for a drag modulation system at Mars computed by using AMAT.

![](https://i.imgur.com/YlMk6Th.png)
**Figure 3.** Altitude-velocity profiles and deceleration profiles for a drag modulation aerocapture system at Mars. Solid and dashed lines correspond to a 450 km circular and 1-sol target orbit respectively. Figure based on the flight system presented in [@Werner2019].  

AMAT can be used to perform Monte Carlo simulations to assess user-defined guidance schemes and system performance with a vehicle design and interplanetary trajectory considering navigation, atmospheric, and aerodynamic uncertainties as shown in Figure 4.

![](https://i.imgur.com/Jefki5T.png)
**Figure 4.** Examples of Monte Carlo simulation results showing post-aerocapture orbit parameters for a drag modulation system at Venus (left) and a lift modulation system at Neptune (right) computed by using AMAT. Based on [@Girija2020b].

AMAT can be used to generate EDL **carpet plots** which are commonly used to assess the trade space for preliminary entry system and mission design as shown in Figure 5.

![](https://i.imgur.com/uDxfzsS.png)
**Figure 5.** Examples of carpet plots for Venus entry (left), and Titan entry (right) created by using AMAT. Based on [@scott2018preliminary].

# AMAT Modules

The key functionality of the AMAT package is organized into three modules. These modules and their descriptions are shown in Table \ref{modules_table}. 

Table: Modules of AMAT \label{modules_table}

| Module        | Description                                                         |
| ------------- | --------------------------------------------------------------------|
| AMAT.planet   | Stores planetary constants, atmosphere models, entry interface etc. |
| AMAT.vehicle  | Stores vehicle parameters, guidance scheme, propogate trajectory    |
| AMAT.launcher | Stores Earth escape performance for a list of launch vehicles       |

# Documentation and Example Jupyter Notebooks

AMAT documentation along with a number of example Jupyter notebooks are available at https://amat.readthedocs.io. AMAT can be installed from the [Python Package Index](https://pypi.org/project/AMAT/) (`pip install AMAT`).

# Limitations

AMAT is intended to be used as a low-to-mid fidelity preliminary analysis tool to perform trade studies and select a baseline concept, which can be analyzed in further detail using high-fidelity tools such as DSENDS [@balaram2002dsends] or POST [@brauer1977capabilities].

AMAT uses publicly available empirical relations to compute the stagnation-point aerothermal environment. While this is sufficient for preliminary mission analysis, detailed studies require propietary higher-fidelity CFD tools such as LAURA [@mazaheri2013laura], DPLR and NEQAIR. While AMAT computes the stagnation-point total load which can be used to roughly estimate the thermal protection system (TPS) mass fraction [@Laub2004], there is no functionality to provide an accurate TPS mass estimate which is prerequisite when comparing an aerocapture system with propulsive insertion, especially for outer planet missions. Incorporation of improved TPS sizing models and addition of interactive visualization capabilities using Blender is planned for future work. 

# Acknowledgements

AMAT was developed at the Advanced Astrodynamics Concepts (AAC) group at Purdue University. Parts of the source code were originally developed in support of contracts between AAC and the Jet Propulsion Laboratory for various aerocapture mission studies. 

# References


tests package
=============

Submodules
----------

tests.test\_craigLyne2005 module
--------------------------------

.. automodule:: tests.test_craigLyne2005
   :members:
   :undoc-members:
   :show-inheritance:

tests.test\_girija2020 module
-----------------------------

.. automodule:: tests.test_girija2020
   :members:
   :undoc-members:
   :show-inheritance:

tests.test\_launcher module
---------------------------

.. automodule:: tests.test_launcher
   :members:
   :undoc-members:
   :show-inheritance:

tests.test\_planet module
-------------------------

.. automodule:: tests.test_planet
   :members:
   :undoc-members:
   :show-inheritance:

tests.test\_vehicle module
--------------------------

.. automodule:: tests.test_vehicle
   :members:
   :undoc-members:
   :show-inheritance:

tests.test\_wernerBraun2019 module
----------------------------------

.. automodule:: tests.test_wernerBraun2019
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: tests
   :members:
   :undoc-members:
   :show-inheritance:
AMAT package
============

Submodules
----------

AMAT.launcher module
--------------------

.. automodule:: AMAT.launcher
   :members:
   :undoc-members:
   :show-inheritance:

AMAT.planet module
------------------

.. automodule:: AMAT.planet
   :members:
   :undoc-members:
   :show-inheritance:

AMAT.vehicle module
-------------------

.. automodule:: AMAT.vehicle
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: AMAT
   :members:
   :undoc-members:
   :show-inheritance:
Example Jupyter Notebooks
==========================

.. toctree::
   :maxdepth: 1


   /examples/example-01-hello-world.ipynb
   /examples/example-02-atmosphere.ipynb
   /examples/example-03-venus-aerocapture-1.ipynb
   /examples/example-04-venus-aerocapture-2.ipynb
   /examples/example-05-titan-aerocapture-1.ipynb
   /examples/example-06-titan-aerocapture-2.ipynb
   /examples/example-07-venus-aerocapture-3.ipynb
   /examples/example-08-venus-aerocapture-4.ipynb
   /examples/example-09-uranus-aerocapture.ipynb
   /examples/example-10-neptune-aerocapture-part-1.ipynb
   /examples/example-11-neptune-aerocapture-part-2a.ipynb
   /examples/example-12-neptune-aerocapture-part-2b.ipynb
   /examples/example-13-venus-aerocapture-5a.ipynb
   /examples/example-14-fire-II-earth.ipynb
   /examples/example-15-apollo-as-201-earth.ipynb
   /examples/example-16-apollo-as-202-earth.ipynb
   /examples/example-17-apollo-4-earth.ipynb
   /examples/example-18-apollo-6-earth.ipynb
   /examples/example-19-paet-earth.ipynb
   /examples/example-20-viking-1-mars.ipynb
   /examples/example-21-viking-2-mars.ipynb
   /examples/example-22-pv-small-north-venus.ipynb
   /examples/example-23-pv-small-night-venus.ipynb
   /examples/example-24-pv-small-day-venus.ipynb
   /examples/example-25-pv-large-venus.ipynb
   /examples/example-26-galileo-jupiter.ipynb
   /examples/example-27-orex-earth.ipynb
   /examples/example-28-pathfinder-mars.ipynb
   /examples/example-29-mirka-earth.ipynb
   /examples/example-30-huygens-titan.ipynb
   /examples/example-31-atm-reentry-demonstrator-earth.ipynb
   /examples/example-32-deep-space-2-mars.ipynb
   /examples/example-33-stardust-earth.ipynb
   /examples/example-34-genesis-earth.ipynb
   /examples/example-35-hayabusa-earth.ipynb
   /examples/example-36-beagle-2-mars.ipynb
   /examples/example-37-opportunity-mars.ipynb
   /examples/example-38-curiosity-mars.ipynb
   /examples/example-39-oceanus-probe-saturn.ipynb
   /examples/example-40-oceanus-probe-uranus.ipynb
   /examples/example-41-hera-probe-saturn.ipynb
   /examples/example-42-dragonfly-titan.ipynb
   /examples/example-43-ice-giant-pre-decadal-uranus-probe.ipynb
   /examples/example-44-ice-giant-pre-decadal-neptune-probe.ipynb
   /examples/example-45-adept-vital-venus.ipynb
   /examples/example-46-uranus-entry-probe-decadal-2010.ipynb
   /examples/example-47-venus-flagship-entry-vehicle-decadal-2010.ipynb
   /examples/example-48-crew-module-atmospheric-re-entry-experiment-earth.ipynb
   /examples/example-49-earth-smallsat-aerocapture-demonstration-part-1.ipynb
   /examples/example-50-earth-smallsat-aerocapture-demonstration-part-2.ipynb
   /examples/example-51-mars-smallsat-aerocapture-demonstration-part-1.ipynb
   /examples/example-52-mars-smallsat-aerocapture-demonstration-part-2.ipynb
   /examples/example-53-mars-smallsat-aerocapture-demonstration-part-3.ipynb
   /examples/example-54-exomars-2016-mars.ipynb
   /examples/example-55-titan-aerocapture-systems-study-part-1.ipynb
   /examples/example-56-venus-entry-tradespace.ipynb
   /examples/example-57-titan-entry-tradespace.ipynbAPI Reference
=============

.. toctree::
   :maxdepth: 4

.. automodule:: AMAT.planet
   :members:
.. automodule:: AMAT.vehicle
   :members:
.. automodule:: AMAT.launcher
   :members:Credits
===========
AMAT was developed at the Advanced Astodynamics Concepts (AAC) research group at Purdue University. Parts of the AMAT source code were originally developed in support of contracts between AAC and the Jet Propulsion Laboratory for various aerocapture mission studies between 2016 and 2020. Samples of atmospheric data from Global Reference Atmospheric Model (GRAM) software is used for illustration purpose only, and was developed by NASA Marshall Space Flight Center. The use of these GRAM models does not imply endorsement by NASA in any way whatsoever. A minimal working set of atmospheric profiles is included with AMAT to run the example notebooks. A minimal working interplanetary trajctory dataset is included with AMAT. The dataset was generated at Purdue University using the STOUR software package by Alec Mudek, and is also derived from trajectories published in the NASA Ice Giants Pre-Decadal Mission Study. The author plans to augment the interplanetary dataset with more publicly available information as it becomes available.

Extras
-------

The AMAT repository includes representative atmospheric profiles for Solar System bodies, an Excel sheet with a comprehensive literature review of aerocapture, sample feasibility charts for aerocapture at all destinations, reference journal articles (by the author), a PDF version of the documentation, and the author's Ph.D. dissertation which provides more details on the methods and algorithms implemented in AMAT.
setup module
============

.. automodule:: setup
   :members:
   :undoc-members:
   :show-inheritance:
JSR Article Notebooks
==========================

Jupyter Notebooks associated with Girija et al. (2022), Quantitative Assessment of Aerocapture and Applications to Future Solar System Exploration, Journal of Spacecraft and Rockets. 

.. toctree::
   :maxdepth: 1

   /jsr-notebooks/01-literature-survey.ipynb
   /jsr-notebooks/02-interplanetary.ipynb
   /jsr-notebooks/03-a-venus-feasibility-lift.ipynb 
   /jsr-notebooks/03-b-venus-feasibility-drag.ipynb
   /jsr-notebooks/04-a-earth-feasibility-lift.ipynb 
   /jsr-notebooks/04-b-earth-feasibility-drag.ipynb
   /jsr-notebooks/05-a-mars-feasibility-lift.ipynb 
   /jsr-notebooks/05-b-mars-feasibility-drag.ipynb
   /jsr-notebooks/06-a-jupiter-feasibility-lift.ipynb 
   /jsr-notebooks/06-b-jupiter-feasibility-drag.ipynb
   /jsr-notebooks/07-a-saturn-feasibility-lift.ipynb 
   /jsr-notebooks/07-b-saturn-feasibility-drag.ipynb
   /jsr-notebooks/08-a-titan-feasibility-lift.ipynb 
   /jsr-notebooks/08-b-titan-feasibility-drag.ipynb
   /jsr-notebooks/09-a-uranus-feasibility-lift.ipynb 
   /jsr-notebooks/09-b-uranus-feasibility-drag.ipynb
   /jsr-notebooks/10-a-neptune-feasibility-lift.ipynb 
   /jsr-notebooks/10-b-neptune-feasibility-drag.ipynb
   /jsr-notebooks/11-mass-benefit-analysis-venus.ipynb
   /jsr-notebooks/12-mass-benefit-analysis-earth.ipynb
   /jsr-notebooks/13-mass-benefit-analysis-mars.ipynb
   /jsr-notebooks/14-mass-benefit-analysis-titan.ipynb
   /jsr-notebooks/15-mass-benefit-analysis-uranus.ipynb
   /jsr-notebooks/16-mass-benefit-analysis-neptune.ipynb
   /jsr-notebooks/17-comparative-studies.ipynb
   




make\_docs module
=================

.. automodule:: make_docs
   :members:
   :undoc-members:
   :show-inheritance:
Capabilities
=============


AMAT allows the user to perform low-fidelity broad sweep parametric studies; as well as high fidelity Monte Carlo simulations to quantify aerocapture performance. AMAT supports analysis for all atmosphere-bearing destinations in the Solar System: Venus, Earth, Mars, Jupiter, Saturn, Titan, Uranus, and Neptune.

Sample results
-----------------

Venus Aerocapture Trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: _images/craig-lyne-2005-higher-res.png
	:width: 400px
	:alt: Alitude vs time for an aerocapture trajectory at Venus
	:align: center


Neptune Aerocapture Feasibility Chart
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: _images/neptune-feasibility.png
	:width: 600px
	:alt: Neptune aerocapture feasibility chart
	:align: center

Monte Carlo simulations - Neptune aerocapture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: _images/prograde-higher-res.png
	:width: 400px
	:alt: Neptune aerocapture Monte Carlo simulation results
	:align: center

Monte Carlo simulations - SmallSat aerocapture at Venus
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. image:: _images/austin-drag-modulation-N1000.png
	:width: 400px
	:alt: Venus SmallSat aerocapture Monte Carlo simulation results
	:align: center



What kind of problems can AMAT solve?
--------------------------------------

AMAT can be used to quickly assess the feasibility of an aerocapture mission concept using aerocapture feasibiility charts, and perform trade studies involving vehicle type, control authority, thermal protection system materials, and useful delivered mass to orbit. AMAT can also be used to set up and run high-fidelity Monte Carlo simulations of aerocapture trajectories considering delivery errors, atmospheric uncertainties, and aerodynamic uncertainties to evaluate system performance under uncertainty.
Installation
=============

Note: AMAT is designed to work with Python 3.0 or greater. You must have a Python 3 installation in your system. AMAT has been tested on Linux/MacOS and Windows machines.

There are three ways to install AMAT. 

Option 1 : Install from pip (recommended)
----------------------------------------------

Note: Python Package Index limits the amount of additional data that can be packaged in the distribution, hence all data cannot be included in the built version. You will need to clone the GitHub repository to get the required data files, examples, and start using AMAT.

**For Linux/MacOS machines:**

  * ``$ pip install AMAT``
  * ``$ git clone https://github.com/athulpg007/AMAT.git``

If you are unable to clone the repository, you can download the repository as a .zip file from GitHub and extract it. Once AMAT is installed, run an example Jupyter notebook to check everything works correctly.

  * ``$ cd AMAT/examples``
  * ``$ jupyter-notebook``

.. note::
   Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ``pip install jupyterlab pandas`` to install Jupyter and pandas if it is not already installed on your system. 

This will display the full list of example Jupyter notebooks included with AMAT.  Open and run the ``example-01-hello-world`` notebook to get started with AMAT.

**For Windows machines:**

  * ``$ pip install AMAT``

(You must have Anaconda installed. Use the pip command from the Anaconda Prompt terminal). Open a Windows Powershell terminal and clone the GitHub repository. You must have Git installed.

  * ``$ git clone https://github.com/athulpg007/AMAT.git``

Run an example Jupyter notebook. From the Anaconda Prompt terminal:

  * ``$ cd AMAT/examples``
  * ``$ jupyter-notebook``

.. note::
   Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ``pip install jupyterlab pandas`` to install Jupyter and pandas if it is not already installed on your system. 


This will display the full list of example Jupyter notebooks included with AMAT. Open and run the ``example-01-hello-world`` notebook to get started with AMAT.


Option 2 : Install from source
-----------------------------------------------

This will clone the repository from GitHub and install AMAT from the source code.

**For Linux/MacOS machines:**

  * ``$ pip install numpy scipy matplotlib``

Clone the GitHub repository and install AMAT.

  * ``$ git clone https://github.com/athulpg007/AMAT.git``
  * ``$ cd AMAT``
  * ``$ python setup.py install``
  * ``$ cd examples``
  * ``$ jupyter-notebook``

.. note::
   Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ``pip install jupyterlab pandas`` to install Jupyter and pandas if it is not already installed on your system. 


**For Windows machines (from the Anaconda Prompt terminal):**

  * ``$ pip install numpy scipy matplotlib``

Open a Windows Powershell terminal, clone the GitHub repository and install AMAT.

  * ``$ git clone https://github.com/athulpg007/AMAT.git``
  * ``$ cd AMAT``
  * ``$ python setup.py install``
  * ``$ cd examples``
  * ``$ jupyter-notebook``

.. note::
   Note that you will need jupyterlab and pandas (for some examples) to run the example notebooks. Use ``pip install jupyterlab pandas`` to install Jupyter and pandas if it is not already installed on your system. 


**To uninstall AMAT:**

1. If you installed AMAT using pip:
  * ``$ pip uninstall AMAT``

2. If you installed AMAT from source, from the main AMAT directory:
  * ``$ python setup.py develop -u``

This will remove the AMAT installation from Python. You may simply delete the root folder where AMAT was installed to completely remove the files.


Option 3 : Install in a virutalenv 
---------------------------------------------------------

If you plan to modifty the source code or add features, the recommended option is to install it in a virtual environment. 

1. Change directory to where you want the virtual environment to be created.
  * ``$ cd home/path``

2. Create a virutal environment and activate it.

**On Linux/MacOS machines:**

  * ``$ python3 -m venv env1``
  * ``$ source env1/bin/activate``

**On Windows machines (from Anaconda Prompt):**

  * ``$ conda create --name env1``
  * ``$ conda activate env1``
  * ``$ conda install pip``

4. Follow the steps outlined in Option #2 (build from source) to clone the repository and install AMAT. If you make changes to the source code, remove the existing installation, update the setup file with a new version number, and re-install:

  * ``$ python setup.py develop -u``
  * ``$ python setup.py install``

If you want to create a new distrubution package:

  * ``$ python3 setup.py sdist bdist_wheel``

To re-make docs if you made changes to the source code (you must have Sphinx installed):

  * ``$ cd ~root/docs``
  * ``$ sphinx-apidoc -f -o source/ ../``
  * ``$ make html``

If you added a new AMAT module, appropriate changes must be made to docs/source/AMAT.rst.

AMAT Usage
------------

  * ``from AMAT.planet import Planet``
  * ``from AMAT.vehicle import Vehicle``
  * ``from AMAT.launcher import Launcher``
AMAT-master-3
=============

.. toctree::
   :maxdepth: 4

   AMAT
   setup
   tests
Contributions
=============

If you wish to contribute or report an issue, feel free to e-mail me here_ or to open an issue_ or create a pull_ request here from the repository_.

If you wish to make a contribution, you can do as follows:

 * fork the GitHub repository 
 * create a feature branch from *master* 
 * add your feature and document it
 * run the tests 
 * open a pull request

.. _here: mailto:athulpg007@gmail.com
.. _issue: https://github.com/athulpg007/AMAT/issues
.. _pull: https://github.com/athulpg007/AMAT/pulls
.. _repository: https://github.com/athulpg007/AMAT.. AMAT documentation master file, created by
   sphinx-quickstart on Thu May  7 10:27:29 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Aerocapture Mission Analysis Tool (AMAT)
========================================

.. image:: _images/AMAT_logo_cropped.png

AMAT_ is designed to provide rapid mission analysis capability for aerocapture mission concepts to the planetary science community. 

.. _AMAT: https://github.com/athulpg007/AMAT

See Jupyter_ notebooks to get started or refer to examples_ in the GitHub repository.

.. _Jupyter: https://amat.readthedocs.io/en/latest/jupyter_link.html
.. _examples: https://github.com/athulpg007/AMAT/tree/master/examples

AMAT allows the user to peform low-fidelity broad sweep parametric studies; as well as high fidelity Monte Carlo simulations to quantify aerocapture performance. AMAT comes with a suite of interplanetary trajectories, planetary atmosphere models, aeroheating correlations, and guidance algorithms for rapid conceptual mission design. AMAT provides aerocapture and Entry, Descent, Landing (EDL) mission analysis capability for Venus, Earth, Mars, Jupiter, Saturn, Titan, Uranus, and Neptune. 

.. image:: _images/planets.png


For sub-routine documentation, see :ref:`modindex`


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   about
   installation
   capabilities
   jupyter_link
   jsr-notebooks
   api_reference
   contributions
   credits
   references
   module_index

Index
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Module Index
=============
:ref:`modindex`About AMAT
===========

Overview
--------

**AMAT** is an open source collection of Python subroutines for rapid conceptual design of aerocapture and atmospheric Entry, Descent, and Landing (EDL) missions in a Jupyter environment.

AMAT comes with a suite of tools to allow end-to-end conceptual design of aerocapture missions: launch vehicle performance calculator, database of interplanetary trajectories, atmosphere models, vehicle control techniques, and aeroheating models. AMAT supports analysis for all atmosphere-bearing Solar System destinations for both lift and drag modulation control techniques.

AMAT has been extensively used in various aerocapture mission studies at the `Advanced Astrodynamics Concepts`_ (AAC) research group at Purdue University in collaboration with the NASA Jet Propulsion Laboratory (JPL). 

.. _Advanced Astrodynamics Concepts: https://engineering.purdue.edu/AAC/

History
-------

The lack of a rapid mission design tool for aerocapture mission concepts was identified by the NASA `Ice Giants Pre Decadal`_ (IGPD) Study led by JPL in 2016. This meant there was no quick way of performing architectural level assessments without resorting to resource intensive, subsystem-level design exercises such as the NASA `Aerocapture Systems Analysis`_ Team studies in 2004. 

.. image:: _images/igpd.png


.. _Ice Giants Pre Decadal: https://www.lpi.usra.edu/icegiants/
.. _Aerocapture Systems Analysis: https://ntrs.nasa.gov/search.jsp?R=20040111217

A team of researchers at Purdue University (Saikia et al.) led the aerocapture assessment studies in support of IGPD. Graduate researchers have since then further developed and extended the methods and tools for other atmosphere-bearing Solar System destinations. The focus was on developing an integrated systems engineering framework to allow mission designers to quickly evaluate the feasibility and performance of aerocapture mission concepts. Ye Lu and Athul P. Girija from the AAC research group conceptualized the **aerocapture feasibility charts**, now a commonly used graphical method for aerocapture mission design. An earlier version of the feasibility charts was presented by Saikia et al. in the IGPD study report.

Athul P. Girija  formulated a systems framework for rapid conceptual design of aerocapture missions for his doctoral thesis. Much of the AMAT source code was originally written in support of his Ph.D. dissertation work. AMAT was first publicly released in November 2019, and has since then been maintained by the author at Purdue University. In the spirit of `open code for open science`_, AMAT is free and open-source to foster universal access to the knowledge, and allow reproducibility of results by other researchers. Sugestions for improvement and potential contributions are greatly welcome.


.. _open code for open science: https://www.cos.io/about/mission

Related software
----------------

There are industrial software tools which offer mission analysis capabilities for aerocapture missions. These offer much higher fidelity, but are also substantially more complex to set up and run. Such fidelity is most often not required at the level conceptual studies. These tools are not in the public domain, and is generally available only to U.S. government affiliated labs, and U.S. persons at academic or industrial institutions. 

* `POST2`_: According to NASA Langley Research Center, "The Program to Optimize Simulated Trajectories II (POST2) is a generalized point mass, discrete parameter targeting and optimization program. POST2 provides the capability to target and optimize point mass trajectories for multiple powered or un-powered vehicles near an arbitrary rotating, oblate planet. POST2 has been used successfully to solve a wide variety of atmospheric ascent and entry problems, as well as exo-atmospheric orbital transfer problems."

* `DSENDS`_: According to NASA Jet Propulsion Laboratory, "DSENDS is a high-fidelity spacecraft simulator for Entry, Descent and Landing (EDL) on planetary and small-bodies. DSENDS (Dynamics Simulator for Entry, Descent and Surface landing) is an EDL-specific extension of a JPL multi-mission simulation toolkit Darts/Dshell which is capable of modeling spacecraft dynamics, devices, and subsystems, and is in use by interplanetary and science-craft missions such as Cassini, Galileo, SIM, and Starlight. DSENDS is currently in use by the JPL Mars Science Laboratory project to provide a high-fidelity testbed for the test of precision landing and hazard avoidance functions for future Mars missions. "


.. _POST2: https://post2.larc.nasa.gov/
.. _DSENDS: https://dartslab.jpl.nasa.gov/DSENDS/index.php


Future ideas
------------

Some things I would like to implement in the future:

* Pairing AMAT with `Blender`_ and `NASA 3D models`_ of planets and spacecraft to produce high resolution renders of aerocapture vehicle trajectories.

* Improved guidance schemes for lift and drag modulation aerocapture such as direct force control.

* Improved support for EDL mission concepts in the areas of precision landing, parachute dynamics, terminal descent and landing phases.

.. _Blender: https://www.blender.org/
.. _NASA 3D models: https://solarsystem.nasa.gov/resources


Note from the original author
------------------------------

I am `Athul P Girija`_, an aerospace engineer with a passion for robotic exploration of our Solar System. At the time of writing, I am a Ph.D. candidate in the Advanced Astrodynamics Concepts (AAC) research group led by Prof. Sarag Saikia and Prof. James Longuski at Purdue University. My interests lie in the general areas of planetary science, mission design and concept formulation, open-source software, and scientific programming.

AMAT is available for general public use under the GNU GPLv3 license.

`Google Scholar`_

.. _Athul P Girija: https://www.linkedin.com/in/athulpg007/
.. _Google Scholar: https://scholar.google.com/citations?hl=en&user=XxLVDPEAAAAJ
References
============

Results from these articles are also used as benchmark examples.

1. Craig, Scott, and James Evans Lyne. Parametric Study of Aerocapture for Missions to Venus. *Journal of Spacecraft and Rockets*, Vol. 42, No. 6, pp. 1035-1038. craiglyne2005_

2. Putnam and Braun, Drag-Modulation Flight-Control System Options for Planetary Aerocapture, *Journal of Spacecraft and Rockets*, Vol. 51, No. 1, 2014. putnamBraun2014_

3. Lu, Ye, and Sarag J. Saikia. Feasibility Assessment of Aerocapture for Future Titan Orbiter Missions. *Journal of Spacecraft and Rockets*, Vol. 55, No. 5, pp. 1125-1135. luSaikia2018_

4. Girija, A. P., Lu, Y., & Saikia, S. J. Feasibility and Mass-Benefit Analysis of Aerocapture for Missions to Venus. *Journal of Spacecraft and Rockets*, Vol. 57, No. 1, pp. 58-73. girijaSaikia2020a_

5. Girija, A. P. et al. Feasibility and Performance Analysis of Neptune Aerocapture Using Heritage Blunt-Body Aeroshells, *Journal of Spacecraft and Rockets*, Vol. 57, No. 6, pp. 1186-1203. girijaSaikia2020b_

.. _craiglyne2005: https://arc.aiaa.org/doi/10.2514/1.2589
.. _putnamBraun2014: https://arc.aiaa.org/doi/10.2514/1.A32589
.. _luSaikia2018: https://arc.aiaa.org/doi/10.2514/1.A34121
.. _girijaSaikia2020a: https://arc.aiaa.org/doi/10.2514/1.A34529
.. _girijaSaikia2020b: https://arc.aiaa.org/doi/10.2514/1.A34719
