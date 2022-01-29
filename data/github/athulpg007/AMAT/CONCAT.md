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


