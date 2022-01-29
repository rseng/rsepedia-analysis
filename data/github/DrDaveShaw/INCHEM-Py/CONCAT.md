---
title: 'INCHEM-Py: An open source Python box model for indoor air chemistry'
tags:
  - Python
  - atmospheric chemistry
  - indoor air
  - air pollution
  - box model
authors:
  - name: David Shaw^[corresponding author]
    orcid: 0000-0001-5542-0334
    affiliation: 1
  - name: Nicola Carslaw
    orcid: 0000-0002-5290-4779
    affiliation: 1
affiliations:
 - name: Department of Environment and Geography, University of York, Wentworth Way, York, YO10 5NG, United Kingdom
   index: 1
date: 12 March 2021
bibliography: paper.bib
---

# Summary

In developed countries people spend over 90% of their time indoors where they are exposed to airborne pollutants that are generated indoors or outdoors, some of which are harmful to human health. Occupant activities indoors such as cooking and cleaning can generate numerous chemical compounds and some of these can undergo further reactions to produce a large range of complex chemicals. Current instrumental techniques are unable to measure many of these compounds at present, so models provide the means to try and understand these processes and their impacts. The INdoor CHEMical model in Python (INCHEM-Py) is an open source and accessible box-model that has been re-factored from the indoor detailed chemical model developed by @Carslaw2007 to give researchers deeper insight into the chemical mechanisms behind indoor air chemistry. 

# Statement of need

Over the last 15 years, the INdoor Detailed Chemical Model (INDCM) has been used to successfully probe the chemistry of indoor air [@Carslaw2007]. However, it relies on proprietary software and requires a high level of chemistry expertise and some Fortran knowledge to edit and use the code. Software tools such as Pybox [@Topping2018], PyCHAM [@OMeara2020] and AtChem2 [@Sommariva2020] facilitate the use of chemical mechanisms to model atmospheric chemistry, but with a focus on chamber studies or ambient conditions. INCHEM-Py has been designed with a unique set of tools for the specific purpose of modelling indoor air chemistry. As well as a detailed gas-phase chemical mechanism, the new model includes gas-to-particle partitioning for three of the commonly encountered terpenes indoors (limonene and alpha- and beta-pinene), novel indoor photolysis parameterisation, indoor-outdoor air exchange and deposition to internal surfaces. INCHEM-Py is open source, has no black box processes and all inputs can be tracked through the model, allowing for complete understanding of the system. It has been designed to be easy to install for use by academics and students of all abilities, and is sufficiently accessible for further development by the wider indoor air community. The functionality embedded within INCHEM-Py will allow for a wide range of uses including in-depth analysis of experimental measurements, development and testing of new chemical mechanisms and probing numerous indoor scenarios, with the impacts on simulated indoor air pollutant concentrations from variations in parameters such as photolysis, ventilation and deposition rates, outdoor pollutant concentrations, time of year, and building location. 

# INCHEM-Py

INCHEM-Py creates and solves a system of coupled Ordinary Differential Equations (ODEs) to progress indoor atmospheric chemical species concentrations through time. It can be used to investigate numerous indoor air chemistry events in world-wide locales providing key insights into the chemical processes involved. INCHEM-Py utilises the Master Chemical Mechanism (MCM) at its core [@Jenkin1997; @Saunders2003], which is near explicit and contains no lumping or use of species surrogates. 

ICHEM-Py does not solve for spatial dimensions and assumes a well mixed atmosphere [@Carslaw2007]. Gas-to-particle partitioning is implemented using absorptive partitioning from @Pankow1994 with relevant reactions between particle and gas phase species included within the mechanism [@Carslaw2012]. Surface deposition is treated as an irreversible process with rates dependant on the simulated surface area to volume ratio and species deposition velocities [@Carslaw2012]. Photolysis rates can be calculated for several indoor lighting sources with different spectral profiles and also for attenuation of sunlight through multiple glass compositions [@Wang2021]. Daylight hours are derived from latitude and time of year.

The additional chemical mechanisms developed specifically for indoor air scenarios contained within the model have been utilised in a numerous published studies. These include: indoor air chemistry following cleaning with terpene based mixtures [@Carslaw2012; @Carslaw2017; @Carslaw2013; @Terry2014]; indoor air chemistry following cleaning with chlorine containing bleach [@Wong2017]; the impact of outdoor vegetation on indoor air chemistry [@Carslaw2015]; the importance of surface interactions for secondary pollutant formation [@Kruza2017]; ranking of harmful volatile organic compounds (VOCs) [@Carslaw2019]; improved model representation of the formation and composition of aerosols [@Kruza2020]; and photolysis of indoor air chemistry following high-concentration hospital/industrial cleaning events [@Wang2020]. INCHEM-Py has already been used to determine production rates and reactivity of indoor radical species, to assess the spatial and temporal scales of variability for indoor air constituents [@Lakey2021], and is currently being used to probe the impact of indoor air chemistry on ambient air, as well as to compare the differential secondary pollutant formation potential for different cleaning formulations.

At publication the current stable release of INCHEM-Py is v1.1.  

# Acknowledgements

The development of this model has been funded by a grant from the Alfred P. Sloan Foundation, grant number 2018-10083. Conclusions reached or positions taken by researchers or other grantees represent the views of the grantees themselves and not those of the Alfred P. Sloan Foundation or its trustees, officers, or staff.

# References
![INCHEM-Py](https://raw.githubusercontent.com/DrDaveShaw/INCHEM-Py/main/INCHEMPY_logo.png "INCHEM-Py logo")

[![status](https://joss.theoj.org/papers/674233ae20533ddbc12e9d89eec7b3bb/status.svg)](https://joss.theoj.org/papers/674233ae20533ddbc12e9d89eec7b3bb)


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5036374.svg)](https://doi.org/10.5281/zenodo.5036374)


# INCHEM-Py

The INdoor CHEMical model in Python (INCHEM-Py) is an open source box model that creates and solves a system of coupled Ordinary Differential Equations (ODEs) to provide predicted concentrations of indoor air pollutants through time. It is a refactor of the indoor detailed chemical model originally developed by [Carslaw (2007)](https://doi.org/10.1016/j.atmosenv.2006.09.038) with improvements in form, function, and accessibility. 

INCHEM-Py takes the [Master Chemical Mechanism (MCM)](http://mcm.leeds.ac.uk/MCM/), a near explicit mechanism developed for atmospheric chemistry, and additional chemical mechanisms that have been developed specifically for indoor air. These include gas-to-particle partitioning for three of the commonly encountered terpenes indoors (limonene and alpha- and beta-pinene), improved photolysis parametrisation, indoor-outdoor air exchange, and deposition to surfaces.

Typical usage of INCHEM-Py is either alongside experiment, where it can be used to gain a deeper insight into the chemistry through its ability to track a vast array of species concentrations or as a stand-alone method of investigating chemical events that occur indoors over a range of conditions. INCHEM-Py is open source, has no black box processes, and all inputs can be tracked through the model allowing for complete understanding of the system.  

A wide array of outputs from the model can be accessed, including species concentrations, species reactivity and production rates, photolysis values, and summations such as the total peroxy radical concentration. Custom reactions and summations can also be added by users to tailor the model to specific indoor scenarios.

A full pdf manual is included within this repository, including a quick start guide, this README is not intended to cover every aspect of INCHEM-Py but should be sufficient for users to get started.

A copy of the GNU General Public License v3.0 under which this project is licenced is included in this repository. 

# Table of contents
1. [Installation](#Installation)
2. [Running the model](#Running-the-model)
	* [Spyder](#Spyder)
	* [Command line](#cmd-line)
3. [Checking model function](#Tests)
4. [Community](#Community)

## Installation<a name="Installation"></a>

If you are new to Python then the "Quick Start Guide" in the INCHEM-Py manual will take you through a step-by-step guide of how to install and run the model.

INCHEM-Py has been developed in Python 3 using the [Anaconda environment](https://www.anaconda.com/products/individual) and as such we recommend its use for installing Python and running the model. The Anaconda environment contains all of the python packages required. The graphical installer from the above link should be used to download and install Anaconda 3 which will contain Python 3 and other useful tools.

If you have chosen not to use the Anaconda environment then any Python 3 installation can be used and the following non-default packages will need to be installed:
* numpy
* numba
* pandas
* tqdm
* scipy
* threadpoolctl
* matplotlib

Downloading the INCHEM-Py repository and extracting it to the location within your computer from which you would like to run the model will complete the installation. Multiple files are produced by INCHEM-Py during a model run. The output folder is named automatically using the current date and time in the format YYYYMMDD_hhmmss, and an additional custom title can be set within the settings.py file. You must have write permission to the directory in which INCHEM-Py is extracted. Further details are included in the user manual.

## Running the model<a name="Running-the-model"></a>

INCHEM-Py is run using the included settings.py script. Within this script are variables that can be adjusted and a full description of how to modify these is discussed in the user manual. The following instructions will be to run the model in its default downloaded state with no editing of files necessary. This process is also included with pictures in the quick start guide within the manual.

### Spyder<a name="Spyder"></a>

INCHEM-Py was written using the IDE (integrated developer environment) Spyder. Spyder can be installed both with the Anaconda install or from the Anaconda Navigator. Instructions on how to both install and run Spyder can be found [here](https://docs.anaconda.com/anaconda/user-guide/getting-started/)

Once Spyder is open you can set the INCHEM-Py directory as your working directory using the folder icon in the top right. Then by clicking the "files" tab at the bottom of that top right window the settings.py file, or any other you wish to edit, can be opened to the window on the left by double clicking on it.

Once you are ready the model can be run by opening the settings.py file in the left hand window and either clicking on the green arrow in the toolbar or by pressing F5 on your keyboard. The model will then run in the bottom right console window.

### Command line<a name="cmd-line"></a>

Assuming you have installed python as recommended via Anaconda it is possible to run INCHEM-Py from the Anaconda prompt from within the Anaconda Navigator on Windows, or from the terminal in MacOS or Linux. Once Anaconda prompt or the terminal is open simply navigate to the INCHEM-Py directory using the change directory command, inputting your file path

    cd C:/Directory/AnotherDirectory/INCHEM-Py/

and run the settings.py script with Python

    python settings.py

## Checking model function<a name="Tests"></a>

Included in the INCHEM-Py module folder is inchem_test.py. This script can be run to test functions of INCHEM-Py that manipulate the input data into useful formats within the model. It uses preset inputs, found within the test_files folder, to check that the model outputs are expected.

When entering new species or chemical mechanisms via the custom_input.txt file, users should be careful that names are correct and reactions are not duplicated. The model does not test for duplicate species or new species as both are valid inputs. Many outputs are produced by INCHEM-Py that can be used to check that the model ODEs are constructed as expected by the user or that custom chemistry has been entered correctly. The master_array can be viewed to validate reactions, and the Jacobian is saved for a similar purpose. Any user entered mechanisms are also saved alongside mechanisms provided by the INCHEM-Py team.

## Community<a name="Community"></a>

We welcome any contributions to INCHEM-Py and look forward to working with a community of people to develop the model further. 

Please contact us if you require any support.

david.shaw@york.ac.uk

nicola.carslaw@york.ac.uk

INCHEM-Py is free software, but since it is a scientific code we also ask that you show professional courtesy when using this code:
* Since you are benefiting from work on INCHEM-Py, we ask that you submit any improvements you make to the code to us by submitting a pull request. Issues should be reported using the issue tracker.
* If you use INCHEM-Py results in a paper or professional publication, we ask that it includes an appropriate reference. It is understood that in most cases if one or more of the INCHEM-Py team are involved in preparing results then they should appear as co-authors.
* The INCHEM-Py logo is included with the model and may optionally be used in any oral or poster presentations.
