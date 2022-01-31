Copyright
=========

fuelcell: A Python package and graphical user interface for electrochemical 
data analysis (fuelcell) Copyright (c) 2021, The Regents of the University 
of California, through Lawrence Berkeley National Laboratory (subject to 
receipt of any required approvals from the U.S. Dept. of Energy) and 
University of California, Berkeley. All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.
[![Documentation Status](https://readthedocs.org/projects/fuelcell/badge/?version=latest)](https://fuelcell.readthedocs.io/en/latest/?badge=latest) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/BSD) [![status](https://joss.theoj.org/papers/79df7852b6fedb417b625979da1f5567/status.svg)](https://joss.theoj.org/papers/79df7852b6fedb417b625979da1f5567) [![DOI](https://zenodo.org/badge/292444002.svg)](https://zenodo.org/badge/latestdoi/292444002)


# fuelcell
`fuelcell` is a Python package which streamlines electrochemical data analysis. 

`fuelcell` includes both a standard Python package which can be used programmatically and incorporated into an existing workflow, as well as a standalone graphical user interface for interactive use.

If you have a feature request or find a bug, please file an [issue](https://github.com/samaygarg/fuelcell/issues) or submit a [pull request](https://github.com/samaygarg/fuelcell/pulls). This is designed to be an open-source tool which the entire electrochemical community can build upon and use.

## Installation
`fuelcell` can be installed from PyPI using `pip`:

```bash
pip install fuelcell
```

## GUI 
The latest version of the standalone GUI can be downlaoded [here](https://fuelcell.readthedocs.io/en/latest/gui.html). There are separate versions for Mac and Windows operating systems.

## Documentation
The documentation can be found at [fuelcell.readthedocs.io](https://fuelcell.readthedocs.io/en/latest/) 

##  Package Structure
- `datums.py`: Data processing functions
- `visuals.py`: Data visualization functions
- `utils.py`: File handling and general auxilliary functions
- `model.py`:  `Datum` class to store electrochemical data along with associated features  and expereimental parameters
- `fuelcell_gui.py`: Graphical user interface for interactive use

## License
[MIT](https://choosealicense.com/licenses/mit/) 
---
title: '`fuelcell`: A Python package and graphical user interface for electrochemical data analysis'
tags:
  - Python
  - electrochemistry
  - chemical engineering
  - fuel cells
  - electrolysis
authors:
  - name: Samay Garg
    affiliation: "1, 2"
    orcid: 0000-0002-9872-6429
  - name: Julie C. Fornaciari
    affiliation: "1, 2"
    orcid: 0000-0002-0473-2298
  - name: Adam Z. Weber
    affiliation: 1
    orcid: 0000-0002-7749-1624
  - name: Nemanja Danilovic
    affiliation: 1
affiliations:
 - name: Energy Storage and Distributed Resources Division, Lawrence Berkeley National Laboratory, Berkeley, CA 94720
   index: 1
 - name: Department of Chemical and Biomolecular Engineering, University of California Berkeley, Berkeley CA 94720
   index: 2
date: 10 November 2020
bibliography: paper.bib
---

# Overview

`fuelcell` is a Python package designed to standardize and streamline the analysis of electrochemical data. This package includes modules for data processing and data visualization, as well as a graphical user interface (GUI) for interactive use.

# Introduction

As the demand for sustainable, carbon-free electricity increases globally, development of electrochemical energy conversion devices is increasing rapidly. These devices include fuel cells, flow batteries, and water electrolysis cells. A wide range of diagnostic experiments is used to assess the performance, durability, and efficiency of electrochemical devices. [@Bard:2001; @Newman:2004]. Among the most commonly used techniques are chronopotentiometry (CP), chronoamperometry (CA), cyclic voltammetry (CV), linear sweep voltammetry (LSV), and electrochemical impedance spectroscopy (EIS) experiments.[@Bard:2001; @Newman:2004; @Orazem:2008; @Wang:2003]. Although these experimental protocols have been well-established in the field of electrochemistry, the protocols for analyzing electrochemical data have not been clearly standardized. Standardizing electrochemical data analysis will also aid in applying machine learning frameworks to extract valuable information from electrochemical data sets. 

# Statement of need

A single electrochemical experiment can generate on the order of ten thousand data points, and several individual experiments are frequently used to assess a single cell. Electrochemical experiments also generate large quantities of raw data, which require extensive preprocessing before the data can be used to assess the performance of an electrochemical device completely. Processing and analyzing the data from a single experiment using conventional methods often is a bottleneck and time consuming. Manually processing this data also introduces unnecessary human error into the results, resulting in increased variation both between individual researchers and between research groups within the electrochemical field [@Agbo:2019]. Therefore, an application that efficiently processes electrochemical data will standardize and expedite the analysis of data generated from electrochemical experiments.

# Functionality

`fuelcell` includes modules for both data processing and visualization to enable a smooth, efficient workflow from raw data to publication-ready figures. These modules can be imported and used programmatically as a standard Python package, and `fuelcell` also includes a standalone GUI that allows users with little or no programming experience to utilize these modules.  `fuelcell` also serves as a platform that can be expanded to facilitate new and more advanced techniques as the needs of the electrochemical community evolve.

## `fuelcell.datums`
Every experiment requires a unique protocol to process the raw experimental data, so the datums module contains experiment-specific functions to read and process data for each experiment. Currently, functions to process CV, CP, CA, LSV, and EIS data are included fuelcell, and new protocols can be added to the project by opening an issue on GitHub. The complexity of these processing protocols varies depending on the experiment, and the specific data processing steps carried out for each experiment are detailed in the documentation. The datums module also includes functions to determine the Tafel slope and exchange current density from LSV data as well as to extract the high-frequency resistance (HFR) value from EIS data [@Orazem:2008; @Agbo:2019].

## `fuelcell.visuals`
The visuals module includes functions to generate visualizations, which are commonly used in the electrochemical community. `fuelcell` currently includes functions to generate polarization curves, cyclic voltammograms, linear sweep voltammograms, and Nyquist plots (Figure 1). This module is built around the matplotlib library, which allows for highly customizable visualizations. The visuals module is designed to integrate both seamlessly with the datums module and to function as an independent module that can be incorporated into an existing workflow.

![Examples of figures created using functions in `fuelcell.visuals`. (a) Polarization curves generated using data from CP experiments. (b) Cyclic voltammograms. (c) LSV data with the Tafel fit overlaid in yellow. (d) EIS data with the HFR value calculated using both a semicircle fit and a linear fit.\label{fig:1}](fig1.png)

## `fuelcell_gui`
The GUI is included in the standard `fuelcell` installation, but it can also be installed independently as a single executable file (Windows and MacOS) that includes all necessary dependencies. The GUI also enables users to interactively create and customize visualizations without being familiar with the ins and outs of the matplotlib library. This GUI has been shown to greatly reduce the time required to process electrochemical data, with researchers using the program reporting that it reduces the time required to process data from testing four cells from close to one hour to about five minutes.

![Data visualization tab of the GUI.\label{fig:2}](fig2.png)

# Acknowledgements
SG acknowledges funding from the Berkeley Lab Undergraduate Research Fellowship. The authors thank Dr. Xiong Peng, Zachary Taie, Eden Tzanetopoulos, and Grace Anderson for helpful discussions and assisting with testing the program.

# References
## Setup

There are two different versions of the program: One for Mac and one for Windows; download the version that will work on your operating system. The program should be able to run as-is without downloading anything else. It may take a minute or two to launch, especially the first time, but this is normal.

## Preprocessing

Very minimal preprocessing of the data is required to ensure the program works correctly. 

### File type

The program is currently able to read  delimited text, Microsoft Excel (xls and xlsx), and CSV files. 

### File naming

Ensure that the abbreviation for the type of experiment which was run appears somewhere in the filename. For example, if the data are from a chronopotentiometry experiment, 'cp' should be in the file name somewhere; if the data were from a linear sweep voltammetry experiment, 'lsv' should be in the file name somewhere; etc.

### File struture

#### Excel and CSV

All data for a particular test should be in a single spreadsheet.

If this text is in the first row of the spreadsheet, delete the first row.

![png](rowone.png) 

If you used EC lab to generate the data, leave the column headings as-is so the program can automatically identify the columns. Otherwise you will need to manually specify the columns. In this case, the column headings and the order of  the columns do not matter. It can be helpful, however, if you maintain your own convention for the order of the columns. Delete any empty columns which may be between columns containing data. Finally, there should only be one column for each variable. This is mainly an issue for CP and CA data, since you might have two time columns; just delete one of them. That's it! 

Example file:

<img src="xlsexample.png" alt="png" style="zoom:50%;" />

#### Text Files

The file should look similar to the example below; again, the order of the columns doesn't matter, and the column headings can remain as is. I've never had to modify text files in any way after I export them from EC Lab

Example file:

![png](txtexample.png)
Graphical User Interface
=========================

The GUI for the fuelcell package can be downloaded as a standalone executable file and can be used on either Mac or Windows operating systems. To use the latest version GUI, download the zip file which contains the GUI appropriate to your operating system from the `github releases page <https://github.com/samaygarg/fuelcell/releases>`_. If you are using MacOS, you will likely have to give permission for the program to run after trying to open in for the first time. This can be done by going to System Preferences -> Security & Privacy -> General and clicking "Open Anyway"
To report any problems that arise while using the GUI, open a `new issue <https://github.com/samaygarg/fuelcell/issues>`_ on github.

New in version 0.5.2
---------------------
* resolved deprecation issues to ensure forward compatibility

Setup
=======

Installation
--------------
The fuelcell library can be installed using `pip <https://pypi.org/project/fuelcell/>`_

.. code-block:: bash

   pip install fuelcell


Updating
---------
The fuelcell library is still under development, so supdates are being made and features are being added regularly. You should ensure that you are updated to the latest version. fuelcell can be updated using pip:

.. code-block:: bash

   pip install --upgrade fuelcell


Graphical User Interface
-------------------------
The graphical user interface (GUI) can be downloaded as a standalone executable file for either Windows or Mac operating systems. The latest version can be downloaded `here <https://github.com/samaygarg/fuelcell/tree/master/gui/fuelcell_gui>`_. See the GUI page for more details.fuelcell Documentation
=====================================

Introduction
-------------
fuelcell is a Python library designed to streamline the process of analyzing data generated from experiments which are commonly used to test electrochemical devices.


Contents
---------

.. toctree::
   :maxdepth: 1

   setup.rst
   gui.rst
   modules/modindex.rst
   
GUI
======
Download the latest version of the GUI `here <https://github.com/samaygarg/fuelcell/tree/master/gui/fuelcell_gui>`_.
fuelcell.visuals
==================

Data visualization methods

.. automodule:: fuelcell.visuals
	:members:fuelcell.utils
================

File handling functions and miscellaneous auxilliary functions

.. automodule:: fuelcell.utils
	:members:Modules
==========

Contents:

.. toctree::
	:maxdepth: 2
	
	datums.rst
	utils.rst
	visuals.rst
fuelcell.datums
=================

Data processing functions

.. automodule:: fuelcell.datums
   :members:
   :undoc-members: