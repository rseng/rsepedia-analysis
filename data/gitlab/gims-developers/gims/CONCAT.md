# CHANGELOG

All changes are tracked here: https://gims-developers.gitlab.io/gims/releaseNotes/changelog.html
# GIMS

GIMS (Graphical Interface for Materials Simulations) is a toolbox for electronic structure codes and supports the generation of input files and the analysis of output files. There are the following elemental apps:
* **Structure Builder:** Visualize, create or modify atomic structures.
* **Control Generator:** Set up the numerical Settings for the calculations.
* **Output Analyzer:** Analyze and visualize the output files.

These elemental apps are connected in two workflows:
* **Simple Calculation:** Prepare everything for a single-point calculation and later analyze the results.
* **Band Structure:** Based on your (periodic) structure the  k-point path is automatically generated and the input files are set up accordingly.

Currently, GIMS supports the following electronic structure codes:
* [FHI-aims](https://aimsclub.fhi-berlin.mpg.de)
* [Exciting](http://exciting-code.org)

## Documentation

For a complete documentation of the whole project please visit: https://gims-developers.gitlab.io/gims/index.html

## Requirements

Python:

| Package | Version |
|---------|---------|
| ase     | 3.20.1  |
| flask   | 1.1.2   |
| spglib  | 1.15.0  |

Javascript:

| Package             | Version  |
|---------------------|----------|
| chart.js            | 3.5.0    |
| parcel              | 2.0.0    |
| regenerator-runtime | 0.13.8   |
| three               | 0.130.1  |

## Quick start commands

To build GIMS two steps are needed: 
1. Build the client (javascript) in the `gims/client`
2. Install the python package in the `gims/app` folder.

After that the development server from flask can be started and you can run the app locally on your machine.

### How to build the client

First, we get into the application `client/` folder:

```
cd gims/client/
```
Please install all dependencies with:
```
npm install --production=false
```
We use the **watch** command from the parcel bunlder, so we watch the source code and rebuild when necessary. This is activated with:

```
npm run dev
```
If you just want to build the client once, please use:
```
npm run build
```

The result of this build process is copied to the `app/gims/static/` folder.

\* WARNING: Depending on the execution environment the *parcel* command parameter `--public-url` and the `Conf.BASE_FOLDER` variable (in the `Conf.js` code file) could need to be changed.


### How to run the application (development)

We use the [Flask](http://flask.pocoo.org/) Python application server. The project [User's guide](http://flask.pocoo.org/docs/1.0/) explains quite well how to set up the server and run a application. In short, if you have build the client go to:
```
cd gims/app/
```
Install the package with:
```
pip install .[dev]
```
and type:
```
export FLASK_APP=gims
flask run
```
Now, open the browser and enter the URL:
```
localhost:5000
```
You should now see the top-level GIMS interface.

### Using Docker

A Dockerfile is provided to build an image that can run a docker container
with both the server and client.
Run the following from the top-level of the repository:

```
docker build -t gims .
docker run -p 5000:5000 gims
```
# CONTRIBUTING

Please visit https://gims-developers.gitlab.io/gims/developerManual/contributing.html for details about how to contribute.
---
title: 'GIMS: Graphical Interface for Materials Simulations'
tags:
  - Python
  - JavaScript
  - Computational Materials Science
  - Electronic Structure Theory
  - Density Functional Theory
authors:
  - name: Sebastian Kokott
    orcid: 0000-0003-1066-6909
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Iker Hurtado
    orcid: 0000-0003-3805-4912
    affiliation: 1
  - name: Christian Vorwerk
    orcid: 0000-0002-2516-9553
    affiliation: 2
  - name: Claudia Draxl
    orcid: 0000-0003-3523-6657
    affiliation: 2
  - name: Volker Blum
    orcid: 0000-0001-8660-7230
    affiliation: 3
  - name: Matthias Scheffler
    orcid: 0000-0002-1280-9873
    affiliation: 1
affiliations:
  - name: The NOMAD Laboratory at the Fritz Haber Institute of the Max Planck Society, Berlin, Germany
    index: 1
  - name: Institut für Physik and IRIS Adlershof, Humboldt-Universität zu Berlin, Berlin, Germany
    index: 2
  - name: Department of Mechanical Engineering and Materials Science, Duke University, Durham, NC, United States of America
    index: 3
date: 17 December 2020
bibliography: paper.bib
---

# Abstract

GIMS (Graphical Interface for Materials Simulations) is an open-source browser-based toolbox for electronic-structure codes. It supports the generation of input files for first-principles electronic-structure calculations and workflows, as well as the automated analysis and visualization of the results. GIMS is deliberately extendable to enable support for any electronic-structure code. Presently, it supports two different software packages: the numerical atom centered orbital package `FHI-aims` and the LAPW code `exciting`.


# Statement of Need

Common workflows for electronic-structure calculations require at least the following steps *1. Generating input files:* This step includes the definition of structural data (e.g. position of atoms) and numerical settings (e.g. basis-set quality, runtime choices, and numerical convergence criteria). *2. Running the calculation(s):* Based on the input files, the electronic-structure code performs the requested calculation. Usually, each calculation produces several output files. *3. Post-Processing:* The output files are parsed, analyzed, and results are finally visualized. Step 1 to 3 can be repeated and connected to workflows.

While step 2 is usually run by an *ab initio* engine on a remote HPC cluster, steps 1 and 3 can be executed on local machines. The electronic-structure community primarily uses the *command line interface* and *text editors* as natural working environments, where workflows are automated using scripts. Each of the above-mentioned steps adds its own technical complexity and, thus, potential barriers for the user.

The objective of *GIMS* is to lower the entry barrier and to provide an easy-to-use, platform-independent, zero-setup toolbox for standard tasks within the framework of first-principles electronic-structure calculations. A running GIMS application can be found here: [https://gims.ms1p.org](https://gims.ms1p.org). GIMS is intentionally written and designed to be easily extendable to any electronic-structure code. At present, it supports the `FHI-aims` [@blum:2009] and `exciting` [@gulans:2014] codes.


# Software Architecture

The application is designed as web client-server system, but can be run entirely on a local machine. The client side is responsible for the user interaction, file parsing, and data visualization. The primary programming language is JavaScript. Conceptually, the web client is designed as a single page application: the client application is loaded at the outset and some of the data it displays is dynamically updated at runtime by the server. The server has no User Interaction (UI) logic nor does it maintain an UI state.

The server part is written in [python](https://www.python.org). The communication between web application and web server is realized by using the web server gateway interface (WSGI) framework [Flask](https://flask.palletsprojects.com/en/1.1.x/) [@flask:2010]. Client requests are interfaced with the [ASE package](https://wiki.fysik.dtu.dk/ase/) [@ase:2017] on the server side. The ASE package provides python objects for the code-independent handling of atomic structures (`Atoms object`) and numerical settings (`Calculator`). Thus, the use of ASE will enable to extend GIMS' functionalities to all codes (that is, by the time of writing more than 40 different codes) supported by ASE in a straightforward way. Moreover, we integrate [`spglib`](https://spglib.github.io/spglib/) [@spglib:2018] to obtain additional symmetry properties for periodic structures.

# Overview of Features

GIMS is structured in terms of three separate *elemental* apps that address step 1 (input generation) and step 3 (post-processing) described above. These apps also serve as building blocks for workflows (see below). The current three elemental apps are:

1. **Structure Builder.** This app allows to import, view, manipulate, taking snapshots of, and export structure files for various file formats. The 3D structure viewer is based on the [threejs](threejs.org) library [@threejs:2010]. The builder enables a user to add, delete, and change properties of atoms, as well as to analyze molecular and periodic structures (e.g. measuring distances between two atoms, angles between three atoms, getting symmetry information for periodic structures). This is a subset of capabilities as found, e.g., in existing visualization and building tools such as [Jmol](http://www.jmol.org/), [Avogadro](https://avogadro.cc) [@avogadro:2012], and [PyMol](https://pymol.org/2/), but intrinsically designed as part of a broader client-server framework in the case of GIMS.

2. **Control Generator.** Another step needed to set up a calculation is the selection of numerical parameters. This process is highly specific to each electronic-structure code. The control generator allows the user to provide the basic parameters for the selected electronic-structure code. The user can select items from a form, where tool-tip help provides code-specific information about the listed keywords. As a final product of this step, the input file for the selected code is created and available for download.

3. **Output Analyzer.** After running the calculation, analysis and post-processing of the output files are needed. The output analyzer facilitates some basic tasks, such as output file identification (that is, automatically identifying the code the output files came from and what kind of output files were provided), file parsing, visualization of the results and numerical convergence of the calculations. Graphs can be interactively modified and downloaded as png picture that can be directly used for a presentation or publication.

*Workflow apps* in GIMS combine different *elemental* apps. For instance, the *band-structure* workflow proceeds as follows: First, the user selects the electronic-structure code with which they want to carry out the corresponding band-structure calculation. Second, based on the provided periodic structure and underlying Bravais lattice defined in the *structure builder*, the correct band path is automatically determined (according to the Setyawan-Curtarolo convention [@setyawan:2010] as implemented in the [ASE package](https://wiki.fysik.dtu.dk/ase/)). Third, mandatory keywords to run a band structure calculation are pre-selected in the *control generator*. Fourth, the band path is incorporated into the corresponding input file. Finally, all resulting output files are processed and visualized by the *output analyzer* in the last step of the workflow.

Parsing and visualizing in- and output files in a browser based framework is also an integral part of other projects that make use of electronic structure data, such as the materials databases [NOMAD](https://nomad-lab.eu), [AFLOW](http://aflowlib.org), and [Materials Project](https://materialsproject.org). However, the GIMS workflow apps focuses also on the *generation* of input files. Workflow apps help to make the user aware of potential pitfalls, e.g., background sanity checks of the input files and cross checks among all input files (a simple example would be to make it mandatory to define a k-grid when a periodic structure is used).

The manual and a detailed description of all features are available at: [https://gims-developers.gitlab.io/gims](https://gims-developers.gitlab.io/gims). Both the client as well as server part are integration tested using [jest](https://jestjs.io) and [pytest](https://docs.pytest.org/en/latest/), respectively.

# Acknowledgements
This work received funding from the European Union’s Horizon 2020 Research and Innovation Programme (grant agreement No. 951786), the NOMAD CoE, MS1P e.V., and ERC:TEC1P (No. 740233).

# References
