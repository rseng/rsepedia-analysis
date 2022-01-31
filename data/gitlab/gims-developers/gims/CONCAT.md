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

Where to find GIMS
==================

Running Versions
----------------

Right now we maintain two servers running the GIMS application.

1. Stable Version
    The latest stable version runs on https://gims.ms1p.org. This server runs version ``1.0.9`` of GIMS.

2. Development Version
    We also try to offer to the user the latest version from the master branch. However, there is no automatic deployment, so there might be some delay involved. You can find the development version right here: https://gims-dev.ms1p.org. Currently, this server runs version ``1.0.9`` of GIMS.

GitLab
------

The GIMS project is hosted on Gitlab: https://gitlab.com/gims-developers/gims.
Feel free to join the group of developers and contribute your own features to GIMS.
.. GIMS documentation master file, created by
   sphinx-quickstart on Tue Feb  4 18:05:07 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Graphical Interface for Materials Simulations
=============================================

GIMS (Graphical Interface for Materials Simulations) is a browser-based toolbox for electronic structure codes and supports the generation of input files for first-principles electronic structure calculations and workflows, as well as the analysis and visualization of the resulting data extracted from the output files.

.. figure:: images/docsOverview.png
   :scale: 20 %


.. toctree::
    about
    whereToFind
    userManual/userManual
    developerManual/developerManual
    releaseNotes/releaseNotes

..
  Indices and tables
  ==================

  * :ref:`genindex`
  * :ref:`modindex`
  * :ref:`search`

About
=====

The objective of *GIMS* is to lower the entry barrier and to provide an easy-to-use, platform-independent, zero-setup toolbox for standard tasks within the framework of first-principles electronic-structure calculations. Information about running GIMS application can be found here: `Where to find GIMS <whereToFind.html>`_. GIMS is intentionally written and designed to be easily extendable to any electronic-structure code. At present, it supports the ``FHI-aims`` and ``exciting`` codes.

Feature List
------------

* Structure visualisation of various different formats: investigate, manipulate, export, and take snapshots of your structure.
* Select your needed input parameters from a list, which includes detailed explanations for all keywords.
* Visualize all of your output files (main output files, density of states, band structure + Brillouin Zone and Path) with ease. Publication ready graphs and pictures can be downloaded.
* Select specific workflows. The workflow will guide you through the needed steps and background checks will help to avoid easy setup errors for your calculation. 

Structure Builder
=================

The *Structure Builder* is an app for viewing and manipulating atomic structures, such as molecules and solids.

Import structure files
----------------------

In the right upper corner you find the *Import* button. The *Structure Builder* first tries to detect the file format, automatically. If it fails, you can choose among the following formats:

+------------------+------------------------------------------+
| Code             | File Detection                           |
+==================+==========================================+
| FHI-aims         | File ends with ``*.in``                  |
+------------------+------------------------------------------+
| Exciting         | File name ``input.xml``                  |
+------------------+------------------------------------------+
| CIF-file         | File ends with ``*.cif``                 |
+------------------+------------------------------------------+
| XYZ-file         | File ends with ``*.xyz``                 |
+------------------+------------------------------------------+
| VASP             | File name ``POSCAR``                     |
+------------------+------------------------------------------+
| quantum espresso | File contains ``&system`` or ``&SYSTEM`` |
+------------------+------------------------------------------+

3D Structure Viewer
-------------------

If the file has been imported successfully, the structure of atoms should be displayed. You can inspect the structure in the following ways:

* Rotate the structure: Hold left mouse and move
* Move the structure: Hold right mouse and move
* Get atom information: Hover over atom

Structure Analysis
""""""""""""""""""

Further, you can easily do some structure analysis:

* Measure atom distance: Hold shift and select two atoms
* Measure angle of three atoms: Hold shift and select three atoms
* Measure torsion angle of four atoms: Hold shift and select four atoms
* De-select atoms: Click right

In the left upper corner the button **View options** allows to further alter the displayed content, e.g. you can switch on and off the display of the unit cell, bonds, lattice vectors and so on.

Changing Species Colors
"""""""""""""""""""""""

The colors of the species and their respective names are shown on in the middle on the left side of the structure viewer. Clicking on the color symbol allows to select a new color for the species. To reset the colors of all elements click on the *circle arrow*.

Structure Manipulation
""""""""""""""""""""""

There are two ways to manipulate your uploaded structure: 1) 3D Structure Viewer 2) Side Panel.

    1. Structure Viewer: To edit the structure click on the button **Edit structure**. You can now perform the following actions:

        * Move atoms: Click on an atom and hold and move atom (drag and drop)
        * Delete atom: Press D and click on atom

    2. Side Panel: To the right  of the 3D structure viewer you find editable fields. You can:

        * Get structure information about your system, such as space group, lattice parameters, occupied Wyckoff positions, etc.
        * Change or delete lattice vectors. When changing lattice vectors you can choose either to scale all atom position with the lattice vectors. To activate this feature click on the checkbox below the lattice vector fields.
        * Create supercell. By filling the fields with integer number in the **Create Supercell** section you can create a supercell with lattice vectors multiples of the initial length.
        * Create the primitive cell, that is, find the smallest possible unit cell (implemented via spglib).
        * Change basis atoms: All properties of the basis atoms can be edited: atom position atom species. Moreover, by **clicking on the species color**, you can **constrain the atom positions** (e.g. for relaxation calculations) and you can set an **initial spin moment** on that atom (if supported by the corresponding code). Clicking on the trash bin deletes the atom from the structure. To add a new atom click on the **New atom** button.

        .. figure:: ../images/basisAtoms.png
           :scale: 40 %

           Click on the species color behind the element symbol field to constrain the element for a relaxation and/or to set the initial moment of the corresponding atom.

You can follow the **number of changes** of your structure in the top center of the 3D structure viewer (below the button **Edit structure**). Click on the arrows to the left and right to undo and redo your last actions, respectively.

3D Brillouin Zone Viewer
------------------------

In case of a periodic structure (structure has three lattice vectors) the corresponding Brillouin zone (BZ) and the band path according to Setjawan/Curtarolo is shown.


Export Structure Files
----------------------

To export your structure click on the **Export** button. The export format will be according your selected code, i.e. if you have selected FHI-aims, a geometry.in file will be downloaded to your download folder.

Please note: Due to browser security restrictions, we can only save files to the download folder. To change the download folder please change it in the settings of your browser.

Reset Structure
---------------

To erase all structural data click on the button "Reset" in the top right corner. Be careful: It will reset all changes and steps before. You will have to confirm this step.
Prerequisites
=============

All main features of GIMS are roughly tested for the following Browsers: Google Chrome, Safari (MacOS), Firefox.
However, some specific features, which we may have not tested explicitly, may not work for all browsers. GIMS is developed to primarily work for the Chrome and Chromium browser. Please consider to use them, if you experience any problems for your favorite browser. Feel free to report any issues in the GIMS repository: https://gitlab.com/gims-developers/gims

Output Analyzer
===============

The *output analyzer* app visualizes output files of implemented electronic structure codes.

Import Files
""""""""""""
Click on the button next to **Import output files** to upload your output files.
You can choose either uploading single files, batches of files or folders.
To upload folders mark the checkbox right to the upload button.

After selection the *output analyzer* figures out the used code and parses the files.
The result of the parsing is shown and organized in different sections.

System Information
""""""""""""""""""

This section gives details about the investigated system, that is, chemical formula and number of atoms. Right to this section the system is visualized. In case of a structural relaxation you can watch an animation, one picture per optimization step. To start the animation, click on the **play** symbol. The optimization progress is given left to the play button, where 0 refers to the input structure.

Results
"""""""

If the calculation has finished as expected, this section summarizes the most important characteristics, e.g., total energy, or information about electron levels.

If the calculation included a relaxation of the structure, the final geometry can be downloaded by clicking on "Download final geometry" in the last line of the Results section.

Summary Calculation
"""""""""""""""""""

Summarizes information about the code itself (e.g. code version, commit number), used computational resources and runtime information (e.g. number of tasks, memory calculation time).

The exit status "Calculation exited regularly" indicates, whether the parsed output file ended correctly and as expected. That means, that the calculation converged to a results, ended within the wall time, and did not crash. If the exit status "Calculation exited regularly" returns **no**, you should manually inspect the output file. Further, for the exit status **no**, no result should be shown, except for the following case: calculations that contain more than one SCF cycle, e.g. optimization of structures. For such a case, the results for the last converged SCF cycle is shown.

Within the results section convergence graphs appear. In case of a relaxation, properties per SCF cycle are visualized, e.g. the difference in total energy or the maximum force component. Moreover, convergence of quantities per SCF iteration of a single SCF cycle are shown as well (the only one in case of a single point calculation).

DOS and Band Structure Graphs, Brillouin Zone Viewer
""""""""""""""""""""""""""""""""""""""""""""""""""""

If the corresponding output files of a Density of states (DOS) or band structure calculation are uploaded, the GIMS output analyzer automatically generates corresponding graphs from this output files. Each graph can be edited by the control elements to the right of the graphs. Editing of the graphs includes:

* Zooming in and out of an energy window
* Changing the colors and thickness of the DOS and band structure lines
* Adapting the axis labels font size and color
* Changing the labels name of the axes

In case of a periodic structure (structure has three lattice vectors) the corresponding Brillouin zone (BZ) and the band path as written in the uploaded in- or/and output files are shown.


Input files
"""""""""""

The code parser also provide the input file for download, in case they have been provided as upload or can be generated from the output file (some codes write the input files in the output file, so we use this to generated them again).

Parsing Errors
""""""""""""""

This section (hopefully) appears, if something went wrong while parsing the output files.

Simple Calculation
==================

The "Simple Calculation" is a workflow to create all needed input files for a single configuration or relaxation calculation.
It connects the following steps:

1. Structure Builder
2. Control Generator
3. Download Input files
4. Output Analyzer

For detailed information about the elemental apps, please visit the corresponding manual page.

Before starting the workflow, please choose your code (top right corner). This can't be reversed once the workflow has started.

The "simple calculation" workflow guides you through the individual steps.
You can navigate forward and backward using the button **Go ahead** and **Go back**, respectively. 
When you reach the download page (*Download the input files and start your calculation*), you can download the prepared input files as a ``input_files.tar`` file.
The file will be saved in your download folder.
(Please note: due to browser security restriction, we cannot let you choose the folder. However, you can change the download folder in the settings of your browser. We recommend to just move the tar file to the right place, manually.)
The tar file contains all needed input files to start the calculation.
Note: GIMS is not capable of running the actual calculation, nor start or submit such a calculation. You have to copy the input files to the computing resources of your choice and run the calculation there.
================
GIMS User Manual
================

This part of the documentation is concerned with the description of the app functionalities. Navigate through the table of contents to learn details about the GIMS application.

.. toctree::

    prerequisites
    mainPage
    structureBuilder
    controlGenerator
    outputAnalyzer
    simpleCalculation
    bandStructure

Control Generator
=================

The *Control Generator* allows to define all necessary numerical settings to run a calculation.

Note: This part is highly Code dependent. Please choose your preferred code in advance in the top right corner (**Choose you code** section). Changing the code after completing the form will delete all progress.

By hovering over the settings text you will get additional information to the right side of the corresponding item.

Downloading control files
-------------------------

If all mandatory fields are filled, you can click on **Generate and download input file(s)**. Then, the control generator creates the input file(s) (input files that contain numerical settings). This also includes corresponding basis sets or basis setting, even if not explicitly requested. (The control generator makes use of code internal defaults as much as possible.)

**FHI-aims:** the ``control.in`` file is generated and the chosen species defaults (light, tight etc.) are pasted into the file as well.

**Exciting:** the ``input.xml`` is generated and all needed species files (``<species>.xml``) are collected, too.

All input files are packed in to a single ``input_files.tar`` file, which is finally saved to your download directory.
Main Page
=========

The figure below shows the main page of the GIMS application. An example can be found here:
https://gims.ms1p.org .

.. image:: ../images/appsMainPage.png



1. **App Selector:** Here, you can select the individual apps. They are described in detail in the subsequent chapters. Hovering over the question mark displays a short description below the app icon. It is further subdivided into two parts: a) workflow apps (*Simple Calculation* and *Band Structure*) and Elemental apps (*Structure Builder*, *Control Generator*, and *Output Analyzer*).

2. **Code Selection:** By clicking on the corresponding code icon the underlying behavior of the apps is adapted to the selected code.

3. **Settings:** The settings menu enables to change display properties, such as the number of shown floating point digits.

4. **GIMS Button:** The GIMS button is the *home button*, that is, clicking on this icon always returns to the apps main page.

5. **Desktop application:** Click on this icon to get a description how to download and install the GIMS offline executable.

6. **Feedback and User Manual:** This *Feedback* icon re-directs to the issue tracker of the GIMS repository. You can leave suggestions, improvements, and bug reports, here. The *User Manual* icon redirects you to the current web site.
Proposal: Phonon Workflow
=========================

This is the sketch for the new phonon workflow. To follow the implementation please checkout the branch:
`28-new-phonon-workflow <https://gitlab.com/gims-developers/gims/-/tree/28-new-phonon-workflow>`_
as well as add your comments in `issue number 28 <https://gitlab.com/gims-developers/gims/issues/28>`_.

Workflow
--------

    1. Structure Builder

        Import/Create structure file. Check, if the generated structure is periodic (i.e. has three lattice vectors). Only allow to proceed if it is periodic.

    2. Control Generator

        Define numerical settings. We will need a new *phonon* section allowing to select:

        * supercell matrix (mandatory); 3x3 Matrix
        * displacement (mandatory, with default 0.01 AA)
        * symmetry threshold (mandatory, with default 1e-5)
        * primitive matrix (optional); radio list: either ``off``, ``auto``, or specify 3x3 matrix. Matrix to obtain the primitive cell. ``off`` switches off this functionality. ``auto`` will use the phonopy auto detection of the primitive cell.

        Moreover, *Calculate Forces* should be set to mandatory.

    3. Download Section

        Provide phonopy information to user. This should contain:

        * Phonopy Version; Citation Reference
        * detected space group and symmetry information
        * Number of atoms in supercell
        * number of needed displacements (thus, the number of needed calculations)
        * Chosen k-grid for supercells
        * A short description how to run the calculations (e.g. we are expecting to preserve the folder structure)

        All the information is contained in the ``phonopy_disp.yaml`` file (generated by ``phonopy``) and can be parsed from there.

    4. Output Analyzer

        The phonon workflow requires the upload of a directory tree, which should be already possible within the current implementation.
        The directory tree of a phonon calculation looks like this::

            phonons
            |--geometry.in
            |--phonopy_disp.yaml
            |--displacement-01
            |  |--aims.out
            |  |--control.in
            |  |--geometry.in
            |--displacement-02
            |  |--aims.out
            |  |--control.in
            |  |--geometry.in
            |--displacement-03
            |  |--aims.out
            |  |--control.in
            |  |--geometry.in
            ...

        It is needed to parse (should be done at the frontend) the forces from the ``aims.out`` in all sub-directories and prepare a special file ``displacement_dataset.json`` with the following structure (from the phonopy documentation)::

          displacement_dataset = {
            'natom': number_of_atoms_in_supercell,
            'first_atoms': [
              {
                'number': atom index of displaced atom (starting with 0),
                'displacement': displacement in Cartesian coordinates,
                'forces': forces on atoms in supercell
              },
                {
                ...
              },
              ...
            ]
          }

        A new section *Phonons* should be added to the site. This section will need some user interaction, since the post-processing has to be initialized. The user has to specify the k-mesh for the post-processing. The specified k-mesh and the ``displacement_dataset.json`` file should be sent to the backend, where the DOS and the Band Structure is calculated. The server answer should contain the DOS and band structure data with the equal data structure as for the electronic band structure and DOS.

        We will need to test the performance of `phonopy` at the backend and limit the calculation time (we will have to set a certain time out and kill phonopy if it exceeds it). Then, notify the user to run the offline app, where this `time out` feature should be deactivated, so the user can do phonopy runs as long as they need.


Backend
-------

For now, we will use directly the phonopy python API. There are two requests needed:

  1. ``phonopy-pre-process``: Generates the displacements and creates a tarball for the download. The response should be the tarball and the phonopy info file.
  2. ``phonopy-post-process``: sends the user settings and the ``displacement_dataset.json`` to the backend. Response should be DOS and Band structure in the json format and phonopy info of the post-processing.

Band Structure
==============

This workflow is similar to the **Simple Calculation** workflow. The only difference is that the needed band path information is automatically generated.
Again, this workflow connects the following steps:

1. Structure Builder
2. Control Generator
3. Download Input files
4. Output Analyzer

This workflow is restricted to work only for periodic systems (for the band structure is a collection of properties of a periodic structure). So make sure you have defined lattice vectors for you structure.

The band path information is automatically generated based on the Bravais lattice from the underlying geometry. More detailed information about the detected Bravais lattice, the special points in the Brillouin Zone and their coordinates are given in the download section.
How to contribute
=================

General git workflow
--------------------

1. `Open a new issue <https://gitlab.com/gims-developers/gims/-/issues>`_ (either to report a bug or to propose a new feature).
2. Go to the issue and below the description of the issue select *Create merge request* (or alternatively choose just *Create branch*). This will create a new branch and open a merge request. The name of the branch is equal to the name of the issue.
3. Update your local repository by typing in the terminal:

    .. code-block:: bash

      git checkout master
      git pull
      git checkout <name-of-new-branch>

4. Now, make your changes. Please commit often and early. Do not commit a huge amount of work at once. Please add a descriptive commit message, so everyone gets an idea about your changes/new implementations. This keeps the process transparent for all contributors. To commit you can, e.g., use:

    .. code-block:: bash

      git add <file-name-1> <file-name-2>
      git commit

5. Push your changes into the remote branch:

    .. code-block:: bash

      git push -u origin <name-of-new-branch>

6. Once you have finished your changes, go back to the browser to the gims GitLab. You will need an approval of your merge request by an eligible user. If your merge request is approved, go to your merge request and select *Merge*.
7. If you have implemented a new or changed an existing feature, please document your changes in the :ref:`changelog`. Further, please also document how to use your new feature in the docs. Please provide information for the user manual.
8. Finally, add tests for the new functionalities for both the client and server side!

Adding support for a new Code
-----------------------------

If you plan to add support for a new electronic structure code, there are only a countable number of steps to do. In what follows gives a rough overview.

Client Side
+++++++++++

1. Add the code icon to ``client/img`` (you'll find examples already there) and add it in the ``index.html`` to the div ``code-selection-box``. For the code name use the format ``<codeName-logo.png>``. Don't use special characters for the code name.
2. Add the definition of the fields needed for the Control Generator. Please use the following name for it: ``Fields_<codeName>.js``. Examples for this file can be found in the folder ``client/src/control-generator-mod``. Import the ``Fields_<codeName>.js`` in the ``AllFields.js`` file.
3. Please add a parser for the output files of your code to the ``output-analyzer-mod``. The prototype is defined as ``Output.js`` in ``client/src/output-analyzer-mod``. Examples for existing supported code parser are in this folder, too. Please add your new parser object to the function ``getOutputInstance`` in ``utils.js`` in the ``output-analyzer-mod`` folder.

Server Side
+++++++++++

1. Please make sure if ASE supports your code. GIMS relies on the auto-detection of the structure-data files. If this is missing, we recommend to add the support for your code first there.
2. Please provide a code class for the new code in the ``app/gims/codes`` directory. A prototype code class and examples for other codes can be found there, too.

Contact
-------

If you experience any problems during the contribution process or if you have any suggestions, feel free to contact me directly via: kokott@fhi-berlin.mpg.de

Developer Manual
================


This part of the documentation is intended to give an overview for developers of GIMS.

.. toctree::

    installation
    codeStructure
    contributing
    deployServer


Deploy the app on a server
==========================

We are assuming that you have a ``ubuntu`` server running with ``apache2``. If you need some help with that we recommend the official documentation for ubuntu here:
`<https://ubuntu.com/tutorials/install-and-configure-apache>`_

The proper configuration of the apache web server might quite specific to your local IT infrastructure, so we cannot describe it here in detail.

Before actually starting you have to install the following in addtion:
::

    sudo apt-get install libapache2-mod-wsgi-py3

Please create folder in ``/var/www/`` named ``gims_gui`` and clone or copy the source code into that folder. This will result in the following directory tree:
::

    /var/www/gims_gui
      |-- gims.wsgi
      |-- gims
            |-- app
                  |-- main.py
                  |-- static

The WSGI file ``gims.wsgi`` should include the following:
::

    #!/usr/bin/python3
    import sys
    import logging
    logging.basicConfig(stream=sys.stderr)
    sys.path.insert(0,"/var/www/gims_gui/gims/app/")

    from main import app as application

Further, in ``/etc/apache2/sites-available`` you should add the following to your web site configuration file:
::

    WSGIScriptAlias / /var/www/gims_gui/gims.wsgi
    <Directory /var/www/gims_gui/gims/app/>
      Order allow,deny
      Allow from all
    </Directory>
    Alias /static /var/www/gims_gui/gims/app/static
    <Directory /var/www/gims_gui/gims/app/static/>
      Order allow,deny
      Allow from all
    </Directory>


Please make sure that you have built the actual application:
::

    cd /var/www/gims_gui/gims/client
    parcel build index.html --out-dir ../app/static --public-url ./

If you have not build the app before, go to the section :ref:`build-app`.
Finally, you have to start the service with:
::

    service apache2 start

If you just want to restart the app after an update, you can simply use:
::

    service apache2 restart

You should now see the web page running. An example for this you can find here: `<https://gims.ms1p.org>`_


Installation
============

If you ended up here to install a local version of GIMS (e.g., to use it offline), make sure that the desktop version provided for download on the GIMS main page is not enough for you. The steps described here require some setup patience. However, you are free to follow them. At the end, you will be able to run the most recent developer version from any branch in the gitlab. So if you can't wait for the next release ... please go ahead!

The installation involves the following three steps:
    1. Download the Code
    2. Build the client in the ``client`` folder (frontend code written in javascript). The result will be copied to ``app/gims/static``.
    3. Install the GIMS python application (backend code written in python) in the folder ``app``.

Obtaining the Code
------------------

To obtain the source code, please clone the repository (make sure you have `deployed you ssh keys on gitlab <https://docs.gitlab.com/ee/ssh/>`_):

.. code-block:: bash

    git clone git@gitlab.com:gims-developers/gims.git

Prerequisites
-------------

We make use of several packages to bundle/build the app and (python) libraries for actually running the app. The most convienient way to quickly install all needed packages is to use the popular package managers ``npm`` and ``pip`` for ``javascript`` and ``python`` packages, respectively. There is a lot of material out there in the so called internet about installing these packages on the OS of your choice. In what follows we assume that they are installed on your system.

A quick overview what is needed to follow the installation instructions:

    * ``node`` (node.js includes the package manager ``npm``)
    * ``python3.6`` (this will include ``pip``)

If you hav install the above, then, the following package managers should be already in place:

    * npm
    * pip

.. _build-app:

Building the Client
-------------------

To bundle and finally build the client side of the app we use `parcel <https://parceljs.org>`_. Please enter the following in your terminal:
::

    npm install -g parcel-bundler

After that we are ready to install the project. Go to the folder ``gims/client/`` and type:
::

    npm install

``npm`` will automtically install the needed modules. The output of the build process is saved to the ``gims/app/static`` folder.

There are two modes availabe:

    1. development mode: This mode uses the ``parcel watch`` option. Parcel rebuilds the app automatically, whenever a change to the javascript files has been made. There is no need to manually recompile the app over and over again while developing. Please type for that:
    ::

        npm run dev

    2. production mode: This mode uses the ``parcel build`` option. It builds the app once and that's it. To use this mode, type:
    ::

        npm run build

To check whether the build process was successfull, please check if the folder ``gims/app/static`` has any content. If it has we are finished with building of the client side of the application.


Install the GIMS python application
-----------------------------------

We now turn to the server part of your local installation. Please change to the folder ``gims/app``. You should see now some python files. To isolate your application from dependencies from the rest of your machine we recommend to use python virtual environment. We will briefly discuss the setup of such an environment, first.

Setup a virtual environment
"""""""""""""""""""""""""""

Please install ``virtualenv`` with:
::

    pip install virtualenv

Now, setup a new virtual environment:
::

    python3 -m venv env

The settings of this virtual environment are saved in the folder `gims/app/env`.

Activate the virtual environment with:
::

    source env/bin/activate

Install GIMS app
""""""""""""""""

GIMS can be installed like any other python package. Just got to ``app`` folder and type :
::

    pip install -e .[dev]

All dependencies should be installed automatically. A list for all GIMS dependencies can be found in the ``setup.py`` file in the ``app`` folder.

Start the flask server
""""""""""""""""""""""

For the use of GIMS on your local machine we will use the development server from Flask. Please go to the folder ``app`` in the `gims` folder and type the following in the terminal:
::

    export FLASK_APP=gims
    flask run

To check if everything works as expected go to your browser and visit (enter in the URL bar):
::

    localhost:5000

That's it! It was a long road to go, though.

Code Structure
==============


GIMS is a web client-server system and so less surprisingly there are two parts of GIMS: The client side is responsible of the user interaction and data visualization and the server part holds the scripts for pre-processing/post-processing information before/after the simulation.

The client side is written in Javascript, while the server side makes use of python. You may ask, whether it would have better to keep the whole application more homogeneous. However, we wanted to have the optimal experience for the end user and making use of already existing libraries in our community, which are mainly written in python. On the one hand using a `javascript-only` framework would have required us to stupidly rewrite big code blocks in javacsript. On the other hand having a minimal frontend solution might have negatively affected the user experience, since communication of big files might significantly reduce the responsiveness of the application. Nowadays, the browser has become a more and more powerful engine by itself: for us a two-component application evolved, naturally.

The following figure may grant you glance on the overall code structure. It also shows the main dependencies of the package on each of the sides.

.. image:: ../images/software_architecture.png

**References to other code projects:**

* Flask: `<https://flask.palletsprojects.com/en/1.1.x/>`_
* ASE: `<https://wiki.fysik.dtu.dk/ase/index.html>`_
* Spglib: `<https://spglib.github.io/spglib/>`_
* threejs: `<https://threejs.org>`_
Release Notes
===============

The following threads will keep you updated on the upcoming GIMS releases.


Release Plans
-------------

This section announces the next planned releases and updates of the running GIMS servers.

October 6, 2020
+++++++++++++++

Release of version 1.0.0 and simultaneous update of the servers https://gims.ms1p.org and https://gims-dev.ms1p.org.
For any future patch or minor update, only the gims-dev.ms1p.org server will be updated first.
Only an official release will trigger an update of the main server gims.ms1p.org.

.. _changelog:

Changelog
---------

The section *Changelog* describes and documents the implementation of new features, fixes and so on, of all releases.

.. toctree::

    changelog
Release 1.0.2
-------------

Official Release for the JOSS Paper! Thanks for all who contributed!
A special thanks to JOSS and the reviewer for the immense helpful comments!

* Restructuring of the backend 
* Improvement of the documentation
* Improvement of the paper
* Adding a docker file
* Adding contribution guidelines

And a few small other things.


Release 1.0.0
-------------

*(From June 03, 2020, to October 06, 2020)*

With this new release we introduce semantic versioning of GIMS. The versioning we had before is deprecated. Once this release is finished, the corresponding commit will be tagged

Changes:

* **Structure Builder**
    * Multiple structure file upload: user can select several structure files. All uploaded structures are accessible via tabs.
    * New Supercell Button: user can create supercell either with simple multiples of lattice vectors (3 parameters) or with a supercell matrix (9 parameters). The new supercell will be created in a new tab.
    * Structure Info field: user gets additional information about structure, e.g. chemical formula, space group, occupied Wyckoff positions, ect. This feature integrate spglib in the backend of GIMS.
    * Primitive Button: New button in the side panel. It creates a new primitive structure (periodic systems only) in a new tab.
    * Brillouin Zone Viewer: If a periodic structure is loaded, the GIMS structure builder will automatically view the Brillouin zone of the Bravais lattice of the loaded structure. The BZ path according to Setyawan/Curtarolo, the cartesian Axis, and the reciprocal unit cell are shown as well.

* **Control Generator**
    * Units (if applicable) are added to the input values
    * Rework input value structure (code structure): Each code has its own list of supported input values. Makes it more easier to include new codes.

* **Output Analyzer**
    * Brillouin Zone Viewer added. BZ viewer displays actual path used in the calculation.

* **GIMS Backend**
    * Restructuring of the backend: Backend is now organized as a python package (which simply allows a ``pip install .`` of the backend).
    * Extending backend tests.

* **GIMS Frontend**
    * Restructuring of the dependencies: Some libraries, where the source code has been included in GIMS, have been removed. Instead, the corresponding libraries have been added as a dependency to the ``package.json`` file.

* **GIMS landing page**
    * Small redesign to adapt the content to the description in the paper.
    * Add linkt to running GIMS versions and to the release notes.


Seventh release (*r7*)
----------------------

*(From March 18, 2020, to April 7, 2020)*

This release improves several aspects of the project like testing, code documentation and adds some enhancements:

- Enhancement: **Consistent browser and in-app navigation**. The usage of the browser native navigation elements is now consistent with the in-app navigation elements for the stand-alone apps as well as for the workflows.

Some important changes were required in the application:

- Support of **multiple instances of modules and nested modules** inside another (e.g. `StructureBuilder` inside a workflow module). The `StructureBuilder` state is not a singleton longer.
- Complete **rewriting of the workflow implementation**: The workflows are now implemented like application modules (`Workflow.js`) and contain UI components (modules or regular/simpler component)
- URL based **navigation point for every step in workflows**: The URL is this way: `#WorkflowName-workflow#StepComponentId`
- Wiki **documentation** piece updated and extended: `Implementation details: workflows, application navigation <https://gitlab.com/gims-developers/gims/-/wikis/Implementation-details:-workflows,-application-navigation>`_. And `Gitlab issue explaining the development process <https://gitlab.com/gims-developers/gims/-/issues/29>`_
- **Extension of the Front-End CI Testing-Suite**. Some parts of the client are tested:
    - ``Output Analyzer`` module testing
    - The **main output file parsing** is tested for **many combination of calculations** (Single-point calculation, relaxation / Spin none, collinear / Band structure+DOS / Periodic, non-periodic / Exciting, FHI-aims).
    - Files parsing. The parsing of **all the files making up the output** is tested.
    - The parsing of the **files containing the DOS and BandStructure data** is tested.
- **File formats generated for exporting** in `Structure Builder`. The functions that generate the different combination of file formats (molecule/periodic, aims/exciting, fractional/cartesian coordinates) to be exported are tested.
- `More details in this Gitlab issue <https://gitlab.com/gims-developers/gims/-/issues/32>`_

- **Code Documentation of Front-End**. All the application code has been documented according to the **JSDoc** style

- Enhancement: **Refactoring Output-Analyzer**. The *OutputAnalyzer* module has been deeply refactored, both the *FHIaims* and *Exciting* code parsers (Sebastian) and the user interface components: The `OutputAnalyzerMod.js` has been refactored to (electronic-structure-)code-independent.

- Some minor UI improvements

Sixth release (*r6*)
--------------------

*(From November 18, 2019, to February 4, 2020)*

The main feature in this release is the addition of *Exciting* code support. On the other hand, we've added many interactivity possibilities on graphs, mainly on the *Band Structure* and *DOS* plots. As this was such important release we took the opportunity, as well, to carry out several meaningful improvements (regarding design and UX) along the application:

Main r6 features:

- **Exciting code support** (major application refactoring and new functionality):
    - Code selector IU component and corresponding new application state and logic
    - Update of the import/export functionality (*StructureBuilder* module)
    - Control generator module adaptation: new form fields, server side development (Sebastian) and integration with client.
    - Workflow adaptation (code election disabled inside, workflow configuration improvement, etc)
- **Advanced interactivity in graphs**:
    - Structure Builder: possibility of color change of the showing species (by clicking at the species circles on the legend)
    - Band Structure and DOS plots:
        - Possibility of taking screenshots of the plots
        - The plot labels are now editable and adjustable in size and color and axes ticks labels in size
        - The lines (by spin) are adjustable in color and thickness
        - Button that returns to the original state of the plot
        - Zooming improvements
- Final release improvements:
    - Explanatory text for the top level page buttons
    - **Back/Forward browser buttons support** (consistent application behavior)
    - **Application Settings**: UI and event system to change characteristics along the application. First case: *Number of decimal digits*
    - **Tooltip system for user help**. First application: field labels in *ControlGenerator*
    - User Feedback button pointing to the issue tracker on *Gitlab*
    - **Application breadcrumbs**: navigation system that shows the user where (module) is in the application
    - Several minor design and interaction improvements along the application

More work done this time:

- New own **2D graphic library** (`Canvas.js` file)
- *InfoCanvas* and DOS/BS plots refactoring (on the new library) and performance optimization.
- Color picker integration (first use: species color change - Structure builder)
- Application refactoring (images imports) and application resources review and cleaning
- *Control Generator* improvement for the *Band Structure* workflow
- Client adaptation to support `tar` file. Input files requested are now bundled (`tar` format) as the response
- Beginning of the formal code documentation

Fifth release (*r5*)
--------------------

*(From September 19 to November 15, 2019)*

The main feature in this release is the addition of a Band Structure and DOS workflow and their visualization possibility on an improved *Output Analyzer*. On the other hand, we have found, tested and implemented a way to bundle and distribute the application as a desktop/offline app.

Main r5 features:

- **Band Structure and DOS data visualization** integrated in the *Output Analyzer*:
    - *Output Analyzer* file import logic and multi-file support (folder source support).
    - Bands and DOS files data parsing and flexible graphical representation
- New **Band Structure & DOS workflow** and improved *Control Generator*:
    - New *Band Structure & DOS workflow*
    - Improved *Control Generator*: More flexible and supporting sections (Accordion user interface)
    - Backend based *control.in* generation (use the *control.in* writer from ASE). (by Sebastian)
- Application **packaging and distribution as a desktop/offline** app research and first implementation:
    - Research into the *PyInstaller* tool and decision of adoption
    - Testing on Linux and MacOS, optimization and documentation
    - The desktop application (versions for Linux and MacOS) is already downloadable (via the web application)
- Top level/dashboard restructuring:
    - Change of workflows, layout and buttons

More work done this time:

- **Error handling** through the application. Refactoring and documentation ([Application error handling](Application-error-handling) doc on the Wiki)
- Structure builder improvements (by Sebastian)
    - Simple supercell handling at the right panel of the structure builder UI
    - Constrain atoms via a checkbox at the basis atoms panel. The material of the corresponding atom changes.
    - The distance between atoms is checked. If they become too close, then, the bond color changes to red and the bond radius is increased.
- User feedback UI component at application level. A way to give simple feedback to the user (errors, warnings and informative tips). It's a text box at the top-center of application layout
- Meaningful modification of FileImporter component
- Simpler console logging/error call (Global utility functions)


Fourth release (*r4*)
---------------------

*(From August 5 to September 18, 2019)*

This time we kept expanding the application: a new module *Control generator* has been added, as well as a new important application feature: workflows support.

Main r4 features:

- ***Control Generator* application module**:
    - It generates *control.in* files from a user form and species defaults files
    - New generic ***Form*** component
    - URL fragment associated to the module (*#ControlGenerator*)
- **Workflows** support:
    - UI header (*AssistantHeader*) enables the workflow interaction
    - Declarative (partially) workflow definition
    - Two workflows (similar, only declarative difference) implemented
    - Modules adaptations to be integrated in the workflows

More work done this time:

- ***ASE* library integration** and server request based geometry importer (*Structure Builder* module)
- Import **geometry file format resolution** user interface
- Some work on the adaptation and integration of the *DOS* and *Band Structure* graphs (from the *NOMAD* project).
- Some refactoring (asynchronous programming based since now on *promises* and *async/await*) and documentation (Workflows feature)


Third release (*r3*)
--------------------

*(From May 2 to July 26, 2019)*

From this version, the application is more than a *Structure builder*. This becomes a module in the app and a new one is added: the *Output Analyzer*.

- Improvements on the previous (r2b) release:
    - Supercell construction based on a supercell matrix (based on Sebastian’s code)
    - Refactoring: ES6 import/export and dependency reduction on the math.js library (goal: get rid of it)
    - New Switch UI component (used for Atom displacement)
    - Add a checkbox to the lattice vector section to enable/disable the scaling of atom positions as lattice vectors change

Main r3 features:

- **Non-periodic systems** support (Structure Builder module): file parsing, Structure viewer and UI adaptation
    - Lattice vectors can be removed or created at any moment and circumstance
    - The lattice vectors creation/removals are supported by the undo/redo system
- Top level **application module: Dashboard**.
- Application **Output Analyzer module**:
    - Output file(s) importer: simple and relaxation calculation output support (integration of Sebastian’s file parsing code)
    - Module page layout (divided in sections)
    - Integration of the StructureViewer and development of structure animations support
    - 2D charts integration and enhancement from Sebastian’s code
    - Input files: geometry.in and control.in popup viewer and downloader
    - URL fragment associated to the module (#OutputAnalyzer)

More work done this time:

- Refactor the **StructureViewer as an app library** (reusable application module)
- Create a changelog document (Gitlab wiki)
- Research and new doc (Gitlab wiki): **Web client architecture and technical decisions**
- Modal popup component
- Some research and discussion about the next big goal: the first complete workflow


Second release - second part (*r2b*)
------------------------------------

*(From March 5 to April 29, 2019)*

- Improvements on the previous release:
    - In an atom displacement, new bonds are calculated and automatically updated on the viewer and the new atom coordinates are shown (atom info area) at mouse release time
    - Fixed performance penalty issue (mathjs library methods) related to bonds calculation
    - Screenshot feature improvement: take a screenshot in one click
    - Unify and make well-proportioned the arrows of the lattice vectors and labels for any cell size.
    - The corresponding atom is highlighted when the mouse pointer hovers the atom row (right panel)

Main features:

- **Measurements** of distance between atoms (2 atoms selected), angle (3 atoms selected) and torsion angle (4 atoms selected). The info is shown on the top-right part of the viewer.

- Possibility of creating **supercells** from the current cell. URL fragments interface: this funtion is launched by adding a fragment to the main URL. e.g. `URL#repeat:2:2:3`

- **Undo/redo feature**: system for going back and forth in the history of changes (both via UI and key combination Ctrl+z Ctrl+y). Changes supported: atom creation, move, species change, removal, lattice vectors modification and supercell transformation.

More work done:

- Important refactoring. The `Reactjs` library was removed from the project and implemented a new client architecture.
- **Transfer the application to the FHI infrastructure**
- **Performance and memory optimization** research, testing and implementation of measures
- We started the **project documentation** on gitlab (wiki section) with three pages: Development roadmap, Performance optimization and Setup WSGI for gui on apache2. We intent to document every meaningful aspect of the project from now
- **Atom selection mechanism** (multi-selection enabled). Now there are 3 atom interaction ways: hovering, selection and drag&drop. Selection: Shift key + mouse click. Clear selection: mouse right-click
- Formalization of the URL fragment (browser feature) use (starting with the hash symbol - #). For now only param repeat supported


Second release - first part (*r2a*)
-----------------------------------

*(From early January to March 1, 2019)*

- Improvements of the previous release:
    - Bonds representation improvement (algorithm to identify bonds across cell boundaries)
    - UI element to switch between fractional and cartesian coordinates
    - More digits (5) in the coordinate fields. Several solutions were considered and visually tested
    - More options (left-top checkboxes) added to the 3D viewer:
        - Option of switching on/off the lattice vectors
        - Option of switching on/off the unit cell box
        - Wrap atoms into cell (If atoms are outside the box, move them to the equivalent positions inside the box)
        - Option of switching on/off the bonds across the cell boundaries

Main features:

- Structure **creation from scratch**
- **Manipulation of imported** (from files) structures. You can (editing text fields):
    - edit the lattice vectors (modify the cell) and the atoms coordinates (move the atoms)
    - change the species (substitution of atoms)
    - add and remove atoms
- **Enabled interaction** with the structure elements **on the 3D viewer**. You can:
    - Move atoms dragging them (mouse click and release)
    - Show textual info (on the viewer) of the atom hovered by the mouse pointer
    - Atom highlights (in 3D view) by hovering the circle in the atom row (side panel)
- **Export** of the editor current state of the structure in 'geometry.in' format (both in fractional and cartesian coordinates)

More work included in this release:

- **Screenshot capture** functionality

- Moving atoms with the mouse enabled via a checkbox (to avoid accidents)

- Side panel **new button bar** (*Import - Export - Reset*) and user interaction

- Default structure when the web app is loaded

First release (r1)
------------------

*(From mid November to the end of December, 2018)*

Main features:

- **Work environment set-up**: `Reactjs` (frontend library) + `Flask` (python backend) + `GoogleAppEngine` (cloud hosting)

- **Structure 3D viewer**: based on `Threejs` library (using *WebGL* browser capability)
  - Basic cell, atoms and lattice vectors and parameters representation
  - Interaction: zoom and rotation
  - Bonds and atoms on cell boundaries
  - Legend and structure visualization control

- **Side panel** showing lattice vectors and atoms alphanumeric data

- **Import of existing structures** feature. Formats supported:
  - FHI-aims format ‘geometry.in’
  - Crystallographic Information File (CIF)
