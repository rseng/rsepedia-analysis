# <img src="https://github.com/SchmollerLab/Cell_ACDC/blob/main/cellacdc/resources/icons/assign-motherbud.svg" width="80" height="80"> Cell-ACDC

### A GUI-based Python framework for **segmentation**, **tracking**, **cell cycle annotations** and **quantification** of microscopy data

*Written in Python 3 by Francesco Padovani and Benedikt Mairhoermann.*

[![build ubuntu](https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-ubuntu.yml/badge.svg)](https://github.com/SchmollerLab/Cell_ACDC/actions)
[![build macos](https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-macos.yml/badge.svg)](https://github.com/SchmollerLab/Cell_ACDC/actions)
[![build windows](https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-windows.yml/badge.svg)](https://github.com/SchmollerLab/Cell_ACDC/actions)
[![Python version](https://img.shields.io/pypi/pyversions/cellacdc)](https://www.python.org/downloads/)
[![pypi version](https://img.shields.io/pypi/v/cellacdc?color=red)](https://pypi.org/project/cellacdc/)
[![Downloads](https://pepy.tech/badge/cellacdc/month)](https://pepy.tech/project/cellacdc)
[![License](https://img.shields.io/badge/license-BSD%203--Clause-brightgreen)](https://github.com/SchmollerLab/Cell_ACDC/blob/main/LICENSE)
[![repo size](https://img.shields.io/github/repo-size/SchmollerLab/Cell_ACDC)](https://github.com/SchmollerLab/Cell_ACDC)
[![DOI](https://img.shields.io/badge/DOI-10.1101%2F2021.09.28.462199-informational)](https://www.biorxiv.org/content/10.1101/2021.09.28.462199v2)

<div align="left">
  <img src="https://github.com/SchmollerLab/Cell_ACDC/blob/main/cellacdc/resources/figures/Fig1.jpg" width="700" alt><br>
    <em>Overview of pipeline and GUI</em>
</div>

## Resources

- [User Manual](https://github.com/SchmollerLab/Cell_ACDC/blob/main/UserManual/Cell-ACDC_User_Manual.pdf) with **detailed instructions**
- [Pre-print](https://www.biorxiv.org/content/10.1101/2021.09.28.462199v2) of Cell-ACDC publication
- [Forum](https://github.com/SchmollerLab/Cell_ACDC/discussions) for discussions (feel free to **ask any question**)
- **Report issues, request a feature or ask questions** by opening a new issue [here](https://github.com/SchmollerLab/Cell_ACDC/issues).
- Twitter [thread](https://twitter.com/frank_pado/status/1443957038841794561?s=20)

## Overview

Let's face it, when dealing with segmentation of microscopy data we often do not have time to check that **everything is correct**, because it is a **tedious** and **very time consuming process**. Cell-ACDC comes to the rescue!
We combined the currently **best available neural network models** (such as [YeaZ](https://www.nature.com/articles/s41467-020-19557-4),
[Cellpose](https://www.nature.com/articles/s41592-020-01018-x), [StarDist](https://github.com/stardist/stardist), and [YeastMate](https://github.com/hoerlteam/YeastMate)) and we complemented them with a **fast and intuitive GUI**.

We developed and implemented several smart functionalities such as **real-time continuous tracking**, **automatic propagation** of error correction, and several tools to facilitate manual correction, from simple yet useful **brush** and **eraser** to more complex flood fill (magic wand) and Random Walker segmentation routines.

See below **how it compares** to other popular tools available (*Table 1 of our [pre-print](https://www.biorxiv.org/content/10.1101/2021.09.28.462199v2)*).

<p align="center">
  <img src="https://github.com/SchmollerLab/Cell_ACDC/blob/main/cellacdc/resources/figures/Table1.jpg" width="600">
</p>


## Is it only about segmentation?

Of course not! Cell-ACDC automatically computes **several single-cell numerical features** such as cell area and cell volume, plus the mean, max, median, sum and quantiles of any additional fluorescent channel's signal. It even performs background correction, to compute the **protein amount and concentration**.

You can load and analyse single **2D images**, **3D data** (3D z-stacks or 2D images over time) and even **4D data** (3D z-stacks over time).

Finally, we provide Jupyter notebooks to **visualize** and interactively **explore** the data produced.

**Do not hesitate to contact me** here on GitHub (by opening an issue) or directly at my email francesco.padovani@helmholtz-muenchen.de for any problem and/or feedback on how to improve the user experience!

## Update v1.2.4
First release that is finally available on PyPi.

Main new feature: custom trackers! You can now add any tracker you want by implementing a simple tracker class. See the [manual](https://github.com/SchmollerLab/Cell_ACDC/blob/main/UserManual/Cell-ACDC_User_Manual.pdf) at the section "**Adding trackers to the pipeline**".

Additionally, this release includes many UI/UX improvements such as color and style customisation, alongside a light/dark mode switch.

## Update v1.2.3

**NOTE: some users had issues installing the environment with this version. Please see this [issue](https://github.com/SchmollerLab/Cell_ACDC/issues/5) for a possible solution**

This release includes new segmentation models:
- Cellpose v0.8.0 with the models cyto2 and omnipose
- StarDist

## Update v1.2.2

This is the first release with **full macOS support**! Additionally, navigating through time-lapse microscopy data is now up to **10x faster** than previous versions.
More details [here](https://github.com/SchmollerLab/Cell_ACDC/releases/tag/v1.2.2)

## Installation using Anaconda (recommended)

*NOTE: If you don't know what Anaconda is or you are not familiar with it, we recommend reading the detailed installation instructions found in manual [here](https://github.com/SchmollerLab/Cell_ACDC/blob/main/UserManual/Cell-ACDC_User_Manual.pdf).*

1. Install [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for **Python 3.9**. *IMPORTANT: For Windows make sure to choose the **64 bit** version*.
2. Update conda with `conda update conda`. Optionally, consider removing unused packages with the command `conda clean --all`
3. Create a virtual environment with the command `conda create -n acdc python=3.9`
4. Activate the environment `conda activate acdc`
5. Install Cell-ACDC with the command `pip install cellacdc`

## Installation using Pip

1. Download and install [Python 3.9](https://www.python.org/downloads/)
2. Upgrade pip with `pip install --updgrade pip`
3. Navigate to a folder where you want to create the virtual environment
4. Create a virtual environment: Windows: `py -m venv acdc`, macOS/Unix `python3 -m venv acdc`
5. Activate the environment: Windows: `.\acdc\Scripts\activate`, macOS/Unix: `source acdc/bin/activate`
6. Install Cell-ACDC with the command `pip install cellacdc`

## Install from source

If you want to contribute or try out experimental features (and, if you have time, maybe report a bug or two :D), you can install the developer version from source as follows:

1. Open a terminal and navigate to a folder where you want to download Cell-ACDC
2. Clone the repo with the command `git clone https://github.com/SchmollerLab/Cell_ACDC.git` (if you are on Windows you need to install `git` first. Install it from [here](https://git-scm.com/download/win))
3. Install [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4. Update conda with `conda update conda`. Optionally, consider removing unused packages with the command `conda clean --all`
5. Create a new conda environment with the command `conda create -n acdc python=3.9`
6. In the terminal, navigate to the `Cell_ACDC` folder that you cloned before and install Cell-ACDC with the command `pip install -e .`.

## Running Cell-ACDC

1. Open a terminal (on Windows use the Anaconda Prompt if you installed with `conda` otherwise we recommend installing and using the [PowerShell 7](https://docs.microsoft.com/en-us/powershell/scripting/install/installing-powershell-on-windows?view=powershell-7.2))
2. Activate the environment (conda: `conda activate acdc`, pip on Windows: `.\env\Scripts\activate`, pip on Unix: `source env/bin/activate`)
3. Run the command `acdc` or `cellacdc`

## Usage

For details about how to use Cell-ACDC please read the User Manual downloadable from [here](https://github.com/SchmollerLab/Cell_ACDC/tree/main/UserManual)
# Custom annotations in the main GUI

Currently we are implementing the following types of annotations:

1. Single time-point
2. Multiple time-points
3. Multiple values class

As of April 2022 we only implemented 1. Single time-point.

## Single time-point annotation

If the user loads time-lapse data then the annotation mode needs to be activated from the `modeComboBox` control. Otherwise it is always active.

### Construction

A custom annotation can be added by clicking on the `addCustomAnnotationAction`
tool button on the left side toolbar. This action is connected (with `triggered` slot) to the `addCustomAnnotation` function.

This function will open the dialog `apps.customAnnotationDialog` (with custom `exec_` function) and the window reference is stored in `addAnnotWin` attribute of the gui.

The `addAnnotWin` window has the following attributes created in its `closeEvent`:

- **symbol**: any of the `pyqtgraph` valid symbols (only the string inside the quotes, e.g., 'o' for circle, see `widgets.pgScatterSymbolsCombobox`)

- **keySequence**: a `QKeySequence` built with valid PyQt shortcut text, see [here](https://doc.qt.io/qt-5/qkeysequence.html#QKeySequence-1). Note that macOS shortcut strings are converted to valid PyQt string using the `widgets.macShortcutToQKeySequence` function.

- **toolTip**: formatted text for `setToolTip` method of the button

- **state**: a dictionary with all the information needed to restore the annotations parameters, such as 'name', 'type', 'symbol'. Note that 'symbolColor' is a `QColor`.

With this info, when the window is closed with the `Ok` button, a tool button is added to the toolbar with the function `addCustomAnnotationButton`. This function adds a `widgets.customAnnotToolButton` (custom `paintEvent`).

The tool button reference is used as a key to create a dictionary inside the `customAnnotDict` (initialized in `__init__` method of the gui). This dictionary has four keys:

1. 'action' with the action linked to the tool button
2. 'state' with the dictionary of `addAnnotWin.state`
3. 'annotatedIDs' with a list of dictionaries, one dictionary per loaded position (same length as `self.data`). Each one of these dictionaries will be populated with the `frame_i` as key and, as value, the list of annotated IDs with that particular button.
4. 'scatterPlotItem', see below.

Next, the parameters of the annotation are saved as a json file to both `cellacdc/temp/custom_annotations.json` path (initialized as a global variable after the imports) and to the `<posData.images_path>/<basename>custom_annot_params.json` path (initialized in `load.loadData.buildPaths`). The saving is performed in the `saveCustomAnnot` function of the gui.

Finally, we create a `pg.ScatterPlotItem` and we add it to the `ax1` (left plot). We also add a column of 0s to the `acdc_df` filled with 0s and with column name = self.addAnnotWin.state['name'].

### Usage

The user clicks on the tool button of the annotation. This tool button is connected to `customAnnotButtonClicked`. Additionally, the tool button also has a right-click context menu with the following actions:

- **sigRemoveAction** --> `removeCustomAnnotButton`
- **sigKeepActiveAction** --> `customAnnotKeepActive`
- **sigModifyAction** --> `customAnnotModify`

At this point, the user can **RIGHT-click** on any segmented object, **only on the left image**. The click can be used also for undoing annotation.

The annotation is performed by the `doCustomAnnotation` function called in the case `elif isCustomAnnot:` of `gui_mousePressEventImg1`.

Note that `doCustomAnnotation` is also called at the end of `updateALLimg` to annotate every time the image changes.

The `doCustomAnnotation` function will check which button is active (from the keys of `customAnnotDict`). The checked button reference is used to access the `customAnnotDict[button]['annotatedIDs']` list of dictionaries. This list is indexed with `self.pos_i`. The resulting dictionary is indexed with `frame_i` to get the list of annotated IDs for the specific frame/position.

If the clicked ID is in `annotIDs_frame_i` then it is removed because the user is asking to undo the annotation, otherwise it is appended.

Next, we get the centroid coordinates of the annotated ID and we draw the scatter symbol. The 'scatterPlotItem' is stored in the `customAnnotDict` dictionary.

Finally, we reset to 0 the column in the acdc_df and we write 1 on the newly annotated IDs.

### Restoring of annotations after GUI is closed and re-opened

Steps for restoring:

1. `load.loadData` loads the attribute `posData.customAnnot` in the `loadOtherFiles` method. This attribute is set to the output of `load.read_json` function. The json file contains a dictionary where the keys are the names of each custom annotation, and the values are the `state` attribute of the `addAnnotWin` saved in `saveCustomAnnot` (with 'symbolColor' converted from QColor to rgb for json saving).
2. At the end of `load.loadData.loadOtherFiles` we call the function `load.loadData.getCustomAnnotatedIDs` that retrieves the annotated IDs from the acdc_df and it adds them to the `customAnnotIDs` dictionary (empty if acdc_df or customAnnot not found). The `customAnnotIDs` dictionary has the name of the saved custom annotations as keys, and, as values, a dictionary with frame_i as keys and the list of annotated IDs for that frame as values.
3. At the end of `gui.loadingDataCompleted` or in `gui.next_pos`/`gui.prev_pos` we call `addCustomAnnotationSavedPos` which iterates the items of `posData.customAnnot` where the keys are the names of the annotation and the values are the state saved in the json files:
    - For each annotation we re-build `symbolColor`, `keySequence` and `toolTip` and we add the tool button with `addCustomAnnotationButton`
    - If the tool button is not already present in the `customAnnotDict` we create it (see the **Construction** section above) and we add the `posData.customAnnotIDs` to the 'annotatedIDs' value `customAnnotDict`. Finaly, We add the scatter plot item with the function `addCustomAnnnotScatterPlot`
    - If the tool button is already present we only add `posData.customAnnotIDs` to the `customAnnotDict['annotatedIDs']` value.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I feel like this repetetive action could be automated [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Informal request
about: If you are not sure you have a bug nor a feature request use this free template
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - Python version: [e.g. Python 3.9]
 - OS: [e.g. Windows 10]
 - Cell-ACDC Version [e.g. v1.2.4, you can find this from the menu "Help-->About Cell-ACDC" on the main Launcher]

**Additional context**
Add any other context about the problem here.
Getting started
---------------
.. image:: images/logo.svg
   :align: left
   :width: 60

Cell-ACDC
=========
|

|BuildUbuntu| |BuildMacOS| |BuildWindows| |PythonVersion| |Licence| |PiPyVersion| |RepoSize|

.. |BuildUbuntu| image:: https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-ubuntu.yml/badge.svg
   :target: https://github.com/SchmollerLab/Cell_ACDC/actions
.. |BuildMacOS| image:: https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-macos.yml/badge.svg
   :target: https://github.com/SchmollerLab/Cell_ACDC/actions
.. |BuildWindows| image:: https://github.com/SchmollerLab/Cell_ACDC/actions/workflows/build-windows.yml/badge.svg
   :target: https://github.com/SchmollerLab/Cell_ACDC/actions
.. |PythonVersion| image:: https://img.shields.io/pypi/pyversions/cellacdc
   :target: https://www.python.org/downloads/
.. |Licence| image:: https://img.shields.io/badge/license-BSD%203--Clause-brightgreen
   :target: https://github.com/SchmollerLab/Cell_ACDC/blob/main/LICENSE
.. |PiPyVersion| image:: https://img.shields.io/pypi/v/cellacdc?color=red
   :target: https://pypi.org/project/cellacdc/
.. |RepoSize| image:: https://img.shields.io/github/repo-size/SchmollerLab/Cell_ACDC
   :target: https://github.com/SchmollerLab/Cell_ACDC

Welcome to Cell-ACDC's documentation!
-------------------------------------

.. note::

   The building of this documentation is under active development. Please, refer to our `User Manual <https://github.com/SchmollerLab/Cell_ACDC/blob/main/UserManual/Cell-ACDC_User_Manual.pdf>`_ in the meanwhile

Cell-ACDC is a GUI-based Python framework for **segmentation**, **tracking**, **cell cycle annotations** and **quantification** of microscopy data.

You can load and analyse **2D, 3D** (either single z-stacks or 2D images over time) and **4D** (3D z-stacks over time) images.

Additionally, you can load **as many additional fluorescent channels** as you wish. Cell-ACDC will then compute many **numerical features** for each segmented cell, such as mean, sum, max, quantiles etc.
It also performs **automatic background correction** and computes **protein amount**.

Other numerical features computed are **cell volume**, and morphological properties of the segmented object.

Resources
---------

* `User Manual`_ with **detailed instructions**
* `Pre-print`_ of Cell-ACDC publication
* `Forum`_ for discussions (feel free to **ask any question**)
* **Report issues, request a feature or ask questions** by opening a `new issue`_
* Twitter `thread`_

.. _User Manual: https://github.com/SchmollerLab/Cell_ACDC/blob/main/UserManual/Cell-ACDC_User_Manual.pdf
.. _Pre-print: https://www.biorxiv.org/content/10.1101/2021.09.28.462199v2
.. _Forum: https://github.com/SchmollerLab/Cell_ACDC/discussions
.. _new issue: https://github.com/SchmollerLab/Cell_ACDC/issues
.. _thread: https://twitter.com/frank_pado/status/1443957038841794561?s=20

Overview
--------

Let's face it, when dealing with segmentation of microscopy data we often do not have time to check that **everything is correct**, because it is a **tedious** and **very time consuming process**. Cell-ACDC comes to the rescue!
We combined the currently **best available neural network models** (such as `YeaZ <https://www.nature.com/articles/s41467-020-19557-4>`_,
`Cellpose <https://www.nature.com/articles/s41592-020-01018-x>`_, `StarDist <https://github.com/stardist/stardist>`_, and `YeastMate <https://github.com/hoerlteam/YeastMate>`_) and we complemented them with a **fast and intuitive GUI**.

We developed and implemented several smart functionalities such as **real-time continuous tracking**, **automatic propagation** of error correction, and several tools to facilitate manual correction, from simple yet useful **brush** and **eraser** to more complex flood fill (magic wand) and Random Walker segmentation routines.


.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   getting-started
Installation
============
