# Chemiscope: interactive structure-property explorer for materials and molecules

![tests](https://github.com/lab-cosmo/chemiscope/workflows/Tests%20&%20Lints/badge.svg)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02117/status.svg)](https://doi.org/10.21105/joss.02117)

Chemiscope is an graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems; as well as a library of
re-usable components useful to create new interfaces.

![Default interface of chemiscope](docs/src/img/screenshot.png)

## Citing chemiscope

Chemiscope is distributed under an open-source license, and you are welcome to
use it and incorporate it into your own research and software projects.
If you find it useful, we would appreciate a citation to the chemiscope
[paper](https://doi.org/10.21105/joss.02117):

> G. Fraux, R. K. Cersonsky, M. Ceriotti, _Chemiscope: Interactive
> Structure-Property Explorer for Materials and Molecules._ **Journal of Open
> Source Software** 5 (51), 2117 (2020)

If you incorporate chemiscope components into a software project, a link back to
the chemiscope homepage (https://chemiscope.org) is the preferred form of
acknowledgement.

## [Documentation](https://chemiscope.org/docs/)

You may be interested in particular about how to [create a visualization of your
own dataset](https://chemiscope.org/docs/tutorial.html#input-file-format-for-chemiscope).

## Getting help for using chemiscope

If you want to get help when using chemiscope either as a JavaScript/TypeScript
library inside your own project; or for creating input files for the default
visualizer at https://chemiscope.org, you can open a [Github
issue](https://github.com/lab-cosmo/chemiscope/issues/new) with your question;
or send an email to the developers (you can find these emails on the lab
webpage: https://www.epfl.ch/labs/cosmo/people/)

## Getting and running the code

```bash
git clone https://github.com/lab-cosmo/chemiscope
cd chemiscope
npm install
npm start

# navigate to localhost:8080
```

## Building the code to use it in other projects

```bash
git clone https://github.com/lab-cosmo/chemiscope
cd chemiscope
npm install
npm run build

# Include dist/chemiscope.min.js or dist/molecule-viewer.min.js
# in your own web page
```

See [app/] or the [documentation](https://chemiscope.org/docs/embedding.html)
for a examples of how to create a webpage using chemiscope.

## License and contributions

If you are interested in contributing to chemiscope, please have a look at our
[contribution guidelines](Contributing.md)

Chemiscope itself is distributed under the 3-Clauses BSD license. By
contributing to this repository, you agree to distribute your contributions
under the same license.
# Contributing to chemfiles

:tada: First off, thanks for taking the time to contribute to chemiscope! :tada:

If you want to contribute but feel a bit lost, do not hesitate to contact us and
ask your questions! We will happily mentor you through your first contributions.

## Area of contributions

The first and best way to contribute to chemiscope is to use it and advertise it
to other potential users. Other than that, you can help with:

-   documentation: correcting typos, making various documentation clearer;
-   bug fixes and improvements to existing code;
-   and many more …

All these contributions are very welcome. We accept contributions via Github
pull request (have a look [here][pr] for Github model of pull request). If you
want to work on the code and pick something easy to get started, have a look at
the [first good issues][easy-issues].

## Bug reports and feature requests

Bug and feature requests should be reported as [Github issue][issue]. For bugs,
you should provide information so that we can reproduce it: what did you try?
What did you expect? What happened instead? Please provide any useful code
snippet or input file with your bug report.

If you want to add a new feature to chemiscope, please create an [issue] so that
we can discuss it, and you have more chances to see your changes incorporated.

### Code contribution check-list

Every item in this list is explained in the next section

-   [ ] Fork chemiscope;
-   [ ] Create a local branch;
-   [ ] Add code / correct typos / ...;
-   [ ] Check that the code passes lint checks;
-   [ ] Push to Github;
-   [ ] Create a Pull Request;
-   [ ] Discuss your changes with the reviewers;
-   [ ] Have your code merged
-   [ ] Celebrate! :tada: :cake: :tada:

### Contribution tutorial

In this small tutorial, you should replace `<angle brackets>` as needed. If
anything is unclear, please ask for clarifications! There are no dumb questions.

---

Start by [forking chemiscope][fork], and then clone and start a development
server for your fork by running:

```bash
git clone https://github.com/<YOUR USERNAME>/chemiscope
cd chemiscope
npm install
npm start
```

Then create a new branch for your changes

```bash
git checkout -b <new-branch>
```

Implement your changes, including the documentation (in `docs`) for new code.
You can then navigate to `http://localhost:8080` to interact with the code and
check that everything works as expected. It is always a good idea to check that
your code is working with multiple browsers (Chrome, Firefox, Edge, Safari).

Run tests and lints (we use [eslint] and [prettier] to ensure a consistent
coding style):

```bash
npm test
```

We suggest that you configure your code editor to automatically re-format the
code when you save the files. There are prettier plugins for most editors, see
the "Editor Support" section of the [prettier] website. If you want to manually
format your files, you can use

```bash
npx prettier --write <path/to/the/files>
```

Finally, you can push your code to Github, and create a [Pull Request][pr] to
the `lab-cosmo/chemiscope` repository.

```bash
git commit  # ask for help if you don't know how to use git
git push -u origin <new-branch>
```

[pr]: https://help.github.com/articles/using-pull-requests/
[easy-issues]: https://github.com/lab-cosmo/chemiscope/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22
[fork]: https://help.github.com/articles/fork-a-repo/
[issue]: https://github.com/lab-cosmo/chemiscope/issues/new
[eslint]: https://eslint.org/
[prettier]: https://prettier.io/
---
title: 'Chemiscope: interactive structure-property explorer for materials and molecules'
tags:
    - TypeScript
    - JavaScript
    - chemistry
    - material science
    - machine learning
authors:
    - name: Guillaume Fraux
      orcid: 0000-0003-4824-6512
      affiliation: 1
    - name: Rose K. Cersonsky
      orcid: 0000-0003-4515-3441
      affiliation: 1
    - name: Michele Ceriotti
      orcid: 0000-0003-2571-2832
      affiliation: 1
affiliations:
    - name: Laboratory of Computational Science and Modeling, IMX, École Polytechnique Fédérale de Lausanne, 1015 Lausanne, Switzerland
      index: 1
date: 30 January 2020
bibliography: paper.bib
---

# Summary

The number of materials or molecules that can be created by combining different
chemical elements in various proportions and spatial arrangements is enormous.
Computational chemistry can be used to generate databases containing billions of
potential structures [@Ruddigkeit2012], and predict some of the associated
properties [@Montavon2013; @Ramakrishnan2014]. Unfortunately, the very large
number of structures makes exploring such database — to understand
structure-property relations or find the _best_ structure for a given
application — a daunting task. In recent years, multiple molecular
_representations_ [@Behler2007; @Bartok2013; @Willatt2019] have been developed
to compute structural similarities between materials or molecules, incorporating
physically-relevant information and symmetries. The features associated with
these representations can be used for unsupervised machine learning
applications, such as clustering or classification of the different structures,
and high-throughput screening of database for specific properties [@Maier2007;
@De2017; @Hautier2019]. Unfortunately, the dimensionality of these features (as
well as most of other descriptors used in chemical and materials informatics) is
very high, which makes the resulting classifications, clustering or mapping very
hard to visualize. Dimensionality reduction algorithms [@Schlkopf1998;
@Ceriotti2011; @McInnes2018] can reduce the number of relevant dimensions to a
handful, creating 2D or 3D maps of the full database.

![The Qm7b database [@Montavon2013] visualized with chemiscope](screenshot.png)

Chemiscope is a graphical tool for the interactive exploration of materials and
molecular databases, correlating local and global structural descriptors with
the physical properties of the different systems. The interface consists of
two panels. The left panel displays a 2D or 3D scatter plot, in which each
point corresponds to a chemical entity. The axes, color, and style of each point
can be set to represent a property or a structural descriptor to visualize
structure-property relations directly. Structural descriptors are not computed
directly by chemiscope, but must be obtained from one of the many codes
implementing general-purpose atomic representation [@librascal; @QUIP] or more specialized descriptors. Since the most common
descriptors can be very high dimensional, it can be convenient to apply a
dimensionality reduction algorithm that maps them to a lower-dimensional space
for easier visualization. For example the sketch-map algorithm [@Ceriotti2011]
was used with the Smooth Overlap of Atomic Positions representation [@Bartok2013] to
generate the visualization in Figure 1. The right panel displays the
three-dimensional structure of the chemical entities, possibly including
periodic repetition for crystals. Visualizing the chemical structure can help
in finding an intuitive rationalization of the layout of the dataset and the
structure-property relations.

Whereas similar tools [@Gong2013; @Gutlein2014; @Probst2017; @ISV] only allow
visualizing maps and structures in which each data point corresponds to a
molecule, or a crystal structure, a distinctive feature of chemiscope is the
possibility of visualizing maps in which points correspond to atom-centred
environments. This is useful, for instance, to rationalize the relationship
between structure and atomic properties such as nuclear chemical shieldings
(Figure 2). This is also useful as a diagnostic tool for the many
machine-learning schemes that decompose properties into atom-centred
contributions [@Behler2007; @Bartok2010].

![Database of chemical shieldings [@Paruzzo2018] in chemiscope demonstrating the use of a 3D plot and highlighting of atomic environments](./screenshot-3d.png)

Chemiscope took strong inspiration from a previous similar graphical software,
the interactive sketch-map visualizer [@ISV]. This previous software was used in
multiple research publication, related to the exploration of large-scale
databases, and the mapping of structure-property relationships [@De2016;
@De2017; @Musil2018].

# Implementation

Chemiscope is implemented using the web platform: HTML5, CSS and WebGL to
display graphical elements, and TypeScript (compiled to JavaScript) for
interactivity. It uses [Plotly.js](https://plot.ly/javascript/) to render and
animate 2D and 3D plots; and the JavaScript version of [Jmol](http://jmol.org/)
to display atomic structures. The visualization is fast enough to be used with
datasets containing up to a million points, reacting to user input within a few
hundred milliseconds in the default 2D mode. More elaborate visualizations are
slower, while still handling 100k points easily.

The use of web technologies makes chemiscope usable from different operating
systems without the need to develop, maintain and package the code for each
operating system. It also means that we can provide an online service at
http://chemiscope.org that allows users to visualize their own dataset without any
local installation. Chemiscope is implemented as a library of re-usable
components linked together via callbacks. This makes it easy to modify the
default interface to generate more elaborate visualizations, for example,
displaying multiple maps generated with different parameters of a dimensionality
reduction algorithm. Chemiscope can also be distributed in a standalone mode,
where the code and a predefined dataset are merged together as a single HTML
file. This standalone mode is useful for archival purposes, for example as
supplementary information for a published article and for use in corporate
environments with sensitive datasets.

# Acknowledgements

The development of chemiscope have been funded by the [NCCR
MARVEL](http://nccr-marvel.ch/), the [MAX](http://max-centre.eu/) European
centre of excellence, and the European Research Council (Horizon 2020 grant
agreement no. 677013-HBMAP).

# References
# Python helpers for chemiscope

This package contains Python code to help generate input files for the
[chemiscope](https://chemiscope.org) default visualizer, and integrate
chemiscope with jupyter notebooks.

## Installation

You should use pip to install this package:

```bash
pip install chemiscope
```

This installs both a `chemiscope-input` command line tool, and the `chemiscope`
package.

## Usage

To create a new chemiscope input file:

```python
import chemiscope
import ase.io

# read frames using ase
frames = ase.io.read("structures.xyz", ":")

# add additional properties to display
properties = {
    "<property name>": {
        target: "atom",
        values: [3, 4, 2, 8, 9, 10],
    }
}

chemiscope.write_input("my-input.json.gz", frames=frames, properties=properties)
```

To display a chemiscope widget inside a jupyter notebook:

```python
import chemiscope
import ase.io

# read frames using ase
frames = ase.io.read("structures.xyz", ":")

# add additional properties to display
properties = {
    "<property name>": [3, 4, 2, 8, 9, 10],
}

chemiscope.show(frames=frames, properties=properties)
```
Chemiscope: interactive structure/property explorer for materials and molecules
===============================================================================

Welcome to the documentation of the `chemiscope`_ visualization tool, an
interactive structure/property explorer for materials and molecules. The goal of
chemiscope is to provide interactive exploration of large databases of materials
and molecules and help researchers to find structure-properties correlations
inside such databases. The screenshot below shows an example of such database
being visualized with chemiscope. The first part of this documentation describes
the default interface of chemiscope and how to use it with your own database.

.. figure:: img/screenshot.png
    :align: center

    Screenshot of the `Qm7b`_ database visualized in the default chemiscope viewer


Chemiscope is built around re-usable components, that can be arranged in
different manners to create visualization adapted to different kinds of data. The
second part of this documentation explains how to build the code and use it in
your own website to create new interfaces.

Citing chemiscope
-----------------

Chemiscope is distributed under an open-source license, and you are welcome to
use it and incorporate it into your own research and software projects.
If you find it useful, we would appreciate a citation to the chemiscope
`paper`_:

G. Fraux, R. K. Cersonsky, M. Ceriotti, *Chemiscope: Interactive
Structure-Property Explorer for Materials and Molecules.* **Journal of Open
Source Software** 5 (51), 2117 (2020)

If you incorporate chemiscope components into a software project, a link
back to the `chemiscope`_ homepage is the preferred form of acknowledgement.


What's in this documentation?
-----------------------------

.. toctree::
    :maxdepth: 2

    tutorial/index
    embedding


.. _chemiscope: https://chemiscope.org
.. _paper: https://doi.org/10.21105/joss.02117
.. _Qm7b: https://doi.org/10.1088/1367-2630/15/9/095003
Chemiscope as a library
=======================

It is possible to use chemiscope as a software library when writing your own
web-based interface. This page document how to get the library and give a few
usage examples.

.. note::
    You should also look at the `API documentation <api/index.html>`_ for a list
    of all public classes, interfaces and functions in chemiscope.

Dependencies
^^^^^^^^^^^^

Chemiscope relies on different external dependencies that you should load in all
the HTML pages using it. You can serve these from your own web server, or use a
CDN to deliver them.

- `Bootstrap <https://getbootstrap.com/>`_ for HTML styling and basic UI;
- Bootstrap and 3Dmol.js rely on the ubiquitous `JQuery <https://jquery.com/>`_

Getting and building the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-built version
-----------------

The easiest way to do so is to download the latest release from `the release
page <https://github.com/lab-cosmo/chemiscope/releases>`_ on GitHub. The main
file is ``chemiscope.min.js``, containing the code required to create the
default visualizer. This file exports a single global object ``Chemiscope``,
which contains references to

- `Chemiscope.DefaultVisualizer <DefaultVisualizer_>`_
- `Chemiscope.addWarningHandler <addWarningHandler_>`_
- `Chemiscope.MetadataPanel <MetadataPanel_>`_
- `Chemiscope.ViewersGrid <ViewersGrid_>`_
- `Chemiscope.PropertiesMap <PropertiesMap_>`_
- `Chemiscope.EnvironmentInfo <EnvironmentInfo_>`_
- `Chemiscope.EnvironmentIndexer <EnvironmentIndexer_>`_

Partial builds are also available, in particular ``molecule-viewer.min.js``
which only contains code related to `Chemiscope.MoleculeViewer
<MoleculeViewer_>`_, making the minified JavaScript file much smaller. Other
partial builds containing only part of chemiscope can be added upon request.

.. _DefaultVisualizer: api/classes/DefaultVisualizer.html
.. _addWarningHandler: api/index.html#addWarningHandler
.. _ViewersGrid: api/classes/ViewersGrid.html
.. _PropertiesMap: api/classes/PropertiesMap.html
.. _EnvironmentInfo: api/classes/EnvironmentInfo.html
.. _MetadataPanel: api/classes/MetadataPanel.html
.. _EnvironmentIndexer: api/classes/EnvironmentIndexer.html
.. _MoleculeViewer: api/classes/MoleculeViewer.html

Build from sources
------------------

Chemiscope is written in TypeScript, a statically typed language which compiles
to JavaScript. It uses the standard JavaScript ecosystem tools for dependency
management, ``npm``. To build chemiscope from sources, you will first need to
get the sources, either as an archive from `the release page
<https://github.com/lab-cosmo/chemiscope/releases>`_, or using git

.. code-block:: bash

    git clone https://github.com/lab-cosmo/chemiscope

You will also need `node.js <https://nodejs.org/en/>`_ and `npm
<https://docs.npmjs.com/cli/npm>`_, which you can install with your favorite
package manager.

.. code-block:: bash

    cd chemiscope
    npm install
    npm run build

This should create a ``dist`` directory, containing all the minified JavaScript
libraries.

Usage inside a JavaScript project
---------------------------------

If you already have a JavaScript project using ``npm`` or ``yarn``, you can use
the version of Chemiscope available on npm. Add it to your project with

.. code-block:: bash

    npm install chemiscope

And import the library inside your own code using

.. code-block:: js

    const Chemiscope = require("chemiscope");

If your are using TypeScript, definition files are also provided with the npm
package, and should give you auto-completion, inline documentation and interface
checking.

Usage example
^^^^^^^^^^^^^

Below is the minimal HTML and JavaScript code needed to load and use the default
chemiscope interface. Additional code would be required to load the dataset,
using JSON files or directly.

.. literalinclude:: embedded-example.html
    :language: html
Jupyter notebooks
=================

.. autofunction:: chemiscope.show
Sharing datasets with collaborators
===================================

Once you have converted your data in the :ref:`format used by chemiscope
<input-format>`, you might want to share it with collaborators. There are
multiple ways to do this, we'll go over them in this section.

Online visualizer at chemiscope.org
-----------------------------------

Uploading datasets
^^^^^^^^^^^^^^^^^^

The simplest way to share chemiscope dataset is to send the corresponding file
(``my-dataset.json`` or ``my-dataset.json.gz``) to your collaborators.

They can then load this file from the main chemiscope website:
https://chemiscope.org/. Here they can use the *Load/Save* menu to upload the
dataset to the website.

Loading from an URL
^^^^^^^^^^^^^^^^^^^

You can also host your dataset in any publicly accessible web server such as
your home page. Given a file hosted at
``https://university.edu/~myself/dataset.json``, you can then load this file
directly in chemiscope by going to
``https://chemiscope.org/?load=https://university.edu/~myself/dataset.json``.
This allow you to create links to directly open a given dataset in the main
chemiscope website, loading such dataset from your own webpage.

In general, you can set the ``load`` GET parameter on ``https://chemiscope.org``
to any url-encoded URL, and chemiscope will try to load the dataset from this
URL.

Saving visualization state
^^^^^^^^^^^^^^^^^^^^^^^^^^

You can save the current visualization settings from the website using the
*Load/Save* menu. This will allow you to download a JSON file, which you can use
later to reset visualization state.

If you are loading a file from an URL by creating a link that looks like
``https://chemiscope.org/?load=https://university.edu/dataset.json.gz``, you can
use the ``settings`` GET parameter to specify the URL of saved settings:
``https://chemiscope.org/?load=https://university.edu/dataset.json.gz&settings=https://university.edu/settings.json``

Standalone offline visualizer
-----------------------------

There are some cases where you do not want to depend on an online tool when
sharing your dataset, such as scientific article supplementation information.
For these use cases, a standalone, offline visualizer exists that uses the same
input file format as the default interface. You can download the latest version
of the standalone viewer at
:download:`https://chemiscope.org/chemiscope_standalone.html`.

This file contains all the required HTML and JavaScript code for chemiscope. You
can then add your own dataset by adding the corresponding JSON file at the end
of the ``chemiscope_standalone.html`` file.

**WARNING:** Only JSON, not compressed JSON (``.json.gz``) files are supported
with the standalone HTML visualizer.

.. code-block:: bash

    cat chemiscope_standalone.html my-dataset.json > my-dataset.html

You can then share the ``my-dataset.html`` file with others, who can open it in
any web browser.

To re-build the ``chemiscope_standalone.html`` file from sources, please follow
the steps below:

.. code-block:: bash

    git clone https://github.com/lab-cosmo/chemiscope
    cd chemiscope
    npm install
    npm run build
    python3 ./utils/generate_standalone.py
Different panels and settings
=============================

The default chemiscope visualizer is organized in three main panels: the map,
the structure viewer and the environment information display. Additionally,
clicking on the dataset title (on top of the map) will display some metadata
about the dataset (description, authors, references). This section will
present each one, as well as the main settings accessible to customize the
display.

The map is a 2D or 3D scatter plot showing properties for all the environments
in the dataset. You can set which properties (structural or physical) should be
used a the x, y, and potentially z axis; as well as for color and size of the
points. Additionally, properties which have string values (an not numeric
values) can be used as category data to set the symbols used for the points. To
open the settings modal window, click on the hamburger menu (the ☰ symbol) on
the left of the dataset title.

.. figure:: ../img/map.png
    :width: 80 %

    The map panel in 2D mode and the related settings

The structure panel is a 3D molecular viewer based on `3Dmol.js`_. The settings are
accessible through the hamburger menu (☰) on the right of the viewer. The
settings are grouped into **representation** (how is the molecule rendered);
**supercell** (how many copies of the unit cell to display); **environments**
(how atom-centered environments are displayed); **camera** (reset the camera in
along one of the given axis); and **trajectory** (playback related settings).

.. figure:: ../img/structure.png
    :width: 80 %

    The structure panel and related settings

Finally, the environments information panel features sliders and text input to
allow for an easy selection of the environment of interest. The play button on
the left of the sliders activates the trajectory playback, looping over the
structures in the datasets or the atoms in a structure. By clicking on the
labels at the top (*structure XXX* and*atom XXX*), one can hide or show the
full property tables. These tables show all properties in the dataset for the
currently selected environment.

.. figure:: ../img/info.png
    :width: 40 %

    The environment information panel fully expanded

.. _3Dmol.js: https://3dmol.csb.pitt.edu/
.. _input-reference:

Input file reference
====================

If you can not or do not want to use the ``chemiscope`` python package to create
your input files, you can also directly write the JSON file conforming to the
schema described here. The input file follows closely the `Dataset`_ typescript
interface used in the library. Using a pseudo-JSON format, the file should
contains the following fields and values:

.. code-block:: javascript

    {
        // metadata of the dataset. `description`, `authors` and `references`
        // will be rendered as markdown.
        "meta": {
            // the name of the dataset
            "name": "this is my name"
            // description of the dataset, OPTIONAL
            "description": "This contains data from ..."
            // authors of the dataset, OPTIONAL
            "authors": ["John Doe", "Mr Green, green@example.com"],
            // references for the dataset, OPTIONAL
            "references": [
                "'A new molecular construction', Journal of Random Words 19 (1923) pp 3333, DOI: 10.0000/0001100",
                "'nice website' http://example.com",
            ],
        },

        // list of properties in this dataset
        "properties": {
            // Each property have at least a name, a target and some values.
            // Optional entries for the units and descriptions can also be added.
            <name>: {
                // the property target: is it defined per atom or for the full
                // structures
                "target": "atom" | "structure",
                // values of the properties can either be numbers or strings.
                // string properties are assumed to represent categories of
                // data.
                "values": [1, 2, 3, ...] | ["first", "second", "first", ...],

                // OPTIONAL: units of the property' value
                "units": "A/fs^2",
                // OPTIONAL: free-form description of the property as a string
                "description": "acceleration of the atoms in the structure ...",
            }
        }

        // list of structures in this dataset
        "structures": [
            {
                // number of atoms in the structure
                "size": 42,
                // names of the atoms in the structure
                "names": ["H", "O", "C", "C", ...],
                // x cartesian coordinate of all the atoms, in Angstroms
                "x": [0, 1.5, 5.2, ...],
                // y cartesian coordinate of all the atoms, in Angstroms
                "y": [5.7, 7, -2.4, ...],
                // z cartesian coordinate of all the atoms, in Angstroms
                "z": [8.1, 2.9, -1.3, ...],
                // OPTIONAL: unit cell of the system, if any.
                //
                // This should be given as [ax ay az bx by bz cx cy cz], where
                // a, b, and c are the unit cell vectors. All values are
                // expressed in Angstroms.
                "cell": [10, 0, 0, 0, 10, 0, 0, 0, 10],
            },
            // other structures as needed
            ...
        ],

        // OPTIONAL: atom-centered environments descriptions
        //
        // If present, there should be one environment for each atom in each
        // structure.
        "environments": [
            {
                // index of the structure in the above structures list
                "structure": 0,
                // index of the central atom in structures
                "center": 8,
                // spherical cutoff radius, expressed in Angstroms
                "cutoff": 3.5,
            },
            // more environments
            ...
        ]

        // OPTIONAL: setting for each panel
        //
        // Adding these values allow to setup how a given dataset should be
        // visualized in chemiscope.
        //
        // Each value inside the settings group is optional
        "settings": {
            // settings related to the map
            "map": {
                // x axis settings
                "x": {
                    // name of the property to use for this axis, this must be
                    // one of the key from the root `properties` table.
                    "property": "<name>",
                    // should the axis use linear or logarithmic scaling
                    "scale": "linear" | "log",
                    // lower bound of the axis
                    "min": -0.23,
                    // upper bound of the axis
                    "max": 1.42,
                },
                // y axis setting, using the the same keys as x axis setting
                "y": {
                    // ...
                },
                // z axis setting, using the the same keys as x axis setting
                "z": {
                    // property can be set to an empty string to get a 2D map
                    "property": "",
                    // ...
                },
                // name of the property to use for markers symbols, this must be
                // one of the key from the root `properties` table. The
                // associated property should have string values
                "symbol": "<name>",
                // point color setting, using the the same keys as x axis setting
                "color": {
                    // property can be set to an empty string for uniform color
                    "property": "",
                    // ...
                },
                // Color palette to use, default to 'inferno'
                "palette": "cividis",
                // settings related to the markers sizes
                "size": {
                    // scaling factor for the axis, between 1 and 100
                    "factor": 55,
                    // mode to scale the markers with respect to the properties
                      // `constant`: all markers are same size, scaled by `factor`
                      // `linear`: markers are directly proportional to the property
                      // `log`: markers are proportional to the logarithm of the property
                      // `sqrt`: markers are proportional to the square root of the property
                      // `inverse`: markers are inversely proportional to the property
                    "mode": "constant" | "linear" | "log" | "sqrt | "inverse"",
                    // name of the property to use for the markers size, this
                    // must be one of the key from the root `properties` table.
                    "property": "<name>",
                    // if false, markers scale from smallest to largest property value
                    // if true, marker scale from largest to smallest property value
                    // in the case of `inverse` scaling, this is reversed.
                    "reverse": false | true,
                },
            },
            // Settings related to the structure viewers grid. This is an array
            // containing the settings for each separate viewer
            "structure": [
                {
                    // show bonds between atoms
                    "bonds": true,
                    //use space filling representation
                    "spaceFilling": false,
                    // show atoms labels
                    "atomLabels": false,
                    // show unit cell information and lines
                    "unitCell": false,
                    // displayed unit cell as a packed cell
                    "packedCell": false,
                    // number of repetitions in the `a/b/c` direction for the supercell
                    "supercell": [2, 2, 3],
                    // make the molecule spin
                    "rotation": false,
                    // which axis system to use
                    "axes": "none" | "xyz" | "abc",
                    // keep the orientation constant when loading a new structure
                    "keepOrientation": false,
                    // options related to atom-centered environments
                    "environments": {
                        // should we display environments & environments options
                        "activated": true,
                        // automatically center the environment when loading it
                        "center": false,
                        // the cutoff value for spherical environments
                        "cutoff": 3.5
                        // which style for atoms not in the environment
                        "bgStyle": "licorice" | "ball-stick" | "hide",
                        // which colors for atoms not in the environment
                        "bgColor": "grey" | "CPK",
                    };
                },
                // ...
            ]
            // List of environments to display (up to 9). These environments
            // will be shown in the structure viewer grid and indicated on
            // the map.
            //
            // This list should containg 0-based indexes of the environment in
            // the root "environments" object; or of the structure in the root
            // "environments" if no environments are present.
            //
            // If both this list and the "structure" settings list above are
            // present, they should have the same size and will be used together
            // (first element of "structure" setting used for the first "pinned"
            // value; and so on).
            //
            // This defaults to [0], i.e. showing only the first
            // environment/structure.
            "pinned": [
                33, 67, 12, 0,
            ]
        }
    }

.. _Dataset: api/interfaces/main.dataset.html

Introduction to structural properties
=====================================

Before we get started, we will introduce a few concepts that underlie the
concept and the usage of chemiscope. Chemiscope is designed to help navigating
*structure-property maps*, i.e. 2D or 3D representations of a set of
atomic scale entities that reflect how structure influences materials
properties.

Chemiscope can work with two kinds of entities: full structures, or
atom-centred environments. A structure consists in a set of atoms, possibly
representing the periodic repeat unit of an infinite structure. An
environment consists in a set of atoms that surround a central atom,
In both cases, these entities are fully defined by the position and nature
of the atoms present in the structure, or in the neighborhood of the
environment center.

For each structure or environment, one may have computed *properties*,
e.g. the cohesive energy of a molecule, or the NMR chemical shielding of
a nucleus, or *structural representations*, i.e. functions of the
spatial arrangement of the atoms that incorporate some fundamental
symmetries to achieve a description of the structure that is as complete
as possible, yet concise. Examples of such representations are for instance
`atom density representationis <soap>`_ or `Behler-Parrinello
symmetry functions <Behler-Parrinello>`_. These representations are usually
high-dimensional vectors, hard to visualize and interpret. For this reason, one
usually applies a dimensionality reduction algorithm, such as `PCA`_, `sketch-map`_,
*etc.*   The interpretation of the resulting  will differ depending on both the
descriptor used to represent the structures or environments and the
dimensionality reduction algorithm applied.

Chemiscope simplifies visualizing the correlations between structural
representations and properties associated with structures and environments,
by representing in an interactive fashion these atomic-scale entities as points
on a map, and by associating these points with an explicit, 3D visualization
of the structure of the material or molecule.

.. figure:: ../img/mol-to-map.*
    :width: 65 %

    Illustration of the process used to create structural properties from a
    molecule.

Chemiscope is completly agnostic with respect to how properties and structural
representations are generated, and do not provide any facilities to generate them.
In the rest of this document, we will refer to properties describing
the structure of an environment or structure as *structural properties*
and other associated properties associated (such as energy, density, ...) as
*physical properties*.

.. _soap: https://doi.org/10.1063/1.5090481
.. _Behler-Parrinello: https://doi.org/10.1103/physrevlett.98.146401
.. _PCA: https://en.wikipedia.org/wiki/Principal_component_analysis
.. _sketch-map: https://doi.org/10.1073/pnas.1108486108
User tutorial
=============

This tutorial will present all then different ways one can use the standard
chemiscope visualization with your own dataset. This section starts by
introducing the concept of structural and physical properties, before describing
how to use the different panels in the standard visualization. It continues by
presenting how you can generate a chemiscope input file to load on
https://chemiscope.org, as well as within a standalone HTML viewer which does
not require internet connectivity. Finally, we present the chemiscope jupyter
extension, which can be used to explore a dataset directly inside a jupyter
notebook.

.. _chemiscope: https://chemiscope.org

.. toctree::
    :maxdepth: 2

    properties
    panels
    input
    input-reference
    sharing
    jupyter
.. _input-format:

Creating chemiscope input files
===============================

When using the default chemiscope interface, all the structures and properties
in a dataset are loaded from a single JSON file. This sections describe how to
generate such JSON file, either using a pre-existing python script that does
most of the work for you, or by writing the JSON file directly. Since the
resulting JSON file can be quite large and thus harder to share with
collaborators, the default chemiscope interface also allows to load JSON files
compressed with gzip.

Tools able to create chemiscope input
-------------------------------------

``chemiscope`` Python module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to create a JSON input file is to use the ``chemiscope`` Python
module. Install the package with ``pip install chemiscope``, and use
:py:func:`chemiscope.write_input` or :py:func:`chemiscope.create_input` in your
own script to generate the JSON file.

If all the properties you want to include into chemiscope are already stored in
a file `ase`_ can read, the ``chemiscope`` python package also install a
`chemiscope-input <chemiscope-input-cli_>`_ command line script.

Note that chemiscope does not compute structural representations or
dimensionality reduction, and you need to do this yourself or use another
package such as ASAP.

``ASAP``
^^^^^^^^

The `ASAP`_ structural analysis package is another tool that can directly
generate an output in chemiscope format.

``chemiscope`` functions reference
----------------------------------

.. autofunction:: chemiscope.write_input

.. autofunction:: chemiscope.create_input

.. autofunction:: chemiscope.all_atomic_environments

.. autofunction:: chemiscope.librascal_atomic_environments

.. _ase: https://wiki.fysik.dtu.dk/ase/index.html
.. _ASAP: https://github.com/BingqingCheng/ASAP


.. _chemiscope-input-cli:

``chemiscope-input`` command line interface
-------------------------------------------

.. sphinx_argparse_cli::
   :module: chemiscope.main
   :func: _chemiscope_input_parser
   :prog: chemiscope-input
   :title:
   :group_title_prefix:
