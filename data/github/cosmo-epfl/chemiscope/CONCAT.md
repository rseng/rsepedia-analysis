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
