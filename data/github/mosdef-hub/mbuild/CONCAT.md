## mBuild: a hierarchical, component based molecule builder

[![Gitter chat](https://badges.gitter.im/mosdef-hub/gitter.svg)](https://gitter.im/mosdef-hub/Lobby)
[![AZP Build Status](https://dev.azure.com/mosdef/mbuild/_apis/build/status/mosdef-hub.mbuild?branchName=master)](https://dev.azure.com/mosdef/mbuild/_build/latest?definitionId=1&branchName=master)
[![Anaconda Badge](https://anaconda.org/conda-forge/mbuild/badges/version.svg)](https://anaconda.org/conda-forge/mbuild)
[![codecov](https://codecov.io/gh/mosdef-hub/mbuild/branch/master/graph/badge.svg)](https://codecov.io/gh/mosdef-hub/mbuild)

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular dynamics simulations.

* mBuild is designed to minimize or even eliminate the need to explicitly translate and
  orient components when building systems: you simply tell it to connect two
  pieces!
* mBuild keeps track of the system's topology so you don't have to
  worry about manually defining bonds when constructing chemically bonded
  structures from smaller components.

To learn more, get started or contribute, check out our [website](http://mbuild.mosdef.org).

#### mBuild within the MoSDeF Ecosystem
<p align="center">
  <img src="docs/images/mosdef.svg?raw=true" alt="mBuild within the MoSDeF Ecosystem" width="500" height="500"/>
</p>


The `mBuild` package is part of the [Molecular Simulation Design Framework (MoSDeF) project](http://mosdef.org/).
Libraries in the MoSDeF ecosystem are designed to provide utilities neccessary to streamline
a researcher's simulation workflow. When setting up simulation studies,
we also recommend users to follow the [TRUE](https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938)
(Transparent, Reproducible, Usable-by-others, and Extensible) standard, which is a set of common
practices meant to improve the reproducibility of computational simulation research.

#### Quick Start with Docker
To use `mbuild` in a jupyter-notebook that runs from a docker container with all the dependencies installed use the following command:

```sh
$ docker pull mosdef/mbuild:latest
$ docker run -it --name mbuild -p 8888:8888 mosdef/mbuild:latest su anaconda -s\
      /bin/sh -l -c "jupyter-notebook --no-browser --ip="0.0.0.0" --notebook-dir\
      /home/anaconda/mbuild-notebooks
```

Alternatively, you can also start a Bourne shell directly:
```sh
$ docker run -it --name mbuild mosdef/mbuild:latest
```

To learn more about using `mBuild` with docker, please refer to the documentation [here](https://mbuild.mosdef.org/en/latest/docker.html).

#### Tutorials

*Interactive tutorials can be found here:*

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/mosdef-hub/mbuild_tutorials/master)

Components in dashed boxes are drawn by hand using, e.g.,
[Avogadro](https://avogadro.cc/) or generated elsewhere. Each
component is wrapped as a simple python class with user defined attachment
sites, or ports. That's the "hard" part! Now mBuild can do the rest. Each component
further down the hierarchy is, again, a simple python class that describes
which piece should connect to which piece.

Ultimately, complex structures can be created with just a line or two
of code. Additionally, this approach seamlessly exposes tunable parameters within
the hierarchy so you can actually create whole families of structures simply
by adjusting a variable:

```python
import mbuild as mb
from mbuild.examples import PMPCLayer

pattern = mb.Random2DPattern(20)  # A random arrangement of 20 pieces on a 2D surface.
pmpc_layer = PMPCLayer(chain_length=20, pattern=pattern, tile_x=3, tile_y=2)
```

![Zwitterionic brushes on beta-cristobalite substrate](docs/images/pmpc.png)
#### Community Recipes
Use case-specific systems can be generated via mBuild recipes.
Some users have graciously contributed recipes for particular systems, including:

* [Graphene slit pores](https://github.com/rmatsum836/Pore-Builder)
* [Nanodroplets on graphene](https://github.com/ftiet/droplet-builder)
* [Coarse-grained DNA](https://github.com/zijiewu3/mbuild_ONA)
* [Lipid bilayers](https://github.com/uppittu11/mbuild_bilayer)

#### Citing mBuild

If you use this package, please cite [our paper](http://dx.doi.org/10.1007/978-981-10-1128-3_5
). The BibTeX reference is
```
@article{Klein2016mBuild,
      author = "Klein, Christoph and Sallai, János and Jones, Trevor J. and Iacovella, Christopher R. and McCabe, Clare and Cummings, Peter T.",
      title = "A Hierarchical, Component Based Approach to Screening Properties of Soft Matter",
      booktitle = "Foundations of Molecular Modeling and Simulation",
      series = "Molecular Modeling and Simulation: Applications and Perspectives",
      year = "2016",
      doi = "http://dx.doi.org/10.1007/978-981-10-1128-3_5"
}
```

#### [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

This material is based upon work supported by the National Science Foundation under grants NSF CBET-1028374 and NSF ACI-1047828. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
# Change Log #

All big and breaking changes for [mBuild](https://mosdef-hub.github.io/mbuild/) will be recorded here.
This project adheres to [Semantic Versioning](http://semver.org/).

## 0.10.4 (2019-12-28)

### Features
- Add `wrap_coords` function for wrapping particles into a box (#643)

### Bug Fixes
- Don't populate empty lists for hoomd snapshot (#647)

### Maintenance
- Add `nbval` as requirement for tests in conda recpe (#656)
- Remove `six` and some Python 2 legacy code (#641)
- Remove examples (from this repository) (#658)
- Fix travis by not installing md5sha1sum (#660)
- Add CHANGELOG (#632)

## 0.10.3 (2019-10-27)

### Features
- Add Cassandra MCF writer (#636)
- Add HOOMD snapshot converter (#622)
- Generate `mb.Compound` from Parmed or MDTraj residues (#526)
- Add `**kwargs` for `write_gsd` (#653)

### Bug Fixes
- Fix unique naming problem in `to_networkx` (#583)
- Default the mBuild PAR-writer to use `IMPROPER` (#648)
- Fix ghost ports attached to removed compound (#593)
- Remove duplicate NP particles in TNP example (#625)
- Only import MCF writer if `networkx` is installed (#654)

### Maintenance
- Add LGTM (#616)
- Add Azure CI testing (#615, #617, #618, #630, #637, #638)
- Update Travis CI (#649)
- Update conda recipe to include nbval (#611, #656)
- Pin `nglview` to 2.7+ (#651, #655)
- Reduce length of some unit tests (#621)
- Rename `AmorphousSilica` to be more descriptive (#630)

## 0.10.1 (2019-10-9)
This is a bugfix release to resolve a potential issue with a foyer dependency with mBuild 0.10.0.

### Bugfixes
* Remove a `foyer` import that was producing a circular dependency (#610)

## 0.10.0 (2019-10-8)
### Breaking Change
- Officially drop Python 2.7 (#573)

### Features
- Load mBuild compounds from existing ParmEd and MDTraj objects (#561)
- Convert mBuild compound to and from JSON format (#581)
- Include testing of notebooks in CI (#590)
- Add NGLView tooltip (#600)
- Initialize mBuild Compounds from SMILES strings (#598)
- Write out parameterized structures to a CHARMM `.par` file (#508)
- Add method to convert to and from Pybel (#555)
- Add POSCAR file writer in an effort to incorporate VASP into mBuild (#468)

### Bug Fixes
- Remove unnecessary Pybel import statement in `mb.load` (#604)
- Change how proxy compounds are named so that `MOL2` files are in the correct format (#605)
- Rename atom names in silica interface example to be compatible with recent mBuild release (#594)
- Improve error handling for Box class (#576)
- Add a `with open` block to manage file open (#585)
- Add `_clone()` function to the Proxy class to properly clone an instance of Proxy (#592)

### Miscellaneous
- Add `compound_pb2.py` file generated from `protoc` compiler to gitignore (#602)
- Improve documentation of various mBuild classes and functions (#577, #578, #579, #580)
- Add additional testing for `foyer_kwargs` (#428)

## 0.9.3 (2019-8-5)
### Breaking Changes
* This is the last release supporting Python 2.7

### Features
* A more descriptive exception is raised when attempting to access a non-existent Port (#544)
* Element capitalization is better enforced in ParmEd conversions (#550)
* The XYZ reader can now act on a passed compound (#567)
* The LAMMPS writer now accurately prints residue IDs (#569)

### Bugfixes
* Visualizing compounds in notebooks no longer returns a duplicate widget (#545)
* Names of custom elements are no longer clobbered during visualization (#563)

### Maintenance
* The image in our gitter link has been updated (#543)
* Some links in tutorials have been corrected (#548)
* Installation documentation has been updated to reflect changes in conda configurations (#558)
* Some other documentation has been refreshed (#560)
* The GSD version is pinned to 1.7 in order to continue Python 2.7 support (#572)

## 0.9.2 (2019-5-27)
### Breaking Changes
* Python 3.5 is no longer officially supported or tested on as part of the development process.

### Features
* mBuild is now tested and packaged on Python 3.7 (#542)

### Maintenance
* MDTraj is no longer pinned to an old version (#542)
* Coveralls is dropped; we have been using codecov for a few months (#542)

## 0.9.1 (2019-5-26)
### Breaking Changes
This is the last release including official support for Python 3.5. It will likely work for some time but mBuild will not be tested on Python 3.5 during development.

### Features
* Residue names can optionally be inferred from compound names in conversion to ParmEd (#475)
* Custom cross-interactions (NBFIXES in ParmEd jargon) can now be written to LAMMPS data files (#456)
* The LAMMPS writer now prints helpful comments to more verbosely describe which atom types are associated with each potential (#535)

### Maintenance
* Some stylistic changes were made as suggested by various linters (#522)
* Appveyor now tests on Python 3.6 (#520)

### Documentation
* Installation docs were updated to explicitly list supported Python versions (#532)
* A comment pointing to the `glozter` Anaconda channel has been updated to point to `conda-forge` (#534)

### Miscellaneous
* Some tests depending on `foyer` are now properly skipped when it it not installed (#521)
* Some examples were updated in accordance with their new structure as internal recipes (#536, #538)

## 0.9.0 (2019-4-11)
### Features
* A plugin or "recipe" architecture has been added to allow external modules to be imported inside of mBuild (#501)
* Python 3.6 is now explicitly supported and tested (#518)

### Maintenance
* A contributor's guide (#500) and `.github` issue & pull request templates (#498) have been added
* A redundant and unused block of code was removed (#515)
* The `glotzer` and `bioconda` channels, which are now obsolete in this scope, have been dropped (#516)

## 0.8.3 (2019-2-23)
### Breaking Changes
* When writing hoomdxml files, units will now be in kJ/mol & nm instead of kcal/mol & ang, so particle positions will differ by a factor of 10.

### Features
* A `to_networkx` function was added to convert the hierarchy of a compound to a graph (#484)
* Packing functions now use XYZ files while running PACKMOL, bypassing some issues with PDB files (#422)
* When saving hoomdxml files, `auto_scale=True` will scale reference units from max forcefield parameters. (#488)

### Maintenance
* Switched to codecov for code coverage testing (#485)
* Some dependencies accidentally missing in earlier PRs were cleaned up (#493)

### Bugfixes
* `update_coordinates` now behaves well when passed an XYZ file or operating on simple hierarchies (#496)
* Internal conversion from ParmEd structures to HOOMDXML files was improved (#463, see above)

## 0.8.2 (2019-1-8)
### Features
* Special Pair Support (1-4 pair information) to GSD writers (#473)
    * GSD files now include 1-4 special pairs for use in OPLS

### Misc and Bugfixes
* Dependency requirements have been updated (#457)
    * A dependency loop between `Foyer` and `mBuild` has been resolved
* Fixed a bug that prevented Appveyor builds from running (#477)
* Temporary PDB files left behind by packing functions are now properly removed (#471)
    * Packing.py uses temporary files which were previously never closed. This sometimes caused the program to reach the limit of open files for a process set by the OS
* `pytest-ignore-flaky` has been replaced in favor of `xfail` (#471)
* Additonal fixes for PACKMOL input files (#474)
    * Input files are now closed by `mBuild` in order to ensure it can be read by PACKMOL
    * Error reporting is now caught when the subprocess returns an error code
* Microsoft VSCode extraneous files are now ignored by git (#478)

## 0.8.1 (2018-11-28)
### Features
* Packing functions can optionally constrain the rotation of `Compounds` when using `fill_box` (#407)
* Additional lammps datafile support (#412)
    * Add functionality for `atomic`, `charge`, and `molecular` atom_styles
    * Fix `atomic` and `molecular` atom-styles
    * Add optional `atom-style` argument to `save` function
    * Add tests to check for correct `Atoms` format
* A `Compound` can be generated from a SMILES string if the user has [Open Babel](http://openbabel.org/wiki/Main_Page) installed (#430)
* The [website](https://mosdef-hub.github.io/mbuild/) was updated with details how to properly cite mBuild (#421)
* OpenMM can now be used for energy minimization (#416)
* A simple xyz file reader was added (#423)
* Defaults in `Compound.visualize` have been improved (#438)
* mBuild boxes can now store angles (#448)
* mBuild boxes can now be passed to various writers (#448)
* A changelog is now included in the root directory (#454)

### Misc and Bugfixes
* Switched from OpenMM to MDTraj to check for element existence (#408)
* Changed bilayer notebook to use `Compound` methods for object manipulation (#406)
* Added test to ensure that users can provide a custom forcefield XML file when applying a forcefield to a `Compound` (#431)
* An error is now generated if the miniconda MD5 sum does not match when performing tests (#409)
* LAMMPS box values are now written appropriately in Angstroms (#459)
* Coordinates in HOOMDXML and GSD files are now correctly written in the range [L /2 to L] (#452)
* A bug in the ordering of some Bravais angles in non-rectangular lattices has been fixed (#450)

## 0.8.0 (2018-01-18)
### Features
* Improved packing API
	* `fill_box` method now supports user-specified densities (#372)
	* Support for non-cubic boxes (#385)
	* Added edge buffer for pseudo-support of periodic boundaries (#385)
	* Allow users to access PACKMOL raw output (#385)
	* Improve documentation and removed repetitious code (#385)
* Proper support for triclinic lattices (#386)
* More intuitive `Port` behavior when adding/removing bonds
	* `Ports` are now added along the bond vector when bonds are removed
	* `Ports` are removed when using `force_overlap` when attribute `add_bond` is `True` (#390)
### Misc and Bugfixes
* Increased precision of PACKMOL overlap tolerance (#367)
* Increased robustness for `packing.py` argument types (#368)
* Continuous Integration (CI) fixes (#369, #375, #392, #402, #405)
* Documentation updates (#371)
* Combining rules for non-bonded interactions can now be specified when saving compounds (#377)
* Remove ambigious data types for `Box` attributes (#384)
* Fixed changing of basis for non-cubic lattices (#386)
* Fixed issue where `Ports` were not aligned properly along user-specified direction (#390)
* Add support for controlling `Foyer` warning verbosity when saving `Compounds` (#391)
* Added more information to `Port.__repr__` (#396)
* Fixed bug in unit conversion for periodicity in `Compound.from_parmed()` (#401)
* Rounding `Lattice.lattice_point` positions to 0. when below a certain threshold (#404)

## 0.7.1 (2017-06-09)

Merge pull request (#352) from summeraz/use_parmed
Parmed loaders by default, remove Mdtraj dependency

## 0.7.0 (2017-06-09)

Bump version to 0.7.0

## 0.6.1 (2017-02-14)

Tag 0.6.1

## 0.6.0 (2017-02-14)

Tag 0.6.0

## 0.5.2 (2015-08-27)

Merge pull request (#131) from ctk3b/dev-package
Fix requirements.
Contributions are welcomed via [pull requests on GitHub](https://github.com/mosdef-hub/mbuild/pulls). Developers and/or
users will review requested changes and make comments. The rest of this file will serve as a set of general guidelines
for contributors.

# Features

## Implement functionality in a general and flexible fashion

mBuild is designed to be general and flexible, not limited to single chemistries, file formats, simulation engines, or
simulation methods. Additions to core features should attempt to provide something that is applicable to a variety of
use-cases and not targeted at only the focus area of your research. However, some specific features targeted toward
a limited use case may be appropriate. Speak to the developers before writing your code and they will help you make design
choices that allow flexibility.

# Version control

We currently use the "standard" Pull Request model. Contributions should be implemented on feature branches of forks.
Please try to keep the `master` branch of your fork up-to-date with the `master` branch of the main repository.

## Propose a single set of related changes

Small changes are preferred over large changes. A major contribution can often be broken down into smaller PRs. Large PRs that
affect many parts of the codebase can be harder to review and are more likely to cause merge conflicts.

# Source code

## Use a consistent style

It is important to have a consistent style throughout the source code. The following criteria are desired:

* Lines wrapped to 80 characters
* Lines are indented with spaces
* Lines do not end with whitespace
* For other details, refer to [PEP8](https://www.python.org/dev/peps/pep-0008)

We use [pre-commit](https://pre-commit.com/) to automatically check our code style. Pre-commit is included in the dev environment and its git hooks can be installed using:

```bash
pre-commit install
```

## Document code with comments

All public-facing functions should have docstrings using the numpy style. This includes concise paragraph-style description
of what the class or function does, relevant limitations and known issues, and descriptions of arguments. Internal functions
can have simple one-liner docstrings.


# Tests

## Write unit tests

All new functionality in mBuild should be tested with automatic unit tests that execute in a few seconds. These tests
should attempt to cover all options that the user can select. All or most of the added lines of source code should be
covered by unit test(s). We currently use [pytest](https://docs.pytest.org/en/latest/), which can be executed simply by calling
`pytest` from the root directory of the package.
Developer Notes / Tools / License
=================================

Assorted notes for developers.

Most of these tools were adapted from MDTraj and are released under the following
 license:

License
-------
Copyright (c) 2012-2015 Stanford University and the Authors
All rights reserved.

Redistribution and use of all files in this folder (devtools) and
(../basesetup.py, ../setup.py) files in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
### PR Summary:

### PR Checklist
------------
 - [ ] Includes appropriate unit test(s)
 - [ ] Appropriate docstring(s) are added/updated
 - [ ] Code is (approximately) PEP8 compliant
 - [ ] Issue(s) raised/addressed?
---
name: Bug report
about: Report a bug in mBuild

---

**Bug summary**

What were you trying to do and what happened instead? Please copy and paste the stack output


**Code to reproduce the behavior**

Please include a code snippet that can be used to reproduce this bug.

```python
# Paste your code here
#
#
```

**Software versions**

- Which version of mBuild are you using? (`python -c "import mbuild as mb; print(mb.__version__)"`)
- Which version of Python (`python --version`)?
- Which operating system?
---
name: Questions/discussions
about: For usage questions, important notes, and other discussions.

---

**A summary of the question or discussion topic.**
---
name: Feature request
about: Suggest an improvement to mBuild

---

**Describe the behavior you would like added to mBuild**
A clear and concise description of what the proposed idea is.

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
mBuild
=======

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

*A hierarchical, component based molecule builder*

With just a few lines of mBuild code, you can assemble reusable components into
complex molecular systems for molecular simulations.


* mBuild is designed to minimize or even eliminate the need to explicitly translate and
  orient components when building systems: you simply tell it to connect two
  pieces!
* mBuild keeps track of the system's topology so you don't have to
  worry about manually defining bonds when constructing chemically bonded
  structures from smaller components.



mBuild is a part of the MoSDeF ecosystem
----------------------------------------

The **mBuild** software, in conjunction with the other `Molecular Simulation Design Framework (MoSDeF) <https://mosdef.org>`_ tools, supports a wide range of
simulation engines, including `Cassandra <https://cassandra.nd.edu>`_, `GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu>`_, `GROMACS <https://www.gromacs.org>`_,
`HOOMD-blue <http://glotzerlab.engin.umich.edu/hoomd-blue/>`_, and
`Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) <https://lammps.sandia.gov>`_.
The **mBuild** and **MoSDeF** tools allow simulation reproducibility
across the various simulation engines, eliminating the need to be an expert user in all
the engines to replicate, continue, or advance the existing research. Additionally,
the software can auto-generate many different systems,
allowing large-scale screening of chemicals and materials using
`Signac <https://signac.io>`_ to manage the simulations and data.

The `MoSDeF <https://mosdef.org>`_ software is comprised the following packages:
    * `mBuild <https://mbuild.mosdef.org/en/stable/>`_ -- A hierarchical, component based molecule builder
    * `foyer <https://foyer.mosdef.org/en/stable/>`_ -- A package for atom-typing as well as applying and disseminating forcefields
    * `GMSO <https://gmso.mosdef.org/en/stable/>`_ -- Flexible storage of chemical topology for molecular simulation


.. toctree::
	:caption: Getting Started
    	:maxdepth: 2

	getting_started/example_system
	getting_started/installation/installation_toc
    	getting_started/quick_start/quick_start
	getting_started/writers/writers
	getting_started/tutorials/tutorials

.. toctree::
	:caption: Topic Guides
    	:maxdepth: 2

	topic_guides/recipe_development
    	topic_guides/data_structures
	topic_guides/load_data
    	topic_guides/coordinate_transforms
    	topic_guides/recipes

.. toctree::
    	:caption: Reference
    	:maxdepth: 2

	reference/units
    	reference/citing_mbuild
    	reference/older_documentation
=============
Citing mBuild
=============

If you use mBuild for your research, please cite `our paper <http://doi.org/10.1007%2F978-981-10-1128-3_5>`_:

**ACS**

    Klein, C.; Sallai, J.; Jones, T. J.; Iacovella, C. R.; McCabe, C.; Cummings, P. T. A Hierarchical, Component Based Approach to Screening Properties of Soft Matter. In *Foundations of Molecular Modeling and Simulation. Molecular Modeling and Simulation (Applications and Perspectives)*; Snurr, R. Q., Adjiman, C. S., Kofke, D. A., Eds.; Springer, Singapore, 2016; pp 79-92.

**BibTeX**

.. code-block:: bibtex

    @Inbook{Klein2016mBuild,
        author      = "Klein, Christoph and Sallai, János and Jones, Trevor J. and Iacovella, Christopher R. and McCabe, Clare and Cummings, Peter T.",
        editor      = "Snurr, Randall Q and Adjiman, Claire S. and Kofke, David A.",
        title       = "A Hierarchical, Component Based Approach to Screening Properties of Soft Matter",
        bookTitle   = "Foundations of Molecular Modeling and Simulation: Select Papers from FOMMS 2015",
        year        = "2016",
        publisher   = "Springer Singapore",
        address     = "Singapore",
        pages       = "79--92",
        isbn        = "978-981-10-1128-3",
        doi         = "10.1007/978-981-10-1128-3_5",
        url         = "https://doi.org/10.1007/978-981-10-1128-3_5"
    }

Download as :download:`BibTeX <../files/mbuild_citation.bib>` or :download:`RIS <../files/mbuild_citation.ris>`
=====
Units
=====

mBuild automatically performs unit conversions in its reader and writer functions.
When working with an :py:class:`mbuild.Compound`, mBuild uses the following units:

+----------+----------+
| Quantity |   Units  |
+==========+==========+
| distance |    nm    |
+----------+----------+
|   angle  | radians* |
+----------+----------+

\* :py:class:`mbuild.Lattice` and :py:class:`mbuild.Box` use degrees.

See also `foyer unit documentation <https://foyer.mosdef.org/en/stable/units.html>`_ and `ele documentation <https://ele-ment.readthedocs.io/en/latest/>`_.
====================
Older Documentation
====================

Up until  `mBuild Version 0.10.4`_, the documentation is
available as pdf files. Please use the following links to download the documentation as a pdf manual:

* `Mbuild Version 0.10.4`_: `Download Here <https://github.com/mosdef-hub/mosdef-hub.github.io/raw/master/old_docs/mbuild.0.10.4.pdf>`__
* `Mbuild Version 0.10.3`_: `Download Here <https://github.com/mosdef-hub/mosdef-hub.github.io/raw/master/old_docs/mbuild.0.10.3.pdf>`__
* `Mbuild Version 0.10.1`_: `Download Here <https://github.com/mosdef-hub/mosdef-hub.github.io/raw/master/old_docs/mbuild.0.10.1.pdf>`__
* `Mbuild Version 0.9.3`_: `Download Here <https://github.com/mosdef-hub/mosdef-hub.github.io/raw/master/old_docs/mbuild.0.9.3.pdf>`__


.. _Mbuild Version 0.10.4: https://github.com/mosdef-hub/mbuild/releases/tag/0.10.4
.. _Mbuild Version 0.10.3: https://github.com/mosdef-hub/mbuild/releases/tag/0.10.3
.. _Mbuild Version 0.10.1: https://github.com/mosdef-hub/mbuild/releases/tag/0.10.1
.. _Mbuild Version 0.9.3: https://github.com/mosdef-hub/mbuild/releases/tag/0.9.3
==================
Recipe Development
==================

There may be cases where your ``Compounds`` and/or building scripts can be generalized to support a broad range of systems.
Such objects would be a valuable resource for many researchers, and might justify development of a Python package that could be distributed to the community.


``mBuild`` has been developed with this in mind, in the form of a plug-in system.
Detailed below are the specifications of this system, how to convert an existing Python project into an mBuild-discoverable plug-in, and an example.

Entry Points
------------

The basis of the plug-in system in mBuild is the `setuptools.entry_points package <https://packaging.python.org/guides/creating-and-discovering-plugins/#using-package-metadata>`_.
This allows other packages to register themselves with the ``entry_point`` group we defined in mBuild, so they are accessible through the ``mbuild.recipes`` location.
Imagine you have a class named ``my_foo`` that inherits from ``mb.Compound``.
It is currently inside of a project ``my_project`` and is accessed via a direct import, i.e. ``from my_project import my_foo``.
You can register this class as an entry point associated with ``mbuild.recipes``.
It will then be accessible from inside mBuild as a plug-in via ``mbuild.recipes.my_foo`` and a direct import will be unncessary.
The call ``import mbuild`` discovers all plug-ins that fit the ``entry_point`` group specification and makes them available under ``mbuild.recipes``.

Registering a Recipe
____________________

Here we consider the case that a user already has a Python project set up with a structure similar to the layout below.

This project can be found `here <https://github.com/justinGilmer/mbuild-fcc>`_.

::

    mbuild_fcc
    ├── LICENSE
    ├── README.md
    ├── mbuild_fcc
    │   ├── mbuild_fcc.py
    │   └── tests
    │       ├── __init__.py
    │       └── test_fcc.py
    └── setup.py


The two important files for the user to convert their ``mBuild`` plug-in to a discoverable plug-in are ``setup.py`` and ``mbuild_fcc.py``.

To begin, lets first inspect the ``mbuild_fcc.py`` file, a shortened snippet is below.

::

    import mbuild


    class FCC(mbuild.Compound):
        """Create a mBuild Compound with a repeating unit of the FCC unit cell.

        ... (shortened for viewability)

        """

        def __init__(self, lattice_spacing=None, compound_to_add=None, x=1, y=1, z=1):
            super(FCC, self).__init__()

            # ... (shortened for viewability)

    if __name__ == "__main__":
        au_fcc_lattice = FCC(lattice_spacing=0.40782,
                             compound_to_add=mbuild.Compound(name="Au"),
                             x=5, y=5, z=1)
        print(au_fcc_lattice)



There are two notable lines in this file that we need to focus on when developing this as a plug-in for mBuild.

The first is the import statement ``import mbuild``.
We must make sure that mbuild is installed since we are inheriting from ``mbuild.Compound``. When you decide to distribute your plug-in,
the dependencies must be listed.

The second is to select the name of the plug-in itself.
It is considered good practice to name it the name of your ``class``.
In this case, we will name the plug-in ``FCC``.

The last step is to edit the ``setup.py`` file such that the plug-in can be registered under the entry_point group ``mbuild.plugins``.

::

    from setuptools import setup

    setup(
        ...
        entry_points={ "mbuild.plugins":[ "FCC = mbuild_fcc.mbuild_fcc:FCC"]},
        ...
    )

The important section is the ``entry_points`` argument. Here we define the entry_point group we want to register with: ``"mbuild.plugins"``.
Finally, we tell Python what name to use when accessing this plug-in.
Earlier, we decided to call it ``FCC``.
This is denoted here by the name before the assignment operator ``FCC =``.
Next, we pass the location of the file with our plug-in: ``mbuild_fcc.mbuild_fcc`` as if we were located at the ``setup.py`` file.
Then, we provide the name of the class within that Python file we want to make discoverable ``:FCC``.

Since the ``setup.py`` file is located in the top folder of the python project, the first ``mbuild_fcc`` is the name of the folder, and the second is the name of the python file. The colon (``:``) is used when accessing the class that is in the python file itself.


Putting it all together
_______________________

Finally, we have ``FCC = mbuild_fcc.mbuild_fcc:FCC``.

To test this feature, you should clone the ``mbuild-fcc`` project listed above.

``git clone https://github.com/justinGilmer/mbuild-fcc``


Make sure you have mBuild installed, then run the command below after changing into the ``mbuild-fcc`` directory.

``cd mbuild-fcc``

``pip install -e .``

Note that this command will install this example from source in an editable format.


Trying it Out
_____________

To test that you set up your plug-in correctly, try importing mBuild:

``import mbuild``

If you do not receive error messages, your plug-in should be discoverable!

``help(mbuild.recipes.FCC)``
`
Recipes
=======

Monolayer
---------

.. autoclass:: mbuild.lib.recipes.monolayer.Monolayer
    :members:

Polymer
-------

.. autoclass:: mbuild.lib.recipes.polymer.Polymer
    :members:

Tiled Compound
--------------

.. autoclass:: mbuild.lib.recipes.tiled_compound.TiledCompound
    :members:

Silica Interface
----------------

.. autoclass:: mbuild.lib.recipes.silica_interface.SilicaInterface
    :members:

Packing
-------
.. automodule:: mbuild.packing
    :members:

Pattern
-------
.. automodule:: mbuild.pattern
    :members:
===============
Data Structures
===============

The primary building blocks in an mBuild hierarchy inherit from the
:py:class:`mbuild.Compound` class.  ``Compounds`` maintain an ordered set of ``children``
which are other ``Compounds``.  In addition, an independent, ordered dictionary
of ``labels`` is maintained through which users can reference any other
``Compound`` in the hierarchy via descriptive strings.  Every ``Compound``
knows its parent ``Compound``, one step up in the hierarchy, and knows which
``Compounds`` reference it in their ``labels``.  :py:class:`mbuild.Port` is a special type
of ``Compound`` which are used internally to connect different ``Compounds``
using the equivalence transformations described below.

``Compounds`` at the bottom of an mBuild hierarchy, the leaves of the tree, are
referred to as ``Particles`` and can be instantiated as ``foo =
mbuild.Particle(name='bar')``. Note however, that this merely serves to illustrate
that this ``Compound`` is at the bottom of the hierarchy; ``Particle`` is simply
an alias for ``Compound`` which can be used to clarify the intended role of an
object you are creating. The method :py:meth:`mbuild.Compound.particles` traverses the
hierarchy to the bottom and yields those ``Compounds``. :py:meth:`mbuild.Compound.root`
returns the compound at the top of the hierarchy.

Compound
--------

.. autoclass:: mbuild.Compound
	:members:


Box
---

.. autoclass:: mbuild.Box
	:members:

Lattice
-------

.. autoclass:: mbuild.Lattice
	:members:


Port
----

.. autoclass:: mbuild.Port
	:members:
==========================
Coordinate Transformations
==========================
The following utility functions provide mechanisms for spatial transformations for mbuild compounds:

.. autofunction:: mbuild.coordinate_transform.force_overlap

.. autofunction:: mbuild.coordinate_transform.translate

.. autofunction:: mbuild.coordinate_transform.translate_to

.. autofunction:: mbuild.coordinate_transform.rotate

.. autofunction:: mbuild.coordinate_transform.spin

.. autofunction:: mbuild.coordinate_transform.x_axis_transform

.. autofunction:: mbuild.coordinate_transform.y_axis_transform

.. autofunction:: mbuild.coordinate_transform.z_axis_transform
.. _loading_data:

===============
Loading Data
===============

Data is input into ``mBuild`` in a few different ways or from different file types.


Load
--------------------

.. automodule:: mbuild.conversion.load
    	:members:


Load CIF
--------------------

.. automodule:: mbuild.lattice.load_cif
    	:members:
Example System
===============


Components in dashed boxes are drawn by hand using, e.g., `Avogadro <https://avogadro.cc>`_ or generated elsewhere.
`mBuild <https://mbuild.mosdef.org/en/stable/>`_ builds up complex systems from simple building blocks through simple attachment sites, called a ``Port`` (i.e., connection points). Each building block is a python class that can be customized or created through the pre-built options in the ``mBuild`` library ( ``mbuild.lib`` ). A hierarchical structure of parents and children is created through these classes, which can be easily parsed or modified.
This allows `mBuild <https://mbuild.mosdef.org/en/stable/>`_ to generate chemical structures in a piecemeal fashion by creating or importing molecular sections, adding ports, and connecting the ports to form bonds.
Together with `Signac <https://signac.io>`_, this functionality enables an automatic and dynamic method for generating chemical systems, allowing large-scale chemical and materials screening with minimal user interaction.

Ultimately, complex systems can be created with just a line or two
of code. Additionally, this approach seamlessly exposes tunable parameters within
the hierarchy so you can actually create whole families of structures
by adjusting a variable or two::

    pattern = Random2DPattern(20)  # A random arrangement of 20 pieces on a 2D surface.
    brush_layer = BrushLayer(chain_lenth=20, pattern=pattern, tile_x=3, tile_y=2)

.. figure:: ../images/pmpc.png
    :width: 100 %
    :align: center

    **Zwitterionic brushes on beta-cristobalite substrate.** Example system that can be created using mBuild.
    Components in dashed boxes are created from some external tool like Avogadro or SMILES strings.
    Components in solid boxes are created from these smaller dashed components and then constructed into larger,
    more complex systems using mBuild functionality.

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.
==============
Installation
==============


.. toctree::
    	installation
    	docker
Using mBuild with Docker
========================

Docker and other containerization technologies allow entire applications
and their dependencies to be packaged and distributed as images. This
simplifies the installation process for the user and substantially
reduces platform dependence (e.g., different compiler versions, libraries,
etc). This section is a how-to guide for using mBuild with docker.

Prerequisites
-------------
A docker installation on your machine. This
`Docker installation documentation <https://docs.docker.com/get-docker/>`_ has instructions to get docker running on your machine.
If you are not familiar with docker, the Internet is full of good tutorials like these from
`Docker curriculum <https://docker-curriculum.com/>`_ and
`YouTube <https://www.youtube.com/watch?v=zJ6WbK9zFpI&feature=youtu.be>`_.

Jupyter Quick Start
-------------------
After you have a working docker installation, use the following command to
start a Jupyter notebook with mBuild and all the required dependencies:

.. code-block:: bash

    $ docker pull mosdef/mbuild:latest
    $ docker run -it --name mbuild -p 8888:8888 mosdef/mbuild:latest

If no command is provided to the container (as above), the container starts a
``jupyter-notebook`` at the (container) location ``/home/anaconda/data``.
To access the notebook, paste the notebook URL into a web browser on your
computer. When you are finished, you can use control-C to exit the notebook
as usual. The docker container will exit upon notebook shutdown.

.. warning::

    Containers by nature are ephemeral, so filesystem changes (e.g., adding
    a new notebook) only persists until the end of the container's lifecycle.
    If the container is removed, any changes or code additions will not persist.
    See the section below for persistent data.

.. note::

    The ``-it`` flags connect your keyboard to the terminal running in the
    container. You may run the prior command without those flags, but be
    aware that the container will not respond to any keyboard input. In
    that case, you would need to use the ``docker ps`` and ``docker kill``
    commands to shut down the container.


Persisting User Volumes
-----------------------
If you are using mBuild from a docker container and need access to data
on your local machine or you wish to save files generated in the container,
you can mount user volumes in the container. User volumes will provide a way
to persist filesystem changes made to a container regardless of the container
lifecycle. For example, you might want to create a directory called
``mbuild-notebooks`` on your local system, which will store all of your mBuild
notebooks/code. In order to make that accessible from within the container
(where the notebooks will be created/edited), use the following steps:

.. code-block:: bash

    $ mkdir mbuild-notebooks
    $ cd mbuild-notebooks/
    $ docker run -it --name mbuild --mount type=bind,source=$(pwd),target=/home/anaconda/data -p 8888:8888  mosdef/mbuild:latest

You can easily mount a different directory from your local machine by changing
``source=$(pwd)`` to ``source=/path/to/my/favorite/directory``.

.. note::

    The ``--mount`` flag mounts a volume into the docker container. Here we
    use a ``bind`` mount to bind the current directory on our local filesystem
    to the ``/home/anaconda/data`` location in the container. The files you see
    in the ``jupyter-notebook`` browser window are those that exist on your
    local machine.

.. warning::

    If you are using the container with jupyter notebooks you should use
    the ``/home/anaconda/data`` location as the mount point inside the container;
    this is the default notebook directory.

Running Python scripts in the container
---------------------------------------
Jupyter notebooks are a great way to explore new software and prototype
code. However, when it comes time for production science, it is often
better to work with python scripts. In order to execute a python script
(``example.py``) that exists in the current working directory of your
local machine, run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/mbuild:latest "python data/test.py"

Note that once again we are ``bind`` mounting the current working directory
to ``/home/anaconda/data``. The command we pass to the container is
``python data/test.py``. Note the prefix ``data/`` to the script; this is because
we enter the container in the home folder (``/home/anaconda``), but our script
is located under ``/home/anaconda/data``.

.. warning::

    Do not bind mount to ``target=/home/anaconda``. This will cause errors.


If you don't require a Jupyter notebook, but just want a Python interpreter,
you can run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/mbuild:latest python

If you don't need access to any local data, you can of course drop the
``--mount`` command:

.. code-block:: bash

    $ docker run mosdef/mbuild:latest python


Different mBuild versions
-------------------------
Instead of using ``latest``, you can use the image ``mosdef/mbuild:stable``
for most recent stable release of mBuild.

Cleaning Up
-----------
You can remove the *container* by using the following command.

.. code-block:: bash

    $ docker container rm mbuild

The *image* will still exist on your machine. See the tutorials at the
top of this page for more information.

.. warning::

    You will not be able to start a second container with the same name
    (e.g., mbuild), until the first container has been removed.

.. note::

    You do not need to name the container `mbuild` as shown in the above
    examples (``--name mbuild``). Docker will give each container a name
    automatically. To see all the containers on your machine, run
    ``docker ps -a``.
============
Installation
============

Install with `conda <https://repo.anaconda.com/miniconda/>`_
-----------------------------------------------------------------
::

    $ conda install -c conda-forge mbuild

Alternatively you can add all the required channels to your ``.condarc``
after which you can simply install without specifying the channels::

    $ conda config --add channels conda-forge
    $ conda install mbuild

.. note::
    The order in which channels are added matters: ``conda-forge`` should be the highest priority as a result of being added last. In your ``.condarc`` file, it should be listed first.

.. note::
    Because ``packmol`` binaries are unavailable for windows from ``conda-forge`` channel, to use mbuild with conda in a Windows system requires the ``omnia`` channel. Use the following command to use ``mbuild`` with conda in a Windows system::

        $ conda install -c conda-forge -c omnia mbuild

.. note::
    The `MDTraj website <http://mdtraj.org/1.9.3/new_to_python.html>`_ makes a
    nice case for using Python and in particular the
    `Anaconda scientific python distribution <https://www.anaconda.com/products/individual>`_
    to manage your numerical and scientific Python packages.

Install an editable version from source
---------------------------------------

To make your life easier, we recommend that you use a pre-packaged Python
distribution like `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
in order to get all of the dependencies::

    $ git clone https://github.com/mosdef-hub/mbuild
    $ cd mbuild
    $ conda env create -f environment-dev.yml
    $ conda activate mbuild-dev
    $ pip install -e .

.. note::
    The above installation is for OSX and Unix. If you are using Windows, use environment-win.yml instead of environment-dev.yml


Install pre-commit
------------------

We use `pre-commit <https://pre-commit.com/>`_ to automatically handle our code formatting and this package is included in the dev environment.
With the ``mbuild-dev`` conda environment active, pre-commit can be installed locally as a git hook by running::

     $ pre-commit install

And (optional) all files can be checked by running::

     $ pre-commit run --all-files


Supported Python Versions
-------------------------

Python 3.6, 3.7 and 3.8 are officially supported, including testing during
development and packaging. Support for Python 2.7 has been dropped as of
August 6, 2019. Other Python versions, such as 3.9 and 3.5 and older, may
successfully build and function but no guarantee is made.

Testing your installation
-------------------------

mBuild uses `py.test <https://docs.pytest.org/en/stable/>`_ for unit testing. To run them simply run the following while in the base directory::

    $ conda install pytest
    $ py.test -v

Building the documentation
--------------------------

mBuild uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::

    $ pip install -r requirements.txt
    $ make html
----------------------------------------------
File Writers
----------------------------------------------

The mBuild library also supports simulation engine-specific file writers.  These writers create a complete set of simulation writers to input files or a partial set of file writers, where the other required files are generated via another means.

mBuild utilizes ParmEd to write ``Compound`` information to a variety of file
formats (e.g. PDB, MOL2, GRO.  The full list of formats supported by ParmEd
can be found at the `ParmEd website <http://parmed.github.io/ParmEd/html/readwrite.html>`_).
Additionally, mBuild features several internal writers for file formats not yet
supported by ParmEd. Information on these internal writers can be found below.

By default, many mBuild functions will only write coordinate and bond information to these files,
i.e. no angles or dihedrals, and no atom typing is performed (atom names are used
as atom types). However, force fields can be applied to Compounds by passing force
field XML files (used by the `Foyer package <https://github.com/mosdef-hub/foyer>`_)
to the ``Compound.save`` function if Foyer is installed. If a force field is applied to a
Compound, the mBuild internal writers will also write angle and dihedral information
to the file in addition to labelling atoms by the atom types specified by the force
field.  The CHARMM-style GOMC writers are the exception to this default rule since
they need a force field to build the files, as these files depend on the force field parameters (Example: charge and MW in the PSF files).

The simulation engine writers that use mBuild or are currently contained in the mBuild library:


* `Cassandra <https://cassandra.nd.edu/>`_
* `GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu/>`_
* `GROMACS <https://www.gromacs.org/>`_
* `HOOMD-blue <http://glotzerlab.engin.umich.edu/hoomd-blue//>`_
* `Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) <https://lammps.sandia.gov/>`_


.. toctree::

	cassandra_file_writers
   	GOMC_file_writers
   	HOOMD_blue_file_writers
	LAMMPS_file_writers
GOMC and NAMD File Writers
======================================================


CHARMM-style PDB, PSF, and Force Field File Writers
--------------------------------------------------------

	.. autoclass:: mbuild.formats.charmm_writer.Charmm
		:special-members: __init__
		:members:




GOMC Control File Writer
--------------------------------------------------------

	.. automodule:: mbuild.formats.gomc_conf_writer
    		:members: write_gomc_control_file



NAMD Control File Writer
--------------------------------------------------------

The NAMD control file writer is not currently available.
LAMMPS File Writers
===================

Write Lammps data
-----------------

	.. automodule:: mbuild.formats.lammpsdata
    		:members:
HOOMD-blue File Writers
===========================


Write GSD (General Simulation Data)
-----------------------------------
Default data file format for HOOMD-blue
+++++++++++++++++++++++++++++++++++++++

    .. automodule:: mbuild.formats.gsdwriter
        :members:

Create HOOMD-blue force field (>= 3.0)
--------------------------------------------------------

    .. automodule:: mbuild.formats.hoomd_forcefield
            :members:

Create HOOMD-blue Simulation (v2.x)
--------------------------------------------------------

    .. automodule:: mbuild.formats.hoomd_simulation
            :members:


HOOMD-blue Snapshot
--------------------------------------------------------

    .. automodule:: mbuild.formats.hoomd_snapshot
            :members:
Cassandra File Writers
===========================

	.. automodule:: mbuild.formats.cassandramcf
    		:members:
.. _QuickStart_Load_files:

Load files
========================


mol2 files
------------------------

Create an ``mbuild.Compound`` (i.e., the "pentane" variable) by loading a molecule from a `mol2 <http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf>`_ file.

Import the required mbuild packages.

.. code:: ipython3

    import mbuild as mb


Load the "pentane.mol2" file from its directory.

.. code:: ipython3

    pentane = mb.load("path_to_mol2_file/pentane.mol2")


CIF files
------------------------

Build an ``mbuild.Compound`` (i.e., the "ETV_triclinic" variable) by loading a `Crystallographic Information File (CIF) <https://www.iucr.org/resources/cif>`_ file and selecting the number of cell units to populate in the x, y, and z-dimensions.


Import the required mbuild packages.

.. code:: ipython3

    import mbuild as mb
    from mbuild.lattice import load_cif



The `CIF <https://www.iucr.org/resources/cif>`_ file is loaded using the ``load_cif`` function. Next, three (3) cell units shall be built for all the x, y, and z-dimensions with the populate function.  Finally, the `CIF <https://www.iucr.org/resources/cif>`_'s residues are named 'ETV'.

.. code:: ipython3

    lattice_cif_ETV_triclinic = load_cif("path_to_cif_file/ETV_triclinic.cif")
    ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=3, y=3, z=3)
    ETV_triclinic.name = 'ETV'


Other file types
------------------------
mBuild also supports :ref:`loading_data` or files via hoomd_snapshot, GSD, SMILES strings, and ParmEd structures.
-----------------------
Quick Start
-----------------------

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: http://opensource.org/licenses/MIT

The `MoSDeF <https://mosdef.org>`_ software is comprised the following packages:
    * `mBuild <https://mbuild.mosdef.org/en/stable/>`_ -- A hierarchical, component based molecule builder
    * `foyer <https://foyer.mosdef.org/en/stable/>`_ -- A package for atom-typing as well as applying and disseminating forcefields
    * `GMSO <https://gmso.mosdef.org/en/stable/>`_ -- Flexible storage of chemical topology for molecular simulation

.. note::
    **foyer** and **GMSO** are used together with **mBuild** to create all the required files to conduct the simulations. Run time parameters for a simulation engine need to be created by the user.

In the following examples, different types of simulation boxes are constructed using the **MoSDeF** software.


Molecular simulations are usually comprised of many molecules contained in a
box (NPT and NVT ensembles), or boxes (GEMC and GCMC ensembles).
The **mBuild** library allows for easy generation of the simulation
box or boxes utilizing only a few lines of python code.


The following tutorials are available either as html or interactive `jupyter <https://jupyter.org/>`_ notebooks.


.. toctree::

   	load_files
   	Box_example
	fill_box_example
   	polymer_example
Fill Box
========


All-Atom (AA) Hexane and Ethanol System
---------------------------------------

.. note::
    `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with ``mBuild`` in
    the following example to demonstrate how the `MoSDeF <https://mosdef.org>`_
    libraries can be used together to generate a simulation box.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct an all-atom (AA) hexane and ethanol using the OPLS-AA force field (FF),
which is shipped as a standard `foyer <https://foyer.mosdef.org/en/stable/>`_ FF.
The hexane and ethanol molecules will be created using `smiles strings <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_.
The hexane and ethanol residues will be named `"HEX"` and `"ETO"`, respectively.
Lastly, the hexane and ethanol molecule's configuration will be energy minimized, properly reorienting the molecule to the specified FF, which is sometimes needed for some simulation engines to ensure the initial configuration energy is not too high.

.. note::
    The energy minimize step requires the `foyer <https://foyer.mosdef.org/en/stable/>`_ package.

.. code:: ipython3

    hexane = mb.load('CCCCCC', smiles=True)
    hexane.name = 'HEX'
    hexane.energy_minimize(forcefield='oplsaa', steps=10**4)


    ethanol = mb.load('CCO', smiles=True)
    ethanol.name = 'ETO'
    ethanol.energy_minimize(forcefield='oplsaa', steps=10**4)

The liquid box is built to a density of 680 kg/m\ :sup:`3`, with a 50/50 mol ratio of hexane and ethanol, and will be in an orthogonal box measuring 5.0 nm in the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound= [hexane, ethanol],
                          density=680,
                          compound_ratio=[0.5, 0.5],
                          box=[5.0, 5.0, 5.0])


United Atom (UA) Methane System
-------------------------------

.. note::
    `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with ``mBuild`` in
    the following example to demonstrate how the `MoSDeF <https://mosdef.org>`_ libraries
    integrate to generate a simulation box.  A subset of the `TraPPE-United Atom <http://trappe.oit.umn.edu>`_
    force field (FF) comes standard with the `foyer <https://foyer.mosdef.org/en/stable/>`_ software package.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct a pseudo-monatomic molecule (united atom (UA) methane), for use with the
`TraPPE <http://trappe.oit.umn.edu>`_ FF.  The UA methane, bead type `"_CH4"`, will be built as a child (``mbuild.Compound.children``), so the parent (``mbuild.Compound``) will
allow a user-selected residue name (``mbuild.Compound.name``). If the methane is built using ``methane = mb.Compound(name="_CH4")``, then the user must keep the residue name `"_CH4"` or `foyer <https://foyer.mosdef.org/en/stable/>`_ will not recognize the bead type when using the standard TraPPE force field XML file.

.. code:: ipython3

    methane = mb.Compound(name="MET")
    methane_child_bead = mb.Compound(name="_CH4")
    methane.add(methane_child_bead, inherit_periodicity=False)

.. note::
    The ``inherit_periodicity`` flag is an optional boolean (default=True), which replaces
    the periodicity of self with the periodicity of the Compound being added.

The orthogonal liquid box contains 1230 methane molecules and measures 4.5 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound=methane,
                          n_compounds=1230,
                          box=[4.5, 4.5, 4.5]
                          )
Box
========================


Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Orthogonal Box
------------------------

Build an empty orthogonal **mBuild** Box (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 90) measuring 4.0 nm in all the x, y, and z-dimensions.

.. note::
    Note: if the angles are not specified, the system will default to an orthogonal box
    (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 90).

.. code:: ipython3

    empty_box = mb.Box(lengths=[4.0, 4.0, 4.0], angles=[90, 90, 90])


Non-Orthogonal Box
------------------------

Build an empty non-orthogonal **mBuild** Box (i.e., the angle in degrees are 𝛼 = 90, 𝛽 = 90, 𝛾 = 120) measuring 4.0 nm in the x and y-dimensions, and 5.0 nm in the z-dimension.

.. code:: ipython3

    empty_box = mb.Box(lengths=[4.0, 4.0, 5.0], angles=[90, 90, 120])
Polymer
========================

Use two (2) different monomer units, A and B, to construct a polymer, capping it with a carboxylic acid and amine end group.


Import the required mbuild packages.

.. code:: ipython3

    import mbuild as mb
    from mbuild.lib.recipes.polymer import Polymer



Create the monomer units `comp_1` and `comp_2` using `SMILES strings <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_.
Set the `chain` as a ``Polymer`` class, adding `comp_1` and `comp_2` as the monomers A and B to the polymer.

.. note::
    Setting the indices identifies which atoms will be removed and have ports created in their place.


.. code:: ipython3

    comp_1 = mb.load('CC', smiles=True) # mBuild compound of the monomer unit
    comp_2 = mb.load('COC', smiles=True) # mBuild compound of the monomer unit
    chain = Polymer()
    chain.add_monomer(compound=comp_1,
                      indices=[2, -2],
                      separation=.15,
                      replace=True)

    chain.add_monomer(compound=comp_2,
                      indices=[3, -1],
                      separation=.15,
                      replace=True)


Select the carboxylic acid and amine end groups that we want to use for the head and tail of the polymer.
Then, build the polymer with three (3) iterations of the AB sequence, and the selected head and tail end groups.


.. code:: ipython3

    chain.add_end_groups(mb.load('C(=O)O',smiles=True),
                         index=3,
                         separation=0.15,
                         duplicate=False,
		         label="head")

    chain.add_end_groups(mb.load('N', smiles=True),
                         index=-1,
		         separation=0.13,
                         duplicate=False,
		         label="tail")

    chain.build(n=3, sequence='AB')
    chain.visualize()


.. figure:: ../../images/polymer_example_image.png
    :width: 60 %
    :align: center

    This **example polymer** is 3 of the AB sequences together with carboxylic acid and amine end groups.
Building a Simple Alkane
========================

The purpose of this tutorial is to demonstrate the construction of an
alkane polymer and provide familiarity with many of the underlying
functions in mBuild. Note that a robust polymer construction recipe
already exists in mBuild, which will also be demonstrated at the end of
the tutorial.

Setting up the monomer
----------------------

The first step is to construct the basic repeat unit for the alkane,
i.e., a :math:`CH_2` group, similar to the construction of the
:math:`CH_3` monomer in the prior methane tutorial. Rather than
importing the coordinates from a pdb file, as in the previous example,
we will instead explicitly define them in the class. Recall that
distance units are nm in mBuild.

.. code:: ipython3

    import mbuild as mb

    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
            # Add carbon
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')

            # Add hydrogens
            self.add(mb.Particle(name='H', pos=[-0.109, 0, 0.0]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0.109, 0, 0.0]), label='HC[$]')

            # Add bonds between the atoms
            self.add_bond((self['C'][0], self['HC'][0]))
            self.add_bond((self['C'][0], self['HC'][1]))

            # Add ports anchored to the carbon
            self.add(mb.Port(anchor=self[0]), label='up')
            self.add(mb.Port(anchor=self[0]), label='down')

            # Move the ports approximately half a C-C bond length away from the carbon
            self['up'].translate([0, -0.154/2, 0])
            self['down'].translate([0, 0.154/2, 0])

    monomer = CH2()
    monomer.visualize(show_ports=True)

This configuration of the monomer is not a particularly realistic
conformation. One could use this monomer to construct a polymer and then
apply an energy minimization scheme, or, as we will demonstrate here, we
can use mBuild’s rotation commands to provide a more realistic starting
point.

Below, we use the same basic script, but now apply a rotation to the
hydrogen atoms. Since the hydrogens start 180° apart and we know they
should be ~109.5° apart, each should be rotated half of the difference
closer to each other around the y-axis. Note that the rotation angle is
given in radians. Similarly, the ports should be rotated around the
x-axis by the same amount so that atoms can be added in a realistic
orientation.

.. code:: ipython3

    import numpy as np
    import mbuild as mb

    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
            # Add carbon
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')

            # Add hydrogens
            self.add(mb.Particle(name='H', pos=[-0.109, 0, 0.0]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0.109, 0, 0.0]), label='HC[$]')

            # Rotate the hydrogens
            theta = 0.5 * (180 - 109.5) * np.pi / 180
            #mb.rotate(self['HC'][0], theta, around=[0, 1, 0])
            #mb.rotate(self['HC'][1], -theta, around=[0, 1, 0])
            self['HC'][0].rotate( theta, around=[0, 1, 0])
            self['HC'][1].rotate(-theta, around=[0, 1, 0])

            # Add bonds between the atoms
            self.add_bond((self['C'][0], self['HC'][0]))
            self.add_bond((self['C'][0], self['HC'][1]))

            # Add the ports and appropriately rotate them
            self.add(mb.Port(anchor=self[0]), label='up')
            self['up'].translate([0, -0.154/2, 0])
            self['up'].rotate(theta, around=[1, 0, 0])

            self.add(mb.Port(anchor=self[0]), label='down')
            self['down'].translate([0, 0.154/2, 0])
            self['down'].rotate(-theta, around=[1, 0, 0])

    monomer = CH2()
    monomer.visualize(show_ports=True)

Defining the polymerization class
---------------------------------

With a basic monomer construct, we can now construct a polymer by
connecting the ports together. Here, we first instantiate one instance
of the CH2 class as ``1ast_monomer``, then use the clone function to
make a copy. The ``force_overlap()`` function is used to connect the
``'up'`` port from ``current_monomer`` to the ``'down'`` port of
``last_mononer``.

.. code:: ipython3

    class AlkanePolymer(mb.Compound):
        def __init__(self):
            super(AlkanePolymer, self).__init__()
            last_monomer = CH2()
            self.add(last_monomer)
            for i in range(3):
                current_monomer = CH2()
                mb.force_overlap(move_this=current_monomer,
                                 from_positions=current_monomer['up'],
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer = current_monomer

    polymer = AlkanePolymer()
    polymer.visualize(show_ports=True)

Visualization of this structure demonstrates a problem; the polymer
curls up on itself. This is a result of the fact that ports not only
define the location in space, but also an orientation. This can be
trivially fixed, by rotating the down port 180° around the y-axis.

We can also add a variable ``chain_length`` both to the for loop and
``init`` that will allow the length of the polymer to be adjusted when
the class is instantiated.

.. code:: ipython3

    import numpy as np
    import mbuild as mb

    class CH2(mb.Compound):
        def __init__(self):
            super(CH2, self).__init__()
             # Add carbons and hydrogens
            self.add(mb.Particle(name='C', pos=[0,0,0]), label='C[$]')
            self.add(mb.Particle(name='H', pos=[-0.109, 0, 0.0]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0.109, 0, 0.0]), label='HC[$]')

            # rotate hydrogens
            theta = 0.5 * (180 - 109.5) * np.pi / 180
            self['HC'][0].rotate(theta, around=[0, 1, 0])
            self['HC'][1].rotate(-theta, around=[0, 1, 0])

            # Add bonds between the atoms
            self.add_bond((self['C'][0], self['HC'][0]))
            self.add_bond((self['C'][0], self['HC'][1]))

            # Add ports
            self.add(mb.Port(anchor=self[0]), label='up')
            self['up'].translate([0, -0.154/2, 0])
            self['up'].rotate(theta, around=[1, 0, 0])

            self.add(mb.Port(anchor=self[0]), label='down')
            self['down'].translate([0, 0.154/2, 0])
            self['down'].rotate(np.pi, [0, 1, 0])
            self['down'].rotate(-theta, around=[1, 0, 0])


    class AlkanePolymer(mb.Compound):
        def __init__(self, chain_length=1):
            super(AlkanePolymer, self).__init__()
            last_monomer = CH2()
            self.add(last_monomer)
            for i in range (chain_length-1):
                current_monomer = CH2()

                mb.force_overlap(move_this=current_monomer,
                                 from_positions=current_monomer['up'],
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer=current_monomer

.. code:: ipython3

    polymer = AlkanePolymer(chain_length=10)
    polymer.visualize(show_ports=True)

Using mBuild’s Polymer Class
----------------------------

``mBuild`` provides a prebuilt class to perform this basic
functionality. Since it is designed to be more general, it takes as an
argument not just the replicates (``n``), ``sequence`` ('A' for a single monomer or 'AB' for two different monomers).
Then, it binds them together by removing atom/bead via specifying its index number (``indices``).
A graphical description of the polymer builder creating ports, then bonding them together is provided below.

.. figure:: ../../images/polymer_image.png
    :width: 100 %
    :align: center

    **Polymer builder class example.** This shows how to define the atoms, which are replaced with ports.  The ports are then bonded together between the monomers.  Additionally, these ports can be utilized for adding different end groups moieties to the polymer.

.. note::
    The port locations may be critical to ensure the molecule is not overlapping when it is built.

Building a Simple Hexane
----------------------------
A simple hexane molecule is built using ``mBuild``'s packaged polymer builder.
This is done by loading a methane molecule via a SMILES string.
The indices are explicitly selected, so the molecule builds out in the proper directions and does not overlap.

.. code:: ipython3

    import mbuild as mb
    from mbuild.lib.recipes.polymer import Polymer

    comp = mb.load('C', smiles=True) # mBuild compound of the monomer unit
    chain = Polymer()

    chain.add_monomer(compound=comp,
                      indices=[1, -2],
                      separation=.15,
                      replace=True)

    chain.build(n=6, sequence='A')


Using Multiple Monomers and Capping the Ends of a Polymer
---------------------------------------------------------
This example uses methyl ether and methane monomers to build a polymer, capping it with fluorinated and alcohol end groups.
The monomers are combined together in the 'AB' sequence two times (n=2), which means the polymer will contain 2 of each monomer (ABAB).
The end groups are added via the ``add_end_groups`` attribute, specifying the atom to use (``index``), the distance of the bond (``separation``),
the location of each end group (``label``), and if the tail end group is duplicated to the head of the polymer (``duplicate``).
The indices are explicitly selected, so the molecule builds out in the proper directions and does not overlap.

.. code:: ipython3

    from mbuild.lib.recipes.polymer import Polymer
    import mbuild as mb

    comp_1 = mb.load('C', smiles=True)
    comp_2 = mb.load('COC', smiles=True)
    chain = Polymer()

    chain.add_monomer(compound=comp_1,
                      indices=[1, -1],
                      separation=.15,
                      replace=True)

    chain.add_monomer(compound=comp_2,
                      indices=[3, -1],
                      separation=.15,
                      replace=True)


    chain.add_end_groups(mb.load('O',smiles=True), # Capping off this polymer with an Alcohol
                         index=1,
                         separation=0.15, label="head", duplicate=False)

    chain.add_end_groups(mb.load('F',smiles=True), # Capping off this polymer with a Fluorine
                         index=1,
                         separation=0.18, label="tail", duplicate=False)


    chain.build(n=2, sequence='AB')
    chain.visualize(show_ports=True)

Building a System of Alkanes
----------------------------

A system of alkanes can be constructed by simply cloning the polymer
constructed above and translating and/or rotating the alkanes in space.
``mBuild`` provides many routines that can be used to create different
patterns, to which the polymers can be shifted.

.. code:: ipython3

    comp = mb.load('C', smiles=True) # mBuild compound of the monomer unit
    polymer = Polymer()

    polymer.add_monomer(compound=comp,
                        indices=[1, -2],
                        separation=.15,
                        replace=True)

    polymer.build(n=10, sequence='A')

    # the pattern we generate puts points in the xy-plane, so we'll rotate the polymer
    # so that it is oriented normal to the xy-plane
    polymer.rotate(np.pi/2, [1, 0, 0])

    # define a compound to hold all the polymers
    system = mb.Compound()

    # create a pattern of points to fill a disk
    # patterns are generated between 0 and 1,
    # and thus need to be scaled to provide appropriate spacing
    pattern_disk = mb.DiskPattern(50)
    pattern_disk.scale(5)

    # now clone the polymer and move it to the points in the pattern
    for pos in pattern_disk:
        current_polymer = mb.clone(polymer)
        current_polymer.translate(pos)
        system.add(current_polymer)

    system.visualize()

Other patterns can be used, e.g., the ``Grid3DPattern``. We can also use
the rotation commands to randomize the orientation.

.. code:: ipython3

    import random

    comp = mb.load('C', smiles=True)
    polymer = Polymer()

    polymer.add_monomer(compound=comp,
                        indices=[1, -2],
                        separation=.15,
                        replace=True)

    polymer.build(n=10, sequence='A')

    system = mb.Compound()
    polymer.rotate(np.pi/2, [1, 0, 0])

    pattern_disk = mb.Grid3DPattern(5, 5, 5)
    pattern_disk.scale(8.0)

    for pos in pattern_disk:
        current_polymer = mb.clone(polymer)
        for around in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]:  # rotate around x, y, and z
            current_polymer.rotate(random.uniform(0, np.pi), around)
        current_polymer.translate(pos)
        system.add(current_polymer)

    system.visualize()


``mBuild`` also provides an interface to ``PACKMOL``, allowing the
creation of a randomized configuration.

.. code:: ipython3

    comp = mb.load('C', smiles=True) # mBuild compound of the monomer unit
    polymer = Polymer()

    polymer.add_monomer(compound=comp,
                        indices=[1, -2],
                        separation=.15,
                        replace=True)

    polymer.build(n=5, sequence='A')

    system = mb.fill_box(polymer, n_compounds=100, overlap=1.5, box=[10,10,10])
    system.visualize()

Variations
----------

Rather than a linear chain, the ``Polymer`` class we wrote can be easily
changed such that small perturbations are given to each port. To avoid
accumulation of deviations from the equilibrium angle, we will clone an
unperturbed monomer each time (i.e., ``monomer_proto``) before applying
a random variation.

We also define a variable ``delta``, which will control the maximum
amount of perturbation. Note that large values of ``delta`` may result
in the chain overlapping itself, as ``mBuild`` does not currently
include routines to exclude such overlaps.

.. code:: ipython3

    import mbuild as mb

    import random

    class AlkanePolymer(mb.Compound):
        def __init__(self, chain_length=1, delta=0):
            super(AlkanePolymer, self).__init__()
            monomer_proto = CH2()
            last_monomer = CH2()
            last_monomer['down'].rotate(random.uniform(-delta,delta), [1, 0, 0])
            last_monomer['down'].rotate(random.uniform(-delta,delta), [0, 1, 0])
            self.add(last_monomer)
            for i in range(chain_length-1):
                current_monomer = mb.clone(monomer_proto)
                current_monomer['down'].rotate(random.uniform(-delta,delta), [1, 0, 0])
                current_monomer['down'].rotate(random.uniform(-delta,delta), [0, 1, 0])
                mb.force_overlap(move_this=current_monomer,
                                 from_positions=current_monomer['up'],
                                 to_positions=last_monomer['down'])
                self.add(current_monomer)
                last_monomer=current_monomer

    polymer = AlkanePolymer(chain_length = 200, delta=0.4)
    polymer.visualize()
Ethane: Reading from files, Ports and coordinate transforms
-----------------------------------------------------------

**Note**: mBuild expects all distance units to be in nanometers.

In this example, we’ll cover reading molecular components from files,
introduce the concept of ``Ports`` and start using some coordinate
transforms.

First, we need to import the mbuild package:

.. code:: ipython3

    import mbuild as mb

As you probably noticed while creating your methane molecule in the last
tutorial, manually adding ``Particles`` and ``Bonds`` to a ``Compound``
is a bit cumbersome. The easiest way to create small, reusable
components, such as methyls, amines or monomers, is to hand draw them
using software like `Avogadro <https://avogadro.cc/>`__ and
export them as either a .pdb or .mol2 file (the file should contain
connectivity information).

Let’s start by reading a methyl group from a ``.pdb`` file:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.visualize()


Now let’s use our first coordinate transform to center the methyl at its
carbon atom:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.translate(-ch3[0].pos)  # Move carbon to origin.

Now we have a methyl group loaded up and centered. In order to connect
``Compounds`` in mBuild, we make use of a special type of ``Compound``:
the ``Port``. A ``Port`` is a ``Compound`` with two sets of four “ghost”
``Particles`` that assist in bond creation. In addition, ``Ports`` have an ``anchor`` attribute which
typically points to a particle that the ``Port`` should be associated
with. In our methyl group, the ``Port`` should be anchored to the carbon
atom so that we can now form bonds to this carbon:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.translate(-ch3[0].pos)  # Move carbon to origin.

    port = mb.Port(anchor=ch3[0])
    ch3.add(port, label='up')

    # Place the port at approximately half a C-C bond length.
    ch3['up'].translate([0, -0.07, 0])

By default, ``Ports`` are never output from the mBuild structure.
However, it can be useful to look at a molecule with the ``Ports`` to
check your work as you go:

.. code:: ipython3

    ch3.visualize(show_ports=True)

Now we wrap the methyl group into a python class, so that we can reuse
it as a component to build more complex molecules later.

.. code:: ipython3

    import mbuild as mb

    class CH3(mb.Compound):
        def __init__(self):
            super(CH3, self).__init__()

            mb.load('ch3.pdb', compound=self)
            self.translate(-self[0].pos)  # Move carbon to origin.

            port = mb.Port(anchor=self[0])
            self.add(port, label='up')
            # Place the port at approximately half a C-C bond length.
            self['up'].translate([0, -0.07, 0])

When two ``Ports`` are connected, they are forced to overlap in space
and their parent ``Compounds`` are rotated and translated by the same
amount.

**Note:** If we tried to connect two of our methyls right now using only
one set of four ghost particles, not only would the ``Ports`` overlap
perfectly, but the carbons and hydrogens would also perfectly overlap -
the 4 ghost atoms in the ``Port`` are arranged identically with respect
to the other atoms. For example, if a ``Port`` and its direction is
indicated by “<-”, forcing the port in <-CH3 to overlap with <-CH3 would
just look like <-CH3 (perfectly overlapping atoms).

To solve this problem, every port contains a second set of 4 ghost atoms
pointing in the opposite direction. When two ``Compounds`` are
connected, the port that places the anchor atoms the farthest away from
each other is chosen automatically to prevent this overlap scenario.

When <->CH3 and <->CH3 are forced to overlap, the CH3<->CH3 is
automatically chosen.

Now the fun part: stick ’em together to create an ethane:

.. code:: ipython3

    ethane = mb.Compound()

    ethane.add(CH3(), label="methyl_1")
    ethane.add(CH3(), label="methyl_2")
    mb.force_overlap(move_this=ethane['methyl_1'],
                             from_positions=ethane['methyl_1']['up'],
                             to_positions=ethane['methyl_2']['up'])

Above, the ``force_overlap()`` function takes a ``Compound`` and then
rotates and translates it such that two other ``Compounds`` overlap.
Typically, as in this case, those two other ``Compounds`` are ``Ports``
- in our case, ``methyl1['up']`` and ``methyl2['up']``.

.. code:: ipython3

    ethane.visualize()

.. code:: ipython3

    ethane.visualize(show_ports=True)

Similarly, if we want to make ethane a reusable component, we need to
wrap it into a python class.

.. code:: ipython3

    import mbuild as mb

    class Ethane(mb.Compound):
        def __init__(self):
            super(Ethane, self).__init__()

            self.add(CH3(), label="methyl_1")
            self.add(CH3(), label="methyl_2")
            mb.force_overlap(move_this=self['methyl_1'],
                             from_positions=self['methyl_1']['up'],
                             to_positions=self['methyl_2']['up'])

.. code:: ipython3

    ethane = Ethane()
    ethane.visualize()

.. code:: ipython3

    # Save to .mol2
    ethane.save('ethane.mol2', overwrite=True)
Methane: Compounds and bonds
----------------------------

**Note**: mBuild expects all distance units to be in nanometers.

The primary building block in mBuild is a ``Compound``. Anything you
construct will inherit from this class. Let’s start with some basic
imports and initialization:

.. code:: ipython3

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()

Any ``Compound`` can contain other ``Compounds`` which can be added
using its ``add()`` method. ``Compounds`` at the bottom of such a
hierarchy are referred to as ``Particles``. Note however, that this is
purely semantic in mBuild to help clearly designate the bottom of a
hierarchy.

.. code:: ipython3

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')

By default a created ``Compound/Particle`` will be placed at ``0, 0, 0``
as indicated by its ``pos`` attribute. The ``Particle`` objects
contained in a ``Compound``, the bottoms of the hierarchy, can be
referenced via the ``particles`` method which returns a generator of all
``Particle`` objects contained below the ``Compound`` in the hierarchy.

**Note:** All positions in mBuild are stored in nanometers.

Any part added to a ``Compound`` can be given an optional, descriptive
string label. If the label ends with the characters ``[$]``, a list will
be created in the labels. Any subsequent parts added to the ``Compound``
with the same label prefix will be appended to the list. In the example
above, we’ve labeled the hydrogen as ``HC[$]``. So this first part, with
the label prefix ``HC``, is now referenceable via ``self['HC'][0]``. The
next part added with the label ``HC[$]`` will be referenceable via
``self['HC'][1]``.

Now let’s use these styles of referencing to connect the carbon to the
hydrogen. Note that for typical use cases, you will almost never have to
explicitly define a bond when using mBuild - this is just to show you
what’s going on under the hood:

.. code:: ipython3

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')

            self.add_bond((self[0], self['HC'][0]))

As you can see, the carbon is placed in the zero index of ``self``. The
hydrogen could be referenced via ``self[1]`` but since we gave it a
fancy label, it’s also referenceable via ``self['HC'][0]``.

Alright now that we’ve got the basics, let’s finish building our
``Methane`` and take a look at it:

.. code:: ipython3

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.1, 0, -0.07])
            self.add(hydrogen, label='HC[$]')

            self.add_bond((self[0], self['HC'][0]))

            self.add(mb.Particle(name='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0, 0.1, 0.07]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0, -0.1, 0.07]), label='HC[$]')

            self.add_bond((self[0], self['HC'][1]))
            self.add_bond((self[0], self['HC'][2]))
            self.add_bond((self[0], self['HC'][3]))

.. code:: ipython3

    methane = Methane()
    methane.visualize()

.. code:: ipython3

    # Save to .mol2
    methane.save('methane.mol2',overwrite=True)
``Lattice`` Tutorial
====================

The following notebook provides a thorough walkthrough to using the
``Lattice`` class to build up crystal systems.

``Lattice`` Functionality
-------------------------

-  **Variable-dimension crystal structures**

   -  ``Lattice`` supports the dimensionality of ``mBuild``, which
      means that the systems can be in 1D, 2D, or 3D. Replace the
      necessary vector components with 0 to emulate the dimensionality
      of interest.

-  **Multicomponent crystals**

   -  ``Lattice`` can support an indefinite amount of lattice points in
      its data structure.

   -  The `repeat` cell can be as large as the user defines useful.

   -  The components that occupy the lattice points are ``mbuild.Compound`` objects.

      -  This allows the user to build any system since compounds are
         only a representation of matter in this design.
      -  Molecular crystals, coarse grained, atomic, even alloy crystal
         structures are all supported

-  **Triclinic Lattices**

   -  With full support for triclinic lattices, any crystal structure
      can technically be developed.
   -  Either the user can provide the lattice parameters, or if they
      know the vectors that span the unit cell, that can also be
      provided.

-  **Generation of lattice structure from crystallographic index file** `(CIF) <https://www.iucr.org/resources/cif/documentation>`_ **formats**

   -  Please also see the :ref:`QuickStart_Load_files` section for other ways to load files.

-  *IN PROGRESS* **Template recipes to generate common crystal
   structures (FCC, BCC, HEX, etc)**

   -  This is currently being developed and will be released relatively
      shortly
   -  To generate these structures currently, the user needs to know the
      lattice parameters or lattice vectors that define these units.

``Lattice`` Data Structure Introduction
---------------------------------------

Below we will explore the relevant data structures that are attributes
of the ``Lattice`` class. This information will be essential to build
desired crystal structures.

To begin, we will call the python ``help()`` method to observe the
parameters and attributes of the ``Lattice`` class.

.. code:: ipython3

    import mbuild
    help(mbuild.Lattice)

As we can see, there are quite a few attributes and parameters that make
up this class. There are also a lot of inline examples as well. If you
ever get stuck, remember to use the python built-in ``help()`` method!

-  **Lattice.lattice_spacing**

   This data structure is a (3,) array that details the lengths of the
   repeat cell for the crystal. You can either use a ``numpy`` array
   object, or simply pass in a list and ``Lattice`` will handle the
   rest. Remember that ``mBuild``\ ’s units of length are in nanometers
   [nm]. You must pass in all three lengths, even if they are all
   equivalent. These are the lattice parameters :math:`a, b, c` when
   viewing crystallographic information.

   For Example:

   .. code:: ipython3

       lattice_spacing = [0.5, 0.5, 0.5]

-  **Lattice.lattice_vectors**

   ``lattice_vectors`` is a 3x3 array that defines the vectors that
   encapsulate the repeat cell. This is an optional value that the user
   can pass in to define the cell. Either this must be passed in, or the
   3 Bravais angles of the cell from the lattice parameters must be
   provided. If neither is passed in, the default value are the vectors
   that encase a cubic lattice.

   .. note::
       Most users will **not** have to use these to build their
       lattice structure of interest. It will usually be easier for the
       users to provide the 3 Bravais angles instead. If the user then wants
       the vectors, the ``Lattice`` object will calculate them for the user.

   For example: Cubic Cell

   .. code:: ipython3

       lattice_vectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

-  **Lattice.angles**

   ``angles`` is a (3,) array that defines the three Bravais angles of
   the lattice. Commonly referred to as :math:`\alpha, \beta, \gamma` in
   the definition of the lattice parameters.

   For example: Cubic Cell

   .. code:: ipython3

       angles = [90, 90, 90]

-  **Lattice.lattice_points**

   ``lattice_points`` can be the most common source of confusion when
   creating a crystal structure. In crystallographic terms, this is the
   minimum basis set of points in space that define where the points in
   the lattice exist. This requires that the user does not over define
   the system.

   .. note::
       MIT's OpenCourseWare has an excellent PDF for more information
       `here <https://ocw.mit.edu/courses/earth-atmospheric-and-planetary-sciences/12-108-structure-of-earth-materials-fall-2004/lecture-notes/lec7.pdf>`_

   The other tricky issue that can come up is the data structure itself.
   ``lattice_points`` is a dictionary where the ``dict.key`` items are
   the ``string`` id’s for each basis point. The ``dict.values`` items
   are a nested list of fractional coordinates of the unique lattice
   points in the cell. If you have the same ``Compound`` at multiple
   lattice_points, it is easier to put all those coordinates in a nested
   list under the same ``key`` value. Two examples will be given below,
   both FCC unit cells, one with all the same id, and one with unique
   ids for each lattice_point.

   For Example: FCC All Unique

   .. code:: ipython3

       lattice_points = {'A' : [[0, 0, 0]],
                         'B' : [[0.5, 0.5, 0]],
                         'C' : [[0.5, 0, 0.5]],
                         'D' : [[0, 0.5, 0.5]]
                         }

   For Example: FCC All Same

   .. code:: ipython3

       lattice_points = {'A' : [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]] }


``Lattice`` Public Methods
--------------------------

The ``Lattice`` class also contains methods that are responsible for
applying ``Compounds`` to the lattice points, with user defined cell
replications in the x, y, and z directions.

-  **Lattice.populate(compound_dict=None, x=1, y=1, z=1)**

   This method uses the ``Lattice`` object to place ``Compounds`` at the
   specified ``lattice_points``. There are 4 optional inputs for this
   class.

   -  ``compound_dict`` inputs another dictionary that
      defines a relationship between the ``lattice_points`` and the
      ``Compounds`` that the user wants to populate the lattice with.
      The ``dict.keys`` of this dictionary must be the same as the
      ``keys`` in the ``lattice_points`` dictionary. However, for the
      ``dict.items`` in this case, the ``Compound`` that the user wants
      to place at that lattice point(s) will be used. An example will
      use the FCC examples from above. They have been copied below:

      For Example: FCC All Unique

          .. code:: ipython3

            lattice_points = {'A' : [[0, 0, 0]],
                              'B' : [[0.5, 0.5, 0]],
                              'C' : [[0.5, 0, 0.5]],
                              'D' : [[0, 0.5, 0.5]]
                              }

            # compound dictionary
            a = mbuild.Compound(name='A')
            b = mbuild.Compound(name='B')
            c = mbuild.Compound(name='C')
            d = mbuild.Compound(name='D')
            compound_dict = {'A' : a, 'B' : b, 'C' : c, 'D' : d}


      For Example: FCC All Same

          .. code:: ipython3

            lattice_points = {'A' : [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]] }

            # compound dictionary
            a = mbuild.Compound(name='A')
            compound_dict = {'A' : a}


Example Lattice Systems
-----------------------

Below contains some examples of homogeneous and heterogeneous 2D and 3D
lattice structures using the ``Lattice`` class.

Simple Cubic (SC)
~~~~~~~~~~~~~~~~~

-  Polonium

.. code:: ipython3

    import mbuild as mb
    import numpy as np

    # define all necessary lattice parameters
    spacings = [0.3359, 0.3359, 0.3359]
    angles = [90, 90, 90]
    points = [[0, 0, 0]]

    # define lattice object
    sc_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'Po' : points})

    # define Polonium Compound
    po = mb.Compound(name='Po')

    # populate lattice with compounds
    po_lattice = sc_lattice.populate(compound_dict={'Po' : po}, x=2, y=2, z=2)

    # visualize
    po_lattice.visualize()

.. figure:: ../../images/lattice_SC_polonium_image.png
    :width: 40 %
    :align: center

    **Polonium simple cubic (SC) structure**

Body-centered Cubic (BCC)
~~~~~~~~~~~~~~~~~~~~~~~~~

-  CsCl

.. code:: ipython3

    import mbuild as mb
    import numpy as np

    # define all necessary lattice parameters
    spacings = [0.4123, 0.4123, 0.4123]
    angles = [90, 90, 90]
    point1 = [[0, 0, 0]]
    point2 = [[0.5, 0.5, 0.5]]

    # define lattice object
    bcc_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : point1, 'B' : point2})

    # define Compounds
    cl = mb.Compound(name='Cl')
    cs = mb.Compound(name='Cs')

    # populate lattice with compounds
    cscl_lattice = bcc_lattice.populate(compound_dict={'A' : cl, 'B' : cs}, x=2, y=2, z=2)

    # visualize
    cscl_lattice.visualize()

.. figure:: ../../images/lattice_BCC_CsCl_image.png
    :width: 40 %
    :align: center

    **CsCl body-centered cubic (BCC) structure**

Face-centered Cubic (FCC)
~~~~~~~~~~~~~~~~~~~~~~~~~

-  Cu

.. code:: ipython3

    import mbuild as mb
    import numpy as np

    # define all necessary lattice parameters
    spacings = [0.36149, 0.36149, 0.36149]
    angles = [90, 90, 90]
    points = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]

    # define lattice object
    fcc_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : points})

    # define Compound
    cu = mb.Compound(name='Cu')

    # populate lattice with compounds
    cu_lattice = fcc_lattice.populate(compound_dict={'A' : cu}, x=2, y=2, z=2)

    # visualize
    cu_lattice.visualize()

.. figure:: ../../images/lattice_FCC_Cu_image.png
    :width: 40 %
    :align: center

    **Cu face-centered cubic (FCC) structure**

Diamond (Cubic)
~~~~~~~~~~~~~~~

-  Si

.. code:: ipython3

    import mbuild as mb
    import numpy as np

    # define all necessary lattice parameters
    spacings = [0.54309, 0.54309, 0.54309]
    angles = [90, 90, 90]
    points = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5],
              [0.25, 0.25, 0.75], [0.25, 0.75, 0.25], [0.75, 0.25, 0.25], [0.75, 0.75, 0.75]]

    # define lattice object
    diamond_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : points})

    # define Compound
    si = mb.Compound(name='Si')

    # populate lattice with compounds
    si_lattice = diamond_lattice.populate(compound_dict={'A' : si}, x=2, y=2, z=2)

    # visualize
    si_lattice.visualize()

.. figure:: ../../images/lattice_Diamond_cubic_Si_image.png
    :width: 40 %
    :align: center

    **Si diamond (Cubic) structure**

Graphene (2D)
~~~~~~~~~~~~~

-  C

.. code:: ipython3

    import mbuild as mb
    import numpy as np


    # define all necessary lattice parameters
    spacings = [0.246, 0.246, 0.335]
    angles = [90, 90, 120]
    points = [[0, 0, 0], [1/3, 2/3, 0]]

    # define lattice object
    graphene_lattice = mb.Lattice(lattice_spacing=spacings, angles=angles, lattice_points={'A' : points})

    # define Compound
    c = mb.Compound(name='C')

    # populate lattice with compounds
    graphene = graphene_lattice.populate(compound_dict={'A' : c}, x=5, y=5, z=1)

    # visualize
    graphene.visualize()

.. figure:: ../../images/lattice_graphene_2D_image.png
    :width: 40 %
    :align: center

    **Graphene (2D) structure**
Monolayer: Complex hierarchies, patterns, tiling and writing to files
---------------------------------------------------------------------

**Note**: mBuild expects all distance units to be in nanometers.

In this example, we’ll cover assembling more complex hierarchies of
components using patterns, tiling and how to output systems to files. To
illustrate these concepts, let’s build an alkane monolayer on a
crystalline substrate.

First, let’s build our monomers and functionalized them with a silane
group which we can then attach to the substrate. The ``Alkane`` example
uses the ``polymer`` tool to combine ``CH2`` and ``CH3`` repeat units.
You also have the option to cap the front and back of the chain or to
leave a ``CH2`` group with a dangling port. The ``Silane`` compound is a
Si(OH)2 group with two ports facing out from the central Si. Lastly, we
combine ``alkane`` with ``silane`` and add a label to ``AlkylSilane``
which points to, ``silane['down']``. This allows us to reference it
later using ``AlkylSilane['down']`` rather than
``AlkylSilane['silane']['down']``.

**Note:** In ``Compounds`` with multiple ``Ports``, by convention, we
try to label every ``Port`` successively as ‘up’, ‘down’, ‘left’,
‘right’, ‘front’, ‘back’ which should roughly correspond to their
relative orientations. This is a bit tricky to enforce because the
system is so flexible so use your best judgement and try to be
consistent! The more components we collect in our library with the same
labeling conventions, the easier it becomes to build ever more complex
structures.

.. code:: ipython3

    import mbuild as mb

    from mbuild.lib.recipes import Alkane
    from mbuild.lib.moieties import Silane


    class AlkylSilane(mb.Compound):
        """A silane functionalized alkane chain with one Port. """
        def __init__(self, chain_length):
            super(AlkylSilane, self).__init__()

            alkane = Alkane(chain_length, cap_end=False)
            self.add(alkane, 'alkane')
            silane = Silane()
            self.add(silane, 'silane')
            mb.force_overlap(self['alkane'], self['alkane']['down'], self['silane']['up'])

            # Hoist silane port to AlkylSilane level.
            self.add(silane['down'], 'down', containment=False)

.. code:: ipython3

    AlkylSilane(5).visualize()

Now let’s create a substrate to which we can later attach our monomers:

.. code:: ipython3

    import mbuild as mb
    from mbuild.lib.surfaces import Betacristobalite

    surface = Betacristobalite()
    tiled_surface = mb.lib.recipes.TiledCompound(surface, n_tiles=(2, 1, 1))

Here we’ve imported a beta-cristobalite surface from our component
library. The ``TiledCompound`` tool allows you replicate any
``Compound`` in the x-, y- and z-directions by any number of times - 2,
1 and 1 for our case.

Next, let’s create our monomer and a hydrogen atom that we’ll place on
unoccupied surface sites:

.. code:: ipython3

    from mbuild.lib.atoms import H
    alkylsilane = AlkylSilane(chain_length=10)
    hydrogen = H()

Then we need to tell mBuild how to arrange the chains on the surface.
This is accomplished with the “pattern” tools. Every pattern is just a
collection of points. There are all kinds of patterns like spherical,
2D, regular, irregular etc. When you use the ``apply_pattern`` command,
you effectively superimpose the pattern onto the host compound, mBuild
figures out what the closest ports are to the pattern points and then
attaches copies of the guest onto the binding sites identified by the
pattern:

.. code:: ipython3

    pattern = mb.Grid2DPattern(8, 8)  # Evenly spaced, 2D grid of points.

    # Attach chains to specified binding sites. Other sites get a hydrogen.
    chains, hydrogens = pattern.apply_to_compound(host=tiled_surface, guest=alkylsilane, backfill=hydrogen)

Also note the ``backfill`` optional argument which allows you to place a
different compound on any unused ports. In this case we want to backfill
with hydrogen atoms on every port without a chain.

.. code:: ipython3

    monolayer = mb.Compound([tiled_surface, chains, hydrogens])
    monolayer.visualize() # Warning: may be slow in IPython notebooks

.. code:: ipython3

    # Save as .mol2 file
    monolayer.save('monolayer.mol2', overwrite=True)

``lib.recipes.monolayer.py`` wraps many these functions into a simple,
general class for generating the monolayers, as shown below:

.. code:: ipython3

    from mbuild.lib.recipes import Monolayer

    monolayer = Monolayer(fractions=[1.0], chains=alkylsilane, backfill=hydrogen,
                          pattern=mb.Grid2DPattern(n=8, m=8),
                          surface=surface, tile_x=2, tile_y=1)
    monolayer.visualize()
Point Particles: Basic system initialization
============================================

**Note**: mBuild expects all distance units to be in nanometers.

This tutorial focuses on the usage of basic system initialization
operations, as applied to simple point particle systems (i.e., generic
Lennard-Jones particles rather than specific atoms).

The code below defines several point particles in a cubic arrangement.
Note, the color and radius associated with a Particle name can be set
and passed to the visualize command. Colors are passed in hex format
(see http://www.color-hex.com/color/bfbfbf).

.. code:: ipython3

    import mbuild as mb

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            lj_particle1 = mb.Particle(name='LJ', pos=[0, 0, 0])
            self.add(lj_particle1)

            lj_particle2 = mb.Particle(name='LJ', pos=[1, 0, 0])
            self.add(lj_particle2)

            lj_particle3 = mb.Particle(name='LJ', pos=[0, 1, 0])
            self.add(lj_particle3)

            lj_particle4 = mb.Particle(name='LJ', pos=[0, 0, 1])
            self.add(lj_particle4)

            lj_particle5 = mb.Particle(name='LJ', pos=[1, 0, 1])
            self.add(lj_particle5)

            lj_particle6 = mb.Particle(name='LJ', pos=[1, 1, 0])
            self.add(lj_particle6)

            lj_particle7 = mb.Particle(name='LJ', pos=[0, 1, 1])
            self.add(lj_particle7)

            lj_particle8 = mb.Particle(name='LJ', pos=[1, 1, 1])
            self.add(lj_particle8)


    monoLJ = MonoLJ()
    monoLJ.visualize()

While this would work for defining a single molecule or very small
system, this would not be efficient for large systems. Instead, the
clone and translate operator can be used to facilitate automation.
Below, we simply define a single prototype particle (lj_proto), which we
then copy and translate about the system.

Note, mBuild provides two different translate operations, “translate”
and “translate_to”. “translate” moves a particle by adding the vector
the original position, whereas “translate_to” move a particle to the
specified location in space. Note, “translate_to” maintains the internal
spatial relationships of a collection of particles by first shifting the
center of mass of the collection of particles to the origin, then
translating to the specified location. Since the lj_proto particle in
this example starts at the origin, these two commands produce identical
behavior.

.. code:: ipython3

    import mbuild as mb

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            for i in range(0,2):
                for j in range(0,2):
                    for k in range(0,2):
                        lj_particle = mb.clone(lj_proto)
                        pos = [i,j,k]
                        lj_particle.translate(pos)
                        self.add(lj_particle)

    monoLJ = MonoLJ()
    monoLJ.visualize()

To simplify this process, mBuild provides several build-in patterning
tools, where for example, Grid3DPattern can be used to perform this same
operation. Grid3DPattern generates a set of points, from 0 to 1, which
get stored in the variable “pattern”. We need only loop over the points
in pattern, cloning, translating, and adding to the system. Note,
because Grid3DPattern defines points between 0 and 1, they must be
scaled based on the desired system size, i.e., pattern.scale(2).

.. code:: ipython3

    import mbuild as mb

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern = mb.Grid3DPattern(2, 2, 2)
            pattern.scale(2)

            for pos in pattern:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)

    monoLJ = MonoLJ()
    monoLJ.visualize()

Larger systems can therefore be easily generated by toggling the values
given to Grid3DPattern. Other patterns can also be generated using the
same basic code, such as a 2D grid pattern:

.. code:: ipython3

    import mbuild as mb

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern = mb.Grid2DPattern(5, 5)
            pattern.scale(5)

            for pos in pattern:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)

    monoLJ = MonoLJ()
    monoLJ.visualize()

Points on a sphere can be generated using SpherePattern. Points on a
disk using DisKPattern, etc.

Note to show both simultaneously, we shift the x-coordinate of Particles
in the sphere by -1 (i.e., pos[0]-=1.0) and +1 for the disk (i.e,
pos[0]+=1.0).

.. code:: ipython3

    import mbuild as mb

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern_sphere = mb.SpherePattern(200)
            pattern_sphere.scale(0.5)

            for pos in pattern_sphere:
                lj_particle = mb.clone(lj_proto)
                pos[0]-=1.0
                lj_particle.translate(pos)
                self.add(lj_particle)

            pattern_disk = mb.DiskPattern(200)
            pattern_disk.scale(0.5)
            for pos in pattern_disk:
                lj_particle = mb.clone(lj_proto)
                pos[0]+=1.0
                lj_particle.translate(pos)
                self.add(lj_particle)

    monoLJ = MonoLJ()
    monoLJ.visualize()

We can also take advantage of the hierachical nature of mBuild to
accomplish the same task more cleanly. Below we create a component that
corresponds to the sphere (class SphereLJ), and one that corresponds to
the disk (class DiskLJ), and then instantiate and shift each of these
individually in the MonoLJ component.

.. code:: ipython3

    import mbuild as mb

    class SphereLJ(mb.Compound):
        def __init__(self):
            super(SphereLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern_sphere = mb.SpherePattern(200)
            pattern_sphere.scale(0.5)

            for pos in pattern_sphere:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)

    class DiskLJ(mb.Compound):
        def __init__(self):
            super(DiskLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern_disk = mb.DiskPattern(200)
            pattern_disk.scale(0.5)
            for pos in pattern_disk:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)


    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()

            sphere = SphereLJ();
            pos=[-1, 0, 0]
            sphere.translate(pos)
            self.add(sphere)

            disk = DiskLJ();
            pos=[1, 0, 0]
            disk.translate(pos)
            self.add(disk)


    monoLJ = MonoLJ()
    monoLJ.visualize()

Again, since mBuild is hierarchical, the pattern functions can be used
to generate large systems of any arbitary component. For example, we can
replicate the SphereLJ component on a regular array.

.. code:: ipython3

    import mbuild as mb

    class SphereLJ(mb.Compound):
        def __init__(self):
            super(SphereLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern_sphere = mb.SpherePattern(13)
            pattern_sphere.scale(0.1)

            for pos in pattern_sphere:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)
    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            sphere = SphereLJ();

            pattern = mb.Grid3DPattern(3, 3, 3)
            pattern.scale(2)

            for pos in pattern:
                lj_sphere = mb.clone(sphere)
                lj_sphere.translate_to(pos)
                #shift the particle so the center of mass
                #of the system is at the origin
                lj_sphere.translate([-5,-5,-5])

                self.add(lj_sphere)

    monoLJ = MonoLJ()
    monoLJ.visualize()

Several functions exist for rotating compounds. For example, the spin
command allows a compound to be rotated, in place, about a specific axis
(i.e., it considers the origin for the rotation to lie at the compound’s
center of mass).

.. code:: ipython3

    import mbuild as mb
    import random
    from numpy import pi


    class CubeLJ(mb.Compound):
        def __init__(self):
            super(CubeLJ, self).__init__()
            lj_proto = mb.Particle(name='LJ', pos=[0, 0, 0])

            pattern = mb.Grid3DPattern(2, 2, 2)
            pattern.scale(0.2)

            for pos in pattern:
                lj_particle = mb.clone(lj_proto)
                lj_particle.translate(pos)
                self.add(lj_particle)

    class MonoLJ(mb.Compound):
        def __init__(self):
            super(MonoLJ, self).__init__()
            cube_proto = CubeLJ();

            pattern = mb.Grid3DPattern(3, 3, 3)
            pattern.scale(2)
            rnd = random.Random()
            rnd.seed(123)

            for pos in pattern:
                lj_cube = mb.clone(cube_proto)
                lj_cube.translate_to(pos)
                #shift the particle so the center of mass
                #of the system is at the origin
                lj_cube.translate([-5,-5,-5])
                lj_cube.spin( rnd.uniform(0, 2 * pi), [1, 0, 0])
                lj_cube.spin(rnd.uniform(0, 2 * pi), [0, 1, 0])
                lj_cube.spin(rnd.uniform(0, 2 * pi), [0, 0, 1])

                self.add(lj_cube)

    monoLJ = MonoLJ()
    monoLJ.visualize()

Configurations can be dumped to file using the save command; this takes
advantage of MDTraj and supports a range of file formats (see
http://MDTraj.org).

.. code:: ipython3

    #save as xyz file
    monoLJ.save('output.xyz')
    #save as mol2
    monoLJ.save('output.mol2')
-----------------------
Tutorials
-----------------------
.. only:: html

    The following tutorials are available either as html or interactive jupyter notebooks.


.. toctree::

    	tutorial_methane
    	tutorial_ethane
    	tutorial_monolayer
    	tutorial_simple_LJ
    	tutorial_polymers
	tutorial_lattice
