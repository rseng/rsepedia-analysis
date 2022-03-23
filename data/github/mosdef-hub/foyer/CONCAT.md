Contributions are welcomed via [pull requests on GitHub](https://github.com/mosdef-hub/foyer/pulls). Developers and/or
users will review requested changes and make comments. The rest of this file will serve as a set of general guidelines
for contributors.

# Features

## Implement functionality in a general and flexible fashion

Foyer is designed to be general and flexible, not limited to single force fields, simulation engines, or  simulation methods.
Additions to core features should attempt to provide something that is applicable to a variety of use-cases and not targeted
at only the focus area of your research. However, some specific features targeted toward a limited use case may be
appropriate. Speak to the developers before writing your code and they will help you make design choices that allow flexibility.

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

All new functionality in Foyer should be tested with automatic unit tests that execute in a few seconds. These tests
should attempt to cover all options that the user can select. All or most of the added lines of source code should be
covered by unit test(s). We currently use [pytest](https://docs.pytest.org/en/latest/), which can be executed simply by calling
`pytest` from the root directory of the package. Additions to force field files should include test molecules that encompass
the added or modified atom types or functionality.
### Foyer: A package for atom-typing as well as applying and disseminating forcefields

[![Gitter chat](https://badges.gitter.im/mosdef-hub/gitter.svg)](https://gitter.im/mosdef-hub/Lobby)
[![AZP Build Status](https://dev.azure.com/mosdef/mosdef/_apis/build/status/mosdef-hub.foyer?branchName=master)](https://dev.azure.com/mosdef/mosdef/_build/latest?definitionId=2&branchName=master)
[![Anaconda Badge](https://anaconda.org/conda-forge/foyer/badges/version.svg)](https://anaconda.org/conda-forge/foyer)
[![codecov](https://codecov.io/gh/mosdef-hub/foyer/branch/master/graph/badge.svg)](https://codecov.io/gh/mosdef-hub/foyer)
[![DOI](https://zenodo.org/badge/34077879.svg)](https://zenodo.org/badge/latestdoi/34077879)


## Overview
Foyer is an open-source Python tool for defining and applying force field atom-typing
rules in a format that is both human- and machine-readable.  It parametrizes chemical topologies,
generating, syntactically correct input files for various simulation engines. Foyer provides a framework for force field
dissemination, helping to eliminate ambiguity in atom-typing and improving reproducibility
(for more information, see [our paper](https://www.sciencedirect.com/science/article/pii/S0927025619303040) or its corresponding [pre-print](https://arxiv.org/pdf/1812.06779.pdf)).

#### Foyer within the MoSDeF Ecosystem
<p align="center">
  <img src="docs/images/mosdef_graphic_foyer.svg?raw=true" alt="Foyer within the MoSDeF Ecosystem" width="500" height="500"/>
</p>

Foyer defines force fields in an XML format, where SMARTS strings are used to define the chemical context
of a particular atom type and “overrides” are used to set rule precedence, rather than a rigid hierarchical scheme.
Foyer builds upon the [OpenMM .xml force field](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields)
file, annotated with SMARTS-based atomtypes, e.g.:

```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
  <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
 </AtomTypes>
</ForceField>
```

Foyer can apply the forcefield to arbitrary chemical topologies. We currently support:

* [OpenMM.Topology](https://github.com/openmm/openmm/blob/7.6.0/wrappers/python/openmm/app/topology.py#L70)
* [ParmEd.Structure](http://parmed.github.io/ParmEd/html/structure.html)
* [mBuild.Compound](http://mosdef-hub.github.io/mbuild/data_structures.html)
* [gmso.Topology](https://gmso.mosdef.org/en/stable/data_structures.html#gmso.Topology)
* [openff.tookit.topology.Topology](https://open-forcefield-toolkit.readthedocs.io/en/0.9.2/api/generated/openff.toolkit.topology.Topology.html#openff-toolkit-topology-topology)

Application of a force field can be as simple as:
```python
from foyer import Forcefield
import parmed as pmd

untyped_ethane = pmd.load_file('ethane.mol2', structure=True)
oplsaa = Forcefield(forcefield_files='oplsaa.xml')
ethane = oplsaa.apply(untyped_ethane)

# Save to any format supported by ParmEd
ethane.save('ethane.top')
ethane.save('ethane.gro')
```

The `Foyer` package is part of the [Molecular Simulation Design Framework (MoSDeF) project](http://mosdef.org/).
Libraries in the MoSDeF ecosystem are designed to provide utilities neccessary to streamline
a researcher's simulation workflow. When setting up simulation studies,
we also recommend users to follow the [TRUE](https://www.tandfonline.com/doi/full/10.1080/00268976.2020.1742938)
(Transparent, Reproducible, Usable-by-others, and Extensible) standard, which is a set of common
practices meant to improve the reproducibility of computational simulation research.

## Getting started

#### Quick Start with Docker
To use `foyer` in a jupyter-notebook that runs from a docker container with all the dependencies installed use the following command:

```sh
$ docker pull mosdef/foyer:latest
$ docker run -it --name foyer -p 8888:8888 mosdef/foyer:latest\
  /opt/conda/envs/foyer-docker/bin/jupyter notebook --ip="*"
```

Alternatively, you can also start a Bourne shell directly:
```sh
$ docker run -it --name foyer mosdef/foyer:latest
```

To learn more about using `foyer` with docker, please refer to the documentation [here](https://foyer.mosdef.org/en/latest/docker.html) .


#### Getting started with SMARTS-based atom-typing
* [SMARTS-based atomtyping](docs/source/topic_guides/smarts.rst)
* [Supported SMARTS Grammar](https://github.com/mosdef-hub/foyer/issues/63)

#### Defining force fields:
* [Defining force field parameters](docs/source/topic_guides/parameter_definitions.rst)
* [Force field file validation](docs/source/reference/validation.rst)


#### Example foyer force field files:
Foyer currently includes a subset of the OPLS AA and TraPPE forcefields, currently part of the source distribution:
* https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields

Additional example force field XML files:
* https://github.com/chrisiacovella/OPLSaa_perfluoroalkanes
* https://github.com/mosdef-hub/forcefield_perfluoroethers
* https://github.com/summeraz/OPLSaa_alkylsilanes

Example template for disseminating force fields:
* https://github.com/mosdef-hub/forcefield_template


#### Using Foyer to perform atom typing:
* [Basic usage examples](docs/source/topic_guides/usage_examples.rst)
* [Detailed Jupyter notebook tutorials, including integration with mBuild](https://github.com/mosdef-hub/foyer_tutorials)
* [Jupyter notebook tutorials](https://github.com/mosdef-hub/foyer_tutorials/tree/master), from [our paper](https://arxiv.org/abs/1812.06779)

### Documentation:
* Documentation website: http://foyer.mosdef.org

### Installation instructions
* [Installation instructions](docs/source/getting_started/install.rst)

### Citing Foyer:
* If you use this package, please cite [our paper](https://www.sciencedirect.com/science/article/pii/S0927025619303040) published in [Computational Materials Science](https://www.journals.elsevier.com/computational-materials-science).
* This manuscript is also available in its pre-print form on [arxiv](https://arxiv.org/pdf/1812.06779.pdf)
* The paper and examples in this work were developed for tag [paper_COMMAT_2019](https://github.com/mosdef-hub/foyer/tree/paper_COMMAT_2019)


* Please also cite the github repository, https://github.com/mosdef-hub/foyer

#### [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.

This material is based upon work supported by the National Science Foundation under grants NSF ACI-1047828 and NSF ACI-1535150. Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
### PR Summary:

### PR Checklist
------------
 - [ ] Includes appropriate unit test(s)
 - [ ] Appropriate docstring(s) are added/updated
 - [ ] Code is (approximately) PEP8 compliant
 - [ ] Issue(s) raised/addressed?
---
name: Feature request
about: Suggest an improvement to Foyer

---

**Describe the behavior you would like added to Foyer**
A clear and concise description of what the proposed idea is.

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Questions/discussions
about: For usage questions, important notes, and other discussions.

---

**A summary of the question or discussion topic.**
---
name: Bug report
about: Report a bug in Foyer

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

- Which version of Foyer are you using? (`python -c "import foyer; print(foyer.__version__)"`)
- Which version of Python (`python --version`)?
- Which operating system?
### Usage Examples

Foyer supports atomtyping of both all-atom and coarse-grained molecular systems, and
also allows for separate force fields to be used to atom-type separate portions of a
molecular system.

#### Creating a box of ethane
Here we use mBuild to construct a box filled with ethane molecules and use Foyer to
atom-type the system, applying the OPLS force field, and save to run-able GROMACS
files.
```python
import mbuild as mb
from mbuild.examples import Ethane

ethane_box = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
ethane_box.save('ethane-box.gro')
ethane_box.save('ethane-box.top', forcefield_name='oplsaa')
```
----

#### Creating a box of coarse-grained ethane
Again we will use mBuild to construct a box filled with ethane molecules.  However,
now we will model ethane using a united-atom description and apply the TraPPE force
field during atom-typing.  Note how particle names are prefixed with an underscore so
that Foyer knows these particles are non-atomistic.
```python
import mbuild as mb

ethane_UA = mb.Compound()
ch3_1 = mb.Particle(name='_CH3', pos=[0, 0, 0])
ch3_2 = mb.Particle(name='_CH3', pos=[0.15, 0, 0])
ethane_UA.add([ch3_1, ch3_2])
ethane_UA.add_bond((ch3_1, ch3_2))

ethane_UA_box = mb.fill_box(ethane_UA, 100, box=[2, 2, 2])
ethane_UA_box.save('ethane-UA-box.gro')
ethane_UA_box.save('ethane-UA-box.top', forcefield_name='trappe-ua')
```
----
#### Combining force fields
In some instances, the use of multiple force fields may be desired to describe a
molecular system.  For example, the user may want to use one force field for a
surface and another for a fluid in the same system.  Foyer supports this
functionality by allowing the user to separately atom-type parts of a system.  In
this example, we take a system featuring bulk united atom ethane above a silica
surface and apply the OPLS force field to the surface and the TraPPE force field to
the ethane.  The two atomtyped Parmed structures are then combined using a simple
'\+' operator and can be saved to Gromacs files.
```python
from foyer import Forcefield
from foyer.tests.utils import get_fn
import mbuild as mb
from mbuild.examples import Ethane
from mbuild.lib.atoms import H
from mbuild.lib.bulk_materials import AmorphousSilica

interface = mb.SilicaInterface(bulk_silica=AmorphousSilica())
interface = mb.Monolayer(surface=interface, chains=H(), guest_port_name='up')

box = mb.Box(mins=[0, 0, max(interface.xyz[:,2])],
             maxs=interface.periodicity + [0, 0, 4])

ethane_box = mb.fill_box(compound=Ethane(), n_compounds=200, box=box)

opls = Forcefield(name='oplsaa')
opls_silica = Forcefield(forcefield_files=get_fn('opls-silica.xml'))
ethane_box = opls.apply(ethane_box)
interface = opls_silica.apply(interface)

system = interface + ethane_box

system.save('ethane-silica.gro')
system.save('ethane-silica.top')
```
Install from conda:
```bash
conda install -c mosdef foyer
```

Install an editable version from source:
```bash
git clone https://github.com/mosdef-hub/foyer.git
cd foyer
pip install -e .
```
### SMARTS-based atomtyping

Foyer allows users to describe atomtypes using a modified version of
[SMARTS](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html)
You may already be familiar with
[SMILES](https://www.wikiwand.com/en/Simplified_molecular-input_line-entry_system)
representations for describing chemical structures. SMARTS is a straightforward
extension of this notation.

#### Basic usage
Consider the following example defining the OPLS-AA atomtypes for a methyl group
carbon and its hydrogen atoms:
```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
  <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
 </AtomTypes>
</ForceField>
```

This `.xml` format is an extension of the [OpenMM force field format](http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields)
The above example utilizes two additional `.xml` attributes supported by foyer:
`def` and `desc`. The atomtype that we are attempting to match is always the
__first__ token in the SMARTS string, in the above example, `[C;X4]` and `H`.
The `opls_135` (methyl group carbon) is defined by a SMARTS
string indicated a carbon with 4 bonds, a carbon neighbor and 3
hydrogen neighbors. The `opls_140` (alkane hydrogen) is defined simply as a
hydrogen atom bonded to a carbon with 4 bonds.


#### Overriding atomtypes
When multiple atomtype definitions can apply to a given atom, we must establish
precedence between those definitions. Many other atomtypers determine rule
precedence by placing more specific rules first and evaluate those in sequence,
breaking out of the loop as soon as a match is found.

While this approach works, it becomes more challenging to maintain the correct
ordering of rules as the number of atomtypes grows. Foyer iteratively runs all
rules on all atoms and each atom maintains a whitelist (rules that apply) and a
blacklist (rules that have been superceded by another rule). The set difference
between the white- and blacklists yields the correct atomtype if the force field
is implemented correctly.

We can add a rule to a blacklist using the `overrides` attribute in the `.xml`
file. To illustrate an example where overriding can be used consider the
following types describing alkenes and benzene:

```xml
<ForceField>
 <AtomTypes>
  <Type name="opls_141" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)C" desc="alkene C (R2-C=)"/>
  <Type name="opls_142" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)H" desc="alkene C (RH-C=)"/>
  <Type name="opls_144" class="HC" element="H" mass="1.00800"  def="[H][C;X3]" desc="alkene H"/>
  <Type name="opls_145" class="CA" element="C" mass="12.01100" def="[C;X3;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1" overrides="opls_141,opls_142"/>
  <Type name="opls_146" class="HA" element="H" mass="1.00800"  def="[H][C;%opls_145]" overrides="opls_144" desc="benzene H"/>
 </AtomTypes>
</ForceField>
```

If we're atomtyping a benzene molecule, the carbon atoms will match the SMARTS
patterns for both `opls_142` and `opls_145`. Without the `overrides` attribute,
foyer will notify you that multiple atomtypes were found for each carbon.
Providing the `overrides` indicates that if the `opls_145` pattern matches, it
should supercede the specified rules.

#### Current Grammar Supported
We currently do not (yet) support all of [SMARTS' features](http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html). [Here](https://github.com/mosdef-hub/foyer/issues/63) we're keeping track of which portions are supported.
### Validation of force field files

Foyer performs several validation steps to help prevent malformed force field
files and SMARTS strings from making it into production code. Our goal is to
provide human readable error messages that users who may not be intimately
familiar with XML or our SMARTS parsing grammar can easily act upon.

However, if you receive any unclear error messages or warnings we strongly
encourage you to
__[submit an issue](https://github.com/mosdef-hub/foyer/issues/new)__
detailing the error message you received and, if possible, attach a minimal
example of the force field file that created the problem.

#### XML schema

As a first line of defense, any force field files loaded by foyer is validated
by this [XML schema definition](../foyer/forcefields/ff.xsd). Here we enforce
which elements (e.g. `HarmonicBondForce`) are valid and how their attributes
should be formatted. Additionally, the schema ensures that atomtypes are not 1)
defined more than once and that 2) atomtypes referenced in other sections
are actually defined in the `<AtomTypes>` element.

#### SMARTS validation
All SMARTS strings used to define atomtypes are parsed. Parsing errors are
captured and re-raised with error messages that allow you to pin point the
location of the problem in the XML file and within the SMARTS string. Wherever
possible, we attempt to provide helpful hints and we welcome any contributions
that help improve the clarity of our error messages.

Additionally, we ensure that any atomtypes referenced using the
[`%type` or `overrides` syntax](smarts.md) are actually defined in the
`<AtomTypes>` element.
# Forcefield Changelog
----
## OPLS-AA

v0.0.1 - April 15, 2021
 - started versioning of forcefield xmls
v0.0.2 - June 24, 2021
 - update SMARTS string for:
    -   opls_182 (from `[C;X4]([O;%opls_180])(H)(H)` to `[C;X4]([O;%opls_180])(H)(H)C`)
    -   opls_282(from `HC[C;%opls_277,%opls_280]` to `HC[C;%opls_277,%opls_280,%opls_465;!%opls_267]`)
    -   opls_468 (from `[C;X4]([O;%opls_467])(H)(H)` to `[C;X4]([O;%opls_467])(H)(H)H`)
    -   opls_469 (from `H[C;%opls_468]` to `H[C;%opls_468,%opls_490]`)
    -   opls_490 (`[C;X4]([O;%opls_467])(H)(H)C`)

 - update overrides for:
    -   opls_279 (from `"opls_144"` to `"opls_185, opls_144"`)
    -   opls_465 (from `""` to `"opls_277"`)
    -   opls_465 (from `""` to `"opls_278"`)

  - update references
    - opls_465 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_466 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_467 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_468 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)
    - opls_469 (from `"10.1021/ja9539195"` to `"10.1002/jcc.1092"`)

v0.0.3 - August 7, 2021
 - update SMARTS string for opls_154 (from `[O;X2]H` to `[O;X2](H)([!H])`)

----
## Trappe-UA

v0.0.1 - April 15, 2021
 - started versioning of forcefield xmls
v0.0.2 - August 9, 2021
 - Updated `combining_rule` to `lorentz`
Developer Notes / Tools / License
=================================

Assorted notes for developers.

Most of these tools were adapted from MDTraj and are released under the following
 license:

License
-------
Copyright (c) 2012-2015 Stanford University and the Authors
All rights reserved.

Redistribution and use of all files in this folder (devtools) and (../basesetup.py,
../setup.py) files in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

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
Foyer
======
*A package for atom-typing as well as applying and disseminating force fields*

|License|
|Citing|
|Anaconda|
|CodeCov|
|Azure|

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
.. |Citing| image:: https://img.shields.io/badge/cite-foyer-green
   :target: reference/citing.html
.. |Anaconda| image:: https://anaconda.org/conda-forge/foyer/badges/version.svg
   :target: https://anaconda.org/conda-forge/foyer
.. |Codecov| image:: https://codecov.io/gh/mosdef-hub/foyer/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/mosdef-hub/foyer
.. |Azure| image:: https://dev.azure.com/mosdef/mosdef/_apis/build/status/mosdef-hub.foyer?branchName=master
   :target: https://dev.azure.com/mosdef/mosdef/_build/latest?definitionId=2&branchName=master


Overview
~~~~~~~~

**Foyer** is an open-source Python tool that provides a framework for the application
and dissemination of classical molecular modeling force fields. Importantly,
it enables users to define and apply atom-typing rules in a format that is
simultaneously human- and machine-readable. A primary goal of **foyer**
is to eliminate ambiguity in the atom-typing and force field application steps of
molecular simulations in order to improve reproducibility. Foyer force fields are
defined in an XML format derived from the
`OpenMM XML specification <http://docs.openmm.org/latest/userguide/application.html#basic-concepts>`_.
`SMARTS strings <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_
are used to define the chemical context of each atom type and "overrides" are used
to define clear precedence of different atom types. **Foyer** is designed to
be compatible with the other tools in the
`Molecular Simulation Design Framework (MoSDeF) ecosystem <https://mosdef.org>`_.


Resources
~~~~~~~~~

* :doc:`Installation guide <getting_started/install>`: Instructions for installing foyer.
* :doc:`Quickstart <getting_started/quickstart>`: A brief introduction to foyer.
* `MoSDeF  <https://mosdef.org>`_: Learn more about the **Mo**\ lecular **S**\ imulation **De**\ sign **F**\ ramework.
* `Foyer paper <https://www.sciencedirect.com/science/article/pii/S0927025619303040>`_: The journal article describing foyer.
* `GitHub repository <https://github.com/mosdef-hub/foyer>`_: Download the source code or contribute to the development of foyer.
* `Issue Tracker <https://github.com/mosdef-hub/foyer/issues>`_: Report issues and request features.


Citation
~~~~~~~~

If you use foyer in your research, please cite the
`foyer paper <https://www.sciencedirect.com/science/article/pii/S0927025619303040>`_.
See :doc:`here <reference/citing>` for details.


Installation
~~~~~~~~~~~~

Complete installation instructions are :doc:`here <getting_started/install>`.
A conda installation is available:

.. code-block:: bash

    conda create --name foyer foyer -c conda-forge

Example
~~~~~~~

Annotate an `OpenMM .xml force
field <http://docs.openmm.org/latest/userguide/application.html#creating-force-fields>`__
file with SMARTS-based atomtypes:

.. code:: xml

    <ForceField>
     <AtomTypes>
      <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
      <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
     </AtomTypes>
    </ForceField>

Apply the forcefield to arbitrary chemical topologies. We currently
support:

-  `OpenMM.Topology <http://docs.openmm.org/latest/api-python/generated/simtk.openmm.app.topology.Topology.html#>`__
-  `ParmEd.Structure <http://parmed.github.io/ParmEd/html/structure.html>`__
-  `mBuild.Compound <http://mosdef-hub.github.io/mbuild/data_structures.html>`__

.. code:: python

    from foyer import Forcefield
    import parmed as pmd

    untyped_ethane = pmd.load_file('ethane.mol2', structure=True)
    oplsaa = Forcefield(forcefield_files='oplsaa.xml')
    ethane = oplsaa.apply(untyped_ethane)

    # Save to any format supported by ParmEd
    ethane.save('ethane.top')
    ethane.save('ethane.gro')

Getting started?
~~~~~~~~~~~~~~~~

Check out our example template for disseminating force fields:
https://github.com/mosdef-hub/forcefield_template


Credits
~~~~~~~

This material is based upon work supported by the National Science
Foundation under grants NSF ACI-1047828 and NSF ACI-1535150. Any
opinions, findings, and conclusions or recommendations expressed in this
material are those of the author(s) and do not necessarily reflect the
views of the National Science Foundation.

Table of Contents
~~~~~~~~~~~~~~~~~

.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    getting_started/install
    getting_started/docker
    getting_started/quickstart

.. toctree::
    :maxdepth: 2
    :caption: Topic Guides

    topic_guides/smarts
    topic_guides/parameter_definitions
    topic_guides/ffapply
    topic_guides/usage_examples
    topic_guides/paper_examples

.. toctree::
    :maxdepth: 2
    :caption: Reference

    reference/units
    reference/validation
    reference/ffclass
    reference/citing
    reference/license
Quickstart
==========

**Foyer** is distributed with partial support for the OPLS-AA
force field. First, we will load a PDB file
(:download:`here <../files/ethane.pdb>`) with ``parmed``
and then load the OPLS-AA force field with foyer.

.. code:: python

    import foyer
    import parmed

    mol = parmed.load("ethane.pdb")
    ff = foyer.forcefield.load_OPLSAA()

    mol_ff = ff.apply(mol)

    mol_ff.save("ethane.gro")
    mol_ff.save("ethane.top")


The first step loads ``ethane.pdb``. ``mol`` is a ``parmed.Structure``.
Next, we load the OPLS force field as a ``foyer.forcefield`` object.
Then, we apply the force field to the molecule. ``mol_ff`` is also
a ``parmed.Structure``, but it now contains all of the force field
parameters (Lennard-Jones parameters, partial charges, bond, angle,
dihedral parameters). Finally, we save the parameterized structure
to input formats for GROMACS. In this example, ``ff`` has the force field parameters
for the OPLS-AA force field. You can view the XML file with all of the
parameters `here <https://github.com/mosdef-hub/foyer/blob/master/foyer/forcefields/xml/oplsaa.xml>`_.
As you can see from the XML file, SMARTS
strings have been added for many of the OPLS-AA atom types. You
should always verify that the atom-typing is performed correctly for
a single molecule before applying the force field to your entire system.
If you wish to add supported molecules to the OPLS-AA force field,
you can view `this issue <https://github.com/mosdef-hub/foyer/issues/314>`_

The real power of **foyer** is to build and distribute XML files with your
own parameter sets. Here we show an example of an XML file for
pentafluoroethane. This is a custom parameter set taken from
`this paper <https://arxiv.org/abs/2103.03208>`_.

.. code:: bash

	<ForceField name="example_custom" version="0.0.1">
	 <AtomTypes>
	  <Type name="C1" class="c3" element="C" mass="12.011" def="C(C)(H)(F)(F)" desc="carbon bonded to 2 Fs, a H, and another carbon"/>
	  <Type name="C2" class="c3" element="C" mass="12.011" def="C(C)(F)(F)(F)" desc="carbon bonded to 3 Fs and another carbon"/>
	  <Type name="F1" class="f" element="F" mass="18.998" def="FC(C)(F)H" desc="F bonded to C1"/>
	  <Type name="F2" class="f" element="F" mass="18.998" def="FC(C)(F)F" desc="F bonded to C2"/>
	  <Type name="H1" class="h2" element="H" mass="1.008" def="H(C)" desc="single H bonded to C1"/>
	 </AtomTypes>
	 <HarmonicBondForce>
	  <Bond class1="c3" class2="c3" length="0.15375" k="251793.12"/>
	  <Bond class1="c3" class2="f" length="0.13497" k="298653.92"/>
	  <Bond class1="c3" class2="h2" length="0.10961" k="277566.56"/>
	 </HarmonicBondForce>
	 <HarmonicAngleForce>
	  <Angle class1="c3" class2="c3" class3="f" angle="1.9065976748786053" k="553.1248"/>
	  <Angle class1="c3" class2="c3" class3="h2" angle="1.9237019015481498" k="386.6016"/>
	  <Angle class1="f" class2="c3" class3="f" angle="1.8737854849411122" k="593.2912"/>
	  <Angle class1="f" class2="c3" class3="h2" angle="1.898743693244631" k="427.6048"/>
	 </HarmonicAngleForce>
	 <PeriodicTorsionForce>
	  <Proper class1="f" class2="c3" class3="c3" class4="f" periodicity1="3" k1="0.0" phase1="0.0" periodicity2="1" k2="5.0208" phase2="3.141592653589793"/>
	  <Proper class1="" class2="c3" class3="c3" class4="" periodicity1="3" k1="0.6508444444444444" phase1="0.0"/>
	 </PeriodicTorsionForce>
	 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
	  <Atom type="C1" charge="0.224067"  sigma="0.371084" epsilon="0.304665"/>
	  <Atom type="C2" charge="0.500886"  sigma="0.393872" epsilon="0.222541"/>
	  <Atom type="F1" charge="-0.167131" sigma="0.298239" epsilon="0.208221"/>
	  <Atom type="F2" charge="-0.170758" sigma="0.276783" epsilon="0.237635"/>
	  <Atom type="H1" charge="0.121583"  sigma="0.264229" epsilon="0.071381"/>
	 </NonbondedForce>
	</ForceField>

The SMARTS strings are defined in the ``<AtomTypes>`` section.
The bond, angle, and dihedral parameters are defined in the following sections.
This time we load a ``.gro`` file for HFC-125
(:download:`here <../files/hfc125.gro>`),
load our custom force field XML file
(:download:`here <../files/ff_custom.xml>`),
and then apply the force field parameters.

.. code:: python

    import foyer
    import parmed

    mol = parmed.load("hfc125.gro")
    ff = foyer.ForceField("ff_custom.xml")

    mol_ff = ff.apply(mol)

	mol_ff.save("hfc125.top")

Foyer can be used to save input files for any simulation engine supported by
``parmed``. If you also install ``mbuild``, then a variety of other simulation
engines are also supported.
Installation
==============

For most users we recommend a conda installation:

.. code:: bash

    conda install -c conda-forge -c omnia foyer


If you wish to install from source, you can use the following commands:

.. code:: bash

    git clone https://github.com/mosdef-hub/foyer.git
    cd foyer
    conda env create -f environment.yml
    conda activate foyer
    pip install .

If you are using windows, you should use ``environment-win.yml`` rather than
``environment.yml``.


If you plan on contributing to the development of foyer, we recommend
you create an editable installation with all the required dependencies:

.. code:: bash

    git clone https://github.com/mosdef-hub/foyer.git
    cd foyer
    conda env create -f environment-dev.yml
    conda activate foyer-dev
    pip install -e .


Install pre-commit
------------------

We use [pre-commit](https://pre-commit.com/) to automatically handle our code formatting and this package is included in the dev environment.
With the ``foyer-dev`` conda environment active, pre-commit can be installed locally as a git hook by running::

    $ pre-commit install

And (optional) all files can be checked by running::

    $ pre-commit run --all-files


Supported Python Versions
-------------------------

Python 3.6 and 3.7 are officially supported, including testing during
development and packaging. Other Python versions, such as 3.8 and 3.5 and
older, may successfully build and function but no guarantee is made.

Testing your installation
-------------------------

foyer uses ``py.test`` for unit testing. To run them simply type run the
following while in the base directory::

    $ conda install pytest
    $ py.test -v

Building the documentation
--------------------------

foyer uses `sphinx <https://www.sphinx-doc.org/en/master/index.html>`_ to build its documentation. To build the docs locally, run the following while in the ``docs`` directory::

    $ conda env create -f docs-env.yml
    $ conda activate foyer-docs
    $ make html
Using foyer with Docker
========================

As much of scientific software development happens in unix platforms, to avoid the quirks of development dependent on system you use, a recommended way is to use docker or other containerization technologies. This section is a how to guide on using ``foyer`` with docker.

Prerequisites
-------------
A docker installation in your machine. Follow this `link <https://docs.docker.com/get-docker/>`_ to get a docker installation working on your machine. If you are not familiar with docker and want to get started with docker, the Internet is full of good tutorials like the ones `here <https://docker-curriculum.com/>`_ and `here <https://www.youtube.com/watch?v=zJ6WbK9zFpI&feature=youtu.be>`_.

Quick Start
-----------
After you have a working docker installation, please use the following command to use run a jupyter-notebook with all the dependencies for `foyer` installed:

.. code-block:: bash

    $ docker pull mosdef/foyer:latest
    $ docker run -it --name foyer -p 8888:8888 mosdef/foyer:latest


If no command is provided to the container (as above), the container starts a ``jupyter-notebook`` at the container location ``/home/anaconda/data``.
Then, the notebook can be accessed by copying and pasting the notebook URL into a web browser on your computer.
When finished with the session, you can use `Ctr`+`C` and follow instruction to exit the notebook as usual.
The docker container will exit upon notebook shutdown.

.. warning::

    Containers by nature are ephemeral, so filesystem changes (e.g., adding a new notebook) only persist until the end of the container's lifecyle.
    If the container is removed, any changes or code addition will not persist.
    See the section below for persistent data.

.. note::
    The ``-it`` flags connect your keyboard to the terminal running in the container.
    You may run the prior command without those flags, but be aware that the container will not respond to any keyboard input.
    In that case, you would need to use the ``docker ps`` and ``docker kill`` commands to shut down the container.


Persisting User Volumes
-----------------------
If you will be using `foyer` from a docker container, a recommended way is to mount what are called user volumes in the container. User volumes will provide a way to persist all filesystem/code additions made to a container regardless of the container lifecycle. For example, you might want to create a directory called `foyer-notebooks` in your local system, which will store all your `foyer` notebooks/code. In order to make that accessible to the container(where the notebooks will be created/edited), use the following steps:

.. code-block:: bash

    $ mkdir -p /path/to/foyer-notebooks
    $ cd /path/to/foyer-notebooks
    $ docker run -it --name foyer --mount type=bind,source=$(pwd),target=/home/anaconda/data -p 8888:8888 mosdef/foyer:latest

You can easily mount a different directory from your local machine by changing ``source=$(pwd)`` to ``source=/path/to/my/favorite/directory``.

.. note::

    The ``--mount`` flag mounts a volume into the docker container.
    Here we use a ``bind`` mount to bind the current directory on our local filesystem to the ``/home/anaconda/data`` location in the container.
    The files you see in the ``jupyter-notebook`` browser window are those that exist on your local machine.

.. warning::

    If you are using the container with jupyter notebooks you should use the ``/home/anaconda/data`` location as the mount point inside the container;
    this is the default notebook directory.

    Running Python scripts in the container
    ---------------------------------------
    Jupyter notebooks are a great way to explore new software and prototype code. However, when it comes time for production sciences, it is often better to work with python scripts.
    In order to execute a python script (``example.py``) that exists in the current working directory of your local machine, run:

.. code-block:: bash

    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/foyer:latest "python data/test.py"

Note that once again we are ``bind`` mounting the current working directory to ``/home/anaconda/data``.
The command we pass to the container is ``python data/test.py``.
Note the prefix ``data/`` to the script; this is because we enter the container in the home folder (``/home/anaconda``), but our script is located under ``/home/anaconda/data``.

.. warning::
    Do not bind mount to ``target=/home/anaconda``. This will cause errors.


If you don't want a Jupyter notebook, but just want a Python interpreter, you can run:

.. code-block:: bash
    $ docker run --mount type=bind,source=$(pwd),target=/home/anaconda/data mosdef/foyer:latest python

If you don't need access to any local data, you can of course drop the ``--mount`` command:

.. code-block:: bash
    $ docker run mosdef/foyer:latest python

Cleaning Up
-----------
You can remove the created container by using the following command:

.. code-block:: bash

    $ docker container rm foyer

.. note::

    Instead of using `latest`, you can use the image :code:`mosdef/foyer:stable` for most recent stable release of `foyer` and run the tutorials.
Examples from Foyer paper
~~~~~~~~~~~~~~~~~~~~~~~~~

Contained below are the toy examples from the *Usage Examples* section of the `foyer paper <https://arxiv.org/pdf/1812.06779.pdf>`__. The source code selections are listed below on this page, there are `Jupyter
Notebooks <https://github.com/mosdef-hub/foyer/tree/master/docs/examples>`__
where you can try these examples yourself. Note that these examples are
meant to showcase the abilities of ``foyer`` through simple examples. If
the user would like to examine more in-depth examples using ``foyer``
with ``mBuild``, please refer to the `tutorial
repository <https://github.com/mosdef-hub/mosdef_tutorials>`__.

Below is *Listing 6* from the paper, a python script to fill a :math:`2x2x2 nm`
box with 100 ethane molecules. The system is then atomtyped using the
OPLS-AA forcefield. There are two approaches to the same problem
detailed below in this listing, the first approach uses the
``forcefield_files`` function argument from
`mBuild <https://github.com/mosdef-hub/mbuild>`__ to atomptype the
system (using foyer under the hood). While the second approach creates a
``foyer`` ``Forcefield`` object, which then calls its ``apply``
function, operating on the ``mBuild`` ``Compound`` to return the
properly atomtyped structure. Note that in all instances when using
``foyer``, the chemical system of interest is converted into a
``ParmEd`` ``Structure``. Even the ``mBuild`` ``Compounds``, when
calling the ``save`` routine, are converted into a ``ParmEd``
``Structure`` before ``foyer`` can atomtype them. The object returned by
``foyer`` after the atomtypes have been applied are ``ParmEd``
``Structures``. This is subject to change in later iterations of
``foyer``.

Homogeneous fluid
^^^^^^^^^^^^^^^^^

.. code:: python

    import mbuild as mb
    from mbuild.lib.molecules import Ethane
    from foyer.examples.utils import example_file_path
    from foyer import Forcefield

    """ Applying a force field while saving from mBuild """
    # Create the chemical topology
    ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
    # Apply and save the topology
    ethane_fluid.save("ethane-box.top", forcefield_files=example_file_path("oplsaa_alkane.xml"))
    ethane_fluid.save("ethane-box.gro")

    """ Applying a force field directly with foyer """
    # Create the chemical topology
    ethane_fluid = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
    # Load the forcefield
    opls_alkane = Forcefield(forcefield_files=example_file_path("oplsaa_alkane.xml"))
    # Apply the forcefield to atom-type
    ethane_fluid = opls_alkane.apply(ethane_fluid)
    # Save the atom-typed system
    ethane_fluid.save("ethane-box.top", overwrite=True)
    ethane_fluid.save("ethane-box.gro", overwrite=True)


---------------------------------------

Fluid on silica substrate
^^^^^^^^^^^^^^^^^^^^^^^^^

The other example listing from the text showcases the ability to create
two separate chemical topologies and applying different forcefield files
to each. The two parameterized systems that are generated are then
combined into a single ``ParmEd`` ``Structure`` and saved to disk.

.. code:: python

    from foyer import Forcefield
    from foyer.examples.utils import example_file_path
    import mbuild as mb
    from mbuild.examples import Ethane
    from mbuild.lib.atoms import H
    from mbuild.lib.bulk_materials import AmorphousSilicaBulk
    from mbuild.lib.recipes import SilicaInterface
    from mbuild.lib.recipes import Monolayer

    # Create a silica substrate, capping surface oxygens with hydrogen
    silica=SilicaInterface(bulk_silica=AmorphousSilicaBulk())
    silica_substrate=Monolayer(surface=silica,chains=H(),guest_port_name="up")
    # Determine the box dimensions dictated by the silica substrate
    box=mb.Box(mins=[0, 0,max(silica.xyz[:,2])],maxs=silica.periodicity+ [0, 0, 4])
    # Fill the box with ethane
    ethane_fluid=mb.fill_box(compound=Ethane(),n_compounds=200,box=box)
    # Load the forcefields
    #opls_silica=Forcefield(forcefield_files=get_example_file("oplsaa_with_silica.xml"))
    opls_silica=Forcefield(forcefield_files=example_file_path("output.xml"))
    opls_alkane=Forcefield(forcefield_files=example_file_path("oplsaa_alkane.xml"))
    # Apply the forcefields
    silica_substrate=opls_silica.apply(silica_substrate)
    ethane_fluid=opls_alkane.apply(ethane_fluid)
    # Merge the two topologies
    system=silica_substrate+ethane_fluid
    # Save the atom-typed system
    system.save("ethane-silica.top")
    system.save("ethane-silica.gro")
SMARTS-based atomtyping
~~~~~~~~~~~~~~~~~~~~~~~

Foyer allows users to describe atomtypes using a modified version of
`SMARTS <http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`__
You may already be familiar with
`SMILES <https://www.wikiwand.com/en/Simplified_molecular-input_line-entry_system>`__
representations for describing chemical structures. SMARTS is a
straightforward extension of this notation.

Basic usage
^^^^^^^^^^^

Consider the following example defining the OPLS-AA atomtypes for a
methyl group carbon and its hydrogen atoms:

.. code:: xml

    <ForceField>
     <AtomTypes>
      <Type name="opls_135" class="CT" element="C" mass="12.01100" def="[C;X4](C)(H)(H)H" desc="alkane CH3"/>
      <Type name="opls_140" class="HC" element="H" mass="1.00800"  def="H[C;X4]" desc="alkane H"/>
     </AtomTypes>
    </ForceField>

This ``.xml`` format is an extension of the `OpenMM force field
format <http://docs.openmm.org/latest/userguide/application.html#creating-force-fields>`__
The above example utilizes two additional ``.xml`` attributes supported
by foyer: ``def`` and ``desc``. The atomtype that we are attempting to
match is always the **first** token in the SMARTS string, in the above
example, ``[C;X4]`` and ``H``. The ``opls_135`` (methyl group carbon) is
defined by a SMARTS string indicated a carbon with 4 bonds, a carbon
neighbor and 3 hydrogen neighbors. The ``opls_140`` (alkane hydrogen) is
defined simply as a hydrogen atom bonded to a carbon with 4 bonds.

Overriding atomtypes
^^^^^^^^^^^^^^^^^^^^

When multiple atomtype definitions can apply to a given atom, we must
establish precedence between those definitions. Many other atomtypers
determine rule precedence by placing more specific rules first and
evaluate those in sequence, breaking out of the loop as soon as a match
is found.

While this approach works, it becomes more challenging to maintain the
correct ordering of rules as the number of atomtypes grows. Foyer
iteratively runs all rules on all atoms and each atom maintains a
whitelist (rules that apply) and a blacklist (rules that have been
superceded by another rule). The set difference between the white- and
blacklists yields the correct atomtype if the force field is implemented
correctly.

We can add a rule to a blacklist using the ``overrides`` attribute in
the ``.xml`` file. To illustrate an example where overriding can be used
consider the following types describing alkenes and benzene:

.. code:: xml

    <ForceField>
     <AtomTypes>
      <Type name="opls_141" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)C" desc="alkene C (R2-C=)"/>
      <Type name="opls_142" class="CM" element="C" mass="12.01100" def="[C;X3](C)(C)H" desc="alkene C (RH-C=)"/>
      <Type name="opls_144" class="HC" element="H" mass="1.00800"  def="[H][C;X3]" desc="alkene H"/>
      <Type name="opls_145" class="CA" element="C" mass="12.01100" def="[C;X3;r6]1[C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6][C;X3;r6]1" overrides="opls_141,opls_142"/>
      <Type name="opls_146" class="HA" element="H" mass="1.00800"  def="[H][C;%opls_145]" overrides="opls_144" desc="benzene H"/>
     </AtomTypes>
    </ForceField>

If we’re atomtyping a benzene molecule, the carbon atoms will match the
SMARTS patterns for both ``opls_142`` and ``opls_145``. Without the
``overrides`` attribute, foyer will notify you that multiple atomtypes
were found for each carbon. Providing the ``overrides`` indicates that
if the ``opls_145`` pattern matches, it should supercede the specified
rules.

Supported SMARTS Grammar
^^^^^^^^^^^^^^^^^^^^^^^^^

We currently do not (yet) support all of `SMARTS’
features <http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`__.
`Here <https://github.com/mosdef-hub/foyer/issues/63>`__ we’re keeping
track of which portions are supported.
Usage Examples
~~~~~~~~~~~~~~

Foyer supports atomtyping of both all-atom and coarse-grained molecular
systems, and also allows for separate force fields to be used to
atom-type separate portions of a molecular system.

Creating a box of ethane
^^^^^^^^^^^^^^^^^^^^^^^^

Here we use mBuild to construct a box filled with ethane molecules and
use Foyer to atom-type the system, applying the OPLS force field, and
save to run-able GROMACS files.

.. code:: python

    import mbuild as mb
    from mbuild.lib.molecules import Ethane

    ethane_box = mb.fill_box(compound=Ethane(), n_compounds=100, box=[2, 2, 2])
    ethane_box.save('ethane-box.gro')
    ethane_box.save('ethane-box.top', forcefield_name='oplsaa')

--------------

Creating a box of coarse-grained ethane
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Again we will use mBuild to construct a box filled with ethane
molecules. However, now we will model ethane using a united-atom
description and apply the TraPPE force field during atom-typing. Note
how particle names are prefixed with an underscore so that Foyer knows
these particles are non-atomistic.

.. code:: python

    import mbuild as mb

    ethane_UA = mb.Compound()
    ch3_1 = mb.Particle(name='_CH3', pos=[0, 0, 0])
    ch3_2 = mb.Particle(name='_CH3', pos=[0.15, 0, 0])
    ethane_UA.add([ch3_1, ch3_2])
    ethane_UA.add_bond((ch3_1, ch3_2))

    ethane_UA_box = mb.fill_box(ethane_UA, 100, box=[2, 2, 2])
    ethane_UA_box.save('ethane-UA-box.gro')
    ethane_UA_box.save('ethane-UA-box.top', forcefield_name='trappe-ua')

--------------

Combining force fields
^^^^^^^^^^^^^^^^^^^^^^

In some instances, the use of multiple force fields may be desired to
describe a molecular system. For example, the user may want to use one
force field for a surface and another for a fluid in the same system.
Foyer supports this functionality by allowing the user to separately
atom-type parts of a system. In this example, we take a system featuring
bulk united atom ethane above a silica surface and apply the OPLS force
field to the surface and the TraPPE force field to the ethane. The two
atomtyped Parmed structures are then combined using a simple ‘+’
operator and can be saved to Gromacs files.

.. code:: python

    from foyer import Forcefield
    from foyer.examples.utils import example_file_path
    import mbuild as mb
    from mbuild.examples import Ethane
    from mbuild.lib.atoms import H
    from mbuild.lib.bulk_materials import AmorphousSilica
    from mbuild.lib.recipes import SilicaInterface
    from mbuild.lib.recipes import Monolayer

    interface = SilicaInterface(bulk_silica=AmorphousSilica())
    interface = Monolayer(surface=interface, chains=H(), guest_port_name='up')

    box = mb.Box(mins=[0, 0, max(interface.xyz[:,2])],
                 maxs=interface.periodicity + [0, 0, 4])

    ethane_box = mb.fill_box(compound=Ethane(), n_compounds=200, box=box)

    opls = Forcefield(name='oplsaa')
    opls_silica = Forcefield(forcefield_files=example_file_path('opls-silica.xml'))
    ethane_box = opls.apply(ethane_box)
    interface = opls_silica.apply(interface)

    system = interface + ethane_box

    system.save('ethane-silica.gro')
    system.save('ethane-silica.top')
Parameter definitions
=====================

Parameter definitions within force field XMLs follow the same
conventions as defined in the `OpenMM
documentation <http://docs.openmm.org/latest/userguide/application.html#creating-force-fields>`__.
Currently, only certain functional forms for molecular forces are
supported, while future developments are expected to allow Foyer to
support any desired functional form, including reactive and tabulated
potentials. The currently supported functional forms for molecular
forces are:

-  **Nonbonded**

   -  `Lennard-Jones
      (12-6) <http://docs.openmm.org/latest/userguide/application.html#nonbondedforce>`__

-  **Bonds**

   -  `Harmonic <http://docs.openmm.org/latest/userguide/application.html#harmonicbondforce>`__

-  **Angles**

   -  `Harmonic <http://docs.openmm.org/latest/userguide/application.html#harmonicangleforce>`__

-  **Torsions (proper)**

   -  `Periodic <http://docs.openmm.org/latest/userguide/application.html#periodictorsionforce>`__
   -  `Ryckaert-Bellemans <http://docs.openmm.org/latest/userguide/application.html#rbtorsionforce>`__

-  **Torsions (improper)**

   -  `Periodic <http://docs.openmm.org/latest/userguide/application.html#periodictorsionforce>`__

Definitions for each molecular force follow the `OpenMM standard <http://docs.openmm.org/latest/userguide/theory.html>`_.

The harmonic bond potential is defined as

.. math::

    E = \frac{1}{2}k(r-r_{0})^{2}

where `k` is the bond coefficient (:math:`\frac{energy}{distance^{2}}`) and `r`\ :sub:`0` is the equilibrium bond distance. Note the factor of :math:`\frac{1}{2}`.

Dihedral potentials reported as a fourier series (e.g., OPLS) can be converted to Ryckaert-Bellemans (RB) torsions as specified in the `GROMACS User Manual <https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#proper-dihedrals-ryckaert-bellemans-function>`_.

Classes vs. Types
-----------------

OpenMM allows users to specify either a ``class`` or a
``type`` (See `Atom Types and Atom Classes
<http://docs.openmm.org/latest/userguide/application.html#atom-types-and-atom-classes>`_),
to define each particle within the force definition. Here, ``type``
refers to a specific atom type (as defined in the ``<AtomTypes>``
section), while ``class`` refers to a more general description that can
apply to multiple atomtypes (i.e. multiple atomtypes may share the same
class). This aids in limiting the number of force definitions required
in a force field XML, as many similar atom types may share force
parameters.

Assigning parameters by specificity
-----------------------------------

Foyer deviates from OpenMM’s convention when matching force definitions
in a force field XML to instances of these forces in a molecular system.
In OpenMM, forces are assigned according to the first matching
definition in a force field XML, even if multiple matching definitions
exist. In contrast, Foyer assigns force parameters based on definition
specificity, where definitions containing more ``type`` attributes are
considered to be more specific.

**Example:**

.. code:: xml

   <RBTorsionForce>
     <Proper class1="CT" class2="CT" class3="CT" class4="CT" c0="2.9288" c1="-1.4644" c2="0.2092" c3="-1.6736" c4="0.0" c5="0.0"/>
     <Proper type1="opls_961" type2="opls_136" type3="opls_136" type4="opls_136" c0="-0.987424" c1="0.08363" c2="-0.08368" c3="-0.401664" c4="1.389088" c5="0.0"/>
   </RBTorsionForce>

Above, two proper torsions are defined, both describing a torsional
force between four tetrahedral carbons. However, the first definition
features four ``class`` attributes and zero ``type`` attributes, as this
describes a general dihedral for all tetrahedral carbons. The second
definition features zero ``class`` attributes and four ``type``
attributes, and describes a more specific dihedral for the case where
one end carbon is of type ``'opls_961'`` (a fluorinated carbon) and the
remaining three carbons are of type ``'opls_136'`` (alkane carbons). Now
consider we want to use a force field containing the above torsion
definitions to assign parameters to some molecular system that features
partially fluorinated alkanes. When assigning torsion parameters to a
quartet of atoms where one end carbon is fluorinated (``'opls_961'``)
and the remaining three are hydrogenated (``'opls_136'``), if using the
OpenMM procedure for parameter assignment the more general
``'CT-CT-CT-CT'`` torsion parameters (the first definition above) would
be assigned because this would be the first matching definition in the
force field XML. However, with Foyer, the second definition will be
auto-detected as more specific, due to the greater number of ``type``
attributes (4 vs. 0) and those parameters will be assigned instead.

It should be noted that if two definitions feature the same specificity
level (i.e. the same number of ``type`` definitions) then automatic
detection of precedence is not possible and parameter assignment will
follow the OpenMM procedure whereby parameters from the first matching
force definition in the XML will be assigned.
Applying a force field
======================

The main method you will use in **foyer** is the
``Forcefield.apply()`` method. There are a few important arguments
you should understand.

The first several are the ``assert_bond_params``
``assert_angle_params``, ``assert_dihedral_params``, and
``assert_improper_params``. These arguments require that the
supplied force field has parameters for every bond, angle, dihedral,
and improper in the system. In most cases, if you get an error,
it means that your force field is missing parameters for one of the
bonds/angles/dihedrals/impropers in the system. This could be
because the parameters are missing or because the atom-typing
(i.e., the SMARTS strings) are incorrect. These arguments are ``True`` by default,
with the exception of ``assert_improper_params``. In all cases, it
is wise to verify that the final files you generate have the expected
number of bonds/angles/dihedrals/impropers for your system.

The other important optional argument is the ``combining_rule`` option,
which is ``"lorentz"`` (Lorentz-Berthelot) by default. The other valid
option is ``"geometric"``, if your force field uses geometric combining rules.
============
Citing foyer
============

If you use foyer for your research, please cite `our paper <https://doi.org/10.1016/j.commatsci.2019.05.026>`_ or
its `pre-print <https://arxiv.org/abs/1812.06779>`_:

**ACS**

   Klein, C.; Summers, A. Z.; Thompson, M. W.; Gilmer, J. B.; Mccabe, C.; Cummings, P. T.; Sallai, J.; Iacovella, C. R. Formalizing Atom-Typing and the Dissemination of Force Fields with Foyer. Computational Materials Science 2019, 167, 215–227.

**BibTeX**

.. code-block:: bibtex

   @article{klein2019,
      title = "Formalizing atom-typing and the dissemination of force fields with foyer",
      journal = "Computational Materials Science",
      volume = "167",
      pages = "215 - 227",
      year = "2019",
      issn = "0927-0256",
      doi = "https://doi.org/10.1016/j.commatsci.2019.05.026",
      url = "http://www.sciencedirect.com/science/article/pii/S0927025619303040",
      author = "Christoph Klein and Andrew Z. Summers and Matthew W. Thompson and Justin B. Gilmer and Clare McCabe and Peter T. Cummings and Janos Sallai and Christopher R. Iacovella",
      keywords = "Molecular simulation, Force fields, Reproducibility, Open-source software",
   }

Download as :download:`BibTeX <../files/foyer_citation.bib>` or :download:`RIS <../files/foyer_citation.ris>`
=====
Units
=====

Foyer forcefield XML files use the following units. These follow
from the units used in
`OpenMM <http://docs.openmm.org/latest/userguide/theory.html#units>`_.

+----------+---------+
| Quantity |  Units  |
+==========+=========+
| distance |    nm   |
+----------+---------+
|   angle  | radians |
+----------+---------+
|   mass   |   amu   |
+----------+---------+
|  energy  |  kJ/mol |
+----------+---------+

Please note that the units of the parameters found in the parameterized
``parmed.Structure`` will follow the
`units used by ParmEd <https://parmed.github.io/ParmEd/html/dimensional_analysis.html#>`_
Forcefield Class
----------------

The primary data structure in foyer is the ``Forcefield`` class, which inherits
from the OpenMM class of the same name. The primary operation on this class is
the ``.apply()`` function, which takes a chemical topology and returns a
parametrized ``ParmEd`` ``Structure``. The user may pass some otions to this
function based on a particular use case.

.. autoclass:: foyer.forcefield.Forcefield
    :members:
Validation of force field files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Foyer performs several validation steps to help prevent malformed force
field files and SMARTS strings from making it into production code. Our
goal is to provide human readable error messages that users who may not
be intimately familiar with XML or our SMARTS parsing grammar can easily
act upon.

However, if you receive any unclear error messages or warnings we
strongly encourage you to `submit an issue <https://github.com/mosdef-hub/foyer/issues/new>`_
detailing the error message you received and, if possible, attach a minimal
example of the force field file that created the problem.

XML schema
^^^^^^^^^^

As a first line of defense, any force field files loaded by foyer is
validated by this `XML schema
definition <../foyer/forcefields/ff.xsd>`__. Here we enforce which
elements (e.g. ``HarmonicBondForce``) are valid and how their attributes
should be formatted. Additionally, the schema ensures that atomtypes are
not 1) defined more than once and that 2) atomtypes referenced in other
sections are actually defined in the ``<AtomTypes>`` element.

SMARTS validation
^^^^^^^^^^^^^^^^^

All SMARTS strings used to define atomtypes are parsed. Parsing errors
are captured and re-raised with error messages that allow you to pin
point the location of the problem in the XML file and within the SMARTS
string. Wherever possible, we attempt to provide helpful hints and we
welcome any contributions that help improve the clarity of our error
messages.

Additionally, we ensure that any atomtypes referenced using the
``%type`` or ``overrides`` `syntax <smarts.html>`__ are actually defined
in the ``<AtomTypes>`` element.
