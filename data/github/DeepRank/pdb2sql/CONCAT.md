# Change Log


All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## 0.5.1
- Updated `compute_CapriClass` conditions

## 0.5.0
- Added atom name selection for lrmsd calculation
- Removed hardcoded chainIDs and added chainID selection in `StructureSimilarity`

## 0.4.0
- Added `many2sql` to support read multiple PDB files
- Added support for `Path` objects of input PDB files
- Added support for `help(pdb2sql)`
- Updated assignment of chain IDs in `StructureSimilarity`

## 0.3.0
- Added `align` to superpose a structure to a specific axis or plane
- Added `superpose` to superpose two structures based on selections
- Added `fetch` to download PDB file with given PDB ID
- Updated `interface` object to take `pdb2sql` object as input# PDB2SQL

[![PyPI](https://img.shields.io/pypi/v/pdb2sql)](https://pypi.org/project/pdb2sql/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3232887.svg)](https://doi.org/10.5281/zenodo.3232887)
[![RSD](https://img.shields.io/badge/RSD-pdb2sql-red)](https://research-software.nl/software/pdb2sql)
![Build_Test](https://github.com/DeepRank/pdb2sql/workflows/Build_Test/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/DeepRank/pdb2sql/badge.svg)](https://coveralls.io/github/DeepRank/pdb2sql)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/36ad228df234488ab70ade6b2a80d54b)](https://www.codacy.com/gh/DeepRank/pdb2sql/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DeepRank/pdb2sql&amp;utm_campaign=Badge_Grade)
[![Documentation Status](https://readthedocs.org/projects/pdb2sql/badge/?version=latest)](https://pdb2sql.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02077/status.svg)](https://doi.org/10.21105/joss.02077)

`pdb2sql` is a Python package that leverage SQL queries to parse, manipulate and process PDB files. It provides:

-   a powerful `pdb2sql` object to convert PDB data in SQL database
-   strcuture transformation functions (rotations, translations...)
-   useful capablities to
    -   calculate structure interface (contact atoms and residues)
    -   calculate structure similarity (iRMSD, lRMSD, FNAT, DockQ...)

## Installation

```
pip install pdb2sql
```

## Documentation
The documentation of the package alongside small tutorial can be found at :
-  <https://pdb2sql.readthedocs.io>

## Quick Example

`pdb2sql` easily allows to load a PDB file in an object. Once loaded, the data can be parsed using SQL queries. To facilitate the adoption of the tool simple methods have been developped to wrap the SQL queries in simple methods. For example obtaining the positions of all carbon, nitrogen and oxygen atoms of chain A from all residues but VAL and LEU, one can use :

```python
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('x,y,z',
                name = ['C','N', 'O'],
                no_resName = ['VAL','LEU'],
                chainID = 'A')
```
---
title: 'The pdb2sql Python Package: Parsing, Manipulation and Analysis of PDB Files Using SQL Queries'
tags:
  - Python
  - Bioinformatics
  - PDB files
authors:
  - name: Nicolas Renaud
    orcid: 0000-0001-9589-2694
    affiliation: 1
  - name: Cunliang Geng
    orcid: 0000-0002-1409-8358
    affiliation: 1
affiliations:
 - name: Netherlands eScience Center, Science Park 140 1098 XG Amsterdam, the Netherlands
   index: 1
date: 13 January 2020
bibliography: paper.bibtex
---

# Summary

The analysis of biomolecular structures is a crucial task for a wide range of applications ranging from drug design to protein engineering. The Protein Data Bank (PDB) file format [@pdb] is the most popular format to describe biomolecular structures such as proteins and nucleic acids. In this text-based format, each line represents a given atom and entails its main properties such as atom name and identifier, residue name and identifier, chain identifier, coordinates, etc. Several solutions have been developed to parse PDB files into dedicated objects that facilitate the analysis and manipulation of biomolecular structures. This is, for example, the case for the ``BioPython``  parser [@biopython,@biopdb] that loads PDB files into a nested dictionary, the structure of which mimics the hierarchical nature of the biomolecular structure. Selecting a given sub-part of the biomolecule can then be done by going through the dictionary and selecting the required atoms. Other packages, such as ``ProDy`` [@prody], ``BioJava`` [@biojava], ``MMTK`` [@mmtk] and ``MDAnalysis`` [@mdanalysis] to cite a few, also offer solutions to parse PDB files. However, these parsers are embedded in large codebases that are sometimes difficult to integrate with new applications and are often geared toward the analysis of molecular dynamics simulations. Lightweight applications such as ``pdb-tools`` [@pdbtools] lack the capabilities to manipulate coordinates.



We present here the Python package ``pdb2sql``, which loads individual PDB files into a relational database. Among different solutions, the Structured Query Language (SQL) is a very popular solution to query a given database. However SQL queries are complex and domain scientists such as bioinformaticians are usually not familiar with them. This represents an important barrier to the adoption of SQL technology in bioinformatics. ``pdb2sql`` exposes complex SQL queries through simple Python methods that are intuitive for end users.  As such, our package leverages the power of SQL queries and removes the barrier that SQL complexity represents. In addition, several advanced modules have also been built, for example, to rotate or translate biomolecular structures, to characterize interface contacts, and to measure structure similarity between two protein complexes. Additional modules can easily be developed following the same scheme. As a consequence, ``pdb2sql`` is a lightweight and versatile PDB tool that is easy to extend and to integrate with new applications.


# Capabilities of ``pdb2sql``

``pdb2sql`` allows a user to query, manipulate, and process PDB files through a series of dedicated classes. We give an overview of these features and illustrate them with snippets of code. More examples can be found in the documentation (https://pdb2sql.readthedocs.io).

## Extracting data from PDB files

``pdb2sql`` allows a user to simply query the database using the ``get(attr, **kwargs)`` method. The attribute ``attr`` here is a list of or a single column name of the ``SQL`` database; see Table 1 for available attributes. The keyword argument ``kwargs`` can then be used to specify a sub-selection of atoms.

Table 1. Atom attributes and associated definitions in ``pdb2sql``

| Attribute | Definition                                |
|-----------|-------------------------------------------|
| serial    | Atom serial number                        |
| name      | Atom name                                 |
| altLoc    | Alternate location indicator              |
| resName   | Residue name                              |
| chainID   | Chain identifier                          |
| resSeq    | Residue sequence number                   |
| iCode     | Code for insertion of residues            |
| x         | Orthogonal coordinates for X in Angstroms |
| y         | Orthogonal coordinates for Y in Angstroms |
| z         | Orthogonal coordinates for Z in Angstroms |
| occ       | Occupancy                                 |
| temp      | Temperature factor                        |
| element   | Element symbol                            |
| model     | Model serial number                       |


Every attribute name can be used to select specific atoms and multiple conditions can be easily combined. For example, let's consider the following example:

```python
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('x,y,z',
                name=['C','H'],
                resName=['VAL','LEU'],
                chainID='A')
```

This snippet extracts the coordinates of the carbon and hydrogen atoms that belong to all the valine and leucine residues of the chain labelled `A` in the PDB file. Atoms can also be excluded from the selection by appending the prefix ``no_`` to the attribute name. This is the case in the following example:

```python
from pdb2sql import pdb2sql
pdb = pdb2sql('1AK4.pdb')
atoms = pdb.get('name, resName',
                no_resName=['GLY', 'PHE'])
```
This snippet extracts the atom and residue names of all atoms except those belonging to the glycine and phenylalanine residues of the structure. Similar combinations of arguments can be designed to obtain complex selection rules that precisely select the desired atom properties.

## Manipulating PDB files

The data contained in the SQL database can also be modified using the ``update(attr, vals, **kwargs)`` method. The attributes and keyword arguments are identical to those in the ``get`` method. The ``vals`` argument should contain a `numpy` array whose dimension should match the selection criteria.  For example:

```python
import numpy as np
from pdb2sql import pdb2sql

pdb = pdb2sql('1AK4.pdb')
xyz = pdb.get('x,y,z', chainID='A', resSeq=1)
xyz = np.array(xyz)
xyz -= np.mean(xyz)
pdb.update('x,y,z', xyz, chainID='A', resSeq=1)
```

This snippet first extracts the coordinates of atoms in the first residue of chain A, then translates this fragment to the origin and updates the coordinate values in the database. ``pdb2sql`` also provides a convenient class ``transform`` to easily translate or rotate structures. For example, to translate the first residue of the structure 5 Å along the Y-axis,

```python
import numpy as np
from pdb2sql import pdb2sql
from pdb2sql import transform

pdb = pdb2sql('1AK4.pdb')
trans_vec = np.array([0,5,0])
transform.translation(pdb, trans_vec, resSeq=1, chainID='A')
```

One can also rotate a given selection around a given axis with the `rotate_axis` method:

```python
angle = np.pi
axis = (1., 0., 0.)
transform.rot_axis(pdb, axis, angle, resSeq=1, chainID='A')
```

## Identifying interface

The ``interface`` class is derived from the ``pdb2sql`` class and offers functionality to identify contact atoms or residues between two different chains with a given contact distance. It is useful for extracting and analysing the interface of, e.g., protein-protein complexes. The following example snippet returns all the atoms and all the residues of the interface of '1AK4.pdb' defined by a contact distance of 6 Å.

```python
from pdb2sql import interface

pdb = interface('1AK4.pdb')
atoms = pdb.get_contact_atoms(cutoff=6.0)
res = pdb.get_contact_residues(cutoff=6.0)
```

It is also possible to directly create an ``interface`` instance with a ``pdb2sql`` instance as input. In this case, all the changes in the ``pdb2sql`` instance before creating the new ``interface`` instance will be kept in the ``interface`` instance; afterwards, the two instances will be independent, which means changes in one will not affect the other.

```python
from pdb2sql import pdb2sql
from pdb2sql import interface

pdb = pdb2sql('1AK4.pdb')
pdbitf = interface(pdb)
atoms = pdbitf.get_contact_atoms(cutoff=6.0)
res = pdbitf.get_contact_residues(cutoff=6.0)
```


## Computing Structure Similarity

The ``StructureSimilarity`` class allows a user to compute similarity measures between two protein-protein complexes. Several popular measures used to classify qualities of protein complex structures in the CAPRI (Critical Assessment of PRedicted Interactions) challenges [@capri] have been implemented: interface rmsd, ligand rmsd,  fraction of native contacts and DockQ [@dockq]. The approach implemented to compute the interface rmsd and ligand rmsd is identical to the well-known package ``ProFit`` [@profit]. All the methods required to superimpose structures have been implemented in the ``transform`` class and therefore this relies on no external dependencies. The following snippet shows how to compute these measures:

```python
from pdb2sql import StructureSimilarity

sim = StructureSimilarity(decoy = '1AK4_model.pdb',
                          ref = '1AK4_xray.pdb')

irmsd = sim.compute_irmsd_fast()
lrmsd = sim.compute_lrmsd_fast()
fnat = sim.compute_fnat_fast()
dockQ = sim.compute_DockQScore(fnat, lrmsd, irmsd)
```


# Application
``psb2sql`` has been used at the Netherlands eScience center for bioinformatics projects. This is, for example, the case of ``iScore`` [@iscore], which uses graph kernels and support vector machines to rank protein-protein interfaces. We illustrate the use of the package here by computing the interface rmsd and ligand rmsd of a series of structural models using the experimental structure as a reference. This is a common task for protein-protein docking, where a large number of docked conformations are generated and have then to be compared to ground truth to identify the best-generated poses. This calculation is usually done using the ProFit software and we, therefore, compare our results with those obtained with ProFit. The code to compute the similarity measure for different decoys is simple:

```python
from pdb2sql import StructureSimilarity

ref = '1AK4.pdb'
decoys = os.listdir('./decoys')
irmsd = {}

for d in decoys:g
    sim = StructureSimilarity(d, ref)
    irmsd[d] = sim.compute_irmsd_fast(method='svd', izone='1AK4.izone')
```

Note that the method will compute the i-zone, i.e., the zone of the proteins that form the interface in a similar way to ProFit. This is done for the first calculations and the i-zone is then reused for the subsequent calculations. The comparison of our interface rmsd values to those given by ProFit is shown in Fig 1.

![Example figure.](sim.png)
Figure 1. Left - Superimposed model (green) and reference (cyan) structures. Right - comparison of interface rmsd values given by `pdb2sql` and by `ProFit`.

# Acknowledgements
We acknowledge contributions from Li Xue, Sonja Georgievska, and Lars Ridder.


# References
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Environment:**
- OS system:
- Version:
- Branch commit ID: 
- Inputs:

**To Reproduce**
Steps/commands to reproduce the behaviour:
 1.  
 2.
 3. 

**Expected Results**
A clear and concise description of what you expected to happen.

**Actual Results or Error Info**
If applicable, add screenshots to help explain your problem.

**Additional Context**
Add any other context about the problem here.
############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

#. you have a question;
#. you think you may have found a bug (including unexpected behavior);
#. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

#. use the search functionality `here <https://github.com/DeepRank/pdb2sql/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here </https://github.com/DeepRank/pdb2sql/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
#. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

#. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
#. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
#. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
#. make sure the existing tests still work by running ``pytest`` from the `test/` folder;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the PDB2SQL repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
###############################################################################
Contributor Covenant Code of Conduct
###############################################################################

Our Pledge
**********

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

Our Standards
*************

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

Our Responsibilities
********************

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

Scope
*****

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

Enforcement
***********

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at n.renaud@esciencecenter.nl. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

Attribution
***********

This Code of Conduct is adapted from the `Contributor Covenant <https://www.contributor-covenant.org>`_, version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
\.setup\_save module
====================

.. automodule:: .setup_save
   :members:
   :undoc-members:
   :show-inheritance:
=========
Utilities
=========
.. currentmodule:: pdb2sql.utils

PDB Tools
~~~~~~~~~
.. autosummary::
    :toctree: api/

    fetch

..
    .. automodule:: pdb2sql.transform
    :members:
    :undoc-members:
    :show-inheritance:
=============================
Superposition and alignement
=============================

Superposition
~~~~~~~~~~~~~~~

.. currentmodule:: pdb2sql.superpose

.. autosummary::
    :toctree: api/

    superpose
    superpose_selection    

..
    .. automodule:: pdb2sql.superpose
    :members:
    :undoc-members:
    :show-inheritance:


Alignement
~~~~~~~~~~~

.. currentmodule:: pdb2sql.align

.. autosummary::
    :toctree: api/

    align
    align_interface

..
    .. automodule:: pdb2sql.align
    :members:
    :undoc-members:
    :show-inheritance:
=========
Interface
=========
.. autoclass:: pdb2sql.interface.interface
.. currentmodule:: pdb2sql.interface.interface

Contact Atoms
~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    get_contact_atoms

Contact Residues
~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    get_contact_residues=======
PDB2SQL
=======

This module is based on :py:mod:`sqlite3`.

.. autoclass:: pdb2sql.pdb2sqlcore.pdb2sql
.. currentmodule:: pdb2sql.pdb2sqlcore.pdb2sql

Process PDB
~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    read_pdb

.. currentmodule:: pdb2sql.pdb2sql_base.pdb2sql_base
.. autosummary::
    :toctree: api/

    exportpdb
    sql2pdb

.. currentmodule:: pdb2sql.pdb2sqlcore.pdb2sql

Get SQL Data
~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    get
    get_colnames
    get_chains
    get_residues
    get_xyz

Set SQL data
~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    update
    add_column
    update_column
    update_xyz

Print SQL data
~~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    print
    print_colnames

..
    .. automodule:: pdb2sql.pdb2sqlcore
    :members:
    :undoc-members:
    :show-inheritance:
====================
Structure Similarity
====================
.. autoclass:: pdb2sql.StructureSimilarity.StructureSimilarity

.. currentmodule:: pdb2sql.StructureSimilarity.StructureSimilarity

i-RMSD (interface RMSD)
~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    compute_irmsd_fast
    compute_irmsd_pdb2sql
    compute_izone

l-RMSD (ligand RMSD)
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    compute_lrmsd_fast
    compute_lrmsd_pdb2sql
    compute_lzone

FNAT (Fraction of native contacts)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    compute_fnat_fast
    compute_fnat_pdb2sql
    compute_residue_pairs_ref

DockQ
~~~~~~
.. autosummary::
    :toctree: api/

    compute_DockQScore

CAPRI classes
~~~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    compute_CapriClass

Clashes
~~~~~~~
.. autosummary::
    :toctree: api/

    compute_clashes

..
    .. automodule:: pdb2sql.StructureSimilarity
        :members:
        :undoc-members:
        :show-inheritance:
=========
Transform
=========
.. currentmodule:: pdb2sql.transform

Rotation
~~~~~~~~
.. autosummary::
    :toctree: api/

    rotate
    rot_mat

    get_rot_axis_angle

    rot_axis
    rot_xyz_around_axis

    rot_euler
    rotation_euler


Translation
~~~~~~~~~~~
.. autosummary::
    :toctree: api/

    translation


..
    .. automodule:: pdb2sql.transform
    :members:
    :undoc-members:
    :show-inheritance:
.. ipython:: python
    :suppress:

    # change working dir to docs/
    import os
    os.chdir('..')

=====================
10 minutes to pdb2sql
=====================

This is a short introduction to pdb2sql.


Download PDB files
------------------

A handy tool `fetch` is provided to download PDB files from `PDB <https://www.rcsb.org>`_ website.

.. ipython:: python

    from pdb2sql import fetch
    fetch('3CRO', './pdb/')
    ls ./pdb

For clear illustration, some `dummy PDB files <https://github.com/DeepRank/pdb2sql/tree/master/docs/pdb>`_
are used in the following examples.

Get and set data
----------------

First, we import as follows:

.. ipython:: python

    from pdb2sql import pdb2sql

Create a SQL database instance:

.. ipython:: python

    db = pdb2sql("./pdb/dummy.pdb")


The ``db`` is a SQL instance that contains one table named *ATOM*.

In this table, each row represents one atom, and columns are atom properties:

.. ipython:: python

    db.print()

Get data
^^^^^^^^

Get chainID, residue number, residue name and atom name of all atoms:

.. ipython:: python

    p = db.get('chainID, resSeq, resName, name')
    p

Get x,y,z coordinates of all atoms:

.. ipython:: python

    p = db.get('x,y,z')
    p

Get x,y,z coordinates of chain A atoms:

.. ipython:: python

    p = db.get('chainID, x,y,z', chainID=['A'])
    p

Get x,y,z coordinates of atoms on residue 1 and 4 of Chain A

.. ipython:: python

    p = db.get('chainID,resSeq,x,y,z', chainID=['A'], resSeq=['1', '4'])
    p

Get data of all atoms except residue MET and GLN atoms

.. ipython:: python

    p = db.get('chainID, resSeq, resName, name', no_resName = ['MET', 'GLN'])
    p

Get data of all atoms except residue MET and GLN atoms or CA (carbon alpha) atoms

.. ipython:: python

    p = db.get('chainID, resSeq, resName, name', no_resName = ['MET', 'GLN'], no_name = ['CA'])
    p


Get all data, a simple way is ``db.get('*')``.

A shortcut to get x,y,z coordinates:

.. ipython:: python

    p = db.get_xyz()
    p

Get chain IDs:

.. ipython:: python

    p = db.get_chains()
    p

Get residue list:

.. ipython:: python

    p = db.get_residues()
    p


Filter the data base
^^^^^^^^^^^^^^^^^^^^^^^^^^

pdb2sql allows to create a new database by filtering the one we jut created

.. ipython:: python

    db_chainA = db(chainID='A')
    db_chainA.print()

In that example `dp_chainA` is a sql database that only includes the atoms from chain A.
All the selection keywords (chainID, resSeq, resName, name) and their negations 
(no_chainID, no_resSeq, no_resName, no_name) can be used and combined to obtain the new database.

Set data
^^^^^^^^

Rename chain B to C:

.. ipython:: python

    num_B_atoms = len(db.get('chainID', chainID=['B']))
    chainC = ['C'] * num_B_atoms
    db.get_chains()
    db.update('chainID', chainC, chainID = ['B'])
    db.get_chains()


Update x,y,z coordinates for structure translatation of [10,10,10]

.. ipython:: python

    xyz_old = db.get_xyz()
    xyz = np.array(xyz_old) + 10.0
    db.update('x,y,z', xyz)
    xyz_new = db.get_xyz()
    print("old:\n", xyz_old)
    print("new:\n", xyz_new)

Update a column using index, e.g. change the x coordinates of the first
10 atoms to 2:

.. ipython:: python

    x = np.ones(10) + 1
    db.update_column('x', values=x, index=list(range(10)))
    db.print('serial, name, x')

Add a new column *type* with value *high*:

.. ipython:: python

    db.add_column('type', value = 'high', coltype = 'str')
    db.print('serial, name, type')


PDB I/O
-------

Read PDB file or data to a list:

.. ipython:: python

    pdb = pdb2sql.read_pdb('./pdb/dummy.pdb')
    pdb

Convert SQL data to PDB-formated data:

.. ipython:: python

    pdb = db.sql2pdb()
    pdb

Write PDB file from SQL database:

.. ipython:: python

    db.exportpdb('./pdb/test.pdb')

    # show the test.pdb file
    ls ./pdb



Interface calculation
---------------------

Create an :class:`~pdb2sql.interface.interface` SQL database instance:

.. ipython:: python

    from pdb2sql import interface

    # use pdb2sql instance as input
    from pdb2sql import pdb2sql
    pdb_db = pdb2sql('./pdb/3CRO.pdb')
    db = interface(pdb_db)

    # or use pdb file as input
    db = interface('./pdb/3CRO.pdb')

Interface atoms
^^^^^^^^^^^^^^^

.. ipython:: python

    itf_atom = db.get_contact_atoms(cutoff = 3)
    itf_atom_pair = db.get_contact_atoms(cutoff = 3, return_contact_pairs=True)
    print("interface atom:\n", itf_atom)
    print("interface atom pairs:\n", itf_atom_pair)


Interface residues
^^^^^^^^^^^^^^^^^^

.. ipython:: python

    itf_residue = db.get_contact_residues(cutoff = 3)
    itf_residue_pair = db.get_contact_residues(cutoff = 3, return_contact_pairs=True)
    itf_residue
    itf_residue_pair


Structure superposition  
--------------------------

pdb2sql allows to superpose two structure on top of each other either using the full structure or with selection keywords.
For example to superpose the chain A of two PDB one can use :

.. ipython:: python 

    from pdb2sql import superpose
    ref = pdb2sql('./pdb/1AK4_5w.pdb')
    decoy = pdb2sql('./pdb/1AK4_10w.pdb')
    superposed_decoy = superpose(decoy, ref, chainID='A', export=True)

This will export a new PDB file containining the structure of the decoy superposed onto the reference.

Structure alignement
---------------------------

pdb2sql allows to align structure along a specific axis

.. ipython:: python

    from pdb2sql import align
    db = pdb2sql('./pdb/1AK4_10w.pdb')
    aligned_db = align(db, axis='z', export=True)

The alignement can  also consider only a subpart of the complex using the selection keywords:

.. ipython:: python

    aligned_db = align(db, axis='z', chainID='A')

There the chain A will be aligned along the z-axis

This will create a new PDB file containing the structure aligned along the z-axis. It is 
also possible aligning an interface in a given plane

.. ipython:: python

    from pdb2sql import align_interface
    db = pdb2sql('./pdb/3CRO.pdb')
    aligned_db = align_interface(db, plane='xy', export=True)

By default the interface formed by chain A and B will be considered. In case multiple chains are present
in the structure it is possible to specify wich interface to consider:

.. ipython:: python

    aligned_db = align_interface(db, plane='xy', chain1='L', chain2='R')


There the interface between chain L and R will be considered. Note that any other selection
keyword can be used to specify which interface to account for.

Structure similarity calculation
--------------------------------

Create a :class:`~pdb2sql.StructureSimilarity.StructureSimilarity` instance:

.. ipython:: python

    from pdb2sql.StructureSimilarity import StructureSimilarity
    sim = StructureSimilarity('./pdb/decoy.pdb', './pdb/ref.pdb')

interface RMSD
^^^^^^^^^^^^^^

.. ipython:: python
    :okwarning:

    irmsd_fast = sim.compute_irmsd_fast()
    irmsd_pdb2sql = sim.compute_irmsd_pdb2sql()
    irmsd_fast
    irmsd_pdb2sql


ligand RMSD
^^^^^^^^^^^

.. ipython:: python
    :okwarning:

    lrmsd_fast = sim.compute_lrmsd_fast()
    lrmsd_pdb2sql = sim.compute_lrmsd_pdb2sql()
    lrmsd_fast
    lrmsd_pdb2sql

FNAT
^^^^

Calculate the fraction of native contacts:

.. ipython:: python
    :okwarning:

    fnat_fast = sim.compute_fnat_fast()
    fnat_pdb2sql = sim.compute_fnat_pdb2sql()
    fnat_fast
    fnat_pdb2sql


DockQ score
^^^^^^^^^^^

.. ipython:: python

    dockQ = sim.compute_DockQScore(fnat_fast, lrmsd_fast, irmsd_fast)
    dockQ


Structure transformation
------------------------

Create SQL instance:

.. ipython:: python

    from pdb2sql import transform
    db = pdb2sql('./pdb/dummy_transform.pdb')

The atom coordinates are:

.. ipython:: python

    db.get_xyz()

Rotations
^^^^^^^^^
Rotate structures 180 degrees along the x-axis:

.. ipython:: python

    angle = np.pi
    axis = (1., 0., 0.)
    transform.rot_axis(db, axis, angle)
    db.get_xyz()

Get random rotation axis and angle:

.. ipython:: python

    axis, angle = transform.get_rot_axis_angle()
    axis
    angle

Translations
^^^^^^^^^^^^

Translate structure 5Å along y-axis:

.. ipython:: python

        trans_vec = np.array([0,5,0])
        transform.translation(db, trans_vec)
        db.get_xyz()*************************************
pdb2sql: Processing PDB data with SQL
*************************************

`pdb2sql`_ is a Python package that allows to use SQL queries to handle `PDB`_ files.
Currently, only 'ATOM' data is parsed, and other items of PDB, e.g. HETATM, are ignored.

Installation:
    ``pip install pdb2sql``

.. _pdb2sql: https://github.com/DeepRank/pdb2sql
.. _PDB: https://www.rcsb.org/

========
Tutorial
========

.. toctree::
    :maxdepth: 1

    10 minutes to pdb2sql <tutorial.rst>

==========
Python API
==========

.. toctree::
    :maxdepth: 1

    PDB2SQL <pdb2sql.pdb2sqlcore>
    Interface <pdb2sql.interface>
    Superposition <pdb2sql.superpose>
    Structure Similarity <pdb2sql.StructureSimilarity>
    Structure Transformation <pdb2sql.transform>
    Utilities <pdb2sql.utils>

.. :caption: Python API
..
    ==================
    Indices and tables
    ==================

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`
{{ name | escape | underline }}

.. currentmodule:: {{ module }}

.. automethod:: {{ fullname }}