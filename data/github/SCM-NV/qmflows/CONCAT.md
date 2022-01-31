# Version 0.11.2 (unreleased)

## New
 * *placeholder*.


# Version 0.11.1 (21/12/2022)

## New
 * Introduced a new template for frequency analyses with CP2K (#278).
 * Allow ``dir()`` to work on result-based generic properties.
 * Added the ``basis`` and ``potential`` generic keywords to CP2K.


# Version 0.11.0 (17/11/2021)

## New
 * Add support for reading CP2K MOs from unrestricted calculations.
 * Add support for reading CP2K >=8.2 MOs.
 * Add a template for (CP2K) cell optimizations: ``qmflows.cell_opt``.
 * Add a generic keyword for the CP2K GAL19 non-bonded forcefield.
 * Add 6 new generic properties to ``qmflows.cp2k`` and ``qmflows.cp2k_mm`` outputs:
   * ``volume``
   * ``forces``
   * ``coordinates``
   * ``temperature``
   * ``lattice``
   * ``pressure``

## Changed
 * Make ``qmflows.Package`` instance more compatible with builtin functions.
 * Remove the unused ``__block_replace`` functionality.
 * Remove the cell parameters from the ``qmflows.cp2k_mm`` templates.
 * Remove the 2-digit restriction from CP2K cell parameters.
 * Check for duplicate keys when parsing .yaml inputs.
 * QMFlows templates are now always copied when getting them (requires Python >= 3.7).
 * Make RDKit an optional dependency (requires Python >= 3.7).

## Fix
 * Fix the ``ResultWrapper`` parameters being ordered incorrectly.
 * Fix ``qmflows.cp2m_mm`` ignoring the ``executable`` key.
 * Fix ``qmflows.InitRestart`` failing on consecutive calls.
 * Fix ``qmflows.CP2KMM_Result`` not inheriting from ``qmflows.CP2K_Result``.
 * Remove usage of the CP2K ``USE_ELEMENT_AS_KIND`` keyword.


# Version 0.10.4 (07/09/2020)
## New
 * Introduced a flag for keeping the Log files

## Fix
 * Improve CP2K error reporting (#209)


# Version 0.10.3 (12/06/2020)

## New
  * Added tests for generating the Sphinx documentation.

## Changed
  * Replaced ``requirements.txt`` with ``.readthedocs.yml``.
  * Fixed the jupyter notebook in the documentation.


# Version 0.10.2 (12/06/2020)

## New
  * Allow other cp2k executable: ``cp2k.sopt``, ``cp2k.psmp``, etc.


# Version 0.10.1 (09/06/2020)

## Changed
  * Exposed ``InitRestart`` to the main QMFlows ``__init__.py`` file.
  * Exchanged ``plams.init()`` / ``plams.finish()`` for ``qmflows.InitRestart`` in the ``qmflows.run()`` function.
  * Store the ``cache.db`` file in the PLAMS working directory.


# Version 0.10.0 (XX/03/2020)

## Added
  * Introduced the ``CP2KMM`` class for classical forcefield calculations with CP2K: [qmflows/pull/150](https://github.com/SCM-NV/qmflows/pull/150).
  * Introduced the ``PackageWrapper`` class: [qmflows/pull/149](https://github.com/SCM-NV/qmflows/pull/149).
  * Introduced updates and code-style improvements to the ``Package`` and ``Result`` classes: [qmflows/pull/146](https://github.com/SCM-NV/qmflows/pull/146)
  * Added workflow for [GitHub Actions](https://github.com/SCM-NV/qmflows/actions)

## Removed
  * [Removed references to Dirac](https://github.com/SCM-NV/qmflows/issues/152)
  * [Removed Pymonad](https://github.com/SCM-NV/qmflows/issues/156)
  * Remove support for [FDE](https://github.com/SCM-NV/qmflows/issues/171)

# Changed
  * Used [Path](https://github.com/SCM-NV/qmflows/issues/153) instead of ``str``.


# Version 0.9.0 (27/11/2019)

## Changed
  * Use autopep to format the code

## Removed
  * Interface to HDF5
  * Turbomol Parser
  * graphviz dependency


# Version 0.8.0 (17/06/2019)

## Changed

 * Used [pyyaml](https://pyyaml.org/wiki/PyYAMLDocumentation) for the [templates](https://github.com/SCM-NV/qmflows/blob/master/src/qmflows/templates/templates.py) instead of *JSON*
 * Updated documentation
 * Test wiht python 3.7


# Version 0.4.0 (25/02/2019)

## Changed

  * Moved all the functionality to build and analysis quantum dot structures to [their own repo](https://github.com/BvB93/CAT)
  * Moved `molkit` functionality to [PLAMS](https://github.com/SCM-NV/PLAMS)

# 08/01/2019

## Removed
*  Quantum Dots builder functionality moved to [CAT](https://github.com/BvB93/CAT)



# 19/10/2018

## Added
 * Quantum Dots builder functionality

## Changed

 * Use [noodles==0.3.0](https://github.com/NLeSC/noodles/releases)
 * Replace [nose](https://nose.readthedocs.io/en/latest/) with [pytest](https://docs.pytest.org/en/latest/)
 * Imported only the core functionality of the library
 * Used [intra-package-references](https://docs.python.org/3/tutorial/modules.html#intra-package-references)
 * Used `__all__` to limit exposed fuctionality

## Removed

 * Dead code related to all noodles API
 * All the `import *`
 * Dead code from components including PES

## Fixed

 * Job manager issue when removing a SCM job



# 02/12/2018

## Added
 * Ligand MOPAC+COSMO-RS property calculation
 * Inter-ligand activation strain analysis have been (UFF)

## Changed

 * Ligand optimization has been overhauled
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


**Additional context**
Add any other context about the problem here.
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

#. use the search functionality `here <https://github.com/SCM-NV/qmflows/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/SCM-NV/qmflows/issues>`__ to see if someone already filed the same issue;
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
#. make sure the existing tests still work by running ``pytest test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the qmflows repository on GitHub;
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

.. image:: https://github.com/SCM-NV/qmflows/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/SCM-NV/qmflows/actions
.. image:: https://codecov.io/gh/SCM-NV/qmflows/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/SCM-NV/qmflows
.. image:: https://readthedocs.org/projects/qmflows/badge/?version=latest
   :target: https://qmflows.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3274284.svg
   :target: https://doi.org/10.5281/zenodo.3274284
.. image:: https://badge.fury.io/py/qmflows.svg
   :target: https://badge.fury.io/py/qmflows
.. image:: qmflows.png

QMFlows
#######
See documentation_ for tutorials and documentation.

Motivation
==========
Research on modern computational quantum chemistry relies on a set of computational
tools to carry out calculations. The complexity of the calculations usually requires
intercommunication between the aforementioned tools, such communication is usually done
through shell scripts that try to automate input/output actions like: launching
the computations in a cluster, reading the resulting output and feeding the relevant
numerical result to another program. Such scripts are difficult to maintain and extend,
requiring a significant programming expertise to work with them. Being then desirable a
set of automatic and extensible tools that allows to perform complex simulations in
heterogeneous hardware platforms.

This library tackles the construction and efficient execution of computational chemistry workflows.
This allows computational chemists to use the emerging massively parallel compute environments in
an easy manner and focus on interpretation of scientific data rather than on tedious job submission
procedures and manual data processing.

Description
===========
This library consists of a set of modules written in Python3 to
automate the following tasks:

 1. Input generation.
 2. Handle tasks dependencies (Noodles_).
 3. Advanced molecular manipulation capabilities with (rdkit_).
 4. Jobs failure detection and recovery.
 5. Numerical data storage (h5py_).

Tutorial and Examples
---------------------
A tutorial written as a jupyter-notebook_ is available from: tutorial-qmflows_. You can
also access direclty more advanced examples_.

Installation
============

- Download miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Install according to: installConda_.

- Create a new virtual environment using the following commands:

  - ``conda create -n qmflows``

- Activate the new virtual environment

  - ``source activate qmflows``

To exit the virtual environment type  ``source deactivate``.


.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate qmflows``

Using the conda environment the following packages should be installed:


- install rdkit_ and h5py_ using conda:

  - ``conda install -y -q -c conda-forge rdkit h5py``

  - Note that ``rdkit`` is optional for Python 3.7 and later.

.. _installation:

Package installation
--------------------
Finally install the package:

- Install **QMFlows** using pip:
  - ``pip install qmflows``

Now you are ready to use *qmflows*.


  **Notes:**

  - Once the libraries and the virtual environment are installed, you only need to type
    ``conda activate qmflows`` each time that you want to use the software.


.. _documentation: https://qmflows.readthedocs.io/en/latest/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
.. _anaconda: https://www.anaconda.com/distribution/#download-section
.. _installConda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _Noodles: http://nlesc.github.io/noodles/
.. _h5py: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _jupyter-notebook: http://jupyter.org/
.. _tutorial-qmflows: https://github.com/SCM-NV/qmflows/tree/master/jupyterNotebooks
.. _examples: https://github.com/SCM-NV/qmflows/tree/master/src/qmflows/examples
.. _PLAMS: https://github.com/SCM-NV/PLAMS
qmflows.packages
================

The :class:`~qmflows.packages.packages.Package` (sub-)classes of QMFlows.

Package-related Functions
-------------------------
.. currentmodule:: qmflows.packages.packages
.. autosummary::
    run

The Package Class
-----------------
.. autosummary::
    Package
    Package.__init__
    Package.__repr__
    Package.__call__
    Package.prerun
    Package.run_job
    Package.postrun
    Package.generic2specific
    Package.handle_special_keywords

Package Subclasses
------------------
.. currentmodule:: qmflows.packages
.. autosummary::
    ~SCM.ADF
    ~SCM.DFTB
    ~cp2k_package.CP2K
    ~cp2k_mm.CP2KMM
    ~orca.ORCA
    ~package_wrapper.PackageWrapper

Package Instances
------------------
.. autosummary::
    ~SCM.adf
    ~SCM.dftb
    ~cp2k_package.cp2k
    ~cp2k_mm.cp2k_mm
    ~orca.orca

API
---
.. autofunction:: qmflows.packages.packages.run

|

.. autoclass:: qmflows.packages.packages.Package
    :members: generic_mapping, result_type, pkg_name

.. automethod:: qmflows.packages.packages.Package.__init__
.. automethod:: qmflows.packages.packages.Package.__repr__
.. automethod:: qmflows.packages.packages.Package.__call__
.. automethod:: qmflows.packages.packages.Package.prerun
.. automethod:: qmflows.packages.packages.Package.run_job
.. automethod:: qmflows.packages.packages.Package.postrun
.. automethod:: qmflows.packages.packages.Package.generic2specific
.. automethod:: qmflows.packages.packages.Package.handle_special_keywords

|

.. autoclass:: qmflows.packages.SCM.ADF

|

.. autoclass:: qmflows.packages.SCM.DFTB

|

.. autoclass:: qmflows.packages.cp2k_package.CP2K

|

.. autoclass:: qmflows.packages.cp2k_mm.CP2KMM

|

.. autoclass:: qmflows.packages.orca.ORCA

|

.. autoclass:: qmflows.packages.package_wrapper.PackageWrapper
    :noindex:

|

.. autodata:: qmflows.packages.SCM.adf
.. autodata:: qmflows.packages.SCM.dftb
.. autodata:: qmflows.packages.cp2k_package.cp2k
.. autodata:: qmflows.packages.cp2k_mm.cp2k_mm
.. autodata:: qmflows.packages.orca.orca
Extracting numerical properties from output files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Quantum packages simulations generate output file in different formats. For examples the SCM_ simulation suite
(:class:`~qmflows.packages.SCM.ADF` and :class:`~qmflows.packages.SCM.DFTB` in *QMFlows*) generate binary outputs,
while other packages like :class:`~qmflows.packages.cp2k_package.CP2K` and :class:`~qmflows.packages.orca.ORCA` generate ascii text files.

*QMFlows* abstract away all the different commmunication protocols with the different output formats, allowing the user to
extract the desire property by using the convention:

.. code:: python

    >>> result = job.property

where job is the simulation perform with a given package and property is the numerical value of interest (scalar or array).


The  *QMFlows* implementation of the aforemention mechanism search in the YAML files located at ``qmflows/data/dictionaries/``
for instructions about how to read that given property from the output file. Nevertheless, *Not all the properties for a given
pacakge are implemented*. **If the property of your interest is not available you can request it in the Qmflows issues page**.


Parsers
~~~~~~~
Internally *QMFlows* uses different mechanism to extract different properties from the output files. In the case of the :class:`~qmflows.packages.SCM.ADF` and
:class:`~qmflows.packages.SCM.DFTB` packages, *QMFlows* take advantages of the python interface to kftools_ files developed by the SCM_. In the case of *XML* output,
*QMFlows* direcltly uses the python built-in xml_ reader. For the output files in text format *Qmflows* uses a mixuture of awk_ and
*parsers*.

Parsers are a robust alternative to regular expressions, parsers are modular and reusable, while
|regex| tends to be abstruse and difficult to reuse. A parser is a function that decomposes a string (or binary) into its syntactic components using some predefined rules or grammar.
The library  pyparsing_ offers all the functionality to parse strings, some detail explanation about the library can be found at docs_.




.. _pyparsing: https://pyparsing.wikispaces.com/

.. _docs: https://pythonhosted.org/pyparsing/

.. |regex| replace:: :mod:`re`

.. _SCM: https://www.scm.com/

.. _KF: https://www.scm.com/doc/Scripting/Commandline_Tools/KF_command_line_utilities.html

.. _xml: https://docs.python.org/3.5/library/xml.etree.elementtree.html

.. _awk: https://www.gnu.org/software/gawk/manual/gawk.html

.. _properties: https://github.com/SCM-NV/qmflows/tree/master/qmflows/data/dictionaries

.. _kftools: https://www.scm.com/doc/plams/scm.html#kf-files
.. _dictionaries:

Dictionaries
~~~~~~~~~~~~
While :mod:`~templates` are used as defaults, the *YAML* files stored in the dictionaries folder are use for two mainly two purposes:
translation of the generic keywords and properties extraction.


to translate the generic keywords provided by the user to specific keywords used for each package.
For instance these files are used by the :class:`~qmflows.packages.packages.Package` class.

.. include:: ../README.rst
Library Documentation
=====================

For a more detailed description of **QMFlows** read the documentation

.. toctree::
   settings
   templates
   dictionaries
   packages
   promise
   parsers
   hdf5
   package_wrapper


.. toctree::
   :hidden:

   _packages

HDF5
====

.. currentmodule:: qmflows.hdf5.quantumHDF5

HDF5_ is a data model to store and represent complex numerical data in a hierarchical way.
HDF5_ is extremely optimized to perform fast I/O operations in large data set and it is
implemented as a library with APIs to Python, C, C++, etc. The python interface is called
pyh5_ and it was designed to run along with Numpy as Scipy. 


.. _HDF5: https://www.hdfgroup.org/HDF5/

.. _pyh5: http://www.h5py.org/
Running a workflow
~~~~~~~~~~~~~~~~~~
A workflow in **QMFlows** consist of a set of computations and the dependencies between them,
explicitly declared by the user in the python script. This dependecies and the relation between
them form an graph (specifically an acyclic direct graph) that represent these relations.

**QMFlows** Builds the aforemention graph in order to realize the workflow evaluation order.
For instance the figure below represent a simulation where firstly a molecular geometry optimization is carried out using the *ADF* package and
some user defined :class:`~qmflows.settings.Settings` for the *ADF* simulation package.
Subsequently, using the optimized molecular geometry from the previous step and
another :class:`~qmflows.settings.Settings` for an orca simulation a job to compute the molecular frequencies is carried out.

.. image:: _images/simple_graph.png

A python script corresponding with this graph can be

.. code:: python

   >>> from plams import Molecule
   >>> from qmflows import (adf, orca, run, Settings)

   >>> inp = Settings(...)
   >>> acetonitrile = Molecule(...)

   # ADF optimization
   >>> optmized_mol_adf = adf(inp, acetonitrile, job_name='acetonitrile_opt')

   # Orca Settings definition
   >>> s2 = Settings()
   >>> s2.specific.orca.main = "freq"
   >>> s2.specific.orca.basis.basis = 'sto_sz'
   >>> s2.specific.orca.method.functional = 'lda'
   >>> s2.specific.orca.method.method = 'dft'

   # Orca Frequencies
   >>> job_freq = orca(s2, optmized_mol_adf)

   # Extract the frequencies from the Orca job
   >>> frequencies = job_freq.frequencies

   # Run the graph
   >>> result = run(frequencies)
   >>> print(result)


Up to the invocation of the :func:`~qmflows.packages.packages.run` function none of the computations have been executed,
it is the :func:`~qmflows.packages.packages.run` function which builds and executes the dependencies.
Since *QMFlows* needs to figure out all the dependecies in the script,
the :func:`~qmflows.packages.packages.run` function takes as argument last dependency (or inner most dependy),
which in this case are the frequencies. The reason behind this, is that from the last dependency it is possible to
retrace all the dependecies.

**QMFlows** uses  the noodles_ library under the hook to takes care of the construction and
execution of the dependecy graph.

.. _noodles: http://nlesc.github.io/noodles/
PackageWrapper
==============

.. automodule:: qmflows.packages.package_wrapper
Templates
---------
The input generations consist of two parts: chosing a template (see templates_)
for the kind of calculation to perform and adding some settings to that template. Notice
that the user can either pick a specific package template or provides only generic
keywords.

.. _templates:

YAML
~~~~
.. currentmodule:: qmflows.templates.templates

The YAML_ markdown format is used together with the :mod:`yaml` module to implement the mechanism to load the templates.

.. _YAML: https://pyyaml.org/wiki/PyYAMLDocumentation

For Example, the default parameter for a single point calculation for several package is:

.. code-block:: python

    specific:
      adf:
         basis:
           type: SZ
         numericalquality:
           normal
         scf:
           converge: 1e-6
           iterations: 100
      ams:
         ams:
           Task: SinglePoint
      dftb:
        dftb:
          resourcesdir:
            "DFTB.org/3ob-3-1"
        task:
          runtype: SP

      cp2k:
        force_eval:
          dft:
            mgrid:
              cutoff: 400
              ngrids: 4
            print:
              mo:
                add_last: numeric
                each:
                  qs_scf: 0
                eigenvalues: ""
                eigenvectors: ""
                filename: "./mo.data"
                ndigits: 36
                occupation_numbers: ""
            qs:
                method: gpw
            scf:
                eps_scf: 1e-06
                max_scf: 200
                scf_guess: restart

          subsys:
            cell:
              periodic: xyz
        global:
          print_level: low
          project: cp2k
          run_type: energy

      orca:
        method:
            method: dft
            functional: lda
        basis:
            basis: sto_sz
.. qmworks documentation master file, created by
   sphinx-quickstart on Wed Dec 16 14:33:52 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to qmflows' documentation!
===================================

Contents:

.. toctree::
   :maxdepth: 2

   includeme
   interactive_tutorial
   Further_reading
   documentation


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Settings
--------

|Settings| is a :class:`dict` subclass implemented in PLAMS_ and modified in *Qmflows*.
This class represents the data in a hierarchical tree-like structure. for example:

.. code:: python

    >>> from qmflows import Settings, templates

    >>> s = Settings()  # (1)

    # generic keyword
    >>> s.basis = "DZP"  #  (2)

    # "specific" allows the user to apply specific keywords for a package
    >>> s.specific.adf.basis.core = "large"  # (3)

    >>> input_settings = templates.singlepoint.overlay(s)  # (4)


The above code snippet shows how to create a |Settings| instance object in **(1)**,
then in **(2)** the generic keyword *basis*  declares that the "DZP" should be used together with the *large* keyword
of *ADF* as shown at **(3)**.
Finally in line **(4)** the user's keywords are merged with the defaults resultin in a input like:


.. code::

    basis
        Core large
        Type DZP
    end

    integration
        accint 4.0
    end

    scf
        converge 1e-06
        iterations 100
    wnd

    xc
        lda
    end


API
~~~
.. autoclass:: qmflows.settings.Settings
    :members:
    :special-members:
    :exclude-members: __weakref__


.. _PLAMS: https://www.scm.com/doc/plams/components/settings.html
.. packages_:

Packages
========

The base class :class:`~qmflows.packages.packages.Package` is the library core, it provides the general scaffold to call a quantum code.
On top of this infrastructure it has been created several subclasses that contain the specific details for each quantum code.
The available interfaces to quantum codes are:

* :class:`~qmflows.packages.SCM.ADF`
* :class:`~qmflows.packages.SCM.DFTB`
* :class:`~qmflows.packages.cp2k_package.CP2K`
* :class:`~qmflows.packages.cp2k_mm.CP2KMM`
* :class:`~qmflows.packages.orca.ORCA`



**This class must not be call directly**, instead the correspoding class for the quantum package should be called or in case that there is not an interface to your quantum code,
you can make a new subclass that implement the following method:

* ``run_job`` -- This methods takes a :class:`~qmflows.settings.Settings` object a molecule and call a function to create the input automatically and takes cares of the bookkeeping associated with creating new folders, calling the package and retrieving and Object-result.


Instead of implementing all the runners and the bookkeeping functions ourselves, we use the plams_ library.

For example to carry out a simulation using ADF we call plams as follows:

.. code:: python

    >>> from scm import plams

    >>> mol = plams.Molecule(...)
    >>> adf_settings = plams.Settings(...)

    >>> result = plams.ADFJob(molecule=mol, settings=adf_settings).run()


.. _plams: https://www.scm.com/doc/plams/index.html
