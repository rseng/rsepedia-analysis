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

#. use the search functionality `here <https://github.com/https://github.com/nlesc-nano/flamingo/issues>`__ to see if someone already filed the same issue;
#. if your issue search did not yield any relevant results, make a new issue;
#. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

#. use the search functionality `here <https://github.com/https://github.com/nlesc-nano/flamingo/issues>`__ to see if someone already filed the same issue;
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
#. make sure the existing tests still work by running ``python setup.py test``;
#. add your own tests (if necessary);
#. update or expand the documentation;
#. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the flamingo repository on GitHub;
#. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


0.3.1 [unreleased]
******************
* *placeholder*.


0.3.0 [03/12/2021]
******************
* Release flamingo on pypi (#64)
* Run all the filters in parallel (#38)
* Filter molecules with a single functional group (#43)
* Add interface to cosmo-rs (#5)


0.2.1 [14/01/2021]
******************
Change
-----
* Use all the available CPU to compute bulkiness with CAT by calling the `imap_unordered Pool's method <https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.imap_unordered>`_.
* Remove the `batch_size` input parameter and fix it to 1000.


0.2.0 [03/11/2020]
******************

Added
-----
* Workflow to compute properties using `CAT <https://github.com/nlesc-nano/CAT>`_


0.1.0 [14/10/2020]
******************

Added
-----
* Move `swan functionality to compute and filter properties<https://github.com/nlesc-nano/swan/issues/44>`_ to **flamingo**
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
reported by contacting the project team at f.zapata@esciencecenter.nl. All
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
.. image:: https://github.com/nlesc-nano/flamingo/workflows/build%20with%20conda/badge.svg
   :target: https://github.com/nlesc-nano/flamingo/actions
.. image:: https://readthedocs.org/projects/flamingo-molecular-properties-calculator/badge/?version=latest
   :target: https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/?badge=latest
.. image:: https://codecov.io/gh/nlesc-nano/flamingo/branch/master/graph/badge.svg?token=N6CU1B82X0
   :target: https://codecov.io/gh/nlesc-nano/flamingo
.. image:: https://zenodo.org/badge/300545275.svg
   :target: https://zenodo.org/badge/latestdoi/300545275
.. image:: https://badge.fury.io/py/nlesc-flamingo.svg
   :target: https://badge.fury.io/py/nlesc-flamingo

########
flamingo
########

Compute and filter molecular properties. See `documentation <https://flamingo-molecular-properties-calculator.readthedocs.io/en/latest/>`_.

Installation
============

- Download miniconda for python3: miniconda_.

- Install according to: installConda_.

- Create a new virtual environment using the following commands:

  - ``conda create -n flamingo``

- Activate the new virtual environment

  - ``conda activate flamingo``

To exit the virtual environment type  ``conda deactivate``.


.. _dependecies:

Dependencies installation
-------------------------

- Type in your terminal:

  ``conda activate flamingo``

Using the conda environment the following packages should be installed:


- install RDKit_ and H5PY_:

  - `conda install -y -q -c conda-forge  h5py rdkit`

.. _installation:

Package installation
--------------------
Finally install the package:

- Install **flamingo** using pip:
  - ``pip install nlesc-flamingo``

Now you are ready to use *flamingo*.


Contributing
************

If you want to contribute to the development of flamingo,
have a look at the `contribution guidelines <CONTRIBUTING.rst>`_.

License
*******

Copyright (c) 2020-2021, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.



Credits
*******

This package was created with `Cookiecutter <https://github.com/audreyr/cookiecutter>`_ and the `NLeSC/python-template <https://github.com/NLeSC/python-template>`_.

.. _installConda: https://conda.io/projects/conda/en/latest/user-guide/install/index.html
.. _RDKit: https://www.rdkit.org
.. _H5PY: https://www.h5py.org/
.. _miniconda: https://docs.conda.io/en/latest/miniconda.html
Smiles Screener
===============

.. automodule:: flamingo.screen
Tutorial in Silico Screening
============================
This tutorial covers how to perform insilico filtering of a set of molecules
represented as smiles_. Table :ref:`smiles-table` contains some smiles
examples representing molecules with different functional groups.

.. _smiles-table:

.. csv-table:: smiles to filter
   :header: "smiles"

   CN1C=NC2=C1C(=O)N(C(=O)N2C)C
   OC(=O)C1CNC2C3C4CC2C1N34
   C1=CC=CC=C1
   OC(=O)C1CNC2COC1(C2)C#C
   CCO
   CCCCCCCCC=CCCCCCCCC(=O)O
   CC(=O)O
   O=C(O)Cc1ccccc1
   CC(C(=O)O)O


The filtering process consists in excluding (or including) a set of
molecules based on structural characteristics like their functional
groups or derived properties like bulkiness.

To run the simulation, the user must provide two files: one containing the
smiles that she wants to filter and another file containing
the values of the properties used as filters. 


Simulation input
****************
The smiles input should be in csv format like ::

  ,smiles
  ,CC(=O)O
  ,CCC(=O)O


The properties specification file to perform the filtering must be a yaml
file following the subsequent schema yaml_::

 smiles_file:
   smiles.csv

 core:
   "Cd68Se55.xyz"

 anchor:
   "O(C=O)[H]"

 batch_size: 1000
    
 filters:
   include_functional_groups:
     groups:
        - "[CX3](=O)[OX2H1]" # Include carboxylic acids
     maximum: 1
   exclude_functional_groups:
     groups:
        - "[NX3]"  # Exclude tertiary amines
        - "C#C"    # Exclude triplet Carbon-Carbon bonds
   scscore:
     lower_than:
       3.0
   bulkiness:
     lower_than: 20


The *smiles_file* entry contains the path to the files containing the smiles. The
other keywords will be explain in the following sections.
	
Available filters
*****************

.. Note:: The filters are run in sequential order, meaning that second filter is applied
   to the set of molecules remaining after applying the first filters, the third
   filter is applied after the second and so on.


1. Include and exclude function groups
--------------------------------------
The *include_functional_groups* and *exclude_functional_groups* as their names suggest
keep and drop molecules based on a list of functional groups represented as
`SMARTS <https://en.wikipedia.org/wiki/SMILES_arbitrary_target_specification>`_.

the *maximum* keyword indicates what is the maximum number of functional groups
that can be present.

2. Synthesizability scores
--------------------------
The scscore_ is a measure of synthetic complexity. It is scaled from 1 to 5
to facilited human interpretation. See the scscore_ paper for further details.


3. Bulkiness
------------
Assuming that a given molecule can be attached to a given surface, the bulkiness
descriptor gives a measure of the volumen occupied by the molecule from the
anchoring point extending outwards as a cone. It requires the *core* keywords
specifying the surface to attach the molecule and the *anchor* functional
group used as attachment.
See the `the CAT bulkiness <https://cat.readthedocs.io/en/latest/4_optional.html?highlight=bulkiness#optional.qd.bulkiness>`_
for further information.

	
Running the filtering script
****************************
To perform the screening you just need to execute the following command ::
  smiles_screener -i path_to_yaml_input.yml


Job distributions and results
*****************************
For a given filter, **Flamingo** will try to compute the molecular properties in parallel since properties
can be computed independently for each molecule. Therefore **Flamingo** split the molecular set
into batches that can be computed in parallel. The `batch_size` keyword is used to control the
size of these batches.

After the computation has finished the filtered molecules are stored in the **results** folder
in the *current work directory*. In that folder you can find a `candidates.csv` file for
each batch containing the final molecules.

.. _smiles: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
.. _yaml: https://yaml.org/
.. _scscore: https://pubs.acs.org/doi/10.1021/acs.jcim.7b00622
Molecular features
==================

.. automodule:: flamingo.features.featurizer
.. flamingo documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to flamingo's documentation!
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   includereadme
   tutorial_screening
   api_molecular_features
   api_filtering
   cat_interface

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
CAT interface
=============

.. automodule:: flamingo.cat_interface
.. include:: ../README.rst
