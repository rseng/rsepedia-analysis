############################
Contributing guidelines
############################

We welcome any kind of contribution to our software, from simple comment or question to a full fledged `pull request <https://help.github.com/articles/about-pull-requests/>`_. Please read and follow our `Code of Conduct <CODE_OF_CONDUCT.rst>`_.

A contribution can be one of the following cases:

1. you have a question;
1. you think you may have found a bug (including unexpected behavior);
1. you want to make some kind of change to the code base (e.g. to fix a bug, to add a new feature, to update documentation).

The sections below outline the steps in each case.

You have a question
*******************

1. use the search functionality `here <https://github.com/nlesc-nano/CAT/issues>`__ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue;
1. apply the "Question" label; apply other labels when relevant.

You think you may have found a bug
**********************************

1. use the search functionality `here <https://github.com/nlesc-nano/CAT/issues>`__ to see if someone already filed the same issue;
1. if your issue search did not yield any relevant results, make a new issue, making sure to provide enough information to the rest of the community to understand the cause and context of the problem. Depending on the issue, you may want to include:
    - the `SHA hashcode <https://help.github.com/articles/autolinked-references-and-urls/#commit-shas>`_ of the commit that is causing your problem;
    - some identifying information (name and version number) for dependencies you're using;
    - information about the operating system;
1. apply relevant labels to the newly created issue.

You want to make some kind of change to the code base
*****************************************************

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue;
1. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
1. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions `here <https://help.github.com/articles/configuring-a-remote-for-a-fork/>`__ and `here <https://help.github.com/articles/syncing-a-fork/>`__);
1. make sure the existing tests still work by running ``python setup.py test``;
1. add your own tests (if necessary);
1. update or expand the documentation;
1. `push <http://rogerdudler.github.io/git-guide/>`_ your feature branch to (your fork of) the Compound Attachment/Analysis Tool repository on GitHub;
1. create the pull request, e.g. following the instructions `here <https://help.github.com/articles/creating-a-pull-request/>`__.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.
##########
Change Log
##########

All notable changes to this project will be documented in this file.
This project adheres to `Semantic Versioning <http://semver.org/>`_.


0.10.6
******
* *placeholder*.


0.10.5
******
* Fix an issue where certain properties could not be extracted from CP2K frequency jobs.
* Add the ``anchor.multi_anchor_filter`` option.


0.10.4
******
* Added the ``qd.dissociate.qd_opt`` keyword.
* Fix thermochemical properties not properly being set to ``nan`` for crashed jobs.


0.10.3
******
* Allow ``core.allignment: "surface"`` for cores with <4 anchor atoms.
* Allow QMFlows-style keywords to be parsed by the BDE workflow.
* Fix the ``qd.dissociate.core_atom`` and `lig_count` keys not being optional.
* Added the ``qd.dissociate.xyn_pre_opt`` keyword.


0.10.2
******
* Added the option to specify dummy atoms as anchor group.
* Added the option to specify the format type of anchor groups (SMILES, SMARTS, *etc.*).
* Fixed an issue wherein anchor-group-parsing could fail for mono-atomic anchors.
* Added ``packaging`` and ``rdkit-pypi`` as dependencies;
  manually installing ``rdkit`` via conda is no longer necessary.


0.10.1
******
* Added the option to manually specify angle offsets the ligand vectors.
* Added the option to manually specify dihedral angle the ligand vectors.
* Added the option to invert the core vectors.
* Added the various ligand anchor parsing options to its core-based counterpart.
* Deprecated usage of ``Molecule.get_formula`` in favor of a PLAMS <=1.5.1-based backport.
* Fixed CAT modifiying global ``logging`` settings.


0.10.0
******
* Add more advanced anchor-parsing options.
* Update for AMS 2021.
* Fixed various issues.


0.9.11
******
* Changed the default ``optional.core.allignment`` option from ``"sphere"`` to ``"surface"``.


0.9.10
******
* Added the new ``optional.core.anchor`` and ``optional.ligand.anchor`` options,
  which are respectively aliases for the old ``optional.core.dummy`` and
  ``optional.ligand.functional_groups`` options.
* Added a new H3O+ .coskf file.
* Added a new template for LogP calculations.
* Enabled tests for Windows.
* Fixed a recent readthedocs failure.


0.9.9
*****
* Improve the conformations of ringed systems;
  ring-substituents are now treated as fragments.


0.9.8
*****
* For some reason the matplotlib logger is going bananas;
  its level has been changed from DEBUG (default) to INFO in order to silence it.


0.9.7
*****
* Update the documentation for https://github.com/nlesc-nano/data-CAT/pull/38.
* Relax the pandas version requirement.


0.9.6
*****
* Updated the documentation for https://github.com/nlesc-nano/data-CAT/pull/25.


0.9.5
*****
* Fixed a bug where rdkit molecules were not properly converted into numpy arrays.
* Only perform `optional.ligand.optimize.job2` if the preceding UFF optimization finishes without crashing.
* Remove a ligand if its optimization fails (_i.e._ an exception is raised).
* Fixed an issue where ligand anchoring groups could not be explicitly specified
  (`docs <https://cat.readthedocs.io/en/latest/3_input_core_ligand.html#indices>`_).
* Improved the conformations of ligands with tri- and penta-valent pnictogen anchors.


0.9.4
*****
* Fixed an issue where certain `AMSJobs` would have duplicate keys.
* Fixed an issue where certain ligands weren't exported when constructing
  multi-ligand quantum dots (see https://github.com/nlesc-nano/CAT/issues/115).
* Re-enabled tests for Python 3.6 (see https://github.com/nlesc-nano/CAT/issues/125).


0.9.3
*****
* Switched from travis to GitHub Actions.
* Added tests using flake8, pydocstyle and doctest.


0.9.2
*****
* Updated the Database documentation (see https://github.com/nlesc-nano/data-CAT/pull/23).
* Upped the minimum Sphinx version to 2.1.
* Removed ``sphinx-autodoc-typehints``.
* Removed ``requirements.txt`` in favor of ``.readthedocs.yml``.


0.9.1
*****
* Added a new conceptual DFT (CDFT) workflow to Nano-CAT
  (https://github.com/nlesc-nano/nano-CAT/pull/57).


0.9.0
*****
* Moved a number of functions to the `nanoutils <https://github.com/nlesc-nano/Nano-Utils>`_ package.


0.8.9
*****
* Added the option to perform ligand geometry optimizations at
  arbitrary levels of theory.


0.8.8
*****
* Fixed a bug where ligands weren't properly rotated when
  using ``optional.ligand.optimize = False``.
* Added tests for ``CAT.utils.SetAttr``.


0.8.7
*****
* Replaced ``print()`` calls with ``logger.warning()`` in all dye-related functions.
* Added a recipe for calculating the synthetic accessibility score (SAS)
  as adapted from `MolGAN <https://github.com/nicola-decao/MolGAN>`_.


0.8.6
*****
* Added documentation for the new Nano-CAT ``multi_ligand_job()`` recipe.


0.8.5
*****
* Version bump.


0.8.4
*****
* Turned the ``dye`` functionality into a recipe in ``CAT.recipes``.


0.8.3
*****
* Merged all features from the ``dye`` branch into the master.
* Fixed an issue where custom forcefield settings are not properly parsed:
  https://github.com/nlesc-nano/CAT/pull/99.
* Added a try/except clause for job hashing in case rerun prevention is disabled:
  https://github.com/nlesc-nano/CAT/pull/98.
* Added new recipes to the documentation:
  https://github.com/nlesc-nano/CAT/pull/95 & https://github.com/nlesc-nano/CAT/pull/96.
* Fixed an issue where creating an object array would unpack a Molecule into Atoms:
  https://github.com/nlesc-nano/CAT/pull/94.
* Raise an Exception when failing to identify any atoms:
  https://github.com/nlesc-nano/CAT/pull/93.


0.8.2
*****
* Added the option to decorate a qd surface with more than one type of ligand.


0.8.1
*****
* Added the ``optional.core.allignment`` keyword for determining how
  ligands should be alligned with the core.
  Accepted values are ``"sphere"`` and ``"surface"``.
* https://github.com/nlesc-nano/CAT/pull/87:
  Ensure that part of the core-surface is accounted for when rotating ligands.
* https://github.com/nlesc-nano/CAT/pull/85 & https://github.com/nlesc-nano/CAT/pull/86:
  Issue a warning when atoms are too close when constructing QDs.
* https://github.com/nlesc-nano/CAT/pull/85 & https://github.com/nlesc-nano/CAT/pull/86:
  Improved warning handling.


0.8.0
*****
* Moved the ``CAT.recipes`` module to Nano-CAT.
* Moved the ``CAT.attachment.qd_opt_ff`` module to Nano-CAT.
* Created the ``CAT.workflow.key_map module`` for storing aliases
  for ``DataFrame()`` columns.
* Cleaned the modules in ``CAT.workflows``.
* Updated tests.


0.7.15
******
* Moved ``test_distribute()`` to it's own module: ``CAT.attachment.distribution_utils``.
* Added the ``brute_uniform_idx()`` for creating uniform/clustered distributions
  in a brute-force manner, *i.e.* by finding the global minimum/maximum within
  the set of all valid atom combinations.
* Generalized the ``array_combinations()`` function, it now accepts any
  array-like object and can generate combinations along any user-specified axis.
* Added the ``get_nearest_neighbors()`` function for finding the ``k``
  nearest-neighbors within a molecule.
* Added a recipe for marking a (sub-)set of surface atoms:
  ``CAT.recipes.mark_surface()``.
* Added a recipe for dissociating specific sets of surface atoms:
  ``CAT.recipes.dissociate_surface()``.
* Update to the general structure of the ``CAT.recipes`` modules.
* Multiple minor documentation adjustments.


0.7.14
******
* Changed the default value of the CP2K ``EI_SCALE14`` keyword from 0.0 to 1.0
  (*i.e.* the CHARMM forcefield default).
* Renamed the CAT ``activation_strain.scale_elstat`` keyword to ``.el_scale14``.
* Renamed the CAT ``activation_strain.scale_lj`` keyword to ``.lj_scale14``.
* Added the CAT ``activation_strain.dump_csv`` keyword for writing the raw
  potential energies to a set of .csv files.
* Added the CAT ``activation_strain.shift_cutoff`` keyword.
  Sets the value of all non-bonded potential to zero at ``activation_strain.distance_upper_bound``.
* A number of consistency improvements to the Schemas.


0.7.13
******
* Small optimization improvements to ``edge_dist()``.
* Moved a number of functions around in the CAT.utils module.
* Added the ``optional.qd.dissociate.lig_pairs`` keyword for the BDE workflow.


0.7.12
******
* Fixed a bug ``qd_opt_ff()`` where the wrong dictionary key was validated.
* Multiple updates to the CP2K MD template.
* Employ a more duck-typing based approach during the ``schema`` validation.
* Fixed a bug in the ``jobs`` module where incorrect ``Results()`` instances
  were returned.
* Multiple documentation updates.


0.7.11
******
* Updated the ``CAT.attachment.qd_opt_ff`` module in preparation for
  https://github.com/nlesc-nano/nano-CAT/pull/26.


0.7.10
******
* The function for applying distance weights during the
  subset-generation process is now configurable.
* The default distance weighting function has been changed to
  ``weight = "np.exp(-x)"``.
  The old p-norm with ``p=-2`` is still accessible via: ``weight = "x**-2"``


0.7.9
*****
* Added the option to interpolate between ``"uniform"`` / ``"cluster"`` and
  ``"random"``.
* The order of the ``p``-norm is now configurable.
* The variable representing the anchor-atom subset size has been changed
  from ``p`` to ``f``.
  ``p`` is now reserved for the order of the ``p-norm``.
* https://github.com/nlesc-nano/CAT/pull/70: Fixed an issue with the
  ``_parse_cluster_size()`` index offset.


0.7.8
*****
* It is now possible to create ``"uniform"`` distributions of clusters,
  the size of each cluster being user-specified.


0.7.7
*****
* The ``"uniform"`` and ``"cluster"`` distributions are now weighted by
  the distance rather than using a, less robust, distance truncation.


0.7.6
*****
* Added the option, when constructing core atom subsets,
  the use a distance matrix representing the shortest paths along the
  edges of a polyhedron, rather than through space.
  Enabling this option will result in more accurate ``"uniform"`` and
  ``"cluster"`` distributions at the cost of increased computational time.
* Updated and improved the ``"uniform"`` and ``"cluster"`` distributions.
* https://github.com/nlesc-nano/CAT/pull/65: Fixed a bug where ``uniform_idx()`` yielded the rolled,
  rather than unshifted, indices.
* https://github.com/nlesc-nano/CAT/pull/64: Bug fix: the subset Schema now checks for instances of
  int ``Or`` float.
* https://github.com/nlesc-nano/CAT/pull/66: Return the identity (rotation) matrix if a ``FloatingPointError`` is
  encountered during the creation of rotation matrices.
  This can occur if a ligand consists of a single atom.
* https://github.com/nlesc-nano/CAT/pull/66: Fixed a bug in the parsing of the mode parameter of ``distribute_idx()``;
  ``"uniform"`` and ``"cluster"`` will now correctly link to ``np.argmax`` and
  ``np.argmin`` instead of the other way around.


0.7.5
*****
* Added the ability to populate only a (random-ish) subset of
  core anchors with ligands.


0.7.4
*****
* The ligand rotation check is now substantially faster:
  a distance cutoff has been implemented for the construction
  of distance matrices.


0.7.3
*****
* Added an option perform an ensemble-averaged QD activation strain
  analyses in Nano-CAT_.
* Removed a number of redundant modules.
* QD optimization now properly respect the ``optional.qd.opt.use_ff`` keyword.


0.7.2
*****
* Minor tweaks to the default forcefield-related CP2K input files.
* Fixed a couple of bugs in the ligand dissociation workflow.
* Reworked the ligand dissociation procedure in Nano-CAT_.


0.7.1
*****
* Bug fix: Added a missing value to the to-be exported ASA columns.


0.7.0
*****
* Finalize the introduction of a new CAT template system (``WorkFlow()``).
* WiP: Implement an acitvation strain workflow with custom MATCH-based
  forcefields in Nano-CAT_.


0.6.5
*****
* Updated Nano-CAT to 0.2.4: https://github.com/nlesc-nano/nano-CAT/pull/20.
* Updated Data-CAT to 0.1.5: https://github.com/nlesc-nano/data-CAT/pull/17.
* Import assertions from AssertionLib_ rather than CAT_.
* Simplified to ``AsArray()`` context manager.
* Added the ``["keep_files"]`` option for quantum dot optimizations.
* Removed ``CRSJob()`` and ``CRSResults()``; import them from PLAMS_ instead.
* WiP: Introduction of a new CAT template system (``WorkFlow()``).


0.6.4
*****
* Moved the ligand bulkiness workflow from the `ligand` to the `qd` block
  in the CAT input. See `nano-CAT`_ 0.2.3.
* Updated the formula for the ligand bulkiness calculation.
  See `nano-CAT`_ 0.2.3.


0.6.3
*****
* Fixed a bug where hypervalent atoms where assigned incorrect atomic charges.


0.6.2
*****
* Added multiple improvements (and bug fixes) to the
  ligand conformation optimizer.
* Added a context manager for the `plams.Molecule.as_array()` method.
* Added an optimizer for the ligand vector.
* Updated the ligand bulkiness workflow in `nano-CAT`_ 0.2.2.


0.6.1
*****
* Added a workflow for calculating ligand bulkiness in `nano-CAT`_ 0.2.1.


0.6.0
*****
* Implemented an interface to MATCH_ (Multipurpose Atom-Typer for CHARMM)
  in Nano-CAT.
* Added a workflow for creating CP2K input files with
  the MATCH-assigned atom types & charges.
* Updated the handling of assertions, see ``CAT.assertions.assertion_manager``.


0.5.5
*****
* Lowered Python version requirement from >=3.7 to >=3.6.


0.5.4
*****
* Minor updates to the logger.
* Cleaned up CAT.jobs.py.
* ``check_sys_var()`` is now only called if an ADF-specific Job is requirest.
* Job hashes are now stored in (and retrieved from) $JN.hash files (plain text).
* Added a permanent Database_ instance to .optional.database.db.
* Parsing of functional group SMILES_ strings is now carried out during the Schema_ validation.
* Updated Data-CAT_ to 0.1.2; changed status from pre-alpha to alpha
  (see https://github.com/nlesc-nano/data-CAT/pull/13).



0.5.3
*****
* Moved Molecule to file exporting (*i.e.* .xyz and .pdb creation) from data-CAT_ to CAT_.
* Molecules can now be exported to .mol and .mol2 formats (in addition to .pdb and .xyz format).
* Increased the clarity of many exceptions (see https://github.com/nlesc-nano/CAT/issues/45).
* Updated the documentation.
* Introduced a proper logger (see https://github.com/nlesc-nano/CAT/issues/46).
* Updated data-CAT_ to 0.1.1 (https://github.com/nlesc-nano/data-CAT/pull/12) and
  nano_CAT_ to 0.1.2 (https://github.com/nlesc-nano/nano-CAT/pull/10).


0.5.2
*****
* Added more tests.
* Added a more explicit error message to ``_smiles_to_rdmol()``.


0.5.1
*****
* Documentation update.
* Updated to the ligand dissociation module in nano-CAT_ (see https://github.com/nlesc-nano/nano-CAT/issues/1).
* Added the ``keep_files`` keyword to the cosmo-rs and ligand dissociation workflows.
  Default value: ``True``.
* See https://github.com/nlesc-nano/nano-CAT/pull/9.


0.5.0
*****
* CAT_ has been split into 3 seperate packages (see https://github.com/nlesc-nano/CAT/issues/39):

  * CAT_: A collection of tools designed for the automatic construction of composite chemical compounds.
  * nano-CAT_: A collection of tools for the analysis of nanocrystals.
  * data-CAT_: A databasing framework for the Compound Attachment Tools package (CAT_).

* Docstrings have been changed into NumPy style.
* Added typehints.
* Added the CAT.SettingsDataFrame and CAT.SettingsSeries classes.
* Added more tests.
* Cleaned up all input-parsing related modules.
* Custom function groups (*i.e.* SMILES_ strings) can now be specified in the input
  under the optional.ligand.functional_groups key (see https://github.com/nlesc-nano/CAT/issues/13).


0.4.6
*****
* Added an interface between MongoDB_ and the CAT.Database_ class (see https://github.com/nlesc-nano/CAT/issues/11).


0.4.5
*****
* All raw input scripts are now stored in the structures.hdf5 file
  (see: https://github.com/nlesc-nano/CAT/issues/36).


0.4.4
*****
* Split CAT_database.py into database.py and database_functions.py.
* Unoptimized starting structures are now exported to the database.
* Added the sphinx autosummary extension.


0.4.3
*****
* Improved interaction between the database and BDE module.
* Cleaned up BDE module.
* HDF5 indices are now always sorted when itneraction with the database.


0.4.2
*****
* Numerous bug fixes.
* A couple of code-style changes.


0.4.1
*****
* COSMO-RS calculations now allow for COSMO-surface construction
  at the DFT level.


0.4.0
*****
* Introduction of the CAT.Database class.
* Central object of CAT has been changed into a dataframe of
  molecules rather than lists molecules.
* Updated a number of tests.


0.3.3
*****
* Changed qmflows template import syntax (see: https://github.com/SCM-NV/qmflows/pull/132).
* Changed yaml loader.


0.3.2
*****
* Further (minor) updates and bug fixes to the database interaction.
* Overhaul of the bond dissociation energy (BDE) module.
* Job settings are now stored in the database.


0.3.0
*****
* Massive overhaul of the CAT database interaction.
* Moved functions related to functiona group recognizition to
  CAT.attachment.ligand_anchoring.py.
* Multiple minor bug fixes.


[Unreleased]
************
* Empty Python project directory structure.


.. _AssertionLib: https://github.com/nlesc-nano/AssertionLib
.. _CAT: https://github.com/nlesc-nano/CAT
.. _CAT.Database: https://cat.readthedocs.io/en/latest/7_database.html
.. _CP2K: https://www.cp2k.org/
.. _data-CAT: https://github.com/nlesc-nano/data-CAT/
.. _Database: https://cat.readthedocs.io/en/latest/7_database.html#class-api
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _MATCH: http://brooks.chem.lsa.umich.edu/index.php?page=match&subdir=articles/resources/software
.. _MongoDB: https://www.mongodb.com/
.. _nano-CAT: https://github.com/nlesc-nano/nano-CAT/
.. _Schema: https://github.com/keleshev/schema
.. _SMILES: https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system
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
reported by contacting the project team at b.f.van.beek@vu.nl. All
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
.. image:: https://github.com/nlesc-nano/CAT/workflows/Build%20with%20Conda/badge.svg
   :target: https://github.com/nlesc-nano/CAT/actions?query=workflow%3A%22Build+with+Conda%22
.. image:: https://readthedocs.org/projects/cat/badge/?version=latest
   :target: https://cat.readthedocs.io/en/latest/
.. image:: https://zenodo.org/badge/169708827.svg
   :target: https://zenodo.org/badge/latestdoi/169708827
.. image:: https://badge.fury.io/py/nlesc-CAT.svg
   :target: https://badge.fury.io/py/nlesc-CAT

|

.. image:: https://img.shields.io/badge/python-3.6-blue.svg
   :target: https://docs.python.org/3.6/
.. image:: https://img.shields.io/badge/python-3.7-blue.svg
   :target: https://docs.python.org/3.7/
.. image:: https://img.shields.io/badge/python-3.8-blue.svg
   :target: https://docs.python.org/3.8/
.. image:: https://img.shields.io/badge/python-3.9-blue.svg
   :target: https://docs.python.org/3.9/
.. image:: https://img.shields.io/badge/python-3.10-blue.svg
   :target: https://docs.python.org/3.10/

###############################
Compound Attachment Tool 0.10.6
###############################

**CAT** is a collection of tools designed for the construction of various chemical compounds.
Further information is provided in the documentation_.

Package installation
--------------------
**CAT** can be installed via pip as following:

- **CAT**: ``pip install nlesc-CAT --upgrade``

Note that, while not strictly necessary, it is recommended to first create a conda environment:

- Download and install miniconda for python3: miniconda_ (also you can install the complete anaconda_ version).

- Create a new virtual environment:  ``conda create --name CAT python``

- Activate the environment:: ``conda activate CAT``

Input files
============

Running **CAT** and can be done with the following command:
``init_cat my_settings.yaml``. The user merely has to provide a yaml_ file
with the job settings, settings which can be tweaked and altered to suit ones
purposes (see example1_). Alternatively, **CAT** can be run like a regular
python script, bypassing the command-line interface
(*i.e.* ``python input.py``, see example2_).

An extensive description of the various available settings is available in
the documentation_.


.. _yaml: https://yaml.org/
.. _documentation: https://cat.readthedocs.io/en/latest/
.. _miniconda: http://conda.pydata.org/miniconda.html
.. _anaconda: https://www.continuum.io/downloads
.. _installConda: https://docs.anaconda.com/anaconda/install/
.. _HDF5: http://www.h5py.org/
.. _here: https://www.python.org/downloads/
.. _rdkit: http://www.rdkit.org
.. _PLAMS: https://github.com/SCM-NV/PLAMS
.. _QMFlows: https://github.com/SCM-NV/qmflows
.. _example1: https://github.com/BvB93/CAT/blob/master/examples/input_settings.yaml
.. _example2: https://github.com/BvB93/CAT/blob/master/examples/input.py
###
CAT
###

~~~~~~~~~
analysis_
~~~~~~~~~

Modules related to the analysis ligands and quantum dots.

~~~~~~~~~~~
attachment_
~~~~~~~~~~~

Modules related to processing ligands,
cores and the subsequent creation of quantum dots.

~~~~~
data_
~~~~~

Various .yaml templates and COSMO-MOPAC .coskf files.

~~~~~~~~~~~~~~
data_handling_
~~~~~~~~~~~~~~

Modules related to the importing, exporting and general handling of data.

~~~~~~~~
base.py_
~~~~~~~~

A module handling the interaction with all other modules,
functioning as recipe.

~~~~~~~~~
utils.py_
~~~~~~~~~

A module with miscellaneous functions.

~~~~~~~~~~~~~
mol_utils.py_
~~~~~~~~~~~~~

A module with misc functions related to manipulating molecules,
atoms and bonds.

.. _analysis: https://github.com/BvB93/CAT/tree/master/CAT/analysis
.. _attachment: https://github.com/BvB93/CAT/tree/master/CAT/attachment
.. _data: https://github.com/BvB93/CAT/tree/master/CAT/data
.. _data_handling: https://github.com/BvB93/CAT/tree/master/CAT/data_handling
.. _base.py: https://github.com/BvB93/CAT/tree/master/CAT/base.py
.. _utils.py: https://github.com/BvB93/CAT/tree/master/CAT/utils.py
.. _mol_utils.py: https://github.com/BvB93/CAT/tree/master/CAT/mol_utils.py
##########
Attachment
##########

~~~~~~~~~~~~~~~~~
ligand_anchor.py_
~~~~~~~~~~~~~~~~~

A module designed for finding and processing functional groups within ligands.

~~~~~~~~~~~~~~~~~
ligand_attach.py_
~~~~~~~~~~~~~~~~~

A module designed for attaching ligands to cores.

~~~~~~~~~~~~~~
ligand_opt.py_
~~~~~~~~~~~~~~

A module designed for optimizing the geometry of ligands.

~~~~~~~~~~
qd_opt.py_
~~~~~~~~~~

A module designed for optimizing the combined ligand & core.

.. _ligand_anchor.py: https://github.com/BvB93/CAT/tree/master/CAT/attachment/ligand_anchor.py
.. _ligand_attach.py: https://github.com/BvB93/CAT/tree/master/CAT/attachment/ligand_attach.py
.. _ligand_opt.py: https://github.com/BvB93/CAT/tree/master/CAT/attachment/ligand_opt.py
.. _qd_opt.py: https://github.com/BvB93/CAT/tree/master/CAT/attachment/qd_opt.py
####
data
####

~~~~~
coskf
~~~~~

Various solvent .coskf files produced at the COSMO-MOPAC(PM7) level of theory.

~~~~~~~~~
templates
~~~~~~~~~
Various templates input scripts and SMILES strings used internally by **CAT**.

~~~~~
CORES
~~~~~

Various .xyz files of example cores.

~~~~~~~
LIGANDS
~~~~~~~
Various .xyz files of example ligands.
#############
data_handling
#############

~~~~~~~~~~~~~~
mol_import.py_
~~~~~~~~~~~~~~

A module related to the importing of molecules.

~~~~~~~~~~~~~~~~~~~
input_sanitizer.py_
~~~~~~~~~~~~~~~~~~~

A module designed for sanitizing and interpreting the input file.

~~~~~~~~~~~~~~~~
input_parser.py_
~~~~~~~~~~~~~~~~

A module designed for parsing the input .yaml file.

.. _mol_import.py: https://github.com/BvB93/CAT/tree/master/CAT/data_handling/mol_import.py
.. _input_sanitizer.py: https://github.com/BvB93/CAT/tree/master/CAT/data_handling/input_sanitizer.py
.. _input_parser.py: https://github.com/BvB93/CAT/tree/master/CAT/data_handling/input_parser.py
########
examples
########

~~~~~~~~~~
input.yaml
~~~~~~~~~~

Example input file settings.

~~~~~~~~
input.py
~~~~~~~~

An example input file.

~~~~~~~~~
addlig.py
~~~~~~~~~

An example input file for creating dyes.

~~~~~~~~~~~~~~
example_xyz.py
~~~~~~~~~~~~~~

Contains a function for accessing example core & ligand .xyz files
.. _Path:

path
====

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    path: null

|

Arguments
~~~~~~~~~

.. attribute:: path

    :Parameter:     * **Type** - :class:`str` or :class:`NoneType`
                    * **Default value** – ``None``

    The path were all working directories are/will be stored.
    To use the current working directory, use one of the following values:
    ``None``, ``"."``, ``""`` or ``"path_to_workdir"``.

    .. note::
        The yaml format uses ``null`` rather than ``None`` as in Python.
.. _HDF5 logging:

HDF5 Access Logging
===================
.. automodule:: dataCAT.hdf5_log
.. automodule:: nanoCAT.recipes.mark_surface
.. _import_qd:

Importing Quantum Dots
======================

*WiP*: Import pre-built quantum dots rather than constructing them from scratch.


Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    input_qd:
        - Cd68Se55_ethoxide.xyz:
            ligand_smiles: '[O-]CC'
            ligand_anchor: '[O-]'

|

Arguments
~~~~~~~~~

.. attribute:: ligand_smiles

    :Parameter:     * **Type** - :class:`str`
                    * **Default value** – ``None``

    A SMILES string representing the ligand.
    The provided SMILES string will be used for identifying the core and all ligands.

    .. warning::
        This argument has no value be default and thus *must* be provided by the user.


.. attribute:: ligand_anchor

    :Parameter:     * **Type** - :class:`str`
                    * **Default value** – ``None``

    A SMILES string representing the achor functional group of the ligand.
    If the provided SMILES string consists of multiple atoms
    (*e.g.* a carboxylate: ``"[O-]C=O"``), than the first atom will be treated as anchor (``"[O-]"``).

    .. warning::
        This argument has no value be default and thus *must* be provided by the user.
.. automodule:: nanoCAT.recipes.bulk
.. _recipes:

Recipes
=======
.. toctree::

    12_1_recipes.rst
    12_2_recipes.rst
    12_3_recipes.rst
    12_4_recipes.rst
    12_5_recipes.rst
    12_6_recipes.rst
    12_7_recipes.rst
    12_8_recipes.rst
    12_9_recipes.rst
    12_10_recipes.rst
    12_11_recipes.rst
.. automodule:: nanoCAT.recipes.coordination_number
.. _HDF5 property storage:

HDF5 Property Storage
=====================
.. automodule:: dataCAT.property_dset
.. _Gettings Started:

General Overview & Getting Started
==================================

A basic recipe for running **CAT**:

1.  Create two directories named ‘core’ and ‘ligand’. The 'core' directory
should contain the input cores & the 'ligand' should contain the input
ligands. The quantum dots will be exported to the 'QD' directory.

2. 	Customize the job settings to your liking, see
CAT/examples/input_settings.yaml_ for an example.
Note: everything under the ``optional`` section does **not** have to be
included in the input settings.
As is implied by the name, everything in ``optional`` is completely optional.

3.  Run **CAT** with the following command:
``init_cat input_settings.yaml``

4.  Congratulations, you just ran
**CAT**!

The default **CAT** settings, at various levels of verbosity, are provided
below.

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O


Verbose default Settings
~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O

    optional:
        database:
            dirname: database
            read: True
            write: True
            overwrite: False
            thread_safe: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            anchor: Cl
            subset: null

        ligand:
            dirname: ligand
            optimize: True
            split: True
            anchor: null
            cosmo-rs: False

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            bulkiness: False
            activation_strain: False
            dissociate: False

Maximum verbose default Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: yaml

    path: None

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O

    optional:
        database:
            dirname: database
            read: (core, ligand, qd)
            write: (core, ligand, qd)
            overwrite: False
            thread_safe: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            anchor: Cl
            subset: null

        ligand:
            dirname: ligand
            split: True
            anchor: null
            cosmo-rs: False
            optimize:
                use_ff: False
                job1: null
                s1: null
                job2: null
                s2: null

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            bulkiness: False
            activation_strain: False
            dissociate:
                core_atom: Cd
                lig_count: 2
                keep_files: True
                core_core_dist: 5.0
                lig_core_dist: 5.0
                topology: {}

                job1: False
                s1: False
                job2: False
                s2: False

.. _input_settings.yaml: https://github.com/BvB93/CAT/blob/devel/examples/input_settings.yaml
.. _Input Cores and Ligands:

input_cores & input_ligands
===========================

Thia section related relates the importing and processing of cores and ligands.
Ligand & cores can be imported from a wide range of different files and files
types, which can roughly be divided into three categories:

1.  Files containing coordinates of a single molecule: .xyz, .pdb & .mol files.
2.  Python objects: :class:`plams.Molecule`, :class:`rdkit.Chem.Mol` & SMILES strings (:class:`str`).
3.  Containers with one or multiple input molecules: directories & .txt files.

In the later case, the container can consist of multiple SMILES strings or
paths to .xyz, .pdb and/or .mol files. If necessary, containers are searched
recursively. Both absolute and relative paths are explored.

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(C)=O
        - OC(CC)=O
        - OC(CCC)=O
        - OC(CCCC)=O

Optional arguments
~~~~~~~~~~~~~~~~~~

.. attribute:: .guess_bonds

    :Parameter:     * **Type** - :class:`bool`
                    * **Default value** – ``False``

    Try to guess bonds and bond orders in a molecule based on the types atoms
    and the relative of atoms. Is set to False by default, with the exception
    of .xyz files.

    .. warning:
        This option works reliably only for geometries representing complete molecules.
        If some atoms are missing (for example, a protein without hydrogens) the resulting set of bonds
        would usually contain more bonds or bonds with higher order than expected.


.. attribute:: .column

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The column containing the to be imported molecules.
    Relevant when importing structures from .txt and .xlsx files with
    multiple columns.
    Relevant for .txt and .csv files.
    Numbering starts from 0.


.. attribute:: .row

    :Parameter:     * **Type** - :class:`int`
                    * **Default value** – ``0``

    The first row in a column which contains a molecule.
    Useful for when, for example, the very first row contains the title of
    aforementioned row, in which case row = 1 would be a sensible choice.
    Relevant for .txt and .csv files.
    Numbering starts from 0.

.. attribute:: .indices

    :Parameter:     * **Type** - :class:`int` or :class:`tuple` [:class:`int`]
                    * **Default value** – ``None``

    The behaviour of this argument depends on whether it is passed to a molecule
    in :attr:`input_cores` or :attr:`input_ligands`:

    .. attribute:: input_cores

        Manually specify the atomic index of one ore more atom(s) in the core that
        will be replaced with ligands. If left empty, all atoms of a user-specified
        element (see :attr:`optional.cores.dummy`) will be replaced with
        ligands.

    .. attribute:: input_ligands

        Manually specify the atomic index of the ligand atom that will be attached
        to core (implying argument_dict: :attr:`optional.ligand.split` = ``False``).
        If two atomic indices are provided (*e.g.* ``(1, 2)``), the bond between atoms ``1`` and
        [``2``] will be broken and the remaining molecule containing atom ``2`` is attached to the core,
        (implying argument_dict: :attr:`.split` = ``True``).
        Serves as an alternative to the functional group based :func:`CAT.find_substructure` function,
        which identifies the to be attached atom based on connectivity patterns
        (*i.e.* functional groups).

    .. note::
        Atom numbering follows the PLAMS [1_, 2_] convention of starting from 1 rather than 0.

.. _1: https://github.com/SCM-NV/PLAMS
.. _2: https://www.scm.com/doc/plams/index.html
nanoCAT.recipes.mol_filter
==========================

.. automodule:: nanoCAT.recipes.mol_filter
.. _PDBContainer:

The PDBContainer Class
======================

.. automodule:: dataCAT.pdb_array
nanoCAT.recipes.multi_ligand_job
================================

.. automodule:: nanoCAT.recipes.multi_ligand_job
identify_surface
================

.. automodule:: nanoCAT.bde.identify_surface

.. include:: ../README.rst
.. automodule:: nanoCAT.recipes.fast_sigma
Context Managers
================

Various context managers for manipulating molecules.

Index
-----
.. currentmodule:: CAT.attachment
.. autosummary::

    as_array.AsArray
    mol_split_cm.SplitMol
    remove_atoms_cm.RemoveAtoms

API
---
.. currentmodule:: CAT.attachment.as_array
.. autoclass:: AsArray

.. currentmodule:: CAT.attachment.mol_split_cm
.. autoclass:: SplitMol

.. currentmodule:: CAT.attachment.remove_atoms_cm
.. autoclass:: RemoveAtoms
.. _dtype:

Data Types
==========
.. automodule:: dataCAT.dtype
CAT Documentation
=================

For a more detailed description of the **CAT** compound builder read the
documentation. The documentation is divided into three parts: The basics,
further details about the input cores & ligands and finally a more detailed
look into the customization of the various jobs.

.. toctree::
    1_get_started
    2_path
    3_input_core_ligand
    4_optional
    5_bde
    6_type_aliases
    7_database
    14_pdb_container
    15_dtype
    16_hdf5_logging
    17_property_dset
    10_context_managers
    11_md_asa
    12_recipes
    13_multi_ligand
    18_distribution.rst

.. toctree::
    :hidden:

    _1_identify_surface
    _2_distribution_brute
    _4_guess_core_core_dist
    8_import_qd.rst
.. _Database:

The Database Class
==================
.. currentmodule:: dataCAT

A Class designed for the storing, retrieval and updating of results.

.. image:: _images/Database.png
    :scale: 35 %
    :align: center

The methods of the Database class can be divided into three categories
accoring to their functionality:

-   Opening & closing the database - these methods serve as context managers
    for loading and unloading parts of the database from the harddrive.

    The context managers can be accessed by calling
    either :attr:`.Database.csv_lig`, :attr:`.Database.csv_qd`,
    or :attr:`.Database.hdf5`, with the option
    of passing additional positional or keyword arguments.

    .. code:: python

        >>> from dataCAT import Database

        >>> database = Database()
        >>> with database.csv_lig(write=False) as db:
        >>>     print(repr(db))
        DFProxy(ndframe=<pandas.core.frame.DataFrame at 0x7ff8e958ce80>)

        >>> with database.hdf5('r') as db:
        >>>     print(type(db))
        <class 'h5py._hl.files.File'>

-   Importing to the database - these methods handle the importing of new data
    from python objects to the Database class:

    ============================  =============================  ================================
    :meth:`~Database.update_csv`  :meth:`~Database.update_hdf5`  :meth:`~Database.update_mongodb`
    ============================  =============================  ================================

-   Exporting from the database - these methods handle the exporting of data
    from the Database class to other python objects or remote locations:

    ==========================  ===========================
    :meth:`~Database.from_csv`  :meth:`~Database.from_hdf5`
    ==========================  ===========================


Index
-----

.. autosummary::

    Database.dirname
    Database.csv_lig
    Database.csv_qd
    Database.hdf5
    Database.mongodb

    Database.update_mongodb
    Database.update_csv
    Database.update_hdf5
    Database.from_csv
    Database.from_hdf5

    DFProxy
    OpenLig
    OpenQD


API
---
.. autoclass:: Database
    :members:

.. autoclass:: DFProxy

.. autoclass:: OpenLig
    :members:
    :inherited-members:

.. autoclass:: OpenQD
    :members:
    :inherited-members:
.. automodule:: nanoCAT.recipes.entropy
.. _Optional:

Optional
========

There are a number of arguments which can be used to modify the
functionality and behavior of the quantum dot builder. Herein an
overview is provided.

Note: Inclusion of this section in the input file is not required,
assuming one is content with the default settings.

Index
~~~~~

========================================= =========================================================================================================
Option                                    Description
========================================= =========================================================================================================
:attr:`optional.database.dirname`         The name of the directory where the database will be stored.
:attr:`optional.database.read`            Attempt to read results from the database before starting calculations.
:attr:`optional.database.write`           Export results to the database.
:attr:`optional.database.overwrite`       Allow previous results in the database to be overwritten.
:attr:`optional.database.thread_safe`     Ensure that the created workdir has a thread-safe name.
:attr:`optional.database.mol_format`      The file format(s) for exporting moleculair structures.
:attr:`optional.database.mongodb`         Options related to the MongoDB format.

:attr:`optional.core.dirname`             The name of the directory where all cores will be stored.
:attr:`optional.core.anchor`              Atomic number of symbol of the core anchor atoms.
:attr:`optional.core.allignment`          How the to-be attached ligands should be alligned with the core.
:attr:`optional.core.subset`              Settings related to the partial replacement of core anchor atoms.

:attr:`optional.ligand.dirname`           The name of the directory where all ligands will be stored.
:attr:`optional.ligand.optimize`          Optimize the geometry of the to-be attached ligands.
:attr:`optional.ligand.anchor`            Manually specify SMILES strings representing functional groups.
:attr:`optional.ligand.split`             If the ligand should be attached in its entirety to the core or not.
:attr:`optional.ligand.cosmo-rs`          Perform a property calculation with COSMO-RS on the ligand.
:attr:`optional.ligand.cdft`              Perform a conceptual DFT calculation with ADF on the ligand.

:attr:`optional.qd.dirname`               The name of the directory where all quantum dots will be stored.
:attr:`optional.qd.construct_qd`          Whether or not the quantum dot should actually be constructed or not.
:attr:`optional.qd.optimize`              Optimize the quantum dot (i.e. core + all ligands).
:attr:`optional.qd.multi_ligand`          A workflow for attaching multiple non-unique ligands to a single quantum dot.
:attr:`optional.qd.bulkiness`             Calculate the :math:`V_{bulk}`, a ligand- and core-sepcific descriptor of a ligands' bulkiness.
:attr:`optional.qd.activation_strain`     Perform an activation strain analyses.
:attr:`optional.qd.dissociate`            Calculate the ligand dissociation energy.
========================================= =========================================================================================================

Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    optional:
        database:
            dirname: database
            read: True
            write: True
            overwrite: False
            thread_safe: False
            mol_format: (pdb, xyz)
            mongodb: False

        core:
            dirname: core
            anchor: Cl
            allignment: surface
            subset: null

        ligand:
            dirname: ligand
            optimize: True
            anchor: null
            split: True
            cosmo-rs: False
            cdft: False

        qd:
            dirname: qd
            construct_qd: True
            optimize: False
            activation_strain: False
            dissociate: False
            bulkiness: False

Arguments
~~~~~~~~~

Database
--------

.. attribute:: optional.database

    All database-related settings.

    .. note::
        For :attr:`optional.database` settings to take effect the `Data-CAT <https://github.com/nlesc-nano/data-CAT>`_ package has to be installed.

    Example:

    .. code:: yaml

        optional:
            database:
                dirname: database
                read: True
                write: True
                overwrite: False
                mol_format: (pdb, xyz)
                mongodb: False

|

    .. attribute:: optional.database.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default Value** - ``"database"``

        The name of the directory where the database will be stored.

        The database directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.database.read

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Attempt to read results from the database before starting calculations.

        Before optimizing a structure, check if a geometry is available from
        previous calculations. If a match is found, use that structure and
        avoid any geometry (re-)optimizations. If one wants more control then the
        boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures will be read only for a
        specific subset.


        .. admonition:: Example

            Example #1:

            .. code:: yaml

                optional:
                    database:
                        read: (core, ligand, qd)  # This is equivalent to read: True

            Example #2:

            .. code:: yaml

                optional:
                    database:
                        read: ligand


    .. attribute:: optional.database.write

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("core", "ligand", "qd")``

        Export results to the database.

        Previous results will **not** be overwritten unless
        :attr:`optional.database.overwrite` = ``True``. If one wants more control then
        the boolean can be substituted for a list of strings (*i.e.* ``"core"``,
        ``"ligand"`` and/or ``"qd"``), meaning that structures written for a specific
        subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.overwrite

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``False``

        Allow previous results in the database to be overwritten.

        Only applicable if :attr:`optional.database.write` = ``True``.
        If one wants more control then the boolean can be substituted for
        a list of strings (*i.e.* ``"core"``, ``"ligand"`` and/or ``"qd"``), meaning
        that structures written for a specific subset.

        See :attr:`optional.database.read` for a similar relevant example.


    .. attribute:: optional.database.thread_safe

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** - ``False``

        Ensure that the created workdir has a thread-safe name.

        Note that this disables the restarting of partially completed jobs.


    .. attribute:: optional.database.mol_format

        :Parameter:     * **Type** - :class:`bool`, :class:`str` or :class:`tuple` [:class:`str`]
                        * **Default value** - ``("pdb", "xyz")``

        The file format(s) for exporting moleculair structures.

        By default all structures are stored in the .hdf5 format as
        (partially) de-serialized .pdb files. Additional formats can be
        requested with this keyword.
        Accepted values: ``"pdb"``, ``"xyz"``, ``"mol"`` and/or ``"mol2"``.


    .. attribute:: optional.database.mongodb

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default Value** – ``False``

        Options related to the MongoDB format.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Database`:.

|

Core
----

.. attribute:: optional.core

    All settings related to the core.

    Example:

    .. code:: yaml

        optional:
            core:
                dirname: core
                anchor: Cl
                allignment: surface
                subset: null

|

    .. attribute:: optional.core.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"core"``

        The name of the directory where all cores will be stored.

        The core directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.core.anchor

        :Parameter:     * **Type** - :class:`str` or :class:`int`
                        * **Default value** – ``17``

        Atomic number of symbol of the core anchor atoms.

        The atomic number or atomic symbol of the atoms in the core which are to be
        replaced with ligands. Alternatively, anchor atoms can be manually specified
        with the core_indices variable.

        Further customization can be achieved by passing a dictionary:

        * :attr:`anchor.group <optional.ligand.anchor.group>`
        * :attr:`anchor.group_idx <optional.ligand.anchor.group_idx>`
        * :attr:`anchor.group_format <optional.ligand.anchor.group_format>`
        * :attr:`anchor.remove <optional.ligand.anchor.remove>`

        .. note::

            .. code:: yaml

                optional:
                    core:
                        anchor:
                            group: "[H]Cl"  # Remove HCl and attach at previous Cl position
                            group_idx: 1
                            group_format: "SMILES"
                            remove: [0, 1]



    .. attribute:: optional.core.allignment

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"surface"``

        How the to-be attached ligands should be alligned with the core.

        Has four allowed values:

        * ``"surface"``: Define the core vectors as those orthogonal to the cores
          surface. Not this option requires at least four core anchor atoms.
          The surface is herein defined by a convex hull constructed from the core.
        * ``"sphere"``: Define the core vectors as those drawn from the core anchor
          atoms to the cores center.
        * ``"surface invert"``/``"surface_invert"``: The same as ``"surface"``,
          except the core vectors are inverted.
        * ``"sphere invert"``/``"sphere_invert"``: The same as ``"sphere"``,
          except the core vectors are inverted.

        Note that for a spherical core both approaches are equivalent.

        .. note::
            An example of a ``"sphere"`` (left) and ``"surface"`` (right) allignment.

            .. image:: _images/allignment.png
                :scale: 15 %
                :align: center



    .. attribute:: optional.core.subset

        :Parameter:     * **Type** - :class:`dict`, optional
                        * **Default value** – ``None``

        Settings related to the partial replacement of core anchor atoms with ligands.

        If not ``None``, has access to six further keywords,
        the first two being the most important:

        * :attr:`subset.f`
        * :attr:`subset.mode`
        * :attr:`subset.follow_edge`
        * :attr:`subset.weight`
        * :attr:`subset.randomness`
        * :attr:`subset.cluster_size`


    .. attribute:: optional.core.subset.f

        :Parameter:     * **Type** - :class:`float`

        The fraction of core anchor atoms that will actually be exchanged for ligands.

        The provided value should satisfy the following condition: :math:`0 < f \le 1`.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.core.subset.mode

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"uniform"``

        Defines how the anchor atom subset, whose size is defined by the fraction :math:`f`, will be generated.

        Accepts one of the following values:

        * ``"uniform"``: A uniform distribution; the nearest-neighbor distances between each
          successive anchor atom and all previous anchor atoms is maximized.
          can be combined with :attr:`subset.cluster_size<optional.core.subset.cluster_size>`
          to create a uniform distribution of clusters of a user-specified size.
        * ``"cluster"``: A clustered distribution; the nearest-neighbor distances between each
          successive anchor atom and all previous anchor atoms is minimized.
        * ``"random"``: A random distribution.

        It should be noted that all three methods converge towards the same set
        as :math:`f` approaches :math:`1.0`.

        If :math:`\boldsymbol{D} \in \mathbb{R}_{+}^{n,n}` is the (symmetric) distance matrix constructed
        from the anchor atom superset and :math:`\boldsymbol{a} \in \mathbb{N}^{m}` is the vector
        of indices which yields the anchor atom subset. The definition of element :math:`a_{i}`
        is defined below for the ``"uniform"`` distribution.
        All elements of :math:`\boldsymbol{a}` are furthermore constrained to be unique.

        .. math::
            :label: 1

            \DeclareMathOperator*{\argmin}{\arg\!\min}
            a_{i} = \begin{cases}
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D_{k, \hat{\imath}} \right) &
                \text{if} & i=0 \\
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{i-1} f \left( D[k, a_{\hat{\imath}}]\ \right) &
                \text{if} & i > 0
            \end{cases} \begin{matrix} & \text{with} & f(x) = e^{-x} \end{matrix}

        For the ``"cluster"`` distribution all :math:`\text{argmin}` operations
        are exchanged for :math:`\text{argmax}`.

        The old default, the p-norm with :math:`p=-2`, is equivalent to:

        .. math::
            :label: 2

            \DeclareMathOperator*{\argmax}{\arg\!\max}
            \begin{matrix}
            \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D_{k, \hat{\imath}} \right) =
            \argmax\limits_{k \in \mathbb{N}} \left( \sum_{\hat{\imath}=0}^{n} | D_{k, \hat{\imath}} |^p \right)^{1/p}
            & \text{if} & f(x) = x^{-2} \end{matrix}

        Note that as the elements of :math:`\boldsymbol{D}` were defined as positive or zero-valued real numbers;
        operating on :math:`\boldsymbol{D}` is thus equivalent to operating on its absolute.

        .. note::
            An example of a ``"uniform"``, ``"cluster"`` and ``"random"`` distribution with :math:`f=1/3`.

            .. image:: _images/distribution.png
                :scale: 15 %
                :align: center

        .. note::
            An example of four different ``"uniform"`` distributions at :math:`f=1/16`,
            :math:`f=1/8`, :math:`f=1/4` and :math:`f=1/2`.

            .. image:: _images/distribution_p_var.png
                :scale: 20 %
                :align: center


    .. attribute:: optional.core.subset.follow_edge

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Construct the anchor atom distance matrix by following the shortest path along the
        edges of a (triangular-faced) polyhedral approximation of the core rather than the
        shortest path through space.

        Enabling this option will result in more accurate ``"uniform"`` and ``"cluster"``
        distributions at the cost of increased computational time.

        Given the matrix of Cartesian coordinates :math:`\boldsymbol{X} \in \mathbb{R}^{n, 3}`,
        the matching edge-distance matrix :math:`\boldsymbol{D}^{\text{edge}} \in \mathbb{R}_{+}^{n, n}`
        and the vector :math:`\boldsymbol{p} \in \mathbb{N}^{m}`, representing a (to-be optimized)
        path as the indices of edge-connected vertices, then element :math:`D_{i,j}^{\text{edge}}`
        is defined as following:

        .. math::
            :label: 3

            D_{i, j}^{\text{edge}} = \min_{\boldsymbol{p} \in \mathbb{N}^{m}; m \in \mathbb{N}}
            \sum_{k=0}^{m-1} || X_{p_{k},:} - X_{p_{k+1},:} ||
            \quad \text{with} \quad p_{0} = i \quad \text{and} \quad p_{m} = j

        The polyhedron edges are constructed, after projecting all vertices on the surface of a sphere,
        using Qhull's :class:`ConvexHull<scipy.spatial.ConvexHull>` algorithm
        (`The Quickhull Algorithm for Convex Hulls <https://doi.org/10.1145/235815.235821>`_).
        The quality of the constructed edges is proportional to the convexness of the core,
        more specifically: how well the vertices can be projected on a spherical surface without
        severely distorting the initial structure.
        For example, spherical, cylindrical or cuboid cores will yield reasonably edges,
        while the edges resulting from torus will be extremely poor.

        .. note::
            An example of a cores' polyhedron-representation; displaying the shortest path
            between points :math:`i` and :math:`j`.

            .. image:: _images/polyhedron.png
                :scale: 15 %
                :align: center


    .. attribute:: optional.core.subset.cluster_size

        :Parameter:     * **Type** - :class:`int` or :class:`Iterable<collections.abc.Iterable>` [:class:`int`]
                        * **Default value** – ``1``

        Allow for the creation of uniformly distributed clusters of size :math:`r`;
        should be used in conjunction with :attr:`subset.mode = "uniform"<optional.core.subset.mode>`.

        The value of :math:`r` can be either
        a single cluster size (*e.g.* :code:`cluster_size = 5`) or an iterable of various
        sizes (*e.g.* :code:`cluster_size = [2, 3, 4]`).
        In the latter case the iterable will be repeated as long as necessary.

        Compared to Eq :eq:`2` the vector of indices :math:`\boldsymbol{a} \in \mathbb{N}^{m}` is,
        for the purpose of book keeping, reshaped into the matrix
        :math:`\boldsymbol{A} \in \mathbb{N}^{q, r} \; \text{with} \; q*r = m`.
        All elements of :math:`\boldsymbol{A}` are, again, constrained to be unique.

        .. math::
            :label: 4

            \DeclareMathOperator*{\argmin}{\arg\!\min}
            A_{i,j} = \begin{cases}
                \argmin\limits_{k \in \mathbb{N}} \sum_{\hat{\imath}=0}^{n} f \left( D[k, \, \hat{\imath}] \right) &
                \text{if} & i=0 & \text{and} & j=0 \\
            \argmin\limits_{k \in \mathbb{N}}
                \sum_{\hat{\imath}=0}^{i-1} \sum_{\hat{\jmath}=0}^{r} f \left( D[k, A_{\hat{\imath}, \, \hat{\jmath}}] \right) &
            \text{if} & i > 0 & \text{and} & j = 0 \\
            \argmin\limits_{k \in \mathbb{N}}
            \dfrac
                { \sum_{\hat{\imath}=0}^{i-1} \sum_{\hat{\jmath}=0}^{r} f \left( D[k, A_{\hat{\imath}, \, \hat{\jmath}}] \right) }
                { \sum_{\hat{\jmath}=0}^{j-1} f \left( D[k, A_{i, \, \hat{\jmath}}] \right) }
            &&& \text{if} & j > 0
            \end{cases}

        |

        .. note::
            An example of various cluster sizes (1, 2, 3 and 4) with :math:`f=1/4`.

            .. image:: _images/cluster_size.png
                :scale: 15 %
                :align: center

        .. note::
            An example of clusters of varying size (:code:`cluster_size = [1, 2, 9, 1]`)
            with :math:`f=1/4`.

            .. image:: _images/cluster_size_variable.png
                :scale: 5 %
                :align: center


    .. attribute:: optional.core.subset.weight

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"numpy.exp(-x)"``

        The function :math:`f(x)` for weighting the distance.; its default value corresponds to: :math:`f(x) = e^{-x}`.

        For the old default, the p-norm with :math:`p=-2`, one can use ``weight = "x**-2"``: :math:`f(x) = x^-2`.

        Custom functions can be specified as long as they satisfy the following constraints:

        * The function must act an variable by the name of ``x``,
          a 2D array of positive and/or zero-valued floats (:math:`x \in \mathbb{R}_{+}^{n, n}`).
        * The function must take a single array as argument and return a new one.
        * The function must be able to handle values of ``numpy.nan`` and ``numpy.inf`` without
          raising exceptions.
        * The shape and data type of the output array should not change with respect to the input.

        Modules specified in the weight function will be imported when required,
        illustrated here with SciPy's :func:`expit<scipy.special.expit>`
        function: ``weight = "scipy.special.expit(x)"`` aka ``weight = "1 / (1 + numpy.exp(-x))"``

        Multi-line statements are allowed: ``weight = "a = x**2; b = 5 * a; numpy.exp(b)"``.
        The last part of the statement is assumed to be the to-be returned value
        (*i.e.* ``return numpy.exp(b)``).


    .. attribute:: optional.core.subset.randomness

        :Parameter:     * **Type** - :class:`float`, optional
                        * **Default value** – ``None``

        The probability that each new core anchor atom will be picked at random.

        Can be used in combination with ``"uniform"`` and ``"cluster"`` to introduce
        a certain degree of randomness (*i.e.* entropy).

        If not ``None``, the provided value should satisfy the following condition:
        :math:`0 \le randomness \le 1`. A value of :math:`0` is equivalent to a
        ``"uniform"`` / ``"cluster"`` distribution while :math:`1` is equivalent
        to ``"random"``.

        .. note::
            A demonstration of the ``randomness`` parameter for a ``"uniform"`` and
            ``"cluster"`` distribution at :math:`f=1/4`.

            The ``randomness`` values are (from left to right) set to :math:`0`,
            :math:`1/4`, :math:`1/2` and :math:`1`.

            .. image:: _images/randomness.png
                :scale: 13 %
                :align: center

|

Ligand
------

.. attribute:: optional.ligand

    All settings related to the ligands.

    Example:

    .. code:: yaml

        optional:
            ligand:
                dirname: ligand
                optimize: True
                anchor: null
                split: True
                cosmo-rs: False
                cdft: False

|

    .. attribute:: optional.ligand.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"ligand"``

        The name of the directory where all ligands will be stored.

        The ligand directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.


    .. attribute:: optional.ligand.optimize

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``True``

        Optimize the geometry of the to-be attached ligands.

        The ligand is split into one or multiple (more or less) linear fragments,
        which are subsequently optimized (RDKit UFF [1_, 2_, 3_]) and reassembled
        while checking for the optimal dihedral angle. The ligand fragments are
        biased towards more linear conformations to minimize inter-ligand
        repulsion once the ligands are attached to the core.

        After the conformation search a final (unconstrained) geometry optimization
        is performed, RDKit UFF again being the default level of theory.
        Custom job types and settings can, respectivelly, be specified with the
        ``job2`` and ``s2`` keys.

        .. note::

            .. code:: yaml

                optional:
                    ligand:
                        optimize:
                            job2: ADFJob


    .. attribute:: optional.ligand.anchor

        :Parameter:     * **Type** - :class:`str`, :class:`Sequence[str] <collections.abc.Sequence>` or :class:`dict[str, Any] <dict>`
                        * **Default value** – ``None``

        Manually specify SMILES strings representing functional groups.

        For example, with :attr:`optional.ligand.anchor` = ``("O[H]", "[N+].[Cl-]")`` all
        ligands will be searched for the presence of hydroxides and ammonium chlorides.

        The first atom in each SMILES string (*i.e.* the "anchor") will be used for attaching the ligand
        to the core, while the last atom (assuming :attr:`optional.ligand.split` = ``True``) will be
        dissociated from the ligand and discarded.

        If not specified, the default functional groups of **CAT** are used.

        This option can alternatively be provided as ``optional.ligand.functional_groups``.

        Further customization can be achieved by passing dictionaries:

        * :attr:`anchor.group`
        * :attr:`anchor.group_idx`
        * :attr:`anchor.group_format`
        * :attr:`anchor.remove`
        * :attr:`anchor.kind`
        * :attr:`anchor.angle_offset`
        * :attr:`anchor.dihedral`
        * :attr:`anchor.multi_anchor_filter`

        .. note::

            .. code:: yaml

                optional:
                    ligand:
                        anchor:
                            - group: "[H]OC(=O)C"  # Remove H and attach at the (formal) oxyanion
                              group_idx: 1
                              remove: 0
                            - group: "[H]OC(=O)C"  # Remove H and attach at the mean position of both oxygens
                              group_idx: [1, 3]
                              remove: 0
                              kind: mean

        .. note::
            This argument has no value be default and will thus default to SMILES strings of the default
            functional groups supported by **CAT**.

        .. note::
            The yaml format uses ``null`` rather than ``None`` as in Python.


    .. attribute:: optional.ligand.anchor.group

        :Parameter:     * **Type** - :class:`str`

        A SMILES string representing the anchoring group.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.ligand.anchor.group_idx

        :Parameter:     * **Type** - :class:`int` or :class:`Sequence[int] <collections.abc.Sequence>`

        The indices of the anchoring atom(s) in :attr:`anchor.group <optional.ligand.anchor.group>`.

        Indices should be 0-based.
        These atoms will be attached to the core, the manner in which is determined by the :attr:`anchor.kind` option.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.ligand.anchor.group_format

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – :data:`"SMILES"`

        The format used for representing :attr:`anchor.group <optional.ligand.anchor.group>`.

        Defaults to the SMILES format.
        The supported formats (and matching RDKit parsers) are as following:

        .. code-block:: python

            >>> import rdkit.Chem

            >>> FASTA      = rdkit.Chem.MolFromFASTA
            >>> HELM       = rdkit.Chem.MolFromHELM
            >>> INCHI      = rdkit.Chem.MolFromInchi
            >>> MOL2       = rdkit.Chem.MolFromMol2Block
            >>> MOL2_FILE  = rdkit.Chem.MolFromMol2File
            >>> MOL        = rdkit.Chem.MolFromMolBlock
            >>> MOL_FILE   = rdkit.Chem.MolFromMolFile
            >>> PDB        = rdkit.Chem.MolFromPDBBlock
            >>> PDB_FILE   = rdkit.Chem.MolFromPDBFile
            >>> PNG        = rdkit.Chem.MolFromPNGString
            >>> PNG_FILE   = rdkit.Chem.MolFromPNGFile
            >>> SVG        = rdkit.Chem.MolFromRDKitSVG
            >>> SEQUENCE   = rdkit.Chem.MolFromSequence
            >>> SMARTS     = rdkit.Chem.MolFromSmarts
            >>> SMILES     = rdkit.Chem.MolFromSmiles
            >>> TPL        = rdkit.Chem.MolFromTPLBlock
            >>> TPL_FILE   = rdkit.Chem.MolFromTPLFile


    .. attribute:: optional.ligand.anchor.remove

        :Parameter:     * **Type** - :data:`None`, :class:`int` or :class:`Sequence[int] <collections.abc.Sequence>`
                        * **Default value** – :data:`None`

        The indices of the to-be removed atoms in :attr:`anchor.group <optional.ligand.anchor.group>`.

        No atoms are removed when set to :data:`None`.
        Indices should be 0-based.
        See also the :attr:`~optional.ligand.split` option.


    .. attribute:: optional.ligand.anchor.kind

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"first"``

        How atoms are to-be attached when multiple anchor atoms are specified in :attr:`anchor.group_idx <optional.ligand.anchor.group_idx>`.

        Accepts one of the following options:

        * ``"first"``: Attach the first atom to the core.
        * ``"mean"``: Attach the mean position of all anchoring atoms to the core.
        * ``"mean_translate"``: Attach the mean position of all anchoring atoms to the core and then translate back to the first atom.


    .. attribute:: optional.ligand.anchor.angle_offset

        :Parameter:     * **Type** - :data:`None`, :class:`float` or :class:`str`
                        * **Default value** – :data:`None`

        Manually offset the angle of the ligand vector by a given number.

        The plane of rotation is defined by the first three indices in :attr:`anchor.group_idx <optional.ligand.anchor.group_idx>`.

        By default the angle unit is assumed to be in degrees,
        but if so desired one can explicitly pass the unit: ``angle_offset: "0.25 rad"``.


    .. attribute:: optional.ligand.anchor.dihedral

        :Parameter:     * **Type** - :data:`None`, :class:`float` or :class:`str`
                        * **Default value** – :data:`None`

        Manually specify the ligands vector dihedral angle, rather than optimizing it w.r.t. the inter-ligand distance.

        The dihedral angle is defined by three vectors:

        * The first two in dices in :attr:`anchor.group_idx <optional.ligand.anchor.group_idx>`.
        * The core vector(s).
        * The Cartesian X-axis as defined by the core.

        By default the angle unit is assumed to be in degrees,
        but if so desired one can explicitly pass the unit: ``dihedral: "0.5 rad"``.


    .. attribute:: optional.ligand.anchor.multi_anchor_filter

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – :data:`"ALL"`

        How ligands with multiple valid anchor sites are to-be treated.

        Accepts one of the following options:

        * ``"all"``: Construct a new ligand for each valid anchor/ligand combination.
        * ``"first"``: Pick only the first valid functional group, all others are ignored.
        * ``"raise"``: Treat a ligand as invalid if it has multiple valid anchoring sites.


    .. attribute:: optional.ligand.split

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        If ``False``: The ligand is to be attached to the core in its entirety .

        =================== ==================
        Before              After
        =================== ==================
        :math:`{NR_4}^+`    :math:`{NR_4}^+`
        :math:`O_2 CR`      :math:`O_2 CR`
        :math:`HO_2 CR`     :math:`HO_2 CR`
        :math:`H_3 CO_2 CR` :math:`H_3 CO_2 CR`
        =================== ==================

        ``True``: A proton, counterion or functional group is to be removed from
        the ligand before attachment to the core.

        ========================= ==================
        Before                    After
        ========================= ==================
        :math:`Cl^- + {NR_4}^+`   :math:`{NR_4}^+`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`Na^+ + {O_2 CR}^-` :math:`{O_2 CR}^-`
        :math:`HO_2 CR`           :math:`{O_2 CR}^-`
        :math:`H_3 CO_2 CR`       :math:`{O_2 CR}^-`
        ========================= ==================


    .. attribute:: optional.ligand.cosmo-rs

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``


        Perform a property calculation with COSMO-RS [4_, 5_, 6_, 7_] on the ligand.

        The COSMO surfaces are by default constructed using ADF MOPAC [8_, 9_, 10_].

        The solvation energy of the ligand and its activity coefficient are
        calculated in the following solvents: acetone, acetonitrile,
        dimethyl formamide (DMF), dimethyl sulfoxide (DMSO), ethyl acetate,
        ethanol, *n*-hexane, toluene and water.


    .. attribute:: optional.ligand.cdft

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``


        Perform a conceptual DFT (CDFT) calculation with `ADF <https://www.scm.com/doc/ADF/Input/Advanced_analysis.html#conceptual-dft>`_ on the ligand.

        All global descriptors are, if installed, stored in the database.
        This includes the following properties:

        * Electronic chemical potential (mu)
        * Electronic chemical potential (mu+)
        * Electronic chemical potential (mu-)
        * Electronegativity (chi=-mu)
        * Hardness (eta)
        * Softness (S)
        * Hyperhardness (gamma)
        * Electrophilicity index (w=omega)
        * Dissocation energy (nucleofuge)
        * Dissociation energy (electrofuge)
        * Electrodonating power (w-)
        * Electroaccepting power(w+)
        * Net Electrophilicity
        * Global Dual Descriptor Deltaf+
        * Global Dual Descriptor Deltaf-

        This block can be furthermore customized with one or more of the following keys:

        * ``"keep_files"``: Whether or not to delete the ADF output afterwards.
        * ``"job1"``: The type of PLAMS Job used for running the calculation.
          The only value that should be supplied here (if any) is ``"ADFJob"``.
        * ``"s1"``: The job Settings used for running the CDFT calculation.
          Can be left blank to use the default template (:data:`nanoCAT.cdft.cdft`).

        .. admonition:: Examples

            .. code:: yaml

                optional:
                    ligand:
                        cdft: True

            .. code:: yaml

                optional:
                    ligand:
                        cdft:
                            job1: ADFJob
                            s1: ...  # Insert custom settings here

|

QD
--

.. attribute:: optional.qd

    All settings related to the quantum dots.

    Example:

    .. code:: yaml

        optional:
            qd:
                dirname: qd
                construct_qd: True
                optimize: False
                bulkiness: False
                activation_strain: False
                dissociate: False

|

    .. attribute:: optional.qd.dirname

        :Parameter:     * **Type** - :class:`str`
                        * **Default value** – ``"qd"``

        The name of the directory where all quantum dots will be stored.

        The quantum dot directory will be created (if it does not yet exist)
        at the path specified in :ref:`Path`.

    .. attribute:: optional.qd.construct_qd

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Whether or not the quantum dot should actually be constructed or not.

        Setting this to ``False`` will still construct ligands and carry out ligand workflows,
        but it will not construct the actual quantum dot itself.


    .. attribute:: optional.qd.optimize

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Optimize the quantum dot (i.e. core + all ligands) .

        By default the calculation is performed with ADF UFF [3_, 11_].
        The geometry of the core and ligand atoms directly attached to the core
        are frozen during this optimization.


    .. attribute:: optional.qd.multi_ligand

        :Parameter:     * **Type** - ``None`` or :class:`dict`
                        * **Default value** – ``None``

        A workflow for attaching multiple non-unique ligands to a single quantum dot.

        Note that this is considered a seperate workflow besides the normal ligand attachment.
        Consequently, these structures will *not* be passed to further workflows.

        See :ref:`Multi-ligand` for more details regarding the available options.

        .. note::
            An example with ``[O-]CCCC`` as main ligand and
            ``[O-]CCCCCCCCCCCCC`` & ``[O-]C`` as additional ligands.

            .. image:: _images/multi_ligand.png
                :scale: 13 %
                :align: center


    .. attribute:: optional.qd.bulkiness

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Calculate the :math:`V_{bulk}`, a ligand- and core-specific descriptor of a ligands' bulkiness.

        Supplying a dictionary grants access to the two additional :attr:`~optional.qd.bulkiness.h_lim`
        and :attr:`~optional.qd.bulkiness.d` sub-keys.

        .. math::
            :label: 5

            V(r_{i}, h_{i}; d, h_{lim}) =
            \sum_{i=1}^{n} e^{r_{i}} (\frac{2 r_{i}}{d} - 1)^{+} (1 - \frac{h_{i}}{h_{lim}})^{+}


    .. attribute:: optional.qd.bulkiness.h_lim

        :Parameter:     * **Type** - :class:`float` or :data:`None`
                        * **Default value** – ``10.0``

        Default value of the :math:`h_{lim}` parameter in :attr:`~optional.qd.bulkiness`.

        Set to :data:`None` to disable the :math:`h_{lim}`-based cutoff.


    .. attribute:: optional.qd.bulkiness.d

        :Parameter:     * **Type** - :class:`float`/:class:`list[float] <list>`, :data:`None` or ``"auto"``
                        * **Default value** – ``"auto"``

        Default value of the :math:`d` parameter in :attr:`~optional.qd.bulkiness`.

        Set to ``"auto"`` to automatically infer this parameters value based on the mean
        nearest-neighbor distance among the core anchor atoms.
        Set to :data:`None` to disable the :math:`d`-based cutoff.
        Supplying multiple floats will compute the bulkiness for all specified values.


    .. attribute:: optional.qd.activation_strain

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Perform an activation strain analysis [12_, 13_, 14_].

        The activation strain analysis (kcal mol\ :sup:`-1`\) is performed
        on the ligands attached to the quantum dot surface with RDKit UFF [1_, 2_, 3_].

        The core is removed during this process; the analysis is thus exclusively
        focused on ligand deformation and inter-ligand interaction.
        Yields three terms:

        1.  d\ *E*\ :sub:`strain`\  : 	The energy required to deform the ligand
        from their equilibrium geometry to the geometry they adopt on the quantum
        dot surface. This term is, by definition, destabilizing. Also known as the
        preparation energy (d\ *E*\ :sub:`prep`\).

        2.  d\ *E*\ :sub:`int`\  :	The mutual interaction between all deformed
        ligands. This term is characterized by the non-covalent interaction between
        ligands (UFF Lennard-Jones potential) and, depending on the inter-ligand
        distances, can be either stabilizing or destabilizing.

        3.  d\ *E* :	The sum of d\ *E*\ :sub:`strain`\  and d\ *E*\ :sub:`int`\ .
        Accounts for both the destabilizing ligand deformation and (de-)stabilizing
        interaction between all ligands in the absence of the core.

        See :ref:`md_asa` for more details.


    .. attribute:: optional.qd.dissociate

        :Parameter:     * **Type** - :class:`bool` or :class:`dict`
                        * **Default value** – ``False``

        Calculate the ligand dissociation energy.

        Calculate the ligand dissociation energy (BDE) of ligands attached to the
        surface of the core. See :ref:`Bond Dissociation Energy` for more details.
        The calculation consists of five distinct steps:

            1.  Dissociate all combinations of |n| ligands (|Y|) and an atom from the core (|X|)
            within a radius *r* from aforementioned core atom.
            The dissociated compound has the general structure of |XYn|.

            2.  Optimize the geometry of |XYn| at the first level of theory
            (:math:`1`). Default: ADF MOPAC [1_, 2_, 3_].

            3.  Calculate the "electronic" contribution to the BDE (|dE|)
            at the first level of theory (:math:`1`): ADF MOPAC [1_, 2_, 3_].
            This step consists of single point calculations of the complete
            quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

            4.  Calculate the thermochemical contribution to the BDE (|ddG|) at the
            second level of theory (:math:`2`). Default: ADF UFF [4_, 5_]. This step
            consists of geometry optimizations and frequency analyses of the same
            compounds used for step 3.

            5.  :math:`\Delta G_{tot} = \Delta E_{1} + \Delta \Delta G_{2} = \Delta E_{1} + (\Delta G_{2} - \Delta E_{2})`.

        .. admonition:: See also

            More extensive options for this argument are provided in :ref:`Bond Dissociation Energy`:.



.. _1: http://www.rdkit.org
.. _2: https://github.com/rdkit/rdkit
.. _3: https://doi.org/10.1021/ja00051a040
.. _4: https://www.scm.com/doc/COSMO-RS/index.html
.. _5: https://doi.org/10.1021/j100007a062
.. _6: https://doi.org/10.1021/jp980017s
.. _7: https://doi.org/10.1139/V09-008
.. _8: https://www.scm.com/doc/MOPAC/Introduction.html
.. _9: http://openmopac.net
.. _10: https://doi.org/10.1007/s00894-012-1667-x
.. _11: https://www.scm.com/doc/UFF/index.html
.. _12: https://doi.org/10.1002/9780470125922.ch1
.. _13: https://doi.org/10.1002/wcms.1221
.. _14: https://doi.org/10.1021/acs.jpcc.5b02987

.. |dE| replace:: :math:`\Delta E`
.. |dE_lvl1| replace:: :math:`\Delta E_{1}`
.. |dE_lvl2| replace:: :math:`\Delta E_{2}`
.. |dG| replace:: :math:`\Delta G_{tot}`
.. |dG_lvl2| replace:: :math:`\Delta G_{2}`
.. |ddG| replace:: :math:`\Delta \Delta G`
.. |ddG_lvl2| replace:: :math:`\Delta \Delta G_{2}`
.. |XYn| replace:: :math:`XY_{n}`
.. |Yn| replace:: :math:`Y_{n}`
.. |n| replace:: :math:`{n}`
.. |X| replace:: :math:`X`
.. |Y| replace:: :math:`Y`
.. automodule:: CAT.dye.addlig
distribution_brute
******************

.. automodule:: CAT.attachment.distribution_brute
.. _Type Aliases:

Type Aliases
============

Aliases are available for a large number of job types,
allowing one to pass a :class:`str` instead of a :class:`type` object, thus simplifying
the input settings for **CAT**. Aliases are insensitive towards capitalization
(or lack thereof).

A comprehensive list of :class:`plams.Job<scm.plams.core.basejob.Job>` subclasses and their respective
aliases (*i.e.* :class:`str`) is presented below.

Aliases
~~~~~~~

-   |ADFJob| = ``"adf"`` = ``"adfjob"``

-   |AMSJob| = ``"ams"`` = ``"amsjob"``

-   |UFFJob| = ``"uff"`` = ``"uffjob"``

-   |BANDJob| = ``"band"`` = ``"bandjob"``

-   |DFTBJob| = ``"dftb"`` = ``"dftbjob"``

-   |MOPACJob| = ``"mopac"`` = ``"mopacjob"``

-   |ReaxFFJob| = ``"reaxff"`` = ``"reaxffjob"``

-   |Cp2kJob| = ``"cp2k"`` = ``"cp2kjob"``

-   |ORCAJob| = ``"orca"`` = ``"orcajob"``

-   |DiracJob| = ``"dirac"`` = ``"diracjob"``

-   |GamessJob| = ``"gamess"`` = ``"gamessjob"``

-   |DFTBPlusJob| = ``"dftbplus"`` = ``"dftbplusjob"``

-   |CRSJob| = ``"crs"`` = ``"cosmo-rs"`` = ``"crsjob"``


.. |ADFJob| replace:: :class:`ADFJob<scm.plams.interfaces.adfsuite.adf.ADFJob>`
.. |AMSJob| replace:: :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`
.. |UFFJob| replace:: :class:`UFFJob<scm.plams.interfaces.adfsuite.uff.UFFJob>`
.. |BANDJob| replace:: :class:`BANDJob<scm.plams.interfaces.adfsuite.band.BANDJob>`
.. |DFTBJob| replace:: :class:`DFTBJob<scm.plams.interfaces.adfsuite.dftb.DFTBJob>`
.. |MOPACJob| replace:: :class:`MOPACJob<scm.plams.interfaces.adfsuite.mopac.MOPACJob>`
.. |ReaxFFJob| replace:: :class:`ReaxFFJob<scm.plams.interfaces.adfsuite.reaxff.ReaxFFJob>`
.. |Cp2kJob| replace:: :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>`
.. |ORCAJob| replace:: :class:`ORCAJob<scm.plams.interfaces.thirdparty.orca.ORCAJob>`
.. |DiracJob| replace:: :class:`DiracJob<scm.plams.interfaces.thirdparty.dirac.DiracJob>`
.. |GamessJob| replace:: :class:`GamessJob<scm.plams.interfaces.thirdparty.gamess.GamessJob>`
.. |DFTBPlusJob| replace:: :class:`DFTBPlusJob<scm.plams.interfaces.thirdparty.dftbplus.DFTBPlusJob>`
.. |CRSJob| replace:: :class:`CRSJob<scm.plams.interfaces.adfsuite.crs.CRSJob>`
.. automodule:: nanoCAT.recipes.dissociation
.. _md_asa:

Ensemble-Averaged Activation Strain Analysis
============================================
.. math::
    :label: 1a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}

Herein we describe an Ensemble-Averaged extension of the
activation/strain  analysis (ASA; also known as the
distortion/interaction model), wherein the ASA is utilized
for the analyses of entire molecular dynamics trajectories.
The implementation utilizes CHARMM-style forcefields for the
calculation of all energy terms.

.. note::
    Throughout this document an overline will be used to distinguish between "normal"
    and ensemble-averaged quantities: *e.g.* :math:`E_{\text{strain}}` versus
    :math:`\overline{E}_{\text{strain}}`.

|

Strain/Distortion
-----------------
The ensemble averaged strain :math:`\Delta \overline{E}_{\text{strain}}`
represents the distortion of all ligands with respect to their equilibrium
geometry.
Given an MD trajectory with :math:`m` iterations and :math:`n` ligands per
quantum dot, the energy is averaged over all :math:`m` MD iterations and
summed over all :math:`n` ligands.

The magnitude of this term is determined by all covalent and non-covalent
intra-ligand interactions.
As this term quantifies the deviation of a ligand from its equilibrium geometry,
it is, by definition, always positive.

.. math::
    :label: 2a

    \Delta E_{\text{strain}} = E_{\text{lig-pert}} - E_{\text{lig-eq}}
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{strain}} = \frac{1}{m} \sum_{i=0}^{m} \sum_{j=0}^{n}
    E_{\text{lig-pert}}(i, j) - E_{\text{lig-eq}}


.. math::
    :label: 3a

    \Delta E_{\text{strain}} = \Delta V_{\text{bond}} + \Delta V_{\text{angle}} +
    \Delta V_{\text{Urey-Bradley}} + \Delta V_{\text{dihedral}} + \Delta V_{\text{improper}} +
    \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`E_{\text{lig-eq}}` is herein the total energy of a (single) ligand at
its equilibrium geometry, while :math:`E_{\text{lig-pert}}(i, j)` is the
total energy of the (perturbed) ligand :math:`j` at MD iteration :math:`i`.

|

Interaction
-----------
The ensemble averaged interaction :math:`\Delta \overline{E}_{\text{int}}`
represents the mutual interaction between all ligands in a molecule.
The interaction is, again, averaged over all MD iterations and
summed over all ligand-pairs.

The magnitude of this term is determined by all non-covalent inter-ligand
interactions and can be either positive (dominated by Pauli and/or
Coulombic repulsion) or negative (dominated by dispersion and/or Coulombic
attraction).

.. math::
    :label: 4a

    \Delta E_{\text{int}} = \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (j, k)
    \quad \Rightarrow \quad
    \Delta \overline{E}_{\text{int}} = \frac{1}{m} \sum_{i=0}^{m}
    \sum_{j=0}^{n} \sum_{k \gt j}^{n} \Delta E_{\text{lig-int}} (i, j, k)

.. math::
    :label: 5a

    \Delta E_{\text{int}} = \Delta V_{\text{Lennard-Jones}} + \Delta V_{\text{elstat}}

:math:`\Delta E_{\text{lig-int}}(i, j, k)` represents the pair-wise
interactions between ligands :math:`j` and :math:`k` at MD iteration :math:`i`.
Double counting is avoided by ensuring that :math:`k > j`.

.. note::
    In order to avoid the substantial Coulombic repulsion between negatively charged ligands,
    its parameters are substituted with those from its neutral (*i.e.* protonated) counterpart.
    This correction is applied, exclusively, for the calculation of :math:`\Delta E_{\text{lig-int}}`.

|

Total Energy
------------
The total (ensemble-averaged) energy is the sum of
:math:`\Delta \overline{E}_{\text{strain}}` and
:math:`\Delta \overline{E}_{\text{int}}`.
Note that the energy is associated with a set of :math:`n` ligands,
*i.e.* the distortion and mutual interaction between all :math:`n` ligands.
Division by :math:`n` will thus yield the averaged energy per ligand
per MD iteration.

.. math::
    :label: 6a

    \Delta \overline{E} = \Delta \overline{E}_{\text{strain}} + \Delta \overline{E}_{\text{int}}
    = \frac{1}{m} \sum_{i=0}^{m} \Delta E_{\text{strain}}(i) + \Delta E_{\text{int}}(i)

|

Examples
--------
An example input script using the ``Cd68Se55`` core and ``OC(=O)CC`` ligand.

The :attr:`activation_strain.md<optional.qd.activation_strain.md>` key enables the MD-ASA procedure;
:attr:`activation_strain.use_ff<optional.qd.activation_strain.use_ff>` ensures
that the user-specified forcefield is used during the construction of the MD trajectory.

.. code:: yaml

    path: ...

    input_cores:
        - Cd68Se55.xyz:
            guess_bonds: False

    input_ligands:
        - OC(=O)CC

    optional:
        core:
            dummy: Cl

        ligand:
            optimize: True
            split: True

        qd:
            activation_strain:
                use_ff: True
                md: True
                job1: Cp2kJob

        forcefield:
            charge:
                keys: [input, force_eval, mm, forcefield, charge]
                Cd: 0.9768
                Se: -0.9768
                O2D2: -0.4704
                C2O3: 0.4524
            epsilon:
                unit: kjmol
                keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                Cd Cd: 0.3101
                Se Se: 0.4266
                Cd Se: 1.5225
                Cd O2D2: 1.8340
                Se O2D2: 1.6135
            sigma:
                unit: nm
                keys: [input, force_eval, mm, forcefield, nonbonded, lennard-jones]
                Cd Cd: 0.1234
                Se Se: 0.4852
                Cd Se: 0.2940
                Cd O2D2: 0.2471
                Se O2D2: 0.3526

|

activation_strain
-----------------
.. attribute:: optional.qd.activation_strain
    :noindex:

    All settings related to the activation strain analyses.

    Example:

    .. code:: yaml

        optional:
            qd:
                activation_strain:
                    use_ff: True
                    md: True
                    iter_start: 500
                    dump_csv: False

                    el_scale14: 1.0
                    lj_scale14: 1.0

                    distance_upper_bound: "inf"
                    k: 20
                    shift_cutoff: True

                    job1: cp2kjob
                    s1: ...


            forcefield:
                ...

|

    .. attribute:: optional.qd.activation_strain.use_ff

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Utilize the parameters supplied in the :attr:`optional.forcefield` block.


    .. attribute:: optional.qd.activation_strain.md

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Perform an ensemble-averaged activation strain analysis.

        If ``True``, perform the analysis along an entire molecular dynamics trajectory.
        If ``False``, only use a single geometry instead.


    .. attribute:: optional.qd.activation_strain.iter_start

        :Parameter:     * **Type** - :class:`int`
                        * **Default value** – ``500``

        The MD iteration at which the ASA will be started.

        All preceding iteration are disgarded, treated as pre-equilibration steps.
        Note that this refers to the iteration is specified in the .xyz file.
        For example, if a geometry is written to the .xyz file very 10 iterations
        (as is the default), then ``iter_start=500`` is equivalent to
        MD iteration 5000.


    .. attribute:: optional.qd.activation_strain.dump_csv

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Dump a set of .csv files containing all potential energies gathered over the course of the MD simulation.

        For each quantum dot two files are created in the ``.../qd/asa/`` directory,
        one containing the potentials over the course of the MD simulation (``.qd.csv``) and
        for the optimized ligand (``.lig.csv``).


    .. attribute:: optional.qd.activation_strain.el_scale14

        :Parameter:     * **Type** - :class:`float`
                        * **Default value** – ``1.0``

        Scaling factor to apply to all 1,4-nonbonded electrostatic interactions.

        Serves the same purpose as the cp2k EI_SCALE14_ keyword.


    .. attribute:: optional.qd.activation_strain.lj_scale14

        :Parameter:     * **Type** - :class:`float`
                        * **Default value** – ``1.0``

        Scaling factor to apply to all 1,4-nonbonded Lennard-Jones interactions.

        Serves the same purpose as the cp2k VDW_SCALE14_ keyword.


    .. attribute:: optional.qd.activation_strain.distance_upper_bound

        :Parameter:     * **Type** - :class:`float` or :class:`str`
                        * **Default value** – ``"inf"``

        Consider only atom-pairs within this distance for calculating inter-ligand interactions.

        Units are in Angstrom.
        Using ``"inf"`` will default to the full, untruncated, distance matrix.


    .. attribute:: optional.qd.activation_strain.k

        :Parameter:     * **Type** - :class:`int`
                        * **Default value** – ``20``

        The (maximum) number of to-be considered distances per atom.

        Only relevant when :attr:`distance_upper_bound != "inf"<optional.qd.activation_strain.distance_upper_bound>`.


    .. attribute:: optional.qd.activation_strain.shift_cutoff

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Add a constant to all electrostatic and Lennard-Jones potentials such that the potential is zero at the :attr:`distance upper bound<optional.qd.activation_strain.distance_upper_bound>`.

        Serves the same purpose as the cp2k SHIFT_CUTOFF_ keyword.
        Only relevant when :attr:`distance_upper_bound != "inf"<optional.qd.activation_strain.distance_upper_bound>`.


    .. attribute:: optional.qd.activation_strain.job1

        :Parameter:     * **Type** - :class:`type` or :class:`str`
                        * **Default value** – :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass,
        used for performing the activation strain analysis.

        Should be set to :class:`Cp2kJob<scm.plams.interfaces.thirdparty.cp2k.Cp2kJob>` if :attr:`activation_strain.md = True<optional.qd.activation_strain.md>`.


    .. attribute:: optional.qd.activation_strain.s1

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s1:
                input:
                    motion:
                        print:
                            trajectory:
                                each:
                                    md: 10
                        md:
                            ensemble: NVT
                            temperature: 300.0
                            timestep: 1.0
                            steps: 15000
                            thermostat:
                                type: CSVR
                                csvr:
                                    timecon: 1250

                    force_eval:
                        method: FIST
                        mm:
                            forcefield:
                                ei_scale14: 1.0
                                vdw_scale14: 1.0
                                ignore_missing_critical_params: ''
                                parmtype: CHM
                                parm_file_name: null
                                do_nonbonded: ''
                                shift_cutoff: .TRUE.
                                spline:
                                    emax_spline: 10e10
                                    r0_nb: 0.2
                            poisson:
                                periodic: NONE
                                ewald:
                                    ewald_type: NONE
                        subsys:
                            cell:
                                abc: '[angstrom] 100.0 100.0 100.0'
                                periodic: NONE
                            topology:
                                conn_file_format: PSF
                                conn_file_name: null
                                coord_file_format: 'OFF'
                                center_coordinates:
                                    center_point: 0.0 0.0 0.0

                    global:
                        print_level: low
                        project: cp2k
                        run_type: MD

        The job settings used for calculating the performing the ASA.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        The default settings above are specifically for the ensemble-averaged ASA
        (:attr:`activation_strain.md = True<optional.qd.activation_strain.md>`.).

.. _EI_SCALE14: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_EI_SCALE14
.. _VDW_SCALE14: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_VDW_SCALE14
.. _SHIFT_CUTOFF: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/MM/FORCEFIELD.html#list_SHIFT_CUTOFF

Welcome to the Compound Attachment/Analysis Tools' documentation!
=================================================================

Contents:

.. toctree::
   :maxdepth: 2

   includeme
   0_documentation
nanoCAT.recipes.cdft_utils
==========================

.. automodule:: nanoCAT.recipes.cdft_utils
.. _Bond Dissociation Energy:

Bond Dissociation Energy
========================

Calculate the bond dissociation energy (BDE) of ligands attached to the
surface of the core. The calculation consists of five distinct steps:

    1.  Dissociate all combinations of |n| ligands (|Y|, see :attr:`optional.qd.dissociate.lig_count`) a
    nd an atom from the core (|X|, see :attr:`optional.qd.dissociate.core_atom`)
    within a radius :math:`r` from aforementioned
    core atom (see :attr:`optional.qd.dissociate.lig_core_dist` and
    :attr:`optional.qd.dissociate.core_core_dist`).
    The dissociated compound has the general structure of |XYn|.

    2.  Optimize the geometry of |XYn| at the first level of theory
    (:math:`1`). Default: ADF MOPAC [1_, 2_, 3_].

    3.  Calculate the "electronic" contribution to the BDE (|dE|)
    at the first level of theory (:math:`1`): ADF MOPAC [1_, 2_, 3_].
    This step consists of single point calculations of the complete
    quantum dot, |XYn| and all |XYn|-dissociated quantum dots.

    4.  Calculate the thermalchemical contribution to the BDE (|ddG|) at the
    second level of theory (:math:`2`). Default: ADF UFF [4_, 5_]. This step
    consists of geometry optimizations and frequency analyses of the same
    compounds used for step 3.

    5.  :math:`\Delta G_{tot} = \Delta E_{1} + \Delta \Delta G_{2} = \Delta E_{1} + (\Delta G_{2} - \Delta E_{2})`.


Default Settings
~~~~~~~~~~~~~~~~

.. code:: yaml

    optional:
        qd:
            dissociate:
                core_atom: Cd
                core_index: null
                lig_count: 2
                core_core_dist: 5.0  # Ångström
                lig_core_dist: 5.0  # Ångström
                lig_core_pairs: 1
                topology: {}

                keep_files: True
                job1: AMSJob
                s1: True
                job2: AMSJob
                s2: True

|

Arguments
~~~~~~~~~

.. attribute:: optional.qd.dissociate
    :noindex:

    .. code:: yaml

        optional:
            qd:
                dissociate:
                    core_atom: Cd
                    core_index: null
                    lig_count: 2
                    lig_pairs: 1
                    core_core_dist: null  # Ångström
                    lig_core_dist: 5.0  # Ångström
                    topology:
                        7: vertice
                        8: edge
                        10: face

|

    .. attribute:: optional.qd.dissociate.core_atom

        :Parameter:     * **Type** - :class:`str` or :class:`int`

        The atomic number or atomic symbol of the core atoms (:math:`X`) which are to be
        dissociated. The core atoms are dissociated in combination with :math:`n` ligands
        (:math:`Y`, see :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`).
        Yields a compound with the general formula |XYn|.

        Atomic indices can also be manually specified with :attr:`dissociate.core_index<optional.qd.dissociate.core_index>`

        If one is interested in dissociating ligands in combination with
        a molecular species (*e.g.* :math:`X = {NR_4}^+`) the atomic number (or symbol)
        can be substituted for a SMILES string represting a poly-atomic ion
        (*e.g.* tetramethyl ammonium: C[N+](C)(C)C).

        If a SMILES string is provided it must satisfy the following 2 requirements:

            1. The SMILES string *must* contain a single charged atom; unpredictable behaviour can occur otherwise.
            2. The provided structure (including its bonds) must be present in the core.

        .. warning::
            This argument has no value be default and thus *must* be provided by the user.


    .. attribute:: optional.qd.dissociate.lig_count

        :Parameter:     * **Type** - :class:`int`

        The number of ligands, :math:`n`, which is to be dissociated in combination
        with a single core atom (:math:`X`, see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`).

        Yields a compound with the general formula |XYn|.

        .. warning::
            This argument has no value be default and thus *must* be provided by the user.


    .. attribute:: optional.qd.dissociate.core_index

        :Parameter:     * **Type** - :class:`int` or :class:`list` [:class:`int`], optional
                        * **Default value** – ``None``

        Alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_core_dist>` and
        :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`.
        Manually specify the indices of all to-be dissociated atoms in the core.
        Core atoms will be dissociated in combination with the :math:`n` closest ligands.

        .. note::
            Atom numbering follows the PLAMS [1_, 2_] convention of starting from 1 rather than 0.

        .. note::
            The yaml format uses ``null`` rather than ``None`` as in Python.


    .. attribute:: optional.qd.dissociate.core_core_dist

        :Parameter:     * **Type** - :class:`float` or :class:`int`, optional
                        * **Default value** – ``None``

        The maximum to be considered distance (Ångström) between atoms in
        :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`.
        Used for determining the topology of the core atom

        (see :attr:`dissociate.topology<optional.qd.dissociate.topology>`) and whether it is exposed to the
        surface of the core or not. It is recommended to use a radius which
        encapsulates a single (complete) shell of neighbours.

        If not specified (or equal to ``0.0``) **CAT** will attempt to guess a suitable value
        based on the cores' radial distribution function.


    .. attribute:: optional.qd.dissociate.lig_core_dist

        :Parameter:     * **Type** - :class:`float` or :class:`int`, optional
                        * **Default value** – ``None``

        Dissociate all combinations of a single core atom (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`)
        and the :math:`n` closests ligands within a user-specified radius.

        Serves as an alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_pairs>`,
        which removes a set number of combinations rather than everything withing a certain radius.

        The number of ligands dissociated in combination with a single core atom is controlled by
        :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`.

        .. image:: _images/BDE_XY2.png
            :scale: 25 %
            :align: center

|


    .. attribute:: optional.qd.dissociate.lig_pairs

        :Parameter:     * **Type** - :class:`int`, optional
                        * **Default value** – ``None``

        Dissociate a user-specified number of combinations of a single core atom (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>`)
        and the :math:`n` closests ligands.

        Serves as an alternative to :attr:`dissociate.lig_core_dist<optional.qd.dissociate.lig_core_dist>`,
        removing a preset number of (closest) pairs rather than all combinations within a certain radius.

        The number of ligands dissociated in combination with a single core atom is controlled by
        :attr:`dissociate.lig_count<optional.qd.dissociate.lig_count>`.


    .. attribute:: optional.qd.dissociate.topology

        :Parameter:     * **Type** - :class:`dict`
                        * **Default value** – ``{}``

        A dictionary which translates the number neighbouring core atoms
        (see :attr:`dissociate.core_atom<optional.qd.dissociate.core_atom>` and
        :attr:`dissociate.core_core_dist<optional.qd.dissociate.core_core_dist>`)
        into a topology. Keys represent the number of neighbours, values represent
        the matching topology.

        .. admonition:: Example

            Given a :attr:`dissociate.core_core_dist<optional.qd.dissociate.core_core_dist>` of ``5.0`` Ångström,
            the following options can be interpreted as following:

            .. code:: yaml

                optional:
                    qd:
                        dissociate:
                            7: vertice
                            8: edge
                            10: face

            Core atoms with ``7`` other neighbouring core atoms (within a radius of ``5.0`` Ångström)
            are marked as ``"vertice"``, the ones with ``8`` neighbours are marked as ``"edge"``
            and the ones with ``10`` neighbours as ``"face"``.


    .. attribute:: optional.qd.dissociate.qd_opt

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``False``

        Whether to optimize the quantum dot and |XYn| -dissociated quantum dot.


|

Arguments - Job Customization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. attribute:: optional.qd.dissociate
    :noindex:

    .. code:: yaml

        optional:
            qd:
                dissociate:
                    keep_files: True
                    job1: AMSJob
                    s1: True
                    job2: AMSJob
                    s2: True

|

    .. attribute:: optional.qd.dissociate.keep_files

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Whether to keep or delete all BDE files after all calculations are finished.

    .. attribute:: optional.qd.dissociate.xyn_pre_opt

        :Parameter:     * **Type** - :class:`bool`
                        * **Default value** – ``True``

        Pre-optimize the |XYn| fragment with UFF.

        .. note::
            Requires AMS.

    .. attribute:: optional.qd.dissociate.job1

        :Parameter:     * **Type** - :class:`type`, :class:`str` or :class:`bool`
                        * **Default value** – :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass, used for calculating the
        "electronic" component (|dE_lvl1|) of the bond dissociation energy.
        Involves single point calculations.

        Alternatively, an alias can be provided for a specific
        job type (see :ref:`Type Aliases`).

        Setting it to ``True`` will default to :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`,
        while ``False`` is equivalent to :attr:`optional.qd.dissociate` = ``False``.


    .. attribute:: optional.qd.dissociate.s1

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s1:
                input:
                    mopac:
                        model: PM7
                    ams:
                        system:
                            charge: 0

        The job settings used for calculating the "electronic" component
        (|dE_lvl1|) of the bond dissociation energy.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        Setting it to ``True`` will default to the ``["MOPAC"]`` block in
        CAT/data/templates/qd.yaml_, while ``False`` is equivalent to
        :attr:`optional.qd.dissociate` = ``False``.


    .. attribute:: optional.qd.dissociate.job2

        :Parameter:     * **Type** - :class:`type`, :class:`str` or :class:`bool`
                        * **Default value** – :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`

        A :class:`type` object of a :class:`Job<scm.plams.core.basejob.Job>` subclass, used for calculating the
        thermal component (|ddG_lvl2|) of the bond dissociation energy.
        Involves a geometry reoptimizations and frequency analyses.

        Alternatively, an alias can be provided for a specific
        job type (see :ref:`Type Aliases`).


        Setting it to ``True`` will default to :class:`AMSJob<scm.plams.interfaces.adfsuite.ams.AMSJob>`,
        while ``False`` will skip the thermochemical analysis completely.


    .. attribute:: optional.qd.dissociate.s2

        :Parameter:     * **Type** - :class:`dict`, :class:`str` or :class:`bool`
                        * **Default value** – See below

        .. code:: yaml

            s2:
                input:
                    uff:
                        library: uff
                    ams:
                        system:
                            charge: 0
                            bondorders:
                                _1: null

        The job settings used for calculating the thermal component (|ddG_lvl2|)
        of the bond dissociation energy.

        Alternatively, a path can be provided to .json or .yaml file
        containing the job settings.

        Setting it to ``True`` will default to the the *MOPAC* block in
        CAT/data/templates/qd.yaml_, while ``False`` will skip the
        thermochemical analysis completely.

Index
-----
.. currentmodule:: nanoCAT.bde.dissociate_xyn
.. autosummary::
    dissociate_ligand
    MolDissociater
    MolDissociater.remove_bulk
    MolDissociater.assign_topology
    MolDissociater.get_pairs_closest
    MolDissociater.get_pairs_distance
    MolDissociater.combinations
    MolDissociater.__call__

API
---
.. autofunction:: dissociate_ligand
.. autoclass:: MolDissociater
.. automethod:: MolDissociater.remove_bulk
.. automethod:: MolDissociater.assign_topology
.. automethod:: MolDissociater.get_pairs_closest
.. automethod:: MolDissociater.get_pairs_distance
.. automethod:: MolDissociater.combinations
.. automethod:: MolDissociater.__call__


.. _1: https://www.scm.com/doc/MOPAC/Introduction.html
.. _2: http://openmopac.net
.. _3: https://doi.org/10.1007/s00894-012-1667-x
.. _4: https://doi.org/10.1021/ja00051a040
.. _5: https://www.scm.com/doc/UFF/index.html
.. _qd.yaml: https://github.com/BvB93/CAT/blob/master/CAT/data/templates/qd.yaml

.. |dE| replace:: :math:`\Delta E`
.. |dE_lvl1| replace:: :math:`\Delta E_{1}`
.. |dE_lvl2| replace:: :math:`\Delta E_{2}`
.. |dG| replace:: :math:`\Delta G_{tot}`
.. |dG_lvl2| replace:: :math:`\Delta G_{2}`
.. |ddG| replace:: :math:`\Delta \Delta G`
.. |ddG_lvl2| replace:: :math:`\Delta \Delta G_{2}`
.. |XYn| replace:: :math:`XY_{n}`
.. |Yn| replace:: :math:`Y_{n}`
.. |n| replace:: :math:`{n}`
.. |X| replace:: :math:`X`
.. |Y| replace:: :math:`Y`
.. _Multi-ligand:

Multi-ligand attachment
=======================

.. attribute:: optional.qd.multi_ligand
    :noindex:

    All settings related to the multi-ligand attachment procedure.

    Example:

    .. code:: yaml

        optional:
            qd:
                multi_ligand:
                    ligands:
                        - OCCC
                        - OCCCCCCC
                        - OCCCCCCCCCCCC
                    anchor:
                        - F
                        - Br
                        - I

|

    .. attribute:: optional.qd.multi_ligand.ligands

        :Parameter:     * **Type** - :class:`list` [:class:`str`]

        SMILES strings of to-be attached ligands.

        Note that these ligands will be attached *in addition* to whichever ligands are
        specified in :ref:`Input Cores and Ligands`.

        .. note::
            This argument has no value be default and must thus be provided by the user.


    .. attribute:: optional.qd.multi_ligand.anchor

        :Parameter:     * **Type** - :class:`list` [:class:`str` or :class:`int`]

        Atomic number of symbol of the core anchor atoms.

        The first anchor atom will be assigned to the first ligand in
        :attr:`multi_ligand.ligands<optional.qd.multi_ligand.ligands>`, the second anchor atom
        to the second ligand, *etc.*.
        The list's length should consequently be of the same length as
        :attr:`multi_ligand.ligands<optional.qd.multi_ligand.ligands>`.

        Works analogous to :attr:`optional.core.anchor`.

        This optiona can alternatively be provided as ``optional.qd.multi_ligand.dummy``.

        .. note::
            This argument has no value be default and must thus be provided by the user.
.. automodule:: nanoCAT.recipes.charges
.. _distribution:

Subset Generation
=================
.. automodule:: CAT.distribution
guess_core_dist
===============

.. automodule:: nanoCAT.bde.guess_core_dist
