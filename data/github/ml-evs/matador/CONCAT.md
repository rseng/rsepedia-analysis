# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

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

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at matador@ml-evs.science. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
---
title: 'matador: a Python library for analysing, curating and performing high-throughput density-functional theory calculations'
tags:
  - Python
  - density-functional theory
  - ab initio
  - crystal structure prediction
  - materials discovery
  - databases
  - castep
  - quantum espresso
  - mongodb
authors:
  - name: Matthew L. Evans
    orcid:  0000-0002-1182-9098 
    affiliation: 1
  - name: Andrew J. Morris
    orcid: 0000-0001-7453-5698
    affiliation: 2
affiliations:
 - name: Theory of Condensed Matter Group, Cavendish Laboratory, University of Cambridge, J. J. Thomson Avenue, Cambridge, CB3 0HE, U.K.
   index: 1
 - name: School of Metallurgy and Materials, University of Birmingham, Edgbaston, Birmingham, B15 2TT, U.K.
   index: 2

date: July 2020
bibliography: paper.bib
---

# Summary

The properties of materials depend heavily on their atomistic structure; knowledge of the possible stable atomic configurations that define a material is required to understand the performance of many technologically and ecologically relevant devices, such as those used for energy storage [@NaP; @CuP]. First-principles crystal structure prediction (CSP) is the art of finding these stable configurations using only quantum mechanics [@JM]. Density-functional theory (DFT) is a ubiquitous theoretical framework for finding approximate solutions to quantum mechanics; calculations using a modern DFT package are sufficiently robust and accurate that insight into real materials can be readily obtained. The computationally intensive work is performed by well-established, low-level software packages, such as CASTEP [@castep] or Quantum Espresso [@qe], which are able to make use of modern high-performance computers. In order to use these codes easily, reliably and reproducibly, many high-level libraries have been developed to create, curate and manipulate the calculations from these low-level workhorses; `matador` is one such framework.

# Statement of need

The purpose of `matador` is fourfold:

- to promote the use of local databases and high-throughput workflows to increase the reproducibility of the computational results, 
- to perform reliable analysis of the stability, structure and properties of materials derived from calculations, 
- to provide tools to create customisable, publication-quality plots of phase diagrams, spectral properties and electrochemistry,
- to make the above functionality available to those with limited programming experience.

# `matador`

`matador` is a Python 3.6+ library and set of command-line tools for performing and analysing high-throughput DFT calculations using the CASTEP [@castep] and Quantum Espresso [@qe] packages. It is well-tested and fully-documented at [ReadTheDocs](https://matador-db.readthedocs.io), and comes with several tutorials and examples. The package is available on PyPI under the name [`matador-db`](https://pypi.org/project/matador-db). As with many projects, `matador` is built on top of the scientific Python ecosystem of NumPy [@numpy], SciPy [@scipy] and matplotlib [@matplotlib].

`matador` has been developed with high-throughput CSP in mind and has found use in the application of CSP to energy storage materials [@NaP; @CuP]; in this use case, a single compositional phase diagram can consist of tens of thousands of structural relaxation calculations. This package is aimed at users of CASTEP or Quantum Espresso who are comfortable with the command-line, yet maybe lack the Python knowledge required to start from scratch with more sophisticated packages. There are many mature packages that provide overlapping functionality with `matador`, the most widespread of which being the Atomic Simulation Environment (ASE) [@ase] and pymatgen [@pymatgen]. A translation layer to and from the structure representation of both of these packages is provided, such that analysis can be reused and combined. 

![Li-Zn-P ternary phase diagram created with matador, plot generated with matplotlib [@matplotlib] and python-ternary [@ternary].\label{fig:hull}](./LiZnP_hull.pdf)

# Overview of functionality

There are two ways of working with `matador`, either from the command-line interface (CLI) or through the Python library directly, with some features that are unique to each. The functionality of `matador` can be broadly split into three categories: 

1. *Creation and curation of databases of the results of first-principles calculations.*

`matador` allows for the creation of [MongoDB](https://mongodb.com) databases of CASTEP (6.0+) geometry optimisations from the command-line, using `matador import`. Only calculations that are deemed "successful", and that have fully-specified inputs are stored, with errors displayed for the rest. The resulting database can be queried with `matador query`, either with Python or through the powerful CLI. The results can be filtered for structural "uniqueness" and written to one of several supported file types and exported for use in other frameworks, such as ASE or pymatgen. Prototyping of structures directly from the database is achieved using `matador swaps`, which uses the same interface as `matador query` to return structure files with "swapped" elements [@NaP].

2. *High-throughput calculations and automated workflows.*

The `run3` executable bundled with `matador` allows for high-throughput calculations to be performed with little setup and no programming knowledge. Specialised support for CASTEP and the post-processing tool OptaDOS [@optados; @optados2] is provided to perform high-throughput geometry optimisations, orbital-projected band structures and densities of states, phonon calculations and elastic properties, however `run3` can also be used to run generic MPI programs concurrently on a set of structures. Sensible defaults for these workflows are provided by leveraging the open-source SeeK-path [@seekpath] and spglib [@spglib] libraries. The bundled `dispersion` script and associated library functionality allows for the creation of publication-quality spectral and vibrational property plots, in a similar fashion to the `sumo` package [@sumo]. The `matador.compute` module behind `run3` also powers the `ilustrado` genetic algorithm code [@ilustrado].

3. *Stability and structural analysis (with an emphasis on battery materials).*

The construction of reliable compositional phase diagrams requires several independent calculations to be performed on different atomic configurations with a compatible set of external parameters. These can be generated from a database query using `matador hull`, which allows the user to filter between different sets of calculations, and, where relevant, `matador voltage` can provide the electrochemical properties of that same phase diagram. Structural fingerprints implemented include pair distribution functions, powder X-ray diffraction patterns, and periodic crystal bond graphs. As more calculations are performed, changes to phase diagrams stored in the local database can be tracked with `matador hulldiff`. Phase diagrams can also be constructed from multiple energy values per structure, for example to show the effects of finite temperature [@CuP], or in the specific case of ensemble-based exchange-correlation functionals like the Bayesian Error Estimate Functional (BEEF) [@beef]. An example of a ternary phase diagram is shown \autoref{fig:hull}.

# Acknowledgements

We acknowledge all the contributors, users and testers of this package, primarily Angela Harper, James Darby, Jordan Dorrell and Matthew Cliffe. M.E. would like to acknowledge the EPSRC Centre for Doctoral Training in Computational Methods for Materials Science for funding under grant number EP/L015552/1. A.J.M. acknowledges funding from EPSRC (EP/P003532/1). The authors acknowledge networking support via the EPSRC Collaborative Computational Projects, CCP9 (EP/M022595/1) and CCP-NC (EP/T026642/1). Much of the development and testing was performed on the Cambridge Service for Data Driven Discovery (CSD3) operated by the University of Cambridge Research Computing Service (http://www.csd3.cam.ac.uk/), provided by Dell EMC and Intel using Tier-2 funding from the Engineering and Physical Sciences Research Council, and DiRAC funding from the Science and Technology Facilities Council (www.dirac.ac.uk).
Date: 13:41 16/04/2020  
Command: matador voltage -int -c NaSnP --db me388_NaPSn --pathways -hc 0.0 --res  
Version: 0.8b1-772-g94db115-voltages-and-volumes  

--------------------------------------------------------------------------------------------------------------------------
```
                  Root                   Pressure  Volume/fu  Hull dist.   Space group     Formula      # fu   Prov.  
                                          (GPa)     (Ang^3)   (meV/atom)  
--------------------------------------------------------------------------------------------------------------------------
Na-OQMD_8535-CollCode44758                   0.04       36.1          0.0  P6_3/mmc         Na          2      ICSD  
Na15Sn4-OQMD_90363-CollCode105167           -0.00      560.7          0.0    I-43d        Na15Sn4       2      ICSD  
NaP-CollCode182164                           0.02       93.9          0.0   P6_3cm         Na3P         6      ICSD  
NaSn-CollCode167672                         -0.01      283.7          0.0   P3_212        Na7Sn3        6      ICSD  
NaPSn-GA-dw2swh-1x103                       -0.03      320.0          0.0     P-1         Na8P4Sn       2       GA   
NaP-CollCode56530                           -0.02      196.1          0.0    C2/m          Na5P4        1      ICSD  
NaPSn-GA-dw2swh-1x79                        -0.06      208.5          0.0     P-1         Na5P3Sn       4       GA   
NaP-CollCode421420                          -0.04       41.9          0.0   P2_1/c          NaP         8      ICSD  
NaSn-CollCode641294                         -0.01       59.5          0.0  I4_1/acd        NaSn         16     ICSD  
NaPSn-GA-npf6il-1x3                         -0.04       78.2          0.0     P1           NaPSn        2       GA   
NaPSn-GA-npf6il-2x47                        -0.04       78.2          0.0     P1           NaPSn        2       GA   
NaP-CollCode60774                           -0.03      251.5          0.0 P2_12_12_1       Na3P7        4      ICSD  
NaSn-NaSn2-Baggeto                           0.02       80.1          0.0    Cmmm          NaSn2        2     SWAPS  
NaP-na3p11collo                              0.01      326.1          0.0    Pbcn         Na3P11        4      ICSD  
NaSn-146907-4243-58                         -0.02      105.6          0.0    R-3m          NaSn3        4     SWAPS  
NaP-NaP7                                     0.03      192.1          0.0  I4_1/acd        NaP7         8     SWAPS  
P-OQMD_5708-CollCode29273                    0.02       23.6          0.0    P2/c            P          84     ICSD  
P3Sn-OQMD_3387-CollCode16293                 0.00       85.9          0.0    R-3m          P3Sn         2      ICSD  
P3Sn4-OQMD_2953-CollCode15014                0.01      165.8          0.0    R-3m          P3Sn4        1      ICSD  
Sn-Collo                                     0.03       28.2          0.0  I4_1/amd         Sn          2      ICSD  
```Date: 11:58 23/04/2020  
Command: matador voltage -c LiCoP -int --db oqmd --spin any --volume --res  
Version: 0.8b1-772-g94db115-voltages-and-volumes  

--------------------------------------------------------------------------------------------------------------------------
```
                  Root                   Pressure  Volume/fu  Hull dist.   Space group     Formula      # fu   Prov.  
                                          (GPa)     (Ang^3)   (meV/atom)  
--------------------------------------------------------------------------------------------------------------------------
OQMD 5154                                    0.00       55.8          0.0   P63/mmc        Li3P         2      ICSD  
OQMD 12657                                   0.00      122.3          0.0    Amm2         Co6LiP4       1      ICSD  
OQMD 8021                                    0.00       22.8          0.0    Pnma           CoP         4      ICSD  
OQMD 6978                                    0.00       38.9          0.0    P21/c         CoP2         4      ICSD  
OQMD 4397                                    0.00       56.9          0.0    Im-3          CoP3         4      ICSD  
OQMD 17007                                   0.00       30.3          0.0    P21/c          LiP         8      ICSD  
OQMD 4845                                    0.00       31.7          0.0    Pnma          Co2P         4      ICSD  
OQMD 11352                                   0.00      192.0          0.0   P212121        Li3P7        4      ICSD  
OQMD 4371                                    0.00      164.2          0.0   I41/acd        LiP7         8      ICSD  
OQMD 24584                                   0.00       23.8          0.0     P-1            P          42     ICSD  
OQMD 594128                                  0.00       10.8          0.0    Fm-3m          Co          1      ICSD  
OQMD 30676                                   0.00       18.2          0.0    R-3m           Li          3      ICSD  
```---
name: Issue template
about: Default issue template
title: ''
labels: ''
assignees: ''

---

If you think you have found a bug, please raise an issue on GitHub, providing information about what you were trying to do, the function/script you ran, the error message/output and your matador version. If you are able to, please try to replicate the problem on the master branch before posting.

If you have a feature request, or you would like to contribute a feature, please raise an issue so that its suitability for this package can be discussed. Ideally, you could do this before writing your implementation, as the feature may already exist! Code contributions must include tests for any new functionality and must pass existing tests and flake8 linting run by the CI. Any new modules must use the black auto-formatter.
# Examples

This folder contains an assortment of examples, some as Jupyter notebooks and
some as tutorials to be followed according to the text files within. 

A full explanation, with links to the appropriate Binder can be found in the
[Examples section](https://docs.matador.science/en/stable/examples_index.html)
of the online documentation.

Examples in the `./interactive/` folder can be run without any external data,
either locally or with Binder. Those in the `./non-interactive/` folder require
external data but provide useful code snippets. 
.. _contributing:

Contributing
============

Contributions and suggestions to this package are very welcome.

If you think you have found a bug, please raise an issue on GitHub, providing information about what you were trying to do, the function/script you ran, the error message/output and your ``matador`` version. If you are able to, please try to replicate the problem on the ``master`` branch before posting.

If you have a feature request, or you would like to contribute a feature, please raise an issue so that its suitability for this package can be discussed. Ideally, you could do this before writing your implementation, as the feature may already exist! Code contributions must include tests for any new functionality and must pass existing tests and ``flake8`` linting run by the CI. Any new modules must use the ``black`` auto-formatter.
.. _changelog:

Changelog
=========

New in release (0.9.11) [03/06/2021]
------------------------------------

- Minor change: allow specification of arbitrary strings for CASTEP pseudopotential libraries (#156).
- Bug fix: ``standardize_cell`` script failing to use default symmetry tolerance value (#157).
- Bug fix: scraping of .cif files with single atoms and no symmetries (#173)
- Bug fix: scraping of Hubbard U values from .castep files, and associated bugs when performing relaxations with Hubbard U (#180)
- Dependency updates and Python 3.6 deprecation warning (#158, #181)

New in release (0.9.10) [23/02/2021]
------------------------------------

- Windows compatibility changes (#149)
- Dependency updates (#146, #148, #149)

New in release (0.9.9) [16/10/2020]
-----------------------------------

- Added support for CASTEP kpoint path ``BREAK`` directive (#107)
- Improvements to magres plotting and magres workflow (#112)
- Added ability to scrape electric field gradient values and compute quadrupolar quantities from NMR calculations (#115)
- Added ability to run all several examples under Binder (#106, #130).
- JOSS paper accepted! (#129)


New in release (0.9.8) [10/08/2020]
-----------------------------------
- Improvements to PDIS functionality (#94).

  - Rasterized scatter points for more efficient exporting and fewer graphical artifacts
  - Made underlying :func:`matador.plotting.spectral_plotting.dos_plot` and :func:`matador.plotting.spectral_plotting.dispersion_plot` more API friendly, and added example notebook.
  - Fixed bug in cell scraping for old ``BS_*`` style keywords.

- Improvements to magres functionality, including scraping of units (#90)
- Example notebooks that do not need external data/databases are now run as part of CI (#91).
- New workflow for NMR calculations and refactoring of old workflows (#96).

  - New workflow performs relaxation and high-quality SCF before NMR calculation.
  - Old workflows refactored and improved to enforce certain required parameters for e.g. checkpointing.
  - Enabled phonon workflow for CASTEP ``PHONON+EFIELD`` task.
  - Made file scrapers less dependent on file type.

- Updated CASTEP parameter list to 20.1 (#97).
- Tweaked spectral plotting defaults, including ``--colours`` flag to dispersion script (#98).


New in release (0.9.7) [29/07/2020]
-----------------------------------
- Bug fixes to problems introduced in 0.9.6.
- Cosmetic fixes to logging and misleading status reports in workflows.


New in release (0.9.6) [28/07/2020]
-----------------------------------
- Improvements to ASE and pymatgen interoperability (#80)
- Fixed bug in :class:`matador.hull.TemperatureDependentHull` which would crash when not provided a list of temperatures (#82).
- Added plotting functions for magres data, and improved its handling inside :class:`matador.crystal.Crystal` (#79).

New in release (0.9.5) [25/06/2020]
-----------------------------------
- This release is mostly to trigger Zenodo archiving.
- Updated README and tests for recent Python versions.


New in release (0.9.4) [08/06/2020]
-----------------------------------
- Fixed flag help strings for ``pxrd_calculator`` (#65)
- Changed default PDF broadening for 3x speedup (#65)
- Reverted ``cpu_count`` to use version that works correctly in most cases, by chance (#66).


New in release (0.9.3) [07/06/2020]
-----------------------------------

- Fixes for the CIF reader: now works with awkward linebreaks and alternative symmetry operation specifications (#61).
- Added several new flags to ``pxrd_calculator`` script (#60 and 61).
- Usability fixes for ``spectral_plotting`` in the case of projected dispersion curves (#59).


New in release (0.9.2) [01/06/2020]
-----------------------------------

- Optimised CIF reader considerably (#50)
- Updated PXRD calculator to allow for partial occupancy, monochromated beam angles and text export, and added ``pxrd_calculator`` script for convenience when handling CIF files.
- Added ability to choose which projectors are plotted with dispersion (#47)
- Various minor fixes and updates:

  - Updates to docs for CLI and configuration.
  - Allow nan-values to be reset inside :class:`matador.crystal.Crystal`.
  - Fixed display ordering of fingerprint-filtered cursors.


New in release (0.9.1) [20/05/2020]
-----------------------------------

- Fixed issue with local pip installs after 0.9 release
- Fixed issue with multi-node MPI tasks by switching to ``proc.communicate()`` after an initial polling stage (#37)
- Fixed issue where bands would be reordered multiple times in spectral plots (#40)
- Tweaked spectral plot defaults (#40)
- Replaced ``multiprocessing.cpu_count()`` calls with ``psutil.cpu_count(logical=False)`` to avoid double-counting hyperthreaded cores


New in release (0.9) [15/05/2020]
---------------------------------

- PyPI release! Can now install with ``pip install matador-db`` (unfortunately ``matador`` was taken, but they are sufficiently orthogonal that the package name ``matador`` is retained here.
- Much improved code structure and many additional classes that wrap the raw calculation dictionaries for e.g. :class:`matador.crystal.Crystal` and spectral classes.
- New module :mod:`matador.orm` containing useful models for data handling.

  - :class:`matador.orm.orm.DataContainer` as a base class that allows for easy
    access to underlying dictionaries.
  - :mod:`matador.orm.spectral` module that contains many useful classes for
    manipulating and plotting e.g. bandstructures, DOS and finite temperature
    properties.

- New features in :mod:`matador.hull` module:

  - Pseudo-ternary phase diagrams (building towards arbitrary n-dimensional phase diagrams).
  - :class:`matador.hull.EnsembleHull` class and submodule to support the Bayesian Error Estimate Functional (BEEF) and finite temperature phase diagrams.
  - Refactoring of hull calculation into light-weight :class:`matador.hull.PhaseDiagram` class.
  - Finite temperature hulls based on :class:`matador.hull.EnsembleHull` with
    :class:`matador.hull.TemperatureDependentHull`.

- Refactored old PDF `similarity` module into new module :mod:`matador.fingerprints`.

  - Added new fingerprint class, :class:`matador.fingerprints.PXRD`, with associated plots (thanks for James Darby for some initial code). Defaults calibrated with GSAS-II.
  - :class:`matador.fingerprints.PDF` sped up by an order of magnitude using `numba`.

- :class:`matador.workflows.castep.CastepSpectralWorkflow` extended to include latest projected dispersion curve developments from OptaDOS, with associated projected dispersion plots (see tutorial).

  - Updated dispersion script to automatically perform naive Gaussian smearing if OptaDOS output not detected.

- Abstracted and simplified :mod:`matador.compute` module to allow for extension to new codes via :mod:`matador.compute.calculators` submodule.

  - Should now be more robust and transferrable, with many HPC environments automatically detected.
  - Added ``--scratch_prefix`` to run3 to allow for temporary files to e.g. be written to faster filesystem with appropriate symlinks to work folder.

- All CASTEP 19 keywords supported, as well as `devel_code` blocks.
- Several new tests: coverage now around 75% when CASTEP is available.
- New tutorials:

  - :ref:`MongoDB setup<mongo>`
  - :ref:`Spectral calculations with run3<run3_spectral>`
  - Example notebooks


New in release (0.8b) [03/08/2018]
----------------------------------

- Wholesale changes, complete refactoring of most of the code.
- Released open source under the MIT license!
- Documentation now hosted on `readthedocs <matador-db.readthedocs.org>`_,
- Workflows: chaining up job steps with run3:

  - spectral and phonons (combined DOS, dispersion calculations) with automated kpoint paths.
  - bulk modulus calculations and EOS fitting.

- New tutorials:

  - :ref:`Geometry optimisations with run3<run3_geom>`

- Temperature-dependent convex hulls (thanks to Angela Harper).
- New per-used configuration that allows changing of plotting styles, colourschemes, database names etc.
- Improvements to compute module:

  - automatically handle walltime constraints for Slurm and PBS.
  - estimate memory usage with CASTEP and skip if exceeds machine capacity,

- All CASTEP 18 keywords supported.
- Better support for electronic structure data, OptaDOS, NMR calculations, CIF files, partial occupancy.


New in version (0.7b) [13/04/2017]
----------------------------------

-  Ternary voltage curves.
-  Similarity/uniqueness filtering with element-projected PDFs.
-  Updated compute engine for remote calculations (see ``compute.py`` and new script ``oddjob``).
-  Improved test suite and full pip compatiblity.
-  Many bugfixes and usability changes.

New in version (0.6b) [01/06/2017]
----------------------------------

-  Intercalation voltage curves, e.g. ``matador voltage -c Li:SnS2``.
-  Ternary phase diagrams with heatmaps for structure prediction sampling, gravimetric capacity and formation enthalpy ``matador hull -c ABC --sampmap --efmap --capmap``.
-  Substructural similarity interface with Can Kocer's code, as proposed by `Yang et al., PRB (2014) <http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.054102>`_.
.. _install:

Installation
============

If you have any issues with installation, feel free to raise an issue on GitHub outlining your approach and any errors you received.


Simple installation with pip
----------------------------

The matador package can be found on PyPI under the name `matador-db <https://pypi.org/project/matador-db>`_ and installed with
``pip install matador-db``, preferably in a fresh virtual environment (see conda instructions below). Extra dependencies may be installed with e.g. ``pip install matador-db[all]``.

Development installation with conda/pip
---------------------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

1. Clone the matador source onto your local machine ``git clone https://github.com/ml-evs/matador.git``.

Optional (but recommended) steps:

2. `Install conda <https://conda.io/miniconda.html>`_, if you have not already. There may be a package available already if you are using a supercomputer (e.g. `anaconda-compute/2.2.0-python3` on ARCHER 30/10/2017).
3. Create a new conda environment to install matador into (``conda create -n matador python=3.7``) and activate it with (``conda activate matador``).
4. Install some of the heavier requirements (e.g. NumPy and SciPy) through conda with ``conda install --yes --file requirements/requirements.txt``.

Required steps:

5. Run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install.
6. You now have a basic matador API installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/`` using either conda or pip. If you wish to just install everything use ``pip install .[all]``.
7. To use matador, you will need to activate the conda environment from step 2, by running ``conda activate matador``. You will also need this in e.g. any job scripts. After installing the test dependencies with ``pip install .[test]``, you can test your installation using ``python -m unittest discover -v -b`` or simply ``py.test``. By default this will look for an MPI-enabled executable called ``castep`` on your ``$PATH`` to run CASTEP tests.

Troubleshooting
---------------

Below are some problems encountered on various machines that may be helpful:

1. (10/09/2019) When installing with ``conda``, if you receive the following error (or
   similar): ``/home/#####/.local/conda/envs/matador/compiler_compat/ld: build/temp.linux-x86_64-3.6/psutil/_psutil_common.o: unable to initialize decompress status for section .debug_info``, then you are using a modern compiler that breaks ``conda``'s attempts to be backwards compatible (in this case it was GCC 9). The simple fix is to rename/remove the copy of ``ld`` inside your conda environment (path in the message above) such that your system ``ld`` is used.
2. (10/10/2017) On some machines (e.g. ARCHER/Thomas) you may receive permissions errors at step 5; if so, try moving matador's `.git` and install again (``mv .git $HOME/matador_git_stash; pip install . ; mv $HOME/matador_git_stash .git``).
3. Some dependencies may not have compiled packages (wheels) for your distribution on PyPI, and may have compilation errors. In this case, you could try finding the relevant package on conda instead.
=======
matador
=======

| |PyPI Version| |GH Actions| |Binder|
| |Documentation Status| |MIT License| |Coverage Status|
| |JOSS| |Zenodo|


matador is an aggregator, manipulator and runner of first-principles calculations, written with a bent towards battery electrode materials.
The source can be found on `GitHub <https://github.com/ml-evs/matador>`_ and online documentation is hosted at `ReadTheDocs <https://docs.matador.science>`_.

Example Jupyter notebooks and tutorials can be found `online <https://docs.matador.science/en/latest/examples_index.html>`_ or in the ``examples/`` folder of the matador source code.

Written & maintained by `Matthew Evans <https://ml-evs.science>`_ (2016-).


.. image:: docs/src/img/lipzn.png
   :name: LiPZn
   :align: center

Installation
------------

In the simplest case (e.g. you already have Python 3.6+ set up), ``pip install matador-db`` is sufficient to get up and running, preferably in a fresh virtual environment.

Upgrading to the latest version should be as simple as ``pip install -U matador-db``.

For an editable development installation, clone the source code from this repository and run ``pip install -e .`` from the matador folder. Tests can be run on your local machine with ``python -m unittest discover -v -b`` or simply with ``py.test`` after test dependencies have been installed with ``pip install .[test]``. In order to test CASTEP-running functionality, the tests will look for an MPI-enabled executable named ``castep`` on your ``$PATH``.

Further instructions can be found in the `Installation instructions <https://docs.matador.science/en/latest/install.html>`_.


Usage
------

``matador`` is primarily a Python *library* that can be used inside Python scripts/modules to create a custom workflow. There are, however, several command-line scripts bundled with ``matador`` itself. All of these scripts are listed under `CLI Usage <https://docs.matador.science/en/latest/cli.html>`_.

For basic command-line usage, please explore the help system for each command. Common workflows can be found inside ``examples/`` and in the `online docs <http://docs.matador.science/en/latest/examples_index.html>`_.

Please consult the full `Python API documentation <http://docs.matador.science/en/latest/modules.html>`_ for programmatic usage.

Core functionality
-------------------

The API has many features that can be explored in the examples and API documentation. As a summary, ``matador`` can be used for:

- Scraping of CASTEP (and Quantum Espresso) input/output files into flexible Python dictionaries/models.
- The creation and curation of MongoDB collections of geometry optimisation calculations, with a powerful querying CLI/API.
- Customisable, publication-ready plots for all models, e.g. phase diagrams, PDF, PXRD, voltage profiles, electronic/vibrational bandstructures etc.
- High-throughput geometry optimisations, electronic and vibrational properties using CASTEP (and Quantum Espresso) with ``run3``. Tested on several supercomputers. ``run3`` is designed primarily for simple workflows and offers little in the way of tools for creating complex workflows; if this is your use case, then check out some of the other codes listed below.
- Creation of phase diagrams and electrochemical voltage profiles from the results of DFT calculations.

This functionality is achieved by interfacing with much of the standard scientific Python stack (`NumPy <https://numpy.org>`_, `SciPy <https://scipy.org>`_, `matplotlib <https://matplotlib.org>`_), some more materials-specific packages (`spglib <https://github.com/atztogo/spglib/>`_, `SeeK-path <https://github.com/giovannipizzi/seekpath>`_, `periodictable <https://github.com/pkienzle/periodictable>`_) and other general packages (`pymongo <https://github.com/mongodb/mongo-python-driver>`_, `python-ternary <https://github.com/marcharper/python-ternary>`_, `numba <https://numba.org>`_).

Similar packages
----------------

This package is by no means unique in its functionality or goals. Below is a list of similar packages and an overview of where they overlap with ``matador``:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_: manipulation of structures, parsing and exporting files, running jobs and local databases with ``ase-db``. An interface is provided to ASE's ``Atoms`` object.
- `pymatgen <https://pymatgen.org>`_: similar to ASE, with a focus on the VASP DFT code. An interface is provided to pymatgen's ``Structure`` object. Local databases can be constructed with the `pymatgen-db <https://github.com/materialsproject/pymatgen-db>`_ add-on and high-throughput workflows are achieved with `Fireworks <https://github.com/materialsproject/fireworks>`_.
- `AiiDA <https://www.aiida.net>`_: high-throughput job running, provenance tracking and database storage for many simulation codes.
- `sumo <https://github.com/SMTG-UCL/sumo>`_ publication quality plotting, primarily for VASP but also with support for other codes.

If you think this list is outdated, incorrect or simply incomplete, then please raise an issue!

Citing matador
--------------

If you use matador in your work, we kindly ask that you cite

    Matthew L. Evans, Andrew J. Morris, *matador: a Python library for analysing, curating and performing high-throughput density-functional theory calculations* Journal of Open Source Software, 5(54), 2563 (2020)
    `10.21105/joss.02563 <https://doi.org/10.21105/joss.02563>`_

Source code archives for all versions above 0.9 can be found at Zenodo `DOI 10.5281/zenodo.3908573 <https://doi.org/10.5281/zenodo.3908573>`_.


.. |PyPI Version| image:: https://img.shields.io/pypi/v/matador-db?label=PyPI&logo=pypi
   :target: https://pypi.org/project/matador-db/
.. |GH Actions| image:: https://img.shields.io/github/workflow/status/ml-evs/matador/Run%20tests/master?label=master&logo=github
   :target: https://github.com/ml-evs/matador/actions?query=branch%3Amaster
.. |MIT License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/ml-evs/matador/blob/master/LICENSE
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/ml-evs/matador/master?logo=codecov
   :target: https://codecov.io/gh/ml-evs/matador
.. |Documentation Status| image:: https://readthedocs.org/projects/matador-db/badge/?version=stable
   :target: https://matador-db.readthedocs.io/en/stable/?badge=stable
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3908573.svg
   :target: https://doi.org/10.5281/zenodo.3908573
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/ml-evs/matador/master?filepath=examples/interactive
.. |JOSS| image:: https://joss.theoj.org/papers/4d0eea9bea4362dec4cb6d62ebccc913/status.svg
   :target: https://joss.theoj.org/papers/4d0eea9bea4362dec4cb6d62ebccc913
.. index:: run3_elastic

.. highlight:: bash

.. _run3_elastic:


Example 4: Bulk moduli with CASTEP and run3
-------------------------------------------

In this tutorial we will calculate the bulk moduli of diamond, silicon and lithium, using CASTEP and run3. 

The general process for calculating the bulk modulus is as follows:

1. Relax a structure to its equilibrium volume :math:`V_0` (ideally relaxing the positions
   too).
2. Run total energy calculations at a series of perturbed volumes :math:`\alpha V_0`
   for :math:`\alpha \in [0.95, 1.05]`.
3. Fit a particular form of equation of state to the resulting :math:`E(V)` curve and extract the bulk modulus.

Note: there is a supported set of scripts for calculating all the elastic constants in CASTEP, which
can be found on `GitHub <https://github.com/andreww/elastic-constants>`_ with an additional tutorial on the `CASTEP website <http://www.castep.org/Tutorials/ElasticConstants>`_.

As with the other tutorials, run3 expects to find a series of structures as ``.res`` files and single ``$seed.cell`` and ``$seed.param`` files. The files for this example can be found in ``examples/run3_elastic``. The ``.cell`` file in this tutorial is basically identical to that in the geometry optimisation tutorial, but the ``.param`` file this time contains the line ``task: bulk_modulus``. This is *not* a valid CASTEP task (as of 2019), but instead run3 will capture this task and use it to spawn a series of single point jobs.::

    $ cat bulk_mod.cell
    kpoints_mp_spacing: 0.05
    snap_to_symmetry
    symmetry_tol: 0.01
    symmetry_generate
    %block species_pot
    QC5
    %endblock species_pot

    $ cat bulk_mod.param
    task: bulk_modulus
    cut_off_energy: 300 eV
    write_bib: false
    write_checkpoint: none
    geom_max_iter: 100
    xc_functional: LDA


Once the files are in place, calling run3 as ``run3 bulk_mod`` will begin to relax
the first structure before deforming it. You can track the calculations
progress in the log files found in ``logs/``.

After the calculations have completed, your completed folder should look
something like this::

   $ ls completed
    C-OQMD_675640-CollCode28857.bands             Li-OQMD_30676-CollCode642106_bulk_mod.castep
    C-OQMD_675640-CollCode28857.castep            Li-OQMD_30676-CollCode642106_bulk_mod.cell
    C-OQMD_675640-CollCode28857.cell              Li-OQMD_30676-CollCode642106_bulk_mod.cst_esp
    C-OQMD_675640-CollCode28857.cst_esp           Li-OQMD_30676-CollCode642106_bulk_mod.param
    C-OQMD_675640-CollCode28857.geom              Li-OQMD_30676-CollCode642106_bulk_mod.png
    C-OQMD_675640-CollCode28857.param             Li-OQMD_30676-CollCode642106_bulk_mod.res
    C-OQMD_675640-CollCode28857.res               Li-OQMD_30676-CollCode642106_bulk_mod.results
    C-OQMD_675640-CollCode28857_bulk_mod.bands    Si-OQMD_5714-CollCode29287.bands
    C-OQMD_675640-CollCode28857_bulk_mod.castep   Si-OQMD_5714-CollCode29287.castep
    C-OQMD_675640-CollCode28857_bulk_mod.cell     Si-OQMD_5714-CollCode29287.cell
    C-OQMD_675640-CollCode28857_bulk_mod.cst_esp  Si-OQMD_5714-CollCode29287.cst_esp
    C-OQMD_675640-CollCode28857_bulk_mod.param    Si-OQMD_5714-CollCode29287.geom
    C-OQMD_675640-CollCode28857_bulk_mod.png      Si-OQMD_5714-CollCode29287.param
    C-OQMD_675640-CollCode28857_bulk_mod.res      Si-OQMD_5714-CollCode29287.res
    C-OQMD_675640-CollCode28857_bulk_mod.results  Si-OQMD_5714-CollCode29287_bulk_mod.bands
    Li-OQMD_30676-CollCode642106.bands            Si-OQMD_5714-CollCode29287_bulk_mod.castep
    Li-OQMD_30676-CollCode642106.castep           Si-OQMD_5714-CollCode29287_bulk_mod.cell
    Li-OQMD_30676-CollCode642106.cell             Si-OQMD_5714-CollCode29287_bulk_mod.cst_esp
    Li-OQMD_30676-CollCode642106.cst_esp          Si-OQMD_5714-CollCode29287_bulk_mod.param
    Li-OQMD_30676-CollCode642106.geom             Si-OQMD_5714-CollCode29287_bulk_mod.png
    Li-OQMD_30676-CollCode642106.param            Si-OQMD_5714-CollCode29287_bulk_mod.res
    Li-OQMD_30676-CollCode642106.res              Si-OQMD_5714-CollCode29287_bulk_mod.results
    Li-OQMD_30676-CollCode642106_bulk_mod.bands

You can see the results of the bulk modulus fits for each structure in
``completed/*_bulk_mod.results``, along with plots of the fits saved as pngs.

Compare these LDA values from the three different fits to the experimental values from Wikipedia:

+----------+--------------------------+------------------------------+
| Material | Expt. Bulk modulus (GPa) | Predicted Bulk Modulus (GPa) |
+==========+==========================+==============================+
| Lithium  | 11                       | 17.6, 17.7, 17.8             |
+----------+--------------------------+------------------------------+
| Silicon  | 98                       | 97.5, 98.0, 97.4             |
+----------+--------------------------+------------------------------+
| Diamond  | 442                      | 503, 496, 496                |
+----------+--------------------------+------------------------------+
.. index:: run3_qe

.. highlight:: bash

.. _run3_elastic:

Example 5: Geometry optimisations with Quantum Espresso and run3
================================================================

This tutorial uses the files found in ``examples/run3_quantum_espresso/vc-relax`` to relax
a folder of res files in a similar way to the CASTEP tutorial. Quantum Espresso (QE) uses a single
input file (as opposed to CASTEP's two), and these must be prepared beforehand (rather than
generated per structure with run3).

First things first, set up your directory so that you have:

* a load of res files
* a template QE input file containing DFT parameters

First, we make the QE input files using the ``shx3pwscf`` script:
``shx3pwscf *.res --template vcr.template --kpoint_spacing 0.03``

Then, to run the optimisations with run3 (either interactively, or at the bottom of a job script), run3 must be called on all the ``*.in`` files:
``run3 --redirect "$seed.out" -nc 4 --mode generic --executable 'pw.x -i $seed.in' *.in``
The ``$seed`` variables will be expanded by run3 to take the values of the ``.res`` file names. In this case, our only res file is called NaP.res, so the only calculation will be called as ``mpirun -n 4 pw.x -i NaP.in > NaP.out`` (or equivalent).
.. index:: run3_spectral

.. _run3_spectral:


Example 2: Spectral calculations with CASTEP and run3
-----------------------------------------------------

In this example, we will go from a crystal structure to a dispersion and DOS plot using run3, CASTEP and `OptaDOS <https://github.com/optados-developers/optados>`_. For this use case, run3 uses the `SeeK-path library <https://github.com/giovannipizzi/seekpath>`_ to generate standardised band paths through reciprocal space to automatically compute a useful bandstructure for all crystal types.

Spectral calculations follow a similar setup to geometry optimisations: run3 expects to find a folder containing .res files with lattice/atomic position data, and one .cell and one .param file specifying the CASTEP options. Standard run3 rules apply: if a ``<seed>.res.lock`` file is found, or if the .res file is listed in ``jobs.txt``, the structure will be skipped. Such a folder can be found in ``examples/bandstructure+dos/simple`` which contains some LiCoO\ :sub:`2` polymorphs. The Jupyter
notebook ``simple_spectral.ipynb`` will also show you exactly how to run a standard BS/DOS calculation and plot the results with the API. The files in ``examples/bandstructure+dos/projected`` will show you how to use matador and OptaDOS to get projected densities of states and bandstructures. Here, we shall run through some more simple cases first.

run3 will follow some simple rules to decide what kind of spectral calculation you want to run. First, it will check the ``task`` and ``spectral_task`` keywords. The ``task`` keyword needs to be set to ``'spectral'`` to trigger a spectral workflow. If ``spectral_task`` is either ``dos`` or ``bandstructure``, then this calculation type will always be included in the workflow, with default parameters if unset. Otherwise, run3 will check the cell file for the
``spectral_kpoints_mp_spacing`` (DOS) and ``spectral_kpoints_path_spacing`` (bandstructure) keywords and will perform either one or both of the corresponding calculation types. The first step will always be to perform an SCF calculation, which is continued from to obtain the DOS or BS (if any check file is found in the folder, run3 will attempt to skip the SCF and restart accordingly).

Both the bandstructure and DOS tasks can be post-processed with OptaDOS; run3 will perform this automatically if it finds a .odi file in alongside the .cell and .param. Similarly, if the ``write_orbitals_file`` keyword is set to True in the param file, ``orbitals2bands`` will be run automatically to reorder the eigenvalues in the .bands file. The output can be plotted using either the ``dispersion`` script bundled with matador, or through the API (see the Jupyter notebook).

Let us now consider some concrete examples.

Example 2.1: Bandstructure calculation with automatic path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can compute just bandstructure by specifying ``spectral_kpoints_path_spacing`` and ``spectral_task = bandstructure`` in the cell and param respectively, e.g. ::

   $ cat LiCoO2.param
   task: spectral
   spectral_task: bandstructure
   xc_functional: LDA
   cut_off_energy: 400 eV

   $ cat LiCoO2.cell
   kpoints_mp_spacing: 0.05
   spectral_kpoints_path_spacing: 0.05

with .res file ::

   $ cat LiCoO2-CollCode29225.res
   TITL LiCoO2-CollCode29225
   CELL 1.0  4.91 4.91 4.91 33.7 33.7 33.7
   LATT -1
   SFAC 	 Co Li O
   Co        1       0.000000        0.000000        0.000000         1.000000
   Li        2       0.500000        0.500000        0.500000         1.000000
   O         3       0.255463        0.255463        0.255463         1.000000
   O         3       0.744538        0.744538        0.744537         1.000000
   END

Simply calling ``run3 LiCoO2`` (with ``-v 4`` if you want to see what's going on) will perform three steps:

   1. Analyse the lattice with SeeKPath, perhaps creating a primitive cell or standardising if need be, and generate the standard Brillouin zone path with desired spacing set by ``spectral_kpoints_path_spacing``.
   2. Perform a self-consistent singlepoint energy calculation to obtain the electron density on the k-point grid specified by ``kpoints_mp_spacing``.
   3. Interpolate the desired bandstructure path, yielding a ``.bands`` file containing the bandstructure.

If you now run ``dispersion LiCoO2-CollCode29225 --png`` from the ``completed/`` folder, you should see this:

.. image:: LiCoO2-CollCode29225_spectral.png
   :name: bandstructure_only
   :align: center

Note that the path labels are generated from the .cell/.res file in the ``completed/`` folder, the .bands file does not contain enough information to make the entire plot portable. As mentioned above, if ``write_orbitals_file: true`` was found in the .param file, ``orbitals2bands`` would have been called to reorder the bands based on Kohn-Sham orbital overlaps.

.. tip::
   Plots can be customised using a custom matplotlib stylesheet. By default, matador plots will use the stylesheet found in ``matador/config/matador.mplstyle``, which can be copied elsewhere and specified by path in your ``.matadorrc`` under the ``plotting.default_style`` option. The default matplotlib styles can also be used directly by name, e.g. "dark_background".

.. tip::
   You can also specify a ``spectral_kpoint_path`` in your .cell file, that will be used for all structures in that
   folder. This can be useful when comparing the electronic structure of e.g. defected cells. The ``dispersion``
   script will try to scrape the labels to use for plotting from the comments at the specification of the kpoint path,
   e.g. ::

    %block spectral_kpoint_path
    0.0 0.0 0.0 ! \Gamma
    0.5 0.0 0.0 ! X
    %endblock spectral_kpoint_path

Example 2.2: A simple density of states (DOS) calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We shall now calculate the density of states in a similar way to the above, using the matador default for ``spectral_kpoints_mp_spacing``, such that the cell file only contains ``kpoints_mp_spacing: 0.05`` and the param file now has ``spectral_task: DOS``.

.. tip::
   If you are starting from the example above, you will need to move the .res file back into the current folder, delete the .txt files and rename the ``completed/`` folder to e.g. ``completed_bs/``. You might also consider moving the .check file to current folder, so that the calculation will re-use the SCF results.

Again, simply running ``run3 LiCoO2`` will do the trick. Eventually, a .bands_dos file will be produced (any existing .bands files will be backed up and replaced at the end of the run). The dispersion script will recognise this as a density of states, and will apply some naive gaussian smearing that can be controlled with the ``-gw/--gaussian_width`` flag. Running ``dispersion LiCoO2-CollCode29225 --png -gw 0.01`` will produce the following:

.. image:: LiCoO2-CollCode29225_spectral_dos.png
   :name: dos_only
   :align: center

.. tip::
   If you have a .bands file remaining in your top directory, ``dispersion`` will try to plot this as a bandstructure alongside your DOS, which may look terrible if .bands contains a DOS calculation! You can plot just the DOS using the ``--dos_only`` flag.


Example 2.3: Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a DOS and bandstructure on the same structure, simply include both ``spectral_kpoints_mp_spacing`` and ``spectral_kpoints_path_spacing`` in your .cell file. Your ``spectral_task`` keyword in the param file will be ignored. This exact example can be found in ``examples/bandstructure+dos/simple``, with an example Jupyter notebook showing how to make plots with the API directly, rather than the dispersion script.

After calling run3 again, the ``completed/`` folder in this case should contain both a .bands and a .bands_dos file which can be plotted alongside one another using ``dispersion LiCoO2-CollCode29225``, to produce the following:

.. image:: LiCoO2-CollCode29225_spectral_both.png
   :name: dos_bs
   :align: center

Example 2.4: Using OptaDOS for post-processing: projected DOS and bandstructures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final piece of the puzzle is `OptaDOS <https://github.com/optados-developers/optados>`_, a package for broadening and projecting densities of states (amongst other things) that comes with CASTEP. By default, run3 will turn on the required CASTEP settings (namely ``pdos_calculate_weights``) required by OptaDOS. In order for OptaDOS to be run automatically by run3, an extra .odi file must be added into our input deck, containing the details of the desired OptaDOS calculation.

.. note::
   This example assumes that the OptaDOS binary is called ``optados`` and resides in your PATH, likewise ``orbitals2bands``. This can altered by setting the ``run3.optados_executable`` setting in your matador config.

.. warning:: By default, OptaDOS will *not* be invoked with ``mpirun`` (i.e., your executable should work for serial runs too). A parallel OptaDOS run can be performed by setting the ``run3.optados_executable`` to e.g. ``mpirun optados.mpi`.

.. warning::
   The projected dispersion curve feature is quite new to OptaDOS and thus is temperamental. Depending on when you are reading this, it may require you to have compiled OptaDOS from the development branch on GitHub.

run3 will try to perform three types of calculation: a simple DOS smearing, a projected density of states (with projectors specified by the ``pdos`` keyword), and a projected bandstructure (with projectors specified by the ``pdispersion`` keyword). If ``pdos``/``pdispersion`` is not found in the .odi, this corresponding task will be skipped. Likewise, if ``broadening`` is not found in the .odi, the standard DOS broadening will not be performed.::

   $ cat LiCoO2.odi
   pdos: species_ang
   pdispersion: species
   adaptive_smearing: 1
   set_efermi_zero: True
   dos_per_volume: True
   broadening: adaptive
   dos_spacing: 0.01

With all these files in place, simply running ``run3 LiCoO2`` and ``dispersion (-interp 3 -scale 25) LiCoO2-CollCode29225`` (optional flags in brackets) should yield the following plot:

.. image:: LiCoO2-CollCode29225_spectral_pdis.png
   :name: full_spectral
   :align: center

.. tip:: Note that the colours of each projectors in these plots is set by your VESTA colour scheme, which is bundled by default inside ``matador/config``.
.. index:: run3_geom

.. _run3_geom:

Example 1: High-throughput geometry optimisations with CASTEP
-------------------------------------------------------------

.. _ex1:


Example 1.1: Using run3 locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In this example, we will suppose that you want to perform a geometry optimisation on several different polymorphs of TiO\ :sub:`2` from the ICSD. The files for this example can be found in ``examples/run3_tutorial``, `here <https://github.com/ml-evs/matador/blob/develop/examples/run3_tutorial>`_.

Setting up the files
^^^^^^^^^^^^^^^^^^^^

By default, run3 expects the following files in the job folder:

* one ``.res``/SHELX file `per structure` to optimise
* one ``$seed.cell`` file that contains the CASTEP CELL keywords common to every structure, e.g. pseudopotentials, kpoint spacing, etc.
* one ``$seed.param`` file that contains the CASTEP PARAM keywords common to every calculation, e.g. ``cut_off_energy``, ``xc_functional``, ``task``, ``geom_max_iter``, etc.
* any required pseudopotential files.

.. tip:: If you have a database set up with structures from the OQMD these could be obtained via ``matador query --db oqmd_1.1 -f TiO2 --icsd --res``.

.. tip:: Alternatively, you can turn many file types into ``.res`` using the various ``<format>3shx`` scripts (shx standing for SHELX), e.g. ``cell3shx *.cell``.

The job folder should look something like this::

    $ ls
    O2Ti-OQMD_112497-CollCode171670.res  O2Ti-OQMD_2575-CollCode9852.res     O2Ti-OQMD_7500-CollCode41493.res
    O2Ti-OQMD_117323-CollCode182578.res  O2Ti-OQMD_3070-CollCode15328.res    O2Ti-OQMD_84685-CollCode97008.res
    O2Ti-OQMD_13527-CollCode75179.res    O2Ti-OQMD_31247-CollCode657748.res  O2Ti-OQMD_97161-CollCode154036.res
    O2Ti-OQMD_19782-CollCode154035.res   O2Ti-OQMD_5979-CollCode31122.res    TiO2.cell
    O2Ti-OQMD_2475-CollCode9161.res      O2Ti-OQMD_7408-CollCode41056.res    TiO2.param

with ``.param`` file containing::

    $ cat TiO2.param
    task                 : geometryoptimization
    xc_functional        : LDA 
    cut_off_energy       : 300.0 eV
    geom_force_tol       : 0.1
    spin_polarized       : false
    fix_occupancy        : false
    max_scf_cycles       : 100
    opt_strategy         : speed
    page_wvfns           : 0
    perc_extra_bands     : 40
    num_dump_cycles      : 0
    backup_interval      : 0
    geom_method          : LBFGS
    geom_max_iter        : 300
    mix_history_length   : 20
    finite_basis_corr    : 0
    fixed_npw            : false
    write_cell_structure : true
    write_checkpoint     : none
    write_bib            : false
    bs_write_eigenvalues : false
    calculate_stress     : true

and ``.cell`` file containing::

    $ cat TiO2.cell
    kpoint_mp_spacing: 0.07
    
    %block species_pot
    QC5
    %endblock species_pot

    symmetry_generate
    symmetry_tol: 0.01
    snap_to_symmetry

.. highlight:: bash


Calling run3
^^^^^^^^^^^^

Once these files are in place, we can begin the geometry optimisations. To run the current host machine, simply call::

    $ run3 TiO2

This will start a single node CASTEP job on the current machine, using all available cores. If you are on a local cluster without a queuing system, and wish to run on several nodes at once (say ``node3``, ``node6`` and ``node8``), the oddjob script can be used as follows::

    $ oddjob 'run3 TiO2' -n 3 6 8

This will start 3 single node CASTEP jobs on the desired nodes. If instead your nodes are called ``cpu00010912``, ``cpu323232`` and ``cpu123123``, the ``--prefix`` flag is needed::

    $ oddjob 'run3 TiO2' --prefix cpu -n 00010912 323232 123123


.. tip:: On a supercomputer with a queuing system, e.g. PBS or slurm, run3 must be called in your submission script. Array jobs are typically an effective way of spreading out over multiple nodes. An example of this kind can be found in `example 1.2 <ex.1.2_>`__.


Monitoring your calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you look at the job folder as run3, er... runs, you will see several files and folders being created. Firstly, 3 ``.txt`` files will be made:

* ``jobs.txt``: this file contains a list of jobs that, at some point, __started__ running.
* ``finished_cleanly.txt``: this file lists jobs that completed without error.
* ``failures.txt``: this file lists jobs that produced an error.

Every structure in progress will have a ``<structure_name>.lock`` file to prevent clashes with other nodes.

Several folders will also be created:

* ``logs/``: log file per structure containing a complete history of the run.
* ``input/``: a backup of the starting configuration as a ``.res`` file.
* ``completed/``: all successful calculations will end up here, usually as a ``.res`` file with the final configuration, a concatenated ``.castep`` file containing every job step, and if requested (via ``write_cell_structure: true``), CASTEP's ``-out.cell`` file.
* ``bad_castep/``: all failed calculations end up here, including all auxiliary files.
* ``<node_name>/``: a folder is created per hostname (e.g. when running on multiple nodes) that contains the interim calculations. On failures/timeouts, all files in here are moved back to the main job folder.

Eventually, all jobs will hopefully be moved to ``completed/``, then you are done!


Example 1.1.1: High-throughput geometry optimisations with CASTEP with per-structure parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a few occasions where you might need a custom ``.param`` file for each structure, for example, if using the implicit nanotube ``%devel_code`` in CASTEP.

These calculations are performed in exactly the same was as above, except a ``<structure_name>.param`` file must be made containing the required DFT parameters AND the nanotube parameters. In this case, run3 must now be called as::

    $ run3 --custom_params TiO2

.. tip:: If you have a .res file that contains a PyAIRSS "REM NT_PROPS" line, this will be ignored.


Example 1.2: High-throughput geometry optimisations with CASTEP on a supercomputer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each HPC facility has its own quirks, so in this example we will try to be as explicit as possible. The set up of the job is exactly the same as in `example 1 <ex1_>`__, but we now must add run3 to our job submission script. The following examples are for the SLURM queuing system on the BlueBear machine at the University of Birmingham and PBS on ARCHER (Tier-1), but run3 has also been tested on CSD3 (Tier-2), HPC Midlands-Plus (Tier-2), Thomas (Tier-2) and several local group-scale clusters.

Example 1.2.1: SLURM on BlueBear
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this job, we will submit a run3 job that performs CASTEP calculations across 2 24-core nodes per structure. Let us presume we have many thousand structures to run. The submission script looks as follows::

    $ cat run3.sub
    #!/bin/bash -l
    
    ###### MACHINE/USER-SPECIFIC OPTIONS ######
    
    #SBATCH --ntasks 48
    #SBATCH --nodes 2-2
    #SBATCH --time 24:00:00
    #SBATCH --qos <REDACTED> 
    ##SBATCH --qos bbshort
    #SBATCH --mail-type ALL
    #SBATCH --account=<REDACTED>

    module purge
    export PATH="$HOME/bin/CASTEP-17.21:$HOME/.conda/bin"
    module load bluebear
    module load mpi/impi/2017.1.132-iccifort-2017.1.132
    unset I_MPI_PMI_LIBRARY
    
    # RUN3 COMMANDS
    # (assuming installation guide followed at
    #  https://matador-db.readthedocs.io/en/latest/install.html)

    source activate matador
    run3 -nc 48 -v 4 --executable castep.mpi --ignore_jobs_file TiO2

Let's unpick a few of the flags used to call run3 here:

* ``-nc/--ncores``: the number of cores to use per structure, per calculation. It is often worth specifying this if more than one node is being used, as the correctness of run3's core counter is queue/machine-specific.
* ``-v 4``: sets the verbosity in the log file to the highest level.
* ``--ignore_jobs_file``: by default run3 will for both ``<structure>.lock`` files and entries in ``jobs.txt`` before running a new structure. It is often worth disabling the ``jobs.txt`` check if it is not expected that all structures complete in one job submission (see below).
  
try to call an executable called simply ``castep``. On many machines, CASTEP is installed as ``castep.mpi``.

Now to submit this script as a 200-node array job (i.e. running a maximum of 100 structures concurrently, depending on the queue), we call the following::

    $ sbatch --array=1-100 run3.job

It may be that this job is not large enough to optimise all structures within the walltime limit. In this case, it can be resubmitted using the same command. Jobs that were running when the walltime ran out should automatically be pushed back into the job folder so that they will be available to the next run3 call. In the event that this does not happen (for example MPI kept control of the Python thread for too long so the queuieng system interrupted run3's clean up), ``<hostname>`` folder
will be left hanging around in the main jobs folder. Jobs must then be manually made restartable by deleting ``<structure>.lock`` (and removing ``<structure>>`` from ``jobs.txt`` if not using ``--ignore_jobs_file``). It may also be that the intermediate CASTEP calculation was not copied over from the ``<hostname>`` folder: in this case, the CASTEP files can be updated by running::
    
    $ cp -u node*/*.castep .

from inside the root job folder.

Example 1.2.2: PBS on ARCHER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ex.1.2:

Instructions are almost identical to the above, but the array job script looks a little different, for the same 100 copies of 2 node jobs (this time 24 cores per node)::

    $ cat run3.job
    #!/bin/bash --login
    # PBS job options (name, compute nodes, job time) # PBS -N is the job name (e.g. Example_MixedMode_Job)
    #PBS -N my_run3_job
    # PBS -l select is the number of nodes requested (e.g. 128 node=3072 cores)
    #PBS -l select=2
    # PBS -l walltime, maximum walltime allowed (e.g. 6 hours)
    #PBS -l walltime=24:00:00
    # Replace [budget code] below with your project code (e.g. t01)
    #PBS -A <REDACTED>
    #PBS -m abe
    #PBS -M <REDACTED>
    #PBS -J 1-100
    #PBS -r y

    # Make sure any symbolic links are resolved to absolute path
    export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

    # Change to the direcotry that the job was submitted from
    # (remember this should be on the /work filesystem)
    cd $PBS_O_WORKDIR

    source $HOME/.bashrc
    module load anaconda-compute/python3
    source activate $HOME/work/.conda/matador

    run3 --archer -v 4 -nc 48 KSnP

Notice here we have specified ``--archer``: again, run3 should be able to detect that ``mpirun`` is missing and thus try ``aprun``, but it can be worth specifying just in case. With PBS, the whole array can be submitted with just::

    $ qsub run3.job
.. index:: mongo

.. highlight:: bash

.. _mongo:

Setting up your own database
============================

``matador`` uses MongoDB to store and operate on results from first principles calculations. In this tutorial, we will set up a MongoDB database containing some dummy test data.

Installing and configuring MongoDB
----------------------------------

Instructions on how to install MongoDB community edition can be found at the `MongoDB website <https://docs.mongodb.com/manual/administration/install-community/>`_. If you are using Linux, you may find that your OS provides a MongoDB package that you can install.

Once installed, you can run your MongoDB server with the ``mongod`` command. An example configuration file can be found in ``matador/config/example_mongod.conf``, you may wish to change the port and default data directory (where the raw database is written to disk). This command will run persistently so you may wish to run it in a persistent terminal (using e.g. GNU screen or tmux), or as an OS-wide service (e.g. using systemctl or your OS's equivalent).

.. warning::

   It is up to the user to ensure that their MongoDB server is secure! This guide will assume access from ``localhost`` only, and thus will not set a password. If you do not wish external parties to have access to your database, make sure your server is running on secure/closed port, or set up authorisation as outlined on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_.

.. tip::

   You should ensure the path found under the ``storage.dbPath`` option in your configuration file exists before running ``mongod``, and that you have access to it!

Importing structures
--------------------

With a MongoDB server running, you can now import some structures. Navigate to ``examples/mongo/castep_files`` and run ``matador import``. If this is your first time using ``matador``, you should be prompted to make a configuration file. Follow the instructions and enter the values you used to set up your database. You will be given the option to test the connection to your database at the end.

If all is successful, you should see that 3 structures are added to your database, along with 4 failures. You can have a look at the spatula output files to see what happened. To finish this tutorial, use ``matador query`` to look at the structures you have just imported, before moving onto the command-line tutorial at :ref:`cli`.

Caring for your database
------------------------

Collections
^^^^^^^^^^^

MongoDB databases are made out of "collections" of data. To create a new collection, simply specify ``--db <collection_name>`` when performing an import. You may wish to create multiple collections of data, or one monolithic database.
   
The default collection can be set in your ``.matadorrc`` with the ``mongo.default_collection`` option. If you want to ensure that this database can only be imported to from a particular path, set the ``mongo.default_collection_file_path`` option.

A summary of the contents of existing collections can be generated using ``matador stats --db <collection_name>``. Similarly, a collection can be deleted using ``matador stats --db <collection_name> --delete``. If you wish to get a list of all existing collections, use ``matador stats --list``.

.. tip:: 
   In our group we combine these two approaches: a watched folder of "finalised" data is automatically scraped every day into a single monolithic collection (a cron job calls ``matador import``), interim calculations are then stored in collections per research project.

.. tip:: 
   Each collection also saves its own "changelog". You can view this changelog and undo changes with the ``matador changes`` interface.

Backing up
^^^^^^^^^^

You should always back up the data used to create a database, but if you wish, you can also directly backup the database itself using the MongoDB tool ``mongodump`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongodump/>`_). Similarly, restoration can be performed using ``mongorestore`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongorestore/>`_). You may also wish to read the general back up and restore tutorial on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/backup-and-restore-tools/>`_.
.. index:: run3_phonon

.. highlight:: bash

.. _run3_phonon:


Example 3: Phonon calculations with CASTEP and run3
-----------------------------------------------------

In this example, we will go from a crystal structure to phonon dispersion and DOS plot using run3 and CASTEP. As with electronic dispersion, run3 uses the `SeeK-path library <https://github.com/giovannipizzi/seekpath>`_ to generate standardised band paths through reciprocal space to automatically compute a useful dispersion for all crystal types.

Phonon calculations follow a similar procedure to the spectral calculations in Example 2: run3 expects to find a folder containing .res files with lattice/atomic position data, and one .cell and one .param file specifying the CASTEP options. Standard run3 rules apply: if a ``<seed>.res.lock`` file is found, or if the .res file is listed in ``jobs.txt``, the structure will be skipped. Such a folder can be found in ``examples/phonons/`` which contain structures of Si, Li and C. The Jupyter notebook ``phonon_results.ipynb`` will also show you exactly how to plot and analyse the results of the phonon calculations with the API. 

run3 will follow some simple rules to decide what kind of phonon calculation you want to run. First, it will check that ``task`` is set to ``phonon`` in the base param file. The workflow will then perform the following steps:

1. Relax the structure to the given ``geom_force_tol`` (which often needs to be very low for phonon calculations ~0.001 eV/A). If a ``.check`` file containing a relaxation is found, then this step is skipped.
2. From the structure/wavefunction found in the ``.check`` file computed in step 1 (or otherwise), compute the dynamical matrix with the given supercell/q-point spacing with CASTEP.
3. If ``phonon_fine_kpoint_path_spacing`` was specified in the base cell file (or a particular set of kpoints), then perform the dispersion calculation with CASTEP's built-in Fourier interpolation.
4. If ``phonon_fine_kpoint_mp_grid_spacing`` was specified in the base cell file, then perform a phonon DOS calculation with CASTEP's Fourier interpolation.

The output files will be cached in such a way that the ``.phonon`` file from the DOS calculation will not be overwritten by the other steps, and vice versa.

The specific example cell and param files will perform all of the above workflow steps: ::

   $ cat phonon.cell
   kpoints_mp_spacing: 0.07
   phonon_kpoints_mp_spacing: 0.1
   phonon_fine_kpoints_path_spacing: 0.05
   phonon_fine_kpoints_mp_spacing: 0.05
   snap_to_symmetry
   symmetry_tol: 0.01
   symmetry_generate
   %block species_pot
   QC5
   %endblock species_pot

   $ cat phonon.param
   task                          : phonon 
   phonon_method                 : finitedisplacement
   phonon_fine_method            : interpolate
   cut_off_energy                : 300.0 eV
   write_bib                     : False
   xc_functional                 : LDA
   geom_force_tol                : 0.05
   finite_basis_corr             : 0
   write_checkpoint              : none

The results can then be plotted with the ``dispersion`` script, if supplied with the ``-ph/--phonons`` flag, e.g. ::
   
   $ dispersion -ph completed/C

.. image:: C-OQMD_675640-CollCode28857.png
   :name: C phonons 
   :align: center
.. _contributing:

Contributing
============

Contributions and suggestions to this package are very welcome.

If you think you have found a bug, please raise an issue on GitHub, providing information about what you were trying to do, the function/script you ran, the error message/output and your ``matador`` version. If you are able to, please try to replicate the problem on the ``master`` branch before posting.

If you have a feature request, or you would like to contribute a feature, please raise an issue so that its suitability for this package can be discussed. Ideally, you could do this before writing your implementation, as the feature may already exist! Code contributions must include tests for any new functionality and must pass existing tests and ``flake8`` linting run by the CI. Any new modules must use the ``black`` auto-formatter.
.. index:: cli

.. _cli:

Command-line usage
==================

Many of ``matador``'s entry points will make use of any ``.matadorrc`` config file they find. Customisation options can be found under :ref:`getting_started`.

Main entry points
-----------------

``matador`` has 3 main CLI entry points that are added to your ``$PATH`` on installation:

- ``matador``: used to interact with databases.
- ``run3``: used to run calculations at high-throughput on the local machine.
- ``dispersion``: used to plot electronic or vibrational properties from local files containing the output of the relevant calculations.

Each of these has many options/sub-commands, which can be explored using ``--help``.

Script entry points
-------------------

There are also several convenience scripts bundled with ``matador``. These scripts are not as well-supported as the main API, but please raise an issue on GitHub if they do not behave as expected. These scripts are also added to your ``$PATH`` on installation and can be found in the ``scripts/`` folder of a development install. If you have a workflow that only requires simple script usage but is not yet supported, or if you have a script you would like to contribute, please raise an issue on GitHub.

Analysis scripts
~~~~~~~~~~~~~~~~

- ``compare_pdfs``: compute PDF overlaps and plot differences from ``.res`` file input.
- ``pxrd_calculator``: compute, plot and export PXRD patterns from ``.cif`` file input.
- ``dbtools``: simple script for dropping MongoDB collections created by ``matador``.
- ``oddjob``: a script for spawning jobs on compute nodes with no queuing and shared filesystems, e.g. ``oddjob 'run3 <seed>' -n 1 3 5 -p cpu-`` will spawn ``run3`` processes on compute nodes with names ``cpu-1``, ``cpu-3`` and ``cpu-5``.
- ``fryan``: a stripped down version of ``matador query`` that operates on folders of ``.res`` files.
- ``plot_convergence``: plots the output of a ``run3`` convergence run.
- ``standardize_cell``: for an input ``.res`` file, standardize the cell (optionally to primitive) with spglib and write out a new file.

File converter scripts
~~~~~~~~~~~~~~~~~~~~~~

These scripts all use ``matador``'s internal file readers/writers to operate on lists of files provided at the CLI: e.g. ``x3y *.x`` will read all ``.x`` files in  a folder and output a series of ``.y`` files. This is not an all-to-all process, these scripts have only been written when needed. If you need quick access to a file converter, either write your own script based on the below, or request the converter script in an issue on GitHub.

- ``atoms3shx``: use ASE's file reader on any arbitrary file, and output a ``.res`` file.
- ``castep3shx``, ``castep3cell``: convert a CASTEP output file into ``.res`` or ``.cell`` respectively.
- ``cell3shx``, ``cif3shx``, ``magres3shx``: convert a CASTEP ``.cell``, ``.cif``, ``.magres`` or Quantum Espresso output file into a ``.res`` file.
- ``shx3cell``, ``shx3cif``, ``shx3pwscf``: convert a ``.res`` file into either ``.cell``, ``.cif``, or Quantum Espresso ``.in`` input file.
.. index:: examples_index

.. _examples_index:

Tutorials and example notebooks
===============================

Real use-cases
--------------

- `Harper et al, Research Data Accompanying "Computational Discovery of Novel Copper Phosphide Conversion Anode for Lithium-Ion Batteries" <https://github.com/harpaf13/data.copper-phosphides>`_. |CuP Binder|

Example Jupyter notebooks
-------------------------

Interactive
~~~~~~~~~~~

These interactive examples can be run under Binder; some features that required system configuration (e.g. fonts when plotting) may not work correctly.

|Binder|

.. toctree::
    :maxdepth: 1

    notebooks/interactive/matador_plot_styles
    notebooks/interactive/voltage_from_res_files
    notebooks/interactive/pxrd_plotting
    notebooks/interactive/magres_plotting
    notebooks/interactive/pymatgen_and_ase_interface

Non-interactive
~~~~~~~~~~~~~~~

These examples require external data, but can be used as example code.

.. toctree::
    :maxdepth: 1

    notebooks/non-interactive/hull
    notebooks/non-interactive/phonon_results
    notebooks/non-interactive/projected_spectral
    notebooks/non-interactive/projected_spectral_as_subplots

Tutorials
---------

These tutorials will guide you through common use cases of the command-line tools bundled with ``matador``.

.. toctree::
    :maxdepth: 2

    tutorials/run3
    tutorials/mongo

.. |Binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/ml-evs/matador/master?filepath=examples/interactive/

.. |CuP Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/harpaf13/data.copper-phosphides/master?filepath=CuP_results.ipynb
=======
matador
=======

| |PyPI Version| |GH Actions| |Binder|
| |Documentation Status| |MIT License| |Coverage Status|
| |JOSS| |Zenodo|


matador is an aggregator, manipulator and runner of first-principles calculations, written with a bent towards battery electrode materials.
The source can be found on `GitHub <https://github.com/ml-evs/matador>`_ and online documentation is hosted at `ReadTheDocs <https://docs.matador.science>`_.

Example Jupyter notebooks and tutorials can be found `online <https://docs.matador.science/en/latest/examples_index.html>`_ or in the ``examples/`` folder of the matador source code.

Written & maintained by `Matthew Evans <https://ml-evs.science>`_ (2016-).


.. image:: docs/src/img/lipzn.png
   :name: LiPZn
   :align: center

Installation
------------

In the simplest case (e.g. you already have Python 3.6+ set up), ``pip install matador-db`` is sufficient to get up and running, preferably in a fresh virtual environment.

Upgrading to the latest version should be as simple as ``pip install -U matador-db``.

For an editable development installation, clone the source code from this repository and run ``pip install -e .`` from the matador folder. Tests can be run on your local machine with ``python -m unittest discover -v -b`` or simply with ``py.test`` after test dependencies have been installed with ``pip install .[test]``. In order to test CASTEP-running functionality, the tests will look for an MPI-enabled executable named ``castep`` on your ``$PATH``.

Further instructions can be found in the `Installation instructions <https://docs.matador.science/en/latest/install.html>`_.


Usage
------

``matador`` is primarily a Python *library* that can be used inside Python scripts/modules to create a custom workflow. There are, however, several command-line scripts bundled with ``matador`` itself. All of these scripts are listed under `CLI Usage <https://docs.matador.science/en/latest/cli.html>`_.

For basic command-line usage, please explore the help system for each command. Common workflows can be found inside ``examples/`` and in the `online docs <http://docs.matador.science/en/latest/examples_index.html>`_.

Please consult the full `Python API documentation <http://docs.matador.science/en/latest/modules.html>`_ for programmatic usage.

Core functionality
-------------------

The API has many features that can be explored in the examples and API documentation. As a summary, ``matador`` can be used for:

- Scraping of CASTEP (and Quantum Espresso) input/output files into flexible Python dictionaries/models.
- The creation and curation of MongoDB collections of geometry optimisation calculations, with a powerful querying CLI/API.
- Customisable, publication-ready plots for all models, e.g. phase diagrams, PDF, PXRD, voltage profiles, electronic/vibrational bandstructures etc.
- High-throughput geometry optimisations, electronic and vibrational properties using CASTEP (and Quantum Espresso) with ``run3``. Tested on several supercomputers. ``run3`` is designed primarily for simple workflows and offers little in the way of tools for creating complex workflows; if this is your use case, then check out some of the other codes listed below.
- Creation of phase diagrams and electrochemical voltage profiles from the results of DFT calculations.

This functionality is achieved by interfacing with much of the standard scientific Python stack (`NumPy <https://numpy.org>`_, `SciPy <https://scipy.org>`_, `matplotlib <https://matplotlib.org>`_), some more materials-specific packages (`spglib <https://github.com/atztogo/spglib/>`_, `SeeK-path <https://github.com/giovannipizzi/seekpath>`_, `periodictable <https://github.com/pkienzle/periodictable>`_) and other general packages (`pymongo <https://github.com/mongodb/mongo-python-driver>`_, `python-ternary <https://github.com/marcharper/python-ternary>`_, `numba <https://numba.org>`_).

Similar packages
----------------

This package is by no means unique in its functionality or goals. Below is a list of similar packages and an overview of where they overlap with ``matador``:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_: manipulation of structures, parsing and exporting files, running jobs and local databases with ``ase-db``. An interface is provided to ASE's ``Atoms`` object.
- `pymatgen <https://pymatgen.org>`_: similar to ASE, with a focus on the VASP DFT code. An interface is provided to pymatgen's ``Structure`` object. Local databases can be constructed with the `pymatgen-db <https://github.com/materialsproject/pymatgen-db>`_ add-on and high-throughput workflows are achieved with `Fireworks <https://github.com/materialsproject/fireworks>`_.
- `AiiDA <https://www.aiida.net>`_: high-throughput job running, provenance tracking and database storage for many simulation codes.
- `sumo <https://github.com/SMTG-UCL/sumo>`_ publication quality plotting, primarily for VASP but also with support for other codes.

If you think this list is outdated, incorrect or simply incomplete, then please raise an issue!

Citing matador
--------------

If you use matador in your work, we kindly ask that you cite

    Matthew L. Evans, Andrew J. Morris, *matador: a Python library for analysing, curating and performing high-throughput density-functional theory calculations* Journal of Open Source Software, 5(54), 2563 (2020)
    `10.21105/joss.02563 <https://doi.org/10.21105/joss.02563>`_

Source code archives for all versions above 0.9 can be found at Zenodo `DOI 10.5281/zenodo.3908573 <https://doi.org/10.5281/zenodo.3908573>`_.


.. |PyPI Version| image:: https://img.shields.io/pypi/v/matador-db?label=PyPI&logo=pypi
   :target: https://pypi.org/project/matador-db/
.. |GH Actions| image:: https://img.shields.io/github/workflow/status/ml-evs/matador/Run%20tests/master?label=master&logo=github
   :target: https://github.com/ml-evs/matador/actions?query=branch%3Amaster
.. |MIT License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://github.com/ml-evs/matador/blob/master/LICENSE
.. |Coverage Status| image:: https://img.shields.io/codecov/c/gh/ml-evs/matador/master?logo=codecov
   :target: https://codecov.io/gh/ml-evs/matador
.. |Documentation Status| image:: https://readthedocs.org/projects/matador-db/badge/?version=stable
   :target: https://matador-db.readthedocs.io/en/stable/?badge=stable
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3908573.svg
   :target: https://doi.org/10.5281/zenodo.3908573
.. |Binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/ml-evs/matador/master?filepath=examples/interactive
.. |JOSS| image:: https://joss.theoj.org/papers/4d0eea9bea4362dec4cb6d62ebccc913/status.svg
   :target: https://joss.theoj.org/papers/4d0eea9bea4362dec4cb6d62ebccc913
Getting Started
===============

After successfully installing from the guide found at :ref:`install`, we can now try to do some science!

When running matador from the command line for the first time, you will be prompted to create a `.matadorrc` config file. This allows you to point to the MongoDB instance you want to use, or set up default plotting options. Alternatively, you can create your own `.matadorrc` based on the default found under ``matador/config/matadorrc.yml`` on `GitHub <https://github.com/ml-evs/matador/blob/master/matador/config/matadorrc.yml>`_.

For API examples, you should now be able to follow/run the Jupyter notebooks found in the ``examples/``. Some of them will require a MongoDB server accessible from your machine, or a CASTEP installation. For some details on how to set this up, please look at the :ref:`MongoDB tutorial <mongo>`. For command-line usage, tutorials should start appearing in due course over at :ref:`cli`.
.. _modules:

Modules
=======

.. toctree::
   :maxdepth: 6

   matador
.. include:: readme.rst

.. toctree::
    :hidden:
    :maxdepth: 4

    install
    getting_started
    changelog
    examples_index
    cli
    Python API <modules>
    contributing
.. _install:

Installation
============

If you have any issues with installation, feel free to raise an issue on GitHub outlining your approach and any errors you received.


Simple installation with pip
----------------------------

The matador package can be found on PyPI under the name `matador-db <https://pypi.org/project/matador-db>`_ and installed with
``pip install matador-db``, preferably in a fresh virtual environment (see conda instructions below). Extra dependencies may be installed with e.g. ``pip install matador-db[all]``.

Development installation with conda/pip
---------------------------------------

The tl;dr way to install matador, on e.g. a computing cluster, is as follows:

1. Clone the matador source onto your local machine ``git clone https://github.com/ml-evs/matador.git``.

Optional (but recommended) steps:

2. `Install conda <https://conda.io/miniconda.html>`_, if you have not already. There may be a package available already if you are using a supercomputer (e.g. `anaconda-compute/2.2.0-python3` on ARCHER 30/10/2017).
3. Create a new conda environment to install matador into (``conda create -n matador python=3.7``) and activate it with (``conda activate matador``).
4. Install some of the heavier requirements (e.g. NumPy and SciPy) through conda with ``conda install --yes --file requirements/requirements.txt``.

Required steps:

5. Run ``pip install .`` from inside the top-level matador directory, or ``pip install -e .`` for an editable developer install.
6. You now have a basic matador API installation, if you wish to use all matador features, install extra dependencies from the other requirements files inside ``requirements/`` using either conda or pip. If you wish to just install everything use ``pip install .[all]``.
7. To use matador, you will need to activate the conda environment from step 2, by running ``conda activate matador``. You will also need this in e.g. any job scripts. After installing the test dependencies with ``pip install .[test]``, you can test your installation using ``python -m unittest discover -v -b`` or simply ``py.test``. By default this will look for an MPI-enabled executable called ``castep`` on your ``$PATH`` to run CASTEP tests.

Troubleshooting
---------------

Below are some problems encountered on various machines that may be helpful:

1. (10/09/2019) When installing with ``conda``, if you receive the following error (or
   similar): ``/home/#####/.local/conda/envs/matador/compiler_compat/ld: build/temp.linux-x86_64-3.6/psutil/_psutil_common.o: unable to initialize decompress status for section .debug_info``, then you are using a modern compiler that breaks ``conda``'s attempts to be backwards compatible (in this case it was GCC 9). The simple fix is to rename/remove the copy of ``ld`` inside your conda environment (path in the message above) such that your system ``ld`` is used.
2. (10/10/2017) On some machines (e.g. ARCHER/Thomas) you may receive permissions errors at step 5; if so, try moving matador's `.git` and install again (``mv .git $HOME/matador_git_stash; pip install . ; mv $HOME/matador_git_stash .git``).
3. Some dependencies may not have compiled packages (wheels) for your distribution on PyPI, and may have compilation errors. In this case, you could try finding the relevant package on conda instead.
.. _changelog:

Changelog
=========

New in release (0.9.11) [03/06/2021]
------------------------------------

- Minor change: allow specification of arbitrary strings for CASTEP pseudopotential libraries (#156).
- Bug fix: ``standardize_cell`` script failing to use default symmetry tolerance value (#157).
- Bug fix: scraping of .cif files with single atoms and no symmetries (#173)
- Bug fix: scraping of Hubbard U values from .castep files, and associated bugs when performing relaxations with Hubbard U (#180)
- Dependency updates and Python 3.6 deprecation warning (#158, #181)

New in release (0.9.10) [23/02/2021]
------------------------------------

- Windows compatibility changes (#149)
- Dependency updates (#146, #148, #149)

New in release (0.9.9) [16/10/2020]
-----------------------------------

- Added support for CASTEP kpoint path ``BREAK`` directive (#107)
- Improvements to magres plotting and magres workflow (#112)
- Added ability to scrape electric field gradient values and compute quadrupolar quantities from NMR calculations (#115)
- Added ability to run all several examples under Binder (#106, #130).
- JOSS paper accepted! (#129)


New in release (0.9.8) [10/08/2020]
-----------------------------------
- Improvements to PDIS functionality (#94).

  - Rasterized scatter points for more efficient exporting and fewer graphical artifacts
  - Made underlying :func:`matador.plotting.spectral_plotting.dos_plot` and :func:`matador.plotting.spectral_plotting.dispersion_plot` more API friendly, and added example notebook.
  - Fixed bug in cell scraping for old ``BS_*`` style keywords.

- Improvements to magres functionality, including scraping of units (#90)
- Example notebooks that do not need external data/databases are now run as part of CI (#91).
- New workflow for NMR calculations and refactoring of old workflows (#96).

  - New workflow performs relaxation and high-quality SCF before NMR calculation.
  - Old workflows refactored and improved to enforce certain required parameters for e.g. checkpointing.
  - Enabled phonon workflow for CASTEP ``PHONON+EFIELD`` task.
  - Made file scrapers less dependent on file type.

- Updated CASTEP parameter list to 20.1 (#97).
- Tweaked spectral plotting defaults, including ``--colours`` flag to dispersion script (#98).


New in release (0.9.7) [29/07/2020]
-----------------------------------
- Bug fixes to problems introduced in 0.9.6.
- Cosmetic fixes to logging and misleading status reports in workflows.


New in release (0.9.6) [28/07/2020]
-----------------------------------
- Improvements to ASE and pymatgen interoperability (#80)
- Fixed bug in :class:`matador.hull.TemperatureDependentHull` which would crash when not provided a list of temperatures (#82).
- Added plotting functions for magres data, and improved its handling inside :class:`matador.crystal.Crystal` (#79).

New in release (0.9.5) [25/06/2020]
-----------------------------------
- This release is mostly to trigger Zenodo archiving.
- Updated README and tests for recent Python versions.


New in release (0.9.4) [08/06/2020]
-----------------------------------
- Fixed flag help strings for ``pxrd_calculator`` (#65)
- Changed default PDF broadening for 3x speedup (#65)
- Reverted ``cpu_count`` to use version that works correctly in most cases, by chance (#66).


New in release (0.9.3) [07/06/2020]
-----------------------------------

- Fixes for the CIF reader: now works with awkward linebreaks and alternative symmetry operation specifications (#61).
- Added several new flags to ``pxrd_calculator`` script (#60 and 61).
- Usability fixes for ``spectral_plotting`` in the case of projected dispersion curves (#59).


New in release (0.9.2) [01/06/2020]
-----------------------------------

- Optimised CIF reader considerably (#50)
- Updated PXRD calculator to allow for partial occupancy, monochromated beam angles and text export, and added ``pxrd_calculator`` script for convenience when handling CIF files.
- Added ability to choose which projectors are plotted with dispersion (#47)
- Various minor fixes and updates:

  - Updates to docs for CLI and configuration.
  - Allow nan-values to be reset inside :class:`matador.crystal.Crystal`.
  - Fixed display ordering of fingerprint-filtered cursors.


New in release (0.9.1) [20/05/2020]
-----------------------------------

- Fixed issue with local pip installs after 0.9 release
- Fixed issue with multi-node MPI tasks by switching to ``proc.communicate()`` after an initial polling stage (#37)
- Fixed issue where bands would be reordered multiple times in spectral plots (#40)
- Tweaked spectral plot defaults (#40)
- Replaced ``multiprocessing.cpu_count()`` calls with ``psutil.cpu_count(logical=False)`` to avoid double-counting hyperthreaded cores


New in release (0.9) [15/05/2020]
---------------------------------

- PyPI release! Can now install with ``pip install matador-db`` (unfortunately ``matador`` was taken, but they are sufficiently orthogonal that the package name ``matador`` is retained here.
- Much improved code structure and many additional classes that wrap the raw calculation dictionaries for e.g. :class:`matador.crystal.Crystal` and spectral classes.
- New module :mod:`matador.orm` containing useful models for data handling.

  - :class:`matador.orm.orm.DataContainer` as a base class that allows for easy
    access to underlying dictionaries.
  - :mod:`matador.orm.spectral` module that contains many useful classes for
    manipulating and plotting e.g. bandstructures, DOS and finite temperature
    properties.

- New features in :mod:`matador.hull` module:

  - Pseudo-ternary phase diagrams (building towards arbitrary n-dimensional phase diagrams).
  - :class:`matador.hull.EnsembleHull` class and submodule to support the Bayesian Error Estimate Functional (BEEF) and finite temperature phase diagrams.
  - Refactoring of hull calculation into light-weight :class:`matador.hull.PhaseDiagram` class.
  - Finite temperature hulls based on :class:`matador.hull.EnsembleHull` with
    :class:`matador.hull.TemperatureDependentHull`.

- Refactored old PDF `similarity` module into new module :mod:`matador.fingerprints`.

  - Added new fingerprint class, :class:`matador.fingerprints.PXRD`, with associated plots (thanks for James Darby for some initial code). Defaults calibrated with GSAS-II.
  - :class:`matador.fingerprints.PDF` sped up by an order of magnitude using `numba`.

- :class:`matador.workflows.castep.CastepSpectralWorkflow` extended to include latest projected dispersion curve developments from OptaDOS, with associated projected dispersion plots (see tutorial).

  - Updated dispersion script to automatically perform naive Gaussian smearing if OptaDOS output not detected.

- Abstracted and simplified :mod:`matador.compute` module to allow for extension to new codes via :mod:`matador.compute.calculators` submodule.

  - Should now be more robust and transferrable, with many HPC environments automatically detected.
  - Added ``--scratch_prefix`` to run3 to allow for temporary files to e.g. be written to faster filesystem with appropriate symlinks to work folder.

- All CASTEP 19 keywords supported, as well as `devel_code` blocks.
- Several new tests: coverage now around 75% when CASTEP is available.
- New tutorials:

  - :ref:`MongoDB setup<mongo>`
  - :ref:`Spectral calculations with run3<run3_spectral>`
  - Example notebooks


New in release (0.8b) [03/08/2018]
----------------------------------

- Wholesale changes, complete refactoring of most of the code.
- Released open source under the MIT license!
- Documentation now hosted on `readthedocs <matador-db.readthedocs.org>`_,
- Workflows: chaining up job steps with run3:

  - spectral and phonons (combined DOS, dispersion calculations) with automated kpoint paths.
  - bulk modulus calculations and EOS fitting.

- New tutorials:

  - :ref:`Geometry optimisations with run3<run3_geom>`

- Temperature-dependent convex hulls (thanks to Angela Harper).
- New per-used configuration that allows changing of plotting styles, colourschemes, database names etc.
- Improvements to compute module:

  - automatically handle walltime constraints for Slurm and PBS.
  - estimate memory usage with CASTEP and skip if exceeds machine capacity,

- All CASTEP 18 keywords supported.
- Better support for electronic structure data, OptaDOS, NMR calculations, CIF files, partial occupancy.


New in version (0.7b) [13/04/2017]
----------------------------------

-  Ternary voltage curves.
-  Similarity/uniqueness filtering with element-projected PDFs.
-  Updated compute engine for remote calculations (see ``compute.py`` and new script ``oddjob``).
-  Improved test suite and full pip compatiblity.
-  Many bugfixes and usability changes.

New in version (0.6b) [01/06/2017]
----------------------------------

-  Intercalation voltage curves, e.g. ``matador voltage -c Li:SnS2``.
-  Ternary phase diagrams with heatmaps for structure prediction sampling, gravimetric capacity and formation enthalpy ``matador hull -c ABC --sampmap --efmap --capmap``.
-  Substructural similarity interface with Can Kocer's code, as proposed by `Yang et al., PRB (2014) <http://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.054102>`_.
.. index:: run3_phonon

.. highlight:: bash

.. _run3_phonon:


Example 3: Phonon calculations with CASTEP and run3
-----------------------------------------------------

In this example, we will go from a crystal structure to phonon dispersion and DOS plot using run3 and CASTEP. As with electronic dispersion, run3 uses the `SeeK-path library <https://github.com/giovannipizzi/seekpath>`_ to generate standardised band paths through reciprocal space to automatically compute a useful dispersion for all crystal types.

Phonon calculations follow a similar procedure to the spectral calculations in Example 2: run3 expects to find a folder containing .res files with lattice/atomic position data, and one .cell and one .param file specifying the CASTEP options. Standard run3 rules apply: if a ``<seed>.res.lock`` file is found, or if the .res file is listed in ``jobs.txt``, the structure will be skipped. Such a folder can be found in ``examples/phonons/`` which contain structures of Si, Li and C. The Jupyter notebook ``phonon_results.ipynb`` will also show you exactly how to plot and analyse the results of the phonon calculations with the API. 

run3 will follow some simple rules to decide what kind of phonon calculation you want to run. First, it will check that ``task`` is set to ``phonon`` in the base param file. The workflow will then perform the following steps:

1. Relax the structure to the given ``geom_force_tol`` (which often needs to be very low for phonon calculations ~0.001 eV/A). If a ``.check`` file containing a relaxation is found, then this step is skipped.
2. From the structure/wavefunction found in the ``.check`` file computed in step 1 (or otherwise), compute the dynamical matrix with the given supercell/q-point spacing with CASTEP.
3. If ``phonon_fine_kpoint_path_spacing`` was specified in the base cell file (or a particular set of kpoints), then perform the dispersion calculation with CASTEP's built-in Fourier interpolation.
4. If ``phonon_fine_kpoint_mp_grid_spacing`` was specified in the base cell file, then perform a phonon DOS calculation with CASTEP's Fourier interpolation.

The output files will be cached in such a way that the ``.phonon`` file from the DOS calculation will not be overwritten by the other steps, and vice versa.

The specific example cell and param files will perform all of the above workflow steps: ::

   $ cat phonon.cell
   kpoints_mp_spacing: 0.07
   phonon_kpoints_mp_spacing: 0.1
   phonon_fine_kpoints_path_spacing: 0.05
   phonon_fine_kpoints_mp_spacing: 0.05
   snap_to_symmetry
   symmetry_tol: 0.01
   symmetry_generate
   %block species_pot
   QC5
   %endblock species_pot

   $ cat phonon.param
   task                          : phonon 
   phonon_method                 : finitedisplacement
   phonon_fine_method            : interpolate
   cut_off_energy                : 300.0 eV
   write_bib                     : False
   xc_functional                 : LDA
   geom_force_tol                : 0.05
   finite_basis_corr             : 0
   write_checkpoint              : none

The results can then be plotted with the ``dispersion`` script, if supplied with the ``-ph/--phonons`` flag, e.g. ::
   
   $ dispersion -ph completed/C

.. image:: C-OQMD_675640-CollCode28857.png
   :name: C phonons 
   :align: center
.. index:: run3_geom

.. _run3_geom:

Example 1: High-throughput geometry optimisations with CASTEP
-------------------------------------------------------------

.. _ex1:


Example 1.1: Using run3 locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


In this example, we will suppose that you want to perform a geometry optimisation on several different polymorphs of TiO\ :sub:`2` from the ICSD. The files for this example can be found in ``examples/run3_tutorial``, `here <https://github.com/ml-evs/matador/blob/develop/examples/run3_tutorial>`_.

Setting up the files
^^^^^^^^^^^^^^^^^^^^

By default, run3 expects the following files in the job folder:

* one ``.res``/SHELX file `per structure` to optimise
* one ``$seed.cell`` file that contains the CASTEP CELL keywords common to every structure, e.g. pseudopotentials, kpoint spacing, etc.
* one ``$seed.param`` file that contains the CASTEP PARAM keywords common to every calculation, e.g. ``cut_off_energy``, ``xc_functional``, ``task``, ``geom_max_iter``, etc.
* any required pseudopotential files.

.. tip:: If you have a database set up with structures from the OQMD these could be obtained via ``matador query --db oqmd_1.1 -f TiO2 --icsd --res``.

.. tip:: Alternatively, you can turn many file types into ``.res`` using the various ``<format>3shx`` scripts (shx standing for SHELX), e.g. ``cell3shx *.cell``.

The job folder should look something like this::

    $ ls
    O2Ti-OQMD_112497-CollCode171670.res  O2Ti-OQMD_2575-CollCode9852.res     O2Ti-OQMD_7500-CollCode41493.res
    O2Ti-OQMD_117323-CollCode182578.res  O2Ti-OQMD_3070-CollCode15328.res    O2Ti-OQMD_84685-CollCode97008.res
    O2Ti-OQMD_13527-CollCode75179.res    O2Ti-OQMD_31247-CollCode657748.res  O2Ti-OQMD_97161-CollCode154036.res
    O2Ti-OQMD_19782-CollCode154035.res   O2Ti-OQMD_5979-CollCode31122.res    TiO2.cell
    O2Ti-OQMD_2475-CollCode9161.res      O2Ti-OQMD_7408-CollCode41056.res    TiO2.param

with ``.param`` file containing::

    $ cat TiO2.param
    task                 : geometryoptimization
    xc_functional        : LDA 
    cut_off_energy       : 300.0 eV
    geom_force_tol       : 0.1
    spin_polarized       : false
    fix_occupancy        : false
    max_scf_cycles       : 100
    opt_strategy         : speed
    page_wvfns           : 0
    perc_extra_bands     : 40
    num_dump_cycles      : 0
    backup_interval      : 0
    geom_method          : LBFGS
    geom_max_iter        : 300
    mix_history_length   : 20
    finite_basis_corr    : 0
    fixed_npw            : false
    write_cell_structure : true
    write_checkpoint     : none
    write_bib            : false
    bs_write_eigenvalues : false
    calculate_stress     : true

and ``.cell`` file containing::

    $ cat TiO2.cell
    kpoint_mp_spacing: 0.07
    
    %block species_pot
    QC5
    %endblock species_pot

    symmetry_generate
    symmetry_tol: 0.01
    snap_to_symmetry

.. highlight:: bash


Calling run3
^^^^^^^^^^^^

Once these files are in place, we can begin the geometry optimisations. To run the current host machine, simply call::

    $ run3 TiO2

This will start a single node CASTEP job on the current machine, using all available cores. If you are on a local cluster without a queuing system, and wish to run on several nodes at once (say ``node3``, ``node6`` and ``node8``), the oddjob script can be used as follows::

    $ oddjob 'run3 TiO2' -n 3 6 8

This will start 3 single node CASTEP jobs on the desired nodes. If instead your nodes are called ``cpu00010912``, ``cpu323232`` and ``cpu123123``, the ``--prefix`` flag is needed::

    $ oddjob 'run3 TiO2' --prefix cpu -n 00010912 323232 123123


.. tip:: On a supercomputer with a queuing system, e.g. PBS or slurm, run3 must be called in your submission script. Array jobs are typically an effective way of spreading out over multiple nodes. An example of this kind can be found in `example 1.2 <ex.1.2_>`__.


Monitoring your calculations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you look at the job folder as run3, er... runs, you will see several files and folders being created. Firstly, 3 ``.txt`` files will be made:

* ``jobs.txt``: this file contains a list of jobs that, at some point, __started__ running.
* ``finished_cleanly.txt``: this file lists jobs that completed without error.
* ``failures.txt``: this file lists jobs that produced an error.

Every structure in progress will have a ``<structure_name>.lock`` file to prevent clashes with other nodes.

Several folders will also be created:

* ``logs/``: log file per structure containing a complete history of the run.
* ``input/``: a backup of the starting configuration as a ``.res`` file.
* ``completed/``: all successful calculations will end up here, usually as a ``.res`` file with the final configuration, a concatenated ``.castep`` file containing every job step, and if requested (via ``write_cell_structure: true``), CASTEP's ``-out.cell`` file.
* ``bad_castep/``: all failed calculations end up here, including all auxiliary files.
* ``<node_name>/``: a folder is created per hostname (e.g. when running on multiple nodes) that contains the interim calculations. On failures/timeouts, all files in here are moved back to the main job folder.

Eventually, all jobs will hopefully be moved to ``completed/``, then you are done!


Example 1.1.1: High-throughput geometry optimisations with CASTEP with per-structure parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are a few occasions where you might need a custom ``.param`` file for each structure, for example, if using the implicit nanotube ``%devel_code`` in CASTEP.

These calculations are performed in exactly the same was as above, except a ``<structure_name>.param`` file must be made containing the required DFT parameters AND the nanotube parameters. In this case, run3 must now be called as::

    $ run3 --custom_params TiO2

.. tip:: If you have a .res file that contains a PyAIRSS "REM NT_PROPS" line, this will be ignored.


Example 1.2: High-throughput geometry optimisations with CASTEP on a supercomputer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each HPC facility has its own quirks, so in this example we will try to be as explicit as possible. The set up of the job is exactly the same as in `example 1 <ex1_>`__, but we now must add run3 to our job submission script. The following examples are for the SLURM queuing system on the BlueBear machine at the University of Birmingham and PBS on ARCHER (Tier-1), but run3 has also been tested on CSD3 (Tier-2), HPC Midlands-Plus (Tier-2), Thomas (Tier-2) and several local group-scale clusters.

Example 1.2.1: SLURM on BlueBear
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this job, we will submit a run3 job that performs CASTEP calculations across 2 24-core nodes per structure. Let us presume we have many thousand structures to run. The submission script looks as follows::

    $ cat run3.sub
    #!/bin/bash -l
    
    ###### MACHINE/USER-SPECIFIC OPTIONS ######
    
    #SBATCH --ntasks 48
    #SBATCH --nodes 2-2
    #SBATCH --time 24:00:00
    #SBATCH --qos <REDACTED> 
    ##SBATCH --qos bbshort
    #SBATCH --mail-type ALL
    #SBATCH --account=<REDACTED>

    module purge
    export PATH="$HOME/bin/CASTEP-17.21:$HOME/.conda/bin"
    module load bluebear
    module load mpi/impi/2017.1.132-iccifort-2017.1.132
    unset I_MPI_PMI_LIBRARY
    
    # RUN3 COMMANDS
    # (assuming installation guide followed at
    #  https://matador-db.readthedocs.io/en/latest/install.html)

    source activate matador
    run3 -nc 48 -v 4 --executable castep.mpi --ignore_jobs_file TiO2

Let's unpick a few of the flags used to call run3 here:

* ``-nc/--ncores``: the number of cores to use per structure, per calculation. It is often worth specifying this if more than one node is being used, as the correctness of run3's core counter is queue/machine-specific.
* ``-v 4``: sets the verbosity in the log file to the highest level.
* ``--ignore_jobs_file``: by default run3 will for both ``<structure>.lock`` files and entries in ``jobs.txt`` before running a new structure. It is often worth disabling the ``jobs.txt`` check if it is not expected that all structures complete in one job submission (see below).
  
try to call an executable called simply ``castep``. On many machines, CASTEP is installed as ``castep.mpi``.

Now to submit this script as a 200-node array job (i.e. running a maximum of 100 structures concurrently, depending on the queue), we call the following::

    $ sbatch --array=1-100 run3.job

It may be that this job is not large enough to optimise all structures within the walltime limit. In this case, it can be resubmitted using the same command. Jobs that were running when the walltime ran out should automatically be pushed back into the job folder so that they will be available to the next run3 call. In the event that this does not happen (for example MPI kept control of the Python thread for too long so the queuieng system interrupted run3's clean up), ``<hostname>`` folder
will be left hanging around in the main jobs folder. Jobs must then be manually made restartable by deleting ``<structure>.lock`` (and removing ``<structure>>`` from ``jobs.txt`` if not using ``--ignore_jobs_file``). It may also be that the intermediate CASTEP calculation was not copied over from the ``<hostname>`` folder: in this case, the CASTEP files can be updated by running::
    
    $ cp -u node*/*.castep .

from inside the root job folder.

Example 1.2.2: PBS on ARCHER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _ex.1.2:

Instructions are almost identical to the above, but the array job script looks a little different, for the same 100 copies of 2 node jobs (this time 24 cores per node)::

    $ cat run3.job
    #!/bin/bash --login
    # PBS job options (name, compute nodes, job time) # PBS -N is the job name (e.g. Example_MixedMode_Job)
    #PBS -N my_run3_job
    # PBS -l select is the number of nodes requested (e.g. 128 node=3072 cores)
    #PBS -l select=2
    # PBS -l walltime, maximum walltime allowed (e.g. 6 hours)
    #PBS -l walltime=24:00:00
    # Replace [budget code] below with your project code (e.g. t01)
    #PBS -A <REDACTED>
    #PBS -m abe
    #PBS -M <REDACTED>
    #PBS -J 1-100
    #PBS -r y

    # Make sure any symbolic links are resolved to absolute path
    export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

    # Change to the direcotry that the job was submitted from
    # (remember this should be on the /work filesystem)
    cd $PBS_O_WORKDIR

    source $HOME/.bashrc
    module load anaconda-compute/python3
    source activate $HOME/work/.conda/matador

    run3 --archer -v 4 -nc 48 KSnP

Notice here we have specified ``--archer``: again, run3 should be able to detect that ``mpirun`` is missing and thus try ``aprun``, but it can be worth specifying just in case. With PBS, the whole array can be submitted with just::

    $ qsub run3.job
.. index:: mongo

.. highlight:: bash

.. _mongo:

Setting up your own database
============================

``matador`` uses MongoDB to store and operate on results from first principles calculations. In this tutorial, we will set up a MongoDB database containing some dummy test data.

Installing and configuring MongoDB
----------------------------------

Instructions on how to install MongoDB community edition can be found at the `MongoDB website <https://docs.mongodb.com/manual/administration/install-community/>`_. If you are using Linux, you may find that your OS provides a MongoDB package that you can install.

Once installed, you can run your MongoDB server with the ``mongod`` command. An example configuration file can be found in ``matador/config/example_mongod.conf``, you may wish to change the port and default data directory (where the raw database is written to disk). This command will run persistently so you may wish to run it in a persistent terminal (using e.g. GNU screen or tmux), or as an OS-wide service (e.g. using systemctl or your OS's equivalent).

.. warning::

   It is up to the user to ensure that their MongoDB server is secure! This guide will assume access from ``localhost`` only, and thus will not set a password. If you do not wish external parties to have access to your database, make sure your server is running on secure/closed port, or set up authorisation as outlined on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/enable-authentication/>`_.

.. tip::

   You should ensure the path found under the ``storage.dbPath`` option in your configuration file exists before running ``mongod``, and that you have access to it!

Importing structures
--------------------

With a MongoDB server running, you can now import some structures. Navigate to ``examples/mongo/castep_files`` and run ``matador import``. If this is your first time using ``matador``, you should be prompted to make a configuration file. Follow the instructions and enter the values you used to set up your database. You will be given the option to test the connection to your database at the end.

If all is successful, you should see that 3 structures are added to your database, along with 4 failures. You can have a look at the spatula output files to see what happened. To finish this tutorial, use ``matador query`` to look at the structures you have just imported, before moving onto the command-line tutorial at :ref:`cli`.

Caring for your database
------------------------

Collections
^^^^^^^^^^^

MongoDB databases are made out of "collections" of data. To create a new collection, simply specify ``--db <collection_name>`` when performing an import. You may wish to create multiple collections of data, or one monolithic database.
   
The default collection can be set in your ``.matadorrc`` with the ``mongo.default_collection`` option. If you want to ensure that this database can only be imported to from a particular path, set the ``mongo.default_collection_file_path`` option.

A summary of the contents of existing collections can be generated using ``matador stats --db <collection_name>``. Similarly, a collection can be deleted using ``matador stats --db <collection_name> --delete``. If you wish to get a list of all existing collections, use ``matador stats --list``.

.. tip:: 
   In our group we combine these two approaches: a watched folder of "finalised" data is automatically scraped every day into a single monolithic collection (a cron job calls ``matador import``), interim calculations are then stored in collections per research project.

.. tip:: 
   Each collection also saves its own "changelog". You can view this changelog and undo changes with the ``matador changes`` interface.

Backing up
^^^^^^^^^^

You should always back up the data used to create a database, but if you wish, you can also directly backup the database itself using the MongoDB tool ``mongodump`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongodump/>`_). Similarly, restoration can be performed using ``mongorestore`` (`documentation <https://docs.mongodb.com/manual/reference/program/mongorestore/>`_). You may also wish to read the general back up and restore tutorial on the `MongoDB website <https://docs.mongodb.com/manual/tutorial/backup-and-restore-tools/>`_.
.. index:: run3_spectral

.. _run3_spectral:


Example 2: Spectral calculations with CASTEP and run3
-----------------------------------------------------

In this example, we will go from a crystal structure to a dispersion and DOS plot using run3, CASTEP and `OptaDOS <https://github.com/optados-developers/optados>`_. For this use case, run3 uses the `SeeK-path library <https://github.com/giovannipizzi/seekpath>`_ to generate standardised band paths through reciprocal space to automatically compute a useful bandstructure for all crystal types.

Spectral calculations follow a similar setup to geometry optimisations: run3 expects to find a folder containing .res files with lattice/atomic position data, and one .cell and one .param file specifying the CASTEP options. Standard run3 rules apply: if a ``<seed>.res.lock`` file is found, or if the .res file is listed in ``jobs.txt``, the structure will be skipped. Such a folder can be found in ``examples/bandstructure+dos/simple`` which contains some LiCoO\ :sub:`2` polymorphs. The Jupyter
notebook ``simple_spectral.ipynb`` will also show you exactly how to run a standard BS/DOS calculation and plot the results with the API. The files in ``examples/bandstructure+dos/projected`` will show you how to use matador and OptaDOS to get projected densities of states and bandstructures. Here, we shall run through some more simple cases first.

run3 will follow some simple rules to decide what kind of spectral calculation you want to run. First, it will check the ``task`` and ``spectral_task`` keywords. The ``task`` keyword needs to be set to ``'spectral'`` to trigger a spectral workflow. If ``spectral_task`` is either ``dos`` or ``bandstructure``, then this calculation type will always be included in the workflow, with default parameters if unset. Otherwise, run3 will check the cell file for the
``spectral_kpoints_mp_spacing`` (DOS) and ``spectral_kpoints_path_spacing`` (bandstructure) keywords and will perform either one or both of the corresponding calculation types. The first step will always be to perform an SCF calculation, which is continued from to obtain the DOS or BS (if any check file is found in the folder, run3 will attempt to skip the SCF and restart accordingly).

Both the bandstructure and DOS tasks can be post-processed with OptaDOS; run3 will perform this automatically if it finds a .odi file in alongside the .cell and .param. Similarly, if the ``write_orbitals_file`` keyword is set to True in the param file, ``orbitals2bands`` will be run automatically to reorder the eigenvalues in the .bands file. The output can be plotted using either the ``dispersion`` script bundled with matador, or through the API (see the Jupyter notebook).

Let us now consider some concrete examples.

Example 2.1: Bandstructure calculation with automatic path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can compute just bandstructure by specifying ``spectral_kpoints_path_spacing`` and ``spectral_task = bandstructure`` in the cell and param respectively, e.g. ::

   $ cat LiCoO2.param
   task: spectral
   spectral_task: bandstructure
   xc_functional: LDA
   cut_off_energy: 400 eV

   $ cat LiCoO2.cell
   kpoints_mp_spacing: 0.05
   spectral_kpoints_path_spacing: 0.05

with .res file ::

   $ cat LiCoO2-CollCode29225.res
   TITL LiCoO2-CollCode29225
   CELL 1.0  4.91 4.91 4.91 33.7 33.7 33.7
   LATT -1
   SFAC 	 Co Li O
   Co        1       0.000000        0.000000        0.000000         1.000000
   Li        2       0.500000        0.500000        0.500000         1.000000
   O         3       0.255463        0.255463        0.255463         1.000000
   O         3       0.744538        0.744538        0.744537         1.000000
   END

Simply calling ``run3 LiCoO2`` (with ``-v 4`` if you want to see what's going on) will perform three steps:

   1. Analyse the lattice with SeeKPath, perhaps creating a primitive cell or standardising if need be, and generate the standard Brillouin zone path with desired spacing set by ``spectral_kpoints_path_spacing``.
   2. Perform a self-consistent singlepoint energy calculation to obtain the electron density on the k-point grid specified by ``kpoints_mp_spacing``.
   3. Interpolate the desired bandstructure path, yielding a ``.bands`` file containing the bandstructure.

If you now run ``dispersion LiCoO2-CollCode29225 --png`` from the ``completed/`` folder, you should see this:

.. image:: LiCoO2-CollCode29225_spectral.png
   :name: bandstructure_only
   :align: center

Note that the path labels are generated from the .cell/.res file in the ``completed/`` folder, the .bands file does not contain enough information to make the entire plot portable. As mentioned above, if ``write_orbitals_file: true`` was found in the .param file, ``orbitals2bands`` would have been called to reorder the bands based on Kohn-Sham orbital overlaps.

.. tip::
   Plots can be customised using a custom matplotlib stylesheet. By default, matador plots will use the stylesheet found in ``matador/config/matador.mplstyle``, which can be copied elsewhere and specified by path in your ``.matadorrc`` under the ``plotting.default_style`` option. The default matplotlib styles can also be used directly by name, e.g. "dark_background".

.. tip::
   You can also specify a ``spectral_kpoint_path`` in your .cell file, that will be used for all structures in that
   folder. This can be useful when comparing the electronic structure of e.g. defected cells. The ``dispersion``
   script will try to scrape the labels to use for plotting from the comments at the specification of the kpoint path,
   e.g. ::

    %block spectral_kpoint_path
    0.0 0.0 0.0 ! \Gamma
    0.5 0.0 0.0 ! X
    %endblock spectral_kpoint_path

Example 2.2: A simple density of states (DOS) calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We shall now calculate the density of states in a similar way to the above, using the matador default for ``spectral_kpoints_mp_spacing``, such that the cell file only contains ``kpoints_mp_spacing: 0.05`` and the param file now has ``spectral_task: DOS``.

.. tip::
   If you are starting from the example above, you will need to move the .res file back into the current folder, delete the .txt files and rename the ``completed/`` folder to e.g. ``completed_bs/``. You might also consider moving the .check file to current folder, so that the calculation will re-use the SCF results.

Again, simply running ``run3 LiCoO2`` will do the trick. Eventually, a .bands_dos file will be produced (any existing .bands files will be backed up and replaced at the end of the run). The dispersion script will recognise this as a density of states, and will apply some naive gaussian smearing that can be controlled with the ``-gw/--gaussian_width`` flag. Running ``dispersion LiCoO2-CollCode29225 --png -gw 0.01`` will produce the following:

.. image:: LiCoO2-CollCode29225_spectral_dos.png
   :name: dos_only
   :align: center

.. tip::
   If you have a .bands file remaining in your top directory, ``dispersion`` will try to plot this as a bandstructure alongside your DOS, which may look terrible if .bands contains a DOS calculation! You can plot just the DOS using the ``--dos_only`` flag.


Example 2.3: Putting it all together
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a DOS and bandstructure on the same structure, simply include both ``spectral_kpoints_mp_spacing`` and ``spectral_kpoints_path_spacing`` in your .cell file. Your ``spectral_task`` keyword in the param file will be ignored. This exact example can be found in ``examples/bandstructure+dos/simple``, with an example Jupyter notebook showing how to make plots with the API directly, rather than the dispersion script.

After calling run3 again, the ``completed/`` folder in this case should contain both a .bands and a .bands_dos file which can be plotted alongside one another using ``dispersion LiCoO2-CollCode29225``, to produce the following:

.. image:: LiCoO2-CollCode29225_spectral_both.png
   :name: dos_bs
   :align: center

Example 2.4: Using OptaDOS for post-processing: projected DOS and bandstructures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The final piece of the puzzle is `OptaDOS <https://github.com/optados-developers/optados>`_, a package for broadening and projecting densities of states (amongst other things) that comes with CASTEP. By default, run3 will turn on the required CASTEP settings (namely ``pdos_calculate_weights``) required by OptaDOS. In order for OptaDOS to be run automatically by run3, an extra .odi file must be added into our input deck, containing the details of the desired OptaDOS calculation.

.. note::
   This example assumes that the OptaDOS binary is called ``optados`` and resides in your PATH, likewise ``orbitals2bands``. This can altered by setting the ``run3.optados_executable`` setting in your matador config.

.. warning:: By default, OptaDOS will *not* be invoked with ``mpirun`` (i.e., your executable should work for serial runs too). A parallel OptaDOS run can be performed by setting the ``run3.optados_executable`` to e.g. ``mpirun optados.mpi`.

.. warning::
   The projected dispersion curve feature is quite new to OptaDOS and thus is temperamental. Depending on when you are reading this, it may require you to have compiled OptaDOS from the development branch on GitHub.

run3 will try to perform three types of calculation: a simple DOS smearing, a projected density of states (with projectors specified by the ``pdos`` keyword), and a projected bandstructure (with projectors specified by the ``pdispersion`` keyword). If ``pdos``/``pdispersion`` is not found in the .odi, this corresponding task will be skipped. Likewise, if ``broadening`` is not found in the .odi, the standard DOS broadening will not be performed.::

   $ cat LiCoO2.odi
   pdos: species_ang
   pdispersion: species
   adaptive_smearing: 1
   set_efermi_zero: True
   dos_per_volume: True
   broadening: adaptive
   dos_spacing: 0.01

With all these files in place, simply running ``run3 LiCoO2`` and ``dispersion (-interp 3 -scale 25) LiCoO2-CollCode29225`` (optional flags in brackets) should yield the following plot:

.. image:: LiCoO2-CollCode29225_spectral_pdis.png
   :name: full_spectral
   :align: center

.. tip:: Note that the colours of each projectors in these plots is set by your VESTA colour scheme, which is bundled by default inside ``matador/config``.
.. index:: run3_qe

.. highlight:: bash

.. _run3_elastic:

Example 5: Geometry optimisations with Quantum Espresso and run3
================================================================

This tutorial uses the files found in ``examples/run3_quantum_espresso/vc-relax`` to relax
a folder of res files in a similar way to the CASTEP tutorial. Quantum Espresso (QE) uses a single
input file (as opposed to CASTEP's two), and these must be prepared beforehand (rather than
generated per structure with run3).

First things first, set up your directory so that you have:

* a load of res files
* a template QE input file containing DFT parameters

First, we make the QE input files using the ``shx3pwscf`` script:
``shx3pwscf *.res --template vcr.template --kpoint_spacing 0.03``

Then, to run the optimisations with run3 (either interactively, or at the bottom of a job script), run3 must be called on all the ``*.in`` files:
``run3 --redirect "$seed.out" -nc 4 --mode generic --executable 'pw.x -i $seed.in' *.in``
The ``$seed`` variables will be expanded by run3 to take the values of the ``.res`` file names. In this case, our only res file is called NaP.res, so the only calculation will be called as ``mpirun -n 4 pw.x -i NaP.in > NaP.out`` (or equivalent).
.. index:: run3_elastic

.. highlight:: bash

.. _run3_elastic:


Example 4: Bulk moduli with CASTEP and run3
-------------------------------------------

In this tutorial we will calculate the bulk moduli of diamond, silicon and lithium, using CASTEP and run3. 

The general process for calculating the bulk modulus is as follows:

1. Relax a structure to its equilibrium volume :math:`V_0` (ideally relaxing the positions
   too).
2. Run total energy calculations at a series of perturbed volumes :math:`\alpha V_0`
   for :math:`\alpha \in [0.95, 1.05]`.
3. Fit a particular form of equation of state to the resulting :math:`E(V)` curve and extract the bulk modulus.

Note: there is a supported set of scripts for calculating all the elastic constants in CASTEP, which
can be found on `GitHub <https://github.com/andreww/elastic-constants>`_ with an additional tutorial on the `CASTEP website <http://www.castep.org/Tutorials/ElasticConstants>`_.

As with the other tutorials, run3 expects to find a series of structures as ``.res`` files and single ``$seed.cell`` and ``$seed.param`` files. The files for this example can be found in ``examples/run3_elastic``. The ``.cell`` file in this tutorial is basically identical to that in the geometry optimisation tutorial, but the ``.param`` file this time contains the line ``task: bulk_modulus``. This is *not* a valid CASTEP task (as of 2019), but instead run3 will capture this task and use it to spawn a series of single point jobs.::

    $ cat bulk_mod.cell
    kpoints_mp_spacing: 0.05
    snap_to_symmetry
    symmetry_tol: 0.01
    symmetry_generate
    %block species_pot
    QC5
    %endblock species_pot

    $ cat bulk_mod.param
    task: bulk_modulus
    cut_off_energy: 300 eV
    write_bib: false
    write_checkpoint: none
    geom_max_iter: 100
    xc_functional: LDA


Once the files are in place, calling run3 as ``run3 bulk_mod`` will begin to relax
the first structure before deforming it. You can track the calculations
progress in the log files found in ``logs/``.

After the calculations have completed, your completed folder should look
something like this::

   $ ls completed
    C-OQMD_675640-CollCode28857.bands             Li-OQMD_30676-CollCode642106_bulk_mod.castep
    C-OQMD_675640-CollCode28857.castep            Li-OQMD_30676-CollCode642106_bulk_mod.cell
    C-OQMD_675640-CollCode28857.cell              Li-OQMD_30676-CollCode642106_bulk_mod.cst_esp
    C-OQMD_675640-CollCode28857.cst_esp           Li-OQMD_30676-CollCode642106_bulk_mod.param
    C-OQMD_675640-CollCode28857.geom              Li-OQMD_30676-CollCode642106_bulk_mod.png
    C-OQMD_675640-CollCode28857.param             Li-OQMD_30676-CollCode642106_bulk_mod.res
    C-OQMD_675640-CollCode28857.res               Li-OQMD_30676-CollCode642106_bulk_mod.results
    C-OQMD_675640-CollCode28857_bulk_mod.bands    Si-OQMD_5714-CollCode29287.bands
    C-OQMD_675640-CollCode28857_bulk_mod.castep   Si-OQMD_5714-CollCode29287.castep
    C-OQMD_675640-CollCode28857_bulk_mod.cell     Si-OQMD_5714-CollCode29287.cell
    C-OQMD_675640-CollCode28857_bulk_mod.cst_esp  Si-OQMD_5714-CollCode29287.cst_esp
    C-OQMD_675640-CollCode28857_bulk_mod.param    Si-OQMD_5714-CollCode29287.geom
    C-OQMD_675640-CollCode28857_bulk_mod.png      Si-OQMD_5714-CollCode29287.param
    C-OQMD_675640-CollCode28857_bulk_mod.res      Si-OQMD_5714-CollCode29287.res
    C-OQMD_675640-CollCode28857_bulk_mod.results  Si-OQMD_5714-CollCode29287_bulk_mod.bands
    Li-OQMD_30676-CollCode642106.bands            Si-OQMD_5714-CollCode29287_bulk_mod.castep
    Li-OQMD_30676-CollCode642106.castep           Si-OQMD_5714-CollCode29287_bulk_mod.cell
    Li-OQMD_30676-CollCode642106.cell             Si-OQMD_5714-CollCode29287_bulk_mod.cst_esp
    Li-OQMD_30676-CollCode642106.cst_esp          Si-OQMD_5714-CollCode29287_bulk_mod.param
    Li-OQMD_30676-CollCode642106.geom             Si-OQMD_5714-CollCode29287_bulk_mod.png
    Li-OQMD_30676-CollCode642106.param            Si-OQMD_5714-CollCode29287_bulk_mod.res
    Li-OQMD_30676-CollCode642106.res              Si-OQMD_5714-CollCode29287_bulk_mod.results
    Li-OQMD_30676-CollCode642106_bulk_mod.bands

You can see the results of the bulk modulus fits for each structure in
``completed/*_bulk_mod.results``, along with plots of the fits saved as pngs.

Compare these LDA values from the three different fits to the experimental values from Wikipedia:

+----------+--------------------------+------------------------------+
| Material | Expt. Bulk modulus (GPa) | Predicted Bulk Modulus (GPa) |
+==========+==========================+==============================+
| Lithium  | 11                       | 17.6, 17.7, 17.8             |
+----------+--------------------------+------------------------------+
| Silicon  | 98                       | 97.5, 98.0, 97.4             |
+----------+--------------------------+------------------------------+
| Diamond  | 442                      | 503, 496, 496                |
+----------+--------------------------+------------------------------+
.. index:: run3

.. highlight:: bash

.. _run3:

Getting started with run3
=========================

run3 is the main entrypoint for performing high-throughput calculations with the matador library. Broadly, the aim of run3 is to perform the same calculation, with the same parameters, on multiple structures concurrently.

.. toctree::
    :maxdepth: 0

    run3_geom
    run3_spectral
    run3_phonon
    run3_elastic
    run3_qe
