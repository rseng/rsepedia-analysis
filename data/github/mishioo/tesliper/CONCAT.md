[![Build status](https://ci.appveyor.com/api/projects/status/vh0t6udj7mnpnfoe?svg=true)](https://ci.appveyor.com/project/mishioo/tesliper-jjshl)
[![Documentation Status](https://readthedocs.org/projects/tesliper/badge/?version=stable)](https://tesliper.readthedocs.io/en/stable/?badge=stable)
[![Coverage Status](https://coveralls.io/repos/github/mishioo/tesliper/badge.svg)](https://coveralls.io/github/mishioo/tesliper)
[![PyPi version](https://badgen.net/pypi/v/tesliper/)](https://pypi.org/project/tesliper)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/tesliper.svg)](https://pypi.python.org/pypi/tesliper/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![License](https://img.shields.io/badge/License-BSD_2--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04164/status.svg)](https://doi.org/10.21105/joss.04164)

<img align="right" width="100" height="100" src="https://raw.githubusercontent.com/mishioo/tesliper/master/tesliper/tesliper.ico">

# tesliper

`tesliper`: Theoretical Spectroscopist's Little Helper is a program for batch processing
of Gaussian output files, focused on calculations of vibrational, electronic, and
scattering spectra from Gaussian-calculated quantum properties of molecule's conformers.

It allows to easily exclude conformers that are not suitable for further analysis:
erroneous, not optimized, of higher energy or lower contribution than a user-given
threshold. It also implements an RMSD sieve, enabling one to filter out similar
structures. It lets you calculate theoretical IR, VCD, UV, ECD, Raman, and ROA spectra
for each conformer or as a population-weighted average and export obtained spectral data
in one of supported file formats: `.txt`, `.csv`, or `.xlsx`. Finally, if allows to
easily setup a next calculations step with batch export to `.gjf` Gaussian input files.

`tesliper` is written in Python 3.6 and makes use of some additional third party
packages (see below or requirements.txt). It may be used as a package or as a
stand-alone application with dedicated GUI.

- [tesliper](#tesliper)
- [Getting Started](#getting-started)
  - [Requirements](#requirements)
  - [Installing to your Python distribution](#installing-to-your-python-distribution)
  - [A standalone application](#a-standalone-application)
- [How to use](#how-to-use)
  - [Primer](#primer)
  - [In Python scripts](#in-python-scripts)
  - [A graphical interface](#a-graphical-interface)
- [License](#license)
- [Contributing to `tesliper`](#contributing-to-tesliper)
  - [Bugs and suggestions](#bugs-and-suggestions)
  - [Participating in code](#participating-in-code)
  - [Roadmap](#roadmap)
  - [Acknowledgements](#acknowledgements)
- [How to cite](#how-to-cite)

# Getting Started

You can use `tesliper` from python or as standalone application with dedicated graphical
user interface. See below for details.

## Requirements

```
Python 3.6+
numpy
openpyxl
matplotlib (optional, for GUI)
```

`tesliper` uses `tkinter` to deliver the graphical interface. It is included in
most Python distributions, but please be aware, that some might miss it. You will
need to install it manually in such case.

## Installing to your Python distribution

`tesliper` is available on PyPI, you can install it to your python distribution by
simply running:

`$ python -m pip install tesliper`

or

`$ python -m pip install tesliper[gui]`

if you would like to be able to use a graphical interface.

## A standalone application

This option is currently available only for Windows users.
To get your copy of `tesliper`, simply download a `Tesliper.exe` file from the
[latest relase](https://github.com/Mishioo/tesliper/releases/latest).
This file is a standalone application, no installation is required.

# How to use
Full documentation is available online: https://tesliper.readthedocs.io/.

## Primer
Conventions that are important to note:

- `tesliper` stores multiple data entries of various types for each conformer. To
prevent confusion with Python's data ``type`` and with data itself, `tesliper` refers to
specific kinds of data as "genres". Genres in code are represented by specific strings,
used as identifiers. To learn about data genres known to `tesliper`, see [Available data
genres](https://tesliper.readthedocs.io/en/stable/genres.html).
- `tesliper` identifies conformers using a stem of an extracted file (i.e. its filename
without extension). When files with identical names are extracted in course of
subsequent `Tesliper.extract` calls or in recursive extraction using
``Tesliper.extract(recursive=True)``, they are treated as data for one conformer.
This enables to join data from subsequent calculations steps, e.g. geometry
optimization, vibrational spectra simulation, and electronic spectra simulation.
Please note that if specific data genre is available from more than one calculation
job, only recently extracted values will be stored.
- `tesliper` was designed to deal with multiple conformers of single molecule and may
not work properly when used to process data concerning different molecules (i.e.
having different number of atoms, different number of degrees of freedom, etc.).
If you want to use it for such purpose anyway, you may set
`Tesliper.conformers.allow_data_inconsistency` to ``True``. `tesliper` will then stop
complaining and try to do its best.

## In Python scripts
`Tesliper` class is the main access point to `tesliper`'s functionality. It allows you
to extract data from specified files, provides a proxy to the trimming
functionality, gives access to data in form of specialized arrays, enables you
to calculate and average desired spectra, and provides an easy way to export data.

Most basic use might look like this:
```python
from tesliper import Tesliper
tslr = Tesliper()
tslr.extract()
tslr.calculate_spectra()
tslr.average_spectra()
tslr.export_averaged()
```
This extracts data from files in the current working directory, calculates
available spectra using standard parameters, averages them using available energy
values, and exports to current working directory in .txt format.

You can customize this process by specifying call parameters for used methods
and modifying `Tesliper`'s configuration attributes:
- to change source directory
or location of exported files instantiate `Tesliper` object with `input_dir`
and `output_dir` parameters specified, respectively. You can also set appropriate
attributes on the instance directly;
- to extract only selected files in `input_dir` use `wanted_files` init parameter.
It should be given an iterable of filenames you want to parse. Again, you can
also directly set an identically named attribute;
- use `Tesliper.conformers.trim...` methods to easily filter out conformers you wish
to ignore in further analysis;
- to change parameters used for calculation of spectra, modify appropriate entries
of `parameters` attribute;
- use other export methods to export more data and specify `fmt` parameter
in method's call to export to other file formats.

```python
tslr = Tesliper(input_dir="./myjob/optimization/", output_dir="./myjob/output/")
tslr.wanted_files = ["one", "two", "three"]  # only files with these names
tslr.extract()  # use tslr.input_dir as source
tslr.extract(path="./myjob/vcd_sim/")  # use other input_dir
tslr.conformers.trim_not_optimized()  # filtering out unwanted conformers
tslr.parameters["vcd"].update({"start": 500, "stop": 2500, "width": 2})
tslr.calculate_spectra(genres=["vcd"])  # we want only VCD spectrum
tslr.average_spectra()
tslr.export_averaged(mode="w")  # overwrite previously exported files
tslr.export_activities(fmt="csv")  # save activities for analysis elsewhere
tslr.output_dir = "./myjob/ecd_sim/"
tslr.export_job_file(  # prepare files for next step of calculations
    route="# td=(singlets,nstates=80) B3LYP/Def2TZVP"
)
```

## A graphical interface
If you are using `tesliper` as a standalone application, simply double click on the
`Tesliper.exe` file to start the application. To invoke it from the command line,
just run `tesliper-gui`. GUI consists of three panels and a number of controls.
The panels are: "Extracted data", "Energies list", and "Spectra view". First two
offer a list of conformers read so far using "Chose files" and "Chose folder" buttons
on the left. The last enables to preview calculated spectra.

![screenshot](https://raw.githubusercontent.com/mishioo/tesliper/master/docs/source/_static/screenshots/16426974012.png)

- "Extracted data" panel shows an identifier of each conformer (a file name) and an
overview of data extracted. Little checkboxes on the left of each conformer may be
clicked to mark this conformer as "kept" or "not kept".
- "Energies list" offers the same set of conformers and checkboxes, but with energies
values listed for each conformer. The view may be changed using "Show" dropdown box
in "Energies and Structure" section of controls, to present difference in energy
between conformers or their percentage contribution in population.
- "Spectra view" tab shows calculated spectra. It may be controlled using "Calculate
spectra" section. After choosing a spectra genre to calculate you may control if it is
simulated using lorentzian or gaussian fitting function, change peak width, spectra
bounds, etc. You may view spectra for one conformer, all of them stacked, or averaged.
You may also load an experimental spectrum (.txt or .csv format) for comparison.

![screenshot](https://raw.githubusercontent.com/mishioo/tesliper/master/docs/source/_static/screenshots/16426984646.png)

Once done with extracting files and tweaking parameters, export selected data to desired
format or save the session for later using buttons in "Session control" section.

A detailed tutorial with screenshots is available in the documentation:
https://tesliper.readthedocs.io/en/stable/gui.html.

# License
This project is licensed with BSD 2-Clause license.
See [LICENSE.txt](https://github.com/mishioo/tesliper/blob/master/LICENSE.txt)
file for details.

# Contributing to `tesliper`
Contributions are welcome! `tesliper` is a growing project and definitely has room
for improvements. 

## Bugs and suggestions
Bug reports are of great value, if you encounter a problem please let me know by
submitting a [new issue](https://github.com/mishioo/tesliper/issues/new).
If you have a suggestion how `tesliper` can be improved, please let me know as well!

## Participating in code
If you'd like to contribute to `tesliper`'s codebase, that's even better! If there is a
specific bug that you know how to fix or a feature you'd like to develop, please let me
know via [issues](https://github.com/mishioo/tesliper/issues). To start coding, get your
working copy of `tesliper`'s source code by cloning the repository and then setup your
environment by installing development dependencies (probably to a virtual environment):
```
python -m pip install .[dev]
```
Please remember to add/update relevant tests along with your code changes.
Make sure the test suite passes by running
```
python -m pytest test
```
or, if you'd like to also generate the coverage report, you can invoke the pytest
with the pytest-cov plugin (included in the "dev" installation), as shown below.
The .coveragerc file includes plugin's configuration (i.e. to omit GUI code).
```
python -m pytest --cov=tesliper --cov-config=.coveragerc test
```
If you contribute code to `tesliper`, the test coverage badge will be updated by the CI
service after your changes are merged.

Although [`mypy`](https://mypy.readthedocs.io/) is not incorporated in `tesliper`'s
development (yet!), I believe that type hints greatly improve the experience of using
the package. Please, add them to the new code you submit. If you're willing to
supplement typing of the existing code, that would be very much welcome!

`tesliper`'s codebase is formatted with [`black`](https://black.readthedocs.io/) and
[`isort`](https://pycqa.github.io/isort/), please use these tools before submitting your
code contribution. To make it easier, you may use a pre-commit configuration available
in this repository. To include it in your workflow, simply run `pre-commit install` in
your copy's root directory. Note, that this configuration also sets up
[`flake8`](https://flake8.pycqa.org/) linter.

To get your change introduced to the codebase, please make a Pull Request to the `fixes`
branch for quick bug fixes or to the `dev` branch for new features and bigger changes.
If at a loss, do not hesitate to reach to me directly! :)

## Roadmap
Ideas for possible future improvements to the software are listed below. Based on the
feedback from the Community, I will decide, which ones are desired and worth working on.

### New Functionality: <!-- omit in toc -->
- command line interface
- support for Jaguar & other packages
- option for using `cclib` for parsing
- spectra comparison feature
- velocity and length electronic spectra comparison
- parsing conformational search files
- GUI: manually choose spectra colour

### UX Improvements: <!-- omit in toc -->
- supplement logging
- auto finding optimal spectra range
- GUI: export spectra as image
- GUI: spectra colour by population
- GUI: drag&drop support
- GUI: display tooltips on hover


## Acknowledgements

Many thanks to the scientists, who advised me on the domain-specific details and helped
to test the software:

[orcid_logo]: https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png

- Joanna Rode [![orcid_id][orcid_logo] 0000-0003-0592-4053](https://orcid.org/0000-0003-0592-4053)
- Magdalena Jawiczuk [![orcid_id][orcid_logo] 0000-0003-2576-4042](https://orcid.org/0000-0003-2576-4042)
- Marcin Górecki [![orcid_id][orcid_logo] 0000-0001-7472-3875](https://orcid.org/0000-0001-7472-3875)

as well as to people, who reviewed the project: [@alejandrogallo](https://github.com/alejandrogallo) and [@arepstein](https://github.com/arepstein).

# How to cite

I'm very happy to announce, that `tesliper` has been published in the [Journal of Open
Source Software](https://joss.theoj.org/papers/10.21105/joss.04164)! If you find
`tesliper` useful in your research, please give it a credit by citing it:

- Więcław, M. M., (2022). tesliper: a theoretical spectroscopist's little helper.
  *Journal of Open Source Software*, 7(72), 4164, DOI:
  [10.21105/joss.04164](https://doi.org/10.21105/joss.04164)

```
@article{tesliper2022,
  doi = {10.21105/joss.04164},
  url = {https://doi.org/10.21105/joss.04164},
  year = {2022},
  publisher = {The Open Journal},
  volume = {7},
  number = {72},
  pages = {4164},
  author = {Michał M. Więcław},
  title = {tesliper: a theoretical spectroscopist's little helper},
  journal = {Journal of Open Source Software}
}
```


---
title: "tesliper: a theoretical spectroscopist's little helper"
tags:
  - Python
  - chemistry
  - spectroscopy
  - Gaussian
  - chemical computing 
  - optical spectroscopy
  - spectral simulations
  - workflow automation
  - batch processing 
authors:
  - name: Michał M. Więcław^[corresponding author]
    orcid: 0000-0001-7884-8982
    affiliation: 1
affiliations:
 - name: Institute of Organic Chemistry, Polish Academy of Sciences
   index: 1
date: 26 January 2022
bibliography: paper.bib
---

# Summary
`tesliper` is a software package for the bulk processing of Gaussian (quantum chemistry
calculations software) output files regarding conformational searches and spectra
simulations. Simulation of a molecule's optical spectra usually requires structure
optimization and calculation of the electronic properties of multiple conformers of the
studied molecule. Gaussian [@gaussian] is one of the most commonly used software
packages to perform such calculations and `tesliper` was created to aid the handling and
analysis of these calculations' outputs.

It allows for easy exclusion of conformers that are not suitable for further analysis:
erroneous, not optimized, of higher energy or lower contribution than a desired
threshold. It also implements a geometry comparison feature: an RMSD sieve, enabling the
filtering out of similar structures. A theoretical IR, VCD, UV, ECD, Raman, and ROA
spectra may be calculated for individual conformers or as a population-weighted average.
Offering a graphical user interface and Python API, it is easily accessible by the
users' preferred method of interaction with the computer: visual or textual.

# Statement of Need
Simulation of optical spectra of organic compounds has become a routine task for
chemical analysts. For example, it is a necessary step in one of the increasingly
popular methods of establishing a compound's absolute configuration through comparison
of recorded and simulated circular dichroism spectra. However, the process of obtaining
a simulated spectrum may be cumbersome, as it usually involves analyzing a large number
of potentially stable conformers of the studied molecule.

Several software packages capable of simulating a spectrum from Gaussian calculations
are already available, e.g. SpecDis [@specdis], CDspecTech [@cdspectech], ComputeVOA
[@computevoa], ChemCraft [@chemcraft], or GaussView [@gaussview]. However, each of these
programs has some severe limitations, for example they are not freely available
(ComputeVOA, ChemCraft, GaussView) or cannot simulate certain types of spectra (SpecDis,
ComputeVOA, ChemCraft). Other packages  lack important features, like comparison of
conformers' geometry (SpecDis, CDspecTech, GaussView) or population-based spectra
averaging (ComputeVOA, ChemCraft, GaussView).

Even with adoption of one or more of the software packages mentioned above, the process is often
suboptimal, incomplete, or unable to be done in an automated fashion. Many research groups
tackle this problem with home-brewed scripting solutions that are usually not easily
adjusted or extended, or with manual work, which may be extremely time-consuming. There
is a clear need for a simple interface to automate tedious parts of the typical spectra
simulation workflow. `tesliper` aims to satisfy this need.

# Acknowledgements
Many thanks to the scientists, who advised me on the domain-specific details and helped
to test the software:

- Joanna Rode [ORCID iD 0000-0003-0592-4053](https://orcid.org/0000-0003-0592-4053)
- Magdalena Jawiczuk [ORCID iD 0000-0003-2576-4042](https://orcid.org/0000-0003-2576-4042)
- Marcin Górecki [ORCID iD 0000-0001-7472-3875](https://orcid.org/0000-0001-7472-3875)

# References

Change Log
==========

v. 0.9.3
--------

GUI:
    - Added button for recursive extraction.

Other Changes:
    - Now warning will be issued after reading abnormally terminated files.
    - Minor corrections in the documentation.

v. 0.9.2
--------

Bug Fixes:
    - Fixed ``__version__`` and other metadata attributes broken in 0.9.1.

v. 0.9.1
--------

Bug Fixes:
    - Fixed ``ImportError`` occurring in Python 3.10.
    - Corrected creation of ``"filanemes"`` pseudo-genre.
    - Corrected ``len()`` behavior with ``Spectra`` instances.

New Features:
    - Added "top-level" temperature setting in both, API and GUI.
    - Allowed ignoring of unexpected keyword arguments in ``Conformers.arrayed()``.

Other Changes:
    - Moved requirements to setup.py file.
    - Added ``tesliper-gui`` entry point.
    - ``Tesliper.get_averaged_spectrum()`` now tries to calculate missing spectra.
    - Minor supplementation to documentation and READEME.

v. 0.9.0
--------

Created online documentation! Available at https://tesliper.readthedocs.io/

Bug Fixes:
    - Fixed error on parsing radical molecules.
    - Corrected ``ArrayProperty`` ignoring it's ``.fill_value``.
    - Fixed infinite recursion error on ``SpectralData.wavelen`` access.
    - Prevented creation of empty files on export of empty data arrays.
    - Prevented intermediate ``.xlsx`` file saving when exporting multiple data genres.
    - Corrected trimming abnormally terminated conformers in GUI.

New Features:
    - ``rmsd_sieve`` and ``Conformers.trim_rmsd`` now allow for arbitrary windows.
    - Added ``datawork.geometry.pyramid_windows`` window strategy function.
    - Extended ``Soxhlet`` to allow use of arbitrary registered parsers.
    - Allowed for automatic instantiation of data arrays for genres that depend on a different genre.
    - Introduced *optimized_geom* genre
    - Added export of generic data arrays.
    - Added parametrization of ``GjfWriter.link0`` commands.

Other Changes:
    - Reviewed and corrected calculation of intensities.
    - Improved automatic scaling of spectra.
    - Renamed ``Parser`` to ``ParserBase`` for consistency with other base classes.
    - Unified base classes' registering mechanism of their subclasses.
    - Cleaned up ``extraction.gaussian_parser``. Changed all data sequences to lists. 
    - Supplemented type hints.
    - Renamed *geometry* genre to *last_read_geom*.
    - Supplemented ``Conformers`` to fully implement ``OrderedDict`` interface.
    - Added storage and serialization of experimental spectra.

GUI:
    - Unified terminology used with the one in code and documentation.

v. 0.8.2
--------

API:
    - Corrected data export when ``Tesliper``'s default genres used.
    - Corrected error when ``Tesliper.calculate_spectra`` called with default values.
    - Corrected default filenames generated for spectral data and activities.
    - Supplemented genres' full names and other metadata.

v. 0.8.1
--------

API:
    - Corrected handling of invalid start, stop, step parameters combination when calculating spectra.
GUI:
    - Fixed incorrect floats' rounding in numeric entries.
    - Added reaction (trim conformers/redraw spectra) to "Enter" key press, when editing a numeric entry.
    - Fixed an error occurring when "show activities" is checked but there are no activities in a plotting range.
    - Added auto-update of energies-related values after trimming.


v. 0.8.0
--------

API:
    - added RMSD-based trimming of conformers with similar geometry
    - added auto scaling and shifting spectra to match reference
    - added support for handling and exporting electronic transitions
    - added export to .gjf files
    - added serialization of ``Tesliper`` class
    - renamed ``Molecules`` class to ``Conformers``
    - significant changes to ``...Writer`` classes
    - significant changes to ``DataArray`` subclasses
    - major code refactoring
    - many smaller changes and improvements
GUI:
    - new application layout
    - added scroll response to numeric fields
    - changed available and default colour schemes
    - supplemented data export options


v. 0.7.4
--------

API:
    - Tesliper's method 'average_spectra' returns reference to dict of averaged spectra
GUI:
    - fixed files export (broken in v. 0.7.3)


v. 0.7.3
--------

API:
    - introduced exceptions.py submodule
    - glassware module turned into package
    - improved mechanism for dealing with inconsistent data sizes
    - added mechanism for trimming conformers with inconsistent data sizes
    - fixed Molecules' trim_incomplete function
    - enhanced Molecules' trim_non_matching_stoichiometry function
    - introduced dict_view classes for iteration through trimmed Molecules 
    - improved Molecules indexing mechanism to return in O(1)
    - removed 'cpu_time' from data extracted by gaussian_parser
    - fixed error on parsing ECD calculations from g.09B 
GUI:
    - fixed problem with stacked spectra drawing 
    - added spectra reversing on demand
    - fixed stacked spectra coloring
    - corrected bars drawing for uv and ecd spectra
    - added option for filtering conformers with inconsistent data sizes
    - split un/check into separate buttons
    - fixed checking/unchecking incomplete entries
    - added checking/unchecking inconsistent sizes
    - other minor changes and fixes


v. 0.7.2
--------

- added support for string 'genres' parameter in Tesliper.calculate_spectra method
- added support for .xy spectra files
- gui: fixed problem with averaged and stacked spectra drawing 
- gui: set "user_home_dir/tesliper/" as default location for tslr_err_log.exe
- other minor fixes and enhancements


v. 0.7.1
--------

- fixed crash on spectra drawing when Matplotlib 3 used
- fixed problem with loading spectra from some txt files
- added support for loading spectra from csv files
- other minor fixes


v. 0.7.0
--------

- graphical user interface redesigned
- significant changes in code architecture
- many fixes


v. 0.6.4
--------

- calculated spectra precision in txt files changed to e-4
- spectra lines width changed
- data trimming features corrected
- spectra plot erasing on session clearing implemented
- inverting x axis for uv and ecd spectra added


v. 0.6.3
--------

- fixed export error when not chosen, but all data were exported
- fixed export error when export occurred after closing popup window
- fixed export error when energies were not exported to separate txt files
- entry validation improved


v. 0.6.2
--------

- solved some problems with corrupted files extraction
- added warning when files from mixed gaussian runs found
- fixed RuntimeError on overlapping actions
- fixed export popup error
- errors description moved to tslr_err_log.txt
- fixed ValueError on empty settings in gui_main.current_settings
- corrected session instantiation from files (unwanted files problem)
- changed energies precision to .6
- added Min. Boltzmann factor in GUI


v. 0.6.1
--------

First beta release


v. 0.6.0 and earlier
--------------------

Early development stagesAdvanced guide
==============

``tesliper`` handles data extracted from the source files in a form of specialized
objects, called :term:`data array`\s. These objects are instances of one of the
:class:`.DataArray` subclasses (hence sometimes referenced as ``DataArray``-like
objects), described here in a greater detail. :class:`.DataArray` base class defines a
basic interface and implements data validation, while its subclasses provided by
``tesliper`` define how certain data :term:`genre`\s should be treated and processed.

.. note::

    Under the hood, :term:`data array`\s, and ``tesliper`` in general, use `numpy
    <https://numpy.org>`_ to provide fast numeric operations on data.

This part of documentation also shows how to take more control over the data export.
:class:`.Tesliper` autotomizes this process quite a bit and exposes only a limited set
of possibilities provided by the underlying writer classes. It will be shown here how to
use these writer classes directly in your code.

Data array classes
------------------

Each ``DataArray``-like object has the following four attributes:

genre
    name of the data genre that *values* represent;
filenames
    sequence of conformers' identifiers as a ``numpy.ndarray(dtype=str)``;
values
    sequence of values of *genre* data genre for each conformer in *filenames*. It is
    also a :class:`numpy.ndarray`, but its ``dtype`` depends on the particular data
    array class;
allow_data_inconsistency
    a flag that controls the process of data validation. More about data inconsistency
    will be said later.

Some data arrays may provide more data. For example, any spectral data values wouldn't
be complete without the information about the band that they corresponds to, so data
arrays that handle this kind of data also provide a *frequencies* or *wavelengths*
attribute.

.. note::

    Attributes that hold a band information are actually *freq* and *wavelen*
    respectively, *frequencies* and *wavelengths* are convenience aliases.


Creating data arrays
''''''''''''''''''''

The easiest way to instantiate the data array of desired data genre is to use
:meth:`.Conformers.arrayed` factory method. It transforms it's stored data into the
``DataArray``-like object associated with a particular data genre, ignoring any
conformer that is :term:`not kept <kept>` or doesn't provide data for the requested
genre. You may force it to ignore any trimming applied by adding ``full=True`` to call
parameters (conformers without data for requested genre still will be ignored).
Moreover, any other keyword parameters provided will be forwarded to the class
constructor, allowing you to override any default values.

.. code-block:: python

    >>> from tesliper import Conformers
    >>> c = Conformers(
    ...     one={"gib":-123.5},
    ...     two={},
    ...     three={"gib": -123.6},
    ...     four={"gib":-123.7}
    ... )
    >>> c.kept = ["one", "three"]
    >>> c.arrayed("gib")
    Energies(genre='gib', filenames=['one' 'three'], values=[-123.5 -123.6], t=298.15)
    >>> c.arrayed("gib", full=True) 
    Energies(genre='gib', filenames=['one' 'three' 'four'], values=[-123.5 -123.6 -123.7], t=298.15)
    >>> c.arrayed("gib", t=1111)             
    Energies(genre='gib', filenames=['one' 'three'], values=[-123.5 -123.6], t=1111)

You can also instantiate any data array directly, providing data by yourself.

.. code-block:: python

    >>> from tesliper import Energies
    >>> Energies(
    ...     genre='gib', 
    ...     filenames=['one' 'three'], 
    ...     values=[-123.5 -123.6]
    ... )
    Energies(genre='gib', filenames=['one' 'three'], values=[-123.5 -123.6], t=298.15)


Data validation
'''''''''''''''

On instantiation of a data array class, *values* provided to its constructor are
transformed to the ``numpy.ndarray`` of the appropriate type. If this cannot be done
due to the incompatibility of type of *values* elements and data array's ``dtype``,
an exception is raised. However, ``tesliper`` will try to convert given values to the
target type, if possible.

.. code-block:: python

    >>> from tesliper import IntegerArray
    >>> arr = IntegerArray(genre="example", filenames=["one"], values=["1"])
    >>> arr
    IntegerArray(genre="example", filenames=["one"], values=[1])
    >>> type(arr.values)
    <class 'numpy.ndarray'>

    >>> IntegerArray(genre="example", filenames=["one"], values=["1.0"])
    Traceback (most recent call last):
    ...
    ValueError: invalid literal for int() with base 10: '1.0'

    >>> IntegerArray(genre="example", filenames=["one"], values=[None])
    Traceback (most recent call last):
    ...
    TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType

Also *values* size is checked: its first dimension must be of the same size, as the
number of entries in the *filenames*, otherwise :exc:`ValueError` is raised.

.. code-block:: python

    >>> IntegerArray(genre="example", filenames=["one"], values=[1, 2])
    Traceback (most recent call last):
    ...
    ValueError: values and filenames must have the same shape up to 1 dimensions. Arrays of shape (2,) and (1,) were given.

:exc:`.InconsistentDataError` exception is raised when *values* are multidimensional,
but provide uneven number of entries for each conformer (*values* are a jagged array).

.. code-block:: python

    >>> IntegerArray(genre="example", filenames=["one", "two"], values=[[1, 2], [3]])
    Traceback (most recent call last):
    ...
    InconsistentDataError: IntegerArray of example genre with unequal number of values for conformer requested.

This behavior may be suppressed, if the instance is initiated with
``allow_data_inconsistency=True`` keyword parameter. In such case no exception is raised
if numbers of entries doesn't match, and jagged arrays will be turned into
``numpy.ma.masked_array`` instead of ``numpy.ndarray``, if it is possible.

.. code-block:: python

    >>> IntegerArray(
    ...     genre="example", 
    ...     filenames=["one"], 
    ...     values=[1, 2],
    ...     allow_data_inconsistency=True
    ... )
    IntegerArray(genre="genre", filenames=["one"], values=[1,2], allow_data_incosistency=True)

    >>> IntegerArray(
    ...     genre="example", 
    ...     filenames=["one", "two"], 
    ...     values=[[1, 2], [3]],
    ...     allow_data_inconsistency=True
    ... )
    IntegerArray(genre='genre', filenames=['one' 'two'], values=[[1 2]
     [3 --]], allow_data_inconsistency=True)

Some data array classes validate also other data provided to its constructor, e.g.
:class:`.Geometry` checks if *atoms* provides an atom specification for each atom in the
conformer.

.. note::
    Each validated field is actually a :class:`.ArrayProperty` or its subclass under the
    hood, which provides the validation mechanism.

Available data arrays
'''''''''''''''''''''

Data arrays provided by ``tesliper`` are listed below in categories, along with a short
description and with a list of data genres that are associated with a particular data
array class. More information about a ``DataArray``-like class of interest may be learn
in the :mod:`API reference <tesliper.glassware.arrays>`.


Generic types
"""""""""""""

Simple data arrays, that hold a data of particular type. They do not provide any
functionality beside initial data validation. They are used by ``tesliper`` for
segregation of simple data an as a base classes for other data arrays (concerns mostly
:class:`.FloatArray`).

:class:`.IntegerArray`
    For handling data of ``int`` type.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - charge
          - multiplicity

:class:`.FloatArray`
    For handling data of ``float`` type.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - zpecorr
          - tencorr
          - entcorr
          - gibcorr

:class:`.BooleanArray`
    For handling data of ``bool`` type.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - normal_termination
          - optimization_completed

:class:`.InfoArray`
    For handling data of ``str`` type.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - command
          - stoichiometry

Spectral data
"""""""""""""

Each data array in this category provides a *freq* or *wavelen* attribute, also
accessible by their convenience aliases *frequencies* and *wavelengths*. These
attributes store an information about frequency or wavelength that the particular
spectral value is associated with (x-axis value of the center of the band).

Activities genres, that are the genres that may be used to simulate the spectrum, also
provide a *calculate_spectra()* method for this purpose (see
:meth:`.VibrationalActivities.calculate_spectra`,
:meth:`.ScatteringActivities.calculate_spectra`, and
:meth:`.ElectronicActivities.calculate_spectra`), as well as a
:attr:`~.SpectralActivities.intensities` property that calculates a theoretical
intensity for each activity value. A convince :attr:`~.SpectralActivities.spectra_name`
property may be used to get the name of spectra pseudo-genre calculated with particular
activities genre.

.. code-block:: python

    >>> act = c["dip"]
    >>> act.spectra_name
    "ir"
    >>> from tesliper import lorentzan
    >>> spc = act.calculate_spectra(
    ...     start=200,  # cm^(-1)
    ...     stop=1800,  # cm^(-1)
    ...     step=1,     # cm^(-1)
    ...     width=5,    # cm^(-1)
    ...     fitting=lorentzan
    )
    >>> type(spc), spc.genre
    (<class 'tesliper.glassware.spectra.Spectra'>, 'ir')

:class:`.VibrationalData`
    For handling vibrational (IR and VCD related) data that is not a spectral activity.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - mass
          - frc
          - emang

:class:`.ScatteringData`
    For handling scattering (Raman and ROA related) data that is not a spectral activity.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - depolarp
          - depolaru
          - depp
          - depu
          - alpha2
        * - beta2
          - alphag
          - gamma2
          - delta2
          - cid1
        * - cid2
          - cid3
          - rc180
          -
          -

:class:`.ElectronicData`
    For handling electronic (UV and ECD related) data that is not a spectral activity.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - eemang

:class:`.VibrationalActivities`
    For handling vibrational (IR and VCD related) spectral activity data.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - iri
          - dip
          - rot

:class:`.ScatteringActivities`
    For handling scattering (Raman and ROA related) spectral activity data.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - ramanactiv
          - ramact
          - raman1
          - roa1
        * - raman2
          - roa2
          - raman3
          - roa3

:class:`.ElectronicActivities`
    For handling electronic (UV and ECD related) spectral activity data.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - vdip
          - ldip
          - vrot
          - lrot
          - vosc
          - losc

Other data arrays
"""""""""""""""""

:class:`.FilenamesArray`
    Special case of :class:`.DataArray`, holds only filenames. *values* property
    returns same as *filenames* and ignores any value given to its setter.
    The only genre associated with this class is *filenames* pseudo-genre.

:class:`.Bands`
    Special kind of data array for band values, to which spectral data or activities
    correspond. Provides an easy way to convert values between their different
    representations: frequency, wavelength, and excitation energy. Also allows to easily
    locate conformers with imaginary frequencies.

    .. code-block:: python

        >>> arr = Bands(
        ...     genre="freq",
        ...     filenames=["one", "two", "three"],
        ...     values=[[-15, -10, 105], [30, 123, 202], [-100, 12, 165]]
        ... )
        >>> arr.imaginary
        array([2, 0, 1])
        >>> arr.find_imaginary()
        {'one': 2, 'three': 1}

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - freq
          - wavelen
          - ex_en

:class:`.Energies`
    For handling data about the energy of conformers. Provides an easy way of
    calculating Boltzmann distribution-based population of conformers *via* a
    :attr:`~.Energies.populations` property.

    .. code-block:: python

        >>> arr = Energies(
        ...     genre="gib", 
        ...     filenames=["one", "two", "three"], 
        ...     values=[-123.505977, -123.505424, -123.506271]
        ... )
        >>> arr.deltas  # difference from lowest energy in kcal/mol
        array([0.18448779, 0.53150055, 0.        ])
        >>> arr.populations
        array([0.34222796, 0.19052561, 0.46724643])

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - scf
          - zpe
          - ten
          - ent
          - gib

:class:`.Transitions`
    For handling information about electronic transitions from ground
    to excited state contributing to each band.

    Data is stored in three attributes: :attr:`~.Transitions.ground`,
    :attr:`~.Transitions.excited`, and :attr:`~.Transitions.values`, which are
    respectively: list of ground state electronic subshells, list of excited state
    electronic subshells, and list of coefficients of transitions from corresponding
    ground to excited subshell. Each of these arrays is of shape (conformers, bands,
    max_transitions), where 'max_transitions' is a highest number of transitions
    contributing to single band across all bands of all conformers.

    Allows to easily calculate contribution of each transition using
    :attr:`~.Transitions.contribution` and to find which transition contributes the most
    to the particular transition with :attr:`~.Transitions.highest_contribution`. 

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - transitions

:class:`.Geometry`
    For handling information about geometry of conformers.

    .. list-table:: Genres associated with this class:
        :width: 100%

        * - last_read_geom
          - input_geom
          - optimized_geom

Writing to disk
---------------

:class:`.Tesliper` object provides an easy, but not necessarily a flexible way of
writing calculated and extracted data to disk. If your process requires more flexibility
in this matter, you may use ``teslier``\s writer objects directly. This will allow you
to adjust how generated files are named and will give you more control over what is
exported.

Writer classes
''''''''''''''

A writer object may be created using a :func:`.writer` factory function. It expects a
string parameter, that specifies a desired format for data export. ``tesliper`` provides
writers for ``"txt"``, ``"csv"``, ``"xlsx"``, and ``"gjf"`` file formats. The second
mandatory parameter is a *destination*: the (existing) directory to which files should
be written. Just like writing methods of :class:`.Tesliper` object, the function also
takes a *mode* parameter that defines what should happen if any file already exists. Any
additional keyword parameters are forwarded to the writer object constructor.

.. code-block:: python

    >>> from tesliper import writer
    >>> wrt = writer("txt", "/path/to/dir")
    >>> type(wrt)
    <class 'tesliper.writing.txt_writer.TxtWriter'>

    >>> wrt = writer("txt", "/doesnt/exists")
    Traceback (most recent call last):
    ...
    FileNotFoundError: Given destination doesn't exist or is not a directory.

.. note::

    :func:`.writer` factory function is used by ``tesliper`` mostly to provide a dynamic
    access to the writer class most recently registered (on class definition) to handle
    a particular format. This is useful when you modify an existing writer class or
    provide a new one.

You can also create any of the writer objects directly, by importing and instantiating
its class. The four available writer classes are listed below with a short comment.
For more information on which methods they implement and how to use them, refer to the
relevant API documentation.

.. code-block:: python

    from tesiper import TxtWriter
    wrt = TxtWriter(destination="/path/to/dir")

:class:`.TxtWriter`
    Generates human-readable text files.

:class:`.CsvWriter`
    Generates files in CSV format with optional headers. Allows for the same level of
    output format customization as Python's :class:`csv.writer` (supports specification
    of *dialect* and other formatting parameters).

:class:`.XlsxWriter`
    Instead of generating multiple files, creates a single .xlsx file and a variable
    number of spreadsheets inside it.

:class:`.GjfWriter`
    Allows to create input files for new calculation job in Gaussian software.


``write()`` and other methods
'''''''''''''''''''''''''''''

Writer objects expect data they receive to be a :class:`.DataArray`-like instances. Each
writer object provides a :meth:`~.WriterBase.write` method for writing arbitrary data
arrays to disk. This method dispatches received data arrays to appropriate writing
methods, based on their type. You are free to use either :meth:`~.WriterBase.write` for
easily writing a number of data genres in batch, or other methods for more control. The
table below lists these methods, along with a brief description and
:class:`.DataArray`-like object, for which the method will be called by writer's
:meth:`~.WriterBase.write` method.

.. list-table:: Methods used to write certain data
    :header-rows: 1

    * - Writer's Method
      - Description
      - Supported arrays
      - Created files
    * - :meth:`~.WriterBase.generic`
      - Generic data: any genre that provides one value for each conformer.
      - :class:`.DataArray`, :class:`.IntegerArray`,
        :class:`.FloatArray`, :class:`.BooleanArray`, :class:`.InfoArray`.
      - one
    * - :meth:`~.WriterBase.overview`
      - General information about conformers: energies, imaginary frequencies,
        stoichiometry.
      - :class:`.Energies`
      - one
    * - :meth:`~.WriterBase.energies`
      - Detailed information about conformers' relative energy,
        including calculated populations
      - :class:`.Energies`
      - for each genre
    * - :meth:`~.WriterBase.single_spectrum`
      - A spectrum - calculated for single conformer or averaged.
      - :class:`.SingleSpectrum`
      - one
    * - :meth:`~.WriterBase.spectral_data`
      - Data related to spectral activity, but not convertible to spectra.
      - :class:`.VibrationalData`, :class:`.ScatteringData`, :class:`.ElectronicData`
      - for each conformer
    * - :meth:`~.WriterBase.spectral_activities`
      - Data that may be used to simulate conformers' spectra.
      - :class:`.VibrationalActivities`, :class:`.ScatteringActivities`,
        :class:`.ElectronicActivities`
      - for each conformer
    * - :meth:`~.WriterBase.spectra`
      - Spectra for multiple conformers.
      - :class:`.Spectra`
      - for each conformer
    * - :meth:`~.WriterBase.transitions`
      - Electronic transitions from ground to excited state, contributing to each band.
      - :class:`.Transitions`
      - for each conformer
    * - :meth:`~.WriterBase.geometry`
      - Geometry (positions of atoms in space) of conformers.
      - :class:`.Geometry`
      - for each conformer

:meth:`~.WriterBase.spectral_data` and :meth:`~.WriterBase.spectral_activities` methods
need some clarification. They will create one file for each conformer in given data
arrays, with data from each provided data array joined in the conformer's file. It's
important to remember, that only *values* from each data array are displayed, contrary
to band values, which are displayed only once, as provided with the *band* parameter.
Consequently, mixing vibrational and scattering data with a custom *name_template* is
fine, but mixing either of those with electronic data in a single call is not possible.

.. warning::

    You need to make sure that data contained in ``DataArray``-like objects cover the
    same set of conformers, when passing multiple data array objects to the
    :meth:`~.WriterBase.write` method or any other writing method. Passing two data
    arrays with data for different sets of conformers may produce files with corrupted
    data or fail silently. :meth:`.Conformers.trim_incomplete` trimming method may be
    helpful in preventing such fails.

Not all writer objects implement each of these writing methods, e.g.
:class:`.GjfWriter`, that allows to create Gaussian input files, only implements
:meth:`~.GjfWriter.geometry` method (because export of, e.g. a calculated spectrum as a
Gaussian input would be pointless). Trying to ``write()`` a data array that should be
written by a method that is not implemented, or calling such method directly, will raise
a :exc:`NotImplementedError`.

Naming files
''''''''''''

Usually, calling any of writing methods will produce multiple files in the *destination*
directory: one for each given genre, each conformer, etc. ``tesliper`` provides a
reasonable naming scheme for these files, but you can modify it, by providing your own
*name_template*\s in place of the default ones. To do this you will need to call
desired writing methods directly, instead of using :meth:`~.WriterBase.write`.

Each writing method uses a value of *name_template* parameter given to the method call
to create a filename for each file it generates. *name_template* should be a string that
contains (zero, one, or more) label identifiers in form of ``${identifier}``. These
identifiers will be substituted to produce a final filename. Available identifiers and
their meaning are as follows:

| ``${ext}`` - appropriate file extension;
| ``${conf}`` - name of the conformer;
| ``${num}`` - number of the file according to internal counter;
| ``${genre}`` - genre of exported data;
| ``${cat}`` - category of produced output;
| ``${det}`` - category-specific detail.

The ``${ext}`` identifier is filled with the value of Writers ``.extension`` attribute,
which value is also used to identify a writer class: ``"txt"``, ``"csv"``, etc. Other
values are provided by the particular writing method.

.. code-block:: python

    from tesiper import Tesliper, writer
    tslr = Tesliper(input_dir="/project/input")
    ...  # data extracted and trimmed
    # tslr.conformers.kept_keys() == {"conf_one", "conf_four"}
    freq, dip, rot = tslr["freq"], tslr["dip"], tslr["rot"]
    wrt = writer("txt", "/project/default")
    wrt.spectral_activities(band=freq, data=[dip, rot])
    wrt = writer("txt", "/project/custom")
    wrt.spectral_activities(
        band=freq, data=[dip, rot],
        name_template="name_${num}_${genre}.xy"
    )

.. code-block:: none
    :caption: contents of ``/project``

    .
    ├───input
    │   └─── ...
    ├───default
    │   ├───conf_one.activities-vibrational.txt
    │   └───conf_four.activities-vibrational.txt
    └───custom
        ├───name_1_freq.xy
        └───name_2_freq.xy

.. _conventions:

Conventions and Terms
=====================

Reading and writing
-------------------

``tesliper`` was designed to deal with multiple conformers of a single molecule. It
identifies conformers using a stem of an extracted file (i.e. its filename without
extension). When files with identical names (save extension) are extracted in course of
subsequent :meth:`.Tesliper.extract` calls (or in recursive extraction, see method's
documentation), they are treated as the same conformer. This enables to join data
from subsequent calculations steps, e.g. geometry optimization, vibrational spectra
simulation, and electronic spectra simulation.

.. note::
    If specific data genre is available from more than one file, only recently
    extracted values will be stored.

Also, writing extracted and calculated data to files is done in batch, as usually
multiple files are produced. Hence, ``tesliper`` will chose names for these files
automatically, only allowing to specify output directory (as
:attr:`.Tesliper.output_dir` attribute). If you need more control over this process,
you will need to use one of the writer objects directly. These are easily available
*via* the :func:`.writer_base.writer` factory function.

Handling data
-------------

``tesliper`` stores multiple data entries of various types for each conformer. To
prevent confusion with Python's data ``type`` and with data itself, ``tesliper`` refers
to specific kinds of data as *genres*. Genres in code are represented by specific
strings, used as identifiers. To learn about data genres known to ``tesliper``, see
documentation for :class:`.GaussianParser`, which lists them.

.. note::
    Given the above, you may wonder why is it *genres* and not just *kinds* of data
    then? The reason is that naming things is hard (one of the only two hard things in
    Computer Science, as `Phil Karlton said
    <https://www.karlton.org/2017/12/naming-things-hard/>`_). As of time of deciding on
    this name, I did not come up with the second one. Hopefully, this small oddity will
    not bother you too much.

``tesliper`` may not work properly when used to process data concerning different
molecules (i.e. having different number of atoms, different number of degrees of
freedom, etc.). If you want to use it for such purpose anyway, you may set
:attr:`Tesliper.conformers.allow_data_inconsistency <
.Conformers.allow_data_inconsistency>` to ``True``. ``tesliper`` will then stop
complaining and try to do its best.

Glossary
--------

.. glossary::

    genre
        A specific kind of data, e.g. SCF energy, dipole strengths, atoms' positions in
        space, or command used for calculations. Represented in code by a short string. Not
        to be confused with Python's data ``type``. See :ref:`available genres`.

    trimming
        Internally marking certain conformers as *not kept*. ``tesliper`` provides an easy
        way to trim conformers to user's needs, see :ref:`filtering conformers`.

    kept
        Conformers may be internally marked as *kept* or *not kept*. *Kept* conformers
        will be normally processed by ``tesliper``, *not kept* conformers will be ignored.
        See :attr:`.Conformers.kept`.

    arrayed
        About data turned into an instance of :class:`.DataArray`-like object, usually by
        :class:`.Conformers`' method of the same name. See :meth:`.Conformers.arrayed`.

    data array
        Type of objects used by ``tesliper`` to handle data read from multiple conformers.
        The same data array class may be used to represent more than one genre. Sometimes
        referred to as :class:`.DataArray`-like classes or objects. See :mod:`.arrays`.

    data inconsistency
        An event of data having non-uniform properties, e.g. when number of values doesn't
        match number of conformers, or when some conformers provide a different number of
        values than other conformers for a particular data genre. See :mod:`.array_base`.
Scripting with ``tesliper``
===========================

This part discusses basics of using Python API. For tutorial on using a Graphical
Interface, see :ref:`gui`.

``tesliper`` provides a :class:`.Tesliper` class as a main entry point to its
functionality. This class allows you to easily perform any typical task: read and write
files, filter data, calculate and average spectra. It is recommended to read
:ref:`conventions` to get a general idea of what to expect. The next paragraphs will
introduce you to basic use cases with examples and explanations. The examples do not use
real data, but simplified mockups, to not obscure the logic presented.

.. code-block:: python

    from tesliper import Tesliper

    # extract data from Gaussian output files
    tslr = Tesliper(input_dir="./opt_and_freq")
    tslr.extract()

    # conditional filtering of conformers
    tslr.conformers.trim_non_normal_termination()
    tslr.conformers.trim_not_optimized()
    tslr.conformers.trim_imaginary_frequencies()
    tslr.conformers.trim_to_range("gib", maximum=10, attribute="deltas")
    tslr.conformers.trim_rmsd(threshold=1.0, window_size=0.5, energy_genre="gib")

    # calculate and average spectra, export data
    tslr.calculate_spectra()
    tslr.average_spectra()
    tslr.export_energies(fmt="txt")
    tslr.export_averaged(fmt="csv")

Reading files
-------------

After importing :class:`.Tesliper` class, we instantiate it with an *input_dir*
parameter, which is a path to the directory containing output files from quantum
chemical calculations software. You may also provide an *output_dir* parameter, defining
where ``tesliper`` should write the files it generates. Both of those parameters are
optional and default to the current working directory, if omitted. You may also provide
a *wanted_files* parameter, which should be a list of filenames that :class:`.Tesliper`
should parse, ignoring any other files present in *input_dir*. Omitting *wanted_files*
means that no file should be ignored.

.. note::
    
    :class:`.Tesliper` accepts also *quantum_software* parameter, which is a hint for
    ``tesliper`` on how it should parse output files it reads. However, only Gaussian
    software is supported out-of-the-box, and ``quantum_software="gaussian"`` is a
    default value. If you wish to use ``tesliper`` to work with another qc package,
    you will need to define a custom parser that subclasses the :class:`.ParserBase`
    class. Refer to its documentation for more information.

You can extract data from the files in *output_dir* using :meth:`.Tesliper.extract`
method. :meth:`.Tesliper.extract` respects *input_dir* and *wanted_files* given to
:class:`.Tesliper`, but *path* and *wanted_files* parameters provided to the method call
will take precedence. If you would like to read files in the whole directory tree, you
may perform a recursive extraction, using ``extract(recursive=True)``. So assuming a
following directory structure::

    project
    ├── optimization
    │   ├── conf_one.out
    │   └── conf_two.out
    │   └── conf_three.out
    └── vibrational
        ├── conf_one.out
        └── conf_two.out
        └── conf_three.out

you could use any of the following to get the same effect.

.. code-block:: python

    # option 1: change *input_dir*
    tslr = Tesliper(input_dir="./project/optimization")
    tslr.extract()
    tslr.input_dir = "./project/vibrational"
    tslr.extract()

    # option 2: override *input_dir* only for one call
    tslr = Tesliper(input_dir="./project/optimization")
    tslr.extract()
    tslr.extract(path="./project/vibrational")

    # option 3: read the whole tree
    tslr = Tesliper(input_dir="./project")
    tslr.extract(recursive=True)


``tesliper`` will try to guess the extension of files it should parse: e.g. Gaussian
output files may have ".out" or ".log" extension. If those are mixed in the source
directory, an exception will be raised. You can prevent this by providing the
*extension* parameter, only files with given extension will be parsed.

.. code-block:: none

    project
    ├── conf_one.out
    └── conf_two.log
    
.. code-block:: python

    tslr = Tesliper(input_dir="./project")
    tslr.extract()  # raises ValueError
    tslr.extract(extension="out")  # ok

.. _filtering conformers:

Filtering conformers
--------------------

:meth:`.Tesliper.extract` will read and parse files it thinks are output files of the
quantum chemical software and update a :attr:`.Tesliper.conformers` internal data
storage. It is a ``dict``-like :class:`.Conformers` instance, that stores data for each
conformer in a form of an ordinary :class:`dict`. This inner dict uses :term:`genre`
names as keys and data as values (the form of which depends on the genre itself).
:class:`.Conformers` provide a number of methods for filtering conformers it knows,
allowing to easily hide data that should excluded from further analysis. ``tesliper``
calls this process a *trimming*. The middle part of the first code snippet are example
of trimming conformers:

.. code-block:: python

    tslr.conformers.trim_non_normal_termination()
    tslr.conformers.trim_not_optimized()
    tslr.conformers.trim_imaginary_frequencies()
    tslr.conformers.trim_to_range("gib", maximum=10, attribute="deltas")
    tslr.conformers.trim_rmsd(threshold=1.0, window_size=0.5, energy_genre="gib")

As you may suspect, :meth:`~.Conformers.trim_non_normal_termination` hides data from
calculations that did not terminate normally, :meth:`~.Conformers.trim_not_optimized`
hides data from conformers that are not optimized, and
:meth:`~.Conformers.trim_imaginary_frequencies` hides data from conformers that have at
least one imaginary frequency. More trimming methods is described :ref:`below
<trimming>`.

Conformers hidden are :term:`not kept <kept>`.
Information about which conformers are *kept* and *not kept* is stored in
:attr:`.Conformers.kept` attribute, which may also be manipulated more directly. More on
this topic will be :ref:`explained later <manipulating kept>`.

As mentioned earlier, :class:`Tesliper.conformers <.Conformers>` is a dict-like
structure, and as such offers a typical functionality of Python's ``dict``\s. However,
checking for presence with ``conf in tslr.conformers`` or requesting a view with
standard :meth:`~dict.keys`, :meth:`~dict.values`, or :meth:`~dict.items` will operate
on the whole data set, ignoring any trimming applied earlier. :class:`.Conformers` class
offers additional :meth:`~.Conformers.kept_keys`, :meth:`~.Conformers.kept_values`, and
:meth:`~.Conformers.kept_items` methods, that return views that acknowledge trimming.

.. _trimming:

Trimming methods
''''''''''''''''

There is a number of those methods available for you, beside those mentioned above.
Below you will find them listed with a short summary and a link to a more comprehensive
explanation in the method's documentation.

:meth:`~.Conformers.trim_incomplete`
    Filters out conformers that doesn't contain data for as many expected genres as
    other conformers.

:meth:`~.Conformers.trim_imaginary_frequencies`
    Filters out conformers that contain imaginary frequencies (any number of negative
    frequency values).

:meth:`~.Conformers.trim_non_matching_stoichiometry`
    Filters out conformers that have different stoichiometry than expected.

:meth:`~.Conformers.trim_not_optimized`
    Filters out conformers that failed structure optimization.

:meth:`~.Conformers.trim_non_normal_termination`
    Filters out conformers, which calculation job did not terminate normally (was
    erroneous or interrupted).

:meth:`~.Conformers.trim_inconsistent_sizes`
    Filters out conformers that have iterable data genres in different size than most
    conformers. Helpful when :exc:`.InconsistentDataError` occurs.

:meth:`~.Conformers.trim_to_range`
    Filters out conformers that have a value of some specific data or property outside
    of the given range, e.g. their calculated population is less than 0.01.

:meth:`~.Conformers.trim_rmsd`
    Filters out conformers that are identical to another conformer, judging by a given
    threshold of the root-mean-square deviation of atomic positions (RMSD).

:meth:`~.Conformers.select_all`
    Marks all conformers as :term:`kept`.

:meth:`~.Conformers.reject_all`
    Marks all conformers as :term:`not kept <kept>`.

.. _manipulating kept:

Manipulating ``Conformers.kept``
''''''''''''''''''''''''''''''''

Information, which conformer is :term:`kept` and which is not, is stored in the
:attr:`Conformers.kept` attribute. It is a list of booleans, one for each conformer
stored, defining which conformers should be processed by ``tesliper``.

.. code-block:: python

    # assuming "conf_two" has imaginary frequencies
    tslr.conformers.trim_imaginary_frequencies()
    tslr.conformers.kept == [True, False, True]  # True
    tslr.export_data(["genres", "to", "export"])
    # only files for "conf_one" and "conf_three" are generated

:attr:`.Conformers.kept` may be modified using trimming methods described :ref:`earlier
<trimming>`, but also more directly: by setting it to a new value. Firstly, it is the
most straightforward to just assign a new list of boolean values to it. This list should
have the same number of elements as the number of conformers contained. A
:exc:`ValueError` is raised if it doesn't.

.. code-block:: python

    >>> tslr.conformers.kept
    [True, True, True]
    >>> tslr.conformers.kept = [False, True, False]
    >>> tslr.conformers.kept
    [False, True, False]
    >>> tslr.conformers.kept = [False, True, False, True]
    Traceback (most recent call last):
    ...
    ValueError: Must provide boolean value for each known conformer.
    4 values provided, 3 excepted.

Secondly, list of filenames of conformers intended to be *kept* may be given. Only these
conformers will be *kept*. If given filename is not in the underlying
:class:`tslr.conformers <.Conformers>`' dictionary, :exc:`KeyError` is raised.

.. code-block:: python

    >>> tslr.conformers.kept = ['conf_one']
    >>> tslr.conformers.kept
    [True, False, False]
    >>>  tslr.conformers.kept = ['conf_two', 'other']
    Traceback (most recent call last):
    ...
    KeyError: Unknown conformers: other.

Thirdly, list of integers representing conformers' indices may be given.
Only conformers with specified indices will be *kept*. If one of given integers
can't be translated to conformer's index, IndexError is raised. Indexing with
negative values is not supported currently.

.. code-block:: python

    >>> tslr.conformers.kept = [1, 2]
    >>> tslr.conformers.kept
    [False, True, True]
    >>> tslr.conformers.kept = [2, 3]
    Traceback (most recent call last):
    ...
    IndexError: Indexes out of bounds: 3.

Fourthly, assigning ``True`` or ``False`` to this attribute will mark all
conformers as *kept* or *not kept* respectively.

.. code-block:: python

    >>> tslr.conformers.kept = False
    >>> tslr.conformers.kept
    [False, False, False]
    >>> tslr.conformers.kept = True
    >>> tslr.conformers.kept
    [True, True, True]

.. warning::

    List of *kept* values may be also modified by setting its elements to ``True`` or
    ``False``. It is advised against, however, as a mistake such as
    ``tslr.conformers.kept[:2] = [True, False, False]`` will break some functionality by
    forcibly changing size of :attr:`tslr.conformers.kept <.Conformers.kept>` list.

Trimming temporarily
''''''''''''''''''''

:class:`.Conformers` provide two convenience context managers for temporarily trimming
its data: :attr:`~.Conformers.untrimmed` and :meth:`~.Conformers.trimmed_to`. The first
one will simply undo any trimming previously done, allowing you to operate on the full
data set or apply new, complex trimming logic. When Python exits
:attr:`~.Conformers.untrimmed` context, previous trimming is restored.

.. code-block:: python

    >>> tslr.conformers.kept = [False, True, False]
    >>> with tslr.conformers.untrimmed:
    >>>     tslr.conformers.kept
    [True, True, True]
    >>> tslr.conformers.kept
    [False, True, False]

The second one temporarily applies an arbitrary trimming, provided as a parameter to the
:meth:`~.Conformers.trimmed_to` call. Any value normally accepted by :attr:`.Conformers.kept`
may be used here.

.. code-block:: python

    >>> tslr.conformers.kept = [True, True, False]
    >>> with tslr.conformers.trimmed_to([1, 2]):
    >>>     tslr.conformers.kept
    [False, True, True]
    >>> tslr.conformers.kept
    [True, True, False]


.. tip::

    To trim conformers temporarily without discarding a currently applied trimming, you
    may use:

    .. code-block:: python

        with tslr.conformers.trimmed_to(tslr.conformers.kept):
            ...  # temporary trimming upon the current one


Simulating spectra
------------------

To calculate a simulated spectra you will need to have spectral activities extracted.
These will most probably come from a *freq* or *td* Gaussian calculation job, depending
on a genre of spectra you would like to simulate. ``tesliper`` can simulate IR, VCD, UV,
ECD, Raman, and ROA spectra, given the calculated values of conformers' optical
activity. When you call :meth:`.Tesliper.calculate_spectra` without any parameters, it
will calculate spectra of all available genres, using default activities genres and
default parameters, and store them in the :attr:`.Tesliper.spectra` dictionary. Aside
form this, the spectra calculated are returned by the method.

You can calculate a specific spectra genres only, by providing a list of their names as
a parameter to the :meth:`.Tesliper.calculate_spectra` call. Also in this case a default
activities genres and default parameters will be used to calculate desired spectra, see
`Activities genres`_ and `Calculation parameters`_ below to learn how this can be
customized. 

.. code-block:: python

    ir_and_uv = tslr.calculate_spectra(["ir", "uv"])
    assert ir_and_uv["ir"] is tslr.spectra["ir"]

Calculation parameters
''''''''''''''''''''''

``tesliper`` uses `Lorentzian <https://en.wikipedia.org/wiki/Cauchy_distribution>`_ or
`Gaussian <https://en.wikipedia.org/wiki/Normal_distribution>`_ fitting function to
simulate spectra from corresponding optical activities values. Both of these require to
specify a desired width of peak, as well as the beginning, end, and step of the abscissa
(x-axis values). If not told otherwise, ``tesliper`` will use a default values for these
parameters and a default fitting function for a given spectra genre. These default
values are available *via* :attr:`.Tesliper.standard_parameters` and are as follows.

.. table:: Default calculation parameters

    +-----------+--------------------------------+--------------------------------+
    | Parameter |  IR, VCD, Raman, ROA           | UV, ECD                        |
    +===========+================================+================================+
    | width     | 6 [:math:`\mathrm{cm}^{-1}`]   | 0.35 [:math:`\mathrm{eV}`]     |
    +-----------+--------------------------------+--------------------------------+
    | start     | 800 [:math:`\mathrm{cm}^{-1}`] | 150 [:math:`\mathrm{nm}`]      |
    +-----------+--------------------------------+--------------------------------+
    | stop      | 2900 [:math:`\mathrm{cm}^{-1}`]| 800 [:math:`\mathrm{nm}`]      |
    +-----------+--------------------------------+--------------------------------+
    | step      | 2 [:math:`\mathrm{cm}^{-1}`]   | 1 [:math:`\mathrm{nm}`]        |
    +-----------+--------------------------------+--------------------------------+
    | fitting   | :func:`.lorentzian`            | :func:`.gaussian`              |
    +-----------+--------------------------------+--------------------------------+

You can change the parameters used for spectra simulation by altering values in the
:attr:`.Tesliper.parameters` dictionary. It stores a ``dict`` of parameters' values for
each of spectra genres ("ir", "vcd", "uv", "ecd", "raman", and "roa"). *start*, *stop*,
and *step* expect its values to by in :math:`\mathrm{cm}^{-1}` units for vibrational and
scattering spectra, and :math:`\mathrm{nm}` units for electronic spectra. *width*
expects its value to be in :math:`\mathrm{cm}^{-1}` units for vibrational and scattering
spectra, and :math:`\mathrm{eV}` units for electronic spectra. *fitting* should be a
callable that may be used to simulate peaks as curves, preferably one of:
:func:`.gaussian` or :func:`.lorentzian`.

.. code-block:: python

    # change parameters' values one by one 
    tslr.parameters["uv"]["step"] = 0.5
    tslr.parameters["uv"]["width"] = 0.5

    tslr.parameters["vcd"].update(  # or with an update
        {"start": 500, "stop": 2500, "width": 2}
    )

    # "fitting" should be a callable
    from tesliper import lorentzian
    tslr.parameters["uv"]["fitting"] = lorentzian


.. table:: Descriptions of parameters

    +-----------+----------------------+-------------------------------------------+
    | Parameter |  ``type``            | Description                               |
    +===========+======================+===========================================+
    | width     | ``float`` or ``int`` | the beginning of the spectral range       |
    +-----------+----------------------+-------------------------------------------+
    | start     | ``float`` or ``int`` | the end of the spectral range             |
    +-----------+----------------------+-------------------------------------------+
    | stop      | ``float`` or ``int`` | step of the abscissa                      |
    +-----------+----------------------+-------------------------------------------+
    | step      | ``float`` or ``int`` | width of the peak                         |
    +-----------+----------------------+-------------------------------------------+
    | fitting   | ``Callable``         | function used to simulate peaks as curves |
    +-----------+----------------------+-------------------------------------------+

.. warning::

    When modifying :attr:`.Tesliper.parameters` be careful to not delete any of the
    key-value pairs. If you need to revert to standard parameters' values, you can just
    reassign them to :attr:`.Tesliper.standard_parameters`.
        
    .. code-block:: python
        
        tslr.parameters["ir"] = {
        ...     "start": 500, "stop": 2500, "width": 2
        ... }  # this will cause problems!
        # revert to default values
        tslr.parameters["ir"] = tslr.standard_parameters["ir"]

Activities genres
'''''''''''''''''

Instead of specifying a spectra genre you'd like to get, you may specify an activities
genre you prefer to use to calculate a corresponding spectrum. The table below summarizes
which spectra genres may be calculated from which activities genres.

.. list-table:: Spectra and corresponding activities genres
    :header-rows: 1

    * - Spectra
      - Default activity
      - Other activities
    * - IR
      - dip
      - iri
    * - VCD
      - rot
      - 
    * - UV
      - vosc
      - losc, vdip, ldip
    * - ECD
      - vrot
      - lrot
    * - Raman
      - raman1
      - ramact, ramanactiv, raman2, raman3
    * - ROA
      - roa1
      - roa2, roa3

.. warning::

    If you provide two different genres that map to the same spectra genre, only one of
    them will be accessible, the other will be thrown away. If you'd like to compare
    results of simulations using different genres, you need to store the return value
    of :meth:`.Tesliper.calculate_spectra` call.

    .. code-block:: python

        >>> out = tslr.calculate_spectra(["vrot", "lrot"])
        >>> list(out.keys())  # only one representation returned
        ["ecd"]
        >>> velo = tslr.calculate_spectra(["vrot"])
        >>> length = tslr.calculate_spectra(["lrot"])
        >>> assert not velo["ecd"] == length["ecd"]  # different

Averaging spectra
'''''''''''''''''

Each possible conformer contributes to the compound's spectrum proportionally to it's
population in the mixture. ``tesliper`` can calculate conformers' population from their
relative energies, using a technique called `Boltzmann distribution
<https://en.wikipedia.org/wiki/Boltzmann_distribution>`_. Assuming that any energies
genre is available (usually at least *scf* energies are), after calculating spectra you
want to simulate, you should call :meth:`.Tesliper.average_spectra` to get the final
simulated spectra.

.. code-block:: python

    >>> averaged = tslr.average_spectra()  # averages available spectra
    >>> assert tslr.averaged is averaged  # just a reference

:meth:`.Tesliper.average_spectra` averages each spectra stored in the
:attr:`.Tesliper.spectra` dictionary, using each available energies genre. Generated
average spectra are stored in :attr:`.Tesliper.averaged` dictionary, using a tuples
of ``("spectra_genre", "energies_genre")`` as keys. :meth:`~.Tesliper.average_spectra`
returns a reference to this attribute.

.. note::

    There is also a :meth:`.Tesliper.get_averaged_spectrum` method for calculating a
    single averaged spectrum using a given spectra genre and energies genre. Value
    returned by this method is not automatically stored.

Temperature of the system
'''''''''''''''''''''''''

Boltzmann distribution depends on the temperature of the system, which is assumed to be
the room temperature, expressed as :math:`298.15\ \mathrm{Kelvin}`
(:math:`25.0^{\circ}\mathrm{C}`). You can change it by setting
:meth:`.Tesliper.temperature` attribute to the desired value. This must be done before
calculation of the average spectrum to have an effect.

.. code-block:: python

    >>> tslr.temperature
    298.15  # default value in Kelvin
    >>> averaged = tslr.average_spectra()  # averages available spectra
    >>> tslr.temperature = 300.0  # in Kelvin
    >>> high_avg = tslr.average_spectra()
    >>> assert not averaged == high_avg  # resulting average is different

.. note::

    :meth:`.Tesliper.temperature` value must be a positive number, absolute zero or
    lower is not allowed, as it would cause problems in calculation of Boltzmann
    distribution. An attempt to set temperature to equal or below :math:`0\ \mathrm{K}`
    will raise a `ValueError`.

Comparing with experiment
-------------------------

The experimental spectrum may be loaded with ``tesliper`` from a text or CSV file. The
software helps you adjust the shift and scale of your simulated spectra to match the
experiment. Unfortunately, ``tesliper`` does not offer broad possibilities when it comes
to mathematical comparison of the simulated spectra and the experimental one. You will
need to use an external library or write your own logic to do that.

Loading experimental spectra
''''''''''''''''''''''''''''

To load an experimental spectrum use :meth:`.Tesliper.load_experimental` method. You
will need to provide a path to the file (absolute or relative to the current
:attr:`.Tesliper.input_dir`) and a genre name of the loaded experimental spectrum. When
the file is read, its content is stored in :attr:`.Tesliper.experimental` dictionary.

.. code-block:: python

    >>> spectrum = tslr.load_experimental("path/to/spectrum.xy", "ir")
    >>> tslr.experimental["ir"] is spectrum
    True

Adjusting calculated spectra
''''''''''''''''''''''''''''

Spectra calculated and loaded from disk with ``tesliper`` are stored as instances of
:class:`.Spectra` or :class:`.SingleSpectrum` classes. Both of them provide a
:meth:`~.SingleSpectrum.scale_to` and :meth:`~.SingleSpectrum.shift_to` methods that
adjust a scale and offset (respectively) to match another spectrum, provided as a
parameter. Parameters found automatically may not be perfect, so you may provide them
yourself, by manually setting :attr:`~.SingleSpectrum.scaling` and
:attr:`~.SingleSpectrum.offset` to desired values.

.. code-block:: python

    >>> spectra = tslr.spectrum["ir"]
    >>> spectra.scaling  # affects spectra.y
    1.0
    >>> spectra.scale_to(tslr.experimental["ir"])
    >>> spectra.scaling
    1.32
    >>> spectra.offset = 50  # bathochromic shift, affects spectra.x

Corrected values may be accessed *via* ``spectra.x`` and ``spectra.y``,
original values may be accessed *via* ``spectra.abscissa`` and ``spectra.values``.


Writing to disk
---------------

Once you have data you care about, either extracted or calculated, you most probably
would like to store it, process it with another software, or visualize it. ``tesliper``
provides a way to save this data in one of the supported formats: CSV, human-readable
text files, or .xlsx spreadsheet files. Moreover, ``tesliper`` may produce Gaussian
input files, allowing you to easily setup a next step of calculations.

Exporting data
''''''''''''''

``tesliper`` provides a convenient shorthands for exporting certain kinds of data:

- :meth:`.Tesliper.export_energies` will export information about conformers' energies
  and their calculated populations;
- :meth:`.Tesliper.export_spectral_data` will export data that is related to spectral
  activity, but cannot be used to simulate spectra;
- :meth:`.Tesliper.export_activities` will export unprocessed spectral activities,
  normally used to simulate spectra;
- :meth:`.Tesliper.export_spectra` will export spectra calculated so far that are stored
  in :attr:`.Tesliper.spectra` dictionary;
- :meth:`.Tesliper.export_averaged` will export each averaged spectrum calculated so far
  that is stored in :attr:`.Tesliper.averaged` dictionary.

Each of these methods take two parameters: *fmt* and *mode*. *fmt* is a file format, to
which data should be exported. It should be one of the following values: ``"txt"`` (the
default), ``"csv"``, ``"xlsx"``. *mode* denotes how files should be opened and should be
one of: ``"a"`` (append to existing file), ``"x"`` (the default, only write if file
doesn't exist yet), ``"w"`` (overwrite file if it already exists).

These export methods will usually produce a number of files in the
:attr:`.Tesliper.output_dir`, which names will be picked automatically, according to the
genre of the exported data and/or the conformer that data relates to.
:meth:`~.Tesliper.export_energies` will produce a file for each available energies genre
and an additional overview file, :meth:`~.Tesliper.export_spectra` will create a file
for spectra genre and each conformer, and so on.

.. note::
    ``"xlsx"`` format is an exception from the above - it will produce only one file,
    named "tesliper-output.xlsx", and create multiple spreadsheets inside this file. The
    appending mode is useful when exporting data to ``"xlsx"`` format, as it allows to
    write multiple kinds of data (with calls to multiple of these methods) to this
    single destination .xlsx file.

There is also a :meth:`.Tesliper.export_data` available, which will export only genres
you specifically request (plus "freq" or "wavelen", if any genre given genre is a
spectral data-related). The same applies here in the context of output format, write
mode, and names of produced files.

.. tip::

    If you would like to customize names of the files produced, you will need to
    directly use one of the writer objects provided by ``tesliper``. Refer to the
    :mod:`.writing` module documentation for more information.

    .. TODO: link to advanced guide when its done

Creating input files
''''''''''''''''''''

You can use :meth:`.Tesliper.export_job_file` to prepare input files for the quantum
chemical calculations software. Apart from the typical *fmt* (only ``"gjf"`` is
supported by default) and *mode* parameters, this method also accepts the
*geometry_genre* and any number of additional keyword parameters, specifying
calculations details. *geometry_genre* should be a name of the data genre, representing
conformers' geometry, that should be used as input geometry. Additional keyword
parameters are passed to the writer object, relevant to the *fmt* requested. Keywords
supported by the :class:`default "gjf"-format writer <.GjfWriter>` are as follows:

    route
        A calculations route: keywords specifying calculations directives for quantum
        chemical calculations software.
    link0
        Dictionary with "link zero" commands, where each key is command's name and each
        value is this command's parameter.
    comment
        Contents of title section, i.e. a comment about the calculations.
    post_spec
        Anything that should be placed after conformer's geometry specification.
        Will be written to the file as given.

``link0`` parameter should be explained in more details. It supports standard link zero
commands used with Gaussian software, like ``Mem``, ``Chk``, or ``NoSave``. Full list of
these commands may be found in the documentation for :attr:`.GjfWriter.link0`. Any
non-parametric ``link0`` command (i.e. ``Save``, ``NoSave``, and ``ErrorSave``),  should
be given a ``True`` value if it should be included in the ``link0`` section.

Path-like commands, e.g. ``Chk`` or ``RWF``, may be parametrized for each conformer. You
can put a placeholder inside a given string path, which will be substituted when writing
to file. The most useful placeholders are probably ``${conf}`` and ``${num}`` that
evaluate to conformer's name and ordinal number respectively. More information about
placeholders may be found in :meth:`.GjfWriter.make_name` documentation.

.. code-block:: python

    >>> list(tslr.conformers.kept_keys())
    ["conf_one", "conf_three"]
    >>> tslr.export_job_file(
    ...     geometry_genre="optimized_geom",
    ...     route="# td=(nstates=80)",
    ...     comment="Example of parametrization in .gjf files",
    ...     link0={
    ...         "Mem": "10MW",
    ...         "Chk": "path/to/${conf}.chk",
    ...         "Save": True,
    ...     },
    ... )
    >>> [file.name for file in tslr.output_dir.iterdir()]
    ["conf_one.gjf", "conf_three.gjf"]

Then contents of "conf_one.gjf" is:

.. code-block:: none

    %Mem=10MW
    %Chk=path/to/conf_one.gjf
    %Save

    # td=(nstates=80)

    Example of parametrization in .gjf files

    [geometry specification...]

Saving session for later
------------------------

If you'd like to come back to the data currently contained within a ``tesliper``
instance, you may serialize it using :meth:`.Tesliper.serialize` method. Provide the
method call with a *filename* parameter, under which filename the session should be
stored inside the current :attr:`~.Tesliper.output_dir`. You may also omit it to use the
default ``".tslr"`` name. All data, extracted and calculated, including current
:term:`kept` status of each conformer, is saved and may be loaded later using
:meth:`.Tesliper.load` class method.

.. code-block:: python

    curr_dir = tslr.output_dir
    tslr.serialize()
    loaded = Tesliper.load(curr_dir / ".tslr")
    assert loaded.conformers == tslr.conformers
    assert loaded.conformers.kept == tslr.conformers.kept
    assert loaded.spectra.keys() == tslr.spectra.keys()
.. _available genres:

Available data genres
=====================

:freq: *list of floats*; available from freq job
    
    harmonic vibrational frequencies (cm^-1)

:mass: *list of floats*; available from freq job
    
    reduced masses (AMU)

:frc: *list of floats*; available from freq job
    
    force constants (mDyne/A)

:iri: *list of floats*; available from freq job
    
    IR intensities (KM/mole)

:dip: *list of floats*; available from freq=VCD job
    
    dipole strengths (10**-40 esu**2-cm**2)

:rot: *list of floats*; available from freq=VCD job
    
    rotational strengths (10**-44 esu**2-cm**2)

:emang: *list of floats*; available from freq=VCD job
    
    E-M angle = Angle between electric and magnetic dipole transition moments (deg)

:depolarp: *list of floats*; available from freq=Raman job
    
    depolarization ratios for plane incident light

:depolaru: *list of floats*; available from freq=Raman job
    
    depolarization ratios for unpolarized incident light

:ramanactiv: *list of floats*; available from freq=Raman job
    
    Raman scattering activities (A**4/AMU)

:ramact: *list of floats*; available from freq=ROA job
    
    Raman scattering activities (A**4/AMU)

:depp: *list of floats*; available from freq=ROA job
    
    depolarization ratios for plane incident light

:depu: *list of floats*; available from freq=ROA job
    
    depolarization ratios for unpolarized incident light

:alpha2: *list of floats*; available from freq=ROA job
    
    Raman invariants Alpha2 = alpha**2 (A**4/AMU)

:beta2: *list of floats*; available from freq=ROA job
    
    Raman invariants Beta2 = beta(alpha)**2 (A**4/AMU)

:alphag: *list of floats*; available from freq=ROA job
    
    ROA invariants AlphaG = alphaG'(10**4 A**5/AMU)

:gamma2: *list of floats*; available from freq=ROA job
    
    ROA invariants Gamma2 = beta(G')**2 (10**4 A**5/AMU)

:delta2: *list of floats*; available from freq=ROA job
    
    ROA invariants Delta2 = beta(A)**2, (10**4 A**5/AMU)

:raman1: *list of floats*; available from freq=ROA job
    
    Far-From-Resonance Raman intensities =ICPu/SCPu(180) (K)

:roa1: *list of floats*; available from freq=ROA job
    
    ROA intensities =ICPu/SCPu(180) (10**4 K)

:cid1: *list of floats*; available from freq=ROA job
    
    CID=(ROA/Raman)*10**4 =ICPu/SCPu(180)

:raman2: *list of floats*; available from freq=ROA job
    
    Far-From-Resonance Raman intensities =ICPd/SCPd(90) (K)

:roa2: *list of floats*; available from freq=ROA job
    
    ROA intensities =ICPd/SCPd(90) (10**4 K)

:cid2: *list of floats*; available from freq=ROA job
    
    CID=(ROA/Raman)*10**4 =ICPd/SCPd(90)

:raman3: *list of floats*; available from freq=ROA job
    
    Far-From-Resonance Raman intensities =DCPI(180) (K)

:roa3: *list of floats*; available from freq=ROA job
    
    ROA intensities =DCPI(180) (10**4 K)

:cid3: *list of floats*; available from freq=ROA job
    
    CID=(ROA/Raman)*10**4 =DCPI(180)

:rc180: *list of floats*; available from freq=ROA job
    
    RC180 = degree of circularity

:wavelen: *list of floats*; available from td job
    
    excitation energies (nm)

:ex_en: *list of floats*; available from td job
    
    excitation energies (eV)

:eemang: *list of floats*; available from td job
    
    E-M angle = Angle between electric and magnetic dipole transition moments (deg)

:vdip: *list of floats*; available from td job
    
    dipole strengths (velocity)

:ldip: *list of floats*; available from td job
    
    dipole strengths (length)

:vrot: *list of floats*; available from td job
    
    rotatory strengths (velocity) in cgs (10**-40 erg-esu-cm/Gauss)

:lrot: *list of floats*; available from td job
    
    rotatory strengths (length) in cgs (10**-40 erg-esu-cm/Gauss)

:vosc: *list of floats*; available from td job
    
    oscillator strengths

:losc: *list of floats*; available from td job
    
    oscillator strengths

:transitions: *list of lists of lists of (int, int, float)*; available from td job
    
    transitions (first to second) and their coefficients (third)

:scf: *float*; always available
    
    SCF energy

:zpe: *float*; available from freq job
    
    Sum of electronic and zero-point Energies (Hartree/Particle)

:ten: *float*; available from freq job
    
    Sum of electronic and thermal Energies (Hartree/Particle)

:ent: *float*; available from freq job
    
    Sum of electronic and thermal Enthalpies (Hartree/Particle)

:gib: *float*; available from freq job
    
    Sum of electronic and thermal Free Energies (Hartree/Particle)

:zpecorr: *float*; available from freq job
    
    Zero-point correction (Hartree/Particle)

:tencorr: *float*; available from freq job
    
    Thermal correction to Energy (Hartree/Particle)

:entcorr: *float*; available from freq job
    
    Thermal correction to Enthalpy (Hartree/Particle)

:gibcorr: *float*; available from freq job
    
    Thermal correction to Gibbs Free Energy (Hartree/Particle)

:command: *str*; always available
    
    command used for calculations

:normal_termination: *bool*; always available
    
    true if Gaussian job seem to exit normally, false otherwise

:optimization_completed: *bool*; available from opt job
    
    true if structure optimization was performed successfully

:version: *str*; always available
    
    version of Gaussian software used

:charge: *int*; always available
    
    molecule's charge

:multiplicity: *int*; always available
    
    molecule's spin multiplicity

:input_atoms: *list of str*; always available
    
    input atoms as a list of atoms' symbols

:input_geom: *list of lists of floats*; always available
    
    input geometry as X, Y, Z coordinates of atoms

:stoichiometry: *str*; always available
    
    molecule's stoichiometry

:last_read_atoms: *list of ints*; always available
    
    molecule's atoms as atomic numbers

:last_read_geom: *list of lists of floats*; always available
    
    molecule's geometry (last one found in file) as X, Y, Z coordinates of atoms

:optimized_atoms: *list of ints*; available from successful opt job
    
    molecule's atoms read from optimized geometry as atomic numbers

:optimized_geom: *list of lists of floats*; available from successful opt job
    
    optimized geometry as X, Y, Z coordinates of atoms

.. image:: /_static/banner.png

|

.. just a vertical space to separate banner and header

Welcome to tesliper's documentation!
====================================

``tesliper`` is a package for batch processing of Gaussian output files, focusing on
extraction and processing of data related to simulation of optical spectra. The software
offers a Python API and a graphical user interface (GUI), allowing for your preferred
style of interaction with the computer: visual or textual. It's main goal is to minimize
time and manual work needed to simulate optical spectrum of investigated compound.

Key features
------------

``tesliper`` was designed for working with multiple conformers of a compound,
represented by a number of files obtained from Gaussian quantum-chemical computations
software. It allows you easily exclude conformers that are not suitable for further
analysis: erroneous, not optimized, of higher energy than a user-given threshold, or
very similar to some other structure in the set. Data parsed from files and data
calculated may be exported to other file formats for storage or further analysis with
other tools. Below is a quick overview of features it provides:

- Batch processing of Gaussian output files regarding structure optimization and
  simulation of spectral properties
- Conditional, property-based filtering of conformers
- Geometry comparison via the RMSD sieve
- Calculation of Boltzmann distribution—based populations of conformers
- Simulation of IR, VCD, UV, ECD, Raman, and ROA spectra from spectral activities
- Export of extracted and calculated data to .txt, .csv, and .xlsx file formats
- Export of .gjf files for further calculations in Gaussian software
- Free & open source (OSI approved BSD 2-Clause license)
- Graphical and programmatic interfaces

Motivation and context
----------------------

Simulation of optical spectra of organic compounds becomes one of the routine tasks
for chemical analysts -- it is a necessary step in one of the increasingly popular
methods of establishing compound's absolute configuration. However, the process
of obtaining a simulated spectrum may be cumbersome, as it usually involves analyzing
a large number of potentially stable conformers of the studied molecule.
``tesliper`` was created to aid in such work.

It should be noted that ``tesliper`` is not the only software that is capable of
providing a simulated spectrum, given output of quantum-chemical computations. The table
below summarizes other available GUI tools and compares features they offer.
Among listed ``tesliper`` is the only one that is open source
and allows to easily filter parsed data.

.. list-table:: How does ``tesliper`` fit into the market?
   :header-rows: 1
   
   *  -
      - Tesliper
      - SpecDis [1]_
      - CDspecTech [2]_
      - ComputeVOA [3]_
      - GaussView [4]_
      - ChemCraft [5]_
   *  - Free
      - |check|        
      - |check|          
      - |check|           
      - |cross|          
      - |cross|           
      - |cross|
   *  - Open Source
      - |check|        
      - |cross|          
      - |cross|           
      - |cross|          
      - |cross|           
      - |cross|
   *  - Batch Processing
      - |check|        
      - |check|          
      - |check|           
      - |cross|     
      - .gjf export
      - .gjf modif.
   *  - Geometry Comparison
      - |check|        
      - |cross|          
      - |cross|           
      - |check|          
      - |cross|           
      - |check|
   *  - Averaging
      - |check|        
      - |check|          
      - |check|           
      - |cross|          
      - |cross|           
      - |cross|
   *  - Conditional Filtering
      - |check|        
      - |cross|          
      - |cross|           
      - |cross|          
      - |cross|           
      - |cross|
   *  - Job File Creation
      - |check|        
      - |cross|          
      - |cross|           
      - |check|          
      - |check|           
      - |check|
   *  - Electronic Spectra
      - |check|        
      - |check|          
      - |check|           
      - |cross|          
      - |check|           
      - |cross|
   *  - Scattering Spectra
      - |check|        
      - |cross|          
      - |check|           
      - |check|          
      - |check|           
      - |check|
   *  - Multi-platform
      - |check|        
      - |check|          
      - |check|           
      - |cross|          
      - |check|       
      - needs wine
   *  - Conformational Search
      - |cross|        
      - |cross|          
      - |cross|           
      - |check|       
      - optional
      - |cross|
   *  - Molecule Visualization
      - |cross|       
      - |cross|          
      - |cross|           
      - |check|          
      - |check|           
      - |check|   


References
----------

.. [1] **SpecDis**: T. Bruhn, A. Schaumlöffel, Y. Hemberger, G. Pescitelli,
   *SpecDis version 1.71*, Berlin, Germany, **2017**, http://specdis-software.jimdo.com
.. [2] **CDspecTech**: C. L. Covington, P. L. Polavarapu, Chirality, **2017**, 29, 5,
   p. 178, DOI: `10.1002/chir.22691 <https://doi.org/10.1002/chir.22691>`_
.. [3] **ComputeVOA**: E. Debie, P. Bultinck, L. A. Nafie, R. K. Dukor, BioTools Inc.,
   Jupiter, FL, **2010**, https://biotools.us/software-2
.. [4] **GaussView**: R. Dennington, T. A. Keith, J. M. Millam, Semichem Inc.,
   *GaussView version 6.1*, Shawnee Mission, KS, **2016**
.. [5] **ChemCraft**: https://www.chemcraftprog.com


.. toctree::
   :hidden:
   :caption: Introduction

   Home page <self>
   Installation <installation>
   Conventions and Terms <glossary>

.. toctree::
   :hidden:
   :caption: Tutorials

   Graphical interface <gui>
   Scripting basics <tutorial>
   Advanced guide <advanced>

.. toctree::
   :hidden:
   :caption: Reference

   Available data genres <genres>
   Math and Algorithms <math>
   API reference <_autosummary/tesliper>

.. toctree::
   :hidden:
   :caption: Appendix

   Change Log <changelog>
   Index <genindex>


.. |check| unicode:: U+2714
.. |cross| unicode:: U+2718
Math and Algorithms
===================

.. |ith| replace:: :math:`i`:superscript:`th`

Simulation of spectra
---------------------

To simulate a spectrum, each band's theoretical signal intensity, derived from quantum
chemical calculations of corresponding optical activity, must be expressed as a
broadened peak, instead of the single scalar value. To simulate peak's shape one of the
curve fitting functions is used. ``tesliper`` implements two, most commonly used, such
functions: gaussian function\ [#gaussian]_ and lorentzian function\ [#lorentzian]_.

For each point on the simulated spectrum's abscissa, the corresponding signal intensity
is calculated by applying the fitting function to all bands of the conformer and summing
resulting values.

.. [#gaussian] https://mathworld.wolfram.com/GaussianFunction.html
.. [#lorentzian] https://mathworld.wolfram.com/LorentzianFunction.html

Gaussian fitting function
'''''''''''''''''''''''''

.. math::

    f(\nu) = \frac{1}{\sigma\sqrt{2\pi}}\sum\limits_i I_i e^{
        -(\nu_i - \nu)^2 / (2\sigma^2)
    }

:math:`\nu`
    Arbitrary point on the :math:`x`-axis, for which the signal intensity is calculated.
:math:`\nu_i`
    Point on the :math:`x`-axis, at which the |ith| band occur.
:math:`I_i`
    Intensity of the |ith| band.
:math:`\sigma = \sqrt{2}\omega`
    Standard derivation, in this context interpreted as equal to :math:`\sqrt{2}`
    times :math:`\omega`.
:math:`\omega`
    Half width of the peak at :math:`\frac{1}{e}` of its maximum value (HW1OeM),
    expressed in the :math:`x`-axis units.

Lorentzian fitting function
'''''''''''''''''''''''''''

.. math::

    f(\nu) = \frac{\gamma}{\pi}\sum\limits_i\frac{I_i}{(\nu_i - \nu)^2 + \gamma^2}

:math:`\nu`
    Arbitrary point on the :math:`x`-axis, for which the signal intensity is calculated.
:math:`\nu_i`
    Point on the :math:`x`-axis, at which the |ith| band occur.
:math:`I_i`
    Intensity of the |ith| band.
:math:`\gamma`
    Half width of the peak at half of its maximum value (HWHM),
    expressed in the :math:`x`-axis units.

Calculation of intensities
--------------------------

Dipole strength and rotator strength is converted to the theoretical intensity as
described by Polavarapu\ [#polavarapu]_. Constants used in below equations are as follows.

:math:`c = 2.99792458 \times 10^{10}\ \mathrm{cm}\cdot\mathrm{s}^{-1}`
    Speed of light.
:math:`h = 6.62606896 \times 10^{-30}\ \mathrm{kg}\cdot\mathrm{cm}^2\cdot\mathrm{s}^{-1}`
    Planck's constant.
:math:`N_A = 6.02214199 \times 10^{23}\ \mathrm{mol}^{-1}`
    Avogadro's constant.
:math:`m_e = 9.10938 \times 10^{-28}\ \mathrm{g}`
    Mass of the electron.
:math:`e = 4.803204 \times 10^{-10}\ \mathrm{esu}`
    The charge on the electron.

.. [#polavarapu] Prasad L. Polavarapu (2017), *Chiroptical Spectroscopy Fundamentals and
    Applications*, CRC Press


Dipole strength to IR intensities
'''''''''''''''''''''''''''''''''

.. math::

    I_k = \frac{100 \cdot 8 \pi^3 N_A}{3 \cdot \ln(10) \cdot hc} D_k \nu_k
    = 0.010886 \cdot D_k \nu_k

:math:`D_k \ [\times 10^{-40}\ \mathrm{esu}^2\mathrm{cm}^2]`
    Dipole strength of :math:`k^{\mathrm{th}}` transition.
:math:`\nu_k \ [\mathrm{cm}^{-1}]`
    Frequency of the :math:`k^{\mathrm{th}}` transition.
:math:`I_k \ [\mathrm{L}\,\mathrm{mol}^{-1}\mathrm{cm}^{-1}]`
    Theoretical ("zero-width") peak intensity in terms of the decadic molar absorption
    coefficient :math:`\epsilon`.

Rotator strength to VCD intensities
'''''''''''''''''''''''''''''''''''

.. math::

    I_k = 4 \cdot \frac{8 \pi^3 N_A}{3 \cdot \ln(10) \cdot hc} R_k \nu_k
    = 0.0435441 \cdot R_k \nu_k

:math:`R_k \ [\times 10^{-44}\ \mathrm{esu}^2\mathrm{cm}^2]`
    Rotator strength of :math:`k^{\mathrm{th}}` transition.
:math:`\nu_k \ [\mathrm{cm}^{-1}]`
    Frequency of the :math:`k^{\mathrm{th}}` transition.
:math:`I_k \ [\mathrm{L}\,\mathrm{mol}^{-1}\mathrm{cm}^{-1}]`
    Theoretical ("zero-width") peak intensity in terms of the decadic molar absorption
    coefficient :math:`\epsilon`.

Dipole strength to UV intensities
'''''''''''''''''''''''''''''''''

.. math::

    I_k = \frac{8 \pi^3 N_A}{3 \cdot \ln(10) \cdot hc} D_k \lambda_k
    = 0.010886 \cdot D_k \lambda_k

:math:`D_k \ [\times 10^{-40}\ \mathrm{esu}^2\mathrm{cm}^2]`
    Dipole strength of :math:`k^{\mathrm{th}}` transition.
:math:`\lambda_k \ [\mathrm{cm}^{-1}]`
    Wavelength of the :math:`k^{\mathrm{th}}` transition.
:math:`I_k \ [\mathrm{L}\,\mathrm{mol}^{-1}\mathrm{cm}^{-1}]`
    Theoretical ("zero-width") peak intensity in terms of the decadic molar absorption
    coefficient :math:`\epsilon`.

Rotator strength to ECD intensities
'''''''''''''''''''''''''''''''''''

.. math::

    I_k = 4 \cdot \frac{8 \pi^3 N_A}{3 \cdot \ln(10) \cdot hc} R_k \lambda_k
    = 0.0435441 \cdot R_k \lambda_k

:math:`R_k \ [\times 10^{-44}\ \mathrm{esu}^2\mathrm{cm}^2]`
    Rotator strength of :math:`k^{\mathrm{th}}` transition.
:math:`\lambda_k \ [\mathrm{cm}^{-1}]`
    Wavelength of the :math:`k^{\mathrm{th}}` transition.
:math:`I_k \ [\mathrm{L}\,\mathrm{mol}^{-1}\mathrm{cm}^{-1}]`
    Theoretical ("zero-width") peak intensity in terms of the decadic molar absorption
    coefficient :math:`\epsilon`.

Oscillator strength to UV intensities
'''''''''''''''''''''''''''''''''''''

Conversion from oscillator strength to signal intensity of UV spectrum is calculated
as described by Gaussian\ [#gaussian_uv]_. 

.. math::

    I_k = \frac{e^2 N_A}{10^3 \ln(10) m_e c^2} f_k
    = 2.315351857 \times 10^8 f_k

:math:`f_k`
    Oscillator strength of :math:`k^{\mathrm{th}}` transition.
:math:`I_k \ [\mathrm{L}\,\mathrm{mol}^{-1}\mathrm{cm}^{-1}]`
    Theoretical ("zero-width") peak intensity in terms of the decadic molar absorption
    coefficient :math:`\epsilon`.

.. [#gaussian_uv] https://gaussian.com/uvvisplot/

Raman/ROA intensities
'''''''''''''''''''''

Gaussian-provided Raman and ROA activities are used without any conversion.

Population of conformers
------------------------

Population of conformers is calculated according to the Boltzmann probability
distribution that "gives the probability that a system will be in a certain state as a
function of that state's energy and the temperature of the system."\ [#boltzmann]_ In
this context each conformer is considered one of the possible states of the system (a
studied molecule).

Firstly, we calculate a Boltzmann factors for each conformer in respect to the most
stable conformer (the one of the lowest energy). Boltzmann factor of two states is
defined as:

.. math::

    B^a_b = \frac{F(state_a)}{F(state_b)} = e^{(E_b - E_a)/kt}

where:

:math:`E_a` and :math:`E_b`
    energies of states :math:`a` and :math:`b`;
:math:`k = 0.0019872041 \: \mathrm{kcal/(mol*K)}`
    Boltzmann constant;
:math:`t` 
    temperature of the system.

Boltzmann factor represents a ratio of probabilities of the two states being occupied.
In other words, it shows how much more likely it is for the molecule to take the form of
one conformer over another conformer. Having a ratio of these probabilities for each
possible conformer in respect to the most stable conformer, we are able to find the
distribution of conformers (probability of taking the form of each conformer):

.. math::

    p_i = \frac{B_0^i}{\sum\limits_j^{states}B_0^j}

assuming that :math:`state_0` is the state of the lowest energy (the most stable
conformer).

.. [#boltzmann] https://en.wikipedia.org/wiki/Boltzmann_distribution

RMSD of conformers
------------------

RMSD, or root-mean-square deviation of atomic positions, is used as a measure of
similarity between two conformers. As its name hints, it is an average distance between
atoms in the two studied conformers: the lower the RMSD value, the more similar are
conformers in question.

Finding minimized value of RMSD
'''''''''''''''''''''''''''''''

In a typical output of the quantum chemical calculations software, molecule is
represented by a number of points (mapping to particular atoms) in a 3-dimensional
space. Usually, orientation and position of the molecule in the coordinate system is
arbitrary and simple overlay of the two conformers may not be the same as their optimal
overlap. To neglect the effect of conformers' rotation and shift on the similarity
measure, we will look for the common reference frame and optimal alignment of atoms.

Zero-centring atomic coordinates
""""""""""""""""""""""""""""""""

To find the common reference frame for two conformers we move both to the origin of the
coordinate system. This is done by calculating a centroid of a conformer and subtracting
it from each point representing an atom. The centroid is given as an arithmetic mean
of all the atoms in the conformer:

.. math::

    a^0_i = a_i - \frac{1}{n}\sum\limits_{j=a}^{n}a_j

where:

:math:`a_i`, :math:`a_j`
    atom's original position in the coordinate system;
:math:`a^0_i`
    atom's centered position in the coordinate system;
:math:`n`
    number of atoms in the molecule.

Rotating with Kabsch algorithm
""""""""""""""""""""""""""""""

Optimal rotation of one conformer onto another is achieved using a Kabsch algorithm\
[#kabsch]_ (also known as Wahba's problem\ [#wanba]_). Interpreting positions of each
conformers' atoms as a matrix, we find the covariance matrix :math:`H` of these
matrices (:math:`P` and :math:`Q`):

.. math::

    H = P^\intercal Q

and then we use the singular value decomposition (SVD)\ [#svd]_ routine to get :math:`U`
and :math:`V` unitary matrices.

.. math::

    H = U \Sigma V ^\intercal

Having these, we can calculate the optimal rotation matrix as:

.. math::

    R = V \begin{pmatrix}1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & d\end{pmatrix} U ^\intercal

where :math:`d = \mathrm{sign}(\mathrm{det}(VU^\intercal))` that allows to ensure a
right-handed coordinate system.

.. note::

    To allow for calculation of th best rotation between sets of molecules and to
    compromise between efficiency and simplicity of implementation, ``tesliper`` uses
    Einstein summation convention\ [#einsum]_ *via* :func:`numpy.einsum` function. The
    implementation is as follows:

    .. literalinclude:: ../../tesliper/datawork/geometry.py
        :language: python
        :pyobject: kabsch_rotate
        :lines: 1-11,23,25-
        :emphasize-lines: 19,32

.. [#kabsch] https://en.wikipedia.org/wiki/Kabsch_algorithm
.. [#wanba] https://en.wikipedia.org/wiki/Wahba%27s_problem
.. [#svd] https://en.wikipedia.org/wiki/Singular_value_decomposition
.. [#einsum] https://en.wikipedia.org/wiki/Einstein_notation

Calculating RMSD of atomic positions
""""""""""""""""""""""""""""""""""""

Once conformers are aligned, the value of RMSD\ [#rmsd]_ is calculated simply by finding
a distance between each equivalent atoms and averaging their squares and finding the
root of this average:

.. math::

    \mathrm{RMSD} = \sqrt{\frac{1}{n}\sum\limits_i^n(p_i - q_i)^2}

where:

:math:`p_i` and :math:`q_i`
    positions of |ith| equivalent atoms in conformers :math:`P` and :math:`Q`;
:math:`n`
    number of atoms in each conformer.

.. [#rmsd] https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions

Comparing conformers
''''''''''''''''''''

To compare conformers as efficiently as possible, the RMSD values are calculated not
in the each-to-each scheme, but inside a rather small moving window. The size of this
window determines how many calculations will be done for the whole collection.

Moving window mechanism
"""""""""""""""""""""""

``tesliper`` provides three types of moving windows: a :func:`fixed
<.geometry.fixed_windows>`, :func:`stretching <.geometry.stretching_windows>`, and
:func:`pyramid <.geometry.pyramid_windows>` windows. The strategy you choose will affect
both the performance and the accuracy of the RMSD sieve, as described below.

:func:`fixed <.geometry.fixed_windows>`
    The most basic sliding window of a fixed size. Provides the most control over the
    performance of the sieve, but is the least accurate.
:func:`stretching <.geometry.stretching_windows>`
    The default, allows to specify the size of the window in the context of some numeric
    property, usually the energy of conformers. The size may differ in the sense of the
    number of conformers in each window, but the difference between maximum and minimum
    values of said property inside a window will not be bigger than the given *size*.
    Provides a best compromise between the performance and the accuracy.
:func:`pyramid <.geometry.pyramid_windows>`
    The first window will contain the whole collection and each consecutive window will
    be smaller by one conformer. Allows to perform a each-to-each comparison, but in
    logarithmic time rather than quadratic time. Best accuracy but worst performance.

.. note::

    The actual windows produced by sliding window functions are iterables of
    :class:`numpy.ndarray`\s of indices (that point to the value in the original
    array of conformers).

The sieve
"""""""""

The :func:`RMSD sieve function <.geometry.rmsd_sieve>` takes care of zero-centring and
finding the best overlap of the conformers, as described previously. Aside form this, it
works as follows: for each window, provided by one of the moving window functions
described above, it takes the first conformer in the window (reference) and calculates
it's minimum RMSD value with respect to each of the other conformers in this window. Any
conformers that have lower RMSD value than a given threshold, will be considered
identical to the reference conformer and internally marked as :term:`not kept <kept>`.
The sieve returns an array of boolean value for each conformer: ``True`` if conformer's
structure is "original" and should be kept, ``False`` if it is a duplicate of other,
"original" structure (at least according to threshold given), and should be discarded.

Spectra transformation
----------------------

Finding best shift
''''''''''''''''''

Optimal offset of two spectra is determined by calculating their cross-correlation\
[#cross-corr]_ (understood as in the signal processing context) and finding its maximum
value. Index of this max value of the discrete cross-correlation array indicates the
position of one spectrum in respect to the other spectrum, in which the overlap of the
two is the greatest.

.. [#cross-corr] https://en.wikipedia.org/wiki/Cross-correlation

Finding optimal scaling
'''''''''''''''''''''''

Optimal scaling factor of spectra is determined by comparing a mean *y* values of target
spectrum and a reference spectrum. Values lower than 1% of maximum absolute *y* value of
each spectrum are ignored.

Installation
============

.. _install-gui:

GUI for Windows
---------------

For Windows users that only intend to use graphical interface, ``tesliper`` is available
as a standalone .exe application, available for download from the `latest release
<https://github.com/mishioo/tesliper/releases/latest/>`_ under the **Assets** section at
the bottom of the page. No installation is required, just double-click the downloaded
**Tesliper.exe** file to run the application.

Unfortunately, a single-file installation is not available for unix-like systems.
Please follow a terminal-based installation instructions below.

Install from terminal
---------------------

``tesliper`` is a Python package `distributed via PyPI <https://pypi.org/project/tesliper/>`_.
You can install it to your python distribution simply by running::

    python -m pip install tesliper

in your terminal. This will download and install ``tesliper`` along with it's essential
dependencies. A graphical interface have an additional dependency, but it may be
easily included in your installation if you use ``python -m pip install tesliper[gui]``
instead. Some users of unix-like systems may also need to instal ``tkinter`` manually,
if it is not included in their distribution by default. Please refer to relevant online
resources on how to do this in your system, if that is your case.

.. note::
    Reminder for zsh users to quote the braces like this ``'tesliper[gui]'``
    or like this ``tesliper\[gui]`` when installing extras. This is necessary, because
    normally zsh uses ``some[thing]`` syntax for pattern matching.

Requirements
------------

This software needs at least Python 3.6 to run. It also uses some additional packages::

    numpy
    openpyxl
    matplotlib (optional, for GUI)

.. note::
    ``tesliper`` uses ``tkinter`` to deliver the graphical interface. It is included in
    most Python distributions, but please be aware, that some might miss it. You will
    need to install it manually in such case... autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

   tesliper

.. This file is a placeholder and will be replaced

Index
=====Graphical Interface
===================

This part discusses the use of the Graphical User Interface (GUI). For tutorial on using
``tesliper`` in Python scripts, see :ref:`tutorial`.

On Windows system you may start the GUI by downloading and double-clicking
``Tesliper.exe`` file available in the `latest release
<https://github.com/mishioo/tesliper/releases/latest/>`_, as described in the
:ref:`Installation section <install-gui>`. Executable files are not available for other
systems, unfortunately, but you may start the GUI from the command line as well:

.. code-block:: bash

    $ python -m pip install tesliper[gui]  # only once
    $ tesliper-gui  # starts GUI

.. note::

    If you'd like to start the graphical interface from the local copy, you may also
    run it as a module with ``python -m tesliper.gui``.

Please note that the first launch may take additional time. After the application
starts, a window like the one bellow will appear. It's actual looks will depend on your
operating system. 

.. image:: _static/screenshots/16426971762.png

The Interface is divided in two parts: controls on the left and views (initially empty)
on the right. Controls panel are further divided into sections, some of which may be
collapsed by clicking on the section title (those with a small arrow on the left).
Each section will be described further in this tutorial, in the appropriate section.

There are tree views available: ``Extracted data`` and ``Energies list`` summarize all
conformers read from files. ``Extracted data`` details what data is available and shows
status of calculations for each conformer. ``Energies list`` shows values of conformers'
energies calculated by quantum chemical software and derived values: Boltzmann factors
and conformers' populations.

Reading files
-------------

``tesliper`` supports reading data from computations performed using Gaussian software.
To load data, use controls in the ``Extract data`` section. ``Choose files`` button
allows you to select individual files to read using the popup dialog. ``Choose folder``
button shows a similar dialog, but allowing you to select a single directory - all
Gaussian output files in this directory (but not subdirectories) will be read. Finally,
using the ``Recursive`` button will also read files form all subdirectories,
recursively.

.. image:: _static/screenshots/16426972102.png
.. image:: _static/screenshots/16426972279.png
.. image:: _static/screenshots/16426972380.png

.. note::

    Make sure you do not have mixed .log and .out files in the directory, when using
    ``Choose folder`` button.

Once you select files or directory and confirm your selection, the process of data
extraction will start and ``Extract data`` view will be updated for each read conformer.
It will show if calculation job terminated normally (``Termination``), if conformer's
structure was optimized and if optimization was successful (``Opt``), if extended set of
energies is available (``Energy``), what spectral data it available (``IR``, ``VCD``,
``UV``, ``ECD``, ``Raman``, ``ROA``), how many imaginary frequencies are reported for
conformer (``Imag``), and what is conformer's stoichiometry (``Stoichiometry``).

.. image:: _static/screenshots/16426972737.png

When data extraction is finished, the status barr att the bottom will show ``Idle``
again. After reading first portion of files, you may tick the ``Ignore unknown
conformers`` option. When this option is ticked, ``tesliper`` will only read files that
correspond to conformers it already knows (judging by the filename).

Trimming conformers
-------------------

Conformers may be marked as :term:`kept` :term:`not kept <kept>` (:term:`trimmed
<trimming>`). Only :term:`kept` conformers are processed by ``tesliper``,
:term:`trimmed <trimming>` ones are ignored. This mechanism allows you to select which
conformers should be included in the final averaged spectrum, etc. Trimmed conformers
are shown in gray.

.. image:: _static/screenshots/16426974012.png

``Kept conformers`` section shows how many conformers contain certain data and allows to
easily keep/trim whole groups of conformers, using ``keep`` and ``trim`` buttons beside
the appropriate group. You may also keep/trim individual conformers by ticking/unticking
checkboxes beside the conformers name (left of ``Filenames`` column).

.. image:: _static/screenshots/16426974044.png

After finished data extraction and after each manual trimming, auto-trimming is
performed to make sure corrupted or invalid conformers are not accidentally
:term:`kept`. Checkboxes in the ``Auto-trim`` subsection, shown below, control which
conformers should be always trimmed.

.. image:: _static/screenshots/16426974255.png

.. tip::

    ``Incomplete entries`` are conformers that miss some data, which other conformers
    include, e.g. those that were left out in one of calculations steps. ``Inconsistent
    data sizes`` indicates that some multi-value data has different number of data
    points than in case of other conformers. This usually suggests that conformer in
    question is not actually a conformer but a different molecule.

Trimming with sieves
--------------------

The ``Energies and structure`` section, described in this part, is related with the
``Energies list`` view. This view shows, as the name suggests, list of energies for each
conformer and energies-derived values.

.. image:: _static/screenshots/16426974610.png

Using a ``Show:`` drop-down menu you may select a different energies-derived data to
show in the view. ``Delta`` is conformer's energy difference to the most stable
(lowest-energy) conformer (in :math:`\mathrm{kcal}/\mathrm{mol}` units), ``Min.
Boltzmann factor`` is conformer's Boltzmann factor in respect to the most stable
conformer (unitless) and ``Popuation`` is population of conformers according to the
Boltzmann distribution (in perecnt). Original ``Energy`` values are shown in Hartree
units.

.. image:: _static/screenshots/16426975018.png

Both types of sieves provided depend on the selected value of the ``Use:`` drop-down
menu. It determines, which energy values are used by the sieves. Only available energies
wil be shown in the list. In case their names are not intuitive enough, here is the
explanation:

| ``Thermal``: sum of electronic and thermal Energies;
| ``Enthalpy``: sum of electronic and thermal Enthalpies;
| ``Gibbs``: sum of electronic and thermal Free Energies;
| ``SCF``: energy calculated with the self-consistent field method;
| ``Zero-Point``: sum of electronic and zero-point Energies.

.. image:: _static/screenshots/16426976160.png

The ``Range sieve`` lets you to trim conformers that have a current ``Show:`` value
outside of the specified range. After you fill the ``Minimum`` and ``Maximum`` fields to
match your needs, click ``Trim to...`` button to perform trimming. The example below
shows trimming of conformers, which Free Energy-derived population is below 1%. Please
note that valuesin the ``Energies list`` are recalculated and ``Minimum`` and
``Maximum`` fields are updated to show real current max and min values.

.. image:: _static/screenshots/16426976574.png
.. image:: _static/screenshots/16426976633.png

The ``RMSD Sieve`` lets you mathematically compare structures of conformers and trim
duplicates and almost-duplicates. RMSD stands for root-mean-square deviation of atomic
positions and is a conformers similarity measure. The sieve calculates the average
distance between atoms of two conformers and trims the less stable (higher-energy)
conformer of the two, if the resulting RMSD value is smaller than value ot the
``Threshold`` field.

Calculating an RMSD value is quite resource-costly. To assure efficient trimming, each
conformer is compared only with conformers inside its energy window, defied by the
``Window size`` filed value. Conformers of energy this much higher or lower are
automatically considered different.

.. image:: _static/screenshots/16426977988.png

Temperature of the system
-------------------------

The ``Energies and structure`` section also allows you to specify the temperature
of the studied system. This parameter is important for calculation of the Boltzmann
distribution of conformers, which is used to estimate conformers' population
and average conformers' spectra. The default value is the room temperature,
expressed as :math:`298.15\ \mathrm{Kelvin}` (:math:`25.0^{\circ}\mathrm{C}`).
Changing this value will trigger automatic recalculation of ``Min. Boltzmann factor``
and ``Population`` values, and average spectra will be redrawn.

.. image:: _static/screenshots/16426977988.png

.. versionadded:: 0.9.1
    The ``Temperature`` entry allowing to change the temperature value.

Spectra simulation
------------------

``Calculate Spectra`` controls section and ``Spectra view`` tab allow to preview the
simulation of selected spectrum type with given parameters.

.. image:: _static/screenshots/16426978290.png

The ``Spectra view`` tab is initially empty, but when you select one of the available
``Spectra type``\s, ``Settings`` subsection will become enabled and the spectrum will be
drawn.

.. tip::

    You can turn off automatic recalculation of the spectrum by unchecking the ``Live
    preview`` box.

.. image:: _static/screenshots/16426978558.png

Beginning and end of the simulated spectral range may be set using ``Start`` and
``Stop`` fields. The view on the right will match these boundaries. Please note that
``Start`` must have lower value than ``Stop``. There is also a ``Step`` field that
allows you to adjust points density in the simulated spectrum.

.. image:: _static/screenshots/16426979133.png

``Width`` field defines a peak width in the simulated spectrum. It exact meaning depends
on the chosen fitting function (see below). For gaussian fitting ``Width`` is
interpreted as **half width of the peak at** :math:`\frac{1}{e}` **of its maximum** value
(HW1OeM). For lorentzian function it is interpreted as **half width at half maximum**
height of the peak (HWHM).

.. tip::

    You may change fields' values with the mouse wheel. Point the field with mouse
    cursor and allow for a small delay before switching form the scroll mode to the
    value-changing mode. Move the mouse cursor away from the field to switch back.

.. image:: _static/screenshots/16426979452.png

Finally, you may choose the fitting function used to simulate the spectrum from the
calculated intensities values - this will have a big impact on simulated peaks' shape.
Two such functions are available: gaussian and lorentzian functions. Usually lorentzian
function is used to simulate vibrational spectra and gaussian function for electronic
spectra.

.. image:: _static/screenshots/16426979810.png

The default spectra preview is a ``Single file`` preview that allows you to see the
simulated spectrum for the selected conformer. You may change the conformer to preview
using the drop-down menu shown in the screenshot below.

.. image:: _static/screenshots/16426980108.png

When in a ``Single file`` preview, spectral activities used to simulate the spectrum are
also shown on the right. You may turn this off by unticking the ``Show activities`` box.

.. image:: _static/screenshots/16426980385.png

You can also preview an population-weighted average spectrum of all :term:`kept`
conformers, by selecting ``Average by energy``. The drop-down menu lets you select the
energies that ``tesliper`` should use to calculate conformers populations.

.. image:: _static/screenshots/16426980675.png

The final option is to show all :term:`kept` conformers at once by selecting ``Stack by
overview`` option. The drop-down menu allows to choose a color scheme for the stacked
spectra lines.

.. image:: _static/screenshots/16426980850.png

Comparing with experiment
-------------------------

It's possible to and an overlay with the experimental spectrum to ``Single file`` and
``Average by energy`` previews. To load an experimental spectrum, use ``Load from file``
button in the ``Experimental spectrum`` subsection. ``tesliper`` can read spectrum
in the .txt (or .xy) file format. Binary .spc formats are not supported.

.. image:: _static/screenshots/16426981100.png

When you choose the experimental spectrum file, it's curve is drown on the right with
respect to the ``Start`` and ``Stop`` bounds. Red color is used for the experiment.
In case of a significant difference in the magnitude of intensity in both spectra,
the second scale will be added to the drawing.

.. image:: _static/screenshots/16426982209.png

The scale of the simulated values may be automatically adjusted to roughly match the
experiment with the ``Auto-scale`` button. It may be also adjusted manually by changing
the value of the ``Scaling`` field.

.. image:: _static/screenshots/16426982533.png
.. image:: _static/screenshots/16426982888.png

Similarly, ``Auto-shift`` button and ``Offset`` field let you to adjust simulated
spectrum's position on the x-axis. Positive ``Offset`` shifts the spectrum
bathochromically, a negative one shifts it hypsochromically. 

.. image:: _static/screenshots/16426982916.png
.. image:: _static/screenshots/16426982961.png

``Scaling`` and ``Offset`` values are remembered for the current spectra type, just like
the other parameters.

.. image:: _static/screenshots/16426984646.png

Data export
-----------

Calculated and extracted data may be exported to disk in three different formats: text
files with ``Export to .txt`` button, csv files with ``Export to .csv`` button and
Excel files with ``Export to .xlsx`` button. Clicking on any of those will bring up
the ``Export...`` dialog.

.. image:: _static/screenshots/16426985398.png

At the top of the ``Export...`` dialog is displayed the path to the currently selected
output directory. It may be changed by clicking on the ``Browse`` button and selecting
a new destination. Files generated by ``tesliper`` will be written to this directory.

.. image:: _static/screenshots/16426986807.png

On the left side of the dialog window you may select what type of data you want to
export by ticking appropriate boxes. Once you hover over the certain category, more
detailed list of available data will be shown on the right. By ticking/unticking
selected boxes you can fine-tune what should be written to disk.

.. image:: _static/screenshots/16426986886.png

In the ``Spectra`` category, beside each available spectra type, there is a note that
informs if calculation parameters were altered by the user. Spectra will be recalculated
with current parameters upon the export confirmation.

.. image:: _static/screenshots/16426986926.png

Creating Gaussian input
-----------------------

Clicking on the ``Create .gjf files...`` will open a dialog window that lets you setup
a next step of calculations to conduct with the Gaussian software.
    
.. image:: _static/screenshots/16426987223.png

Similarly to the previews one, this dialog also features a ``Path`` field that specifies
the output directory, which may be changed by clicking on the ``Browse`` button. Bellow
it is the ``Geometry type`` drop-down menu that allows you to select, which geometry
specification should be used in the new input files. ``Input`` is the geometry used as
an input in the extracted .log/.out files, ``Last read`` is the one that was lastly
encountered in these files. ``Optimized`` is the geometry marked as optimized by
Gaussian, but it is only available from the successful optimization calculations.
You also need to specify the ``Charge`` and the ``Multiplicity`` of the molecule.

.. image:: _static/screenshots/16426987882.png

Below are the ``Route`` and ``Comment`` fields. The first one specifies the calculation
directives for the Gaussian software. The second one is a title section required by
GAussian.

.. image:: _static/screenshots/16426988924.png

Further below is the expandable ``Link0 commands`` panel that allows to specify Link 0
directives, which define location of scratch files, memory usage, etc. Select a command
name from the drop-down menu, filed on right will show a hint about its purpose.

.. image:: _static/screenshots/16426989388.png

Provide a value in the input filed and click a ``+`` button to add a command. It will be
added to the list below. You can update the selected command by providing a new value
and clicking the ``+`` button again or remove it by clicking the ``-`` button.

.. image:: _static/screenshots/16426989597.png

Path-like commands may be parametrized: ``${conf}`` will be substituted with the name of
conformer and ``${num}`` will be substituted with the sequential number.

.. image:: _static/screenshots/16426990467.png

Finally, you can add a post-geometry specification. It will be written to the end of
each .gjf file.

.. image:: _static/screenshots/16426993446.png

Saving session
--------------

You can save a session (all data, along with current trimmed and parameters) with a
``Save session`` button. A popup dialog will be opened, where you can specify a target
session file location.

.. image:: _static/screenshots/16426994047.png

To load previously saved session use the ``Load session`` button. You can also discard
all currently held data by clicking the ``Clear session`` button.

.. warning::

    Loading and clearing session cannot be undone! A confirmation dialog will be
    displayed for those actions.
{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}
  :members:
  :member-order: bysource
  :inherited-members:

  {% block attributes %}
  {% if attributes %}
  .. rubric:: Module Attributes

  .. autosummary::
  {% for item in attributes %}
    {{ item }}
  {%- endfor %}
  {% endif %}
  {% endblock %}

  {% block functions %}
  {% if functions %}
  .. rubric:: {{ _('Functions') }}

  .. autosummary::
  {% for item in functions %}
    {{ item }}
  {%- endfor %}
  {% endif %}
  {% endblock %}

  {% block classes %}
  {% if classes %}
  .. rubric:: {{ _('Classes') }}

  .. autosummary::
  {% for item in classes %}
    {{ item }}
  {%- endfor %}
  {% endif %}
  {% endblock %}

  {% block exceptions %}
  {% if exceptions %}
  .. rubric:: {{ _('Exceptions') }}

  .. autosummary::
  {% for item in exceptions %}
    {{ item }}
  {%- endfor %}
  {% endif %}
  {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}