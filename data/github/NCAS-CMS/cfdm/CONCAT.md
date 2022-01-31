To update the master (publicly viewable) documentation run, from
within the `cfdm/docs` directory:

  ``make html .``

To create a development version of the documentation on your local
machine run, from within the `cfdm/docs` directory:

  ``make html dev``

To create an archive of a particular version's development on your
local machine run, from within the `cfdm/docs`
directory, for example:

  ``make html 1.7.1``
cfdm
====

A Python reference implementation of the CF data model.

[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/NCAS-CMS/cfdm?color=000000&label=latest%20version)](https://ncas-cms.github.io/cfdm/Changelog.html)
[![PyPI](https://img.shields.io/pypi/v/cfdm?color=000000)](https://pypi.org/project/cfdm/)
[![Conda](https://img.shields.io/conda/v/ncas/cfdm?color=000000)](https://ncas-cms.github.io/cfdm/installation.html#conda)

[![Conda](https://img.shields.io/conda/pn/ncas/cfdm?color=2d8659)](https://ncas-cms.github.io/cfdm/installation.html#operating-systems) [![Website](https://img.shields.io/website?color=2d8659&down_message=online&label=documentation&up_message=online&url=https%3A%2F%2Fncas-cms.github.io%2Fcfdm%2F)](https://ncas-cms.github.io/cfdm/index.html) [![GitHub](https://img.shields.io/github/license/NCAS-CMS/cfdm?color=2d8659)](https://github.com/NCAS-CMS/cfdm/blob/master/LICENSE)

[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/NCAS-CMS/cfdm/Run%20test%20suite?color=006666&label=test%20suite%20workflow)](https://github.com/NCAS-CMS/cfdm/actions) [![Codecov](https://img.shields.io/codecov/c/github/NCAS-CMS/cfdm?color=006666)](https://codecov.io/gh/NCAS-CMS/cfdm)

#### References

[![Website](https://img.shields.io/website?down_color=264d73&down_message=10.21105%2Fjoss.02717&label=JOSS&up_color=264d73&up_message=10.21105%2Fjoss.02717&url=https:%2F%2Fjoss.theoj.org%2Fpapers%2F10.21105%2Fjoss.02717%2Fstatus.svg)](https://doi.org/10.21105/joss.02717) [![Website](https://img.shields.io/website?color=264d73&down_message=10.5281%2Fzenodo.3894524&label=DOI&up_message=10.5281%2Fzenodo.3894524&url=https%3A%2F%2Fzenodo.org%2Frecord%2F3894524%23.Xuf2uXVKjeQ)](https://doi.org/10.5281/zenodo.3894524) [![Website](https://img.shields.io/website?down_color=264d73&down_message=10.5194%2Fgmd-10-4619-2017&label=GMD&up_color=264d73&up_message=10.5194%2Fgmd-10-4619-2017&url=https%3A%2F%2Fwww.geosci-model-dev.net%2F10%2F4619%2F2017%2F)](https://www.geosci-model-dev.net/10/4619/2017/)

#### Compliance with [FAIR principles](https://fair-software.eu/about/)

[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)


Documentation
=============

https://ncas-cms.github.io/cfdm

Tutorial
========

https://ncas-cms.github.io/cfdm/tutorial

Installation
============

https://ncas-cms.github.io/cfdm/installation

Functionality
=============

The ``cfdm`` package implements the CF data model
(https://doi.org/10.5194/gmd-10-4619-2017) for its internal data
structures and so is able to process any CF-compliant dataset. It is
not strict about CF-compliance, however, so that partially conformant
datasets may be ingested from existing datasets and written to new
datasets. This is so that datasets which are partially conformant may
nonetheless be modified in memory.

The central elements defined by the CF data model are the **field
construct**, which corresponds to CF-netCDF data variable with all of
its metadata; and the **domain contruct**, which may be the domain of
a field construct or corresponds to a CF-netCDF domain variable with
all of its metadata.

A simple example of reading a field construct from a file and
inspecting it:

    >>> import cfdm
    >>> f = cfdm.read('file.nc')
    >>> f
    [<Field: air_temperature(time(12), latitude(64), longitude(128)) K>]
    >>> print(f[0])
    Field: air_temperature (ncvar%tas)
    ----------------------------------
    Data            : air_temperature(time(12), latitude(64), longitude(128)) K
    Cell methods    : time(12): mean (interval: 1.0 month)
    Dimension coords: time(12) = [0450-11-16 00:00:00, ..., 0451-10-16 12:00:00] noleap
                    : latitude(64) = [-87.8638, ..., 87.8638] degrees_north
                    : longitude(128) = [0.0, ..., 357.1875] degrees_east
                    : height(1) = [2.0] m

The ``cfdm`` package can:

* read field and domain constructs from netCDF and CDL datasets,
* create new field and domain constructs in memory,
* write and append field and domain constructs to netCDF datasets on disk,
* read, write, and create coordinates defined by geometry cells,
* read and write netCDF4 string data-type variables,
* read, write, and create netCDF and CDL datasets containing
  hierarchical groups,
* inspect field and domain constructs,
* test whether two constructs are the same,
* modify field and domain construct metadata and data,
* create subspaces of field and domain constructs,
* incorporate, and create, metadata stored in external files, and
* read, write, and create data that have been compressed by convention
  (i.e. ragged or gathered arrays), whilst presenting a view of the
  data in its uncompressed form.

Command line utility
====================

During installation the `cfdump` command line tool is also installed,
which generates text descriptions of the field constructs contained in
a netCDF dataset:

    $ cfdump file.nc
    Field: air_temperature (ncvar%tas)
    ----------------------------------
    Data            : air_temperature(time(12), latitude(64), longitude(128)) K
    Cell methods    : time(12): mean (interval: 1.0 month)
    Dimension coords: time(12) = [0450-11-16 00:00:00, ..., 0451-10-16 12:00:00] noleap
                    : latitude(64) = [-87.8638, ..., 87.8638] degrees_north
                    : longitude(128) = [0.0, ..., 357.1875] degrees_east
                    : height(1) = [2.0] m

Tests
=====

Tests are run from within the ``cfdm/test`` directory:

    $ python run_tests.py

Citation
========

If you use cfdm, either as a stand-alone application or to provide a CF
data model implementation to another software library, please consider
including the reference:

Hassell et al., (2020). cfdm: A Python reference implementation of the
CF data model. Journal of Open Source Software, 5(54), 2717,
https://doi.org/10.21105/joss.02717

```
@article{Hassell2020,
  doi = {10.21105/joss.02717},
  url = {https://doi.org/10.21105/joss.02717},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2717},
  author = {David Hassell and Sadie L. Bartholomew},
  title = {cfdm: A Python reference implementation of the CF data model},
  journal = {Journal of Open Source Software}
}
```
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, sex
characteristics, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance,
race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive
environment include:

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

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable
behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that
they deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies within all project spaces, and it also
applies when an individual is representing the project or its
community in public spaces. Examples of representing a project or
community include posting via an official social media account, or
acting as an appointed representative at an online or offline
event. Representation of a project may be further defined and
clarified by project maintainers.

## Attribution

This Code of Conduct is adapted from the [Contributor
Covenant][homepage], version 1.4, available at
https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
Thank you for taking the time to consider making a contribution to the
cfdm package.

Please consult the
[dedicated guidance page](https://ncas-cms.github.io/cfdm/contributing.html)
in the documentation for detailed guidance on contributing to cfdm.

Overall, general questions, suggestions for enhancements, and reports
of bugs are appreciated and should be reported via the
[GitHub issue tracker](https://github.com/NCAS-CMS/cfdm/issues).

Code-change contributions to cfdm are also very welcome, though to
ensure the work is in line with plans for development of the library, please
always discuss any intended changes with the core development team in the
first instance. The standard way to do so is also through the
issue tracker:

* if there is an existing issue in the tracker than you would like
  to address, comment on that issue to indicate you would like to work on
  it;
* conversely, if there is not an issue corresponding to your ideas for
  contribution, please raise a new issue outlining the idea.

If you are not sure about any aspect related to contributing after
reading the guidance page, do not hesitate to get in touch by posting
a question in an issue.
* Change the version and date in `cfdm/core/__init__.py`
  (`__version__` and `__date__` variables); and in the `codemeta.json`
  file.

* Ensure that the requirements on dependencies and their versions are
  up-to-date and consistent in both the `requirements.txt` file and in
  `docs/source/installation.rst`; and in the `_requires` list and
  `LooseVersion` checks in `cfdm/core/__init__.py` and
  `cfdm/__init__.py`.

* If required, change the CF conventions version in
  `cfdm/core/__init__.py` (`__cf_version__` variable)

* Make sure that `README.md` is up to date.

* Make sure that the `long_description` in `setup.py` is up to date.

* Make sure that `Changelog.rst` is up to date.

* Check that the documentation API coverage is complete:

  ```bash
  ./check_docs_api_coverage
  ```

  * If it is not complete, add any undocumented attributes, methods,
    functions and keyword arguments (e.g. as listed in the change log)
    to the `.rst` files in `docs/source/class/`.

* Check external links to the CF conventions are up to date in
  `docs/source/tutorial.rst`

* Create a link to the new documentation in `docs/source/releases.rst`

* Test tutorial code:

  ```bash
  export PYTHONPATH=$PWD:$PYTHONPATH
  ./test_tutorial_code
  ```

* Build a development copy of the documentation using to check API
  pages for any new methods are present & correct, & that the overall
  formatting has not been adversely affected for comprehension by any
  updates in the latest Sphinx or theme etc. (Do not manually commit
  the dev build.)

  ```bash
  ./release_docs <vn> dev-clean # E.g. ./release_docs 1.8.7.0 dev-clean
  ```

* Check that no typos or spelling mistakes have been introduced to the
  documentation:

  * Run a dummy build of the documentation to detect invalid words:

     ```console
     $ cd docs
     $ make spelling build
     ```

  * If there are words raised with 'Spell check' warnings for the dummy
    build, such as:

    ```bash
    /home/sadie/cf-python/docs/source/class/cf.NetCDFArray.rst:18: Spell check: isw: element in the sequence isw the name of the group in which.
    Writing /home/sadie/cf-python/docs/spelling/class/cf.NetCDFArray.spelling
    /home/sadie/cf-python/docs/source/class/cf.Query.rst:3: Spell check: encapulates:  object encapulates a condition, such as.
    ```

    they may or may not be typos or mis-spellings. Address all the warnings
    (except those relating to files under `docs/source/class/`,
    `/attribute` or `/function` which will be fixed along with the origin
    docstrings after a 'latest' build) as follows:

    * If there are words that are in fact valid, add the valid words to
      the list of false positives for the spelling checker extension,
      `docs/source/spelling_false_positives.txt`.
    * Correct any words that are not valid in the codebase under `cfdm` or
      in the `docs/source` content files.

  * Note that, in the case there are many words raised as warnings, it
    helps to automate the above steps. The following commands are a means
    to do this processing:

    1. Copy all 'spell check' warnings (there will be 'Writing to ...' lines
       interspersed which can be removed by command so can be copied here too)
       output to STDOUT during the build to a file (here we use
       `spellings-file-1` as an example name).
    2. Cut all 'Writing to ...' lines interspersed with the warnings by
       running `sed -i '/^riting/d' spellings-file-1`.
    3. Cut all of the invalid words detected from the warning messages via
       `cat spellings-file-1 | cut -d':' -f 4 > spellings-file-2`
    4. Sift through these new words and remove any words that are true
       positives i.e. typos or mis-spellings. Correct them in the
       docstrings or documentation source files. If there are many
       instances across the docs, it helps to do a substitution of all
       occurences, e.g. via `find . -type f | xargs sed -i 's/<typo>/<correction>/g'`,
       though take care to have spaces surrounding words which may be
       part of other words, e.g. use
       `find . -type f | xargs sed -i 's/ ot / to /g'` to correct `ot` to `to`.
    5. Remove the leading whitespace character on each line and add
       all the new words to the current list of false positives:
       `sed 's/^.//' spellings-file-2 >> docs/source/spelling_false_positives.txt`
    6. Remove duplicate words and sort alphabetically via:
       `sort -u -o docs/source/spelling_false_positives.txt docs/source/spelling_false_positives.txt`

* Create an archived copy of the documentation:

  ```bash
  ./release_docs <vn> archive # E.g. ./release_docs 1.8.7.0 archive
  ```

* Update the latest documentation:

  ```bash
  ./release_docs <vn> latest # E.g. ./release_docs 1.8.7.0 latest
  ```

* Create a source tarball:

  ```bash
  python setup.py sdist
  ```

* Test the tarball release using

  ```bash
  ./test_release <vn> # E.g. ./test_release 1.8.7.0
  ```

* Push recent commits using

  ```bash
  git push origin master
  ```
  
* Tag the release:

  ```bash
  ./tag <vn> # E.g. ./tag 1.8.7.0
  ```
  
* Upload the source tarball to PyPi. Note this requires the `twine`
  library (which can be installed via `pip`) and relevant project
  privileges on PyPi.

  ```bash
  ./upload_to_pypi <vn> # E.g. ./upload_to_pypi 1.8.7.0
  ```

* Update the GitHub releases page for the new version:
  https://github.com/NCAS-CMS/cfdm/releases

* Upload the new release to Zenodo: https://zenodo.org/record/5521505
---
name: Bug report
about: Report something that is not working
title: 'Bug: '
labels: bug
assignees: ''

---

To report a bug, please provide:

* The version of the software and the environment in which you are encountering an issue. The output of `cfdm.environment(paths=False)` is useful for this.

* A description of the issue with, if possible,
  - what you expected to happen,
  - the steps needed to reproduce it,
  - any traceback information.
---
name: Feature request
about: Suggest an enhancement
title: 'Request: '
labels: enhancement
assignees: ''

---


---
name: Question
about: 'Ask a question, such as "how can I do this?", "why does it behave like that?",
  "how can I make it faster?", etc. '
title: 'Question: '
labels: question
assignees: ''

---


---
title: 'cfdm: A Python reference implementation of the CF data model'
tags:
  - CF
  - Python
  - metadata
  - climate
  - meteorology
  - oceanography
authors:
  - name: David Hassell
    orcid: 0000-0001-5106-7502
    affiliation: "1, 2"
  - name: Sadie L. Bartholomew
    orcid: 0000-0002-6180-3603
    affiliation: "1, 2" 
affiliations:
 - name: National Centre for Atmospheric Science, UK
   index: 1
 - name: University of Reading, UK
   index: 2
date: 24 July 2020
bibliography: paper.bib
---

# Summary

The `cfdm` open-source Python library [@Hassell:2020] implements the
data model [@Hassell:2017] of the CF (Climate and Forecast) metadata
conventions [@Eaton:2020] and so should be able to represent and
manipulate all existing and conceivable CF-compliant datasets.

The CF conventions are designed to promote the creation, processing,
and sharing of climate and forecasting data using Network Common Data
Form (netCDF) files and libraries [@Rew:1990; @Rew:2006]. They cater
for data from model simulations as well as from observations, made in situ
or by remote sensing platforms, of the planetary surface, ocean, and
atmosphere. For a netCDF data variable, they provide a description of
the physical meaning of data and of its spatial, temporal, and other
dimensional properties. The CF data model is an abstract
interpretation of the CF conventions that is independent of the netCDF
encoding.

The `cfdm` library has been designed as a stand-alone application,
e.g. as deployed in the pre-publication checks for the CMIP6 data request
[@Juckes:2020; @Eyring:2016], and also to provide a CF data model
implementation to other software libraries, such as
`cf-python` [@Hassell2:2020].

# Statement of need

The complexity of scientific datasets tends to increase with
improvements in scientific capabilities and it is essential that
software interfaces are able to understand new research outputs. To
the authors' knowledge, `cfdm` and software built on it are currently
the only libraries that can understand all CF-netCDF datasets, made
possible by the complete implementation of the CF data model. All
others omit facets that are not currently of interest to their
particular user communities.

# Functionality

NetCDF variables can be stored in a variety of representations
(including the use of compression techniques) but the CF data model,
and therefore `cfdm`, transcends the netCDF encoding to retain only the
logical structure. A key feature of `cfdm` is that the in-memory
representation and user-facing API are unaffected by the particular
choices made during dataset creation, which are often outside of the
user's control.

The latest version of the CF conventions (CF-1.8) is fully represented
by `cfdm`, including the recent additions of simple geometries
[@iso19125:2004] and netCDF group hierarchies.

The central element of the CF data model is the "field construct"
that encapsulates all of the data and metadata for a single
variable. The `cfdm` library can create field constructs ab initio, or
read them from netCDF files; inspect, subspace, and modify in memory;
and write them to CF-netCDF dataset files. As long as it can interpret
the data, `cfdm` does not enforce CF-compliance, allowing non-compliant
datasets to be read, processed, corrected, and rewritten.

This represents a limited functionality in comparison to other
software libraries used for analysis, which often include higher-level
functions such as those for regridding and statistical analysis, etc.
The decision to restrict the functionality was made for the following
reasons:

* The controlled functionality is sufficient for dataset inspection
  and creation, as well as for modifying non-CF-compliant datasets,
  activities that are an important part of both archive curation and
  data analysis workflows.

* An extended functionality could complicate the implementation,
  making it harder to update the library as the CF data model evolves.

* The anticipation is that other libraries will build on `cfdm`,
  inheriting its knowledge of the CF conventions and extending the API
  to add more sophisticated functions that are appropriate to their
  users (notably `cf-python`).

# Example usage

In this example, a netCDF dataset is read from disk and the resulting
field construct is inspected. The field construct is then subspaced,
has its standard name property changed, and finally is
re-inspected and written to a new dataset on disk:

```python
>>> import cfdm
>>> f = cfdm.read('file.nc')[0]
>>> print(f)
Field: specific_humidity (ncvar%q)
----------------------------------
Data            : specific_humidity(latitude(5), longitude(8)) 1
Cell methods    : area: mean
Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                : longitude(8) = [22.5, ..., 337.5] degrees_east
                : time(1) = [2019-01-01 00:00:00]
>>> g = f[0, 2:6]
>>> g.set_property('standard_name', 'relative humidity')
>>> print(g)
Field: relative humidity (ncvar%q)
----------------------------------
Data            : relative humidity(latitude(1), longitude(4)) 1
Cell methods    : area: mean
Dimension coords: latitude(1) = [-75.0] degrees_north
                : longitude(4) = [112.5, ..., 247.5] degrees_east
                : time(1) = [2019-01-01 00:00:00]
>>> cfdm.write(g, 'new_file.nc')
```	

# Evolution

The CF data model will evolve in line with the CF conventions and the
`cfdm` library will need to respond to such changes. To facilitate this,
there is a core implementation (`cfdm.core`) that defines an in-memory
representation of a field construct, with no further features. The
implementation of an enhancement to the CF data model would proceed
first with an independent update to the core implementation, then with
an update, outside of the inherited core implementation, to the
functionality for dataset interaction and further field construct
modification.

# Extensibility

To encourage other libraries to build on `cfdm`, it has been designed
to be subclassable so that the CF data model representation is easily
importable into third-party software. An important part of this
framework is the ability to inherit the mapping of CF data model
constructs to, and from, netCDF datasets. This is made possible by
use of the bridge design pattern [@Gamma:1995] that decouples the
implementation of the CF data model from the netCDF encoding so that
the two can vary independently. Such an inheritance is employed by the
`cf-python` library, which adds many metadata-aware analytical
capabilities and employs a more sophisticated data class. By
preserving the API of the `cfdm` data class, the `cf-python` data
class can be used within the inherited `cfdm` code base with almost no
modifications.

# Acknowledgements

We acknowledge Bryan Lawrence and Jonathan Gregory for advice on the
API and comments that greatly improved this manuscript; Allyn
Treshansky for suggesting improvements on the use of `cfdm` in other
libraries; and the CF community for their work on the CF conventions.

This work has received funding from the core budget of the UK National
Centre for Atmospheric Science, the European Commission Horizon 2020
programme (project "IS-ENES3", number 824084), the European Research
Council (project "Couplet", number 786427), and Research Councils
UK (project "UKFAFMIP", number NE/R000727/1).

# References
Version 1.9.?.?
---------------

**2021-??-??**

* Fixed bug that caused a `cfdm.write` failure when a vertical
  coordinate reference construct has no coordinates
  (https://github.com/NCAS-CMS/cfdm/issues/164)
* Fixed bug that caused a failure when downstream `identities` methods
  return an `itertools.chain` object
  (https://github.com/NCAS-CMS/cfdm/issues/170)

----
  
Version 1.9.0.1
---------------

**2021-10-12**

* Fixed bug that prevented some geometry coordinates being written to
  netCDF CLASSIC files (https://github.com/NCAS-CMS/cfdm/issues/140)
* Fixed bug that a caused segmentation fault when appending a string
  data type to netCDF files
  (https://github.com/NCAS-CMS/cfdm/issues/155)
* Fixed bug in `cf.Field.get_domain` when there are climatological
  time axes (https://github.com/NCAS-CMS/cfdm/issues/159)

----
  
Version 1.9.0.0
---------------

**2021-09-21**

* Python 3.6 support removed
  (https://github.com/NCAS-CMS/cfdm/issues/139)
* Conversion of `cfdm.Domain` to a non-abstract that may be read from
  and written to a netCDF dataset
  (https://github.com/NCAS-CMS/cfdm/issues/111)
* New method: `cfdm.Domain.creation_commands`
* New method: `cfdm.Domain.climatological_time_axes`
* New method: `cfdm.AuxiliaryCoordinate.del_climatology`
* New method: `cfdm.AuxiliaryCoordinate.get_climatology`
* New method: `cfdm.AuxiliaryCoordinate.is_climatology`
* New method: `cfdm.AuxiliaryCoordinate.set_climatology`
* New method: `cfdm.DimensionCoordinate.del_climatology`
* New method: `cfdm.DimensionCoordinate.get_climatology`
* New method: `cfdm.DimensionCoordinate.is_climatology`
* New method: `cfdm.DimensionCoordinate.set_climatology`
* New function: `cfdm.unique_constructs`
* New function: `cfdm.example_fields`
* Construct access API changes from 1.8.9.0 applied to `Field.convert`
* Improved error message for invalid inputs to `Field.convert`
* Raise exception when attempting to write multiply defined coordinate
  reference parameters (https://github.com/NCAS-CMS/cfdm/issues/148)
* Interpret format specifiers for size 1 `cfdm.Data` arrays
  (https://github.com/NCAS-CMS/cfdm/issues/152)
* Fix file name expansions in `cfdm.write`
  (https://github.com/NCAS-CMS/cfdm/issues/157)
  
----

Version 1.8.9.0
---------------

**2021-05-25**

* Construct access API changes
  (https://github.com/NCAS-CMS/cfdm/issues/124,
  https://github.com/NCAS-CMS/cfdm/issues/130,
  https://github.com/NCAS-CMS/cfdm/issues/132,
  https://github.com/NCAS-CMS/cfdm/issues/137)
* Performance enhancements
  (https://github.com/NCAS-CMS/cfdm/issues/124,
  https://github.com/NCAS-CMS/cfdm/issues/130)
* New write mode ``mode='a'`` for appending to, rather than over-writing,
  a netCDF file on disk (https://github.com/NCAS-CMS/cfdm/issues/143)
* Better error message in the case of a `numpy.ma.core.MaskError` occurring
  upon reading of CDL files with only header or coordinate information
  (https://github.com/NCAS-CMS/cfdm/issues/128)
* Fix for zero-sized unlimited dimensions when read from a grouped
  netCDF file (https://github.com/NCAS-CMS/cfdm/issues/113)
* Fix bug causing occasional non-symmetric `equals` operations
  (https://github.com/NCAS-CMS/cfdm/issues/133)
* Changed dependency: ``cftime>=1.5.0``
* Changed dependency: ``netCDF4>=1.5.4``

----

Version 1.8.8.0
---------------

**2020-12-18**

* The setting of global constants can now be controlled by a context
  manager (https://github.com/NCAS-CMS/cfdm/issues/100)
* Fixed bug that caused a failure when writing a dataset that contains
  a scalar domain ancillary construct
  (https://github.com/NCAS-CMS/cfdm/issues/98)
* Changed dependency: ``cftime>=1.3.0``

----

Version 1.8.7.0
---------------

**2020-10-09**

* Python 3.5 support deprecated (3.5 was retired on 2020-09-13)
* New method: `cfdm.Field.creation_commands`
* New method: `cfdm.Data.creation_commands`
* New method: `cfdm.Field._docstring_special_substitutions`
* New method: `cfdm.Field._docstring_substitutions`
* New method: `cfdm.Field._docstring_package_depth`
* New method: `cfdm.Field._docstring_method_exclusions`
* New method: `cfdm.Data.filled`
* New keyword parameter to `cfdm.Field.set_data`: ``inplace``
* New keyword parameter to `cfdm.write`: ``coordinates``
  (https://github.com/NCAS-CMS/cfdm/issues/81)
* New class: `cfdm.core.DocstringRewriteMeta`
* Comprehensive documentation coverage of class methods.
* Improved documentation following JOSS review.
* Enabled "creation commands" methods
  (https://github.com/NCAS-CMS/cfdm/issues/53)
* Fixed bug that caused failures when reading or writing a dataset
  that contains multiple geometry containers
  (https://github.com/NCAS-CMS/cfdm/issues/65)
* Fixed bug that prevented the writing of multiple fields to netCDF when
  at least one dimension was shared between some of the fields.

----

Version 1.8.6.0
---------------

**2020-07-24**

* Removed Python 2.7 support
  (https://github.com/NCAS-CMS/cfdm/issues/55)
* Implemented the reading and writing of netCDF4 group hierarchies for
  CF-1.8 (https://github.com/NCAS-CMS/cfdm/issues/13)
* Renamed to lower-case (but otherwise identical) names all functions
  which get and set global constants: `cfdm.atol`, `cfdm.rtol`,
  `cfdm.log_level`. The old names e.g. `cfdm.ATOL` remain functional
  as aliases.
* New function: `cfdm.configuration`
* New method: `cfdm.Field.nc_variable_groups`
* New method: `cfdm.Field.nc_set_variable_groups`
* New method: `cfdm.Field.nc_clear_variable_groups`
* New method: `cfdm.Field.nc_group_attributes`
* New method: `cfdm.Field.nc_set_group_attribute`
* New method: `cfdm.Field.nc_set_group_attributes`
* New method: `cfdm.Field.nc_clear_group_attributes`
* New method: `cfdm.Field.nc_geometry_variable_groups`
* New method: `cfdm.Field.nc_set_geometry_variable_groups`
* New method: `cfdm.Field.nc_clear_geometry_variable_groups`
* New method: `cfdm.DomainAxis.nc_dimension_groups`
* New method: `cfdm.DomainAxis.nc_set_dimension_groups`
* New method: `cfdm.DomainAxis.nc_clear_dimension_groups`
* New method: `cfdm.AuxiliaryCoordinate.del_interior_ring`
* New keyword parameter to `cfdm.write`: ``group``
* Keyword parameter ``verbose`` to multiple methods now accepts named
  strings, not just the equivalent integer levels, to set verbosity.
* Added test to check that cell bounds have more dimensions than the
  data.
* Added test to check that dimension coordinate construct data is
  1-dimensional.
* Fixed bug in `cfdm.CompressedArray.to_memory`.
* Fixed bug that caused an error when a coordinate bounds variable is
  missing from a dataset (https://github.com/NCAS-CMS/cfdm/issues/63)
* New dependency: ``netcdf_flattener>=1.2.0``
* Changed dependency: ``cftime>=1.2.1``
* Removed dependency: ``future``

----

Version 1.8.5
-------------

**2020-06-10**

* Fixed bug that prevented the reading of certain netCDF files, such
  as those with at least one external variable.

----

Version 1.8.4
-------------

**2020-06-08**

* Added new example field ``7`` to `cfdm.example_field`.
* Enabled configuration of the extent and nature of informational and
  warning messages output by `cfdm` using a logging framework (see
  points below) (https://github.com/NCAS-CMS/cfdm/issues/31)
* New function `cfdm.LOG_LEVEL` to set the minimum log level for which
  messages are displayed globally, i.e. to change the project-wide
  verbosity (https://github.com/NCAS-CMS/cfdm/issues/35).
* Changed behaviour and default of `verbose` keyword argument when
  available to a function/method so it interfaces with the new logging
  functionality (https://github.com/NCAS-CMS/cfdm/issues/35).
* Changed dependency: ``cftime>=1.1.3``
* Fixed bug the wouldn't allow the reading of a netCDF file which
  specifies Conventions other than CF
  (https://github.com/NCAS-CMS/cfdm/issues/36).

----

Version 1.8.3
-------------

**2020-04-30**

* `cfdm.Field.apply_masking` now masks metadata constructs.
* New method: `cfdm.Field.get_filenames`
* New method: `cfdm.Data.get_filenames`
* New function: `cfdm.abspath`
* New keyword parameter to `cfdm.read`: ``warn_valid``
  (https://github.com/NCAS-CMS/cfdm/issues/30)
* New keyword parameter to `cfdm.write`: ``warn_valid``
  (https://github.com/NCAS-CMS/cfdm/issues/30)
  

----

Version 1.8.2
-------------

**2020-04-24**

* Added time coordinate bounds to the polygon geometry example field
  ``6`` returned by `cfdm.example_field`.
* New method: `cfdm.Field.apply_masking`
* New method: `cfdm.Data.apply_masking`
* New keyword parameter to `cfdm.read`: ``mask``
* New keyword parameter to `cfdm.Field.nc_global_attributes`:
  ``values``
* Fixed bug in `cfdm.write` that caused (what are effectively)
  string-valued scalar auxiliary coordinates to not be written to disk
  as such, or even an exception to be raised.
  
----

Version 1.8.1
-------------

**2020-04-16**

* Improved source code highlighting in links from the documentation
  (https://github.com/NCAS-CMS/cfdm/issues/21).
* Fixed bug that erroneously required netCDF geometry container
  variables to have a ``geometry_dimension`` netCDF attribute.

----

Version 1.8.0
-------------

**2020-03-23**

* First release for CF-1.8 (does not include netCDF hierarchical
  groups functionality).
* Implementation of simple geometries for CF-1.8
  (https://github.com/NCAS-CMS/cfdm/issues/11).
* Implementing of string data-types for CF-1.8
  (https://github.com/NCAS-CMS/cfdm/issues/12).
* New function: `cfdm.example_field`
  (https://github.com/NCAS-CMS/cfdm/issues/18)
* New attributes: `cfdm.Field.dtype`, `cfdm.Field.ndim`,
  `cfdm.Field.shape`, `cfdm.Field.size`
* New method: `cfdm.Data.any`
* New ``paths`` keyword parameter to `cfdm.environment`
* Changed dependency: ``netCDF4>=1.5.3``
* Changed dependency: ``cftime>=1.1.1``
* Fixed bug that prevented the writing of ``'NETCDF3_64BIT_OFFSET'``
  and ``'NETCDF3_64BIT_DATA'`` format files
  (https://github.com/NCAS-CMS/cfdm/issues/9).
* Fixed bug that caused a failure when a "_FillValue" or
  "missing_value" property is set and data type conversions are
  specified with the ``datatype`` keyword to `cfdm.write`
  (https://github.com/NCAS-CMS/cfdm/issues/16).
* Fixed bug whereby `cfdm.Field.has_construct` would try to delete the
  construct rather than check whether it existed.

----

Version 1.7.11
--------------

**2019-11-27**

* New methods: `cfdm.Field.compress`, `cfdm.Field.uncompress`
* New methods: `cfdm.Data.flatten`, `cfdm.Data.uncompress`
* New  ``dtype`` and ``mask`` keyword parameters to `cfdm.Data`
* Changed the default value of the ``ignore_compression`` parameter to
  `True`.

----

Version 1.7.10
--------------

**2019-11-14**

* New method: `cfdm.Field.nc_set_global_attributes`.
* Fixed bug relating to the reading of some CDL files
  (https://github.com/NCAS-CMS/cfdm/issues/5).
* Fixed bug relating numpy warning when printing a field with masked
  reference time values (https://github.com/NCAS-CMS/cfdm/issues/8).

----

Version 1.7.9
-------------

**2019-11-07**

* Fixed bug relating to setting of parameters on datum and coordinate
  conversion objects of coordinate conversion constructs
  (https://github.com/NCAS-CMS/cfdm/issues/6).

----

Version 1.7.8
-------------

**2019-10-04**

* During writing to netCDF files, ensured that _FillValue and
  missing_value have the same data type as the data.
* Fixed bug during construct equality testing that didn't recognise
  equal cell method constructs in transposed, but otherwise equal
  field constructs.
* Bounds netCDF dimension name is now saved, and can be set. The
  saved/set value is written out to disk.
* Now reads CDL files (https://github.com/NCAS-CMS/cfdm/issues/5)

----

Version 1.7.7
-------------

**2019-06-13**

* Don't set the fill mode for a `netCDF4.Dataset` open for writing to
  `off`, to prevent incorrect reading of some netCDF4 files
  (https://github.com/NCAS-CMS/cfdm/issues/4).
* Updated documentation
  
----

Version 1.7.6
-------------

**2019-06-05**

* Added attributes `_ATOL` and `_RTOL` to facilitate subclassing.
* Fixed bug in `cfdm.Field.convert`.
* Fixed bug in `cfdm.core.constructs.new_identifier`.
  
----

Version 1.7.5
-------------

**2019-05-15**

* New methods: `Datum.nc_has_variable`, `Datum.nc_get_variable`,
  `Datum.nc_has_variable`, `Datum.nc_set_variable`
  (https://github.com/NCAS-CMS/cfdm/issues/3).
  
----

Version 1.7.4
-------------

**2019-05-14**

* Changed behaviour of `cfdm.Constructs.filter_by_axis`.
* New methods: `cfdm.Data.has_units`, `cfdm.Data.has_calendar`,
  `cfdm.Data.has_fill_value`.
* New ``constructs`` keyword parameter to `Field.transpose`.
* Keyword parameter ``axes`` to `cfdm.Field.set_data` is now optional.
* Added the 'has_bounds' method to constructs that have data but can't
  have bounds.
* New methods: `cfdm.DomainAxis.nc_is_unlimited`,
  `cfdm.DomainAxis.nc_set_unlimited`.
* Made Data a virtual subclass of Array.   
* Deprecated methods: `cfdm.Field.nc_unlimited`,
  `cfdm.Field.nc_clear_unlimited`, `cfdm.Field.nc_clear_unlimited`.
* Fixed bug when writing new horizontal coordinate reference for the
  vertical datum.
* Fixed bug in `del_data` methods.
* Fixed bug with in-place operations.
* Fixed bug with position in some `insert_dimension` methods.
* Fixed bug that sometimes made duplicate netCDF dimensions when
  writing to a file.
* Added _shape keyword to `cfdm.Field.set_data_axes` to allow the data
  shape to be checked prior to insertion.
* Added the '_custom' attribute to facilitate subclassing.
* New class `cfdm.mixin.NetCDFUnlimitedDimension` replaces
  `cfdm.mixin.NetCDFUnlimitedDimensions`, which is deprecated.
* New method `cfdm.CFDMImplementation.nc_is_unlimited_axis` replaces
  `cfdm.CFDMImplementation.nc_get_unlimited_axes`, which is
  deprecated.
* New method `cfdm.CFDMImplementation.nc_set_unlimited_axis` replaces
  `cfdm.CFDMImplementation.nc_set_unlimited_dimensions`, which is
  deprecated.
  
----

Version 1.7.3
-------------

**2019-04-24**

* New method: `cfdm.Constructs.filter_by_size`.
* New method: `cfdm.Data.uncompress`.
* Changed the default behaviours of the
  `cfdm.Construct.filter_by_axis`, `cfdm.Construct.filter_by_size`,
  `cfdm.Construct.filter_by_naxes`,
  `cfdm.Construct.filter_by_property`,
  `cfdm.Construct.filter_by_ncvar`, `cfdm.Construct.filter_by_ncdim`,
  `cfdm.Construct.filter_by_method`,
  `cfdm.Construct.filter_by_measure` methods in the case when no
  arguments are provided: Now returns all possible constructs that
  *could* have the feature, with any values.
* Renamed the "underlying_array" methods to "source"
* Added _field_data_axes attribute to `Constructs` instances.
* Added _units and _fill_value arguments to get_data method.
* Moved contents of cfdm/read_write/constants.py to `NetCDFRead` and
  `NetCDFWrite`.
* Fixed bug in `cfdm.CoordinateReference.clear_coordinates`.
* Fixed bug in `cfdm.Field.convert` (which omitted domain ancillaries
  in the result).
* Added ``kwargs`` parameter to
  `cfdm.CFDMImplementation.initialise_Data`, to facilitate
  subclassing.
* Added `NetCDFRead._customize_read_vars` to facilitate subclassing.
* Added `NetCDFWrite._transform_strings` to facilitate subclassing.

----

Version 1.7.2
-------------

**2019-04-05**

* New ``mode`` parameter options to `cfdm.Constructs.filter_by_axis`:
  ``'exact'``, ``'subset'``, ``'superset'``.
* Enabled setting of HDF5 chunksizes.
* Fixed bug that caused coordinate bounds to be not sliced during
  subspacing (https://github.com/NCAS-CMS/cfdm/issues/1).

----

Version 1.7.1
-------------

**2019-04-02**

* New methods `cfdm.Constructs.clear_filters_applied`,
  `cfdm.Constructs.filter_by_naxes`.
* Changed behaviour of `cfdm.Constructs.unfilter` and
  `cfdm.Constructs.inverse_filters`: added depth keyword and changed
  default.

----

Version 1.7.0
-------------

**2019-04-02**

* First release for CF-1.7

----
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autofunction:: {{ fullname }}
:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. automethod:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autoattribute:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. auto{{ objtype }}:: {{ fullname }}
.. currentmodule:: cfdm
.. default-role:: obj

.. _Contributing:

**Contributing**
================

----

Version |release| for version |version| of the CF conventions.

.. contents::
   :local:
   :backlinks: entry

**Reporting bugs**
------------------

Please report bugs via a new issue in issue tracker
(https://github.com/NCAS-CMS/cfdm/issues), using the **Bug report**
issue template.

----

**Feature requests and suggested improvements**
-----------------------------------------------

Suggestions for new features and any improvements to the
functionality, API, documentation and infrastructure can be submitted
via a new issue in issue tracker
(https://github.com/NCAS-CMS/cfdm/issues), using the **Feature
request** issue template.

----

**Questions**
-------------

Questions, such as "how can I do this?", "why does it behave like
that?", "how can I make it faster?", etc., can be raised via a new
issue in issue tracker (https://github.com/NCAS-CMS/cfdm/issues),
using the **Question** issue template.

----

**Preparing pull requests**
---------------------------

Pull requests should follow on from a discussion in the issue tracker
(https://github.com/NCAS-CMS/cfdm/issues).

Fork the cfdm GitHub repository (https://github.com/NCAS-CMS/cfdm).

Clone your fork locally and create a branch:

.. code-block:: console
	  
    $ git clone git@github.com:<YOUR GITHUB USERNAME>/cfdm.git
    $ cd cfdm
    $ git checkout -b <your-bugfix-feature-branch-name master>

Break your edits up into reasonably-sized commits, each representing
a single logical change:

.. code-block:: console
	  
    $ git commit -a -m "<COMMIT MESSAGE>"

Create a new changelog entry in ``Changelog.rst``. The entry should be
written (where ``<description>`` should be a *brief* description of
the change) as:

.. code-block:: rst

   * <description> (https://github.com/NCAS-CMS/cfdm/issues/<issue number>)

Run the test suite to make sure the tests all pass:
	
.. code-block:: console

   $ cd cfdm/test
   $ python run_tests.py

Add your name to the list of contributors list at
``docs/source/contributing.rst``.

Finally, make sure all commits have been pushed to the remote copy of
your fork and submit the pull request via the GitHub website, to the
``master`` branch of the ``NCAS-CMS/cfdm`` repository. Make sure to
reference the original issue in the pull request's description.

Note that you can create the pull request while you're working on
this, as it will automatically update as you add more commits. If it is
a work in progress, you can mark it initially as a draft pull request.

----

**Contributors**
----------------

We would like to acknowledge and thank all those who have contributed
ideas, code, and documentation to the cfdm library:

* Alan Iwi
* Allyn Treshansky
* Bryan Lawrence
* David Hassell
* Jonathan Gregory
* Martin Juckes
* Riley Brady  
* Sadie Bartholomew  
.. currentmodule:: cfdm
.. default-role:: obj

**Change log**
==============

.. contents::
   :local:
   :backlinks: entry

.. include:: ../../Changelog.rst
.. currentmodule:: cfdm
.. default-role:: obj

.. raw:: html

    <style> .small {font-size:small} </style>

.. role:: small

.. _CF-data-model:

**CF data model**
=================

----

The CF (Climate and Forecast) metadata conventions
(http://cfconventions.org) provide a description of the physical
meaning of data and of their spatial and temporal properties and are
designed to promote the creation, processing, and sharing of climate
and forecasting data using netCDF files and libraries
(https://www.unidata.ucar.edu/software/netcdf).

`The CF data model
<https://cfconventions.org/cf-conventions/cf-conventions.html#appendix-CF-data-model>`_
identifies the fundamental elements ("constructs") of the CF
conventions and shows how they relate to each other, independently of
the netCDF encoding.

The CF data model defines a **field construct** for storing data with
all of its metadata. It is defined as follows:

.. glossary::

  field construct
    corresponds to a CF-netCDF data variable with all of its
    metadata. It consists of

    - descriptive properties that apply to field construct as a whole
      (e.g. the standard name),
    
    - a data array,

    - a **domain construct** that describes the locations of each cell
      of the data array (i.e. the "domain"),
    
    - **metadata constructs** that describe the physical nature of the
      data array, defined by

      .. glossary::
         
        field ancillary constructs
          corresponding to CF-netCDF ancillary variables
      
        cell method constructs
          corresponding to a CF-netCDF cell_methods attribute of data
          variable

  domain construct
    that describes the locations of each cell of the domain. It may
    exist independently of a **field construct** and consists of

    - descriptive properties that apply to domain construct as a whole,
    
    - **metadata constructs** that describe the locations of each cell
      of the domain, defined by 
    
    .. glossary::
         
      domain axis constructs
        corresponding to CF-netCDF dimensions or scalar coordinate
        variables
    
      dimension coordinate constructs
        corresponding to CF-netCDF coordinate variables or numeric
        scalar coordinate variables
    
      auxiliary coordinate constructs
        corresponding to CF-netCDF auxiliary coordinate variables and
        non-numeric scalar coordinate variables
    
      coordinate reference constructs
        corresponding to CF-netCDF grid mapping variables or the
        formula_terms attribute of a coordinate variable
    
      domain ancillary constructs
        corresponding to CF-netCDF variables named by the
        formula_terms attribute of a coordinate variable
    
      cell measure constructs
        corresponding to CF-netCDF cell measure variables

----

|

.. figure:: images/cfdm_field.svg
   :scale: 35 %

   *The constructs of the CF data model described using UML. The field construct corresponds to a CF-netCDF data variable. The domain construct provides the linkage between the field construct and the constructs which describe measurement locations and cell properties. It is useful to define an abstract generic coordinate construct that can be used to refer to coordinates when the their type (dimension or auxiliary coordinate construct) is not an issue.*

----
.. currentmodule:: cfdm
.. default-role:: obj


.. _Performance:

**Performance**
===============
----

Version |release| for version |version| of the CF conventions.


.. contents::
   :local:
   :backlinks: entry

.. _Memory:

**Memory**
----------

When a dataset is read using `cfdm.read`, `lazy loading
<https://en.wikipedia.org/wiki/Lazy_loading>`_ is employed for all
data arrays, which means that no data is read into memory until the
data is required for inspection or to modify the array contents. This
maximises the number of :term:`field constructs <field construct>`
that may be read within a session, and makes the read operation
fast. If a :ref:`subspace <Subspacing>` of data still in the file is
requested then only that subspace is read into memory. These
behaviours are inherited from the `netCDF4 python package
<http://unidata.github.io/netcdf4-python/netCDF4/index.html>`_.

When an instance is copied with its `!copy` method, all data are
copied with a `copy-on-write
<https://en.wikipedia.org/wiki/Copy-on-write>`_ technique. This means
that a copy takes up very little memory, even when the original data
comprises a very large array in memory, and the copy operation is
fast.

----

.. _In-place-operations:

**In-place operations**
-----------------------

Some methods that create new a instance have an option to perform the
operation in-place, rather than creating a new independent object. The
in-place operation can be considerably faster. These methods have the
``inplace`` keyword parameter, such as the `~Field.squeeze`,
`~Field.transpose`, `~Field.insert_dimension`, `~Field.compress`, and
`~Field.uncompress` methods of a field construct.
  
For example, in one particular test, transposing the data dimensions
of the field construct was ~10 times faster when done in-place,
compared with creating a new independent field construct:

.. code-block:: python
   :caption: *Calculate the speed-up of performing the "transpose"
             operation in-place.
      
   >>> import timeit
   >>> import cfdm
   >>> f = cfdm.example_field(0)
   >>> print(f)
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   >>> min(timeit.repeat('g = f.transpose()',
   ...                   globals=globals(), number=1000))
   1.2819487630004005
   >>> min(timeit.repeat('f.transpose(inplace=True)',
   ...                   globals=globals(), number=1000))
   0.13453567200122052

.. currentmodule:: cfdm
.. default-role:: obj


**API reference**
=================

----

Version |release| for version |version| of the CF conventions.

* **Construct classes**
  
  * :ref:`cfdm.Field <cfdm-Field>`
  * :ref:`cfdm.Domain <cfdm-Domain>`
  * :ref:`cfdm.AuxiliaryCoordinate <cfdm-AuxiliaryCoordinate>`
  * :ref:`cfdm.CellMeasure <cfdm-CellMeasure>`
  * :ref:`cfdm.CellMethod <cfdm-CellMethod>`
  * :ref:`cfdm.Coordinatereference <cfdm-Coordinatereference>`
  * :ref:`cfdm.DimensionCoordinate <cfdm-DimensionCoordinate>`
  * :ref:`cfdm.DomainAncillary <cfdm-DomainAncillary>`
  * :ref:`cfdm.DomainAxis <cfdm-DomainAxis>`
  * :ref:`cfdm.FieldAncillary <cfdm-FieldAncillary>`
    
.. toctree::
   :maxdepth: 2

   functions
   class
   constant
   class_core
.. currentmodule:: cfdm
.. default-role:: obj

.. _Sample-datasets:

**Sample datasets**
-------------------

This tutorial uses a number of small sample datasets, all of which can
be found in the zip file ``cfdm_tutorial_files.zip``
(:download:`download <../source/sample_files/cfdm_tutorial_files.zip>`,
8kB):
		    
.. code-block:: console
   :caption: *Unpack the sample datasets.*
		
   $ unzip -q cfdm_tutorial_files.zip
   $ ls -1
   cfdm_tutorial_files.zip
   contiguous.nc
   external.nc
   file.nc
   gathered.nc
   geometry.nc
   parent.nc
   
The tutorial examples assume that the Python session is being run from
the directory that contains the zip file and its unpacked contents,
and no other files.
   
----
.. currentmodule:: cfdm
.. default-role:: obj
 
.. _Releases:

**Releases**
============
----

Documentation for all versions of cfdm.

.. contents::
   :local:
   :backlinks: entry

**CF-1.9**
----------

* `Version 1.9.0.1 <https://ncas-cms.github.io/cfdm/1.9.0.1>`_ (2021-10-12)
* `Version 1.9.0.0 <https://ncas-cms.github.io/cfdm/1.9.0.0>`_ (2021-09-21)

----

**CF-1.8**
----------

* `Version 1.8.9.0 <https://ncas-cms.github.io/cfdm/1.8.9.0>`_ (2021-05-25)
* `Version 1.8.8.0 <https://ncas-cms.github.io/cfdm/1.8.8.0>`_ (2020-12-18)
* `Version 1.8.7.0 <https://ncas-cms.github.io/cfdm/1.8.7.0>`_ (2020-10-09)
* `Version 1.8.6.0 <https://ncas-cms.github.io/cfdm/1.8.6.0>`_ (2020-07-24)
* `Version 1.8.5 <https://ncas-cms.github.io/cfdm/1.8.5>`_ (2020-06-10)
* `Version 1.8.4 <https://ncas-cms.github.io/cfdm/1.8.4>`_ (2020-06-08)
* `Version 1.8.3 <https://ncas-cms.github.io/cfdm/1.8.3>`_ (2020-04-30)
* `Version 1.8.2 <https://ncas-cms.github.io/cfdm/1.8.2>`_ (2020-04-24)
* `Version 1.8.1 <https://ncas-cms.github.io/cfdm/1.8.1>`_ (2020-04-14)
* `Version 1.8.0 <https://ncas-cms.github.io/cfdm/1.8.0>`_ (2020-03-23)

----

**CF-1.7**
----------

* `Version 1.7.11 <https://ncas-cms.github.io/cfdm/1.7.11>`_ (2019-11-27)
* `Version 1.7.10 <https://ncas-cms.github.io/cfdm/1.7.10>`_ (2019-11-14)
* `Version 1.7.9 <https://ncas-cms.github.io/cfdm/1.7.9>`_ (2019-11-07)
* `Version 1.7.8 <https://ncas-cms.github.io/cfdm/1.7.8>`_ (2019-10-04)
* `Version 1.7.7 <https://ncas-cms.github.io/cfdm/1.7.7>`_ (2019-06-13)
* `Version 1.7.6 <https://ncas-cms.github.io/cfdm/1.7.6>`_ (2019-06-05)
* `Version 1.7.5 <https://ncas-cms.github.io/cfdm/1.7.5>`_ (2019-05-15)
* `Version 1.7.4 <https://ncas-cms.github.io/cfdm/1.7.4>`_ (2019-05-14)
* `Version 1.7.3 <https://ncas-cms.github.io/cfdm/1.7.3>`_ (2019-04-24)
* `Version 1.7.2 <https://ncas-cms.github.io/cfdm/1.7.2>`_ (2019-04-05)
* `Version 1.7.1 <https://ncas-cms.github.io/cfdm/1.7.1>`_ (2019-04-02)
* `Version 1.7.0 <https://ncas-cms.github.io/cfdm/1.7.0>`_ (2019-04-02)

----

.. _Versioning:

**Versioning**
--------------

Finding versions
^^^^^^^^^^^^^^^^

Version |release| for version |version| of the CF conventions.

The version of the CF conventions and the CF data model being used may
be found with the `cfdm.CF` function:

.. code-block:: python
   :caption: *Retrieve the version of the CF conventions.*
	     
   >>> import cfdm
   >>> cfdm.CF()
   '1.9'

This indicates which version of the CF conventions are represented by
this release of the cfdm package, and therefore the version can not be
changed.

The version identifier of the cfdm package is based on the version of
the CF conventions to which it applies, with the addition of extra
integer values for updates that apply to the same version of CF:

.. code-block:: python
   :caption: *Retrieve the version of the cfdm package.*
	     	     
   >>> cfdm.__version__
   '1.9.0.0'

The next section outlines the scheme used to set version identifiers.

Versioning strategy
^^^^^^^^^^^^^^^^^^^

A ``CF.major.minor`` numeric version scheme is used, where ``CF`` is
the version of the CF conventions (e.g. ``1.9``) to which a particular
version of cfdm applies.

**Major** changes comprise:

* changes to the API, such as:

  * changing the name of an existing function or method;
  * changing the behaviour of an existing function or method;
  * changing the name of an existing keyword parameter;
  * changing the default value of an existing keyword parameter;
  * changing the meaning of a value of an existing keyword parameter.
  * introducing a new function or method;
  * introducing a new keyword parameter;
  * introducing a new permitted value of a keyword parameter;

* changes to required versions of the dependencies.

**Minor** changes comprise:

* bug fixes that do not change the API;
* changes to the documentation;
* code tidying.
.. currentmodule:: cfdm
.. default-role:: obj

.. _Extensions:

**Extensions**
==============

----

Version |release| for version |version| of the CF conventions.

.. contents::
   :local:
   :backlinks: entry

The cfdm package has been designed to be subclassed, so that the
creation of a new implementation of the CF data model, based on cfdm,
is straight forward. For example:

.. code-block:: python
   :caption: *Create a new implementation with a new field construct
             class.*

   import cfdm
   
   class my_Field(cfdm.Field):
       def info(self):
           return 'I am a {!r} instance'.format(
               self.__class__.__name__)
   
The interpretation of CF-netCDF files that is encoded within the
`cfdm.read` and `cfdm.write` functions is also inheritable, so that an
extended data model implementation need *not* recreate the complicated
mapping of CF data model constructs to, and from, CF-netCDF
elements. This is made possible by the `bridge design pattern
<https://en.wikipedia.org/wiki/Bridge_pattern>`_, that decouples the
implementation of the CF data model from the CF-netCDF encoding so
that the two can vary independently.

.. code-block:: python
   :caption: *Define an implementation that is the same as cfdm, but
              which uses the my_Field class to represent field
              constructs*   
	      
   >>> my_implementation = cfdm.implementation()
   >>> my_implementation.set_class('Field', my_Field)

.. code-block:: python
   :caption: *Define new functions that can read into my_Field
              instances, and write write my_Field instances to
              datasets on disk.*
   
   >>> import functools
   >>> my_read = functools.partial(cfdm.read,
   ...                             _implementation=my_implementation)
   >>> my_write = functools.partial(cfdm.write,
   ...                              _implementation=my_implementation)

.. code-block:: python
   :caption: *Read my_field constructs from 'file.nc', the netCDF file
              used in the tutorial, using the new my_read
              function. Inspect a field construct read from the
              dataset, demonstrating that it is a my_Field instance
              from the new implementation that has the inherited
              functionality of a cfdm.Field instance.*

   >>> q, t = my_read('file.nc')
   >>> print(type(q))
   <class '__main__.my_Field'>  
   >>> print(q.info())
   I am a 'my_Field' instance
   >>> print(repr(q))
   <my_Field: specific_humidity(latitude(5), longitude(8)) 1>
   >>> print(q.data.array)
   [[0.007 0.034 0.003 0.014 0.018 0.037 0.024 0.029]
    [0.023 0.036 0.045 0.062 0.046 0.073 0.006 0.066]
    [0.11  0.131 0.124 0.146 0.087 0.103 0.057 0.011]
    [0.029 0.059 0.039 0.07  0.058 0.072 0.009 0.017]
    [0.006 0.036 0.019 0.035 0.018 0.037 0.034 0.013]]

.. code-block:: python
   :caption: *Write the my_field constructs to a netCDF file using the
             new my_write function.*
	     
   >>> my_write([q, t], 'new_file.nc')

Note that, so far, we have only replaced the field construct class in
the new implementation, and not any of the metadata constructs or
other component classes:
    
.. code-block:: python
   :caption: *Demonstrate that the metadata construct classes of
             within the my_Field instance are still cfdm classes.*

   >>> print(type(q))
   <class '__main__.my_Field'>  
   >>> print(type(q.construct('latitude')))
   <class 'cfdm.dimensioncoordinate.DimensionCoordinate'>

If the API of the new implementation is changed such that a given cfdm
functionality has a different API in the new implementation, then the
new read-from-disk and write-to-disk functions defined above can still
be used provided that the new implementation is created from a
subclass of `cfdm.CFDMImplementation`, with the new API being applied
in overridden methods.

.. code-block:: python
   :caption: *Create an implementation with a different API.*
   
   class my_Field_2(cfdm.Field):
      def my_coordinates(self):
          """Get coordinate constructs with a different API."""
          c = self.coordinates
          if not c:
              return {}
          return c
   
   class my_CFDMImplementation(cfdm.CFDMImplementation):
      def get_coordinates(self, field):
          """Get coordinate constructs from a my_Field_2 instance,
          using its different API.
          """
          return field.my_coordinates()
   
   my_implementation_2 = my_CFDMImplementation(
      cf_version=cfdm.CF(),
      
      Field=my_Field_2,
      
      AuxiliaryCoordinate=cfdm.AuxiliaryCoordinate,
      CellMeasure=cfdm.CellMeasure,
      CellMethod=cfdm.CellMethod,
      CoordinateReference=cfdm.CoordinateReference,
      DimensionCoordinate=cfdm.DimensionCoordinate,
      DomainAncillary=cfdm.DomainAncillary,
      DomainAxis=cfdm.DomainAxis,
      FieldAncillary=cfdm.FieldAncillary,
      
      Bounds=cfdm.Bounds,
      InteriorRing=cfdm.InteriorRing,
      CoordinateConversion=cfdm.CoordinateConversion,
      Datum=cfdm.Datum,
      
      List=cfdm.List,
      Index=cfdm.Index,
      Count=cfdm.Count,
      NodeCountProperties=cfdm.NodeCountProperties,
      PartNodeCountProperties=cfdm.PartNodeCountProperties,
      
      Data=cfdm.Data,
      GatheredArray=cfdm.GatheredArray,
      NetCDFArray=cfdm.NetCDFArray,
      RaggedContiguousArray=cfdm.RaggedContiguousArray,
      RaggedIndexedArray=cfdm.RaggedIndexedArray,
      RaggedIndexedContiguousArray=cfdm.RaggedIndexedContiguousArray,
   )

As all classes are required for the initialisation of the new
implementation class, this demonstrates explicitly that, in the absence
of subclasses of the other classes, the cfdm classes may be used.

.. code-block:: python
   :caption: *Read the file into 'my_Field_2' instances.*
	 
   >>> my_read_2 = functools.partial(cfdm.read,
   ...                               _implementation=my_implementation2)
   >>> q, t = my_read_2('file.nc')
   >>> print(repr(q))
   <my_Field_2: specific_humidity(latitude(5), longitude(8)) 1>

Finally, the mapping of CF data model constructs from CF-netCDF
elements, and vice versa, may be modified where desired, leaving all
other aspects it unchanged

.. code-block:: python
   :caption: *Modify the mapping of netCDF elements to CF data model
              instances.*
	 
   class my_NetCDFRead(cfdm.read_write.netcdf.NetCDFRead):
       def read(self, filename):
           """Read my fields from a netCDF file on disk or from
           an OPeNDAP server location, using my modified mapping
           from netCDF to the CF data model.
           """
           print("Reading dataset using my modified mapping")
           return super().read(filename)

.. code-block:: python
   :caption: *Create a new read-from-disk function that uses the
              modified mapping.*
	     
   my_netcdf = my_NetCDFRead(my_implementation_2)
   def my_read_3(filename, ):
       """Read my field constructs from a dataset."""
       return my_netcdf.read(filename)
   
.. code-block:: python
   :caption: *Read the file from disk into 'my_Field_2' instances,
              demonstrating that the modified mapping is being used.*
	   
   >>> q, t = my_read_3('~/cfdm/docs/source/sample_files/file.nc')
   Reading dataset using my modified mapping
   >>> print(repr(q))
   <my_Field_2: specific_humidity(latitude(5), longitude(8)) 1>

In the same manner, `cfdm.read_write.netcdf.NetCDFWrite` may be
subclassed, and a new write-to-disk function defined, to override
aspects of the mapping from CF data model constructs to netCDF
elements in a dataset.

**The _custom dictionary**
--------------------------

All cfdm classes have a `_custom` attribute that contains a dictionary
meant for use in external subclasses.

It is intended for the storage of extra objects that are required by
an external subclass, yet can be transferred to copied instances using
the inherited cfdm infrastructure. The `_custom` dictionary is shallow
copied, rather than deep copied, when using the standard cfdm deep
copy method techniques (i.e. the `!copy` method, initialisation with
the *source* parameter, or applying the `copy.deepcopy` function) so
that subclasses of cfdm are not committed to potentially expensive
deep copies of the dictionary values, of which cfdm has no
knowledge. Note that calling `copy.deepcopy` on a cfdm (sub)class
simply invokes its `!copy` method. The cfdm library itself does not
use the `_custom` dictionary, other than to pass on a shallow copy of
it to copied instances.

The consequence of this shallow-copy behaviour is that if an external
subclass stores a mutable object within its custom dictionary then, by
default, a deep copy will contain the identical mutable object, to
which in-place changes will affect both the original and copied
instances.

To account for this, the external subclass can either simply commit to
never updating such mutables in-place (which is can be acceptable for
private quantities which are tightly controlled); or else include
extra code that does deep copy such mutables when any deep copy (or
equivalent) operation is called. The latter approach should be
implemented in the subclass's `__init__` method, similarly to this:

.. code-block:: python
   :caption: *If desired, ensure that the _custom dictionary 'x'
             value is deep copied when a deep copy of an instance is
             requested.*

   import copy
   
   class my_Field_3(cfdm.Field):
       def __init__(self, properties=None, source=None, copy=True,
                    _use_data=True):
           super().__init__(properties=properties, source=source,
                            copy=copy, _use_data=_use_data)
           if source and copy:
   	       # Deep copy the custom 'x' value
               try:
    	           self._custom['x'] = copy.deepcopy(source._custom['x'])
               except (AttributeError, KeyError):
                   pass  
	   
**Documentation**
-----------------

The cfdm package uses a "docstring rewriter" that allows commonly used
parts of class and class method docstrings to be written once in a
central location, and then inserted into each class at import time. In
addition, parts of a docstring are modified to reflect the appropriate
package and class names. This functionality extends to subclasses of
cfdm classes. New docstring substitutions may also be defined for the
subclasses.

See `cfdm.core.meta.DocstringRewriteMeta` for details on how to add to
create new docstring substitutions for extensions, and how to modify
the substitutions defined in the cfdm package.


**A complete example**
----------------------

See `cf-python <https://github.com/NCAS-CMS/cf-python>`_ for a
complete example of extending the cfdm package in the manner described
above.

cf-python adds more flexible inspection, reading and writing; and
provides metadata-aware analytical processing capabilities such as
regridding and statistical calculations.

It also has a more sophisticated data class that subclasses
`cfdm.Data`, but allows for larger-than-memory manipulations and
parallel processing.

cf-python strictly extends the cfdm API, so that a cfdm command will
always work on its cf-python counterpart.
.. currentmodule:: cfdm
.. default-role:: obj

.. em dash, trimming surrounding whitespace
.. |---| unicode:: U+2014  
   :trim:
      
.. _philosophy:

**Philosophy**
==============

----

Version |release| for version |version| of the CF conventions.

.. contents::
   :local:
   :backlinks: entry

**Two levels of implementation**
--------------------------------
	       
The basic requirement of the reference implementation is to represent
the logical :ref:`CF data model <CF-data-model>` in memory with a
package of Python classes, with no further features. However, in order
to be useful the implementation must also have the practical
functionality to read and write netCDF datasets, and inspect :ref:`CF
data model constructs <CF-data-model>`.

In order to satisfy both needs there is a stand-alone core
implementation, the :ref:`cfdm.core <class_core>` package, that
includes no functionality beyond that mandated by the CF data model
(and therefore excludes any information about the netCDF encoding of
constructs). The core implementation provides the basis (via
inheritance) for the :ref:`cfdm <class_extended>` package that allows,
in addition, the reading and writing of netCDF datasets, as well as
comprehensive inspection capabilities and extra field and domain
construct modification capabilities.

----

.. _CF-conventions:

**CF conventions**
------------------

The CF data model does not enforce the CF conventions. CF-compliance
is the responsibility of the user. For example, a "units" property
whose value is not a valid `UDUNITS
<https://www.unidata.ucar.edu/software/udunits>`_ string is not
CF-compliant, but is allowed by the CF data model. This is also true,
in general, for the cfdm package. The few exceptions to this occur
when field and domain constructs are read from, or written to, a
netCDF file: it may not be possible to parse a non-CF-compliant netCDF
variable or attribute to create an unambiguous CF data model
construct; or create an unambiguous netCDF variable or attribute from
a non-CF-compliant CF data model construct.

----

**Functionality**
-----------------

The cfdm package has, with few exceptions, only the functionality
required to read and write datasets, and to create, modify and inspect
field and domain constructs in memory.

The cfdm package is not, and is not meant to be, a general analysis
package. Therefore it can't, for example, regrid field constructs to
new domains, perform statistical collapses, combine field constructs
arithmetically, etc. It has, however, been designed to be
:ref:`extensible <Extensions>` to facilitate the creation of other
packages that build on this cfdm implementation whilst also adding
extra, higher level functionality.

The `cf-python <https://ncas-cms.github.io/cf-python>`_ and `cf-plot
<http://ajheaps.github.io/cf-plot/>`_ packages, that are built on top
of the cfdm package, include much more higher level functionality.

----

**API**
-------

The design of an application programming interface (API) needs to
strike a balance between being verbose and terse. A verbose API is
easier to understand and is more memorable, but usually involves more
typing; whilst a terse API is more efficient for the experienced
user. The cfdm package has aimed for an API that is more at the
verbose end of the spectrum: in general it does not use abbreviations
for method and parameter names, and each method performs a sole
function.
.. currentmodule:: cfdm
.. default-role:: obj

.. _Support:

**Support**
===========

----

Version |release| for version |version| of the CF conventions.

General questions, suggestions for enhancements, and reports of bugs
may be reported in the GitHub issue tracker:
https://github.com/NCAS-CMS/cfdm/issues

----

.. _Contributing:

**Contributing**
----------------

Contributions to cfdm are welcome and should be raised as issues in
the GitHub issue tracker: https://github.com/NCAS-CMS/cfdm/issues,
prior to submitting a pull request containing the new code made from a
fork of the cfdm GitHub code repository:
https://github.com/NCAS-CMS/cfdm.
.. currentmodule:: cfdm
.. default-role:: obj

.. _Installation:

**Installation**
================

----

Version |release| for version |version| of the CF conventions.

.. contents::
   :local:
   :backlinks: entry

.. note:: The latest version to be released and the newest versions
          available from the Python package index (PyPI) and conda are
          confirmed at `the top of the README document
          <https://github.com/NCAS-CMS/cfdm#cfdm>`_.

.. _Python-versions:

**Operating systems**
---------------------

cfdm works for Linux, Mac and Windows operating systems.

**Python versions**
-------------------

cfdm works for Python versions 3.7 or newer.

----

.. _pip:
  
**pip**
-------

To install cfdm and all of its :ref:`dependencies <Dependencies>` run,
for example:

.. code-block:: console
   :caption: *Install as root, with any missing dependencies.*
	     
   $ pip install cfdm

.. code-block:: console
   :caption: *Install as a user, with any missing dependencies.*
	     
   $ pip install cfdm --user

To install cfdm without any of its dependencies then run, for example:

.. code-block:: console
   :caption: *Install as root without installing any of the
             dependencies.*
	     
   $ pip install cfdm --no-deps

See the `documentation for pip install
<https://pip.pypa.io/en/stable/reference/pip_install/>`_ for further
options.

----

.. _Source:

**Source**
----------

To install from source:

1. Download the cfdm package from https://pypi.org/project/cfdm

2. Unpack the library (replacing ``<version>`` with the version that
   you want to install, e.g. ``1.9.0.0``):

   .. code:: console
	 
      $ tar zxvf cfdm-<version>.tar.gz
      $ cd cfdm-<version>

3. Install the package:
  
  * To install the cfdm package to a central location:

    .. code:: console
	 
       $ python setup.py install

  * To install the cfdm package locally to the user in the default
    location:

    .. code:: console

       $ python setup.py install --user

  * To install the cfdm package in the ``<directory>`` of your choice:

    .. code:: console

       $ python setup.py install --home=<directory>

----

.. _cfdump-utility:

**cfdump utility**
------------------

During installation the :ref:`cfdump command line utility <cfdump>` is
also installed, which generates text descriptions of the :term:`field
constructs <field construct>` contained in a netCDF dataset.

----

.. _Tests:

**Tests**
---------

Tests are run from within the ``cfdm/test`` directory:

.. code:: console
 
   $ python run_tests.py
       
----

.. _Dependencies:

**Dependencies**
----------------

The cfdm package requires:

* `Python <https://www.python.org/>`_, version 3.7 or newer,

* `numpy <http://www.numpy.org/>`_, version 1.15 or newer,

* `netCDF4 <https://pypi.org/project/netCDF4/>`_, version 1.5.4 or
  newer,

* `cftime <https://pypi.org/project/cftime/>`_, version 1.5.0 or
  newer,

* `netcdf_flattener <https://pypi.org/project/netcdf-flattener/>`_,
  version 1.2.0 or newer.

----

.. _Code-repository:

**Code repository**
-------------------

The source code is available at https://github.com/NCAS-CMS/cfdm

.. .. rubric:: Footnotes

   .. [#installfiles] The ``requirements.txt`` file contains

     .. include:: ../../requirements.txt
        :literal:
     
.. currentmodule:: cfdm
.. default-role:: obj
		  
**cfdm constants**
==================

----

Version |release| for version |version| of the CF conventions.

.. data:: cfdm.masked

    A constant that allows data values to be masked by direct
    assignment. This is consistent with the :ref:`behaviour of numpy
    masked arrays <numpy:maskedarray.generic.constructing>`.

    **Examples:**

    Masking every element of a field construct's data could be done as
    follows:

    >>> f[...] = cfdm.masked
    
.. currentmodule:: cfdm
.. default-role:: obj

.. _function:

**cfdm functions**
==================

----

Version |release| for version |version| of the CF conventions.

Reading and writing
-------------------

.. autosummary::
   :nosignatures:
   :toctree: function/
   :template: function.rst

   cfdm.read 
   cfdm.write

Constants
---------

.. autosummary::
   :nosignatures:
   :toctree: function/
   :template: function.rst

   cfdm.CF
   cfdm.atol
   cfdm.rtol
   cfdm.log_level
   cfdm.configuration
   cfdm.ATOL
   cfdm.RTOL
   cfdm.LOG_LEVEL

Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: function/
   :template: function.rst

   cfdm.abspath
   cfdm.environment
   cfdm.example_field
   cfdm.example_fields
   cfdm.example_domain
   cfdm.implementation
   cfdm.unique_constructs
.. Currentmodule:: cfdm
.. default-role:: obj

.. _Tutorial:


**Tutorial**
============

----

Version |release| for version |version| of the CF conventions.

All of the Python code in this tutorial is available in an executable
script (:download:`download <../source/tutorial.py>`, 12kB).

.. contents::
   :local:
   :backlinks: entry

.. include:: sample_datasets.rst

.. _Import:

**Import**
----------

The cfdm package is imported as follows:

.. code-block:: python
   :caption: *Import the cfdm package.*

   >>> import cfdm

.. tip:: It is possible to change the extent to which cfdm outputs
         feedback messages and it may be instructive to increase the
         verbosity whilst working through this tutorial to see and
         learn more about what cfdm is doing under the hood and about
         the nature of the dataset being operated on. This can be done
         for example by running:

         .. code-block:: python
            :caption: *Increase the verbosity of cfdm from the
                      default.*

            >>> cfdm.log_level('INFO')

         See :ref:`the section on 'Logging' <Logging>` for more
         information.

.. _CF-version:

CF version
^^^^^^^^^^

The version of the `CF conventions <http://cfconventions.org>`_ and
the :ref:`CF data model <CF-data-model>` being used may be found with
the `cfdm.CF` function:

.. code-block:: python
   :caption: *Retrieve the version of the CF conventions.*
      
   >>> cfdm.CF()
   '1.9'

This indicates which version of the CF conventions are represented by
this release of the cfdm package, and therefore the version can not be
changed.

Note, however, that datasets of different versions may also be
:ref:`read <Reading-datasets>` from, or :ref:`written
<Writing-to-disk>` to, disk.

----

.. _Field-and-domain-constructs:

**Field and domain constructs**
-------------------------------

The central constructs of CF are the :term:`field construct` and
:term:`domain construct`.

The field construct, that corresponds to a CF-netCDF data variable,
includes all of the metadata to describe it:

    * descriptive properties that apply to field construct as a whole
      (e.g. the standard name),
    * a data array, and
    * "metadata constructs" that describe the locations of each cell
      (i.e. the "domain") of the data array, and the physical nature
      of each cell's datum.

Likewise, the domain construct, that corresponds to a CF-netCDF domain
variable or to the domain of a field construct, includes all of the
metadata to describe it:

    * descriptive properties that apply to field construct as a whole
      (e.g. the long name), and
    * metadata constructs that describe the locations of each cell of
      the domain.

A field construct or domain construct is stored in a `cfdm.Field`
instance or `cfdm.Domain` instance respectively. Henceforth the phrase
"field construct" will be assumed to mean "`cfdm.Field` instance", and
"domain construct" will be assumed to mean "`cfdm.Domain` instance".

----

.. _Reading-datasets:

**Reading field or domain constructs from datasets**
----------------------------------------------------

The `cfdm.read` function reads a `netCDF
<https://www.unidata.ucar.edu/software/netcdf/>`_ file from disk, or
from an `OPeNDAP <https://www.opendap.org/>`_ URL [#dap]_, and by
default returns the contents as a Python list of zero or more field
constructs. This list contains a field construct to represent each of
the CF-netCDF data variables in the file.

Datasets of any version of CF up to and including CF-|version| can be
read [#caveat]_.

All formats of netCDF3 and netCDF4 files can be read.

The file name may describe relative paths, and standard tilde and
shell parameter expansions are applied to it.

The following file types can be read:

* All formats of netCDF3 and netCDF4 files can be read, containing
  datasets for versions of CF up to and including CF-|version|.

..

* Files in `CDL format
  <https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_working_with_netcdf_files.html#netcdf_utilities>`_,
  with or without the data array values.

For example, to read the file ``file.nc`` (found in the :ref:`sample
datasets <Sample-datasets>`), which contains two field constructs:

.. code-block:: python
   :caption: *Read file.nc and show that the result is a two-element
             list.*
		
   >>> x = cfdm.read('file.nc')
   >>> print(type(x))
   <type 'list'>
   >>> len(x)
   2

Descriptive properties are always read into memory, but `lazy loading
<https://en.wikipedia.org/wiki/Lazy_loading>`_ is employed for all
data arrays, which means that no data is read into memory until the
data is required for inspection or to modify the array contents. This
maximises the number of field constructs that may be read within a
session, and makes the read operation fast.

Note that when reading netCDF4 files that contain :ref:`hierachical
groups <Hierarchical-groups>`, the group structure is saved via the
:ref:`netCDF interface <NetCDF-interface>` so that it may be re-used,
or modified, if the field constructs are written to back to disk.

The `cfdm.read` function has optional parameters to

* allow the user to provide files that contain :ref:`external
  variables <External-variables>`;

* request :ref:`extra field constructs to be created from "metadata"
  netCDF variables <Creation-by-reading>`, i.e. those that are
  referenced from CF-netCDF data variables, but which are not regarded
  by default as data variables in their own right;

* return only domain constructs derived from CF-netCDF domain
  variables;

* request that masking is *not* applied by convention to data elements
  (see :ref:`data masking <Data-mask>`);

* issue warnings when ``valid_min``, ``valid_max`` and ``valid_range``
  attributes are present (see :ref:`data masking <Data-mask>`); and

* display information and issue warnings about the mapping of the
  netCDF file contents to CF data model constructs.

.. _CF-compliance:

CF-compliance
^^^^^^^^^^^^^
  
If the dataset is partially CF-compliant to the extent that it is not
possible to unambiguously map an element of the netCDF dataset to an
element of the CF data model, then a field construct is still
returned, but may be incomplete. This is so that datasets which are
partially conformant may nonetheless be modified in memory and written
to new datasets. Such "structural" non-compliance would occur, for
example, if the ``coordinates`` attribute of a CF-netCDF data variable
refers to another variable that does not exist, or refers to a
variable that spans a netCDF dimension that does not apply to the data
variable. Other types of non-compliance are not checked, such whether
or not controlled vocabularies have been adhered to. The structural
compliance of the dataset may be checked with the
`~cfdm.Field.dataset_compliance` method of the field construct, as
well as optionally displayed when the dataset is read.

----

.. _Inspection:

**Inspection**
--------------

The contents of a field construct may be inspected at three different
levels of detail.

.. _Minimal-detail:

Minimal detail
^^^^^^^^^^^^^^

The built-in `repr` function returns a short, one-line description:

.. code-block:: python
   :caption: *Inspect the contents of the two field constructs from
             the dataset and create a Python variable for each of
             them.*
      
   >>> x
   [<Field: specific_humidity(latitude(5), longitude(8)) 1>,
    <Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>]
   >>> q = x[0]
   >>> t = x[1]
   >>> q
   <Field: specific_humidity(latitude(5), longitude(8)) 1>
   
This gives the identity of the field construct
(e.g. "specific_humidity"), the identities and sizes of the dimensions
spanned by the data array ("latitude" and "longitude" with sizes 5 and
8 respectively) and the units of the data ("1").

.. _Medium-detail:

Medium detail
^^^^^^^^^^^^^

The built-in `str` function returns similar information as the
one-line output, along with short descriptions of the metadata
constructs, which include the first and last values of their data
arrays:

.. code-block:: python
   :caption: *Inspect the contents of the two field constructs with
             medium detail.*
   
   >>> print(q)
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: time(1) = [2019-01-01 00:00:00]
                   : latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
      
   >>> print(t)
   Field: air_temperature (ncvar%ta)
   ---------------------------------
   Data            : air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K
   Cell methods    : grid_latitude(10): grid_longitude(9): mean where land (interval: 0.1 degrees) time(1): maximum
   Field ancils    : air_temperature standard_error(grid_latitude(10), grid_longitude(9)) = [[0.81, ..., 0.78]] K
   Dimension coords: time(1) = [2019-01-01 00:00:00]
                   : atmosphere_hybrid_height_coordinate(1) = [1.5]
                   : grid_latitude(10) = [2.2, ..., -1.76] degrees
                   : grid_longitude(9) = [-4.7, ..., -1.18] degrees
   Auxiliary coords: latitude(grid_latitude(10), grid_longitude(9)) = [[53.941, ..., 50.225]] degrees_N
                   : longitude(grid_longitude(9), grid_latitude(10)) = [[2.004, ..., 8.156]] degrees_E
                   : long_name=Grid latitude name(grid_latitude(10)) = [--, ..., kappa]
   Cell measures   : measure:area(grid_longitude(9), grid_latitude(10)) = [[2391.9657, ..., 2392.6009]] km2
   Coord references: atmosphere_hybrid_height_coordinate
                   : rotated_latitude_longitude
   Domain ancils   : ncvar%a(atmosphere_hybrid_height_coordinate(1)) = [10.0] m
                   : ncvar%b(atmosphere_hybrid_height_coordinate(1)) = [20.0]
                   : surface_altitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 270.0]] m

Note that :ref:`time values <Time>` are converted to date-times with
the `cftime package <https://unidata.github.io/cftime/>`_.
		   
.. _Full-detail:

Full detail
^^^^^^^^^^^

The `~cfdm.Field.dump` method of the field construct gives all
properties of all constructs, including metadata constructs and their
components, and shows the first and last values of all data arrays:

.. code-block:: python
   :caption: *Inspect the contents of the two field constructs with
             full detail.*

   >>> q.dump()
   ----------------------------------
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Conventions = 'CF-1.9'
   project = 'research'
   standard_name = 'specific_humidity'
   units = '1'
   
   Data(latitude(5), longitude(8)) = [[0.003, ..., 0.032]] 1
   
   Cell Method: area: mean
   
   Domain Axis: latitude(5)
   Domain Axis: longitude(8)
   Domain Axis: time(1)
   
   Dimension coordinate: latitude
       standard_name = 'latitude'
       units = 'degrees_north'
       Data(latitude(5)) = [-75.0, ..., 75.0] degrees_north
       Bounds:Data(latitude(5), 2) = [[-90.0, ..., 90.0]]
   
   Dimension coordinate: longitude
       standard_name = 'longitude'
       units = 'degrees_east'
       Data(longitude(8)) = [22.5, ..., 337.5] degrees_east
       Bounds:Data(longitude(8), 2) = [[0.0, ..., 360.0]]
   
   Dimension coordinate: time
       standard_name = 'time'
       units = 'days since 2018-12-01'
       Data(time(1)) = [2019-01-01 00:00:00]
  
   >>> t.dump()
   ---------------------------------
   Field: air_temperature (ncvar%ta)
   ---------------------------------
   Conventions = 'CF-1.9'
   project = 'research'
   standard_name = 'air_temperature'
   units = 'K'
   
   Data(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) = [[[0.0, ..., 89.0]]] K
   
   Cell Method: grid_latitude(10): grid_longitude(9): mean where land (interval: 0.1 degrees)
   Cell Method: time(1): maximum
   
   Field Ancillary: air_temperature standard_error
       standard_name = 'air_temperature standard_error'
       units = 'K'
       Data(grid_latitude(10), grid_longitude(9)) = [[0.81, ..., 0.78]] K
   
   Domain Axis: atmosphere_hybrid_height_coordinate(1)
   Domain Axis: grid_latitude(10)
   Domain Axis: grid_longitude(9)
   Domain Axis: time(1)
   
   Dimension coordinate: atmosphere_hybrid_height_coordinate
       computed_standard_name = 'altitude'
       standard_name = 'atmosphere_hybrid_height_coordinate'
       Data(atmosphere_hybrid_height_coordinate(1)) = [1.5]
       Bounds:Data(atmosphere_hybrid_height_coordinate(1), 2) = [[1.0, 2.0]]
   
   Dimension coordinate: grid_latitude
       standard_name = 'grid_latitude'
       units = 'degrees'
       Data(grid_latitude(10)) = [2.2, ..., -1.76] degrees
       Bounds:Data(grid_latitude(10), 2) = [[2.42, ..., -1.98]]
   
   Dimension coordinate: grid_longitude
       standard_name = 'grid_longitude'
       units = 'degrees'
       Data(grid_longitude(9)) = [-4.7, ..., -1.18] degrees
       Bounds:Data(grid_longitude(9), 2) = [[-4.92, ..., -0.96]]
   
   Dimension coordinate: time
       standard_name = 'time'
       units = 'days since 2018-12-01'
       Data(time(1)) = [2019-01-01 00:00:00]
   
   Auxiliary coordinate: latitude
       standard_name = 'latitude'
       units = 'degrees_N'
       Data(grid_latitude(10), grid_longitude(9)) = [[53.941, ..., 50.225]] degrees_N
   
   Auxiliary coordinate: longitude
       standard_name = 'longitude'
       units = 'degrees_E'
       Data(grid_longitude(9), grid_latitude(10)) = [[2.004, ..., 8.156]] degrees_E
   
   Auxiliary coordinate: long_name=Grid latitude name
       long_name = 'Grid latitude name'
       Data(grid_latitude(10)) = [--, ..., kappa]
   
   Domain ancillary: ncvar%a
       units = 'm'
       Data(atmosphere_hybrid_height_coordinate(1)) = [10.0] m
       Bounds:Data(atmosphere_hybrid_height_coordinate(1), 2) = [[5.0, 15.0]]
   
   Domain ancillary: ncvar%b
       Data(atmosphere_hybrid_height_coordinate(1)) = [20.0]
       Bounds:Data(atmosphere_hybrid_height_coordinate(1), 2) = [[14.0, 26.0]]
   
   Domain ancillary: surface_altitude
       standard_name = 'surface_altitude'
       units = 'm'
       Data(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 270.0]] m
   
   Coordinate reference: atmosphere_hybrid_height_coordinate
       Coordinate conversion:computed_standard_name = altitude
       Coordinate conversion:standard_name = atmosphere_hybrid_height_coordinate
       Coordinate conversion:a = Domain Ancillary: ncvar%a
       Coordinate conversion:b = Domain Ancillary: ncvar%b
       Coordinate conversion:orog = Domain Ancillary: surface_altitude
       Datum:earth_radius = 6371007
       Dimension Coordinate: atmosphere_hybrid_height_coordinate
   
   Coordinate reference: rotated_latitude_longitude
       Coordinate conversion:grid_mapping_name = rotated_latitude_longitude
       Coordinate conversion:grid_north_pole_latitude = 38.0
       Coordinate conversion:grid_north_pole_longitude = 190.0
       Datum:earth_radius = 6371007
       Dimension Coordinate: grid_longitude
       Dimension Coordinate: grid_latitude
       Auxiliary Coordinate: longitude
       Auxiliary Coordinate: latitude
   
   Cell measure: measure:area
       units = 'km2'
       Data(grid_longitude(9), grid_latitude(10)) = [[2391.9657, ..., 2392.6009]] km2

  
.. _cfdump:
       
cfdump
^^^^^^

The description for every field construct in a file can also be
generated from the command line, with minimal, medium or full detail,
by using the ``cfdump`` tool, for example:

.. code-block:: console
   :caption: *Use cfdump on the command line to inspect the field
             constructs contained in a dataset. The "-s" option
             requests short, minimal detail as output.*

   $ cfdump
   USAGE: cfdump [-s] [-c] [-e file [-e file] ...] [-h] file
     [-s]      Display short, one-line descriptions
     [-c]      Display complete descriptions
     [-e file] External files
     [-h]      Display the full man page
     file      Name of netCDF file (or URL if DAP access enabled)
   $ cfdump -s file.nc
   Field: specific_humidity(latitude(5), longitude(8)) 1
   Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K

``cfdump`` may also be used with :ref:`external files
<External-variables-with-cfdump>`.

----

.. _Properties:

**Properties**
--------------

Descriptive properties that apply to field construct as a whole may be
retrieved with the `~Field.properties` method:

.. code-block:: python
   :caption: *Retrieve all of the descriptive properties*
	     
   >>> t.properties()
   {'Conventions': 'CF-1.9',
    'project': 'research',
    'standard_name': 'air_temperature',
    'units': 'K'}
   
Individual properties may be accessed and modified with the
`~Field.del_property`, `~Field.get_property`, `~Field.has_property`,
and `~Field.set_property` methods:

.. code-block:: python
   :caption: *Check is a property exists, retrieve its value, delete
             it and then set it to a new value.*
      
   >>> t.has_property('standard_name')
   True
   >>> t.get_property('standard_name')
   'air_temperature'
   >>> t.del_property('standard_name')
   'air_temperature'
   >>> t.get_property('standard_name', default='not set')
   'not set'
   >>> t.set_property('standard_name', value='air_temperature')
   >>> t.get_property('standard_name', default='not set')
   'air_temperature'

A collection of properties may be set at the same time with the
`~Field.set_properties` method of the field construct, and all
properties may be completely removed with the
`~Field.clear_properties` method.

.. code-block:: python
   :caption: *Update the properties with a collection, delete all of
             the properties, and reinstate the original properties.*
	     
   >>> original = t.properties()
   >>> original
   {'Conventions': 'CF-1.9',
    'project': 'research',
    'standard_name': 'air_temperature',
    'units': 'K'}
   >>> t.set_properties({'foo': 'bar', 'units': 'K'})
   >>> t.properties()
   {'Conventions': 'CF-1.9',
    'foo': 'bar',
    'project': 'research',
    'standard_name': 'air_temperature',
    'units': 'K'}
   >>> t.clear_properties()
    {'Conventions': 'CF-1.9',
    'foo': 'bar',
    'project': 'research',
    'standard_name': 'air_temperature',
    'units': 'K'}
   >>> t.properties()
   {}
   >>> t.set_properties(original)
   >>> t.properties()
   {'Conventions': 'CF-1.9',
    'project': 'research',
    'standard_name': 'air_temperature',
    'units': 'K'}

All of the methods related to the properties are listed :ref:`here
<Field-Properties>`.

----

.. _Metadata-constructs:

**Metadata constructs**
-----------------------

The metadata constructs describe the field construct that contains
them. Each :ref:`CF data model metadata construct <CF-data-model>` has
a corresponding cfdm class:

=====================  ==============================================================  ==============================
Class                  CF data model construct                                         Description                     
=====================  ==============================================================  ==============================
`DomainAxis`           :term:`Domain axis <domain axis constructs>`                    Independent axes of the domain
`DimensionCoordinate`  :term:`Dimension coordinate <dimension coordinate constructs>`  Domain cell locations         
`AuxiliaryCoordinate`  :term:`Auxiliary coordinate <auxiliary coordinate constructs>`  Domain cell locations         
`CoordinateReference`  :term:`Coordinate reference <coordinate reference constructs>`  Domain coordinate systems     
`DomainAncillary`      :term:`Domain ancillary <domain ancillary constructs>`          Cell locations in alternative 
                                                                                       coordinate systems	       
`CellMeasure`          :term:`Cell measure <cell measure constructs>`                  Domain cell size or shape     
`FieldAncillary`       :term:`Field ancillary <field ancillary constructs>`            Ancillary metadata which vary 
                                                                                       within the domain	       
`CellMethod`           :term:`Cell method <cell method constructs>`                    Describes how data represent  
                                                                                       variation within cells	       
=====================  ==============================================================  ==============================

Metadata constructs of a particular type can be retrieved with the
following methods of the field construct:

==============================  =====================  
Method                          Metadata constructs    
==============================  =====================  
`~Field.domain_axes`            Domain axes            
`~Field.dimension_coordinates`  Dimension coordinates  
`~Field.auxiliary_coordinates`  Auxiliary coordinates  
`~Field.coordinate_references`  Coordinate references  
`~Field.domain_ancillaries`     Domain ancillaries     
`~Field.cell_measures`          Cell measures          
`~Field.field_ancillaries`      Field ancillaries      
`~Field.cell_methods`           Cell methods                               
==============================  =====================  

Each of these attributes returns a `Constructs` class instance that
maps metadata constructs to unique identifiers called "construct
keys". A `~Constructs` instance has methods for selecting constructs
that meet particular criteria (see the section on :ref:`filtering
metadata constructs <Filtering-metadata-constructs>`). It also behaves
like a "read-only" Python dictionary, in that it has
`~Constructs.items`, `~Constructs.keys` and `~Constructs.values`
methods that work exactly like their corresponding `dict` methods. It
also has a `~Constructs.get` method and indexing like a Python
dictionary (see the section on :ref:`metadata construct access
<Metadata-construct-access>` for details).

.. Each of these methods returns a dictionary whose values are the
   metadata constructs of one type, keyed by a unique identifier
   called a "construct key":

.. code-block:: python
   :caption: *Retrieve the field construct's coordinate reference
             constructs, and access them using dictionary methods.*
      
   >>> t.coordinate_references()
   <Constructs: coordinate_reference(2)>
   >>> print(t.coordinate_references())
   Constructs:
   {'coordinatereference0': <CoordinateReference: atmosphere_hybrid_height_coordinate>,
    'coordinatereference1': <CoordinateReference: rotated_latitude_longitude>}
   >>> list(t.coordinate_references().keys())
   ['coordinatereference0', 'coordinatereference1']
   >>> for key, value in t.coordinate_references().items():
   ...     print(key, repr(value))
   ...
   coordinatereference1 <CoordinateReference: rotated_latitude_longitude>
   coordinatereference0 <CoordinateReference: atmosphere_hybrid_height_coordinate>

.. code-block:: python
   :caption: *Retrieve the field construct's dimension coordinate and
             domain axis constructs.*
      
   >>> print(t.dimension_coordinates())
   Constructs:
   {'dimensioncoordinate0': <DimensionCoordinate: atmosphere_hybrid_height_coordinate(1) >,
    'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>,
    'dimensioncoordinate3': <DimensionCoordinate: time(1) days since 2018-12-01 >}
   >>> print(t.domain_axes())
   Constructs:
   {'domainaxis0': <DomainAxis: size(1)>,
    'domainaxis1': <DomainAxis: size(10)>,
    'domainaxis2': <DomainAxis: size(9)>,
    'domainaxis3': <DomainAxis: size(1)>}

The construct keys (e.g. ``'domainaxis1'``) are usually generated
internally and are unique within the field construct. However,
construct keys may be different for equivalent metadata constructs
from different field constructs, and for different Python sessions.

Metadata constructs of all types may be returned by the
`~Field.constructs` attribute of the field construct:

.. code-block:: python
   :caption: *Retrieve all of the field construct's metadata
             constructs.*

   >>> q.constructs
   <Constructs: cell_method(1), dimension_coordinate(3), domain_axis(3)>
   >>> print(q.constructs)
   Constructs:
   {'cellmethod0': <CellMethod: area: mean>,
    'dimensioncoordinate0': <DimensionCoordinate: latitude(5) degrees_north>,
    'dimensioncoordinate1': <DimensionCoordinate: longitude(8) degrees_east>,
    'dimensioncoordinate2': <DimensionCoordinate: time(1) days since 2018-12-01 >,
    'domainaxis0': <DomainAxis: size(5)>,
    'domainaxis1': <DomainAxis: size(8)>,
    'domainaxis2': <DomainAxis: size(1)>}
   >>> t.constructs
   <Constructs: auxiliary_coordinate(3), cell_measure(1), cell_method(2), coordinate_reference(2), dimension_coordinate(4), domain_ancillary(3), domain_axis(4), field_ancillary(1)>
   >>> print(t.constructs)
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>,
    'auxiliarycoordinate1': <AuxiliaryCoordinate: longitude(9, 10) degrees_E>,
    'auxiliarycoordinate2': <AuxiliaryCoordinate: long_name=Grid latitude name(10) >,
    'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>,
    'cellmethod0': <CellMethod: domainaxis1: domainaxis2: mean where land (interval: 0.1 degrees)>,
    'cellmethod1': <CellMethod: domainaxis3: maximum>,
    'coordinatereference0': <CoordinateReference: atmosphere_hybrid_height_coordinate>,
    'coordinatereference1': <CoordinateReference: rotated_latitude_longitude>,
    'dimensioncoordinate0': <DimensionCoordinate: atmosphere_hybrid_height_coordinate(1) >,
    'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>,
    'dimensioncoordinate3': <DimensionCoordinate: time(1) days since 2018-12-01 >,
    'domainancillary0': <DomainAncillary: ncvar%a(1) m>,
    'domainancillary1': <DomainAncillary: ncvar%b(1) >,
    'domainancillary2': <DomainAncillary: surface_altitude(10, 9) m>,
    'domainaxis0': <DomainAxis: size(1)>,
    'domainaxis1': <DomainAxis: size(10)>,
    'domainaxis2': <DomainAxis: size(9)>,
    'domainaxis3': <DomainAxis: size(1)>,
    'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}

----

.. _Data:

**Data**
--------

The field construct's data is stored in a `Data` class instance that
is accessed with the `~Field.data` attribute of the field construct.

.. code-block:: python
   :caption: *Retrieve the data and inspect it, showing the shape and
             some illustrative values.*
		
   >>> t.data
   <Data(1, 10, 9): [[[262.8, ..., 269.7]]] K>

The `Data` instance provides access to the full array of values, as
well as attributes to describe the array and methods for describing
any :ref:`data compression <Compression>`.

The `Data` instance provides access to the full array of values, as
well as attributes to describe the array and methods for describing
any data compression. However, the field construct (and any other
construct that contains data) also provides attributes for direct
access.

.. code-block:: python
   :caption: *Retrieve a numpy array of the data.*
      
   >>> print(t.data.array)
   [[[262.8 270.5 279.8 269.5 260.9 265.0 263.5 278.9 269.2]
     [272.7 268.4 279.5 278.9 263.8 263.3 274.2 265.7 279.5]
     [269.7 279.1 273.4 274.2 279.6 270.2 280.0 272.5 263.7]
     [261.7 260.6 270.8 260.3 265.6 279.4 276.9 267.6 260.6]
     [264.2 275.9 262.5 264.9 264.7 270.2 270.4 268.6 275.3]
     [263.9 263.8 272.1 263.7 272.2 264.2 260.0 263.5 270.2]
     [273.8 273.1 268.5 272.3 264.3 278.7 270.6 273.0 270.6]
     [267.9 273.5 279.8 260.3 261.2 275.3 271.2 260.8 268.9]
     [270.9 278.7 273.2 261.7 271.6 265.8 273.0 278.5 266.4]
     [276.4 264.2 276.3 266.1 276.1 268.1 277.0 273.4 269.7]]]

.. code-block:: python
   :caption: *Inspect the data type, number of dimensions, dimension
             sizes and number of elements of the data.*
	     
   >>> t.dtype
   dtype('float64')
   >>> t.ndim
   3
   >>> t.shape
   (1, 10, 9)
   >>> t.size
   90
   >>> t.data.size
   90

Note it is preferable to access the data type, number of dimensions,
dimension sizes and number of elements of the data via the parent
construct, rather than from the `Data` instance, as there are
:ref:`particular circumstances <Geometry-cells>` when there is no
`Data` instance, but the construct nonetheless has data descriptors.
   
The field construct also has a `~Field.get_data` method as an
alternative means of retrieving the data instance, which allows for a
default to be returned if no data have been set; as well as a
`~Field.del_data` method for removing the data.

All of the methods and attributes related to the data are listed
:ref:`here <Field-Data>`.

.. _Data-axes:

Data axes
^^^^^^^^^

The data array of the field construct spans all the :term:`domain axis
constructs` with the possible exception of size one domain axis
constructs. The domain axis constructs spanned by the field
construct's data are found with the `~Field.get_data_axes` method of
the field construct. For example, the data of the field construct
``t`` does not span the size one domain axis construct with key
``'domainaxis3'``.

.. code-block:: python
   :caption: *Show which data axis constructs are spanned by the field
             construct's data.*
	    
   >>> print(t.domain_axes())
   Constructs:
   {'domainaxis0': <DomainAxis: size(1)>,
    'domainaxis1': <DomainAxis: size(10)>,
    'domainaxis2': <DomainAxis: size(9)>,
    'domainaxis3': <DomainAxis: size(1)>}
   >>> t
   <Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>
   >>> t.shape
   (1, 10, 9)
   >>> t.get_data_axes()
   ('domainaxis0', 'domainaxis1', 'domainaxis2')

The data may be set with the `~Field.set_data` method of the field
construct. The domain axis constructs spanned by the data may also be
set by explicitly providing them via their construct keys. In any
case, the data axes may be set at any time with the
`~Field.set_data_axes` method of the field construct.

.. code-block:: python
   :caption: *Delete the data and then reinstate it, using the
             existing data axes.*
	    
   >>> data = t.del_data()
   >>> t.has_data()
   False
   >>> t.set_data(data, axes=None)
   >>> t.data
   <Data(1, 10, 9): [[[262.8, ..., 269.7]]] K>

See the section :ref:`field construct creation
<Field-creation-in-memory>` for more examples.

.. _Date-time:

Date-time
^^^^^^^^^

Data representing date-times is defined as elapsed times since a
reference date-time in a particular calendar (Gregorian, by
default). The `~cfdm.Data.array` attribute of the `Data` instance
(and any construct that contains it) returns the elapsed times, and
the `~cfdm.Data.datetime_array` (and any construct that contains it)
returns the data as an array of date-time objects.

.. code-block:: python
   :caption: *View date-times aas elapsed time or as date-time
             objects.*
	     
   >>> d = cfdm.Data([1, 2, 3], units='days since 2004-2-28')
   >>> print(d.array)   
   [1 2 3]
   >>> print(d.datetime_array)
   [cftime.DatetimeGregorian(2004-02-29 00:00:00)
    cftime.DatetimeGregorian(2004-03-01 00:00:00)
    cftime.DatetimeGregorian(2004-03-02 00:00:00)]
   >>> e = cfdm.Data([1, 2, 3], units='days since 2004-2-28',
   ...               calendar='360_day')
   >>> print(e.array)   
   [1 2 3]
   >>> print(e.datetime_array)
   [cftime.Datetime360Day(2004-02-29 00:00:00)
    cftime.Datetime360Day(2004-02-30 00:00:00)
    cftime.Datetime360Day(2004-03-01 00:00:00)]

    
.. _Manipulating-dimensions:

Manipulating dimensions
^^^^^^^^^^^^^^^^^^^^^^^

The dimensions of a field construct's data may be reordered, have size
one dimensions removed and have new new size one dimensions included
by using the following field construct methods:

=========================  ===========================================
Method                     Description
=========================  ===========================================
`~Field.insert_dimension`  Insert a new size one data dimension. The
                           new dimension must correspond to an
                           existing size one domain axis construct.

`~Field.squeeze`           Remove size one data dimensions
	   
`~Field.transpose`         Reorder data dimensions
=========================  ===========================================

.. code-block:: python
   :caption: *Remove all size one dimensions from the data, noting
             that metadata constructs which span the corresponding
             domain axis construct are not affected.*

   >>> q, t = cfdm.read('file.nc')
   >>> t
   <CF Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>
   >>> t2 = t.squeeze()
   >>> t2
   <CF Field: air_temperature(grid_latitude(10), grid_longitude(9)) K>
   >>> print(t2.dimension_coordinates())
   Constructs:
   {'dimensioncoordinate0': <CF DimensionCoordinate: atmosphere_hybrid_height_coordinate(1) >,
    'dimensioncoordinate1': <CF DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <CF DimensionCoordinate: grid_longitude(9) degrees>,
    'dimensioncoordinate3': <CF DimensionCoordinate: time(1) days since 2018-12-01 >}

.. code-block:: python
   :caption: *Insert a new size one dimension, corresponding to a size
             one domain axis construct, and then reorder the
             dimensions.*

   >>> t3 = t2.insert_dimension(axis='domainaxis3', position=1)
   >>> t3
   <CF Field: air_temperature(grid_latitude(10), time(1), grid_longitude(9)) K>
   >>> t3.transpose([2, 0, 1])
   <CF Field: air_temperature(grid_longitude(9), grid_latitude(10), time(1)) K>

When transposing the data dimensions, the dimensions of metadata
construct data are, by default, unchanged. It is also possible to
permute the data dimensions of the metadata constructs so that they
have the same relative order as the field construct:

.. code-block:: python
   :caption: *Also permute the data dimension of metadata constructs
             using the 'constructs' keyword.*

   >>> t4 = t.transpose([0, 2, 1], constructs=True)

.. _Data-mask:
   
Data mask
^^^^^^^^^

There is always a data mask, which may be thought of as a separate
data array of Booleans with the same shape as the original data. The
data mask is `False` where the the data has values, and `True` where
the data is missing. The data mask may be inspected with the
`~cfdm.Data.mask` attribute of the data instance, which returns the
data mask in a field construct with the same metadata constructs as
the original field construct.


.. code-block:: python
   :caption: *Inspect the data mask of a field construct.*

   >>> print(q)
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   >>> print(q.data.mask)
   <Data(5, 8): [[False, ..., False]]>
   >>> print(q.data.mask.array)
   [[False False False False False False False False]
    [False False False False False False False False]
    [False False False False False False False False]
    [False False False False False False False False]
    [False False False False False False False False]]

.. code-block:: python
   :caption: *Mask the polar rows (see the "Assignment by index"
             section) and inspect the new data mask.*
	  
   >>> q.data[[0, 4], :] = cfdm.masked            
   >>> print(q.data.mask.array)
   [[ True  True  True  True  True  True  True  True]
    [False False False False False False False False]
    [False False False False False False False False]
    [False False False False False False False False]
    [ True  True  True  True  True  True  True  True]]

The ``_FillValue`` and ``missing_value`` attributes of the field
construct are *not* stored as values of the field construct's
data. They are only used when :ref:`writing the data to a netCDF
dataset <Writing-to-a-netCDF-dataset>`. Therefore testing for missing
data by testing for equality to one of these property values will
produce incorrect results; the `~Data.any` method of the `Data`
instance should be used instead.

.. code-block:: python
   :caption: *See if any data points are masked.*
	     
   >>> q.data.mask.any()
   True

The mask of a netCDF dataset array is implied by array values that
meet the criteria implied by the ``missing_value``, ``_FillValue``,
``valid_min``, ``valid_max``, and ``valid_range`` properties, and is
usually applied automatically by `cfdm.read`. NetCDF data elements
that equal the values of the ``missing_value`` and ``_FillValue``
properties are masked, as are data elements that exceed the value of
the ``valid_max`` property, subceed the value of the ``valid_min``
property, or lie outside of the range defined by the ``valid_range``
property.

However, this automatic masking may be bypassed by setting the *mask*
keyword of the `cfdm.read` function to `False`. The mask, as defined
in the dataset, may subsequently be applied manually with the
`~Field.apply_masking` method of the field construct.

.. code-block:: python
   :caption: *Read a dataset from disk without automatic masking, and
             then manually apply the mask*

   >>> cfdm.write(q, 'masked_q.nc')
   >>> no_mask_q = cfdm.read('masked_q.nc', mask=False)[0]
   >>> print(no_mask_q.data.array)
   [9.96920997e+36, 9.96920997e+36, 9.96920997e+36, 9.96920997e+36,
    9.96920997e+36, 9.96920997e+36, 9.96920997e+36, 9.96920997e+36],
    [0.023 0.036 0.045 0.062 0.046 0.073 0.006 0.066]
    [0.11  0.131 0.124 0.146 0.087 0.103 0.057 0.011]
    [0.029 0.059 0.039 0.07  0.058 0.072 0.009 0.017]
   [9.96920997e+36, 9.96920997e+36, 9.96920997e+36, 9.96920997e+36,
    9.96920997e+36, 9.96920997e+36, 9.96920997e+36, 9.96920997e+36]])
   >>> masked_q = no_mask_q.apply_masking()
   >>> print(masked_q.data.array)
   [[   --    --    --    --    --    --    --    --]
    [0.023 0.036 0.045 0.062 0.046 0.073 0.006 0.066]
    [0.11  0.131 0.124 0.146 0.087 0.103 0.057 0.011]
    [0.029 0.059 0.039 0.07  0.058 0.072 0.009 0.017]
    [   --    --    --    --    --    --    --    --]]

The `~Field.apply_masking` method of the field construct utilises as
many of the ``missing_value``, ``_FillValue``, ``valid_min``,
``valid_max``, and ``valid_range`` properties as are present and may
be used on any construct, not just those that have been read from
datasets.
    
.. _Indexing:

Indexing
^^^^^^^^

When a `Data` instance is indexed a new instance is created for the
part of the data defined by the indices. Indexing follows rules that
are very similar to the `numpy indexing rules
<https://numpy.org/doc/stable/user/basics.indexing.html>`_,
the only differences being:

* An integer index *i* specified for a dimension reduces the size of
  this dimension to unity, taking just the *i*\ -th element, but keeps
  the dimension itself, so that the rank of the array is not reduced.

..

* When two or more dimensions' indices are sequences of integers then
  these indices work independently along each dimension (similar to
  the way vector subscripts work in Fortran). This is the same
  indexing behaviour as on a ``Variable`` object of the `netCDF4
  package <http://unidata.github.io/netcdf4-python>`_.

.. code-block:: python
   :caption: *Create new data by indexing and show the shape
             corresponding to the indices.*
	     
   >>> data = t.data
   >>> data.shape
   (1, 10, 9)
   >>> data[:, :, 1].shape
   (1, 10, 1)
   >>> data[:, 0].shape
   (1, 1, 9)
   >>> data[..., 6:3:-1, 3:6].shape
   (1, 3, 3)
   >>> data[0, [2, 9], [4, 8]].shape
   (1, 2, 2)
   >>> data[0, :, -2].shape
   (1, 10, 1)

.. _Assignment:

Assignment
^^^^^^^^^^

Values can be changed by assigning to elements selected by indices of
the `Data` instance using rules that are very similar to the `numpy
indexing rules
<https://numpy.org/doc/stable/user/basics.indexing.html>`_,
the only difference being:

* When two or more dimensions' indices are sequences of integers then
  these indices work independently along each dimension (similar to
  the way vector subscripts work in Fortran). This is the same
  indexing behaviour as on a ``Variable`` object of the `netCDF4
  package <http://unidata.github.io/netcdf4-python>`_.

A single value may be assigned to any number of elements.
  
.. code-block:: python
   :caption: *Set a single element to -1, a "column" of elements
             to -2 and a "square" of elements to -3.*
	     
   >>> import numpy
   >>> t.data[:, 0, 0] = -1
   >>> t.data[:, :, 1] = -2
   >>> t.data[..., 6:3:-1, 3:6] = -3
   >>> print(t.data.array)
   [[[ -1.0  -2.0 279.8 269.5 260.9 265.0 263.5 278.9 269.2]
     [272.7  -2.0 279.5 278.9 263.8 263.3 274.2 265.7 279.5]
     [269.7  -2.0 273.4 274.2 279.6 270.2 280.0 272.5 263.7]
     [261.7  -2.0 270.8 260.3 265.6 279.4 276.9 267.6 260.6]
     [264.2  -2.0 262.5  -3.0  -3.0  -3.0 270.4 268.6 275.3]
     [263.9  -2.0 272.1  -3.0  -3.0  -3.0 260.0 263.5 270.2]
     [273.8  -2.0 268.5  -3.0  -3.0  -3.0 270.6 273.0 270.6]
     [267.9  -2.0 279.8 260.3 261.2 275.3 271.2 260.8 268.9]
     [270.9  -2.0 273.2 261.7 271.6 265.8 273.0 278.5 266.4]
     [276.4  -2.0 276.3 266.1 276.1 268.1 277.0 273.4 269.7]]]

An array of values can be assigned, as long as it is broadcastable to
the shape defined by the indices, using the `numpy broadcasting rules
<https://numpy.org/doc/stable/user/basics.broadcasting.html>`_.

.. code-block:: python
   :caption: *Overwrite the square of values of -3 with the numbers 0
             to 8, and set the corners of a different square to be
             either -4 or -5.*
	     
   >>> t.data[..., 6:3:-1, 3:6] = numpy.arange(9).reshape(3, 3)
   >>> t.data[0, [2, 9], [4, 8]] =  cfdm.Data([[-4, -5]])
   >>> print(t.data.array)
   [[[ -1.0  -2.0 279.8 269.5 260.9 265.0 263.5 278.9 269.2]
     [272.7  -2.0 279.5 278.9 263.8 263.3 274.2 265.7 279.5]
     [269.7  -2.0 273.4 274.2  -4.0 270.2 280.0 272.5  -5.0]
     [261.7  -2.0 270.8 260.3 265.6 279.4 276.9 267.6 260.6]
     [264.2  -2.0 262.5   6.0   7.0   8.0 270.4 268.6 275.3]
     [263.9  -2.0 272.1   3.0   4.0   5.0 260.0 263.5 270.2]
     [273.8  -2.0 268.5   0.0   1.0   2.0 270.6 273.0 270.6]
     [267.9  -2.0 279.8 260.3 261.2 275.3 271.2 260.8 268.9]
     [270.9  -2.0 273.2 261.7 271.6 265.8 273.0 278.5 266.4]
     [276.4  -2.0 276.3 266.1  -4.0 268.1 277.0 273.4  -5.0]]]

Data array elements may be set to missing values by assigning them to
the `cfdm.masked` constant. Missing values may be unmasked by
assigning them to any other value.

.. code-block:: python
   :caption: *Set a column of elements to missing values, and then
             change one of them back to a non-missing value.*
	     
   >>> t.data[0, :, -2] = cfdm.masked
   >>> t.data[0, 5, -2] = -6
   >>> print(t.data.array)
   [[[ -1.0  -2.0 279.8 269.5 260.9 265.0 263.5    -- 269.2]
     [272.7  -2.0 279.5 278.9 263.8 263.3 274.2    -- 279.5]
     [269.7  -2.0 273.4 274.2  -4.0 270.2 280.0    --  -5.0]
     [261.7  -2.0 270.8 260.3 265.6 279.4 276.9    -- 260.6]
     [264.2  -2.0 262.5   6.0   7.0   8.0 270.4    -- 275.3]
     [263.9  -2.0 272.1   3.0   4.0   5.0 260.0  -6.0 270.2]
     [273.8  -2.0 268.5   0.0   1.0   2.0 270.6    -- 270.6]
     [267.9  -2.0 279.8 260.3 261.2 275.3 271.2    -- 268.9]
     [270.9  -2.0 273.2 261.7 271.6 265.8 273.0    -- 266.4]
     [276.4  -2.0 276.3 266.1  -4.0 268.1 277.0    --  -5.0]]]

----

.. _Subspacing:

**Subspacing**
--------------

Creation of a new field construct which spans a subspace of the domain
of an existing field construct is achieved by indexing the field
itself, rather than its `Data` instance. This is because the operation
must also subspace any metadata constructs of the field construct
(e.g. coordinate metadata constructs) which span any of the
:term:`domain axis constructs` that are affected. The new field
construct is created with the same properties as the original
field. Subspacing uses the same :ref:`cfdm indexing rules <Indexing>`
that apply to the `Data` class.

.. code-block:: python
   :caption: *Create a new field whose domain spans the first longitude
             of the original, and with a reversed latitude axis.*

   >>> print(q)
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: time(1) = [2019-01-01 00:00:00]
                   : latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east

   >>> new = q[::-1, 0]
   >>> print(new)
   Field: specific_humidity (ncvar%q)
   ----------------------------------
   Data            : specific_humidity(latitude(5), longitude(1)) 1
   Cell methods    : area: mean
   Dimension coords: time(1) = [2019-01-01 00:00:00]
                   : latitude(5) = [75.0, ..., -75.0] degrees_north
                   : longitude(1) = [22.5] degrees_east

----

.. _Selecting-metadata-constructs:

**Selecting metadata constructs**
---------------------------------

A `Constructs` instance has filtering methods for selecting constructs
that meet various criteria:

================================  ==========================================================================  
Method                            Filter criteria                                                             
================================  ==========================================================================  
`~Constructs.filter`              General purpose interface to all other filter methods
`~Constructs.filter_by_identity`  Metadata construct identity                
`~Constructs.filter_by_type`      Metadata construct type                       
`~Constructs.filter_by_property`  Property values                                     
`~Constructs.filter_by_axis`      The :term:`domain axis constructs` spanned by the data
`~Constructs.filter_by_naxes`     The number of :term:`domain axis constructs` spanned by the data
`~Constructs.filter_by_size`      The size :term:`domain axis constructs`
`~Constructs.filter_by_measure`   Measure value (for cell measure constructs)
`~Constructs.filter_by_method`    Method value (for cell method constructs)	
`~Constructs.filter_by_data`      Whether or not there could be be data.
`~Constructs.filter_by_key`       Construct key			
`~Constructs.filter_by_ncvar`     NetCDF variable name (see the :ref:`netCDF interface <NetCDF-interface>`)
`~Constructs.filter_by_ncdim`     NetCDF dimension name (see the :ref:`netCDF interface <NetCDF-interface>`)
================================  ==========================================================================  

The `~Constructs.filter` method of a `Constructs` instance allows
these filters to be chained together in a single call.

Each of these methods returns a new `Constructs` instance by default
that contains the selected constructs.

.. code-block:: python
   :caption: *Get constructs by their type*.
	  
   >>> print(t.constructs.filter_by_type('dimension_coordinate'))
   Constructs:
   {'dimensioncoordinate0': <DimensionCoordinate: atmosphere_hybrid_height_coordinate(1) >,
    'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>,
    'dimensioncoordinate3': <DimensionCoordinate: time(1) days since 2018-12-01 >}
   >>> print(t.constructs.filter_by_type('cell_method', 'field_ancillary'))
   Constructs:
   {'cellmethod0': <CellMethod: domainaxis1: domainaxis2: mean where land (interval: 0.1 degrees)>,
    'cellmethod1': <CellMethod: domainaxis3: maximum>,
    'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}

.. code-block:: python
   :caption: *Get constructs by their properties*.

   >>> print(t.constructs.filter_by_property(
   ...             standard_name='air_temperature standard_error'))
   Constructs:
   {'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}
   >>> print(t.constructs.filter_by_property(
   ...             standard_name='air_temperature standard_error',
   ...             units='K'))
   Constructs:
   {'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}
   >>> print(t.constructs.filter_by_property(
   ...             'or',
   ...	           standard_name='air_temperature standard_error',
   ...             units='m'))
   Constructs:
   {'domainancillary0': <DomainAncillary: ncvar%a(1) m>,
    'domainancillary2': <DomainAncillary: surface_altitude(10, 9) m>,
    'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}
   
.. code-block:: python
   :caption: *Get constructs whose data span at least one of the
             'grid_latitude' and 'grid_longitude' domain axis constructs.*

   >>> print(t.constructs.filter_by_axis('grid_latitude', 'grid_longitude',
   ...                                   axis_mode='or'))
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>,
    'auxiliarycoordinate1': <AuxiliaryCoordinate: longitude(9, 10) degrees_E>,
    'auxiliarycoordinate2': <AuxiliaryCoordinate: long_name=Grid latitude name(10) >,
    'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>,
    'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'domainancillary2': <DomainAncillary: surface_altitude(10, 9) m>,
    'fieldancillary0': <FieldAncillary: air_temperature standard_error(10, 9) K>}

.. code-block:: python
   :caption: *Get cell measure constructs by their "measure".*
	     
   >>> print(t.constructs.filter_by_measure('area'))
   Constructs:
   {'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>}

.. code-block:: python
   :caption: *Get cell method constructs by their "method".*
	     
   >>> print(t.constructs.filter_by_method('maximum'))
   Constructs:
   {'cellmethod1': <CellMethod: domainaxis3: maximum>}

As each of these methods returns a `Constructs` instance by default,
it is easy to perform further filters on their results:
   
.. code-block:: python
   :caption: *Make selections from previous selections.*
	     
   >>> print(
   ...     t.constructs.filter_by_type('auxiliary_coordinate').filter_by_axis('domainaxis2')
   ... )
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>,
    'auxiliarycoordinate1': <AuxiliaryCoordinate: longitude(9, 10) degrees_E>}
   >>> c = t.constructs.filter_by_type('dimension_coordinate')
   >>> d = c.filter_by_property(units='degrees')
   >>> print(d)
   Constructs:
   {'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>}

Filters can also be chained with the `~Constructs.filter` method of a
`Constructs` instance
   
.. code-block:: python
   :caption: *Make a chain of selections.*
	  
   >>> c = t.constructs.filter(filter_by_type=('dimension_coordinate',),
   ...                         filter_by_property={'units': 'degrees'})
   >>> print(c)
   Constructs:
   {'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>}

If the result are only required as a `dict`, rather than a
`Constructs` instance, the "todict" parameter can be used to give
faster performance:
    
.. code-block:: python
   :caption: *Make a chain of selections and return as a dictionary.*
	  
   >>> d = t.constructs.filter(filter_by_type=('dimension_coordinate',),
   ...                         filter_by_property={'units': 'degrees'},
   ...                         todict=True)
   >>> type(d)
   dict
   >>> print(d)
   {'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>}

Filters can also be chained with the `~Constructs.filter` method of a
`Constructs` instance
   
.. code-block:: python
   :caption: *Make a chain of selections.*
	  
   >>> c = t.constructs.filter(filter_by_type=('dimension_coordinate',),
   ...                         filter_by_property={'units': 'degrees'})
   >>> print(c)
   Constructs:
   {'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>}

Another method of selection is by metadata construct "identity".
Construct identities are used to describe constructs when they are
inspected, and so it is often convenient to copy these identities
when selecting metadata constructs. For example, the :ref:`three
auxiliary coordinate constructs <Medium-detail>` in the field
construct ``t`` have identities ``'latitude'``, ``'longitude'`` and
``'long_name=Grid latitude name'``.

A construct's identity may be any one of the following

* The value of the ``standard_name`` property,
  e.g. ``'air_temperature'``,
* The value of any property, preceded by the property name and an
  equals, e.g. ``'long_name=Air Temperature'``, ``'axis=X'``,
  ``'foo=bar'``, etc.,
* The cell measure, preceded by "measure:",
  e.g. ``'measure:volume'``
* The cell method, preceded by "method:", e.g. ``'method:maximum'``
* The netCDF variable name, preceded by "ncvar%",
  e.g. ``'ncvar%tas'`` (see the :ref:`netCDF interface
  <NetCDF-interface>`),
* The netCDF dimension name, preceded by "ncdim%" e.g. ``'ncdim%z'``
  (see the :ref:`netCDF interface <NetCDF-interface>`), and 
* The construct key, preceded by "key%"
  e.g. ``'key%auxiliarycoordinate2'``.

.. code-block:: python
   :caption: *Get constructs by their identity.*
	
   >>> print(t)
   Field: air_temperature (ncvar%ta)
   ---------------------------------
   Data            : air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K
   Cell methods    : grid_latitude(10): grid_longitude(9): mean where land (interval: 0.1 degrees) time(1): maximum
   Field ancils    : air_temperature standard_error(grid_latitude(10), grid_longitude(9)) = [[0.81, ..., 0.78]] K
   Dimension coords: time(1) = [2019-01-01 00:00:00]
                   : atmosphere_hybrid_height_coordinate(1) = [1.5]
                   : grid_latitude(10) = [2.2, ..., -1.76] degrees
                   : grid_longitude(9) = [-4.7, ..., -1.18] degrees
   Auxiliary coords: latitude(grid_latitude(10), grid_longitude(9)) = [[53.941, ..., 50.225]] degrees_N
                   : longitude(grid_longitude(9), grid_latitude(10)) = [[2.004, ..., 8.156]] degrees_E
                   : long_name=Grid latitude name(grid_latitude(10)) = [--, ..., kappa]
   Cell measures   : measure:area(grid_longitude(9), grid_latitude(10)) = [[2391.9657, ..., 2392.6009]] km2
   Coord references: atmosphere_hybrid_height_coordinate
                   : rotated_latitude_longitude
   Domain ancils   : ncvar%a(atmosphere_hybrid_height_coordinate(1)) = [10.0] m
                   : ncvar%b(atmosphere_hybrid_height_coordinate(1)) = [20.0]
                   : surface_altitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 270.0]] m
   >>> print(t.constructs.filter_by_identity('latitude'))
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>}
   >>> print(t.constructs.filter_by_identity('long_name=Grid latitude name'))
   Constructs:
   {'auxiliarycoordinate2': <AuxiliaryCoordinate: long_name=Grid latitude name(10) >}
   >>> print(t.constructs.filter_by_identity('measure:area'))
   Constructs:
   {'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>}
   >>> print(t.constructs.filter_by_identity('ncvar%b'))
   Constructs:
   {'domainancillary1': <DomainAncillary: ncvar%b(1) >}

Each construct has an `!identity` method that, by default, returns the
least ambiguous identity (defined in the documentation of a
construct's `!identity` method); and an `!identities` method that
returns a list of all of the identities that would select the
construct.

As a further convenience, selection by construct identity is also
possible by providing identities to a call of a `Constructs` instance
itself, and this technique for selecting constructs by identity will be
used in the rest of this tutorial:

.. code-block:: python
   :caption: *Construct selection by identity is possible with via the
             "filter_by_identity" method, or directly from the
             "Constructs" instance.*

   >>> print(t.constructs.filter_by_identity('latitude'))
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>}
   >>> print(t.constructs('latitude'))
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>}

Selection by construct key is useful for systematic metadata construct
access, or for when a metadata construct is not identifiable by other
means:

.. code-block:: python
   :caption: *Get constructs by construct key.*

   >>> print(t.constructs.filter_by_key('domainancillary2'))
   Constructs:
   {'domainancillary2': <DomainAncillary: surface_altitude(10, 9) m>}
   >>> print(t.constructs.filter_by_key('cellmethod1'))
   Constructs:
   {'cellmethod1': <CellMethod: domainaxis3: maximum>}
   >>> print(t.constructs.filter_by_key('auxiliarycoordinate2', 'cellmeasure0'))
   Constructs:
   {'auxiliarycoordinate2': <AuxiliaryCoordinate: long_name=Grid latitude name(10) >,
    'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>}

If no constructs match the given criteria, then an "empty"
`Constructs` instance is returned:
   
.. code-block:: python
   :caption: *If no constructs meet the criteria then an empty
             "Constructs" object is returned.*

   >>> c = t.constructs('radiation_wavelength')
   >>> c
   <Constructs: >
   >>> print(c)
   Constructs:
   {}
   >>> len(c)
   0

The constructs that were *not* selected by a filter may be returned by
the `~Constructs.inverse_filter` method applied to the results of
filters:

.. code-block:: python
   :caption: *Get the constructs that were not selected by a filter.*

   >>> c = t.constructs.filter_by_type('auxiliary_coordinate')
   >>> c
   <Constructs: auxiliary_coordinate(3)>
   >>> c.inverse_filter()
   <Constructs: cell_measure(1), cell_method(2), coordinate_reference(2), dimension_coordinate(4), domain_ancillary(3), domain_axis(4), field_ancillary(1)>
  
Note that selection by construct type is equivalent to using the
particular method of the field construct for retrieving that type of
metadata construct:

.. code-block:: python
   :caption: *The bespoke methods for retrieving constructs by type
             are equivalent to a selection on all of the metadata
             constructs.*
		
   >>> print(t.constructs.filter_by_type('cell_measure'))
   Constructs:
   {'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>}
   >>> print(t.cell_measures())
   Constructs:
   {'cellmeasure0': <CellMeasure: measure:area(9, 10) km2>}
   
----
   
.. _Metadata-construct-access:

**Metadata construct access**
-----------------------------

An individual metadata construct may be returned, without its
construct key, by any of the following techniques:

* with the `~Field.construct` method of a field construct,

.. code-block:: python
   :caption: *Get the "latitude" metadata construct with its construct
             identity.*
	     
   >>> t.construct('latitude')
   <AuxiliaryCoordinate: latitude(10, 9) degrees_N>

* with the `~Field.construct_key` method of a field construct:

.. code-block:: python
   :caption: *Get the "latitude" metadata construct key with its construct
             identity and use the key to get the construct itself*
	     
   >>> key = t.construct_key('latitude')
   >>> t.construct(key)
   <AuxiliaryCoordinate: latitude(10, 9) degrees_N>

* with the `~Field.construct_item` method of a field construct:

.. code-block:: python
   :caption: *Get the "latitude" metadata construct and its identifier
             via its construct identity.*
      
   >>> key, lat = t.construct_item('latitude')
   ('auxiliarycoordinate0', <AuxiliaryCoordinate: latitude(10, 9) degrees_N>)

* by indexing a `Constructs` instance with  a construct key.

.. code-block:: python
   :caption: *Get the "latitude" metadata construct via its construct
             key and indexing*
	     
   >>> key = t.construct_key('latitude')
   >>> t.constructs[key]
   <AuxiliaryCoordinate: latitude(10, 9) degrees_N>

* with the `~Constructs.get` method of a `Constructs` instance, or

.. code-block:: python
   :caption: *Get the "latitude" metadata construct via its construct
             key and the 'get' method.*
	     
   >>> key = t.construct_key('latitude')
   >>> c = t.constructs.get(key)
   <AuxiliaryCoordinate: latitude(10, 9) degrees_N>

The `~Field.construct` method of the field construct and the
`~Constructs.value` method of the `Constructs` instance will raise an
exception of there is not a unique metadata construct to return, but
this may be replaced with returning a default value or raising a
customised exception:
   
.. code-block:: python
   :caption: *By default an exception is raised if there is not a
             unique construct that meets the criteria. Alternatively,
             the value of the "default" parameter is returned.*

   >>> t.construct('measure:volume')                # Raises Exception
   Traceback (most recent call last):
      ...
   ValueError: Can't return zero constructs
   >>> t.construct('measure:volume', default=False)
   False
   >>> t.construct('measure:volume', default=Exception("my error"))  # Raises Exception
   Traceback (most recent call last):
      ...
   Exception: my error
   >>> c = t.constructs.filter_by_measure("volume")
   >>> len(c)
   0
   >>> d = t.constructs("units=degrees")
   >>> len(d)
   2
   >>> t.construct("units=degrees")  # Raises Exception
   Traceback (most recent call last):
      ...
   ValueError: Field.construct() can't return 2 constructs
   >>> print(t.construct("units=degrees", default=None))
   None

.. _Metadata-construct-properties:

Metadata construct properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Metadata constructs share the :ref:`same API as the field construct
<Properties>` for accessing their properties:

.. code-block:: python
   :caption: *Retrieve the "longitude" metadata construct, set a new
             property, and then inspect all of the properties.*

   >>> lon = q.construct('longitude')
   >>> lon
   <DimensionCoordinate: longitude(8) degrees_east>
   >>> lon.set_property('long_name', 'Longitude')
   >>> lon.properties()
   {'units': 'degrees_east',
    'long_name': 'Longitude',
    'standard_name': 'longitude'}

.. code-block:: python
   :caption: *Get the metadata construct with units of "km2", find its
             canonical identity, and all of its valid identities, that
             may be used for selection by the "filter_by_identity"
             method*

   >>> area = t.constructs.filter_by_property(units='km2').value()
   >>> area
   <CellMeasure: measure:area(9, 10) km2>
   >>> area.identity()
   'measure:area'
   >>> area.identities()
   ['measure:area', 'units=km2', 'ncvar%cell_measure']

.. _Metadata-construct-data:

Metadata construct data
^^^^^^^^^^^^^^^^^^^^^^^

Metadata constructs share the :ref:`a similar API as the field
construct <Data>` as the field construct for accessing their data:

.. code-block:: python
   :caption: *Retrieve the "longitude" metadata construct, inspect its
             data, change the third element of the array, and get the
             data as a numpy array.*
	     
   >>> lon = q.constructs('longitude').value()
   >>> lon
   <DimensionCoordinate: longitude(8) degrees_east>
   >>> lon.data
   <Data(8): [22.5, ..., 337.5] degrees_east>
   >>> lon.data[2]
   <Data(1): [112.5] degrees_east>
   >>> lon.data[2] = 133.33
   >>> print(lon.data.array)
   [22.5 67.5 133.33 157.5 202.5 247.5 292.5 337.5]

The :term:`domain axis constructs` spanned by a particular metadata
construct's data are found with the `~Constructs.get_data_axes` method
of the field construct:

.. code-block:: python
   :caption: *Find the construct keys of the domain axis constructs
             spanned by the data of each metadata construct.*

   >>> key = t.construct_key('latitude')
   >>> key
   'auxiliarycoordinate0'
   >>> t.get_data_axes(key=key)
   ('domainaxis1', 'domainaxis2')
    
The domain axis constructs spanned by all the data of all metadata
construct may be found with the `~Constructs.data_axes` method of the
field construct's `Constructs` instance:

.. code-block:: python
   :caption: *Find the construct keys of the domain axis constructs
             spanned by the data of each metadata construct.*

   >>> t.constructs.data_axes()
   {'auxiliarycoordinate0': ('domainaxis1', 'domainaxis2'),
    'auxiliarycoordinate1': ('domainaxis2', 'domainaxis1'),
    'auxiliarycoordinate2': ('domainaxis1',),
    'cellmeasure0': ('domainaxis2', 'domainaxis1'),
    'dimensioncoordinate0': ('domainaxis0',),
    'dimensioncoordinate1': ('domainaxis1',),
    'dimensioncoordinate2': ('domainaxis2',),
    'dimensioncoordinate3': ('domainaxis3',),
    'domainancillary0': ('domainaxis0',),
    'domainancillary1': ('domainaxis0',),
    'domainancillary2': ('domainaxis1', 'domainaxis2'),
    'fieldancillary0': ('domainaxis1', 'domainaxis2')}

A size one domain axis construct that is *not* spanned by the field
construct's data may still be spanned by the data of metadata
constructs. For example, the data of the field construct ``t``
:ref:`does not span the size one domain axis construct <Data-axes>`
with key ``'domainaxis3'``, but this domain axis construct is spanned
by a "time" dimension coordinate construct (with key
``'dimensioncoordinate3'``). Such a dimension coordinate (i.e. one
that applies to a domain axis construct that is not spanned by the
field construct's data) corresponds to a CF-netCDF scalar coordinate
variable.

----

.. _Time:

**Time**
--------

Constructs representing elapsed time (identified by the presence of
"reference time" units) have data array values that represent elapsed
time since a reference date. These values may be converted into the
date-time objects of the `cftime package
<https://unidata.github.io/cftime/>`_ with the `~Data.datetime_array`
method of the `Data` instance.

.. code-block:: python
   :caption: *Inspect the the values of a "time" construct as elapsed
             times and as date-times.*

   >>> time = q.construct('time')
   >>> time
   <DimensionCoordinate: time(1) days since 2018-12-01 >
   >>> time.get_property('units')
   'days since 2018-12-01'
   >>> time.get_property('calendar', default='standard')
   'standard'
   >>> print(time.data.array)
   [ 31.]
   >>> print(time.data.datetime_array)
   [cftime.DatetimeGregorian(2019, 1, 1, 0, 0, 0, 0, 1, 1)]

----

.. _Domain:

**Domain**
----------

The :ref:`domain of the CF data model <CF-data-model>` is defined
collectively by various other metadata constructs. It is represented
by the `Domain` class. A domain construct may exist independently, or
is accessed from a field construct with its `~Field.domain` attribute,
or `~Field.get_domain` method.

.. code-block:: python
   :caption: *Get the domain, and inspect it.*

   >>> domain = t.domain
   >>> domain
   <Domain: {1, 1, 9, 10}>
   >>> print(domain)
   Dimension coords: atmosphere_hybrid_height_coordinate(1) = [1.5]
                   : grid_latitude(10) = [2.2, ..., -1.76] degrees
                   : grid_longitude(9) = [-4.7, ..., -1.18] degrees
                   : time(1) = [2019-01-01 00:00:00]
   Auxiliary coords: latitude(grid_latitude(10), grid_longitude(9)) = [[53.941, ..., 50.225]] degrees_N
                   : longitude(grid_longitude(9), grid_latitude(10)) = [[2.004, ..., 8.156]] degrees_E
                   : long_name=Grid latitude name(grid_latitude(10)) = [--, ..., kappa]
   Cell measures   : measure:area(grid_longitude(9), grid_latitude(10)) = [[2391.9657, ..., 2392.6009]] km2
   Coord references: atmosphere_hybrid_height_coordinate
                   : rotated_latitude_longitude
   Domain ancils   : ncvar%a(atmosphere_hybrid_height_coordinate(1)) = [10.0] m
                   : ncvar%b(atmosphere_hybrid_height_coordinate(1)) = [20.0]
                   : surface_altitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 270.0]] m
   >>> description = domain.dump(display=False)

The domain construct returned by a field construct is not independent
of its parent field instance, i.e. changes to domain construct are
seen by the field construct, and vice versa. This is because, in this
case, the domain instance is a "view" of the relevant metadata
constructs contained in the field construct.

.. code-block:: python
   :caption: *Change a property of a metadata construct of the domain
             and show that this change appears in the same metadata
             data construct of the parent field, and vice versa.*

   >>> domain_latitude = t.domain.constructs('latitude').value()
   >>> field_latitude = t.constructs('latitude').value()
   >>> domain_latitude.set_property('test', 'set by domain')
   >>> print(field_latitude.get_property('test'))
   set by domain
   >>> field_latitude.set_property('test', 'set by field')
   >>> print(domain_latitude.get_property('test'))
   set by field
   >>> domain_latitude.del_property('test')
   'set by field'
   >>> field_latitude.has_property('test')
   False

----

.. _Domain-axes:

**Domain axes**
---------------

A domain axis metadata construct specifies the number of points along
an independent axis of the field construct's domain and is stored in a
`~cfdm.DomainAxis` instance. The size of the axis is retrieved with
the `~cfdm.DomainAxis.get_size` method of the domain axis construct.

.. code-block:: python
   :caption: *Get the size of a domain axis construct.*

   >>> print(q.domain_axes())
   Constructs:
   {'domainaxis0': <DomainAxis: size(5)>,
    'domainaxis1': <DomainAxis: size(8)>,
    'domainaxis2': <DomainAxis: size(1)>}
   >>> d = q.domain_axes().get('domainaxis1')
   >>> d
   <DomainAxis: size(8)>
   >>> d.get_size()
   8

----

.. _Coordinates:
		
**Coordinates**
---------------

There are two types of coordinate construct, :term:`dimension
<dimension coordinate constructs>` and :term:`auxiliary coordinate
constructs`, which can be retrieved together with the
`~cfdm.Field.coordinates` method of the field construct, as well as
individually with the `~cfdm.Field.auxiliary_coordinates` and
`~cfdm.Field.dimension_coordinates` methods.

.. code-block:: python
   :caption: *Retrieve both types of coordinate constructs.*
      
   >>> print(t.coordinates())
   Constructs:
   {'auxiliarycoordinate0': <AuxiliaryCoordinate: latitude(10, 9) degrees_N>,
    'auxiliarycoordinate1': <AuxiliaryCoordinate: longitude(9, 10) degrees_E>,
    'auxiliarycoordinate2': <AuxiliaryCoordinate: long_name=Grid latitude name(10) >,
    'dimensioncoordinate0': <DimensionCoordinate: atmosphere_hybrid_height_coordinate(1) >,
    'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>,
    'dimensioncoordinate2': <DimensionCoordinate: grid_longitude(9) degrees>,
    'dimensioncoordinate3': <DimensionCoordinate: time(1) days since 2018-12-01 >}

.. _Bounds:

Bounds
^^^^^^

A coordinate construct may contain an array of cell bounds that
provides the extent of each cell by defining the locations of the cell
vertices. This is in addition to the main coordinate data array that
contains a representative grid point location for each cell. The cell
bounds are stored in a `Bounds` class instance that is accessed with
the `~Coordinate.bounds` attribute, or `~Coordinate.get_bounds`
method, of the coordinate construct.

A `Bounds` instance shares the :ref:`the same API as the field
construct <Data>` for accessing its data.

.. code-block:: python
   :caption: *Get the Bounds instance of a coordinate construct and
             inspect its data.*
      
   >>> lon = t.constructs('grid_longitude').value()
   >>> bounds = lon.bounds
   >>> bounds
   <Bounds: grid_longitude(9, 2) >
   >>> bounds.data
   <Data(9, 2): [[-4.92, ..., -0.96]]>
   >>> print(bounds.data.array)
   [[-4.92 -4.48]
    [-4.48 -4.04]
    [-4.04 -3.6 ]
    [-3.6  -3.16]
    [-3.16 -2.72]
    [-2.72 -2.28]
    [-2.28 -1.84]
    [-1.84 -1.4 ]
    [-1.4  -0.96]]

The `Bounds` instance inherits the descriptive properties from its
parent coordinate construct, but it may also have its own properties
(although setting these is not recommended).

.. code-block:: python
   :caption: *Inspect the inherited and bespoke properties of a Bounds
             instance.*
      
   >>> bounds.inherited_properties()
   {'standard_name': 'grid_longitude',
    'units': 'degrees'}  
   >>> bounds.properties()
   {}

.. _Geometry-cells:   

Geometry cells
^^^^^^^^^^^^^^

For many geospatial applications, cell bounds can not be represented
by a simple line or polygon, and different cells may have different
numbers of nodes describing their bounds. For example, if each cell
describes the areal extent of a watershed, then it is likely that some
watersheds will require more nodes than others. Such cells are called
`geometries`_.

If a coordinate construct represents geometries then it will have a
"geometry" attribute (not a :ref:`CF property
<Metadata-construct-properties>`) with one of the values ``'point'``,
``'line'`` or ``'polygon'``.

This is illustrated with the file ``geometry.nc`` (found in the
:ref:`sample datasets <Sample-datasets>`):

.. code-block:: python
   :caption: *Read and inspect a dataset containing geometry cell
             bounds.*

   >>> f = cfdm.read('geometry.nc')[0]
   >>> print(f)
   Field: precipitation_amount (ncvar%pr)
   --------------------------------------
   Data            : precipitation_amount(cf_role=timeseries_id(2), time(4))
   Dimension coords: time(4) = [2000-01-02 00:00:00, ..., 2000-01-05 00:00:00]
   Auxiliary coords: latitude(cf_role=timeseries_id(2)) = [25.0, 7.0] degrees_north
                   : longitude(cf_role=timeseries_id(2)) = [10.0, 40.0] degrees_east
                   : altitude(cf_role=timeseries_id(2)) = [5000.0, 20.0] m
                   : cf_role=timeseries_id(cf_role=timeseries_id(2)) = [b'x1', b'y2']
   Coord references: grid_mapping_name:latitude_longitude
   >>> lon = f.construct('longitude')
   >>> lon.dump()                     
   Auxiliary coordinate: longitude
      standard_name = 'longitude'
      units = 'degrees_east'
      Data(2) = [10.0, 40.0] degrees_east
      Geometry: polygon
      Bounds:axis = 'X'
      Bounds:standard_name = 'longitude'
      Bounds:units = 'degrees_east'
      Bounds:Data(2, 3, 4) = [[[20.0, ..., --]]] degrees_east
      Interior Ring:Data(2, 3) = [[0, ..., --]]
   >>> lon.get_geometry()
   'polygon'

Bounds for geometry cells are also stored in a `Bounds` instance, but
one that always has *two* extra trailing dimensions (rather than
one). The fist trailing dimension indexes the distinct parts of a
geometry, and the second indexes the nodes of each part. When a part
has fewer nodes than another, its nodes dimension is padded with
missing data.


.. code-block:: python
   :caption: *Inspect the geometry nodes.*
 
   >>> print(lon.bounds.data.array)
   [[20.0 10.0  0.0   --]
    [ 5.0 10.0 15.0 10.0]
    [20.0 10.0  0.0   --]]

   [[50.0 40.0 30.0   --]
    [  --   --   --   --]
    [  --   --   --   --]]]

If a cell is composed of multiple polygon parts, an individual polygon
may define an "interior ring", i.e. a region that is to be omitted
from, as opposed to included in, the cell extent. Such cells also have
an interior ring array that spans the same domain axes as the
coordinate cells, with the addition of one extra dimension that
indexes the parts for each cell. This array records whether each
polygon is to be included or excluded from the cell, with values of
``1`` or ``0`` respectively.

.. code-block:: python
   :caption: *Inspect the interior ring information.*
 
   >>> print(lon.get_interior_ring().data.array)
   [[0  1  0]
    [0 -- --]]

Note that it is preferable to access the data type, number of
dimensions, dimension sizes and number of elements of the coordinate
construct via the construct's attributes, rather than the attributes
of the `Data` instance that provides representative values for each
cell. This is because the representative cell values for geometries
are optional, and if they are missing then the construct attributes
are able to infer these attributes from the bounds.
  
When a field construct containing geometries is written to disk, a
CF-netCDF geometry container variable is automatically created, and
the cells encoded with the :ref:`compression <Compression>` techniques
defined in the CF conventions.

----

.. _Domain-ancillaries:
		
**Domain ancillaries**
----------------------

A :term:`domain ancillary <domain ancillary constructs>` construct
provides information which is needed for computing the location of
cells in an alternative :ref:`coordinate system
<Coordinate-systems>`. If a domain ancillary construct provides extra
coordinates then it may contain cell bounds in addition to its main
data array.

.. code-block:: python
   :caption: *Get the data and bounds data of a domain ancillary
             construct.*
      
   >>> a = t.constructs.get('domainancillary0')
   >>> print(a.data.array)
   [10.]
   >>> bounds = a.bounds
   >>> bounds
   <Bounds: ncvar%a_bounds(1, 2) >
   >>> print(bounds.data.array)
   [[  5.  15.]]

----

.. _Coordinate-systems:

**Coordinate systems**
----------------------

A field construct may contain various coordinate systems. Each
coordinate system is either defined by a :term:`coordinate reference
construct <coordinate reference constructs>` that relates dimension
coordinate, auxiliary coordinate and domain ancillary constructs (as
is the case for the field construct ``t``), or is inferred from
dimension and auxiliary coordinate constructs alone (as is the case
for the field construct ``q``).

A coordinate reference construct contains

* references (by construct keys) to the dimension and auxiliary
  coordinate constructs to which it applies, accessed with the
  `~CoordinateReference.coordinates` method of the coordinate
  reference construct;

..

* the zeroes of the dimension and auxiliary coordinate constructs
  which define the coordinate system, stored in a `Datum` instance,
  which is accessed with the `~CoordinateReference.datum` attribute,
  or `~CoordinateReference.get_datum` method, of the coordinate
  reference construct; and

..

* a formula for converting coordinate values taken from the dimension
  or auxiliary coordinate constructs to a different coordinate system,
  stored in a `CoordinateConversion` class instance, which is accessed
  with the `~CoordinateReference.coordinate_conversion` attribute, or
  `~CoordinateReference.get_coordinate_conversion` method, of the
  coordinate reference construct.

.. code-block:: python
   :caption: *Select the vertical coordinate system construct and
             inspect its coordinate constructs.*
     
   >>> crs = t.constructs('standard_name:atmosphere_hybrid_height_coordinate').value()
   >>> crs
   <CoordinateReference: atmosphere_hybrid_height_coordinate>
   >>> crs.dump()
   Coordinate Reference: atmosphere_hybrid_height_coordinate
       Coordinate conversion:computed_standard_name = altitude
       Coordinate conversion:standard_name = atmosphere_hybrid_height_coordinate
       Coordinate conversion:a = domainancillary0
       Coordinate conversion:b = domainancillary1
       Coordinate conversion:orog = domainancillary2
       Datum:earth_radius = 6371007
       Coordinate: dimensioncoordinate0
   >>> crs.coordinates()
   {'dimensioncoordinate0'}

.. code-block:: python
   :caption: *Get the datum and inspect its parameters.*
	     
   >>> crs.datum
   <Datum: Parameters: earth_radius>
   >>> crs.datum.parameters()
   {'earth_radius': 6371007}


.. code-block:: python
   :caption: *Get the coordinate conversion and inspect its parameters
             and referenced domain ancillary constructs.*
	     
   >>> crs.coordinate_conversion
   <CoordinateConversion: Parameters: computed_standard_name, standard_name; Ancillaries: a, b, orog>
   >>> crs.coordinate_conversion.parameters()
   {'computed_standard_name': 'altitude',
    'standard_name': 'atmosphere_hybrid_height_coordinate'}
   >>> crs.coordinate_conversion.domain_ancillaries()
   {'a': 'domainancillary0',
    'b': 'domainancillary1',
    'orog': 'domainancillary2'}    

----

.. _Cell-methods:
   
**Cell methods**
----------------

A cell method construct describes how the data represent the variation
of the physical quantity within the cells of the domain and is stored
in a `~cfdm.CellMethod` instance. A field constructs allows multiple
cell method constructs to be recorded.

.. code-block:: python
   :caption: *Inspect the cell methods. The description follows the CF
             conventions for cell_method attribute strings, apart from
             the use of construct keys instead of netCDF variable
             names for cell method axes identification.*
	     
   >>> print(t.cell_methods())
   Constructs:
   {'cellmethod0': <CellMethod: domainaxis1: domainaxis2: mean where land (interval: 0.1 degrees)>,
    'cellmethod1': <CellMethod: domainaxis3: maximum>}

The application of cell methods is not commutative (e.g. a mean of
variances is generally not the same as a variance of means), and the
cell methods are assumed to have been applied in the order in which
they were added to the field construct during :ref:`field construct
creation <Field-creation-in-memory>`.

The axes to which the method applies, the method itself, and any
qualifying properties are accessed with the
`~cfdm.CellMethod.get_axes`, `~cfdm.CellMethod.get_method`, ,
`~cfdm.CellMethod.get_qualifier` and `~cfdm.CellMethod.qualifiers`
methods of the cell method construct.

.. code-block:: python
   :caption: *Get the domain axes constructs to which the cell method
             construct applies, and the method and other properties.*
     
   >>> cm = t.constructs('method:mean').value()
   >>> cm
   <CellMethod: domainaxis1: domainaxis2: mean where land (interval: 0.1 degrees)>)
   >>> cm.get_axes()
   ('domainaxis1', 'domainaxis2')
   >>> cm.get_method()
   'mean'
   >>> cm.qualifiers()
   {'interval': [<Data(): 0.1 degrees>], 'where': 'land'}
   >>> cm.get_qualifier('where')
   'land'

----

.. _Field-ancillaries:
		
**Field ancillaries**
---------------------

A :term:`field ancillary construct <field ancillary constructs>`
provides metadata which are distributed over the same domain as the
field construct itself. For example, if a field construct holds a data
retrieved from a satellite instrument, a field ancillary construct
might provide the uncertainty estimates for those retrievals (varying
over the same spatiotemporal domain).

.. code-block:: python
   :caption: *Get the properties and data of a field ancillary
             construct.*

   >>> a = t.get_construct('fieldancillary0')
   >>> a
   <FieldAncillary: air_temperature standard_error(10, 9) K>
   >>> a.properties()
   {'standard_name': 'air_temperature standard_error',
    'units': 'K'}
   >>> a.data
   <Data(10, 9): [[0.76, ..., 0.32]] K>

----

.. _Field-creation-in-memory:

**Field creation in memory**
----------------------------

There are various methods for creating a field construct in memory:

* :ref:`Ab initio creation <Ab-initio-creation>`: Instantiate
  instances of field and metadata construct classes and manually
  provide the connections between them.
..

* :ref:`Command modification <Command-modification>`: Produce the
  commands that would create an already existing field construct, and
  then modify them.

..

* :ref:`Creation by conversion <Creation-by-conversion>`: Convert a
  single metadata construct already in memory to an independent field
  construct

..

* :ref:`Creation by reading <Creation-by-reading>`: Create field
  constructs from the netCDF variables in a dataset.

.. _Ab-initio-creation:

Ab initio creation
^^^^^^^^^^^^^^^^^^

Ab initio creation of a field construct has three stages:

**Stage 1:** The field construct is created without metadata
constructs.

..
   
**Stage 2:** Metadata constructs are created independently.

..

**Stage 3:** The metadata constructs are inserted into the field
construct with cross-references to other, related metadata constructs
if required. For example, an auxiliary coordinate construct is related
to an ordered list of the :term:`domain axis constructs` which correspond to
its data array dimensions.

There are two equivalent approaches to **stages 1** and **2**.

Either as much of the content as possible is specified during object
instantiation:

.. code-block:: python
   :caption: *Create a field construct with a "standard_name"
             property. Create dimension coordinate and field ancillary
             constructs, both with properties and data.*
	     
   >>> p = cfdm.Field(properties={'standard_name': 'precipitation_flux'})
   >>> p
   <Field: precipitation_flux>
   >>> dc = cfdm.DimensionCoordinate(properties={'long_name': 'Longitude'},
   ...                               data=cfdm.Data([0, 1, 2.]))
   >>> dc
   <DimensionCoordinate: long_name=Longitude(3) >
   >>> fa = cfdm.FieldAncillary(
   ...        properties={'standard_name': 'precipitation_flux status_flag'},
   ...        data=cfdm.Data(numpy.array([0, 0, 2], dtype='int8')))
   >>> fa
   <FieldAncillary: precipitation_flux status_flag(3) >

or else some or all content is added after instantiation via object
methods:

.. code-block:: python
   :caption: *Create empty constructs and provide them with properties
             and data after instantiation.*
	     
   >>> p = cfdm.Field()
   >>> p
   <Field: >
   >>> p.set_property('standard_name', 'precipitation_flux')
   >>> p
   <Field: precipitation_flux>
   >>> dc = cfdm.DimensionCoordinate()
   >>> dc
   <DimensionCoordinate:  >
   >>> dc.set_property('long_name', 'Longitude')
   >>> dc.set_data(cfdm.Data([1, 2, 3.]))
   >>> dc
   <DimensionCoordinate: long_name=Longitude(3) >
   >>> fa = cfdm.FieldAncillary(
   ...        data=cfdm.Data(numpy.array([0, 0, 2], dtype='int8')))
   >>> fa
   <FieldAncillary: (3) >
   >>> fa.set_property('standard_name', 'precipitation_flux status_flag')
   >>> fa
   <FieldAncillary: precipitation_flux status_flag(3) >

For **stage 3**, the `~cfdm.Field.set_construct` method of the field
construct is used for setting metadata constructs and mapping data
array dimensions to domain axis constructs. This method returns the
construct key for the metadata construct which can be used when other
metadata constructs are added to the field (e.g. to specify which
domain axis constructs correspond to a data array), or when other
metadata constructs are created (e.g. to identify the domain ancillary
constructs forming part of a coordinate reference construct):

.. code-block:: python
   :caption: *Set a domain axis construct and use its construct key
             when setting the dimension coordinate construct. Also
             create a cell method construct that applies to the domain
             axis construct.*
	     
   >>> longitude_axis = p.set_construct(cfdm.DomainAxis(3))
   >>> longitude_axis
   'domainaxis0'
   >>> key = p.set_construct(dc, axes=longitude_axis)
   >>> key
   'dimensioncoordinate0'
   >>> cm = cfdm.CellMethod(axes=longitude_axis, method='minimum')
   >>> p.set_construct(cm)
   'cellmethod0'
   
In general, the order in which metadata constructs are added to the
field does not matter, except when one metadata construct is required
by another, in which case the former must be added to the field first
so that its construct key is available to the latter. Cell method
constructs must, however, be set in the relative order in which their
methods were applied to the data.

The domain axis constructs spanned by a metadata construct's data may
be changed after insertion with the `~Field.set_data_axes` method of
the field construct.

.. Code Block Start 1
   
.. code-block:: python
   :caption: *Create a field construct with properties; data; and
             domain axis, cell method and dimension coordinate
             metadata constructs (data arrays have been generated with
             dummy values using numpy.arange).*

   import numpy
   import cfdm

   # Initialise the field construct with properties
   Q = cfdm.Field(
	     properties={'project': 'research',
                              'standard_name': 'specific_humidity',
                              'units': '1'})
			      
   # Create the domain axis constructs
   domain_axisT = cfdm.DomainAxis(1)
   domain_axisY = cfdm.DomainAxis(5)
   domain_axisX = cfdm.DomainAxis(8)

   # Insert the domain axis constructs into the field. The
   # set_construct method returns the domain axis construct key that
   # will be used later to specify which domain axis corresponds to
   # which dimension coordinate construct.
   axisT = Q.set_construct(domain_axisT)
   axisY = Q.set_construct(domain_axisY)
   axisX = Q.set_construct(domain_axisX)

   # Create and insert the field construct data
   data = cfdm.Data(numpy.arange(40.).reshape(5, 8))
   Q.set_data(data, axes=[axisY, axisX])

   # Create the cell method constructs
   cell_method1 = cfdm.CellMethod(axes='area', method='mean')

   cell_method2 = cfdm.CellMethod()
   cell_method2.set_axes(axisT)
   cell_method2.set_method('maximum')

   # Insert the cell method constructs into the field in the same
   # order that their methods were applied to the data
   Q.set_construct(cell_method1)
   Q.set_construct(cell_method2)

   # Create a "time" dimension coordinate construct, with coordinate
   # bounds
   dimT = cfdm.DimensionCoordinate(
                               properties={'standard_name': 'time',
                                           'units': 'days since 2018-12-01'},
                               data=cfdm.Data([15.5]),
                               bounds=cfdm.Bounds(data=cfdm.Data([[0,31.]])))

   # Create a "longitude" dimension coordinate construct, without
   # coordinate bounds
   dimX = cfdm.DimensionCoordinate(data=cfdm.Data(numpy.arange(8.)))
   dimX.set_properties({'standard_name': 'longitude',
                        'units': 'degrees_east'})

   # Create a "longitude" dimension coordinate construct
   dimY = cfdm.DimensionCoordinate(properties={'standard_name': 'latitude',
		                               'units': 'degrees_north'})
   array = numpy.arange(5.)
   dimY.set_data(cfdm.Data(array))

   # Create and insert the latitude coordinate bounds
   bounds_array = numpy.empty((5, 2))
   bounds_array[:, 0] = array - 0.5
   bounds_array[:, 1] = array + 0.5
   bounds = cfdm.Bounds(data=cfdm.Data(bounds_array))
   dimY.set_bounds(bounds)

   # Insert the dimension coordinate constructs into the field,
   # specifying to # which domain axis each one corresponds
   Q.set_construct(dimT, axes=axisT)
   Q.set_construct(dimY, axes=axisY)
   Q.set_construct(dimX, axes=axisX)

.. Code Block End 1

.. code-block:: python
   :caption: *Inspect the new field construct.* 
	  
   >>> Q.dump()
   ------------------------
   Field: specific_humidity
   ------------------------
   project = 'research'
   standard_name = 'specific_humidity'
   units = '1'
   
   Data(latitude(5), longitude(8)) = [[0.0, ..., 39.0]] 1
   
   Cell Method: area: mean
   Cell Method: time(1): maximum
   
   Domain Axis: latitude(5)
   Domain Axis: longitude(8)
   Domain Axis: time(1)
   
   Dimension coordinate: time
       standard_name = 'time'
       units = 'days since 2018-12-01'
       Data(time(1)) = [2018-12-16 12:00:00]
       Bounds:Data(time(1), 2) = [[2018-12-01 00:00:00, 2019-01-01 00:00:00]]
   
   Dimension coordinate: latitude
       standard_name = 'latitude'
       units = 'degrees_north'
       Data(latitude(5)) = [0.0, ..., 4.0] degrees_north
       Bounds:Data(latitude(5), 2) = [[-0.5, ..., 4.5]] degrees_north
   
   Dimension coordinate: longitude
       standard_name = 'longitude'
       units = 'degrees_east'
       Data(longitude(8)) = [0.0, ..., 7.0] degrees_east

The ``Conventions`` property does not need to be set because it is
automatically included in output files as a netCDF global
``Conventions`` attribute, either as the CF version of the cfdm
package (as returned by the `cfdm.CF` function), or else specified via
the *Conventions* keyword of the `cfdm.write` function. See the
section on :ref:`Writing-to-disk` for details on how to specify
additional conventions.

If this field were to be written to a netCDF dataset then, in the
absence of predefined names, default netCDF variable and dimension
names would be automatically generated (based on standard names where
they exist). The setting of bespoke netCDF names is, however, easily
done with the :ref:`netCDF interface <NetCDF-interface>`.

.. code-block:: python
   :caption: *Set netCDF variable and dimension names for the field
             and metadata constructs.*

   Q.nc_set_variable('q')

   domain_axisT.nc_set_dimension('time')
   domain_axisY.nc_set_dimension('lat')
   domain_axisX.nc_set_dimension('lon')

   dimT.nc_set_variable('time')
   dimY.nc_set_variable('lat')
   dimX.nc_set_variable('lon')

Here is a more complete example which creates a field construct that
contains every type of metadata construct (again, data arrays have
been generated with dummy values using `numpy.arange`):

.. Code Block Start 2
   
.. code-block:: python
   :caption: *Create a field construct that contains at least one
             instance of each type of metadata construct.*

   import numpy
   import cfdm
   
   # Initialize the field construct
   tas = cfdm.Field(
       properties={'project': 'research',
                   'standard_name': 'air_temperature',
                   'units': 'K'})
   
   # Create and set domain axis constructs
   axis_T = tas.set_construct(cfdm.DomainAxis(1))
   axis_Z = tas.set_construct(cfdm.DomainAxis(1))
   axis_Y = tas.set_construct(cfdm.DomainAxis(10))
   axis_X = tas.set_construct(cfdm.DomainAxis(9))
   
   # Set the field construct data
   tas.set_data(cfdm.Data(numpy.arange(90.).reshape(10, 9)),
                axes=[axis_Y, axis_X])
   
   # Create and set the cell method constructs
   cell_method1 = cfdm.CellMethod(
             axes=[axis_Y, axis_X],
	     method='mean',
             qualifiers={'where': 'land',
                         'interval': [cfdm.Data(0.1, units='degrees')]})
   
   cell_method2 = cfdm.CellMethod(axes=axis_T, method='maximum')
   
   tas.set_construct(cell_method1)
   tas.set_construct(cell_method2)
   
   # Create and set the field ancillary constructs
   field_ancillary = cfdm.FieldAncillary(
                properties={'standard_name': 'air_temperature standard_error',
                             'units': 'K'},
                data=cfdm.Data(numpy.arange(90.).reshape(10, 9)))
   
   tas.set_construct(field_ancillary, axes=[axis_Y, axis_X])
   
   # Create and set the dimension coordinate constructs
   dimension_coordinate_T = cfdm.DimensionCoordinate(
                              properties={'standard_name': 'time',
                                          'units': 'days since 2018-12-01'},
                              data=cfdm.Data([15.5]),
                              bounds=cfdm.Bounds(data=cfdm.Data([[0., 31]])))
   
   dimension_coordinate_Z = cfdm.DimensionCoordinate(
           properties={'computed_standard_name': 'altitude',
                       'standard_name': 'atmosphere_hybrid_height_coordinate'},
           data = cfdm.Data([1.5]),
           bounds=cfdm.Bounds(data=cfdm.Data([[1.0, 2.0]])))
   
   dimension_coordinate_Y = cfdm.DimensionCoordinate(
           properties={'standard_name': 'grid_latitude',
                       'units': 'degrees'},
           data=cfdm.Data(numpy.arange(10.)),
           bounds=cfdm.Bounds(data=cfdm.Data(numpy.arange(20).reshape(10, 2))))
   
   dimension_coordinate_X = cfdm.DimensionCoordinate(
           properties={'standard_name': 'grid_longitude',
                       'units': 'degrees'},
       data=cfdm.Data(numpy.arange(9.)),
       bounds=cfdm.Bounds(data=cfdm.Data(numpy.arange(18).reshape(9, 2))))
   
   dim_T = tas.set_construct(dimension_coordinate_T, axes=axis_T)
   dim_Z = tas.set_construct(dimension_coordinate_Z, axes=axis_Z)
   dim_Y = tas.set_construct(dimension_coordinate_Y, axes=axis_Y)
   dim_X = tas.set_construct(dimension_coordinate_X, axes=axis_X)
   
   # Create and set the auxiliary coordinate constructs
   auxiliary_coordinate_lat = cfdm.AuxiliaryCoordinate(
                         properties={'standard_name': 'latitude',
                                     'units': 'degrees_north'},
                         data=cfdm.Data(numpy.arange(90.).reshape(10, 9)))
   
   auxiliary_coordinate_lon = cfdm.AuxiliaryCoordinate(
                     properties={'standard_name': 'longitude',
                                 'units': 'degrees_east'},
                     data=cfdm.Data(numpy.arange(90.).reshape(9, 10)))
   
   array = numpy.ma.array(list('abcdefghij'))
   array[0] = numpy.ma.masked
   auxiliary_coordinate_name = cfdm.AuxiliaryCoordinate(
                          properties={'long_name': 'Grid latitude name'},
                          data=cfdm.Data(array))
   
   aux_LAT  = tas.set_construct(auxiliary_coordinate_lat, axes=[axis_Y, axis_X])
   aux_LON  = tas.set_construct(auxiliary_coordinate_lon, axes=[axis_X, axis_Y])
   aux_NAME = tas.set_construct(auxiliary_coordinate_name, axes=[axis_Y])
   
   # Create and set domain ancillary constructs
   domain_ancillary_a = cfdm.DomainAncillary(
                      properties={'units': 'm'},
                      data=cfdm.Data([10.]),
                      bounds=cfdm.Bounds(data=cfdm.Data([[5., 15.]])))
   
   domain_ancillary_b = cfdm.DomainAncillary(
                          properties={'units': '1'},
                          data=cfdm.Data([20.]),
                          bounds=cfdm.Bounds(data=cfdm.Data([[14, 26.]])))
   
   domain_ancillary_orog = cfdm.DomainAncillary(
                             properties={'standard_name': 'surface_altitude',
                                         'units': 'm'},
                             data=cfdm.Data(numpy.arange(90.).reshape(10, 9)))
   
   domain_anc_A    = tas.set_construct(domain_ancillary_a, axes=axis_Z)
   domain_anc_B    = tas.set_construct(domain_ancillary_b, axes=axis_Z)
   domain_anc_OROG = tas.set_construct(domain_ancillary_orog,
                                       axes=[axis_Y, axis_X])

   # Create the datum for the coordinate reference constructs
   datum = cfdm.Datum(parameters={'earth_radius': 6371007.})

   # Create the coordinate conversion for the horizontal coordinate
   # reference construct
   coordinate_conversion_h = cfdm.CoordinateConversion(
                 parameters={'grid_mapping_name': 'rotated_latitude_longitude',
                             'grid_north_pole_latitude': 38.0,
                             'grid_north_pole_longitude': 190.0})
   
   # Create the coordinate conversion for the vertical coordinate
   # reference construct
   coordinate_conversion_v = cfdm.CoordinateConversion(
            parameters={'standard_name': 'atmosphere_hybrid_height_coordinate',
                        'computed_standard_name': 'altitude'},
            domain_ancillaries={'a': domain_anc_A,
                                'b': domain_anc_B,
                                'orog': domain_anc_OROG})
   
   # Create the vertical coordinate reference construct
   horizontal_crs = cfdm.CoordinateReference(
                      datum=datum,
                      coordinate_conversion=coordinate_conversion_h,
                      coordinates=[dim_X,
                                   dim_Y,
                                   aux_LAT,
                                   aux_LON])

   # Create the vertical coordinate reference construct
   vertical_crs = cfdm.CoordinateReference(
                    datum=datum,
                    coordinate_conversion=coordinate_conversion_v,
                    coordinates=[dim_Z])

   # Set the coordinate reference constructs
   tas.set_construct(horizontal_crs)
   tas.set_construct(vertical_crs)
   
   # Create and set the cell measure constructs
   cell_measure = cfdm.CellMeasure(measure='area',
                    properties={'units': 'km2'},
                    data=cfdm.Data(numpy.arange(90.).reshape(9, 10)))
   
   tas.set_construct(cell_measure, axes=[axis_X, axis_Y])

.. Code Block End 2

The new field construct may now be inspected:

.. code-block:: python
   :caption: *Inspect the new field construct.*

   >>> print(tas)
   Field: air_temperature
   ----------------------
   Data            : air_temperature(grid_latitude(10), grid_longitude(9)) K
   Cell methods    : grid_latitude(10): grid_longitude(9): mean where land (interval: 0.1 degrees) time(1): maximum
   Field ancils    : air_temperature standard_error(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 89.0]] K
   Dimension coords: time(1) = [2018-12-16 12:00:00]
                   : atmosphere_hybrid_height_coordinate(1) = [1.5]
                   : grid_latitude(10) = [0.0, ..., 9.0] degrees
                   : grid_longitude(9) = [0.0, ..., 8.0] degrees
   Auxiliary coords: latitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 89.0]] degrees_north
                   : longitude(grid_longitude(9), grid_latitude(10)) = [[0.0, ..., 89.0]] degrees_east
                   : long_name=Grid latitude name(grid_latitude(10)) = [--, ..., j]
   Cell measures   : measure:area(grid_longitude(9), grid_latitude(10)) = [[0.0, ..., 89.0]] km2
   Coord references: atmosphere_hybrid_height_coordinate
                   : rotated_latitude_longitude
   Domain ancils   : key%domainancillary0(atmosphere_hybrid_height_coordinate(1)) = [10.0] m
                   : key%domainancillary1(atmosphere_hybrid_height_coordinate(1)) = [20.0] 1
                   : surface_altitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 89.0]] m
	
.. _Command-modification:

Command modification
^^^^^^^^^^^^^^^^^^^^

It is sometimes convenient to produce the commands that would create
an already existing field construct, and then modify them to create
the desired field construct. The commands are produced by the
`~Field.creation_commands` method of the existing field construct.

.. code-block:: python
   :caption: *Produce the commands that would create an existing field
             construct.*
	
   >>> q, t = cfdm.read('file.nc')
   >>> print(q.creation_commands())
   #
   # field: specific_humidity
   field = cfdm.Field()
   field.set_properties({'Conventions': 'CF-1.9', 'project': 'research', 'standard_name': 'specific_humidity', 'units': '1'})
   field.nc_set_variable('q')
   data = cfdm.Data([[0.007, 0.034, 0.003, 0.014, 0.018, 0.037, 0.024, 0.029], [0.023, 0.036, 0.045, 0.062, 0.046, 0.073, 0.006, 0.066], [0.11, 0.131, 0.124, 0.146, 0.087, 0.103, 0.057, 0.011], [0.029, 0.059, 0.039, 0.07, 0.058, 0.072, 0.009, 0.017], [0.006, 0.036, 0.019, 0.035, 0.018, 0.037, 0.034, 0.013]], units='1', dtype='f8')
   field.set_data(data)
   #
   # domain_axis: ncdim%lat
   c = cfdm.DomainAxis()
   c.set_size(5)
   c.nc_set_dimension('lat')
   field.set_construct(c, key='domainaxis0', copy=False)
   #
   # domain_axis: ncdim%lon
   c = cfdm.DomainAxis()
   c.set_size(8)
   c.nc_set_dimension('lon')
   field.set_construct(c, key='domainaxis1', copy=False)
   #
   # domain_axis:
   c = cfdm.DomainAxis()
   c.set_size(1)
   field.set_construct(c, key='domainaxis2', copy=False)
   #
   # dimension_coordinate: latitude
   c = cfdm.DimensionCoordinate()
   c.set_properties({'units': 'degrees_north', 'standard_name': 'latitude'})
   c.nc_set_variable('lat')
   data = cfdm.Data([-75.0, -45.0, 0.0, 45.0, 75.0], units='degrees_north', dtype='f8')
   c.set_data(data)
   b = cfdm.Bounds()
   b.nc_set_variable('lat_bnds')
   data = cfdm.Data([[-90.0, -60.0], [-60.0, -30.0], [-30.0, 30.0], [30.0, 60.0], [60.0, 90.0]], units='degrees_north', dtype='f8')
   b.set_data(data)
   c.set_bounds(b)
   field.set_construct(c, axes=('domainaxis0',), key='dimensioncoordinate0', copy=False)
   #
   # dimension_coordinate: longitude
   c = cfdm.DimensionCoordinate()
   c.set_properties({'units': 'degrees_east', 'standard_name': 'longitude'})
   c.nc_set_variable('lon')
   data = cfdm.Data([22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5], units='degrees_east', dtype='f8')
   c.set_data(data)
   b = cfdm.Bounds()
   b.nc_set_variable('lon_bnds')
   data = cfdm.Data([[0.0, 45.0], [45.0, 90.0], [90.0, 135.0], [135.0, 180.0], [180.0, 225.0], [225.0, 270.0], [270.0, 315.0], [315.0, 360.0]], units='degrees_east', dtype='f8')
   b.set_data(data)
   c.set_bounds(b)
   field.set_construct(c, axes=('domainaxis1',), key='dimensioncoordinate1', copy=False)
   #
   # dimension_coordinate: time
   c = cfdm.DimensionCoordinate()
   c.set_properties({'units': 'days since 2018-12-01', 'standard_name': 'time'})
   c.nc_set_variable('time')
   data = cfdm.Data([31.0], units='days since 2018-12-01', dtype='f8')
   c.set_data(data)
   field.set_construct(c, axes=('domainaxis2',), key='dimensioncoordinate2', copy=False)
   #
   # cell_method: mean
   c = cfdm.CellMethod()
   c.set_method('mean')
   c.set_axes(('area',))
   field.set_construct(c)
   #
   # field data axes
   field.set_data_axes(('domainaxis0', 'domainaxis1'))

Some example fields are always available from the `cfdm.example_field`
function.
	  
.. _Creating-data-from-an-array-on-disk:

Creating data from an array on disk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All the of above examples use arrays in memory to construct the data
instances for the field and metadata constructs. It is, however,
possible to create data from arrays that reside on disk. The
`cfdm.read` function creates data in this manner. A pointer to an
array in a netCDF file can be stored in a `~cfdm.NetCDFArray`
instance, which is is used to initialize a `~cfdm.Data` instance.

.. code-block:: python
   :caption: *Define a variable from a dataset with the netCDF package
             and use it to create a NetCDFArray instance with which to
             initialize a Data instance.*
		
   >>> import netCDF4
   >>> nc = netCDF4.Dataset('file.nc', 'r')
   >>> v = nc.variables['ta']
   >>> netcdf_array = cfdm.NetCDFArray(filename='file.nc', ncvar='ta',
   ...	                               dtype=v.dtype, ndim=v.ndim,
   ...	     		  	       shape=v.shape, size=v.size)
   >>> data_disk = cfdm.Data(netcdf_array)

  
.. code-block:: python
   :caption: *Read the netCDF variable's data into memory and
             initialise another Data instance with it. Compare the
             values of the two data instances.*

   >>> numpy_array = v[...]
   >>> data_memory = cfdm.Data(numpy_array)
   >>> data_disk.equals(data_memory)
   True

Note that data type, number of dimensions, dimension sizes and number
of elements of the array on disk that are used to initialize the
`~cfdm.NetCDFArray` instance are those expected by the CF data model,
which may be different to those of the netCDF variable in the file
(although they are the same in the above example). For example, a
netCDF character array of shape ``(12, 9)`` is viewed in cfdm as a
one-dimensional string array of shape ``(12,)``.

.. _Creation-by-conversion:

Creation by conversion
^^^^^^^^^^^^^^^^^^^^^^

An independent field construct may be created from an existing
metadata construct using `~Field.convert` method of the field
construct, which identifies a unique metadata construct and returns a
new field construct based on its properties and data. The new field
construct always has :term:`domain axis constructs` corresponding to
the data, and (by default) any other metadata constructs that further
define its domain.

.. code-block:: python
   :caption: *Create an independent field construct from the "surface
             altitude" metadata construct.*

   >>> key = tas.construct_key('surface_altitude')
   >>> orog = tas.convert(key)
   >>> print(orog)
   Field: surface_altitude
   -----------------------
   Data            : surface_altitude(grid_latitude(10), grid_longitude(9)) m
   Dimension coords: grid_latitude(10) = [0.0, ..., 9.0] degrees
                   : grid_longitude(9) = [0.0, ..., 8.0] degrees
   Auxiliary coords: latitude(grid_latitude(10), grid_longitude(9)) = [[0.0, ..., 89.0]] degrees_north
                   : longitude(grid_longitude(9), grid_latitude(10)) = [[0.0, ..., 89.0]] degrees_east
                   : long_name=Grid latitude name(grid_latitude(10)) = [--, ..., j]
   Cell measures   : measure:area(grid_longitude(9), grid_latitude(10)) = [[0.0, ..., 89.0]] km2
   Coord references: rotated_latitude_longitude

The `~Field.convert` method has an option to only include domain axis
constructs in the new field construct, with no other metadata
constructs.

.. code-block:: python
   :caption: *Create an independent field construct from the "surface
             altitude" metadata construct, but without a complete
             domain.*

   >>> orog1 = tas.convert(key, full_domain=False) 
   >>> print(orog1)
   Field: surface_altitude
   -----------------------
   Data            : surface_altitude(key%domainaxis2(10), key%domainaxis3(9)) m

.. _Creation-by-reading:

Creation by reading
^^^^^^^^^^^^^^^^^^^

The `cfdm.read` function :ref:`reads a netCDF dataset
<Reading-datasets>` and returns the contents as a list of zero or more
field constructs, each one corresponding to a unique CF-netCDF data
variable in the dataset. For example, the field construct ``tas`` that
was created manually can be :ref:`written to a netCDF dataset
<Writing-to-disk>` and then read back into memory:

.. code-block:: python
   :caption: *Write the field construct that was created manually to
             disk, and then read it back into a new field construct.*

   >>> cfdm.write(tas, 'tas.nc')
   >>> f = cfdm.read('tas.nc')
   >>> f
   [<Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>]

The `cfdm.read` function also allows field constructs to be derived
directly from the netCDF variables that correspond to particular types
metadata constructs. In this case, the new field constructs will have
a domain limited to that which can be inferred from the corresponding
netCDF variable, but without the connections that are defined by the
parent netCDF data variable. This will often result in a new field
construct that has fewer metadata constructs than one created with the
`~Field.convert` method.

.. code-block:: python
   :caption: *Read the file, treating formula terms netCDF variables
             (which map to domain ancillary constructs) as additional
             CF-netCDF data variables.*

   >>> fields = cfdm.read('tas.nc', extra='domain_ancillary')
   >>> fields
   [<Field: ncvar%a(atmosphere_hybrid_height_coordinate(1)) m>,
    <Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>,
    <Field: ncvar%b(atmosphere_hybrid_height_coordinate(1)) 1>,
    <Field: surface_altitude(grid_latitude(10), grid_longitude(9)) m>]
   >>> orog_from_file = fields[3]
   >>> print(orog_from_file)
   Field: surface_altitude (ncvar%surface_altitude)
   ------------------------------------------------
   Data            : surface_altitude(grid_latitude(10), grid_longitude(9)) m
   Dimension coords: grid_latitude(10) = [0.0, ..., 9.0] degrees
                   : grid_longitude(9) = [0.0, ..., 8.0] degrees

Comparing the field constructs ``orog_from_file`` (created with
`cfdm.read`) and ``orog`` (created with the `~Field.convert` method of
the ``tas`` field construct), the former lacks the auxiliary
coordinate, cell measure and coordinate reference constructs of the
latter. This is because the surface altitude netCDF variable in
``tas.nc`` does not have the ``coordinates``, ``cell_measures`` nor
``grid_mapping`` netCDF attributes that would link it to auxiliary
coordinate, cell measure and grid mapping netCDF variables.


Creating compressed constructs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
----

.. _Copying:

**Copying**
-----------

A field construct may be copied with its `~Field.copy` method. This
produces a "deep copy", i.e. the new field construct is completely
independent of the original field.

.. code-block:: python
   :caption: *Copy a field construct and change elements of the copy,
             showing that the original field construct has not been
             altered.*
     
   >>> u = t.copy()
   >>> u.data[0, 0, 0] = -1e30
   >>> u.data[0, 0, 0]
   <Data(1, 1, 1): [[[-1e+30]]] K>
   >>> t.data[0, 0, 0]
   <Data(1, 1, 1): [[[-1.0]]] K>
   >>> u.del_construct('grid_latitude')
   <DimensionCoordinate: grid_latitude(10) degrees>
   >>> u.constructs('grid_latitude')
   Constructs:
   {}
   >>> t.constructs('grid_latitude')
   Constructs:
   {'dimensioncoordinate1': <DimensionCoordinate: grid_latitude(10) degrees>}

Equivalently, the `copy.deepcopy` function may be used:

.. code-block:: python
   :caption: *Copy a field construct with the built-in copy module.*
	    
   >>> import copy
   >>> u = copy.deepcopy(t)

Metadata constructs may be copied individually in the same manner:

.. code-block:: python
   :caption: *Copy a metadata construct.*

   >>> orog = t.constructs('surface_altitude').value().copy()

Arrays within `Data` instances are copied with a `copy-on-write
<https://en.wikipedia.org/wiki/Copy-on-write>`_ technique. This means
that a copy takes up very little memory, even when the original
constructs contain very large data arrays, and the copy operation is
fast.

----

.. _Equality:

**Equality**
------------

Whether or not two field constructs are the same is tested with either
field construct's `~Field.equals` method.

.. code-block:: python
   :caption: *A field construct is always equal to itself, a copy of
             itself and a complete subspace of itself. The "verbose"
             keyword will give some (but not necessarily all) of the
             reasons why two field constructs are not the same.*
	     
   >>> t.equals(t)
   True
   >>> t.equals(t.copy())
   True
   >>> t.equals(t[...])
   True
   >>> t.equals(q)
   False
   >>> t.equals(q, verbose=True)
   Field: Different units: 'K', '1'
   Field: Different properties
   False

Equality is strict by default. This means that for two field
constructs to be considered equal they must have corresponding
metadata constructs and for each pair of constructs:

* the descriptive properties must be the same (with the exception of
  the field construct's ``Conventions`` property, which is never
  checked), and vector-valued properties must have same the size and
  be element-wise equal, and
  
* if there are data arrays then they must have same shape, data type
  and be element-wise equal.

Two real numbers :math:`x` and :math:`y` are considered equal if
:math:`|x - y| \le a_{tol} + r_{tol}|y|`, where :math:`a_{tol}` (the
tolerance on absolute differences) and :math:`r_{tol}` (the tolerance
on relative differences) are positive, typically very small
numbers. By default both are set to the system epsilon (the difference
between 1 and the least value greater than 1 that is representable as
a float). Their values may be inspected and changed with the
`cfdm.atol` and `cfdm.rtol` functions.

Note that the above equation is not symmetric in :math:`x` and
:math:`y`, so that for two fields ``f1`` and ``f2``, ``f1.equals(f2)``
may be different from ``f2.equals(f1)`` in some rare cases.
   
.. code-block:: python
   :caption: *The atol and rtol functions allow the numerical equality
             tolerances to be inspected and changed.*
      
   >>> print(cfdm.atol())
   2.220446049250313e-16
   >>> print(cfdm.rtol())
   2.220446049250313e-16
   >>> original = cfdm.rtol(0.00001)
   >>> print(cfdm.rtol())
   1e-05
   >>> print(cfdm.rtol(original))
   1e-05
   >>> print(cfdm.rtol())
   2.220446049250313e-16

The :math:`a_{tol}` and :math:`r_{tol}` constants may be set for a
runtime context established using a `with` statement.

.. code-block:: python
   :caption: *Create a runtime contenxt with a different value of
             'atol'.*
	     
   >>> print(cfdm.atol())
   2.220446049250313e-16
   >>> with cfdm.atol(1e-5):
   ...     print(cfdm.atol())
   ...     
   1e-05
   >>> print(cfdm.atol())
   2.220446049250313e-16
   
NetCDF elements, such as netCDF variable and dimension names, do not
constitute part of the CF data model and so are not checked on any
construct.

The `~Field.equals` method has optional parameters for modifying the
criteria for considering two fields to be equal:

* named properties may be omitted from the comparison,

* fill value and missing data value properties may be ignored,

* the data type of data arrays may be ignored, and

* the tolerances on absolute and relative differences for numerical
  comparisons may be temporarily changed, without changing the default
  settings.

Metadata constructs may also be tested for equality:

.. code-block:: python
   :caption: *Metadata constructs also have an equals method, that
             behaves in a similar manner.*
	  
   >>> orog = t.constructs('surface_altitude').value()
   >>> orog.equals(orog.copy())
   True

----

.. _NetCDF-interface:

**NetCDF interface**
--------------------

The logical CF data model is independent of netCDF, but the CF
conventions are designed to enable the processing and sharing of
datasets stored in netCDF files. Therefore, the cfdm package includes
methods for recording and editing netCDF elements that are not part of
the CF model, but are nonetheless often required to interpret and
create CF-netCDF datasets. See the section on :ref:`philosophy
<philosophy>` for a further discussion.

When a netCDF dataset is read, netCDF elements (such as dimension and
variable names, and some attribute values) that do not have a place in
the CF data model are, nevertheless, stored within the appropriate
cfdm constructs. This allows them to be used when writing field
constructs to a new netCDF dataset, and also makes them accessible as
filters to a `Constructs` instance:

.. code-block:: python
   :caption: *Retrieve metadata constructs based on their netCDF
             names.*
	  
   >>> print(t.constructs.filter_by_ncvar('b'))
   Constructs:
   {'domainancillary1': <DomainAncillary: ncvar%b(1) >}
   >>> t.constructs('ncvar%x').value()
   <DimensionCoordinate: grid_longitude(9) degrees>
   >>> t.constructs('ncdim%x')
   <Constructs: domain_axis(1)>
     
Each construct has methods to access the netCDF elements which it
requires. For example, the field construct has the following methods:

===================================================  ======================================
Method                                               Description
===================================================  ======================================
`~Field.nc_get_variable`                             Return the netCDF variable name
`~Field.nc_set_variable`                             Set the netCDF variable name
`~Field.nc_del_variable`                             Remove the netCDF variable name
				                     
`~Field.nc_has_variable`                             Whether the netCDF variable name has
                                                     been set
				                     
`~Field.nc_global_attributes`                        Return the selection of properties to 
                                                     be written as netCDF global attributes
				                     
`~Field.nc_set_global_attribute`                     Set a property to be written as a
                                                     netCDF global attribute
					             
`~Field.nc_set_global_attributes`                    Set properties to be written as
                                                     netCDF global attributes
					             
`~Field.nc_clear_global_attributes`                  Clear the selection of properties
                                                     to be written as netCDF global
                                                     attributes
					             
`~Field.nc_group_attributes`                         Return the selection of properties to 
                                                     be written as netCDF group attributes
				                     
`~Field.nc_set_group_attribute`                      Set a property to be written as a
                                                     netCDF group attribute
					             
`~Field.nc_set_group_attributes`                     Set properties to be written as
                                                     netCDF group attributes
					             
`~Field.nc_clear_group_attributes`                   Clear the selection of properties
                                                     to be written as netCDF group
                                                     attributes
					             
`~Field.nc_variable_groups`                          Return the netCDF group structure
					             
`~Field.nc_set_variable_groups`                      Set the netCDF group structure
					             
`~Field.nc_clear_variable_groups`                    Remove the netCDF group structure
					             
`~Field.nc_geometry_variable_groups`                 Return the netCDF geometry
                                                     variable group structure
					             
`~Field.nc_set_geometry_variable_groups`             Set the netCDF geometry
                                                     variable group structure
					             
`~Field.nc_clear_geometry_variable_groups`           Remove the netCDF geometry
                                                     variable group structure
					             
`~Field.nc_del_component_variable`                   Remove the netCDF variable name for
                                                     all components of the given type.

`~Field.nc_set_component_variable`                   Set the netCDF variable name for all
                                                     components of the given type.

`~Field.nc_set_component_variable_groups`            Set the netCDF variable groups
                                                     hierarchy for all components of
						     the given type.

`~Field.nc_clear_component_variable_groups`          Remove the netCDF variable groups
                                                     hierarchy for all components of the
						     given type.

`~Field.nc_del_component_dimension`                  Remove the netCDF dimension name for
                                                     all components of the given type.

`~Field.nc_set_component_dimension`                  Set the netCDF dimension name for all
                                                     components of the given type.

`~Field.nc_set_component_dimension_groups`           Set the netCDF dimension groups
                                                     hierarchy for all components of the
						     given type.

`~Field.nc_clear_component_dimension_groups`         Remove the netCDF dimension groups
                                                     hierarchy for all components of the
						     given type.

`~Field.nc_del_component_sample_dimension`           Remove the netCDF sample dimension
                                                     name for all components of the given type.

`~Field.nc_set_component_sample_dimension`           Set the netCDF sample dimension name
                                                     for all components of the given type.

`~Field.nc_set_component_sample_dimension_groups`    Set the netCDF sample dimension
                                                     groups hierarchy for all components
						     of the given type.

`~Field.nc_clear_component_sample_dimension_groups`  Remove the netCDF sample dimension
                                                     groups hierarchy for all components
						     of the given type.
===================================================  ======================================

.. code-block:: python
   :caption: *Access netCDF elements associated with the field and
             metadata constructs.*

   >>> q.nc_get_variable()
   'q'
   >>> q.nc_global_attributes()
   {'project': None, 'Conventions': None}
   >>> q.nc_set_variable('humidity')
   >>> q.nc_get_variable()
   'humidity'
   >>> q.constructs('latitude').value().nc_get_variable()
   'lat'

The complete collection of netCDF interface methods is:

=============================================  =======================================  =====================================
Method                                         Classes                                  NetCDF element
=============================================  =======================================  =====================================
`!nc_del_variable`                             `Field`, `DimensionCoordinate`,          Variable name
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`,  `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       				
`!nc_get_variable`                             `Field`, `DimensionCoordinate`,          Variable name
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_has_variable`                             `Field`, `DimensionCoordinate`,          Variable name
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_set_variable`                             `Field`, `DimensionCoordinate`,          Variable name
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_variable_groups`                          `Field`, `DimensionCoordinate`,          Group hierarchy
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_set_variable_groups`                      `Field`, `DimensionCoordinate`,          Group hierarchy
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_clear_variable_groups`                    `Field`, `DimensionCoordinate`,          Group hierarchy
                                               `AuxiliaryCoordinate`, `CellMeasure`,
                                               `DomainAncillary`, `FieldAncillary`,
                                               `CoordinateReference`, `Bounds`,
			                       `Datum`, `Count`, `Index`, `List`
			                       
`!nc_del_dimension`                            `DomainAxis`, `Bounds`, `Count`,         Dimension name
                                               `Index`
			                       
`!nc_get_dimension`	                       `DomainAxis`, `Bounds`, `Count`,         Dimension name
                                               `Index`
			                       			                    
`!nc_has_dimension`	                       `DomainAxis`, `Bounds`, `Count`,         Dimension name
                                               `Index`
			                       			                    
`!nc_set_dimension`	                       `DomainAxis`, `Bounds`, `Count`,         Dimension name
                                               `Index`
			                       
`!nc_dimension_groups`                         `DomainAxis`, `Bounds`, `Count`,         Group hierarchy
                                               `Index`
			                       
`!nc_set_dimension_groups`	               `DomainAxis`, `Bounds`, `Count`,         Group hierarchy
                                               `Index`
			                       			                    
`!nc_clear_dimension_groups`	               `DomainAxis`, `Bounds`, `Count`,         Group hierarchy
                                               `Index`
				               
`!nc_is_unlimited`                             `DomainAxis`                             Unlimited dimension
				               
`!nc_set_unlimited` 	                       `DomainAxis`   	                        Unlimited dimension
				               
`!nc_global_attributes`	                       `Field`                                  Global attributes
			                       
`!nc_set_global_attribute`                     `Field`                                  Global attributes
			                       
`!nc_set_global_attributes`                    `Field`                                  Global attributes
			                       
`!nc_clear_global_attributes`                  `Field`                                  Global attributes
				               
`!nc_variable_groups`                          `Field`                                  Group hierarchy
 				               
`!nc_set_variable_groups`                      `Field`                                  Group hierarchy
 				               
`!nc_clear_variable_groups`                    `Field`                                  Group hierarchy
				               
`!nc_geometry_variable_groups`                 `Field`                                  Group hierarchy
 				               
`!nc_set_geometry_variable_groups`             `Field`                                  Group hierarchy
 				               
`!nc_clear_geometry_variable_groups`           `Field`                                  Group hierarchy
				               
`!nc_group_attributes`	                       `Field`                                  Group attributes
			                       
`!nc_set_group_attribute`                      `Field`                                  Group attributes
			                       
`!nc_set_group_attributes`                     `Field`                                  Group attributes
			                       
`!nc_clear_group_attributes`                   `Field`                                  Group attributes
			                       
`!nc_del_component_variable`                   `Field`                                  Component common netCDF properties

`!nc_set_component_variable`                   `Field`                                  Component common netCDF properties
					       
`!nc_set_component_variable_groups`            `Field`                                  Component common netCDF properties
					       
`!nc_clear_component_variable_groups`          `Field`                                  Component common netCDF properties
					       
`!nc_del_component_dimension`                  `Field`                                  Component common netCDF properties
					       
`!nc_set_component_dimension`                  `Field`                                  Component common netCDF properties
					       
`!nc_set_component_dimension_groups`           `Field`                                  Component common netCDF properties
					       
`!nc_clear_component_dimension_groups`         `Field`                                  Component common netCDF properties

`!nc_del_component_sample_dimension`           `Field`                                  Component common netCDF properties

`!nc_set_component_sample_dimension`           `Field`                                  Component common netCDF properties

`!nc_set_component_sample_dimension_groups`    `Field`                                  Component common netCDF properties

`!nc_clear_component_sample_dimension_groups`  `Field`                                  Component common netCDF properties

`!nc_get_external`                             `CellMeasure`                            External variable status
				               
`!nc_set_external`                             `CellMeasure`                            External variable status
			                       
`!nc_del_sample_dimension`                     `Count`, `Index`                         Sample dimension name
			                       
`!nc_get_sample_dimension`                     `Count`, `Index`                         Sample dimension name
    			                       
`!nc_has_sample_dimension`                     `Count`, `Index`                         Sample dimension name
			                       
`!nc_set_sample_dimension`                     `Count`, `Index`                         Sample dimension name
				               
`!nc_sample_dimension_groups`                  `Count`                                  Group hierarchy
 				               
`!nc_set_sample_dimension_groups`              `Count`                                  Group hierarchy
 				               
`!nc_clear_sample_dimension_groups`            `Count`                                  Group hierarchy

=============================================  =======================================  =====================================

----

.. _Writing-to-disk:
   
**Writing to disk**
-------------------

The `cfdm.write` function writes a field construct, or a sequence of
field constructs, to a netCDF file on disk:

.. code-block:: python
   :caption: *Write a field construct to a netCDF dataset on disk.*

   >>> print(q)
   Field: specific_humidity (ncvar%humidity)
   -----------------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   >>> cfdm.write(q, 'q_file.nc')

The file name may describe relative paths, and standard tilde and
shell parameter expansions are applied to it.

The new dataset is structured as follows:

.. code-block:: console
   :caption: *Inspect the new dataset with the ncdump command line
             tool.*

   $ ncdump -h q_file.nc
   netcdf q_file {
   dimensions:
   	lat = 5 ;
   	bounds2 = 2 ;
   	lon = 8 ;
   variables:
   	double lat_bnds(lat, bounds2) ;
   	double lat(lat) ;
   		lat:units = "degrees_north" ;
   		lat:standard_name = "latitude" ;
   		lat:bounds = "lat_bnds" ;
   	double lon_bnds(lon, bounds2) ;
   	double lon(lon) ;
   		lon:units = "degrees_east" ;
   		lon:standard_name = "longitude" ;
   		lon:bounds = "lon_bnds" ;
   	double time ;
   		time:units = "days since 2018-12-01" ;
   		time:standard_name = "time" ;
   	double humidity(lat, lon) ;
   		humidity:standard_name = "specific_humidity" ;
   		humidity:cell_methods = "area: mean" ;
   		humidity:units = "1" ;
   		humidity:coordinates = "time" ;
   
   // global attributes:
   		:Conventions = "CF-1.9" ;
   		:project = "research" ;
   }

A sequence of field constructs is written in exactly the same way:
   
.. code-block:: python
   :caption: *Write multiple field constructs to a netCDF dataset on
             disk.*
	     
   >>> x
   [<Field: specific_humidity(latitude(5), longitude(8)) 1>,
    <Field: air_temperature(atmosphere_hybrid_height_coordinate(1), grid_latitude(10), grid_longitude(9)) K>]
   >>> cfdm.write(x, 'new_file.nc')

By default the output file will be for CF-|version|.

The `cfdm.write` function has optional parameters to

* set the output netCDF format (all netCDF3 and netCDF4 formats are
  possible);

* append to the netCDF file rather than over-writing it by default;

* specify which field construct properties should become netCDF data
  variable attributes and which should become netCDF global
  attributes;
  
* set extra netCDF global attributes;
  
* create :ref:`external variables <External-variables>` in an external
  file;

* specify the version of the CF conventions (from CF-1.6 up to
  CF-|version|), and of any other conventions that the file adheres
  to;

* change the data type of output data arrays;
  
* apply netCDF compression and packing;

* set the endian-ness of the output data; and 

* specify whether or not :ref:`netCDF string arrays <Strings>` are to
  be used.

For example, to use the `mode` parameter to append a new field, or fields,
to a netCDF file whilst preserving the field or fields already contained
in that file:

.. code-block:: python
   :caption: *Append field constructs to a netCDF dataset on
             disk.*

   >>> g = cfdm.example_field(2)
   >>> cfdm.write(g, 'append-example-file.nc')
   >>> cfdm.read('append-example-file.nc')
   [<Field: air_potential_temperature(time(36), latitude(5), longitude(8)) K>]
   >>> h = cfdm.example_field(0)
   >>> h
   <Field: specific_humidity(latitude(5), longitude(8)) 1>
   >>> cfdm.write(h, 'append-example-file.nc', mode='a')
   >>> cfdm.read('append-example-file.nc')
   [<Field: air_potential_temperature(time(36), latitude(5), longitude(8)) K>,
    <Field: specific_humidity(latitude(5), longitude(8)) 1>]

Output netCDF variable and dimension names read from a netCDF dataset
are stored in the resulting field constructs, and may also be set
manually with the `!nc_set_variable`, `nc_set_dimension` and
`nc_set_sample_dimension` methods. If a name has not been set then one
will be generated internally (usually based on the standard name if it
exists).

It is possible to create netCDF unlimited dimensions using the
`~DomainAxis.nc_set_unlimited` method of the domain axis construct.

A field construct is not transformed through being written to a file
on disk and subsequently read back from that file.

.. code-block:: python
   :caption: *Read a file that has been created by writing a field
             construct, and compare the result with the original field
             construct in memory.*
	     
   >>> f = cfdm.read('q_file.nc')[0]
   >>> q.equals(f)
   True


.. _Global-attributes:

Global attributes
^^^^^^^^^^^^^^^^^

The field construct properties that correspond to the standardised
description-of-file-contents attributes are automatically written as
netCDF global attributes. Other attributes may also be written as
netCDF global attributes if they have been identified as such with the
*global_attributes* keyword, or via the
`~Field.nc_set_global_attribute` or `~Field.nc_set_global_attributes`
method of the field constructs. In either case, the creation of a
netCDF global attribute depends on the corresponding property values
being identical across all of the field constructs being written to
the file. If they are all equal then the property will be written as a
netCDF global attribute and not as an attribute of any netCDF data
variable; if any differ then the property is written only to each
netCDF data variable.

.. code-block:: python
   :caption: *Request that the "model" property is written as a netCDF
             global attribute, using the "global_attributes" keyword.*
	     
   >>> f.set_property('model', 'model_A')
   >>> cfdm.write(f, 'f_file.nc', global_attributes='model')

.. code-block:: python
   :caption: *Request that the "model" property is written as a netCDF
             global attribute, using the "nc_set_global_attribute"
             method.*
	     
   >>> f.nc_global_attributes()
   {'Conventions': None, 'project': None}
   >>> f.nc_set_global_attribute('model')
   >>> f.nc_global_attributes()
   {'Conventions': None, 'project': None, 'model': None}
   >>> f.nc_global_attributes(values=True)
   {'Conventions': 'CF-1.9', 'project': 'research', 'model': 'model_A'}
   >>> cfdm.write(f, 'f_file.nc')

It is possible to create both a netCDF global attribute and a netCDF
data variable attribute with the same name, but with different
values. This may be done by assigning the global value to the property
name with the `~Field.nc_set_global_attribute` or
`~Field.nc_set_global_attributes` method, or by via the
*file_descriptors* keyword. For the former technique, any
inconsistencies arising from multiple field constructs being written
to the same file will be resolved by omitting the netCDF global
attribute from the file.

.. code-block:: python
   :caption: *Request that the "information" property is written as
             netCDF global and data variable attributes, with
             different values, using the "nc_set_global_attribute"
             method.*
	     
   >>> f.set_property('information', 'variable information')
   >>> f.properties()
   {'Conventions': 'CF-1.9',
    'project': 'research',
    'standard_name': 'specific_humidity',
    'units': '1',
    'model': 'model_A',
    'information': 'variable information'}
   >>> f.nc_set_global_attribute('information', 'global information')
   >>> f.nc_global_attributes()
   {'Conventions': None,
    'project': None,
    'model': None,
    'information': 'global information'}
   >>> cfdm.write(f, 'f_file.nc')

NetCDF global attributes defined with the *file_descriptors* keyword
of the `cfdm.write` function will always be written as requested,
independently of the netCDF data variable attributes, and superseding
any global attributes that may have been defined with the
*global_attributes* keyword, or set on the individual field
constructs.

.. code-block:: python
   :caption: *Insist that the "history" property is written as netCDF
             a global attribute, with the "file_descriptors" keyword.*
	     
   >>> cfdm.write(f, 'f_file.nc', file_descriptors={'history': 'created in 2021'})
   >>> f_file = cfdm.read('f_file.nc')[0]
   >>> f_file.properties()
   {'Conventions': 'CF-1.9',
    'history': 'created in 2021',
    'model': 'model_A',
    'project': 'research',
    'information': 'variable information',
    'standard_name': 'specific_humidity',
    'units': '1'}
   >>> f_file.nc_global_attributes()
   {'Conventions': None,
    'history': None,
    'model': None,
    'project': None,
    'information': 'global information'}

.. _Conventions:

Conventions
^^^^^^^^^^^

The ``Conventions`` netCDF global attribute containing the version of
the CF conventions is always automatically created. If the version of
the CF conventions has been set as a field property, or with the
*Conventions* keyword of the `cfdm.write` function, then it is
ignored. However, other conventions that may apply can be set with
either technique.

.. code-block:: python
   :caption: *Two ways to add additional conventions to the
             "Conventions" netCDF global attribute.*
	     f_file_
   >>> f_file.set_property('Conventions', 'UGRID1.0')
   >>> cfdm.write(f, 'f_file.nc', Conventions='UGRID1.0')   

   
.. _Scalar-coordinate-variables:

Scalar coordinate variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A CF-netCDF scalar (i.e. zero-dimensional) coordinate variable is
created from a size one dimension coordinate construct that spans a
domain axis construct which is not spanned by the field construct's
data, nor the data of any other metadata construct. This occurs for
the field construct ``q``, for which the "time" dimension coordinate
construct was to the file ``q_file.nc`` as a scalar coordinate
variable.

To change this so that the "time" dimension coordinate construct is
written as a CF-netCDF size one coordinate variable, the field
construct's data must be expanded to span the corresponding size one
domain axis construct, by using the `~Field.insert_dimension` method
of the field construct.

.. code-block:: python
   :caption: *Write the "time" dimension coordinate construct to a
             (non-scalar) CF-netCDF coordinate variable by inserting
             the corresponding dimension into the field construct's
             data.*
		   
   >>> print(q)
   Field: specific_humidity (ncvar%humidity)
   -----------------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   <Field: specific_humidity(latitude(5), longitude(8)) 1>
   >>> key = q.construct_key('time')
   >>> axes = q.get_data_axes(key)
   >>> axes
   ('domainaxis2',)
   >>> q2 = q.insert_dimension(axis=axes[0])
   >>> q2
   <Field: specific_humidity(time(1), latitude(5), longitude(8)) 1>
   >>> cfdm.write(q2, 'q2_file.nc')

The new dataset is structured as follows (note, relative to file
``q_file.nc``, the existence of the "time" dimension and the lack of a
``coordinates`` attribute on the, now three-dimensional, data
variable):
   
.. code-block:: console
   :caption: *Inspect the new dataset with the ncdump command line
             tool.*

   $ ncdump -h q2_file.nc
   netcdf q2_file {
   dimensions:
   	lat = 5 ;
   	bounds2 = 2 ;
   	lon = 8 ;
   	time = 1 ;
   variables:
   	double lat_bnds(lat, bounds2) ;
   	double lat(lat) ;
   		lat:units = "degrees_north" ;
   		lat:standard_name = "latitude" ;
   		lat:bounds = "lat_bnds" ;
   	double lon_bnds(lon, bounds2) ;
   	double lon(lon) ;
   		lon:units = "degrees_east" ;
   		lon:standard_name = "longitude" ;
   		lon:bounds = "lon_bnds" ;
   	double time(time) ;
   		time:units = "days since 2018-12-01" ;
   		time:standard_name = "time" ;
   	double humidity(time, lat, lon) ;
   		humidity:units = "1" ;
   		humidity:standard_name = "specific_humidity" ;
   		humidity:cell_methods = "area: mean" ;
   
   // global attributes:
   		:Conventions = "CF-1.9" ;
   		:project = "research" ;
   }

.. _Compressed-constructs:

Compressed constructs
^^^^^^^^^^^^^^^^^^^^^

Constructs that contain compressed data will be automatically written
to a dataset with the correct compression encoding. See the section on
:ref:`compression <Compression>` for details.
   
.. _Strings:
  
Strings
^^^^^^^

String-valued data may be written to netCDF files either as netCDF
character arrays with a trailing dimension large enough to contain the
longest value, or as netCDF4 string arrays. The former is allowed for
all formats of netCDF3 and netCDF4 format files; but string arrays may
only be written to netCDF4 format files (note that string arrays can
not be written to netCDF4 classic format files).

By default, netCDF string arrays will be created whenever possible,
and in all other cases netCDF character arrays will be
used. Alternatively, netCDF character arrays can be used in all cases
by setting the *string* keyword of the `cfdm.write` function.

Groups
^^^^^^

NetCDF4 files with hierarchical groups may be created if a group
structure is defined by the netCDF variable and dimension names,
accessed via the :ref:`netCDF interface <NetCDF-interface>`.  See the
section on :ref:`hierarchical groups <Hierarchical-groups>` for
details.

----
      
.. _Hierarchical-groups:

**Hierarchical groups**
-----------------------

`Hierarchical groups`_ provide a mechanism to structure variables
within netCDF4 datasets, with well defined rules for resolving
references to out-of-group netCDF variables and dimensions.

A group structure that may be applied when writing to disk can be
created ab initio with the :ref:`netCDF interface
<NetCDF-interface>`. For example, the data variable and a coordinate
construct may be moved to a sub-group that has its own group
attribute, and a coordinate construct may be moved to a different
sub-group:

.. code-block:: python
   :caption: *Create a group structure and write it to disk.*

   >>> q, t = cfdm.read('file.nc')
   >>> print(q)
   Field: specific_humidity (ncvar%/forecast/model/q)
   --------------------------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   >>> q.set_property('comment', 'comment')
   >>> q.nc_set_group_attribute('comment', 'group comment')
   >>> q.nc_set_variable_groups(['forecast', 'model'])
   ()
   >>> q.construct('time').nc_set_variable_groups(['forecast'])
   ()
   >>> cfdm.write(q, 'grouped.nc')

.. code-block:: console
   :caption: *Inspect the new grouped dataset with the ncdump command
             line tool.*
   
   $ ncdump -h grouped.nc
   netcdf grouped {
   dimensions:
   	   lat = 5 ;
   	   bounds2 = 2 ;
   	   lon = 8 ;
   variables:
   	   double lat_bnds(lat, bounds2) ;
   	   double lat(lat) ;
   	   	   lat:units = "degrees_north" ;
   	   	   lat:standard_name = "latitude" ;
   	   	   lat:bounds = "lat_bnds" ;
   	   double lon_bnds(lon, bounds2) ;
   	   double lon(lon) ;
   	   	   lon:units = "degrees_east" ;
   	   	   lon:standard_name = "longitude" ;
   	   	   lon:bounds = "lon_bnds" ;
   
   // global attributes:
   		   :Conventions = "CF-1.8" ;
   		   :comment = "comment" ;
   
   group: forecast {
     variables:
     	   double time ;
  		   time:units = "days since 2018-12-01" ;
  		   time:standard_name = "time" ;

     group: model {
       variables:
       	   double q(lat, lon) ;
       		   q:project = "research" ;
       		   q:standard_name = "specific_humidity" ;
       		   q:units = "1" ;
       		   q:coordinates = "time" ;
       		   q:cell_methods = "area: mean" ;
   
       // group attributes:
       		   :comment = "group comment" ;
       } // group model
     } // group forecast
   }

When reading a netCDF dataset, the group structure and groups
attributes are recorded and are made accessible via the :ref:`netCDF
interface <NetCDF-interface>`.

.. code-block:: python
   :caption: *Read the grouped file and inspect its group structure.*

   >>> g = cfdm.read('grouped.nc')[0]
   >>> print(g)
   Field: specific_humidity (ncvar%/forecast/q)
   --------------------------------------------
   Data            : specific_humidity(latitude(5), longitude(8)) 1
   Cell methods    : area: mean
   Dimension coords: latitude(5) = [-75.0, ..., 75.0] degrees_north
                   : longitude(8) = [22.5, ..., 337.5] degrees_east
                   : time(1) = [2019-01-01 00:00:00]
   >>> g.nc_get_variable()
   '/forecast/model/q'
   >>> g.nc_variable_groups()
   ('forecast', 'model')
   >>> g.nc_group_attributes(values=True)
   {'comment': 'group comment'}
   >>> g.construct('latitude').nc_get_variable()
   'lat'
 
By default field constructs are written out to a dataset with their
groups struct (if any) intact. It is always possible, however, to
create a "flat" dataset, i.e. one without any sub-groups. This does
not require the removal of the group structure from the field
construct and all of its components (although that is possible), as it
can be done by simply by overriding the existing group structure by
setting the *group* keyword to `cfdm.write` to `False`.
   
.. code-block:: python
   :caption: *Write the field construct to a file with the same group
             structure, and also to a flat file.*

   >>> cfdm.write(g, 'flat.nc', group=False)

NetCDF variables in the flattened output file will inherit any netCDF
group attributes, providing that they are not superceded by variable
attributes. The output netCDF variable and dimension names will be
taken as the basenames of any that have been pre-defined. This is the
case in file ``flat.nc``, for which the netCDF variable ``q`` has
inherited the ``comment`` attribute that was originally set on the
``/forecast/model`` group. NetCDF group attributes may be set and
accessed via the :ref:`netCDF interface <NetCDF-interface>`, for both
netCDF variable and netCDF dimensions.

.. code-block:: console
   :caption: *Inspect the flat version of the dataset with the ncdump
             command line tool.*
   
   $ ncdump -h flat_out.nc
   netcdf flat {
   dimensions:
   	   lat = 5 ;
   	   bounds2 = 2 ;
   	   lon = 8 ;
   variables:
   	   double lat_bnds(lat, bounds2) ;
   	   double lat(lat) ;
   	   	   lat:units = "degrees_north" ;
   	   	   lat:standard_name = "latitude" ;
   	   	   lat:bounds = "lat_bnds" ;
   	   double lon_bnds(lon, bounds2) ;
   	   double lon(lon) ;
   	   	   lon:units = "degrees_east" ;
   	   	   lon:standard_name = "longitude" ;
   	   	   lon:bounds = "lon_bnds" ;
   	   double time ;
   	   	   time:units = "days since 2018-12-01" ;
   	   	   time:standard_name = "time" ;
   	   double q(lat, lon) ;
   	   	   q:comment = "group comment" ;
		   q:project = "research" ;
   	   	   q:standard_name = "specific_humidity" ;
   	   	   q:units = "1" ;
   	   	   q:coordinates = "time" ;
   	   	   q:cell_methods = "area: mean" ;
   		   
   // global attributes:
   		   :Conventions = "CF-1.8" ;
   		   :comment = "comment" ;
   }

The fields constructs read from a grouped file are identical to those
read from the flat version of the file:
   
.. code-block:: python
   :caption: *Demonstrate that the field constructs are independent of
             the dataset structure.*

   >>> f = cfdm.read('flat.nc')[0]
   >>> f.equals(g)
   True

----
   
.. _External-variables:

**External variables**
----------------------

`External variables`_ are those in a netCDF file that are referred to,
but which are not present in it. Instead, such variables are stored in
other netCDF files known as "external files". External variables may,
however, be incorporated into the field constructs of the dataset, as
if they had actually been stored in the same file, simply by providing
the external file names to the `cfdm.read` function.

An external variables file name may describe relative paths, and
standard tilde and shell parameter expansions are applied to it.

This is illustrated with the files ``parent.nc`` (found in the
:ref:`sample datasets <Sample-datasets>`):

.. code-block:: console
   :caption: *Inspect the parent dataset with the ncdump command line
             tool.*
   
   $ ncdump -h parent.nc
   netcdf parent {
   dimensions:
   	latitude = 10 ;
   	longitude = 9 ;
   variables:
   	double latitude(latitude) ;
   		latitude:units = "degrees_north" ;
   		latitude:standard_name = "latitude" ;
   	double longitude(longitude) ;
   		longitude:units = "degrees_east" ;
   		longitude:standard_name = "longitude" ;
   	double eastward_wind(latitude, longitude) ;
   		eastward_wind:units = "m s-1" ;
   		eastward_wind:standard_name = "eastward_wind" ;
   		eastward_wind:cell_measures = "area: areacella" ;
   
   // global attributes:
   		:Conventions = "CF-1.7" ;
   		:external_variables = "areacella" ;
   }

and ``external.nc`` (found in the :ref:`sample datasets
<Sample-datasets>`):

.. code-block:: console
   :caption: *Inspect the external dataset with the ncdump command
             line tool.*

   $ ncdump -h external.nc 
   netcdf external {
   dimensions:
   	latitude = 10 ;
   	longitude = 9 ;
   variables:
   	double areacella(longitude, latitude) ;
   		areacella:units = "m2" ;
   		areacella:standard_name = "cell_area" ;
   
   // global attributes:
   		:Conventions = "CF-1.7" ;
   }

The dataset in ``parent.nc`` may be read *without* specifying the
external file ``external.nc``. In this case a cell measure construct
is still created, but one without any metadata or data:

.. code-block:: python
   :caption: *Read the parent dataset without specifying the location
             of any external datasets.*

   >>> u = cfdm.read('parent.nc')[0]
   >>> print(u)
   Field: eastward_wind (ncvar%eastward_wind)
   ------------------------------------------
   Data            : eastward_wind(latitude(10), longitude(9)) m s-1
   Dimension coords: latitude(10) = [0.0, ..., 9.0] degrees
                   : longitude(9) = [0.0, ..., 8.0] degrees
   Cell measures   : measure:area (external variable: ncvar%areacella)

   >>> area = u.constructs('measure:area').value()
   >>> area
   <CellMeasure: measure:area >
   >>> area.nc_get_external()
   True
   >>> area.nc_get_variable()
   'areacella'
   >>> area.properties()
   {}
   >>> area.has_data()
   False

If this field construct were to be written to disk using `cfdm.write`,
then the output file would be identical to the original ``parent.nc``
file, i.e. the netCDF variable name of the cell measure construct
(``areacella``) would be listed by the ``external_variables`` global
attribute.

However, the dataset may also be read *with* the external file. In
this case a cell measure construct is created with all of the metadata
and data from the external file, as if the netCDF cell measure
variable had been present in the parent dataset:

.. code-block:: python
   :caption: *Read the parent dataset whilst providing the external
             dataset containing the external variables.*
   
   >>> g = cfdm.read('parent.nc', external='external.nc')[0]
   >>> print(g)
   Field: eastward_wind (ncvar%eastward_wind)
   ------------------------------------------
   Data            : eastward_wind(latitude(10), longitude(9)) m s-1
   Dimension coords: latitude(10) = [0.0, ..., 9.0] degrees
                   : longitude(9) = [0.0, ..., 8.0] degrees
   Cell measures   : cell_area(longitude(9), latitude(10)) = [[100000.5, ..., 100089.5]] m2
   >>> area = g.constructs('measure:area').value()
   >>> area
   <CellMeasure: cell_area(9, 10) m2>
   >>> area.nc_get_external()
   False
   >>> area.nc_get_variable()
   'areacella'
   >>> area.properties()
   {'standard_name': 'cell_area', 'units': 'm2'}
   >>> area.data
   <Data(9, 10): [[100000.5, ..., 100089.5]] m2>
   
If this field construct were to be written to disk using `cfdm.write`
then by default the cell measure construct, with all of its metadata
and data, would be written to the named output file, along with all of
the other constructs. There would be no ``external_variables`` global
attribute.

To create a reference to an external variable in an output netCDF
file, set the status of the cell measure construct to "external" with
its `~CellMeasure.nc_set_external` method.

.. code-block:: python
   :caption: *Flag the cell measure as external and write the field
             construct to a new file.*

   >>> area.nc_set_external(True)
   >>> cfdm.write(g, 'new_parent.nc')

To create a reference to an external variable in the an output netCDF
file and simultaneously create an external file containing the
variable, set the status of the cell measure construct to "external"
and provide an external file name to the `cfdm.write` function:

.. code-block:: python
   :caption: *Write the field construct to a new file and the cell
             measure construct to an external file.*

   >>> cfdm.write(g, 'new_parent.nc', external='new_external.nc')

.. _External-variables-with-cfdump:

External files with cfdump
^^^^^^^^^^^^^^^^^^^^^^^^^^

One or more external files may also be included with :ref:`cfdump
<cfdump>`.

.. code-block:: console
   :caption: *Use cfdump to describe the parent file without resolving
             the external variable reference.*
	     
   $ cfdump parent.nc 
   Field: eastward_wind (ncvar%eastward_wind)
   ------------------------------------------
   Data            : eastward_wind(latitude(10), longitude(9)) m s-1
   Dimension coords: latitude(10) = [0.0, ..., 9.0] degrees_north
                   : longitude(9) = [0.0, ..., 8.0] degrees_east
   Cell measures   : measure:area (external variable: ncvar%areacella)

.. code-block:: console
   :caption: *Providing an external file with the "-e" option allows
             the reference to be resolved.*
	     
   $ cfdump -e external.nc parent.nc 
   Field: eastward_wind (ncvar%eastward_wind)
   ------------------------------------------
   Data            : eastward_wind(latitude(10), longitude(9)) m s-1
   Dimension coords: latitude(10) = [0.0, ..., 9.0] degrees_north
                   : longitude(9) = [0.0, ..., 8.0] degrees_east
   Cell measures   : measure:area(longitude(9), latitude(10)) = [[100000.5, ..., 100089.5]] m2

----

.. _Compression:
   
**Compression**
---------------

The CF conventions have support for saving space by identifying
unwanted missing data.  Such compression techniques store the data
more efficiently and result in no precision loss. The CF data model,
however, views compressed arrays in their uncompressed form.

Therefore, the field construct contains :term:`domain axis constructs`
for the compressed dimensions and presents a view of compressed data
in its uncompressed form, even though the "underlying array" (i.e. the
actual array on disk or in memory that is contained in a `Data`
instance) is compressed. This means that the cfdm package includes
algorithms for uncompressing each type of compressed array.

There are two basic types of compression supported by the CF
conventions: ragged arrays (as used by :ref:`discrete sampling
geometries <Discrete-sampling-geometries>`) and :ref:`compression by
gathering <Gathering>`, each of which has particular implementation
details, but the following access patterns and behaviours apply to
both:

* Whether or not the data are compressed is tested with the
  `~Data.get_compression_type` method of the `Data` instance.

..

* Accessing the data via the `~Data.array` attribute of a `Data`
  instance returns a numpy array that is uncompressed. The underlying
  array will, however, remain in its compressed form. The compressed
  underlying array may be retrieved as a numpy array with the
  `~Data.compressed_array` attribute of the `Data` instance.

..

* A :ref:`subspace <Subspacing>` of a field construct is created with
  indices of the uncompressed form of the data. The new subspace will
  no longer be compressed, i.e. its underlying arrays will be
  uncompressed, but the original data will remain compressed. It
  follows that to uncompress all of the data in a field construct,
  index the field construct with (indices equivalent to) `Ellipsis`.
  
..

* If data elements are modified by :ref:`assigning <Assignment>` to
  indices of the uncompressed form of the data, then the compressed
  underlying array is replaced by its uncompressed form.

..

* An uncompressed field construct can be compressed, prior to being
  written to a dataset, with its `~Field.compress` method, which also
  compresses the metadata constructs as required.

..

* An compressed field construct can be uncompressed with its
  `~Field.uncompress` method, which also uncompresses the metadata
  constructs as required.

..

* If an underlying array is compressed at the time of writing to disk
  with the `cfdm.write` function, then it is written to the file as a
  compressed array, along with the supplementary netCDF variables and
  attributes that are required for the encoding. This means that if a
  dataset using compression is read from disk then it will be written
  back to disk with the same compression, unless data elements have
  been modified by assignment.

Examples of all of the above may be found in the sections on
:ref:`discrete sampling geometries <Discrete-sampling-geometries>` and
:ref:`gathering <Gathering>`.

.. _Discrete-sampling-geometries:
   
Discrete sampling geometries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Discrete sampling geometry (DSG)`_ features may be compressed by
combining them using one of three ragged array representations:
`contiguous`_, `indexed`_ or `indexed contiguous`_.

The count variable that is required to uncompress a contiguous, or
indexed contiguous, ragged array is stored in a `Count` instance and
is accessed with the `~Data.get_count` method of the `Data` instance.

The index variable that is required to uncompress an indexed, or
indexed contiguous, ragged array is stored in an `Index` instance and
is accessed with the `~Data.get_index` method of the `Data` instance.

The contiguous case is is illustrated with the file ``contiguous.nc``
(found in the :ref:`sample datasets <Sample-datasets>`):
     
.. code-block:: console
   :caption: *Inspect the compressed dataset with the ncdump command
             line tool.*
   
   $ ncdump -h contiguous.nc
   dimensions:
   	station = 4 ;
   	obs = 24 ;
   	strlen8 = 8 ;
   variables:
   	int row_size(station) ;
   		row_size:long_name = "number of observations for this station" ;
   		row_size:sample_dimension = "obs" ;
   	double time(obs) ;
   		time:units = "days since 1970-01-01 00:00:00" ;
   		time:standard_name = "time" ;
   	double lat(station) ;
   		lat:units = "degrees_north" ;
   		lat:standard_name = "latitude" ;
   	double lon(station) ;
   		lon:units = "degrees_east" ;
   		lon:standard_name = "longitude" ;
   	double alt(station) ;
   		alt:units = "m" ;
   		alt:positive = "up" ;
   		alt:standard_name = "height" ;
   		alt:axis = "Z" ;
   	char station_name(station, strlen8) ;
   		station_name:long_name = "station name" ;
   		station_name:cf_role = "timeseries_id" ;
   	double humidity(obs) ;
   		humidity:standard_name = "specific_humidity" ;
   		humidity:coordinates = "time lat lon alt station_name" ;
   		humidity:_FillValue = -999.9 ;
   
   // global attributes:
   		:Conventions = "CF-1.7" ;
   		:featureType = "timeSeries" ;
   }

Reading and inspecting this file shows the data presented in
two-dimensional uncompressed form, whilst the underlying array is
still in the one-dimension ragged representation described in the
file:

.. code-block:: python
   :caption: *Read a field construct from a dataset that has been
             compressed with contiguous ragged arrays, and inspect its
             data in uncompressed form.*
   
   >>> h = cfdm.read('contiguous.nc')[0]
   >>> print(h)
   Field: specific_humidity (ncvar%humidity)
   -----------------------------------------
   Data            : specific_humidity(ncdim%station(4), ncdim%timeseries(9))
   Dimension coords: 
   Auxiliary coords: time(ncdim%station(4), ncdim%timeseries(9)) = [[1969-12-29 00:00:00, ..., 1970-01-07 00:00:00]]
                   : latitude(ncdim%station(4)) = [-9.0, ..., 78.0] degrees_north
                   : longitude(ncdim%station(4)) = [-23.0, ..., 178.0] degrees_east
                   : height(ncdim%station(4)) = [0.5, ..., 345.0] m
                   : cf_role:timeseries_id(ncdim%station(4)) = [station1, ..., station4]
   >>> print(h.data.array)
   [[0.12 0.05 0.18   --   --   --   --   --   --]
    [0.05 0.11 0.2  0.15 0.08 0.04 0.06   --   --]
    [0.15 0.19 0.15 0.17 0.07   --   --   --   --]
    [0.11 0.03 0.14 0.16 0.02 0.09 0.1  0.04 0.11]]

.. code-block:: python
   :caption: *Inspect the underlying compressed array and the count
             variable that defines how to uncompress the data.*
	     
   >>> h.data.get_compression_type()
   'ragged contiguous'
   >>> print(h.data.compressed_array)
   [0.12 0.05 0.18 0.05 0.11 0.2 0.15 0.08 0.04 0.06 0.15 0.19 0.15 0.17 0.07
    0.11 0.03 0.14 0.16 0.02 0.09 0.1 0.04 0.11]
   >>> count_variable = h.data.get_count()
   >>> count_variable
   <Count: long_name=number of observations for this station(4) >
   >>> print(count_variable.data.array)
   [3 7 5 9]

The timeseries for the second station is easily selected by indexing
the "station" axis of the field construct:

.. code-block:: python
   :caption: *Get the data for the second station.*
	  
   >>> station2 = h[1]
   >>> station2
   <Field: specific_humidity(ncdim%station(1), ncdim%timeseries(9))>
   >>> print(station2.data.array)
   [[0.05 0.11 0.2 0.15 0.08 0.04 0.06 -- --]]

The underlying array of original data remains in compressed form until
data array elements are modified:
   
.. code-block:: python
   :caption: *Change an element of the data and show that the
             underlying array is no longer compressed.*

   >>> h.data.get_compression_type()
   'ragged contiguous'
   >>> h.data[1, 2] = -9
   >>> print(h.data.array)
   [[0.12 0.05 0.18   --   --   --   --   --   --]
    [0.05 0.11 -9.0 0.15 0.08 0.04 0.06   --   --]
    [0.15 0.19 0.15 0.17 0.07   --   --   --   --]
    [0.11 0.03 0.14 0.16 0.02 0.09 0.1  0.04 0.11]]
   >>> h.data.get_compression_type()
   ''

The easiest way to create a compressed field construct is to create
the equivalent uncompressed field construct and then compress it with
its `~Field.compress` method, which also compresses the metadata
constructs as required.
   
.. Code Block Start 3

.. code-block:: python
   :caption: *Create a field construct and then compress it.*

   import numpy
   import cfdm
   
   # Define the array values
   data = cfdm.Data([[280.0,   -99,   -99,   -99],
	             [281.0, 279.0, 278.0, 279.5]],
      	            mask=[[0, 1, 1, 1],
                          [0, 0, 0, 0]])
   	     
   # Create the field construct
   T = cfdm.Field()
   T.set_properties({'standard_name': 'air_temperature',
                     'units': 'K',
                     'featureType': 'timeSeries'})
   
   # Create the domain axis constructs
   X = T.set_construct(cfdm.DomainAxis(4))
   Y = T.set_construct(cfdm.DomainAxis(2))
   
   # Set the data for the field
   T.set_data(data, axes=[Y, X])
 
   # Compress the data 
   T.compress('contiguous',
              count_properties={'long_name': 'number of obs for this timeseries'},
              inplace=True)
		
.. Code Block End 3

The new field construct can now be inspected and written to a netCDF file:

.. code-block:: python
   :caption: *Inspect the new field construct and write it to disk.*
   
   >>> T
   <Field: air_temperature(key%domainaxis1(2), key%domainaxis0(4)) K>
   >>> print(T.data.array)
   [[280.0    --    --    --]
    [281.0 279.0 278.0 279.5]]
   >>> T.data.get_compression_type()
   'ragged contiguous'
   >>> print(T.data.compressed_array)
   [280.  281.  279.  278.  279.5]
   >>> count_variable = T.data.get_count()
   >>> count_variable
   <Count: long_name=number of obs for this timeseries(2) >
   >>> print(count_variable.data.array)
   [1 4]
   >>> cfdm.write(T, 'T_contiguous.nc')

The content of the new file is:
  
.. code-block:: console
   :caption: *Inspect the new compressed dataset with the ncdump
             command line tool.*   

   $ ncdump T_contiguous.nc
   netcdf T_contiguous {
   dimensions:
   	dim = 2 ;
   	element = 5 ;
   variables:
   	int64 count(dim) ;
   		count:long_name = "number of obs for this timeseries" ;
   		count:sample_dimension = "element" ;
   	float air_temperature(element) ;
   		air_temperature:units = "K" ;
   		air_temperature:standard_name = "air_temperature" ;
   
   // global attributes:
		:Conventions = "CF-1.7" ;
		:featureType = "timeSeries" ;
   data:
   
    count = 1, 4 ;
   
    air_temperature = 280, 281, 279, 278, 279.5 ;
   }

Exactly the same field construct may be also created explicitly with
underlying compressed data. A construct with an underlying ragged
array is created by initialising a `Data` instance with a ragged
array that is stored in one of three special array objects:
`RaggedContiguousArray`, `RaggedIndexedArray` or
`RaggedIndexedContiguousArray`.

.. Code Block Start 4

.. code-block:: python
   :caption: *Create a field construct with compressed data.*

   import numpy
   import cfdm
   
   # Define the ragged array values
   ragged_array = cfdm.Data([280, 281, 279, 278, 279.5])

   # Define the count array values
   count_array = [1, 4]

   # Create the count variable
   count_variable = cfdm.Count(data=cfdm.Data(count_array))
   count_variable.set_property('long_name', 'number of obs for this timeseries')

   # Create the contiguous ragged array object, specifying the
   # uncompressed shape
   array = cfdm.RaggedContiguousArray(
                    compressed_array=ragged_array,
                    shape=(2, 4), size=8, ndim=2,
                    count_variable=count_variable)

   # Create the field construct with the domain axes and the ragged
   # array
   T = cfdm.Field()
   T.set_properties({'standard_name': 'air_temperature',
                     'units': 'K',
                     'featureType': 'timeSeries'})
   
   # Create the domain axis constructs for the uncompressed array
   X = T.set_construct(cfdm.DomainAxis(4))
   Y = T.set_construct(cfdm.DomainAxis(2))
   
   # Set the data for the field
   T.set_data(cfdm.Data(array), axes=[Y, X])

.. Code Block End 4
   
.. _Gathering:

Gathering
^^^^^^^^^

`Compression by gathering`_ combines axes of a multidimensional array
into a new, discrete axis whilst omitting the missing values and thus
reducing the number of values that need to be stored.

The list variable that is required to uncompress a gathered array is
stored in a `List` object and is retrieved with the `~Data.get_list`
method of the `Data` instance.

This is illustrated with the file ``gathered.nc`` (found in the
:ref:`sample datasets <Sample-datasets>`):

.. code-block:: console
   :caption: *Inspect the compressed dataset with the ncdump command
             line tool.*
      
   $ ncdump -h gathered.nc
   netcdf gathered {
   dimensions:
   	time = 2 ;
   	lat = 4 ;
   	lon = 5 ;
   	landpoint = 7 ;
   variables:
   	double time(time) ;
   		time:standard_name = "time" ;
   		time:units = "days since 2000-1-1" ;
   	double lat(lat) ;
   		lat:standard_name = "latitude" ;
   		lat:units = "degrees_north" ;
   	double lon(lon) ;
   		lon:standard_name = "longitude" ;
   		lon:units = "degrees_east" ;
   	int landpoint(landpoint) ;
   		landpoint:compress = "lat lon" ;
   	double pr(time, landpoint) ;
   		pr:standard_name = "precipitation_flux" ;
   		pr:units = "kg m2 s-1" ;
   
   // global attributes:
   		:Conventions = "CF-1.7" ;
   }

Reading and inspecting this file shows the data presented in
three-dimensional uncompressed form, whilst the underlying array is
still in the two-dimensional gathered representation described in the
file:

.. code-block:: python
   :caption: *Read a field construct from a dataset that has been
             compressed by gathering, and inspect its data in
             uncompressed form.*

   >>> p = cfdm.read('gathered.nc')[0]
   >>> print(p)
   Field: precipitation_flux (ncvar%pr)
   ------------------------------------
   Data            : precipitation_flux(time(2), latitude(4), longitude(5)) kg m2 s-1
   Dimension coords: time(2) = [2000-02-01 00:00:00, 2000-03-01 00:00:00]
                   : latitude(4) = [-90.0, ..., -75.0] degrees_north
                   : longitude(5) = [0.0, ..., 40.0] degrees_east
   >>> print(p.data.array)
   [[[--       0.000122 0.0008   --       --      ]
     [0.000177 --       0.000175 0.00058  --      ]
     [--       --       --       --       --      ]
     [--       0.000206 --       0.0007   --      ]]
					  	 
    [[--       0.000202 0.000174 --       --      ]
     [0.00084  --       0.000201 0.0057   --      ]
     [--       --       --       --       --      ]
     [--       0.000223 --       0.000102 --      ]]]

.. code-block:: python
   :caption: *Inspect the underlying compressed array and the list
             variable that defines how to uncompress the data.*
	     
   >>> p.data.get_compression_type()
   'gathered'
   >>> print(p.data.compressed_array)
   [[0.000122 0.0008   0.000177 0.000175 0.00058 0.000206 0.0007  ]
    [0.000202 0.000174 0.00084  0.000201 0.0057  0.000223 0.000102]]
   >>> list_variable = p.data.get_list()
   >>> list_variable
   <List: ncvar%landpoint(7) >
   >>> print(list_variable.data.array)
   [1 2 5 7 8 16 18]

Subspaces based on the uncompressed axes of the field construct are
easily created:

.. code-block:: python
   :caption: *Get subspaces based on indices of the uncompressed
             data.*
	  
   >>> p[0]
   <Field: precipitation_flux(time(1), latitude(4), longitude(5)) kg m2 s-1>
   >>> p[1, :, 3:5]
   <Field: precipitation_flux(time(1), latitude(4), longitude(2)) kg m2 s-1>

The underlying array of original data remains in compressed form until
data array elements are modified:
   
.. code-block:: python
   :caption: *Change an element of the data and show that the
             underlying array is no longer compressed.*

   >>> p.data.get_compression_type()
   'gathered'
   >>> p.data[1] = -9
   >>> p.data.get_compression_type()
   ''
   
A construct with an underlying gathered array is created by
initializing a `Data` instance with a gathered array that is stored in
the special `GatheredArray` array object. The following code creates a
simple field construct with an underlying gathered array:

.. Code Block Start 5

.. code-block:: python
   :caption: *Create a field construct with compressed data.*

   import numpy	  
   import cfdm

   # Define the gathered values
   gathered_array = cfdm.Data([[2, 1, 3], [4, 0, 5]])

   # Define the list array values
   list_array = [1, 4, 5]

   # Create the list variable
   list_variable = cfdm.List(data=cfdm.Data(list_array))

   # Create the gathered array object, specifying the uncompressed
   # shape
   array = cfdm.GatheredArray(
                    compressed_array=gathered_array,
		    compressed_dimension=1,
                    shape=(2, 3, 2), size=12, ndim=3,
                    list_variable=list_variable)

   # Create the field construct with the domain axes and the gathered
   # array
   P = cfdm.Field(properties={'standard_name': 'precipitation_flux',
                              'units': 'kg m-2 s-1'})

   # Create the domain axis constructs for the uncompressed array
   T = P.set_construct(cfdm.DomainAxis(2))
   Y = P.set_construct(cfdm.DomainAxis(3))
   X = P.set_construct(cfdm.DomainAxis(2))

   # Set the data for the field
   P.set_data(cfdm.Data(array), axes=[T, Y, X])

.. Code Block End 5

Note that, because compression by gathering acts on a subset of the
array dimensions, it is necessary to state the position of the
compressed dimension in the compressed array (with the
``compressed_dimension`` parameter of the `GatheredArray`
initialisation).

The new field construct can now be inspected and written a netCDF file:

.. code-block:: python
   :caption: *Inspect the new field construct and write it to disk.*
   
   >>> P
   <Field: precipitation_flux(key%domainaxis0(2), key%domainaxis1(3), key%domainaxis2(2)) kg m-2 s-1>
   >>> print(P.data.array)
   [[[ -- 2.0]
     [ --  --]
     [1.0 3.0]]

    [[ -- 4.0]
     [ --  --]
     [0.0 5.0]]]
   >>> P.data.get_compression_type()
   'gathered'
   >>> print(P.data.compressed_array)
   [[2. 1. 3.]
    [4. 0. 5.]]
   >>> list_variable = P.data.get_list()
   >>> list_variable 
   <List: (3) >
   >>> print(list_variable.data.array)
   [1 4 5]
   >>> cfdm.write(P, 'P_gathered.nc')

The content of the new file is:
   
.. code-block:: console
   :caption: *Inspect new the compressed dataset with the ncdump
             command line tool.*
   
   $ ncdump P_gathered.nc
   netcdf P_gathered {
   dimensions:
   	dim = 2 ;
   	dim_1 = 3 ;
   	dim_2 = 2 ;
   	list = 3 ;
   variables:
   	int64 list(list) ;
   		list:compress = "dim_1 dim_2" ;
   	float precipitation_flux(dim, list) ;
   		precipitation_flux:units = "kg m-2 s-1" ;
   		precipitation_flux:standard_name = "precipitation_flux" ;
   
   // global attributes:
   		:Conventions = "CF-1.7" ;
   data:
   
    list = 1, 4, 5 ;
   
    precipitation_flux =
     2, 1, 3,
     4, 0, 5 ;
   }


.. _Coordinate-subampling:

Coordinate subsampling
^^^^^^^^^^^^^^^^^^^^^^

`Lossy compression by coordinate subsampling`_ was introduced into the
CF conventions at CF-1.9, but is not yet available in cfdm. It will be
ready in a future 1.9.x.0 release.

----

.. _Controlling-output-messages:

**Controlling output messages**
-------------------------------

cfdm will produce messages upon the execution of operations, to
provide feedback about:

* the progress of, and under-the-hood steps involved in, the
  operations it is performing;
* the events that emerge during these operations;
* the nature of the dataset being operated on, including CF compliance
  issues that may be encountered during the operation.

This feedback may be purely informational, or may convey warning(s)
about dataset issues or the potential for future error(s).

It is possible to configure the extent to which messages are output at
runtime, i.e. the verbosity of cfdm, so that less serious and/or more
detailed messages can be filtered out.

There are two means to do this, which are covered in more detail in
the sub-sections below. Namely, you may configure the extent of
messaging:

* **globally** i.e. for all cfdm operations, by setting the
  `cfdm.log_level` which controls the project-wide logging;
  
* **for a specific function only** (for many functions) by setting
  that function's *verbose* keyword (which overrides the global
  setting for the duration of the function call).

Both possibilities use a consistent level-based cut-off system, as
detailed below.

.. _Logging:

Logging
^^^^^^^

All messages from cfdm, excluding exceptions which are always raised
in error cases, are incorporated into a logging system which assigns
to each a level based on the relative seriousness and/or
verbosity. From highest to lowest on this scale, these levels are:

* ``'WARNING'``: conveys a warning;
* ``'INFO'``: provides information concisely, in a few lines or so;
* ``'DETAIL'``: provides information in a more detailed manner than
  ``'INFO'``;
* ``'DEBUG'``: produces highly-verbose information intended mainly for
  the purposes of debugging and cfdm library development.

The function `cfdm.log_level` sets the minimum of these levels for
which messages are displayed. Any message marked as being of any lower
level will be filtered out. Note it sets the verbosity *globally*, for
*all* cfdm library operations (unless these are overridden for
individual functions, as covered below).

As well as the named log levels above, `cfdm.log_level` accepts a
further identifier, ``'DISABLE'``. Each of these potential settings
has a numerical value that is treated interchangeably and may instead
be set (as this may be easier to recall and write, if less
explicit). The resulting behaviour in each case is as follows:

=======================  ============  =========================================
Log level                Integer code  Result when set as the log severity level
=======================  ============  =========================================
``'DISABLE'``            ``0``         *Disable all* logging messages. Note this
                                       does not include exception messages
                                       raised by errors.
				       
``'WARNING'`` (default)  ``1``         *Only show* logging messages that are
                                       *warnings* (those labelled as
				       ``'WARNING'``).
				       
``'INFO'``               ``2``         *Only show* logging messages that are
                                       *warnings or concise informational
                                       messages* (marked as ``'WARNING'`` or
				       ``'INFO'`` respectively).
				       
``'DETAIL'``             ``3``         *Enable all* logging messages *except 
                                       for debugging messages*. In other words,
                                       show logging messages labelled
                                       ``'WARNING'``, ``'INFO'`` and
				       ``'DETAIL'``, but not ``'DEBUG'``.
				       
``'DEBUG'``              ``-1``        *Enable all* logging messages,
                                       *including debugging messages*
				       (labelled as ``'DEBUG'``).
=======================  ============  =========================================

Note ``'DEBUG'`` is intended as a special case for debugging, which
should not be required in general usage of cfdm, hence its equivalence
to ``-1`` rather than ``4`` which would follow the increasing integer
code pattern.  ``-1`` reflects that it is the final value in the
sequence, as with Python indexing.

The default value for `cfdm.log_level` is ``'WARNING'`` (``1``).
However, whilst completing this tutorial, it may be instructive to set
the log level` to a higher verbosity level such as ``'INFO'`` to gain
insight into the internal workings of cfdm calls.


.. _Function-verbosity:

Function verbosity
^^^^^^^^^^^^^^^^^^

Functions and methods that involve a particularly high number of steps
or especially complex processing, for example the `cfdm.read` and
`cfdm.write` functions, accept a keyword argument *verbose*. This be
set to change the minimum log level at which messages are displayed
for the function/method call only, without being influenced by, or
influencing, the global `cfdm.log_level` value.

A *verbose* value effectively overrides the value of `cfdm.log_level`
for the function/method along with any functions/methods it calls in
turn, until the origin function/method completes.

The *verbose* argument accepts the same levels as `cfdm.log_level`
(including ``0`` for ``'DISABLE'``), as described in :ref:`the logging
section <logging>`, namely either an integer or a corresponding string
for example ``verbose=2`` or equivalently ``verbose='INFO'``
(or ``verbose='info'`` since case is ignored).

By default, *verbose* is set to `None`, in which case the value of the
`cfdm.log_level` setting is used to determine which messages,
if any, are filtered out.


.. rubric:: Footnotes

----

.. [#dap] Requires the netCDF4 python package to have been built with
          OPeNDAP support enabled. See
          http://unidata.github.io/netcdf4-python for details.

.. [#caveat] `Lossy compression by coordinate subsampling`_ was
             introduced into the CF conventions at CF-1.9, but is not
             yet available in cfdm. It will be ready in a future
             1.9.x.0 release.

.. External links to the CF conventions (will need updating with new versions of CF)
   
.. _External variables:                          http://cfconventions.org/cf-conventions/cf-conventions.html#external-variables
.. _Discrete sampling geometry (DSG):            http://cfconventions.org/cf-conventions/cf-conventions.html#discrete-sampling-geometries
.. _incomplete multidimensional form:            http://cfconventions.org/cf-conventions/cf-conventions.html#_incomplete_multidimensional_array_representation
.. _Compression by gathering:                    http://cfconventions.org/cf-conventions/cf-conventions.html#compression-by-gathering
.. _contiguous:                                  http://cfconventions.org/cf-conventions/cf-conventions.html#_contiguous_ragged_array_representation
.. _indexed:                                     http://cfconventions.org/cf-conventions/cf-conventions.html#_indexed_ragged_array_representation
.. _indexed contiguous:                          http://cfconventions.org/cf-conventions/cf-conventions.html#_ragged_array_representation_of_time_series_profiles
.. _geometries:                                  http://cfconventions.org/cf-conventions/cf-conventions.html#geometries
.. _Hierarchical groups:                         http://cfconventions.org/cf-conventions/cf-conventions.html#groups
.. _Lossy compression by coordinate subsampling: http://cfconventions.org/cf-conventions/cf-conventions.html#compression-by-coordinate-subsampling
.. currentmodule:: cfdm
.. default-role:: obj

################
**cfdm package**
################

----

.. include:: introduction.rst

**Contents**
============

----

.. toctree::
   :maxdepth: 1

   introduction
   cf_data_model
   installation
   contributing
   tutorial
   api_reference
   philosophy
   performance
   extensions
   releases
   Changelog
   
**Index and search**
====================

----

* :ref:`genindex`
* :ref:`Search <search>`
.. currentmodule:: cfdm
.. default-role:: obj

.. _class_extended:

**cfdm classes**
================

----

Version |release| for version |version| of the CF conventions.


Field construct class
---------------------

.. autosummary::
   :nosignatures:
   :toctree: class/
		 
   cfdm.Field

Domain construct class
----------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Domain

Metadata construct classes
--------------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.AuxiliaryCoordinate
   cfdm.CellMeasure
   cfdm.CellMethod
   cfdm.CoordinateReference
   cfdm.DimensionCoordinate
   cfdm.DomainAncillary
   cfdm.DomainAxis
   cfdm.FieldAncillary
  
Constructs class
----------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Constructs

Coordinate component classes
----------------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Bounds
   cfdm.CoordinateConversion
   cfdm.Datum
   cfdm.InteriorRing

Data classes
------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Data
   cfdm.NetCDFArray
   cfdm.NumpyArray
   cfdm.Array

Data compression classes
------------------------

Classes that support the creation and storage of compressed arrays.

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Count
   cfdm.Index
   cfdm.List
   cfdm.GatheredArray
   cfdm.RaggedContiguousArray
   cfdm.RaggedIndexedArray
   cfdm.RaggedIndexedContiguousArray
   cfdm.CompressedArray


Miscellaneous classes
---------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.Constant
   cfdm.Configuration

.. currentmodule:: cfdm
.. default-role:: obj

.. raw:: html

    <style> .small {font-size:small} </style>

.. role:: small

**Introduction**
================

----

Version |release| for version |version| of the CF conventions.

.. contents::
   :local:
   :backlinks: entry

The cfdm library implements the data model of the CF (Climate and
Forecast) metadata conventions (http://cfconventions.org) and so
should be able to represent and manipulate all existing and
conceivable CF-compliant datasets.

The CF conventions are designed to promote the creation, processing,
and sharing of climate and forecasting data using Network Common Data
Form (netCDF) files and libraries
(https://www.unidata.ucar.edu/software/netcdf). They cater for data
from model simulations as well as from observations, made in situ or
by remote sensing platforms, of the planetary surface, ocean, and
atmosphere. For a netCDF data variable, they provide a description of
the physical meaning of data and of its spatial, temporal, and other
dimensional properties. The CF data model is an abstract
interpretation of the CF conventions that is independent of the netCDF
encoding.

For more details see *cfdm: A Python reference implementation of the
CF data model* in the Journal of Open Source Software:
https://doi.org/10.21105/joss.02717

----
    
**Functionality**
-----------------

The cfdm library can create field constructs ab initio, or read them
from netCDF files, inspect, subspace and modify in memory, and write
them to CF-netCDF dataset files. As long as it can interpret the data,
cfdm does not enforce CF-compliance, allowing non-compliant datasets
to be read, processed, corrected and rewritten.

It does not contain higher-level analysis functions (such as
regridding) because the expectation is that other libraries will build
on cfdm, inheriting its comprehensive knowledge of the CF conventions,
to add more sophisticated methods.

.. code-block:: python
   :caption: *A simple example of reading a field construct from a
             file and inspecting it.*

   >>> import cfdm
   >>> f = cfdm.read('file.nc')
   >>> f
   [<Field: air_temperature(time(12), latitude(64), longitude(128)) K>]
   >>> print(f[0])
   Field: air_temperature (ncvar%tas)
   ----------------------------------
   Data            : air_temperature(time(12), latitude(64), longitude(128)) K
   Cell methods    : time(12): mean (interval: 1.0 month)
   Dimension coords: time(12) = [0450-11-16 00:00:00, ..., 0451-10-16 12:00:00] noleap
                   : latitude(64) = [-87.8638, ..., 87.8638] degrees_north
                   : longitude(128) = [0.0, ..., 357.1875] degrees_east
                   : height(1) = [2.0] m

The cfdm package can

* read :term:`field constructs <field construct>` and :term:`domain
  constructs <domain construct>` from netCDF and CDL datasets,

* create new field and domain constructs in memory,

* write field and domain constructs to netCDF datasets on disk,

* read, write, and create coordinates defined by geometry cells,

* read and write netCDF4 string data-type variables,

* read, write, and create netCDF and CDL datasets containing
  hierarchical groups,

* inspect field and domain constructs,

* test whether two constructs are the same,

* modify field and domain construct metadata and data,

* create subspaces of field and domain constructs,

* incorporate, and create, metadata stored in external files, and

* read, write, and create data that have been compressed by convention
  (i.e. ragged or gathered arrays), whilst presenting a view of the
  data in its uncompressed form.

Note that the cfdm package enables the representation and creation of
CF field constructs, but it is largely :ref:`up to the user to use
them in a CF-compliant way <CF-conventions>`.

A command line tool is provided that allows inspection of datasets
outside of a Python environment:

.. code-block:: console
   :caption: *Inspect a dataset from the command line.*

   $ cfdump file.nc
   Field: air_temperature (ncvar%tas)
   ----------------------------------
   Data            : air_temperature(time(12), latitude(64), longitude(128)) K
   Cell methods    : time(12): mean (interval: 1.0 month)
   Dimension coords: time(12) = [0450-11-16 00:00:00, ..., 0451-10-16 12:00:00] noleap
                   : latitude(64) = [-87.8638, ..., 87.8638] degrees_north
                   : longitude(128) = [0.0, ..., 357.1875] degrees_east
                   : height(1) = [2.0] m

----

**Related packages**
--------------------

The `cf-python <https://ncas-cms.github.io/cf-python>`_ package, which
is built as an extension to cfdm, includes higher-level functionality,
such as regridding, and statistical operations. In turn, the `cf-plot
<http://ajheaps.github.io/cf-plot/>`_ package provides comprehensive
visualisation of field constructs created by cf-python.

----

**Citation**
------------

If you use cfdm, either as a stand-alone application or to provide a
CF data model implementation to another software library, please
consider including the reference:

Hassell, D., and Bartholomew, S. L. (2020). cfdm: A Python reference
  implementation of the CF data model. Journal of Open Source
  Software, 5(54), 2717, https://doi.org/10.21105/joss.02717

.. code-block:: bibtex
   
   @article{Hassell2020,
     doi = {10.21105/joss.02717},
     url = {https://doi.org/10.21105/joss.02717},
     year = {2020},
     publisher = {The Open Journal},
     volume = {5},
     number = {54},
     pages = {2717},
     author = {David Hassell and Sadie L. Bartholomew},
     title = {cfdm: A Python reference implementation of the CF data model},
     journal = {Journal of Open Source Software}
   }

----

**References**
--------------

Eaton, B., Gregory, J., Drach, B., Taylor, K., Hankin, S., Caron, J.,
  Signell, R., et al. (2020). NetCDF Climate and Forecast (CF)
  Metadata Conventions. CF Conventions Committee. Retrieved from
  https://cfconventions.org/cf-conventions/cf-conventions.html

Hassell, D., and Bartholomew, S. L. (2020). cfdm: A Python reference
  implementation of the CF data model. Journal of Open Source
  Software, 5(54), 2717, https://doi.org/10.21105/joss.02717

Hassell, D., Gregory, J., Blower, J., Lawrence, B. N., and
  Taylor, K. E. (2017). A data model of the Climate and Forecast
  metadata conventions (CF-1.6) with a software implementation
  (cf-python v2.1), Geosci. Model Dev., 10, 4619-4646,
  https://doi.org/10.5194/gmd-10-4619-2017

Rew, R., and Davis, G. (1990). NetCDF: An Interface for Scientific
  Data Access. IEEE Computer Graphics and Applications, 10(4),
  76–82. https://doi.org/10.1109/38.56302

Rew, R., Hartnett, E., and Caron, J. (2006). NetCDF-4: Software
  Implementing an Enhanced Data Model for the Geosciences. In 22nd
  International Conference on Interactive Information Processing
  Systems for Meteorology, Oceanography, and Hydrology. AMS. Retrieved
  from
  https://www.unidata.ucar.edu/software/netcdf/papers/2006-ams.pdf
.. currentmodule:: cfdm
.. default-role:: obj

.. _class_core:

**cfdm.core classes**
=====================

----

Version |release| for version |version| of the CF conventions.


Field construct class
---------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Field

Domain construct class
----------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Domain
   
Metadata construct classes
--------------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.AuxiliaryCoordinate
   cfdm.core.CellMeasure
   cfdm.core.CellMethod
   cfdm.core.CoordinateReference
   cfdm.core.DimensionCoordinate
   cfdm.core.DomainAncillary
   cfdm.core.DomainAxis
   cfdm.core.FieldAncillary

Constructs class
----------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Constructs

Coordinate component classes
----------------------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Bounds
   cfdm.core.CoordinateConversion
   cfdm.core.Datum
   cfdm.core.InteriorRing

Data classes
------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Data
   cfdm.core.NumpyArray
   cfdm.core.Array

Abstract base classes
---------------------

Abstract base classes that provide the basis for constructs and
construct components.

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.Container
   cfdm.core.Properties
   cfdm.core.PropertiesData
   cfdm.core.PropertiesDataBounds
   cfdm.core.Coordinate
   cfdm.core.Parameters
   cfdm.core.ParametersDomainAncillaries

Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: class/

   cfdm.core.DocstringRewriteMeta
   
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.RaggedIndexedArray
=======================

----

.. autoclass:: cfdm.RaggedIndexedArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.RaggedIndexedArray.get_compressed_axes
   ~cfdm.RaggedIndexedArray.get_compressed_dimension
   ~cfdm.RaggedIndexedArray.get_compression_type
   ~cfdm.RaggedIndexedArray.get_index

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.RaggedIndexedArray.array
   ~cfdm.RaggedIndexedArray.compressed_array
   ~cfdm.RaggedIndexedArray.dtype
   ~cfdm.RaggedIndexedArray.ndim
   ~cfdm.RaggedIndexedArray.shape
   ~cfdm.RaggedIndexedArray.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedIndexedArray.copy
   ~cfdm.RaggedIndexedArray.get_subspace
   ~cfdm.RaggedIndexedArray.source
   ~cfdm.RaggedIndexedArray.to_memory

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedIndexedArray.__array__
   ~cfdm.RaggedIndexedArray.__deepcopy__
   ~cfdm.RaggedIndexedArray.__getitem__
   ~cfdm.RaggedIndexedArray.__repr__
   ~cfdm.RaggedIndexedArray.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.RaggedIndexedArray._docstring_special_substitutions
   ~cfdm.RaggedIndexedArray._docstring_substitutions        
   ~cfdm.RaggedIndexedArray._docstring_package_depth        
   ~cfdm.RaggedIndexedArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-CellMeasure:

cfdm.CellMeasure
================

----

.. autoclass:: cfdm.CellMeasure
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.dump
   ~cfdm.CellMeasure.identity
   ~cfdm.CellMeasure.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.CellMeasure.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst


   ~cfdm.CellMeasure.del_measure
   ~cfdm.CellMeasure.get_measure
   ~cfdm.CellMeasure.has_measure
   ~cfdm.CellMeasure.set_measure
   ~cfdm.CellMeasure.del_property
   ~cfdm.CellMeasure.get_property
   ~cfdm.CellMeasure.has_property
   ~cfdm.CellMeasure.set_property
   ~cfdm.CellMeasure.properties
   ~cfdm.CellMeasure.clear_properties
   ~cfdm.CellMeasure.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.apply_masking
   ~cfdm.CellMeasure.del_data
   ~cfdm.CellMeasure.get_data
   ~cfdm.CellMeasure.has_data
   ~cfdm.CellMeasure.set_data
   ~cfdm.CellMeasure.insert_dimension
   ~cfdm.CellMeasure.squeeze
   ~cfdm.CellMeasure.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.CellMeasure.data
   ~cfdm.CellMeasure.dtype
   ~cfdm.CellMeasure.ndim
   ~cfdm.CellMeasure.shape
   ~cfdm.CellMeasure.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.copy
   ~cfdm.CellMeasure.creation_commands
   ~cfdm.CellMeasure.equals
   ~cfdm.CellMeasure.has_bounds
   ~cfdm.CellMeasure.uncompress
   ~cfdm.CellMeasure.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.nc_del_variable
   ~cfdm.CellMeasure.nc_get_variable
   ~cfdm.CellMeasure.nc_has_variable
   ~cfdm.CellMeasure.nc_set_variable
   ~cfdm.CellMeasure.nc_get_external
   ~cfdm.CellMeasure.nc_set_external

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.nc_variable_groups
   ~cfdm.CellMeasure.nc_clear_variable_groups
   ~cfdm.CellMeasure.nc_set_variable_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMeasure.__deepcopy__
   ~cfdm.CellMeasure.__getitem__
   ~cfdm.CellMeasure.__repr__
   ~cfdm.CellMeasure.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.CellMeasure._docstring_special_substitutions
   ~cfdm.CellMeasure._docstring_substitutions        
   ~cfdm.CellMeasure._docstring_package_depth        
   ~cfdm.CellMeasure._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Constructs
===============

----

.. autoclass:: cfdm.Constructs
   :no-members:
   :no-inherited-members:

Filtering
---------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.filter_by_identity
   ~cfdm.Constructs.filter_by_property
   ~cfdm.Constructs.filter_by_measure
   ~cfdm.Constructs.filter_by_method
   ~cfdm.Constructs.filter_by_axis
   ~cfdm.Constructs.filter_by_naxes
   ~cfdm.Constructs.filter_by_size
   ~cfdm.Constructs.filter_by_data
   ~cfdm.Constructs.filter_by_type
   ~cfdm.Constructs.filter_by_key
   ~cfdm.Constructs.filter_by_ncdim
   ~cfdm.Constructs.filter_by_ncvar
   ~cfdm.Constructs.filter
   ~cfdm.Constructs.domain_axes
   ~cfdm.Constructs.filters_applied
   ~cfdm.Constructs.clear_filters_applied
   ~cfdm.Constructs.inverse_filter
   ~cfdm.Constructs.unfilter

Constructs and identifiers
--------------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.key
   ~cfdm.Constructs.value
   ~cfdm.Constructs.construct_type
   ~cfdm.Constructs.construct_types
   ~cfdm.Constructs.domain_axis_identity
   ~cfdm.Constructs.new_identifier

Data axes
---------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.data_axes
   ~cfdm.Constructs.get_data_axes

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.copy
   ~cfdm.Constructs.shallow_copy
   ~cfdm.Constructs.equals
   ~cfdm.Constructs.replace
   ~cfdm.Constructs.todict
   ~cfdm.Constructs.ordered

Dictionary-access methods
-------------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.get
   ~cfdm.Constructs.items
   ~cfdm.Constructs.keys
   ~cfdm.Constructs.values
   ~cfdm.Constructs.__getitem__

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constructs.__call__
   ~cfdm.Constructs.__contains__
   ~cfdm.Constructs.__copy__
   ~cfdm.Constructs.__deepcopy__
   ~cfdm.Constructs.__getitem__
   ~cfdm.Constructs.__iter__
   ~cfdm.Constructs.__len__
   ~cfdm.Constructs.__repr__
   ~cfdm.Constructs.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Constructs._docstring_special_substitutions
   ~cfdm.Constructs._docstring_substitutions        
   ~cfdm.Constructs._docstring_package_depth        
   ~cfdm.Constructs._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-DomainAncillary:

cfdm.DomainAncillary
====================

----

.. autoclass:: cfdm.DomainAncillary
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.dump
   ~cfdm.DomainAncillary.identity  
   ~cfdm.DomainAncillary.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DomainAncillary.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.del_property
   ~cfdm.DomainAncillary.get_property
   ~cfdm.DomainAncillary.has_property
   ~cfdm.DomainAncillary.set_property
   ~cfdm.DomainAncillary.properties
   ~cfdm.DomainAncillary.clear_properties
   ~cfdm.DomainAncillary.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.apply_masking
   ~cfdm.DomainAncillary.del_data
   ~cfdm.DomainAncillary.get_data
   ~cfdm.DomainAncillary.has_data
   ~cfdm.DomainAncillary.set_data
   ~cfdm.DomainAncillary.insert_dimension
   ~cfdm.DomainAncillary.squeeze
   ~cfdm.DomainAncillary.transpose
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DomainAncillary.data
   ~cfdm.DomainAncillary.dtype
   ~cfdm.DomainAncillary.ndim
   ~cfdm.DomainAncillary.shape
   ~cfdm.DomainAncillary.size

Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.del_bounds
   ~cfdm.DomainAncillary.get_bounds
   ~cfdm.DomainAncillary.has_bounds
   ~cfdm.DomainAncillary.set_bounds
   ~cfdm.DomainAncillary.get_bounds_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DomainAncillary.bounds

Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.del_geometry
   ~cfdm.DomainAncillary.get_geometry
   ~cfdm.DomainAncillary.has_geometry
   ~cfdm.DomainAncillary.set_geometry
   ~cfdm.DomainAncillary.del_interior_ring
   ~cfdm.DomainAncillary.get_interior_ring
   ~cfdm.DomainAncillary.has_interior_ring
   ~cfdm.DomainAncillary.set_interior_ring
   ~cfdm.DomainAncillary.del_node_count
   ~cfdm.DomainAncillary.get_node_count
   ~cfdm.DomainAncillary.has_node_count
   ~cfdm.DomainAncillary.set_node_count
   ~cfdm.DomainAncillary.del_part_node_count
   ~cfdm.DomainAncillary.get_part_node_count
   ~cfdm.DomainAncillary.has_part_node_count
   ~cfdm.DomainAncillary.set_part_node_count
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DomainAncillary.interior_ring

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.copy
   ~cfdm.DomainAncillary.creation_commands
   ~cfdm.DomainAncillary.equals
   ~cfdm.DomainAncillary.uncompress
   ~cfdm.DomainAncillary.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.nc_del_variable
   ~cfdm.DomainAncillary.nc_get_variable
   ~cfdm.DomainAncillary.nc_has_variable
   ~cfdm.DomainAncillary.nc_set_variable

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.nc_variable_groups
   ~cfdm.DomainAncillary.nc_clear_variable_groups
   ~cfdm.DomainAncillary.nc_set_variable_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAncillary.__deepcopy__
   ~cfdm.DomainAncillary.__getitem__
   ~cfdm.DomainAncillary.__repr__
   ~cfdm.DomainAncillary.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.DomainAncillary._docstring_special_substitutions
   ~cfdm.DomainAncillary._docstring_substitutions        
   ~cfdm.DomainAncillary._docstring_package_depth        
   ~cfdm.DomainAncillary._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.NumpyArray
====================

----

.. autoclass:: cfdm.core.NumpyArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst
   
   ~cfdm.core.NumpyArray.array
   ~cfdm.core.NumpyArray.dtype
   ~cfdm.core.NumpyArray.ndim
   ~cfdm.core.NumpyArray.shape
   ~cfdm.core.NumpyArray.size
      
Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.core.NumpyArray.copy

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.NumpyArray._docstring_special_substitutions
   ~cfdm.core.NumpyArray._docstring_substitutions        
   ~cfdm.core.NumpyArray._docstring_package_depth        
   ~cfdm.core.NumpyArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.PropertiesData
========================

----

.. autoclass:: cfdm.core.PropertiesData
   :no-members:
   :no-inherited-members:

   
Properties
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesData.del_property
   ~cfdm.core.PropertiesData.get_property
   ~cfdm.core.PropertiesData.has_property
   ~cfdm.core.PropertiesData.set_property
   ~cfdm.core.PropertiesData.properties
   ~cfdm.core.PropertiesData.clear_properties
   ~cfdm.core.PropertiesData.set_properties

Data
----

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesData.del_data
   ~cfdm.core.PropertiesData.get_data
   ~cfdm.core.PropertiesData.has_data
   ~cfdm.core.PropertiesData.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.PropertiesData.data

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesData.copy
   ~cfdm.core.PropertiesData.has_bounds

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesData.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.PropertiesData._docstring_special_substitutions
   ~cfdm.core.PropertiesData._docstring_substitutions        
   ~cfdm.core.PropertiesData._docstring_package_depth        
   ~cfdm.core.PropertiesData._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.CoordinateConversion
=========================

----

.. autoclass:: cfdm.CoordinateConversion
   :no-members:
   :no-inherited-members:
   
Parameter terms
---------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateConversion.del_parameter
   ~cfdm.CoordinateConversion.get_parameter
   ~cfdm.CoordinateConversion.has_parameter
   ~cfdm.CoordinateConversion.set_parameter
   ~cfdm.CoordinateConversion.parameters
   ~cfdm.CoordinateConversion.clear_parameters
   ~cfdm.CoordinateConversion.set_parameters
    
Domain ancillary terms
----------------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateConversion.del_domain_ancillary
   ~cfdm.CoordinateConversion.get_domain_ancillary
   ~cfdm.CoordinateConversion.has_domain_ancillary
   ~cfdm.CoordinateConversion.set_domain_ancillary
   ~cfdm.CoordinateConversion.domain_ancillaries
   ~cfdm.CoordinateConversion.clear_domain_ancillaries
   ~cfdm.CoordinateConversion.set_domain_ancillaries

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateConversion.copy
   ~cfdm.CoordinateConversion.equals

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateConversion.__bool__
   ~cfdm.CoordinateConversion.__deepcopy__
   ~cfdm.CoordinateConversion.__nonzero__
   ~cfdm.CoordinateConversion.__repr__
   ~cfdm.CoordinateConversion.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.CoordinateConversion._docstring_special_substitutions
   ~cfdm.CoordinateConversion._docstring_substitutions        
   ~cfdm.CoordinateConversion._docstring_package_depth        
   ~cfdm.CoordinateConversion._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.CoordinateConversion
==============================

----

.. autoclass:: cfdm.core.CoordinateConversion
   :no-members:
   :no-inherited-members:

Parameter terms
---------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateConversion.del_parameter
   ~cfdm.core.CoordinateConversion.get_parameter
   ~cfdm.core.CoordinateConversion.has_parameter
   ~cfdm.core.CoordinateConversion.set_parameter
   ~cfdm.core.CoordinateConversion.parameters
   ~cfdm.core.CoordinateConversion.clear_parameters
   ~cfdm.core.CoordinateConversion.set_parameters
   
Domain ancillary terms
----------------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateConversion.del_domain_ancillary
   ~cfdm.core.CoordinateConversion.get_domain_ancillary
   ~cfdm.core.CoordinateConversion.has_domain_ancillary
   ~cfdm.core.CoordinateConversion.set_domain_ancillary
   ~cfdm.core.CoordinateConversion.domain_ancillaries
   ~cfdm.core.CoordinateConversion.clear_domain_ancillaries
   ~cfdm.core.CoordinateConversion.set_domain_ancillaries

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateConversion.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateConversion.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.CoordinateConversion._docstring_special_substitutions
   ~cfdm.core.CoordinateConversion._docstring_substitutions        
   ~cfdm.core.CoordinateConversion._docstring_package_depth        
   ~cfdm.core.CoordinateConversion._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.DomainAxis
====================

----

.. autoclass:: cfdm.core.DomainAxis
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAxis.construct_type

Size
----

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.core.DomainAxis.del_size
   ~cfdm.core.DomainAxis.get_size
   ~cfdm.core.DomainAxis.has_size
   ~cfdm.core.DomainAxis.set_size

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAxis.construct_type

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAxis.copy
   
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAxis.__deepcopy__
   
Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.DomainAxis._docstring_special_substitutions
   ~cfdm.core.DomainAxis._docstring_substitutions        
   ~cfdm.core.DomainAxis._docstring_package_depth        
   ~cfdm.core.DomainAxis._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Constructs
====================

----

.. autoclass:: cfdm.core.Constructs
   :no-members:
   :no-inherited-members:

Filtering
---------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.filter_by_type

Constructs and identifiers
--------------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.key
   ~cfdm.core.Constructs.value
   ~cfdm.core.Constructs.construct_type
   ~cfdm.core.Constructs.construct_types
   ~cfdm.core.Constructs.new_identifier
   ~cfdm.core.Constructs.replace

Data axes
---------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.data_axes
   ~cfdm.core.Constructs.get_data_axes

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.copy
   ~cfdm.core.Constructs.shallow_copy
   ~cfdm.core.Constructs.ordered
   ~cfdm.core.Constructs.todict

Dictionary-access methods
-------------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.get
   ~cfdm.core.Constructs.items
   ~cfdm.core.Constructs.keys
   ~cfdm.core.Constructs.values
   ~cfdm.core.Constructs.__getitem__

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Constructs.__contains__
   ~cfdm.core.Constructs.__copy__
   ~cfdm.core.Constructs.__deepcopy__
   ~cfdm.core.Constructs.__getitem__
   ~cfdm.core.Constructs.__iter__
   ~cfdm.core.Constructs.__len__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Constructs._docstring_special_substitutions
   ~cfdm.core.Constructs._docstring_substitutions        
   ~cfdm.core.Constructs._docstring_package_depth        
   ~cfdm.core.Constructs._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.PropertiesDataBounds
==============================

----

.. autoclass:: cfdm.core.PropertiesDataBounds
   :no-members:
   :no-inherited-members:

   
Properties
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.core.PropertiesDataBounds.del_property
   ~cfdm.core.PropertiesDataBounds.get_property
   ~cfdm.core.PropertiesDataBounds.has_property
   ~cfdm.core.PropertiesDataBounds.set_property
   ~cfdm.core.PropertiesDataBounds.properties
   ~cfdm.core.PropertiesDataBounds.clear_properties
   ~cfdm.core.PropertiesDataBounds.set_properties

Data
----

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.del_data
   ~cfdm.core.PropertiesDataBounds.get_data
   ~cfdm.core.PropertiesDataBounds.has_data
   ~cfdm.core.PropertiesDataBounds.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.PropertiesDataBounds.data

Bounds
------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.del_bounds
   ~cfdm.core.PropertiesDataBounds.get_bounds
   ~cfdm.core.PropertiesDataBounds.has_bounds
   ~cfdm.core.PropertiesDataBounds.set_bounds
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.PropertiesDataBounds.bounds
   Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.del_geometry
   ~cfdm.core.PropertiesDataBounds.get_geometry
   ~cfdm.core.PropertiesDataBounds.has_geometry
   ~cfdm.core.PropertiesDataBounds.set_geometry
   ~cfdm.core.PropertiesDataBounds.del_interior_ring
   ~cfdm.core.PropertiesDataBounds.get_interior_ring
   ~cfdm.core.PropertiesDataBounds.has_interior_ring
   ~cfdm.core.PropertiesDataBounds.set_interior_ring

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.PropertiesDataBounds.interior_ring

Modification
------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.del_data
   ~cfdm.core.PropertiesDataBounds.set_data
   ~cfdm.core.PropertiesDataBounds.del_bounds
   ~cfdm.core.PropertiesDataBounds.set_bounds
      
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.PropertiesDataBounds.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.PropertiesDataBounds._docstring_special_substitutions
   ~cfdm.core.PropertiesDataBounds._docstring_substitutions        
   ~cfdm.core.PropertiesDataBounds._docstring_package_depth        
   ~cfdm.core.PropertiesDataBounds._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.ParametersDomainAncillaries
=====================================

----

.. autoclass:: cfdm.core.ParametersDomainAncillaries
   :no-members:
   :no-inherited-members:

Parameter terms
---------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.ParametersDomainAncillaries.del_parameter
   ~cfdm.core.ParametersDomainAncillaries.get_parameter
   ~cfdm.core.ParametersDomainAncillaries.has_parameter
   ~cfdm.core.ParametersDomainAncillaries.set_parameter
   ~cfdm.core.ParametersDomainAncillaries.parameters
   ~cfdm.core.ParametersDomainAncillaries.clear_parameters
   ~cfdm.core.ParametersDomainAncillaries.set_parameters

Domain ancillary terms
----------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.ParametersDomainAncillaries.del_domain_ancillary
   ~cfdm.core.ParametersDomainAncillaries.get_domain_ancillary
   ~cfdm.core.ParametersDomainAncillaries.has_domain_ancillary
   ~cfdm.core.ParametersDomainAncillaries.set_domain_ancillary
   ~cfdm.core.ParametersDomainAncillaries.domain_ancillaries
   ~cfdm.core.ParametersDomainAncillaries.clear_domain_ancillaries
   ~cfdm.core.ParametersDomainAncillaries.set_domain_ancillaries
   
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.ParametersDomainAncillaries.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.ParametersDomainAncillaries.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.ParametersDomainAncilliaries._docstring_special_substitutions
   ~cfdm.core.ParametersDomainAncilliaries._docstring_substitutions        
   ~cfdm.core.ParametersDomainAncilliaries._docstring_package_depth        
   ~cfdm.core.ParametersDomainAncilliaries._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-Field:

cfdm.Field
==========

----

.. autoclass:: cfdm.Field
   :no-members:
   :no-inherited-members:

.. _Field-Inspection:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.dump
   ~cfdm.Field.identity  
   ~cfdm.Field.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Field.construct_type

.. _Field-Properties:

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.del_property
   ~cfdm.Field.get_property
   ~cfdm.Field.has_property
   ~cfdm.Field.set_property
   ~cfdm.Field.properties
   ~cfdm.Field.clear_properties
   ~cfdm.Field.set_properties

.. _Field-Data:

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.apply_masking
   ~cfdm.Field.del_data
   ~cfdm.Field.get_data
   ~cfdm.Field.has_data
   ~cfdm.Field.set_data
   ~cfdm.Field.del_data_axes
   ~cfdm.Field.get_data_axes
   ~cfdm.Field.has_data_axes
   ~cfdm.Field.set_data_axes
   ~cfdm.Field.insert_dimension
   ~cfdm.Field.squeeze
   ~cfdm.Field.transpose
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Field.data
   ~cfdm.Field.dtype
   ~cfdm.Field.ndim
   ~cfdm.Field.shape
   ~cfdm.Field.size
   
.. _Field-Metadata-constructs:   
   
Metadata constructs
-------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.construct
   ~cfdm.Field.construct_key
   ~cfdm.Field.construct_item
   ~cfdm.Field.del_construct
   ~cfdm.Field.get_construct
   ~cfdm.Field.has_construct
   ~cfdm.Field.set_construct
   ~cfdm.Field.del_data_axes
   ~cfdm.Field.get_data_axes
   ~cfdm.Field.has_data_axes
   ~cfdm.Field.set_data_axes
   ~cfdm.Field.domain_axis_key
   ~cfdm.Field.auxiliary_coordinates
   ~cfdm.Field.cell_measures
   ~cfdm.Field.cell_methods
   ~cfdm.Field.coordinates
   ~cfdm.Field.coordinate_references
   ~cfdm.Field.dimension_coordinates
   ~cfdm.Field.domain_ancillaries
   ~cfdm.Field.domain_axes
   ~cfdm.Field.field_ancillaries
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Field.constructs

.. _Field-Domain:

Domain
------


.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.get_domain
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Field.domain

.. _Field-Miscellaneous:

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.climatological_time_axes
   ~cfdm.Field.compress
   ~cfdm.Field.copy
   ~cfdm.Field.creation_commands
   ~cfdm.Field.equals
   ~cfdm.Field.convert
   ~cfdm.Field.has_bounds
   ~cfdm.Field.has_geometry
   ~cfdm.Field.uncompress
   ~cfdm.Field.get_filenames

.. _Field-NetCDF:
   
NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.nc_del_variable
   ~cfdm.Field.nc_get_variable
   ~cfdm.Field.nc_has_variable
   ~cfdm.Field.nc_set_variable 
   ~cfdm.Field.nc_global_attributes
   ~cfdm.Field.nc_clear_global_attributes
   ~cfdm.Field.nc_set_global_attribute
   ~cfdm.Field.nc_set_global_attributes

Groups
^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Field.nc_variable_groups
   ~cfdm.Field.nc_set_variable_groups
   ~cfdm.Field.nc_clear_variable_groups
   ~cfdm.Field.nc_group_attributes
   ~cfdm.Field.nc_clear_group_attributes
   ~cfdm.Field.nc_set_group_attribute
   ~cfdm.Field.nc_set_group_attributes
  
Geometries
^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Field.nc_del_geometry_variable
   ~cfdm.Field.nc_get_geometry_variable
   ~cfdm.Field.nc_has_geometry_variable
   ~cfdm.Field.nc_set_geometry_variable 
   ~cfdm.Field.nc_geometry_variable_groups
   ~cfdm.Field.nc_set_geometry_variable_groups
   ~cfdm.Field.nc_clear_geometry_variable_groups

Components
^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Field.nc_del_component_variable
   ~cfdm.Field.nc_set_component_variable
   ~cfdm.Field.nc_set_component_variable_groups
   ~cfdm.Field.nc_clear_component_variable_groups      
   ~cfdm.Field.nc_del_component_dimension
   ~cfdm.Field.nc_set_component_dimension
   ~cfdm.Field.nc_set_component_dimension_groups
   ~cfdm.Field.nc_clear_component_dimension_groups
   ~cfdm.Field.nc_del_component_sample_dimension
   ~cfdm.Field.nc_set_component_sample_dimension   
   ~cfdm.Field.nc_set_component_sample_dimension_groups
   ~cfdm.Field.nc_clear_component_sample_dimension_groups

Dataset compliance
^^^^^^^^^^^^^^^^^^

.. rubric:: Methods


.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.dataset_compliance

.. _Field-Special:

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Field.__deepcopy__
   ~cfdm.Field.__getitem__
   ~cfdm.Field.__repr__
   ~cfdm.Field.__str__

Docstring substitutions
-----------------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.Field._docstring_special_substitutions
   ~cfdm.Field._docstring_substitutions
   ~cfdm.Field._docstring_package_depth
   ~cfdm.Field._docstring_method_exclusions
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.DimensionCoordinate
=============================

----

.. autoclass:: cfdm.core.DimensionCoordinate
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DimensionCoordinate.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.del_property
   ~cfdm.core.DimensionCoordinate.get_property
   ~cfdm.core.DimensionCoordinate.has_property
   ~cfdm.core.DimensionCoordinate.set_property
   ~cfdm.core.DimensionCoordinate.properties
   ~cfdm.core.DimensionCoordinate.clear_properties
   ~cfdm.core.DimensionCoordinate.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.del_data
   ~cfdm.core.DimensionCoordinate.get_data
   ~cfdm.core.DimensionCoordinate.has_data
   ~cfdm.core.DimensionCoordinate.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DimensionCoordinate.data

Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.del_bounds
   ~cfdm.core.DimensionCoordinate.get_bounds
   ~cfdm.core.DimensionCoordinate.has_bounds
   ~cfdm.core.DimensionCoordinate.set_bounds
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DimensionCoordinate.bounds

Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.del_geometry
   ~cfdm.core.DimensionCoordinate.get_geometry
   ~cfdm.core.DimensionCoordinate.has_geometry
   ~cfdm.core.DimensionCoordinate.set_geometry
   ~cfdm.core.DimensionCoordinate.del_interior_ring
   ~cfdm.core.DimensionCoordinate.get_interior_ring
   ~cfdm.core.DimensionCoordinate.has_interior_ring
   ~cfdm.core.DimensionCoordinate.set_interior_ring

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DimensionCoordinate.interior_ring

Climatology
^^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.del_climatology
   ~cfdm.core.DimensionCoordinate.get_climatology
   ~cfdm.core.DimensionCoordinate.is_climatology
   ~cfdm.core.DimensionCoordinate.set_climatology

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DimensionCoordinate.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.DimensionCoordinate._docstring_special_substitutions
   ~cfdm.core.DimensionCoordinate._docstring_substitutions        
   ~cfdm.core.DimensionCoordinate._docstring_package_depth        
   ~cfdm.core.DimensionCoordinate._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.NumpyArray
===============

----

.. autoclass:: cfdm.NumpyArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NumpyArray.get_compression_type
   ~cfdm.NumpyArray.get_subspace
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst
   
   ~cfdm.NumpyArray.array
   ~cfdm.NumpyArray.dtype
   ~cfdm.NumpyArray.ndim
   ~cfdm.NumpyArray.shape
   ~cfdm.NumpyArray.size

Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NumpyArray.copy
   ~cfdm.NumpyArray.to_memory
   
Special
-------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NumpyArray.__getitem__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.NumpyArray._docstring_special_substitutions
   ~cfdm.NumpyArray._docstring_substitutions        
   ~cfdm.NumpyArray._docstring_package_depth        
   ~cfdm.NumpyArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.RaggedIndexedContiguousArray
=================================

----

.. autoclass:: cfdm.RaggedIndexedContiguousArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.RaggedIndexedContiguousArray.get_compressed_axes
   ~cfdm.RaggedIndexedContiguousArray.get_compressed_dimension
   ~cfdm.RaggedIndexedContiguousArray.get_compression_type
   ~cfdm.RaggedIndexedContiguousArray.get_count
   ~cfdm.RaggedIndexedContiguousArray.get_index
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.RaggedIndexedContiguousArray.array
   ~cfdm.RaggedIndexedContiguousArray.compressed_array
   ~cfdm.RaggedIndexedContiguousArray.dtype
   ~cfdm.RaggedIndexedContiguousArray.ndim
   ~cfdm.RaggedIndexedContiguousArray.shape
   ~cfdm.RaggedIndexedContiguousArray.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedIndexedContiguousArray.copy
   ~cfdm.RaggedIndexedContiguousArray.get_subspace
   ~cfdm.RaggedIndexedContiguousArray.source
   ~cfdm.RaggedIndexedContiguousArray.to_memory

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedIndexedContiguousArray.__array__
   ~cfdm.RaggedIndexedContiguousArray.__deepcopy__
   ~cfdm.RaggedIndexedContiguousArray.__getitem__
   ~cfdm.RaggedIndexedContiguousArray.__repr__
   ~cfdm.RaggedIndexedContiguousArray.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.RaggedIndexedContiguousArray._docstring_special_substitutions
   ~cfdm.RaggedIndexedContiguousArray._docstring_substitutions        
   ~cfdm.RaggedIndexedContiguousArray._docstring_package_depth        
   ~cfdm.RaggedIndexedContiguousArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Count
==========

----

.. autoclass:: cfdm.Count
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.dump
   ~cfdm.Count.identity  
   ~cfdm.Count.identities
  
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.del_property
   ~cfdm.Count.get_property
   ~cfdm.Count.has_property
   ~cfdm.Count.set_property
   ~cfdm.Count.properties
   ~cfdm.Count.clear_properties
   ~cfdm.Count.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.apply_masking
   ~cfdm.Count.del_data
   ~cfdm.Count.get_data
   ~cfdm.Count.has_data
   ~cfdm.Count.set_data  
   ~cfdm.Count.insert_dimension
   ~cfdm.Count.squeeze
   ~cfdm.Count.transpose
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Count.data
   ~cfdm.Count.dtype
   ~cfdm.Count.ndim
   ~cfdm.Count.shape
   ~cfdm.Count.size 

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.copy
   ~cfdm.Count.creation_commands
   ~cfdm.Count.equals
   ~cfdm.Count.get_filenames
   ~cfdm.Count.has_bounds
   ~cfdm.Count.uncompress

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.nc_del_variable
   ~cfdm.Count.nc_get_variable
   ~cfdm.Count.nc_has_variable
   ~cfdm.Count.nc_set_variable
   ~cfdm.Count.nc_del_dimension
   ~cfdm.Count.nc_get_dimension
   ~cfdm.Count.nc_has_dimension
   ~cfdm.Count.nc_set_dimension

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.nc_variable_groups
   ~cfdm.Count.nc_clear_variable_groups
   ~cfdm.Count.nc_set_variable_groups
   ~cfdm.Count.nc_dimension_groups
   ~cfdm.Count.nc_clear_dimension_groups
   ~cfdm.Count.nc_set_dimension_groups
   
Sample dimension
^^^^^^^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.nc_del_sample_dimension
   ~cfdm.Count.nc_get_sample_dimension
   ~cfdm.Count.nc_has_sample_dimension
   ~cfdm.Count.nc_set_sample_dimension
   ~cfdm.Count.nc_sample_dimension_groups
   ~cfdm.Count.nc_clear_sample_dimension_groups
   ~cfdm.Count.nc_set_sample_dimension_groups
   
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Count.__deepcopy__
   ~cfdm.Count.__getitem__
   ~cfdm.Count.__repr__
   ~cfdm.Count.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Count._docstring_special_substitutions
   ~cfdm.Count._docstring_substitutions        
   ~cfdm.Count._docstring_package_depth        
   ~cfdm.Count._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-DimensionCoordinate:

cfdm.DimensionCoordinate
========================

----

.. autoclass:: cfdm.DimensionCoordinate
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.dump
   ~cfdm.DimensionCoordinate.identity  
   ~cfdm.DimensionCoordinate.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DimensionCoordinate.construct_type
   
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.del_property
   ~cfdm.DimensionCoordinate.get_property
   ~cfdm.DimensionCoordinate.has_property
   ~cfdm.DimensionCoordinate.set_property
   ~cfdm.DimensionCoordinate.properties
   ~cfdm.DimensionCoordinate.clear_properties
   ~cfdm.DimensionCoordinate.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.apply_masking
   ~cfdm.DimensionCoordinate.del_data
   ~cfdm.DimensionCoordinate.get_data
   ~cfdm.DimensionCoordinate.has_data
   ~cfdm.DimensionCoordinate.set_data
   ~cfdm.DimensionCoordinate.insert_dimension
   ~cfdm.DimensionCoordinate.squeeze
   ~cfdm.DimensionCoordinate.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DimensionCoordinate.data
   ~cfdm.DimensionCoordinate.dtype
   ~cfdm.DimensionCoordinate.ndim
   ~cfdm.DimensionCoordinate.shape
   ~cfdm.DimensionCoordinate.size
   
Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.del_bounds
   ~cfdm.DimensionCoordinate.get_bounds
   ~cfdm.DimensionCoordinate.has_bounds
   ~cfdm.DimensionCoordinate.set_bounds
   ~cfdm.DimensionCoordinate.get_bounds_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DimensionCoordinate.bounds

Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.del_geometry
   ~cfdm.DimensionCoordinate.get_geometry
   ~cfdm.DimensionCoordinate.has_geometry
   ~cfdm.DimensionCoordinate.set_geometry
   ~cfdm.DimensionCoordinate.del_interior_ring
   ~cfdm.DimensionCoordinate.get_interior_ring
   ~cfdm.DimensionCoordinate.has_interior_ring
   ~cfdm.DimensionCoordinate.set_interior_ring
   ~cfdm.DimensionCoordinate.del_node_count
   ~cfdm.DimensionCoordinate.get_node_count
   ~cfdm.DimensionCoordinate.has_node_count
   ~cfdm.DimensionCoordinate.set_node_count
   ~cfdm.DimensionCoordinate.del_part_node_count
   ~cfdm.DimensionCoordinate.get_part_node_count
   ~cfdm.DimensionCoordinate.has_part_node_count
   ~cfdm.DimensionCoordinate.set_part_node_count

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DimensionCoordinate.interior_ring

Climatology
^^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.del_climatology
   ~cfdm.DimensionCoordinate.get_climatology
   ~cfdm.DimensionCoordinate.is_climatology
   ~cfdm.DimensionCoordinate.set_climatology

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.copy
   ~cfdm.DimensionCoordinate.creation_commands
   ~cfdm.DimensionCoordinate.equals
   ~cfdm.DimensionCoordinate.uncompress
   ~cfdm.DimensionCoordinate.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.nc_del_variable
   ~cfdm.DimensionCoordinate.nc_get_variable
   ~cfdm.DimensionCoordinate.nc_has_variable
   ~cfdm.DimensionCoordinate.nc_set_variable

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.nc_variable_groups
   ~cfdm.DimensionCoordinate.nc_clear_variable_groups
   ~cfdm.DimensionCoordinate.nc_set_variable_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DimensionCoordinate.__deepcopy__
   ~cfdm.DimensionCoordinate.__getitem__
   ~cfdm.DimensionCoordinate.__repr__
   ~cfdm.DimensionCoordinate.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.DimensionCoordinate._docstring_special_substitutions
   ~cfdm.DimensionCoordinate._docstring_substitutions        
   ~cfdm.DimensionCoordinate._docstring_package_depth        
   ~cfdm.DimensionCoordinate._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Data
=========

.. autoclass:: cfdm.Data
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
	    
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Data.array
   ~cfdm.Data.dtype
   ~cfdm.Data.ndim
   ~cfdm.Data.shape
   ~cfdm.Data.size
   
Units
-----

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.del_units
   ~cfdm.Data.get_units
   ~cfdm.Data.has_units
   ~cfdm.Data.set_units
   ~cfdm.Data.del_calendar
   ~cfdm.Data.get_calendar
   ~cfdm.Data.has_calendar
   ~cfdm.Data.set_calendar

Data creation routines
----------------------

Ones and zeros
^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.empty

From existing data
^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.copy

Data manipulation routines
--------------------------

Changing data shape
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.flatten


Transpose-like operations
^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.transpose

Changing number of dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.insert_dimension
   ~cfdm.Data.squeeze

Adding and removing elements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.unique

Date-time support
-----------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.del_calendar
   ~cfdm.Data.get_calendar
   ~cfdm.Data.has_calendar
   ~cfdm.Data.set_calendar

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Data.datetime_array
   ~cfdm.Data.datetime_as_string
 
Indexing routines
-----------------

Single value selection
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.first_element
   ~cfdm.Data.second_element
   ~cfdm.Data.last_element

Logic functions
---------------

Truth value testing
^^^^^^^^^^^^^^^^^^^
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.any

Comparison
^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.equals

Mask support
------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.apply_masking
   ~cfdm.Data.filled
   ~cfdm.Data.del_fill_value
   ~cfdm.Data.get_fill_value
   ~cfdm.Data.has_fill_value
   ~cfdm.Data.set_fill_value
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Data.mask

Mathematical functions
----------------------

Sums, products, differences
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.sum

Set routines
-------------

Making proper sets
^^^^^^^^^^^^^^^^^^    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.unique
	    
Sorting, searching, and counting
--------------------------------

Statistics
----------

Order statistics
^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.maximum
   ~cfdm.Data.minimum
   ~cfdm.Data.max
   ~cfdm.Data.min

Sums
^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.sum

Compression by convention
-------------------------
   
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.get_compressed_axes
   ~cfdm.Data.get_compressed_dimension
   ~cfdm.Data.get_compression_type
   ~cfdm.Data.get_count
   ~cfdm.Data.get_index
   ~cfdm.Data.get_list
   ~cfdm.Data.uncompress

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Data.compressed_array

Miscellaneous
-------------
   
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.creation_commands
   ~cfdm.Data.get_filenames
   ~cfdm.Data.source

Performance
-----------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Data.nc_clear_hdf5_chunksizes
   ~cfdm.Data.nc_hdf5_chunksizes
   ~cfdm.Data.nc_set_hdf5_chunksizes
   ~cfdm.Data.to_memory
 
Special
-------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Data.__array__
   ~cfdm.Data.__deepcopy__
   ~cfdm.Data.__getitem__ 
   ~cfdm.Data.__int__
   ~cfdm.Data.__iter__ 
   ~cfdm.Data.__repr__
   ~cfdm.Data.__setitem__ 
   ~cfdm.Data.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Data._docstring_special_substitutions
   ~cfdm.Data._docstring_substitutions        
   ~cfdm.Data._docstring_package_depth        
   ~cfdm.Data._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.CoordinateReference
=============================

----

.. autoclass:: cfdm.core.CoordinateReference
   :no-members:
   :no-inherited-members:
   
Inspection
----------

.. rubric:: Attributes
 
.. autosummary::
   :toctree: ../attribute/
   :template: attribute.rst
	      
   ~cfdm.core.CoordinateReference.construct_type

Coordinates
-----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateReference.del_coordinate
   ~cfdm.core.CoordinateReference.has_coordinate
   ~cfdm.core.CoordinateReference.set_coordinate
   ~cfdm.core.CoordinateReference.coordinates
   ~cfdm.core.CoordinateReference.clear_coordinates
   ~cfdm.core.CoordinateReference.set_coordinates

Datum
-----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateReference.del_datum
   ~cfdm.core.CoordinateReference.get_datum
   ~cfdm.core.CoordinateReference.set_datum

.. rubric:: Attributes
	    
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.CoordinateReference.datum

Coordinate conversion
---------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateReference.del_coordinate_conversion
   ~cfdm.core.CoordinateReference.get_coordinate_conversion
   ~cfdm.core.CoordinateReference.set_coordinate_conversion

.. rubric:: Attributes
	    
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.CoordinateReference.coordinate_conversion

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateReference.copy
      
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CoordinateReference.__deepcopy__

Docstring substitutions
-----------------------                           
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.CoordinateReference._docstring_special_substitutions
   ~cfdm.core.CoordinateReference._docstring_substitutions        
   ~cfdm.core.CoordinateReference._docstring_package_depth        
   ~cfdm.core.CoordinateReference._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Bounds
===========

----

.. autoclass:: cfdm.Bounds
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.dump
   ~cfdm.Bounds.identity
   ~cfdm.Bounds.identities
   
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.del_property
   ~cfdm.Bounds.get_property
   ~cfdm.Bounds.has_property
   ~cfdm.Bounds.set_property
   ~cfdm.Bounds.properties
   ~cfdm.Bounds.clear_properties
   ~cfdm.Bounds.set_properties
   ~cfdm.Bounds.inherited_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.apply_masking
   ~cfdm.Bounds.del_data
   ~cfdm.Bounds.get_data
   ~cfdm.Bounds.has_data
   ~cfdm.Bounds.set_data
   ~cfdm.Bounds.insert_dimension
   ~cfdm.Bounds.squeeze
   ~cfdm.Bounds.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Bounds.data
   ~cfdm.Bounds.dtype
   ~cfdm.Bounds.ndim
   ~cfdm.Bounds.shape
   ~cfdm.Bounds.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.copy
   ~cfdm.Bounds.creation_commands
   ~cfdm.Bounds.equals
   ~cfdm.Bounds.has_bounds
   ~cfdm.Bounds.uncompress
   ~cfdm.Bounds.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.nc_del_variable
   ~cfdm.Bounds.nc_get_variable
   ~cfdm.Bounds.nc_has_variable
   ~cfdm.Bounds.nc_set_variable 
   ~cfdm.Bounds.nc_del_dimension
   ~cfdm.Bounds.nc_get_dimension
   ~cfdm.Bounds.nc_has_dimension
   ~cfdm.Bounds.nc_set_dimension

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.nc_variable_groups
   ~cfdm.Bounds.nc_clear_variable_groups
   ~cfdm.Bounds.nc_set_variable_groups
   ~cfdm.Bounds.nc_dimension_groups
   ~cfdm.Bounds.nc_clear_dimension_groups
   ~cfdm.Bounds.nc_set_dimension_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Bounds.__deepcopy__
   ~cfdm.Bounds.__getitem__
   ~cfdm.Bounds.__repr__
   ~cfdm.Bounds.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Bounds._docstring_special_substitutions
   ~cfdm.Bounds._docstring_substitutions        
   ~cfdm.Bounds._docstring_package_depth        
   ~cfdm.Bounds._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.NetCDFArray
================

----

.. autoclass:: cfdm.NetCDFArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   

   ~cfdm.NetCDFArray.get_ncvar
   ~cfdm.NetCDFArray.get_varid
   ~cfdm.NetCDFArray.get_compression_type
   ~cfdm.NetCDFArray.get_subspace
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst
   
   ~cfdm.NetCDFArray.array
   ~cfdm.NetCDFArray.dtype
   ~cfdm.NetCDFArray.ndim
   ~cfdm.NetCDFArray.shape
   ~cfdm.NetCDFArray.size

File
----
   
.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NetCDFArray.close
   ~cfdm.NetCDFArray.open
   ~cfdm.NetCDFArray.get_filename
   ~cfdm.NetCDFArray.get_group
   ~cfdm.NetCDFArray.get_mask
   
Miscellaneous
-------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NetCDFArray.copy
   ~cfdm.NetCDFArray.to_memory
   
Special
-------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.NetCDFArray.__getitem__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.NetCDFArray._docstring_special_substitutions
   ~cfdm.NetCDFArray._docstring_substitutions        
   ~cfdm.NetCDFArray._docstring_package_depth        
   ~cfdm.NetCDFArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.DomainAncillary
=========================

----

.. autoclass:: cfdm.core.DomainAncillary
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAncillary.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.del_property
   ~cfdm.core.DomainAncillary.get_property
   ~cfdm.core.DomainAncillary.has_property
   ~cfdm.core.DomainAncillary.set_property
   ~cfdm.core.DomainAncillary.properties
   ~cfdm.core.DomainAncillary.clear_properties
   ~cfdm.core.DomainAncillary.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.del_data
   ~cfdm.core.DomainAncillary.get_data
   ~cfdm.core.DomainAncillary.has_data
   ~cfdm.core.DomainAncillary.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAncillary.data

Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.del_bounds
   ~cfdm.core.DomainAncillary.get_bounds
   ~cfdm.core.DomainAncillary.has_bounds
   ~cfdm.core.DomainAncillary.set_bounds
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAncillary.bounds

Geometries
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.del_geometry
   ~cfdm.core.DomainAncillary.get_geometry
   ~cfdm.core.DomainAncillary.has_geometry
   ~cfdm.core.DomainAncillary.set_geometry
   ~cfdm.core.DomainAncillary.del_interior_ring
   ~cfdm.core.DomainAncillary.get_interior_ring
   ~cfdm.core.DomainAncillary.has_interior_ring
   ~cfdm.core.DomainAncillary.set_interior_ring
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.DomainAncillary.interior_ring

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.DomainAncillary.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.DomainAncillary._docstring_special_substitutions
   ~cfdm.core.DomainAncillary._docstring_substitutions        
   ~cfdm.core.DomainAncillary._docstring_package_depth        
   ~cfdm.core.DomainAncillary._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.RaggedContiguousArray
==========================

----

.. autoclass:: cfdm.RaggedContiguousArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
    ~cfdm.RaggedContiguousArray.get_compressed_axes
   ~cfdm.RaggedContiguousArray.get_compressed_dimension
   ~cfdm.RaggedContiguousArray.get_compression_type
   ~cfdm.RaggedContiguousArray.get_count
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.RaggedContiguousArray.array
   ~cfdm.RaggedContiguousArray.compressed_array
   ~cfdm.RaggedContiguousArray.dtype
   ~cfdm.RaggedContiguousArray.ndim
   ~cfdm.RaggedContiguousArray.shape
   ~cfdm.RaggedContiguousArray.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedContiguousArray.copy
   ~cfdm.RaggedContiguousArray.get_subspace
   ~cfdm.RaggedContiguousArray.source
   ~cfdm.RaggedContiguousArray.to_memory

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.RaggedContiguousArray.__array__
   ~cfdm.RaggedContiguousArray.__deepcopy__
   ~cfdm.RaggedContiguousArray.__getitem__
   ~cfdm.RaggedContiguousArray.__repr__
   ~cfdm.RaggedContiguousArray.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.RaggedContiguousArray._docstring_special_substitutions
   ~cfdm.RaggedContiguousArray._docstring_substitutions        
   ~cfdm.RaggedContiguousArray._docstring_package_depth        
   ~cfdm.RaggedContiguousArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.InteriorRing
=================

----

.. autoclass:: cfdm.InteriorRing
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.dump
   ~cfdm.InteriorRing.identity  
   ~cfdm.InteriorRing.identities

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.del_property
   ~cfdm.InteriorRing.get_property
   ~cfdm.InteriorRing.has_property
   ~cfdm.InteriorRing.set_property
   ~cfdm.InteriorRing.properties
   ~cfdm.InteriorRing.clear_properties
   ~cfdm.InteriorRing.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.apply_masking
   ~cfdm.InteriorRing.del_data
   ~cfdm.InteriorRing.get_data
   ~cfdm.InteriorRing.has_data
   ~cfdm.InteriorRing.set_data
   ~cfdm.InteriorRing.insert_dimension
   ~cfdm.InteriorRing.squeeze
   ~cfdm.InteriorRing.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.InteriorRing.data
   ~cfdm.InteriorRing.dtype
   ~cfdm.InteriorRing.ndim
   ~cfdm.InteriorRing.shape
   ~cfdm.InteriorRing.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.copy
   ~cfdm.InteriorRing.creation_commands
   ~cfdm.InteriorRing.equals
   ~cfdm.InteriorRing.has_bounds
   ~cfdm.InteriorRing.uncompress
   ~cfdm.InteriorRing.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.nc_del_variable
   ~cfdm.InteriorRing.nc_get_variable
   ~cfdm.InteriorRing.nc_has_variable
   ~cfdm.InteriorRing.nc_set_variable
   ~cfdm.InteriorRing.nc_del_dimension
   ~cfdm.InteriorRing.nc_get_dimension
   ~cfdm.InteriorRing.nc_has_dimension
   ~cfdm.InteriorRing.nc_set_dimension
   
Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.nc_variable_groups
   ~cfdm.InteriorRing.nc_clear_variable_groups
   ~cfdm.InteriorRing.nc_set_variable_groups
   ~cfdm.InteriorRing.nc_dimension_groups
   ~cfdm.InteriorRing.nc_clear_dimension_groups
   ~cfdm.InteriorRing.nc_set_dimension_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.InteriorRing.__deepcopy__
   ~cfdm.InteriorRing.__getitem__
   ~cfdm.InteriorRing.__repr__
   ~cfdm.InteriorRing.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.InteriorRing._docstring_special_substitutions
   ~cfdm.InteriorRing._docstring_substitutions        
   ~cfdm.InteriorRing._docstring_package_depth        
   ~cfdm.InteriorRing._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-FieldAncillary:

cfdm.FieldAncillary
===================

----

.. autoclass:: cfdm.FieldAncillary
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.dump
   ~cfdm.FieldAncillary.identity  
   ~cfdm.FieldAncillary.identities
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.FieldAncillary.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.del_property
   ~cfdm.FieldAncillary.get_property
   ~cfdm.FieldAncillary.has_property
   ~cfdm.FieldAncillary.set_property
   ~cfdm.FieldAncillary.properties
   ~cfdm.FieldAncillary.clear_properties
   ~cfdm.FieldAncillary.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.apply_masking
   ~cfdm.FieldAncillary.del_data
   ~cfdm.FieldAncillary.get_data
   ~cfdm.FieldAncillary.has_data
   ~cfdm.FieldAncillary.set_data
   ~cfdm.FieldAncillary.insert_dimension
   ~cfdm.FieldAncillary.squeeze
   ~cfdm.FieldAncillary.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.FieldAncillary.data
   ~cfdm.FieldAncillary.dtype
   ~cfdm.FieldAncillary.ndim
   ~cfdm.FieldAncillary.shape
   ~cfdm.FieldAncillary.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.copy
   ~cfdm.FieldAncillary.creation_commands
   ~cfdm.FieldAncillary.equals
   ~cfdm.FieldAncillary.has_bounds
   ~cfdm.FieldAncillary.uncompress
   ~cfdm.FieldAncillary.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.nc_del_variable
   ~cfdm.FieldAncillary.nc_get_variable
   ~cfdm.FieldAncillary.nc_has_variable
   ~cfdm.FieldAncillary.nc_set_variable
   
Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.nc_variable_groups
   ~cfdm.FieldAncillary.nc_clear_variable_groups
   ~cfdm.FieldAncillary.nc_set_variable_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.FieldAncillary.__deepcopy__
   ~cfdm.FieldAncillary.__getitem__
   ~cfdm.FieldAncillary.__repr__
   ~cfdm.FieldAncillary.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.FieldAncillary._docstring_special_substitutions
   ~cfdm.FieldAncillary._docstring_substitutions        
   ~cfdm.FieldAncillary._docstring_package_depth        
   ~cfdm.FieldAncillary._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.FieldAncillary
========================

----

.. autoclass:: cfdm.core.FieldAncillary
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.FieldAncillary.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.FieldAncillary.del_property
   ~cfdm.core.FieldAncillary.get_property
   ~cfdm.core.FieldAncillary.has_property
   ~cfdm.core.FieldAncillary.set_property
   ~cfdm.core.FieldAncillary.properties
   ~cfdm.core.FieldAncillary.clear_properties
   ~cfdm.core.FieldAncillary.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.FieldAncillary.del_data
   ~cfdm.core.FieldAncillary.get_data
   ~cfdm.core.FieldAncillary.has_data
   ~cfdm.core.FieldAncillary.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.FieldAncillary.data

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.FieldAncillary.copy
   ~cfdm.core.FieldAncillary.has_bounds


Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.FieldAncillary.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.FieldAncillary._docstring_special_substitutions
   ~cfdm.core.FieldAncillary._docstring_substitutions        
   ~cfdm.core.FieldAncillary._docstring_package_depth        
   ~cfdm.core.FieldAncillary._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Index
==========

----

.. autoclass:: cfdm.Index
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.dump
   ~cfdm.Index.identity  
   ~cfdm.Index.identities
  
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.del_property
   ~cfdm.Index.get_property
   ~cfdm.Index.has_property
   ~cfdm.Index.set_property
   ~cfdm.Index.properties
   ~cfdm.Index.clear_properties
   ~cfdm.Index.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.apply_masking
   ~cfdm.Index.del_data
   ~cfdm.Index.get_data
   ~cfdm.Index.has_data
   ~cfdm.Index.set_data
   ~cfdm.Index.insert_dimension
   ~cfdm.Index.squeeze
   ~cfdm.Index.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Index.data
   ~cfdm.Index.dtype
   ~cfdm.Index.ndim 
   ~cfdm.Index.shape
   ~cfdm.Index.size
   
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.copy
   ~cfdm.Index.creation_commands
   ~cfdm.Index.equals
   ~cfdm.Index.get_filenames
   ~cfdm.Index.has_bounds
   ~cfdm.Index.uncompress

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.nc_del_variable
   ~cfdm.Index.nc_get_variable
   ~cfdm.Index.nc_has_variable
   ~cfdm.Index.nc_set_variable 
   ~cfdm.Index.nc_del_dimension
   ~cfdm.Index.nc_get_dimension
   ~cfdm.Index.nc_has_dimension
   ~cfdm.Index.nc_set_dimension

Groups
^^^^^^   

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.nc_variable_groups
   ~cfdm.Index.nc_clear_variable_groups
   ~cfdm.Index.nc_set_variable_groups
   ~cfdm.Index.nc_dimension_groups
   ~cfdm.Index.nc_clear_dimension_groups
   ~cfdm.Index.nc_set_dimension_groups

Sample dimension
^^^^^^^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.nc_del_sample_dimension
   ~cfdm.Index.nc_get_sample_dimension
   ~cfdm.Index.nc_has_sample_dimension
   ~cfdm.Index.nc_set_sample_dimension 
   ~cfdm.Index.nc_sample_dimension_groups
   ~cfdm.Index.nc_clear_sample_dimension_groups
   ~cfdm.Index.nc_set_sample_dimension_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Index.__deepcopy__
   ~cfdm.Index.__getitem__
   ~cfdm.Index.__repr__
   ~cfdm.Index.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Index._docstring_special_substitutions
   ~cfdm.Index._docstring_substitutions        
   ~cfdm.Index._docstring_package_depth        
   ~cfdm.Index._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Properties
====================

----

.. autoclass:: cfdm.core.Properties
   :no-members:
   :no-inherited-members:

Properties
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Properties.del_property
   ~cfdm.core.Properties.get_property
   ~cfdm.core.Properties.has_property
   ~cfdm.core.Properties.set_property
   ~cfdm.core.Properties.properties
   ~cfdm.core.Properties.clear_properties
   ~cfdm.core.Properties.set_properties
  
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Properties.copy
   ~cfdm.core.Properties.has_bounds
   ~cfdm.core.Properties.has_data

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Properties.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.Properties._docstring_special_substitutions
   ~cfdm.core.Properties._docstring_substitutions        
   ~cfdm.core.Properties._docstring_package_depth        
   ~cfdm.core.Properties._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Configuration
==================

----

.. autoclass:: cfdm.Configuration
   :no-members:
   :no-inherited-members:

Copying
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Configuration.copy

Dictionary functionality
------------------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Configuration.clear
   ~cfdm.Configuration.fromkeys
   ~cfdm.Configuration.get
   ~cfdm.Configuration.items
   ~cfdm.Configuration.keys
   ~cfdm.Configuration.pop
   ~cfdm.Configuration.popitem
   ~cfdm.Configuration.setdefault
   ~cfdm.Configuration.update
   ~cfdm.Configuration.values

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Configuration.__enter__
   ~cfdm.Configuration.__exit__
   ~cfdm.Configuration.__repr__

.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Datum
==========

----

.. autoclass:: cfdm.Datum
   :no-members:
   :no-inherited-members:

Parameter terms
---------------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Datum.del_parameter
   ~cfdm.Datum.get_parameter
   ~cfdm.Datum.has_parameter
   ~cfdm.Datum.set_parameter
   ~cfdm.Datum.parameters
   ~cfdm.Datum.set_parameters
   ~cfdm.Datum.clear_parameters

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Datum.copy
   ~cfdm.Datum.equals

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Datum.nc_del_variable
   ~cfdm.Datum.nc_get_variable
   ~cfdm.Datum.nc_has_variable
   ~cfdm.Datum.nc_set_variable

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Datum.nc_variable_groups
   ~cfdm.Datum.nc_set_variable_groups
   ~cfdm.Datum.nc_clear_variable_groups
      
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Datum.__bool__
   ~cfdm.Datum.__deepcopy__
   ~cfdm.Datum.__nonzero__
   ~cfdm.Datum.__repr__
   ~cfdm.Datum.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Datum._docstring_special_substitutions
   ~cfdm.Datum._docstring_substitutions        
   ~cfdm.Datum._docstring_package_depth        
   ~cfdm.Datum._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Data
==============

.. autoclass:: cfdm.core.Data
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
	    
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Data.array
   ~cfdm.core.Data.dtype
   ~cfdm.core.Data.ndim
   ~cfdm.core.Data.shape
   ~cfdm.core.Data.size
   
Units
-----

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.del_units
   ~cfdm.core.Data.get_units
   ~cfdm.core.Data.has_units
   ~cfdm.core.Data.set_units
   ~cfdm.core.Data.del_calendar
   ~cfdm.core.Data.get_calendar
   ~cfdm.core.Data.has_calendar
   ~cfdm.core.Data.set_calendar

Data creation routines
----------------------

From existing data
^^^^^^^^^^^^^^^^^^

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.copy

Date-time support
-----------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.del_calendar
   ~cfdm.core.Data.get_calendar
   ~cfdm.core.Data.has_calendar
   ~cfdm.core.Data.set_calendar
 
Mask support
------------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.del_fill_value
   ~cfdm.core.Data.get_fill_value
   ~cfdm.core.Data.has_fill_value
   ~cfdm.core.Data.set_fill_value

Miscellaneous
-------------
   
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.source

Special
-------

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Data.__deepcopy__
   ~cfdm.core.Data.__repr__
   ~cfdm.core.Data.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.core.Data._docstring_special_substitutions
   ~cfdm.core.Data._docstring_substitutions        
   ~cfdm.core.Data._docstring_package_depth        
   ~cfdm.core.Data._docstring_method_exclusions
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.CompressedArray
====================

----

.. autoclass:: cfdm.CompressedArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.CompressedArray.array
   ~cfdm.CompressedArray.compressed_array
   ~cfdm.CompressedArray.dtype
   ~cfdm.CompressedArray.ndim
   ~cfdm.CompressedArray.shape
   ~cfdm.CompressedArray.size
 
Compression
-----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.CompressedArray.get_compressed_axes
   ~cfdm.CompressedArray.get_compressed_dimension
   ~cfdm.CompressedArray.get_compression_type

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CompressedArray.copy
   ~cfdm.CompressedArray.get_subspace
   ~cfdm.CompressedArray.source
   ~cfdm.CompressedArray.to_memory

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CompressedArray.__array__
   ~cfdm.CompressedArray.__deepcopy__
   ~cfdm.CompressedArray.__getitem__
   ~cfdm.CompressedArray.__repr__
   ~cfdm.CompressedArray.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.CompressedArray._docstring_special_substitutions
   ~cfdm.CompressedArray._docstring_substitutions        
   ~cfdm.CompressedArray._docstring_package_depth        
   ~cfdm.CompressedArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.GatheredArray
==================

----

.. autoclass:: cfdm.GatheredArray
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.GatheredArray.get_compressed_axes
   ~cfdm.GatheredArray.get_compressed_dimension
   ~cfdm.GatheredArray.get_compression_type
   ~cfdm.GatheredArray.get_list
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.GatheredArray.array
   ~cfdm.GatheredArray.compressed_array
   ~cfdm.GatheredArray.dtype
   ~cfdm.GatheredArray.ndim
   ~cfdm.GatheredArray.shape
   ~cfdm.GatheredArray.size

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.GatheredArray.copy
   ~cfdm.GatheredArray.get_subspace
   ~cfdm.GatheredArray.source
   ~cfdm.GatheredArray.to_memory

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.GatheredArray.__array__
   ~cfdm.GatheredArray.__deepcopy__
   ~cfdm.GatheredArray.__getitem__
   ~cfdm.GatheredArray.__repr__
   ~cfdm.GatheredArray.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.GatherArray._docstring_special_substitutions
   ~cfdm.GatherArray._docstring_substitutions        
   ~cfdm.GatherArray._docstring_package_depth        
   ~cfdm.GatherArray._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Parameters
====================

----

.. autoclass:: cfdm.core.Parameters
   :no-members:
   :no-inherited-members:

Parameter terms
---------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Parameters.del_parameter
   ~cfdm.core.Parameters.get_parameter
   ~cfdm.core.Parameters.has_parameter
   ~cfdm.core.Parameters.set_parameter
   ~cfdm.core.Parameters.parameters
   ~cfdm.core.Parameters.clear_parameters
   ~cfdm.core.Parameters.set_parameters
   
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Parameters.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Parameters.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.Parameters._docstring_special_substitutions
   ~cfdm.core.Parameters._docstring_substitutions        
   ~cfdm.core.Parameters._docstring_package_depth        
   ~cfdm.core.Parameters._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Array
===============

----

.. autoclass:: cfdm.core.Array
   :no-members:
   :no-inherited-members:

Inspection
----------
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Array.dtype
   ~cfdm.core.Array.ndim
   ~cfdm.core.Array.shape
   ~cfdm.core.Array.size
   ~cfdm.core.Array.array

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Array.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Array.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Array._docstring_special_substitutions
   ~cfdm.core.Array._docstring_substitutions        
   ~cfdm.core.Array._docstring_package_depth        
   ~cfdm.core.Array._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Domain
================

----

.. autoclass:: cfdm.core.Domain
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Domain.construct_type
   
Metadata constructs
-------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Domain.del_construct
   ~cfdm.core.Domain.get_construct
   ~cfdm.core.Domain.has_construct
   ~cfdm.core.Domain.set_construct
   ~cfdm.core.Domain.del_data_axes
   ~cfdm.core.Domain.get_data_axes
   ~cfdm.core.Domain.has_data_axes
   ~cfdm.core.Domain.set_data_axes

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Domain.constructs
   
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Domain.del_property
   ~cfdm.core.Domain.get_property
   ~cfdm.core.Domain.has_property
   ~cfdm.core.Domain.set_property
   ~cfdm.core.Domain.properties
   ~cfdm.core.Domain.clear_properties
   ~cfdm.core.Domain.set_properties

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Domain.copy
   ~cfdm.core.Domain.fromconstructs
   ~cfdm.core.Domain.has_bounds
   ~cfdm.core.Domain.has_data

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Domain.__deepcopy__

Docstring substitutions
-----------------------                           
                                                  
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Domain._docstring_special_substitutions
   ~cfdm.core.Domain._docstring_substitutions        
   ~cfdm.core.Domain._docstring_package_depth        
   ~cfdm.core.Domain._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.InteriorRing
======================

----

.. autoclass:: cfdm.core.InteriorRing
   :no-members:
   :no-inherited-members:

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.InteriorRing.del_property
   ~cfdm.core.InteriorRing.get_property
   ~cfdm.core.InteriorRing.has_property
   ~cfdm.core.InteriorRing.set_property
   ~cfdm.core.InteriorRing.properties
   ~cfdm.core.InteriorRing.clear_properties
   ~cfdm.core.InteriorRing.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.InteriorRing.del_data
   ~cfdm.core.InteriorRing.get_data
   ~cfdm.core.InteriorRing.has_data
   ~cfdm.core.InteriorRing.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.InteriorRing.data

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.InteriorRing.copy
   ~cfdm.core.InteriorRing.has_bounds

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.InteriorRing.__deepcopy__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.core.InteriorRing._docstring_special_substitutions
   ~cfdm.core.InteriorRing._docstring_substitutions        
   ~cfdm.core.InteriorRing._docstring_package_depth        
   ~cfdm.core.InteriorRing._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-CellMethod:

cfdm.CellMethod
===============

----

.. autoclass:: cfdm.CellMethod
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMethod.dump
   ~cfdm.CellMethod.identity  
   ~cfdm.CellMethod.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.CellMethod.construct_type

Properties
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMethod.del_axes
   ~cfdm.CellMethod.get_axes 
   ~cfdm.CellMethod.has_axes
   ~cfdm.CellMethod.set_axes
   ~cfdm.CellMethod.del_method
   ~cfdm.CellMethod.get_method
   ~cfdm.CellMethod.has_method
   ~cfdm.CellMethod.set_method
   ~cfdm.CellMethod.del_qualifier
   ~cfdm.CellMethod.get_qualifier
   ~cfdm.CellMethod.has_qualifier
   ~cfdm.CellMethod.set_qualifier
   ~cfdm.CellMethod.qualifiers

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMethod.copy
   ~cfdm.CellMethod.creation_commands
   ~cfdm.CellMethod.equals
   ~cfdm.CellMethod.sorted

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CellMethod.__deepcopy__
   ~cfdm.CellMethod.__repr__
   ~cfdm.CellMethod.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.CellMethod._docstring_special_substitutions
   ~cfdm.CellMethod._docstring_substitutions        
   ~cfdm.CellMethod._docstring_package_depth        
   ~cfdm.CellMethod._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-AuxiliaryCoordinate:

cfdm.AuxiliaryCoordinate
========================

----

.. autoclass:: cfdm.AuxiliaryCoordinate
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.dump
   ~cfdm.AuxiliaryCoordinate.identity
   ~cfdm.AuxiliaryCoordinate.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.AuxiliaryCoordinate.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.del_property
   ~cfdm.AuxiliaryCoordinate.get_property
   ~cfdm.AuxiliaryCoordinate.has_property
   ~cfdm.AuxiliaryCoordinate.set_property
   ~cfdm.AuxiliaryCoordinate.properties
   ~cfdm.AuxiliaryCoordinate.clear_properties
   ~cfdm.AuxiliaryCoordinate.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.apply_masking
   ~cfdm.AuxiliaryCoordinate.del_data
   ~cfdm.AuxiliaryCoordinate.get_data
   ~cfdm.AuxiliaryCoordinate.has_data
   ~cfdm.AuxiliaryCoordinate.set_data
   ~cfdm.AuxiliaryCoordinate.insert_dimension
   ~cfdm.AuxiliaryCoordinate.squeeze
   ~cfdm.AuxiliaryCoordinate.transpose

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.AuxiliaryCoordinate.data
   ~cfdm.AuxiliaryCoordinate.dtype
   ~cfdm.AuxiliaryCoordinate.ndim
   ~cfdm.AuxiliaryCoordinate.shape
   ~cfdm.AuxiliaryCoordinate.size

Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.del_bounds
   ~cfdm.AuxiliaryCoordinate.get_bounds
   ~cfdm.AuxiliaryCoordinate.has_bounds
   ~cfdm.AuxiliaryCoordinate.set_bounds
   ~cfdm.AuxiliaryCoordinate.get_bounds_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.AuxiliaryCoordinate.bounds

Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.del_geometry
   ~cfdm.AuxiliaryCoordinate.get_geometry
   ~cfdm.AuxiliaryCoordinate.has_geometry
   ~cfdm.AuxiliaryCoordinate.set_geometry
   ~cfdm.AuxiliaryCoordinate.del_interior_ring
   ~cfdm.AuxiliaryCoordinate.get_interior_ring
   ~cfdm.AuxiliaryCoordinate.has_interior_ring
   ~cfdm.AuxiliaryCoordinate.set_interior_ring
   ~cfdm.AuxiliaryCoordinate.del_node_count
   ~cfdm.AuxiliaryCoordinate.get_node_count
   ~cfdm.AuxiliaryCoordinate.has_node_count
   ~cfdm.AuxiliaryCoordinate.set_node_count
   ~cfdm.AuxiliaryCoordinate.del_part_node_count
   ~cfdm.AuxiliaryCoordinate.get_part_node_count
   ~cfdm.AuxiliaryCoordinate.has_part_node_count
   ~cfdm.AuxiliaryCoordinate.set_part_node_count
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.AuxiliaryCoordinate.interior_ring

Climatology
^^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.del_climatology
   ~cfdm.AuxiliaryCoordinate.get_climatology
   ~cfdm.AuxiliaryCoordinate.is_climatology
   ~cfdm.AuxiliaryCoordinate.set_climatology

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst	     

   ~cfdm.AuxiliaryCoordinate.copy
   ~cfdm.AuxiliaryCoordinate.creation_commands
   ~cfdm.AuxiliaryCoordinate.equals
   ~cfdm.AuxiliaryCoordinate.uncompress
   ~cfdm.AuxiliaryCoordinate.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.nc_del_variable
   ~cfdm.AuxiliaryCoordinate.nc_get_variable
   ~cfdm.AuxiliaryCoordinate.nc_has_variable
   ~cfdm.AuxiliaryCoordinate.nc_set_variable

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.nc_variable_groups
   ~cfdm.AuxiliaryCoordinate.nc_clear_variable_groups
   ~cfdm.AuxiliaryCoordinate.nc_set_variable_groups

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.AuxiliaryCoordinate.__deepcopy__
   ~cfdm.AuxiliaryCoordinate.__getitem__
   ~cfdm.AuxiliaryCoordinate.__repr__
   ~cfdm.AuxiliaryCoordinate.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.AuxiliaryCoordinate._docstring_special_substitutions
   ~cfdm.AuxiliaryCoordinate._docstring_substitutions        
   ~cfdm.AuxiliaryCoordinate._docstring_package_depth        
   ~cfdm.AuxiliaryCoordinate._docstring_method_exclusions
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Bounds
================

----

.. autoclass:: cfdm.core.Bounds
   :no-members:
   :no-inherited-members:

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Bounds.del_property
   ~cfdm.core.Bounds.get_property
   ~cfdm.core.Bounds.has_property
   ~cfdm.core.Bounds.set_property
   ~cfdm.core.Bounds.properties
   ~cfdm.core.Bounds.clear_properties
   ~cfdm.core.Bounds.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Bounds.del_data
   ~cfdm.core.Bounds.get_data
   ~cfdm.core.Bounds.has_data
   ~cfdm.core.Bounds.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Bounds.data

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Bounds.copy
   ~cfdm.core.Bounds.has_bounds

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Bounds.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Bounds_special_substitutions
   ~cfdm.core.Bounds_substitutions        
   ~cfdm.core.Bounds_package_depth        
   ~cfdm.core.Bounds_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.CellMeasure
=====================

----

.. autoclass:: cfdm.core.CellMeasure
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.CellMeasure.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMeasure.del_measure  
   ~cfdm.core.CellMeasure.get_measure
   ~cfdm.core.CellMeasure.has_measure
   ~cfdm.core.CellMeasure.set_measure
   ~cfdm.core.CellMeasure.del_property
   ~cfdm.core.CellMeasure.get_property
   ~cfdm.core.CellMeasure.has_property
   ~cfdm.core.CellMeasure.set_property
   ~cfdm.core.CellMeasure.properties
   ~cfdm.core.CellMeasure.clear_properties
   ~cfdm.core.CellMeasure.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMeasure.del_data
   ~cfdm.core.CellMeasure.get_data
   ~cfdm.core.CellMeasure.has_data
   ~cfdm.core.CellMeasure.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.CellMeasure.data

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMeasure.copy
   ~cfdm.core.CellMeasure.has_bounds

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMeasure.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.CellMeasure._docstring_special_substitutions
   ~cfdm.core.CellMeasure._docstring_substitutions        
   ~cfdm.core.CellMeasure._docstring_package_depth        
   ~cfdm.core.CellMeasure._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-Domain:

cfdm.Domain
===========

----

.. autoclass:: cfdm.Domain
   :no-members:
   :no-inherited-members:

.. rubric:: Methods

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.dump
   ~cfdm.Domain.identity  
   ~cfdm.Domain.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Domain.construct_type

Metadata constructs
-------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.construct
   ~cfdm.Domain.construct_key
   ~cfdm.Domain.construct_item
   ~cfdm.Domain.del_construct
   ~cfdm.Domain.get_construct
   ~cfdm.Domain.has_construct
   ~cfdm.Domain.set_construct
   ~cfdm.Domain.del_data_axes
   ~cfdm.Domain.get_data_axes
   ~cfdm.Domain.has_data_axes
   ~cfdm.Domain.set_data_axes
   ~cfdm.Domain.domain_axis_key
   ~cfdm.Domain.auxiliary_coordinates
   ~cfdm.Domain.cell_measures
   ~cfdm.Domain.coordinates
   ~cfdm.Domain.coordinate_references
   ~cfdm.Domain.dimension_coordinates
   ~cfdm.Domain.domain_ancillaries
   ~cfdm.Domain.domain_axes
   ~cfdm.Domain.climatological_time_axes

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Domain.constructs

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.del_property
   ~cfdm.Domain.get_property
   ~cfdm.Domain.has_property
   ~cfdm.Domain.set_property
   ~cfdm.Domain.properties
   ~cfdm.Domain.clear_properties
   ~cfdm.Domain.set_properties

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.apply_masking
   ~cfdm.Domain.climatological_time_axes
   ~cfdm.Domain.copy
   ~cfdm.Domain.creation_commands
   ~cfdm.Domain.equals
   ~cfdm.Domain.fromconstructs
   ~cfdm.Domain.get_filenames
   ~cfdm.Domain.has_bounds
   ~cfdm.Domain.has_data
   ~cfdm.Domain.has_geometry
   ~cfdm.Domain.apply_masking   
   ~cfdm.Domain.get_filenames

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.nc_del_variable
   ~cfdm.Domain.nc_get_variable
   ~cfdm.Domain.nc_has_variable
   ~cfdm.Domain.nc_set_variable 
   ~cfdm.Domain.nc_global_attributes
   ~cfdm.Domain.nc_clear_global_attributes
   ~cfdm.Domain.nc_set_global_attribute
   ~cfdm.Domain.nc_set_global_attributes
   
Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.nc_variable_groups
   ~cfdm.Domain.nc_set_variable_groups
   ~cfdm.Domain.nc_clear_variable_groups
   ~cfdm.Domain.nc_group_attributes
   ~cfdm.Domain.nc_clear_group_attributes
   ~cfdm.Domain.nc_set_group_attribute
   ~cfdm.Domain.nc_set_group_attributes
  
Geometries
^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Domain.nc_del_geometry_variable
   ~cfdm.Domain.nc_get_geometry_variable
   ~cfdm.Domain.nc_has_geometry_variable
   ~cfdm.Domain.nc_set_geometry_variable 
   ~cfdm.Domain.nc_geometry_variable_groups
   ~cfdm.Domain.nc_set_geometry_variable_groups
   ~cfdm.Domain.nc_clear_geometry_variable_groups

Components
^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.Domain.nc_del_component_variable
   ~cfdm.Domain.nc_set_component_variable
   ~cfdm.Domain.nc_set_component_variable_groups
   ~cfdm.Domain.nc_clear_component_variable_groups      
   ~cfdm.Domain.nc_del_component_dimension
   ~cfdm.Domain.nc_set_component_dimension
   ~cfdm.Domain.nc_set_component_dimension_groups
   ~cfdm.Domain.nc_clear_component_dimension_groups
   ~cfdm.Domain.nc_del_component_sample_dimension
   ~cfdm.Domain.nc_set_component_sample_dimension   
   ~cfdm.Domain.nc_set_component_sample_dimension_groups
   ~cfdm.Domain.nc_clear_component_sample_dimension_groups

Dataset compliance
^^^^^^^^^^^^^^^^^^

.. rubric:: Methods


.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.dataset_compliance
   
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Domain.__deepcopy__
   ~cfdm.Domain.__repr__
   ~cfdm.Domain.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Domain._docstring_special_substitutions
   ~cfdm.Domain._docstring_substitutions        
   ~cfdm.Domain._docstring_package_depth        
   ~cfdm.Domain._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-CoordinateReference:

cfdm.CoordinateReference
========================

----

.. autoclass:: cfdm.CoordinateReference
   :no-members:
   :no-inherited-members:
   
Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.dump
   ~cfdm.CoordinateReference.identity  
   ~cfdm.CoordinateReference.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.CoordinateReference.construct_type

Coordinates
-----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.del_coordinate
   ~cfdm.CoordinateReference.has_coordinate
   ~cfdm.CoordinateReference.set_coordinate
   ~cfdm.CoordinateReference.coordinates
   ~cfdm.CoordinateReference.clear_coordinates
   ~cfdm.CoordinateReference.set_coordinates

Datum
-----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.del_datum
   ~cfdm.CoordinateReference.get_datum
   ~cfdm.CoordinateReference.set_datum

.. rubric:: Attributes
	    
.. autosummary::
   :toctree: ../attribute/
   :template: attribute.rst
	      
   ~cfdm.CoordinateReference.datum

Coordinate conversion
---------------------
 
.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.del_coordinate_conversion
   ~cfdm.CoordinateReference.get_coordinate_conversion
   ~cfdm.CoordinateReference.set_coordinate_conversion

.. rubric:: Attributes
	    
.. autosummary::
   :toctree: ../attribute/
   :template: attribute.rst
	      
   ~cfdm.CoordinateReference.coordinate_conversion

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.copy
   ~cfdm.CoordinateReference.creation_commands
   ~cfdm.CoordinateReference.equals

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.nc_del_variable
   ~cfdm.CoordinateReference.nc_get_variable
   ~cfdm.CoordinateReference.nc_has_variable
   ~cfdm.CoordinateReference.nc_set_variable

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.nc_variable_groups
   ~cfdm.CoordinateReference.nc_set_variable_groups
   ~cfdm.CoordinateReference.nc_clear_variable_groups
      
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.CoordinateReference.__deepcopy__
   ~cfdm.CoordinateReference.__repr__
   ~cfdm.CoordinateReference.__str__

Docstring substitutions
-----------------------                      
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.CoordinateReference._docstring_special_substitutions
   ~cfdm.CoordinateReference._docstring_substitutions        
   ~cfdm.CoordinateReference._docstring_package_depth        
   ~cfdm.CoordinateReference._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Constant
=============

----

.. autoclass:: cfdm.Constant
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst
   
   ~cfdm.Constant.value

Copying
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.Constant.copy
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst
   
   ~cfdm.Constant._func

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Constant.__deepcopy__
   ~cfdm.Constant.__enter__
   ~cfdm.Constant.__exit__
   ~cfdm.Constant.__repr__
   ~cfdm.Constant.__str__

.. currentmodule:: cfdm
.. default-role:: obj

cfdm.List
=========

----

.. autoclass:: cfdm.List
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.dump
   ~cfdm.List.identity  
   ~cfdm.List.identities

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.del_property
   ~cfdm.List.get_property
   ~cfdm.List.has_property
   ~cfdm.List.set_property
   ~cfdm.List.properties
   ~cfdm.List.clear_properties
   ~cfdm.List.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.apply_masking
   ~cfdm.List.del_data
   ~cfdm.List.get_data
   ~cfdm.List.has_data
   ~cfdm.List.set_data   
   ~cfdm.List.insert_dimension
   ~cfdm.List.squeeze
   ~cfdm.List.transpose
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.List.data
   ~cfdm.List.dtype
   ~cfdm.List.ndim    
   ~cfdm.List.shape
   ~cfdm.List.size    

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.copy
   ~cfdm.List.creation_commands
   ~cfdm.List.equals
   ~cfdm.List.get_filenames
   ~cfdm.List.has_bounds
   ~cfdm.List.uncompress

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.nc_del_variable
   ~cfdm.List.nc_get_variable
   ~cfdm.List.nc_has_variable
   ~cfdm.List.nc_set_variable 

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.List.nc_variable_groups
   ~cfdm.List.nc_clear_variable_groups
   ~cfdm.List.nc_set_variable_groups
   
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.List.__deepcopy__
   ~cfdm.List.__getitem__
   ~cfdm.List.__repr__
   ~cfdm.List.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.List._docstring_special_substitutions
   ~cfdm.List._docstring_substitutions        
   ~cfdm.List._docstring_package_depth        
   ~cfdm.List._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.AuxiliaryCoordinate
=============================

----

.. autoclass:: cfdm.core.AuxiliaryCoordinate
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.AuxiliaryCoordinate.construct_type

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.del_property
   ~cfdm.core.AuxiliaryCoordinate.get_property
   ~cfdm.core.AuxiliaryCoordinate.has_property
   ~cfdm.core.AuxiliaryCoordinate.set_property
   ~cfdm.core.AuxiliaryCoordinate.properties
   ~cfdm.core.AuxiliaryCoordinate.clear_properties
   ~cfdm.core.AuxiliaryCoordinate.set_properties

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.del_data
   ~cfdm.core.AuxiliaryCoordinate.get_data
   ~cfdm.core.AuxiliaryCoordinate.has_data
   ~cfdm.core.AuxiliaryCoordinate.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.AuxiliaryCoordinate.data

Bounds
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.del_bounds
   ~cfdm.core.AuxiliaryCoordinate.get_bounds
   ~cfdm.core.AuxiliaryCoordinate.has_bounds
   ~cfdm.core.AuxiliaryCoordinate.set_bounds
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.AuxiliaryCoordinate.bounds

Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.del_geometry
   ~cfdm.core.AuxiliaryCoordinate.get_geometry
   ~cfdm.core.AuxiliaryCoordinate.has_geometry
   ~cfdm.core.AuxiliaryCoordinate.set_geometry
   ~cfdm.core.AuxiliaryCoordinate.del_interior_ring
   ~cfdm.core.AuxiliaryCoordinate.get_interior_ring
   ~cfdm.core.AuxiliaryCoordinate.has_interior_ring
   ~cfdm.core.AuxiliaryCoordinate.set_interior_ring

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.AuxiliaryCoordinate.interior_ring

Climatology
^^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.del_climatology
   ~cfdm.core.AuxiliaryCoordinate.get_climatology
   ~cfdm.core.AuxiliaryCoordinate.is_climatology
   ~cfdm.core.AuxiliaryCoordinate.set_climatology

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.AuxiliaryCoordinate.__deepcopy__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.core.AuxiliaryCoordinate._docstring_special_substitutions
   ~cfdm.core.AuxiliaryCoordinate._docstring_substitutions        
   ~cfdm.core.AuxiliaryCoordinate._docstring_package_depth        
   ~cfdm.core.AuxiliaryCoordinate._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Container
===================

----

.. autoclass:: cfdm.core.Container
   :no-members:
   :no-inherited-members:

Miscellaneous
-------------

.. rubric:: Methods
	   
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Container.copy

Private
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Container._del_component
   ~cfdm.core.Container._get_component
   ~cfdm.core.Container._has_component
   ~cfdm.core.Container._set_component

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Container.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.Container._docstring_special_substitutions
   ~cfdm.core.Container._docstring_substitutions        
   ~cfdm.core.Container._docstring_package_depth        
   ~cfdm.core.Container._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Coordinate
====================

----

.. autoclass:: cfdm.core.Coordinate
   :no-members:
   :no-inherited-members:
   
Properties
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
      
   ~cfdm.core.Coordinate.del_property
   ~cfdm.core.Coordinate.get_property
   ~cfdm.core.Coordinate.has_property
   ~cfdm.core.Coordinate.set_property
   ~cfdm.core.Coordinate.properties
   ~cfdm.core.Coordinate.clear_properties
   ~cfdm.core.Coordinate.set_properties

Data
----

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.del_data
   ~cfdm.core.Coordinate.get_data
   ~cfdm.core.Coordinate.has_data
   ~cfdm.core.Coordinate.set_data
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Coordinate.data

Bounds
------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.del_bounds
   ~cfdm.core.Coordinate.get_bounds
   ~cfdm.core.Coordinate.has_bounds
   ~cfdm.core.Coordinate.set_bounds
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Coordinate.bounds
   
Geometries
^^^^^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.del_geometry
   ~cfdm.core.Coordinate.get_geometry
   ~cfdm.core.Coordinate.has_geometry
   ~cfdm.core.Coordinate.set_geometry
   ~cfdm.core.Coordinate.del_interior_ring
   ~cfdm.core.Coordinate.get_interior_ring
   ~cfdm.core.Coordinate.has_interior_ring
   ~cfdm.core.Coordinate.set_interior_ring

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Coordinate.interior_ring

Climatology
^^^^^^^^^^^

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.del_climatology
   ~cfdm.core.Coordinate.get_climatology
   ~cfdm.core.Coordinate.is_climatology
   ~cfdm.core.Coordinate.set_climatology

Modification
------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.del_data
   ~cfdm.core.Coordinate.set_data
   ~cfdm.core.Coordinate.del_bounds
   ~cfdm.core.Coordinate.set_bounds
      
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Coordinate.__deepcopy__

Docstring substitutions
-----------------------                                 
                                                        
.. rubric:: Methods                                     
                                                        
.. autosummary::                                        
   :nosignatures:                                       
   :toctree: ../method/                                 
   :template: method.rst                                
                                                        
   ~cfdm.core.Coordinate._docstring_special_substitutions
   ~cfdm.core.Coordinate._docstring_substitutions        
   ~cfdm.core.Coordinate._docstring_package_depth        
   ~cfdm.core.Coordinate._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.DocstringRewriteMeta
==============================

----

.. autoclass:: cfdm.core.DocstringRewriteMeta
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.core.DocstringRewriteMeta.mro
   ~cfdm.core.DocstringRewriteMeta._docstring_special_substitutions
   ~cfdm.core.DocstringRewriteMeta._docstring_substitutions
   ~cfdm.core.DocstringRewriteMeta._docstring_package_depth
   ~cfdm.core.DocstringRewriteMeta._docstring_method_exclusions

.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Field
===============

----

.. autoclass:: cfdm.core.Field
   :no-members:
   :no-inherited-members:

.. _core-Field-Inspection:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Field.construct_type

.. _core-Field-Properties:

Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.del_property
   ~cfdm.core.Field.get_property
   ~cfdm.core.Field.has_property
   ~cfdm.core.Field.set_property
   ~cfdm.core.Field.properties
   ~cfdm.core.Field.clear_properties
   ~cfdm.core.Field.set_properties

.. _core-Field-Data:

Data
----

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.del_data
   ~cfdm.core.Field.get_data
   ~cfdm.core.Field.has_data
   ~cfdm.core.Field.set_data
   ~cfdm.core.Field.del_data_axes
   ~cfdm.core.Field.get_data_axes
   ~cfdm.core.Field.has_data_axes
   ~cfdm.core.Field.set_data_axes
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Field.data

.. _core-Field-Metadata-constructs:   
   
Metadata constructs
-------------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.del_construct
   ~cfdm.core.Field.get_construct
   ~cfdm.core.Field.has_construct
   ~cfdm.core.Field.set_construct
   ~cfdm.core.Field.del_data_axes
   ~cfdm.core.Field.get_data_axes
   ~cfdm.core.Field.has_data_axes
   ~cfdm.core.Field.set_data_axes

.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Field.constructs

.. _core-Field-Domain:

Domain
------


.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.get_domain
   
.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.Field.domain

.. _core-Field-Miscellaneous:

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.copy
   ~cfdm.core.Field.has_bounds

.. _core-Field-Special:

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Field.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Field._docstring_special_substitutions
   ~cfdm.core.Field._docstring_substitutions        
   ~cfdm.core.Field._docstring_package_depth        
   ~cfdm.core.Field._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.Array
==========

----

.. autoclass:: cfdm.Array
   :no-members:
   :no-inherited-members:

Inspection
----------
   
.. rubric:: Attributes

.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.Array.array
   ~cfdm.Array.dtype
   ~cfdm.Array.ndim
   ~cfdm.Array.shape
   ~cfdm.Array.size

Compression
-----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Array.get_compression_type
   
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Array.copy
   ~cfdm.Array.get_subspace

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.Array.__array__
   ~cfdm.Array.__deepcopy__
   ~cfdm.Array.__getitem__
   ~cfdm.Array.__repr__
   ~cfdm.Array.__str__

Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.Array._docstring_special_substitutions
   ~cfdm.Array._docstring_substitutions        
   ~cfdm.Array._docstring_package_depth        
   ~cfdm.Array._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.CellMethod
====================

----

.. autoclass:: cfdm.core.CellMethod
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.core.CellMethod.construct_type
 
Properties
----------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMethod.del_axes
   ~cfdm.core.CellMethod.get_axes 
   ~cfdm.core.CellMethod.has_axes
   ~cfdm.core.CellMethod.set_axes
   ~cfdm.core.CellMethod.del_method
   ~cfdm.core.CellMethod.get_method
   ~cfdm.core.CellMethod.has_method
   ~cfdm.core.CellMethod.set_method
   ~cfdm.core.CellMethod.del_qualifier
   ~cfdm.core.CellMethod.get_qualifier
   ~cfdm.core.CellMethod.has_qualifier
   ~cfdm.core.CellMethod.set_qualifier
   ~cfdm.core.CellMethod.qualifiers

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMethod.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.CellMethod.__deepcopy__

Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.CellMethod._docstring_special_substitutions
   ~cfdm.core.CellMethod._docstring_substitutions        
   ~cfdm.core.CellMethod._docstring_package_depth        
   ~cfdm.core.CellMethod._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

cfdm.core.Datum
===============

----

.. autoclass:: cfdm.core.Datum
   :no-members:
   :no-inherited-members:

Parameters
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Datum.del_parameter
   ~cfdm.core.Datum.get_parameter
   ~cfdm.core.Datum.has_parameter
   ~cfdm.core.Datum.set_parameter
   ~cfdm.core.Datum.parameters
   ~cfdm.core.Datum.clear_parameters
   ~cfdm.core.Datum.set_parameters

Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Datum.copy

Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.core.Datum.__deepcopy__
   
Docstring substitutions
-----------------------                        
                                               
.. rubric:: Methods                            
                                               
.. autosummary::                               
   :nosignatures:                              
   :toctree: ../method/                        
   :template: method.rst                       
                                               
   ~cfdm.core.Datum._docstring_special_substitutions
   ~cfdm.core.Datum._docstring_substitutions        
   ~cfdm.core.Datum._docstring_package_depth        
   ~cfdm.core.Datum._docstring_method_exclusions    
.. currentmodule:: cfdm
.. default-role:: obj

.. _cfdm-DomainAxis:

cfdm.DomainAxis
===============

----

.. autoclass:: cfdm.DomainAxis
   :no-members:
   :no-inherited-members:

Inspection
----------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.DomainAxis.identity  
   ~cfdm.DomainAxis.identities

.. rubric:: Attributes
   
.. autosummary::
   :nosignatures:
   :toctree: ../attribute/
   :template: attribute.rst

   ~cfdm.DomainAxis.construct_type

Size
----

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
   
   ~cfdm.DomainAxis.del_size
   ~cfdm.DomainAxis.has_size
   ~cfdm.DomainAxis.get_size
   ~cfdm.DomainAxis.set_size
   
Miscellaneous
-------------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAxis.copy
   ~cfdm.DomainAxis.creation_commands
   ~cfdm.DomainAxis.equals

NetCDF
------

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAxis.nc_del_dimension
   ~cfdm.DomainAxis.nc_get_dimension
   ~cfdm.DomainAxis.nc_has_dimension
   ~cfdm.DomainAxis.nc_set_dimension 
   ~cfdm.DomainAxis.nc_is_unlimited
   ~cfdm.DomainAxis.nc_set_unlimited

Groups
^^^^^^

.. rubric:: Methods
	    
.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst
	      
   ~cfdm.DomainAxis.nc_dimension_groups
   ~cfdm.DomainAxis.nc_clear_dimension_groups
   ~cfdm.DomainAxis.nc_set_dimension_groups
   
Special
-------

.. rubric:: Methods

.. autosummary::
   :nosignatures:
   :toctree: ../method/
   :template: method.rst

   ~cfdm.DomainAxis.__deepcopy__
   ~cfdm.DomainAxis.__repr__
   ~cfdm.DomainAxis.__str__
   
Docstring substitutions
-----------------------                   
                                          
.. rubric:: Methods                       
                                          
.. autosummary::                          
   :nosignatures:                         
   :toctree: ../method/                   
   :template: method.rst                  
                                          
   ~cfdm.DomainAxis._docstring_special_substitutions
   ~cfdm.DomainAxis._docstring_substitutions        
   ~cfdm.DomainAxis._docstring_package_depth        
   ~cfdm.DomainAxis._docstring_method_exclusions    
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autofunction:: {{ fullname }}
:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. automethod:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autoattribute:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. auto{{ objtype }}:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autofunction:: {{ fullname }}
:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. automethod:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autoattribute:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. auto{{ objtype }}:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autofunction:: {{ fullname }}
:orphan:

{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. automethod:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. currentmodule:: cfdm
.. default-role:: obj

.. autoattribute:: {{ fullname }}
{{ fullname }}
{{ underline }}

.. auto{{ objtype }}:: {{ fullname }}
