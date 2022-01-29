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
