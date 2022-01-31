---
title: 'ObsPlus: A Pandas-centric ObsPy expansion pack'
tags:
  - Python
  - seismology
  - pandas
authors:
  - name: Derrick J.A. Chambers
    orcid: 0000-0003-3656-6607
    affiliation: 1
  - name: M. Shawn Boltz
    affiliation: 1
  - name: Calum J. Chamberlain
    orcid: 0000-0003-2317-2609
    affiliation: 2
affiliations:
 - name: National Institute for Occupational Safety and Health, Spokane Mining Research Division
   index: 1
 - name: School of Geography, Environment and Earth Sciences, Victoria University of Wellington, New Zealand
   index: 2
date: 03 March 2020
bibliography: paper.bib
---

# Summary

Over the past decade, ``ObsPy``, a python framework for seismology [@Krischer:2015], has become an integral part of many seismology research workflows. ``ObsPy`` provides parsers for most seismological data formats, clients for accessing data-centers, common signal processing routines, and event, station, and waveform data models.

``ObsPlus`` significantly expands ``ObsPy``’s functionality by providing simple data management abstractions and conversions between ``ObsPy`` classes and the ubiquitous ``pandas`` ``DataFrame`` [@mckinney-proc-scipy-2010].

# Statement of Need

ObsPlus benefits researchers by 1) simplifying seismological data access patterns and management of local data, 2) providing means to extract facets of seismic catalog hierarchies to tabular forms, and 3) enabling the packaging and distribution of complete seismological datasets.

# Functionality and Features

1. **A data retrieval interface**
``ObsPlus`` provides a unified data retrieval interface for in-memory, on disk, and remote seismological data. This is enabled by in-process databases which provide a simple mechanism to index and access local seismological data stored in directories of arbitrary organization. Importantly, the classes implement a superset of the interface already provided by ``ObsPy``’s remote clients, making it straight-forward to write data-source agnostic code.

2. **Alternative data structures**
While ``ObsPy``'s data structures are quite powerful, they are not always the most convenient. For example, the canonical event data representation is a deeply nested tree-like structure based on the QuakeML standard [@Schorlemmer:2011]. Working with events in ``ObsPy`` often necessitates deeply nested recursive code which can become difficult to understand and maintain. ``ObsPlus`` provides functionality to flatten desired components of these tree structures into ``DataFrame``s which are simpler and more efficient when the full complexity of QuakeML isn’t merrited for the task at hand.

3. **Datasets**
``ObsPlus`` provides a simple mechanism to bundle, distribute, download, and interact with complete seismological datasets. This is done by creating a simple python package (for which we provide a cookie cutter template) which is published to PyPI. Each package includes small files and instructions to download large files. Optionally, a list of files and their corresponding hashes can be used to validate downloaded data. Datasets are then discovered and loaded through python’s plugin system, and downloaded when needed.

4. **Utilities**
``ObsPlus``’s list of utilities is quite long and more are being added regularly. Many are focused around validating and manipulating event data.

``ObsPlus`` has become an important part of the National Institute for Occupational Safety and Health (NIOSH)’s data processing and management workflows and has enabled rapid prototyping of new ideas while managing complexity through appropriate abstractions. It is our hope that `ObsPlus` will provide similar benefits to the broader seismology community.

# References
* obsplus version:
* Python version:
* Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
---
name: 'Development discussion '
about: Suggest a feature or discuss future development
title: ''
labels: develop
assignees: ''

---

**Is this discussion related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
---
name: Bug report
about: Report a bug or unexpected behavior
title: ''
labels: bug
assignees: ''

---

**Description**
A clear and concise description of the bug.

**To Reproduce**
Include a short code example which reproduces the issue.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Versions (please complete the following information):**
 - OS: [e.g. iOS]
 - ObsPlus Version [e.g. 0.1.0]
 - Python Version [e.g. 3.8]
---
name: Question
about: Ask a question
title: ''
labels: question
assignees: ''

---
.. image:: https://user-images.githubusercontent.com/11671536/51515070-22241c00-1dc7-11e9-90d5-832f3caa5c70.png

"A Pandas-Centric ObsPy_ Expansion Pack"

|Coverage| |Supported Versions| |PyPI| |Conda| |DOI| |Licence|

Documentation_

Code_

Change_Log_

Contributing_

Support_

License_

About_

Code_of_Conduct_

If you find ObsPlus useful consider citing it:

Chambers, D. J., Boltz, M. S., & Chamberlain, C. J. (2021).
ObsPlus: A Pandas-centric ObsPy expansion pack.
Journal of Open Source Software, 6(60), 2696.


.. _About: https://github.com/niosh-mining/about
.. _ObsPy: https://github.com/obspy/obspy
.. _Documentation: https://niosh-mining.github.io/obsplus/versions/latest/index.html
.. _Support: https://niosh-mining.github.io/obsplus/versions/latest/notebooks/support.html
.. _Code: https://github.com/niosh-mining/obsplus
.. _Change_Log: https://github.com/niosh-mining/obsplus/CHANGELOG.txt
.. _License: https://choosealicense.com/licenses/lgpl-3.0/
.. _Code_of_Conduct: https://github.com/niosh-mining/obsplus/blob/master/.github/CODE_OF_CONDUCT.md
.. _Contributing: https://niosh-mining.github.io/obsplus/versions/latest/notebooks/contributing.html

.. |Coverage| image:: https://codecov.io/gh/niosh-mining/obsplus/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/niosh-mining/obsplus

.. |Supported Versions| image:: https://img.shields.io/pypi/pyversions/obsplus.svg
   :target: https://pypi.python.org/pypi/obsplus

.. |Licence| image:: https://www.gnu.org/graphics/lgplv3-88x31.png
   :target: https://www.gnu.org/licenses/lgpl.html

.. |PyPI| image:: https://pepy.tech/badge/obsplus
   :target: https://pepy.tech/project/obsplus

.. |Conda| image:: https://img.shields.io/conda/dn/conda-forge/obsplus?label=conda%20downloads
   :target: https://github.com/conda-forge/obsplus-feedstock

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4544008.svg
   :target: https://doi.org/10.5281/zenodo.4544008
.. image:: images/obsplus.png

ObsPlus
=======

"A Pandas-Centric ObsPy Expansion Pack"

`Introduction <notebooks/intro.ipynb>`_

`Installation <notebooks/installation.ipynb>`_

`Code <https://github.com/niosh-mining/obsplus>`_

`Contributing <notebooks/contributing.ipynb>`_

`Support Bug Reporting <notebooks/support.ipynb>`_

`Quick Reference <quickref/index.rst>`_

`Full API Documentation <api/obsplus.rst>`_

:ref:`genindex`

:ref:`modindex`

:ref:`search`
:orphan:
.. _quickref:

Quick References
################

A brief glance at ObsPlus' main features.

Events
======

.. autosummary::
     :toctree: stubs

     obsplus.get_event_client
     obsplus.EventBank
     obsplus.events_to_df
     obsplus.picks_to_df
     obsplus.json_to_cat
     obsplus.cat_to_json


Stations
========

.. autosummary::
     :toctree: stubs

     obsplus.get_station_client
     obsplus.stations_to_df



Waveforms
=========

.. autosummary::
     :toctree: stubs

     obsplus.get_waveform_client
     obsplus.WaveBank


Datasets
========

.. autosummary::
     :toctree: stubs

     obsplus.datasets.dataset.DataSet
     obsplus.load_dataset
     obsplus.Fetcher


Utils
=====

.. autosummary::
     :toctree: stubs

     obsplus.DataFrameExtractor
     obsplus.utils.yield_obj_parent_attr
     obsplus.utils.time.to_utc
     obsplus.utils.time.to_datetime64
