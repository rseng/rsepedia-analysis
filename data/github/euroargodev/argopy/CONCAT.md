|<img src="https://raw.githubusercontent.com/euroargodev/argopy/master/docs/_static/argopy_logo_long.png" alt="argopy logo" width="400"/>|``argopy`` is a python library that aims to ease Argo data access, visualisation and manipulation for regular users as well as Argo experts and operators|
|:---------:|:-------|
|Documentation|[![JOSS](https://img.shields.io/badge/DOI-10.21105%2Fjoss.02425-brightgreen)](//dx.doi.org/10.21105/joss.02425) <br>[![Documentation](https://img.shields.io/static/v1?label=&message=Read%20the%20documentation&color=blue&logo=read-the-docs&logoColor=white)](https://argopy.readthedocs.io) <br>[![Gitter](https://badges.gitter.im/Argo-floats/argopy.svg)](https://gitter.im/Argo-floats/argopy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)|
|Usage|![License](https://img.shields.io/github/license/euroargodev/argopy) [![Python version](https://img.shields.io/pypi/pyversions/argopy)](//pypi.org/project/argopy/)<br>[![pypi dwn](https://img.shields.io/pypi/dm/argopy?label=Pypi%20downloads)](//pypi.org/project/argopy/) [![conda dwn](https://img.shields.io/conda/dn/conda-forge/argopy?label=Conda%20downloads)](//anaconda.org/conda-forge/argopy)|
|Release|[![](https://img.shields.io/github/release-date/euroargodev/argopy)](//github.com/euroargodev/argopy/releases) [![PyPI](https://img.shields.io/pypi/v/argopy)](//pypi.org/project/argopy/) [![Conda](https://anaconda.org/conda-forge/argopy/badges/version.svg)](//anaconda.org/conda-forge/argopy)|
|Development|![Github Action Status](https://github.com/euroargodev/argopy/workflows/tests/badge.svg?branch=master) [![Documentation Status](https://readthedocs.org/projects/argopy/badge/?version=latest)](https://argopy.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/euroargodev/argopy/branch/master/graph/badge.svg)](https://codecov.io/gh/euroargodev/argopy)<br>[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)|
|Data resources|![Erddap status](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/euroargodev/argopy-status/master/argopy_api_status_erddap.json) ![Argovis status](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/euroargodev/argopy-status/master/argopy_api_status_argovis.json) <br>![Profile count](https://img.shields.io/endpoint?label=Number%20of%20Argo%20profiles%3A&style=social&url=https%3A%2F%2Fapi.ifremer.fr%2Fargopy%2Fdata%2FARGO-FULL.json) <br>[![Statuspage](https://img.shields.io/static/v1?label=&message=Check%20all%20monitors&color=blue&logo=statuspage&logoColor=white)](https://argopy.statuspage.io)|

## Install


Install the last release with conda:
```bash
conda install -c conda-forge argopy
```
or pip:
```bash
pip install argopy
```
or conda:
```bash
conda install -c conda-forge argopy
```

But since this is a young library in active development, use direct install from this repo to benefit from the latest version:

```bash
pip install git+http://github.com/euroargodev/argopy.git@master
```

The ``argopy`` library is tested to work under most OS (Linux, Mac) and with python versions 3.7 and 3.8.

## Usage

[![badge](https://img.shields.io/static/v1.svg?logo=Jupyter&label=Pangeo+Binder&message=Click+here+to+try+argopy+online+!&color=blue&style=for-the-badge)](https://mybinder.org/v2/gh/euroargodev/argopy/master?labpath=docs/tryit.ipynb)

### Fetching Argo Data

Import the data fetcher:
```python
from argopy import DataFetcher as ArgoDataFetcher
```
and then, set it up to request data for a **specific space/time domain**:
```python
argo_loader = ArgoDataFetcher().region([-85,-45,10.,20.,0,10.])
argo_loader = ArgoDataFetcher().region([-85,-45,10.,20.,0,1000.,'2012-01','2012-12'])
```
for **profiles of a given float**: 
```python
argo_loader = ArgoDataFetcher().profile(6902746, 34)
argo_loader = ArgoDataFetcher().profile(6902746, np.arange(12,45))
argo_loader = ArgoDataFetcher().profile(6902746, [1,12])
```
or for **one or a collection of floats**:
```python
argo_loader = ArgoDataFetcher().float(6902746)
argo_loader = ArgoDataFetcher().float([6902746, 6902747, 6902757, 6902766])
```

Once your fetcher is initialized you can trigger fetch/load data like this:
```python
ds = argo_loader.to_xarray()  # or:
ds = argo_loader.load().data
```
By default fetched data are returned in memory as [xarray.DataSet](http://xarray.pydata.org/en/stable/data-structures.html#dataset). 
From there, it is easy to convert it to other formats like a [Pandas dataframe](https://pandas.pydata.org/pandas-docs/stable/getting_started/dsintro.html#dataframe):
```python
df = ArgoDataFetcher().profile(6902746, 34).load().data.to_dataframe()
```

or to export it to files:
```python
ds = ArgoDataFetcher().region([-85,-45,10.,20.,0,100.]).to_xarray()
ds.to_netcdf('my_selection.nc')
# or by profiles:
ds.argo.point2profile().to_netcdf('my_selection.nc')
```


### Fetching only Argo index
Argo index are returned as pandas dataframe. Index fetchers works similarly to data fetchers.

Load the Argo index fetcher:
```python
    from argopy import IndexFetcher as ArgoIndexFetcher
```
then, set it up to request index for a **specific space/time domain**:
```python
    index_loader = ArgoIndexFetcher().region([-85,-45,10.,20.])
    index_loader = ArgoIndexFetcher().region([-85,-45,10.,20.,'2012-01','2014-12'])
```
or for **one or a collection of floats**:
```python
    index_loader = ArgoIndexFetcher().float(6902746)
    index_loader = ArgoIndexFetcher().float([6902746, 6902747, 6902757, 6902766])   
```
Once your fetcher is initialized you can trigger fetch/load index like this:
```python
    df = index_loader.to_dataframe()  # or
    df = index_loader.load().index
```

Note that like the data fetcher, the index fetcher can use different sources, a local copy of the GDAC ftp for instance:
```python
    index_fetcher = ArgoIndexFetcher(src='localftp', path_ftp='/path/to/your/argo/ftp/', index_file='ar_index_global_prof.txt')
```

### Visualisation
For plottings methods, you'll need `matplotlib` and possibly `cartopy` and `seaborn` installed.
Argo Data and Index fetchers provide direct plotting methods, for instance:
```python    
    ArgoDataFetcher().float([6902745, 6902746]).plot('trajectory')    
```
![index_traj](https://github.com/euroargodev/argopy/raw/master/docs/_static/trajectory_sample.png)

See the [documentation page for more examples](https://argopy.readthedocs.io/en/latest/visualisation.html).

## Development roadmap

See milestone here: https://github.com/euroargodev/argopy/milestone/3# Contributor Covenant Code of Conduct

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
reported by contacting the project team at  euroargo@ifremer.fr. All
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
https://www.contributor-covenant.org/faqArgopy's contributor guidelines [can be found in our online documentation](https://argopy.readthedocs.io/en/latest/contributing.html)
1. [ ] Create a new branch and PR to prepare for release:

     ```git checkout -b releaseX.Y```

   1. [ ] Run codespell from repo root and fix errors:

         ```codespell -q 2```

   2. [ ] Make sure that all CI tests are passed with **free* environments

   3. [ ] Update ``./requirements.txt`` and ``./docs/requirements.txt`` with CI free environments dependencies versions

   4. [ ] Update ``./ci/requirements/py*-dev.yml`` with last free environments dependencies versions

   5. [ ] Make sure that all CI tests are passed with **dev* environments

   6. [ ] Increase release version in ``./setup.py`` file

   7. [ ] Update date and release version in ``./docs/whats-new.rst``

   8. [ ] Merge PR to master:
   
        ```git checkout master```
   
        ```git merge releaseX.Y```

2. [ ] On the master branch, commit the release in git:

      ```git commit -a -m 'Release v0.X.Y'```

3. [ ] Tag the release:

      ```git tag -a v0.X.Y -m 'v0.X.Y'```

4. [ ] Push it online:

      ```git push origin v0.X.Y```

5. [ ] Issue the release on GitHub. Click on "Draft a new release" at
    https://github.com/euroargodev/argopy/releases. Type in the version number ``v0.X.Y``, and click on the auto-generate a release note button.
    
    This will trigger the publish Github action that will push the release on Pypi.---
name: Fix bug
about: 'bug'
title: 'Fix a bug or issue'
labels: 'bug'
assignees: ''
---

<!-- Summarize in one line the bug being fixed in this PR -->

This PR will close: #  

#### Changes proposed in this pull request

-
-

#### Task list
- [] CI tests added or updated---
name: New Feature
about: 'API'
title: 'Implement a new feature'
labels: 'enhancement'
assignees: ''
---

<!-- Description of the new feature -->


#### New feature Code Sample
<!-- Describe the new feature API here -->

```python
# Your code here

```

#### Changes proposed in this pull request

-
-

#### Task list
- [] Closes #. 
- [] CI tests added
- [] If necessary, update dependencies in ``requirements.txt`` that are required for this new feature
- [] Fully documented, including ``whats-new.rst`` for all changes and ``api.rst`` for new API---
name: Prepare next release
about: 'release'
title: 'Prepare next release'
labels: 'release'
assignees: ''
---

# Prepare release

- [ ] Make sure that all [CI tests are passed with *free* environments](https://github.com/euroargodev/argopy/actions/workflows/pythonFREEtests.yml?query=event%3Apull_request)
- [ ] Update ``./requirements.txt`` and ``./docs/requirements.txt`` with CI free environments dependencies versions 
- [ ] Update ``./ci/requirements/py*-dev.yml`` with last free environments dependencies versions
- [ ] Make sure that all [CI tests are passed with *dev* environments](https://github.com/euroargodev/argopy/actions/workflows/pythontests.yml?query=event%3Apull_request)
- [ ] Increase release version in ``./setup.py`` file
- [ ] Update date and release version in ``./docs/whats-new.rst``
- [ ] Merge this PR to master

# Publish release

- [ ] On the master branch, commit the release in git:
      ```git commit -a -m 'Release v0.X.Y'```
- [ ] Tag the release:
      ```git tag -a v0.X.Y -m 'v0.X.Y'```
- [ ] Push it online:
       ```git push origin v0.X.Y```
- [ ] Issue the release on GitHub. Click on "Draft a new release" at
     https://github.com/euroargodev/argopy/releases. Type in the version number, but
     don't bother to describe it -- we maintain that on the docs instead.
      
     This will trigger the [publish Github action](https://github.com/euroargodev/argopy/blob/master/.github/workflows/pythonpublish.yml) that will push the release on Pypi.---
name: Bug report / Feature request
about: 'Post a problem or idea'
title: ''
labels: ''
assignees: ''

---

<!-- A short summary of the issue, if appropriate -->


#### MCVE Code Sample
<!-- In order for the maintainers to efficiently understand and prioritize issues, we ask you post a "Minimal, Complete and Verifiable Example" (MCVE): http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports -->

```python
# Your code here

```

#### Expected Output


#### Problem Description
<!-- this should explain why the current behavior is a problem and why the expected output is a better solution -->


#### Versions

<details><summary>Output of `argopy.show_versions()`</summary>

<!-- Paste the output here argopy.show_versions() here -->


</details>
---
title: 'argopy: A Python library for Argo ocean data analysis'
tags:
  - Python
  - ocean
  - oceanography
  - observation
authors:
  - name: Guillaume Maze
    orcid: 0000-0001-7231-2095
    affiliation: 1
  - name: Kevin Balem
    orcid: 0000-0002-4956-8698
    affiliation: 1
affiliations:
 - name: Univ Brest, Ifremer, CNRS, IRD, LOPS, F‐29280 Plouzané, France
   index: 1
date: 24 June 2020
bibliography: paper.bib

---

# Summary

Argo is a real-time global ocean *in situ* observing system. It provides thousands of highly accurate ocean measurements 
every day. The Argo dataset has now accumulated more than 2.3 million vertical ocean profiles and accessing it for scientific 
analysis remains a challenge.

The Argo expert community, focused on delivering a curated dataset of the best scientific quality possible, has never provided 
its user base with a Python software package to easily access and manipulate Argo measurements: the **argopy** software aims 
to fill this gap. The **argopy** software can be used to easily fetch and manipulate measurements from Argo floats. 
It is dedicated to scientists without knowledge of the Argo data management system but is also designed to accommodate expert 
requirements.

# Introduction

The ocean is a key component of the Earth's climate system. It therefore needs continuous real-time monitoring to help scientists 
better understand its dynamics and to predict its evolution. All around the world, oceanographers have managed to join their
efforts and set up a [Global Ocean Observing System](https://www.goosocean.org/) among which *Argo* is a key component. 

Argo is a global network of nearly 4000 autonomous probes measuring pressure, temperature and salinity from the surface 
to 2000m depth every 10 days. The localisation of these probes is nearly random between the $60^o$ parallels ([see live 
coverage here](http://map.argo-france.fr)). Data from the probes are collected by satellite in real-time, processed by several 
data centers, merged in a single dataset (comprising of more than 2.3 million vertical profiles as of 
June 2020) and made freely available to anyone through an [ftp server](ftp://ftp.ifremer.fr/ifremer/argo) or [monthly zip 
snapshots](http://dx.doi.org/10.17882/42182).

The Argo international observation array was initiated in 1999 and soon revolutionised our 
perspective on the large scale structure and variability of the ocean by providing seasonally and regionally unbiased 
in situ temperature/salinity measurements of the ocean interior, key information that satellites can't provide [@riser-2016]. 
The Argo array reached its full global coverage (of 1 profile per month and per 3x3 degree horizontal area) in 2007, and 
pursues its evolution to fulfil new scientific requirements [@roemmich-2019]. Argo data have been used in more than 4000 scientific publications.

This [online figure](http://map.argo-france.fr) shows the current coverage of the network. It now extends to higher latitudes than the 
original $\pm60^o$ and some of the floats are able to profile down to 4000m and 6000m. New floats are also equipped 
with biogeochemical sensors, measuring oxygen and chlorophyll for instance. All these evolutions of the network increase 
the total number of floats to nearly 4000. Argo is thus providing a deluge of in situ data: more than 400 profiles per day.

Each Argo probe is an autonomous, free drifting, profiling float, i.e. a probe that can't control its trajectory but 
is able to control its buoyancy and thus to move up and down the water column as it wishes. Argo floats continuously 
operate the same program, or cycle, illustrated \autoref{fig:argofloat}. After 9 to 10 days of free drift at a parking 
depth of about 1000m, a typical Argo float dives down to 2000m and then rises back to the surface while profiling - measuring pressure, 
temperature and salinity. Once it reaches the surface, the float sends by satellite its measurements to a data center, 
where they are processed in real time and made freely available on the web in less than 24 hours.

![Typical 10 days program, cycle, of an Argo float.\label{fig:argofloat}](_static/argofloats_cycle.png)


# Why **argopy** ?

For non-expert users of the Argo dataset, it is rather complicated to get access to Argo measurements. Even though data are
made freely available on the web, the Argo dataset consists of thousands of files organised using jargon, 
tens of different variables and many reference tables. The exhaustive Argo [user manual](http://dx.doi.org/10.13155/29825) 
is more than 100 pages long, which can be rather intimidating to go through for new users.

This complexity arises from the fact that Argo operates many different models of floats and sensors, quality control 
of *in situ* measurements from autonomous platforms requires a lot of complementary information (meta-data), and the 
Argo data management workflow is distributed between more than 10 Data Assembly Centers all around the world. The Argo 
data management is a model for other ocean observing systems and constantly ensures the highest quality of scientific 
measurements for the community [@wong-2020].

The result of this tremendous success in data management -- in developing good practices and well calibrated 
procedures ([see all the Argo Data Management Team documentation here](http://www.argodatamgt.org/Documentation)) -- is 
a very complex Argo dataset: the **argopy** software aims to help users navigate this complex realm.

Since the Argo community focuses on delivering a curated dataset for science, software packages exist for Argo data operators to decode and quality control the data [e.g. @scoop]. However, no open source software are available for scientists, who therefore must develop their own machinery to download and manipulate the data.

Python is becoming widely used by the scientific community and beyond: worldwide, and is the most popular and fastest growing language in the last 5 years (20%, source: http://pypl.github.io/PYPL.html). It offers a modern, powerful and open
source framework to work with. Since, up to this point, no Python based software has been dedicated to the Argo dataset, it made sense to develop **argopy**.

# Key features of **argopy**

**argopy** is a python software package that simplifies the process of accessing and manipulating Argo data.
The two key features of **argopy** are its trivial fetching API of Argo data and its 
ability to provide data formatted for both beginner and expert users of Argo.

## Data fetching

**argopy** provides a trivial fetching API of Argo data through a simple call to one of the 3 different ways to 
look at Argo data: over a space/time domain (with the *region* access point), for one or a list of specific floats (given 
their unique [WMO number](https://www.wmo.int/pages/prog/amp/mmop/wmo-number-rules.html) with the *float* access point) 
or for one or a list of float profiles (with the *profile* access point). This is as simple as:
```python
from argopy import DataFetcher as ArgoDataFetcher
fetcher = ArgoDataFetcher().region([-75, -45, 20, 30, 0, 100, '2011', '2012'])
ds = fetcher.to_xarray()
```
Here we used **argopy** to fetch data between 75/45W, 20/30N, from 0 to 100db and for the entire year 2011.
Once the user has defined what they need (the ``fetcher`` class instance in the example above), **argopy** will fetch data online and manage 
internally all the complicated processing of formatting the web request and creating a workable in memory data 
structure (the ``to_xarray()`` call above). By default, **argopy** uses the [xarray data model](http://xarray.pydata.org);
*xarray* is an open source Python package to easily work with labelled multi-dimensional arrays.

## Data formatting

**argopy** aims to thrive in providing Argo data to non-experts. One key feature of **argopy** is the option for selecting
a *user mode* that is either ``standard`` or ``expert``. Standard users are those who want to focus on the measurements 
for scientific analysis; those who do not know, or don't want to be bothered with all the Argo jargon and multitude of 
variables and parameters. 

For standard users (the default mode), **argopy** internally runs a series of processes that
curate the raw data and provide a simplified and science focused dataset. For expert users, **argopy** will apply its 
data model to raw fetched data and return Argo variables that experts users are already used to.

## And more

**argopy** has features to manipulate Argo data, for instance:

- the possibility to transform data from a collection of measurements to a collection of vertical profiles, and vice-versa; 
- the possibility to interpolate irregularly sampled measurements onto standard pressure levels.
 
Another feature is the ability to cache fetched data, so that requests provide users with data much more rapidly, 
saving bandwidth and time. 

Two last important features of **argopy** to mention here are: 

- the possibility to fetch data locally, from a user copy of the entire or subset of the Argo database,
- the possibility to fetch only meta data (organised in *index* lookup tables), which allows the user to determine the regional Argo sampling, 
for instance.
These more advanced features may be more of interest for ``expert`` users, since they require more knowledge of the Argo dataset.

# Conclusion

**argopy** is filling an important gap in the ocean science community by providing an easy way to access a large and complex dataset that has proved to be very important in oceanographic studies. For information on all the features available with **argopy**, the reader is referred to the complete software documentation at [https://argopy.readthedocs.io](https://argopy.readthedocs.io).

# Acknowledgements

We acknowledge support from the [Euro-Argo ERIC community](https://www.euro-argo.eu/) during the genesis of this project.
This software was created with support from the EARISE project, a European Union’s Horizon 2020 research and 
innovation programme under grant agreement no 824131. Call INFRADEV-03-2018-2019: Individual support to ESFRI and other 
world-class research infrastructures.

# References
**********************
Contributing to argopy
**********************

.. contents:: Table of contents:
   :local:

First off, thanks for taking the time to contribute!

.. note::

  Large parts of this document came from the `Xarray <http://xarray.pydata.org/en/stable/contributing.html>`_
  and `Pandas <http://pandas.pydata.org/pandas*docs/stable/contributing.html>`_ contributing guides.

If you seek **support** for your argopy usage or if you don't want to read
this whole thing and just have a question: `visit the chat room at gitter <https://gitter.im/Argo-floats/argopy>`_.

Where to start?
===============

All contributions, bug reports, bug fixes, documentation improvements,
enhancements, and ideas are welcome.

We will complete this document for guidelines with regard to each of these contributions over time.

If you are brand new to *argopy* or open source development, we recommend going
through the `GitHub "issues" tab <https://github.com/euroargodev/argopy/issues>`_
to find issues that interest you. There are a number of issues listed under
`Documentation <https://github.com/euroargodev/argopy/issues?q=is%3Aissue+is%3Aopen+label%3Adocumentation>`_
and `good first issue
<https://github.com/euroargodev/argopy/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22>`_
where you could start out. Once you've found an interesting issue, you can
return here to get your development environment setup.

Please don't file an issue to ask a question, instead `visit the chat room at gitter <https://gitter.im/Argo-floats/argopy>`_.

.. _contributing.bug_reports:

Bug reports and enhancement requests
====================================

Bug reports are an important part of making *argopy* more stable. Having a complete bug
report will allow others to reproduce the bug and provide insight into fixing. See
`this stackoverflow article <https://stackoverflow.com/help/mcve>`_ for tips on
writing a good bug report.

Trying the bug producing code out on the *master* branch is often a worthwhile exercise
to confirm the bug still exists. It is also worth searching existing bug reports and
pull requests to see if the issue has already been reported and/or fixed.

Bug reports must:

#. Include a short, self contained Python snippet reproducing the problem.
   You can format the code nicely by using `GitHub Flavored Markdown
   <http://github.github.com/github*flavored*markdown/>`_::

      ```python
      >>> import argopy as ar
      >>> ds = ar.DataFetcher(backend='erddap').float(5903248).to_xarray()
      ...
      ```

#. Include the full version string of *argopy* and its dependencies. You can use the
   built in function::

      >>> import argopy
      >>> argopy.show_versions()

#. Explain why the current behavior is wrong/not desired and what you expect instead.

The issue will then show up to the :mod:`argopy` community and be open to comments/ideas
from others.

`Click here to open an issue with the specific bug reporting template <https://github.com/euroargodev/argopy/issues/new?template=bug_report.md>`_


.. _contributing.documentation:

Contributing to the documentation
=================================

If you're not the developer type, contributing to the documentation is still of
huge value. You don't even have to be an expert on *argopy* to do so! In fact,
there are sections of the docs that are worse off after being written by
experts. If something in the docs doesn't make sense to you, updating the
relevant section after you figure it out is a great way to ensure it will help
the next person.

.. contents:: Documentation:
   :local:


About the *argopy* documentation
--------------------------------

The documentation is written in **reStructuredText**, which is almost like writing
in plain English, and built using `Sphinx <http://sphinx-doc.org/>`__. The
Sphinx Documentation has an excellent `introduction to reST
<http://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`__. Review the Sphinx docs to perform more
complex changes to the documentation as well.

Some other important things to know about the docs:

- The *argopy* documentation consists of two parts: the docstrings in the code
  itself and the docs in this folder ``argopy/docs/``.

  The docstrings are meant to provide a clear explanation of the usage of the
  individual functions, while the documentation in this folder consists of
  tutorial-like overviews per topic together with some other information
  (what's new, installation, etc).

- The docstrings follow the **Numpy Docstring Standard**, which is used widely
  in the Scientific Python community. This standard specifies the format of
  the different sections of the docstring. See `this document
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_
  for a detailed explanation, or look at some of the existing functions to
  extend it in a similar manner.

- The tutorials make heavy use of the `ipython directive
  <http://matplotlib.org/sampledoc/ipython_directive.html>`_ sphinx extension.
  This directive lets you put code in the documentation which will be run
  during the doc build. For example:

  .. code:: rst

      .. ipython:: python

          x = 2
          x ** 3

  will be rendered as::

      In [1]: x = 2

      In [2]: x ** 3
      Out[2]: 8

  Almost all code examples in the docs are run (and the output saved) during the
  doc build. This approach means that code examples will always be up to date,
  but it does make the doc building a bit more complex.

- Our API documentation in ``docs/api.rst`` houses the auto-generated
  documentation from the docstrings. For classes, there are a few subtleties
  around controlling which methods and attributes have pages auto-generated.

  Every method should be included in a ``toctree`` in ``api.rst``, else Sphinx
  will emit a warning.


How to build the *argopy* documentation
---------------------------------------

Requirements
^^^^^^^^^^^^
Make sure to follow the instructions on :ref:`creating a development environment below <contributing.dev_env>`, but
to build the docs you need to use the specific file ``docs/requirements.txt``:

.. code-block:: bash

    $ conda create --yes -n argopy-docs python=3.6 xarray dask numpy pytest future gsw sphinx sphinx_rtd_theme
    $ conda activate argopy-docs
    $ pip install argopy
    $ pip install -r docs/requirements.txt

Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Navigate to your local ``argopy/docs/`` directory in the console and run:

.. code-block:: bash

    make html

Then you can find the HTML output in the folder ``argopy/docs/_build/html/``.

The first time you build the docs, it will take quite a while because it has to run
all the code examples and build all the generated docstring pages. In subsequent
evocations, sphinx will try to only build the pages that have been modified.

If you want to do a full clean build, do:

.. code-block:: bash

    make clean
    make html


.. _working.code:

Working with the code
=====================

Development workflow
--------------------

Anyone interested in helping to develop argopy needs to create their own fork
of our `git repository`. (Follow the github `forking instructions`_. You
will need a github account.)

.. _git repository: https://github.com/euroargodev/argopy
.. _forking instructions: https://help.github.com/articles/fork-a-repo/

Clone your fork on your local machine.

.. code-block:: bash

    $ git clone git@github.com:USERNAME/argopy

(In the above, replace USERNAME with your github user name.)

Then set your fork to track the upstream argopy repo.

.. code-block:: bash

    $ cd argopy
    $ git remote add upstream git://github.com/euroargodev/argopy.git

You will want to periodically sync your master branch with the upstream master.

.. code-block:: bash

    $ git fetch upstream
    $ git rebase upstream/master

**Never make any commits on your local master branch**. Instead open a feature
branch for every new development task.

.. code-block:: bash

    $ git checkout -b cool_new_feature

(Replace `cool_new_feature` with an appropriate description of your feature.)
At this point you work on your new feature, using `git add` to add your
changes. When your feature is complete and well tested, commit your changes

.. code-block:: bash

    $ git commit -m 'did a bunch of great work'

and push your branch to github.

.. code-block:: bash

    $ git push origin cool_new_feature

At this point, you go find your fork on github.com and create a `pull
request`_. Clearly describe what you have done in the comments. If your
pull request fixes an issue or adds a useful new feature, the team will
gladly merge it.

.. _pull request: https://help.github.com/articles/using-pull-requests/

After your pull request is merged, you can switch back to the master branch,
rebase, and delete your feature branch. You will find your new feature
incorporated into argopy.

.. code-block:: bash

    $ git checkout master
    $ git fetch upstream
    $ git rebase upstream/master
    $ git branch -d cool_new_feature

.. _contributing.dev_env:

Virtual environment
-------------------

This is how to create a virtual environment into which to test-install argopy,
install it, check the version, and tear down the virtual environment.

.. code-block:: bash

    $ conda create --yes -n argopy-tests python=3.6 xarray dask numpy pytest future gsw
    $ conda activate argopy-tests
    $ pip install argopy
    $ python -c 'import argopy; print(argopy.__version__);'
    $ conda deactivate
    $ conda env remove --yes -n argopy-tests


Code standards
--------------

Writing good code is not just about what you write. It is also about *how* you
write it. During Continuous Integration testing, several
tools will be run to check your code for stylistic errors.
Generating any warnings will cause the test to fail.
Thus, good style is a requirement for submitting code to *argopy*.

Code Formatting
---------------

*argopy* uses several tools to ensure a consistent code format throughout the project:

* `Flake8 <http://flake8.pycqa.org/en/latest/>`_ for general code quality

``pip``::

   pip install flake8

and then run from the root of the argopy repository::

   flake8

to qualify your code.


.. _contributing.code:

Contributing to the code base
=============================

.. contents:: Code Base:
   :local:

.. _data_fetchers:

Data fetchers
-------------

Introduction
^^^^^^^^^^^^
If you want to add your own data fetcher for a new service, then, keep in mind that:

* Data fetchers are responsible for:

  * loading all available data from a given source and providing at least a :func:`to_xarray()` method
  * making data compliant to Argo standards (data type, variable name, attributes, etc ...)

* Data fetchers must:

  * inherit from the :class:`argopy.data_fetchers.proto.ArgoDataFetcherProto`
  * provide parameters:

    *  ``access_points``, eg: ['wmo', 'box']
    *  ``exit_formats``, eg: ['xarray']
    *  ``dataset_ids``, eg: ['phy', 'ref', 'bgc']

  * provides the facade API (:class:`argopy.fetchers.ArgoDataFetcher`) methods to filter data
    according to user level or requests. These must includes:

    *  :func:`filter_data_mode`
    *  :func:`filter_qc`
    *  :func:`filter_variables`


It is the responsibility of the facade API (:class:`argopy.fetchers.ArgoDataFetcher`) to run
filters according to user level or requests, not the data fetcher.

Detailed guideline
^^^^^^^^^^^^^^^^^^

A new data fetcher must comply with:

Inheritance
"""""""""""

Inherit from the :class:`argopy.data_fetchers.proto.ArgoDataFetcherProto`.
This enforces minimal internal design compliance.

Auto-discovery of fetcher properties
""""""""""""""""""""""""""""""""""""

The new fetcher must come with the ``access_points``, ``exit_formats`` and ``dataset_ids`` properties at the top of the
file, e.g.:

.. code-block:: python

    access_points = ['wmo' ,'box']
    exit_formats = ['xarray']
    dataset_ids = ['phy', 'bgc']  # First is default

Values depend on what the new access point can return and what you want to
implement. A good start is with the ``wmo`` access point and the
``phy`` dataset ID. The ``xarray`` data format is the minimum
required. These variables are used by the facade
to auto-discover the fetcher capabilities. The ``dataset_ids``
property is used to determine which variables can be retrieved.

Auto-discovery of fetcher access points
"""""""""""""""""""""""""""""""""""""""

The new fetcher must come at least with a ``Fetch_box`` or
``Fetch_wmo`` class, basically one for each of the ``access_points``
listed as properties. More generally we may have a main class that
provides the key functionality to retrieve data from the source,
and then classes for each of the ``access_points`` of your fetcher.
This pattern could look like this:

.. code-block:: python

    class NewDataFetcher(ArgoDataFetcherProto)
    class Fetch_wmo(NewDataFetcher)
    class Fetch_box(NewDataFetcher)

It could also be like:

.. code-block:: python

    class Fetch_wmo(ArgoDataFetcherProto)
    class Fetch_box(ArgoDataFetcherProto)

Note that the class names ``Fetch_wmo`` and ``Fetch_box`` must not
change, this is also used by the facade to auto-discover the fetcher
capabilities.

**Fetch\_wmo** is used to retrieve platforms and eventually profiles
data. It must take in the ``__init__()`` method a ``WMO`` and a ``CYC``
as first and second options. ``WMO`` is always passed, ``CYC`` is
optional. These are passed by the facade to implement the
``fetcher.float`` and ``fetcher.profile`` methods. When a float is requested, the ``CYC`` option is
not passed by the facade. Last, ``WMO`` and ``CYC`` are either a single
integer or a list of integers: this means that ``Fetch_wmo`` must be
able to handle more than one float/platform retrieval.

**Fetch\_box** is used to retrieve a rectangular domain in space and
time. It must take in the ``__init__()`` method a ``BOX`` as first
option that is passed a list(lon\_min: float, lon\_max: float, lat\_min:
float, lat\_max: float, pres\_min: float, pres\_max: float, date\_min:
str, date\_max: str) from the facade. The two bounding dates [date\_min
and date\_max] should be optional (if not specified, the entire time
series is requested by the user).

Internal File systems
"""""""""""""""""""""

All http requests must go through the internal
``httpstore``, an internal wrapper around fsspec that allows to
manage request caching very easily. You can simply use it this way
for json requests:

.. code-block:: python

    from argopy.stores import httpstore
    with httpstore(timeout=120).open("https://argovis.colorado.edu/catalog/profiles/5904797_12") as of:
       profile = json.load(of)

Output data format
""""""""""""""""""

Last but not least, about the output data. In **argopy**, we want
to provide data for both expert and standard users. This is explained
and illustrated in the `documentation
here <https://argopy.readthedocs.io/en/latest/user_mode.html>`__.
This means for a new data fetcher that the data content
should be curated and clean of any internal/jargon variables that is
not part of the Argo ADMT vocabulary. For instance,
variables like: ``bgcMeasKeys`` or ``geoLocation`` are not allowed. This will ensure
that whatever the data source set by users, the output xarray or
dataframe will be formatted and contain the same variables. This will
also ensure that other argopy features can be used on the new fetcher
output, like plotting or xarray data manipulation.

.. currentmodule:: argopy

What's New
==========

v0.1.9 (19 Jan. 2022)
---------------------

**Features and front-end API**

- **New method to preprocess data for OWC software**. This method can preprocessed Argo data and possibly create float_source/<WMO>.mat files to be used as inputs for OWC implementations in `Matlab <https://github.com/ArgoDMQC/matlab_owc>`_ and `Python <https://github.com/euroargodev/argodmqc_owc>`_. See the :ref:`Salinity calibration` documentation page for more. (:pr:`142`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    ds = ArgoDataFetcher(mode='expert').float(6902766).load().data
    ds.argo.create_float_source("float_source")
    ds.argo.create_float_source("float_source", force='raw')
    ds_source = ds.argo.create_float_source()


.. currentmodule:: xarray

This new method comes with others methods and improvements:

    - A new :meth:`Dataset.argo.filter_scalib_pres` method to filter variables according to OWC salinity calibration software requirements,
    - A new :meth:`Dataset.argo.groupby_pressure_bins` method to subsample a dataset down to one value by pressure bins (a perfect alternative to interpolation on standard depth levels to precisely avoid interpolation...), see :ref:`Pressure levels: Group-by bins` for more help,
    - An improved :meth:`Dataset.argo.filter_qc` method to select which fields to consider (new option ``QC_fields``),
    - Add conductivity (``CNDC``) to the possible output of the ``TEOS10`` method.

.. currentmodule:: argopy

- **New dataset properties** accessible from the `argo` xarray accessor: ``N_POINTS``, ``N_LEVELS``, ``N_PROF``. Note that depending on the format of the dataset (a collection of points or of profiles) these values do or do not take into account NaN. These information are also visible by a simple print of the accessor. (:pr:`142`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    ds = ArgoDataFetcher(mode='expert').float(6902766).load().data
    ds.argo.N_POINTS
    ds.argo.N_LEVELS
    ds.argo.N_PROF
    ds.argo
    

- **New plotter function** :meth:`argopy.plotters.open_sat_altim_report` to insert the CLS Satellite Altimeter Report figure in a notebook cell. (:pr:`159`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy.plotters import open_sat_altim_report
    open_sat_altim_report(6902766)
    open_sat_altim_report([6902766, 6902772, 6902914])
    open_sat_altim_report([6902766, 6902772, 6902914], embed='dropdown')  # Default
    open_sat_altim_report([6902766, 6902772, 6902914], embed='slide')
    open_sat_altim_report([6902766, 6902772, 6902914], embed='list')
    open_sat_altim_report([6902766, 6902772, 6902914], embed=None)

    from argopy import DataFetcher
    from argopy import IndexFetcher
    DataFetcher().float([6902745, 6902746]).plot('qc_altimetry')
    IndexFetcher().float([6902745, 6902746]).plot('qc_altimetry')


- **New utility method to retrieve topography**. The :class:`argopy.TopoFetcher` will load the `GEBCO topography <https://coastwatch.pfeg.noaa.gov/erddap/griddap/GEBCO_2020.html>`_ for a given region. (:pr:`150`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import TopoFetcher
    box = [-75, -45, 20, 30]
    ds = TopoFetcher(box).to_xarray()
    ds = TopoFetcher(box, ds='gebco', stride=[10, 10], cache=True).to_xarray()


For convenience we also added a new property to the data fetcher that return the domain covered by the dataset.

.. code-block:: python

    loader = ArgoDataFetcher().float(2901623)
    loader.domain  # Returns [89.093, 96.036, -0.278, 4.16, 15.0, 2026.0, numpy.datetime64('2010-05-14T03:35:00.000000000'),  numpy.datetime64('2013-01-01T01:45:00.000000000')]

- Update the documentation with a new section about :ref:`data_qc`.

**Internals**

- Uses a new API endpoint for the ``argovis`` data source when fetching a ``region``. `More on this issue here <https://github.com/donatagiglio/Argovis/issues/3>`_. (:pr:`158`) by `G. Maze <http://www.github.com/gmaze>`_.

- Update documentation theme, and pages now use the `xarray accessor sphinx extension <https://github.com/xarray-contrib/sphinx-autosummary-accessors>`_. (:pr:`104`) by `G. Maze <http://www.github.com/gmaze>`_.

- Update Binder links to work without the deprecated Pangeo-Binder service. (:pr:`164`) by `G. Maze <http://www.github.com/gmaze>`_.

v0.1.8 (2 Nov. 2021)
---------------------

**Features and front-end API**

- Improve plotting functions. All functions are now available for both the index and data fetchers. See the :ref:`data_viz` page for more details. Reduced plotting dependencies to `Matplotlib <https://matplotlib.org/>`_ only. **Argopy** will use `Seaborn <seaborn.pydata.org/>`_ and/or `Cartopy <https://scitools.org.uk/cartopy>`_ if available. (:pr:`56`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import IndexFetcher as ArgoIndexFetcher
    from argopy import DataFetcher as ArgoDataFetcher
    obj = ArgoIndexFetcher().float([6902766, 6902772, 6902914, 6902746])
    # OR
    obj = ArgoDataFetcher().float([6902766, 6902772, 6902914, 6902746])

    fig, ax = obj.plot()
    fig, ax = obj.plot('trajectory')
    fig, ax = obj.plot('trajectory', style='white', palette='Set1', figsize=(10,6))
    fig, ax = obj.plot('dac')
    fig, ax = obj.plot('institution')
    fig, ax = obj.plot('profiler')


- New methods and properties for data and index fetchers. (:pr:`56`) by `G. Maze <http://www.github.com/gmaze>`_. The :meth:`argopy.DataFetcher.load` and :meth:`argopy.IndexFetcher.load` methods internally call on the `to_xarray()` methods and store results in the fetcher instance. The :meth:`argopy.DataFetcher.to_xarray` will trigger a fetch on every call, while the :meth:`argopy.DataFetcher.load` will not.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    loader = ArgoDataFetcher().float([6902766, 6902772, 6902914, 6902746])
    loader.load()
    loader.data
    loader.index
    loader.to_index()

.. code-block:: python

    from argopy import IndexFetcher as ArgoIndexFetcher
    indexer = ArgoIndexFetcher().float([6902766, 6902772])
    indexer.load()
    indexer.index

- Add optional speed of sound computation to xarray accessor teos10 method. (:pr:`90`) by `G. Maze <http://www.github.com/gmaze>`_.

- Code spell fixes (:pr:`89`) by `K. Schwehr <https://github.com/schwehr>`_.

**Internals**

- Check validity of access points options (WMO and box) in the facade, no checks at the fetcher level. (:pr:`92`) by `G. Maze <http://www.github.com/gmaze>`_.

- More general options. Fix :issue:`91`. (:pr:`102`) by `G. Maze <http://www.github.com/gmaze>`_.

    - ``trust_env`` to allow for local environment variables to be used by fsspec to connect to the internet. Useful for those using a proxy.

- Documentation on `Read The Docs` now uses a pip environment and get rid of memory eager conda. (:pr:`103`) by `G. Maze <http://www.github.com/gmaze>`_.

- :class:`xarray.Dataset` argopy accessor ``argo`` has a clean documentation.

**Breaking changes with previous versions**

- Drop support for python 3.6 and older. Lock range of dependencies version support.

- In the plotters module, the ``plot_dac`` and ``plot_profilerType`` functions have been replaced by ``bar_plot``. (:pr:`56`) by `G. Maze <http://www.github.com/gmaze>`_.

**Internals**

- Internal logging available and upgrade dependencies version support (:pr:`56`) by `G. Maze <http://www.github.com/gmaze>`_. To see internal logs, you can set-up your application like this:

.. code-block:: python

    import logging
    DEBUGFORMATTER = '%(asctime)s [%(levelname)s] [%(name)s] %(filename)s:%(lineno)d: %(message)s'
    logging.basicConfig(
        level=logging.DEBUG,
        format=DEBUGFORMATTER,
        datefmt='%m/%d/%Y %I:%M:%S %p',
        handlers=[logging.FileHandler("argopy.log", mode='w')]
    )

v0.1.7 (4 Jan. 2021)
-----------------------

Long due release !

**Features and front-end API**

- Live monitor for the status (availability) of data sources. See documentation page on :ref:`api-status`. (:pr:`36`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    import argopy
    argopy.status()
    # or
    argopy.status(refresh=15)

.. image:: _static/status_monitor.png
  :width: 350

- Optimise large data fetching with parallelization, for all data fetchers (erddap, localftp and argovis). See documentation page on :ref:`parallel`. Two parallel methods are available: multi-threading or multi-processing. (:pr:`28`) by `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    loader = ArgoDataFetcher(parallel=True)
    loader.float([6902766, 6902772, 6902914, 6902746]).to_xarray()
    loader.region([-85,-45,10.,20.,0,1000.,'2012-01','2012-02']).to_xarray()


**Breaking changes with previous versions**

- In the teos10 xarray accessor, the ``standard_name`` attribute will now be populated using values from the `CF Standard Name table <https://cfconventions.org/Data/cf-standard-names/76/build/cf-standard-name-table.html>`_ if one exists.
  The previous values of ``standard_name`` have been moved to the ``long_name`` attribute.
  (:pr:`74`) by `A. Barna <https://github.com/docotak>`_.
  
- The unique resource identifier property is now named ``uri`` for all data fetchers, it is always a list of strings.

**Internals**

- New ``open_mfdataset`` and ``open_mfjson`` methods in Argo stores. These can be used to open, pre-process and concatenate a collection of paths both in sequential or parallel order. (:pr:`28`) by `G. Maze <http://www.github.com/gmaze>`_.

- Unit testing is now done on a controlled conda environment. This allows to more easily identify errors coming from development vs errors due to dependencies update. (:pr:`65`) by `G. Maze <http://www.github.com/gmaze>`_.


v0.1.6 (31 Aug. 2020)
---------------------

- **JOSS paper published**. You can now cite argopy with a clean reference. (:pr:`30`) by `G. Maze <http://www.github.com/gmaze>`_ and `K. Balem <http://www.github.com/quai20>`_.

Maze G. and Balem K. (2020). argopy: A Python library for Argo ocean data analysis. *Journal of Open Source Software*, 5(52), 2425 doi: `10.21105/joss.02425 <http://dx.doi.org/10.21105/joss.02425>`_.


v0.1.5 (10 July 2020)
---------------------

**Features and front-end API**

- A new data source with the **argovis** data fetcher, all access points available (:pr:`24`). By `T. Tucker <https://github.com/tylertucker202>`_ and `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    loader = ArgoDataFetcher(src='argovis')
    loader.float(6902746).to_xarray()
    loader.profile(6902746, 12).to_xarray()
    loader.region([-85,-45,10.,20.,0,1000.,'2012-01','2012-02']).to_xarray()

- Easily compute `TEOS-10 <http://teos-10.org/>`_ variables with new argo accessor function **teos10**. This needs `gsw <https://github.com/TEOS-10/GSW-Python>`_ to be installed. (:pr:`37`) By `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    ds = ArgoDataFetcher().region([-85,-45,10.,20.,0,1000.,'2012-01','2012-02']).to_xarray()
    ds = ds.argo.teos10()
    ds = ds.argo.teos10(['PV'])
    ds_teos10 = ds.argo.teos10(['SA', 'CT'], inplace=False)

- **argopy** can now be installed with conda (:pr:`29`, :pr:`31`, :pr:`32`). By `F. Fernandes <https://github.com/ocefpaf>`_.

.. code-block:: text

    conda install -c conda-forge argopy


**Breaking changes with previous versions**

- The ``local_ftp`` option of the ``localftp`` data source must now points to the folder where the ``dac`` directory is found. This breaks compatibility with rsynced local FTP copy because rsync does not give a ``dac`` folder (e.g. :issue:`33`). An instructive error message is raised to notify users if any of the DAC name is found at the n-1 path level. (:pr:`34`).

**Internals**

- Implement a webAPI availability check in unit testing. This allows for more robust ``erddap`` and ``argovis`` tests that are not only based on internet connectivity only. (:commit:`5a46a39a3368431c6652608ee7241888802f334f`).


v0.1.4 (24 June 2020)
---------------------

**Features and front-end API**

- Standard levels interpolation method available in **standard** user mode (:pr:`23`). By `K. Balem <http://www.github.com/quai20>`_.

.. code-block:: python

    ds = ArgoDataFetcher().region([-85,-45,10.,20.,0,1000.,'2012-01','2012-12']).to_xarray()
    ds = ds.argo.point2profile()
    ds_interp = ds.argo.interp_std_levels(np.arange(0,900,50))

- Insert in a Jupyter notebook cell the `Euro-Argo fleet monitoring <https://fleetmonitoring.euro-argo.eu>`_ dashboard page, possibly for a specific float (:pr:`20`). By `G. Maze <http://www.github.com/gmaze>`_.

.. code-block:: python

    import argopy
    argopy.dashboard()
    # or
    argopy.dashboard(wmo=6902746)

- The ``localftp`` index and data fetcher now have the ``region`` and ``profile`` access points available (:pr:`25`). By `G. Maze <http://www.github.com/gmaze>`_.

**Breaking changes with previous versions**

[None]

**Internals**

- Now uses `fsspec <https://filesystem-spec.readthedocs.io>`_ as file system for caching as well as accessing local and remote files (:pr:`19`). This closes issues :issue:`12`, :issue:`15` and :issue:`17`. **argopy** fetchers must now use (or implement if necessary) one of the internal file systems available in the new module ``argopy.stores``. By `G. Maze <http://www.github.com/gmaze>`_.

- Erddap fetcher now uses netcdf format to retrieve data (:pr:`19`).

v0.1.3 (15 May 2020)
--------------------

**Features and front-end API**

- New ``index`` fetcher to explore and work with meta-data (:pr:`6`). By `K. Balem <http://www.github.com/quai20>`_.

.. code-block:: python

    from argopy import IndexFetcher as ArgoIndexFetcher
    idx = ArgoIndexFetcher().float(6902746)
    idx.to_dataframe()
    idx.plot('trajectory')

The ``index`` fetcher can manage caching and works with both Erddap and localftp data sources. It is basically the same as the data fetcher, but do not load measurements, only meta-data. This can be very useful when looking for regional sampling or trajectories.

.. tip::

  **Performance**: we recommend to use the ``localftp`` data source when working this ``index`` fetcher because the ``erddap`` data source currently suffers from poor performances. This is linked to :issue:`16` and is being addressed by Ifremer.

The ``index`` fetcher comes with basic plotting functionalities with the :func:`argopy.IndexFetcher.plot` method to rapidly visualise measurement distributions by DAC, latitude/longitude and floats type.

.. warning::

  The design of plotting and visualisation features in ``argopy`` is constantly evolving, so this may change in future releases.

- Real documentation written and published (:pr:`13`). By `G. Maze <http://www.github.com/gmaze>`_.

- The :class:`argopy.DataFetcher` now has a :func:`argopy.DataFetcher.to_dataframe` method to return a :class:`pandas.DataFrame`.

- Started a draft for `JOSS <https://joss.theoj.org/>`_ (:commit:`1e37df44073261df2af486a2da014be8f59bc4cd`).

- New utilities function: :func:`argopy.utilities.open_etopo1`, :func:`argopy.show_versions`.

**Breaking changes with previous versions**

- The ``backend`` option in data fetchers and the global option ``datasrc`` have been renamed to ``src``. This makes the code more coherent (:commit:`ec6b32e94b78b2510985cfda49025c10ba97ecab`).

**Code management**

- Add Pypi automatic release publishing with github actions (:commit:`c4307885622709881e34909fd42e43f16a6a7cf4`)

- Remove Travis CI, fully adopt Github actions (:commit:`c4557425718f700b4aee760292b20b0642181dc6`)

- Improved unit testing (:commit:`e9555d1e6e90d3d1e75183cec0c4e14f7f19c17c`, :commit:`4b60ede844e37df86b32e4e2a2008335472a8cc1`, :commit:`34abf4913cb8bec027f88301c5504ebe594b3eae`)

v0.1.2 (15 May 2020)
--------------------

We didn't like this one this morning, so we move one to the next one !

v0.1.1 (3 Apr. 2020)
---------------------

**Features and front-end API**

- Added new data fetcher backend ``localftp`` in DataFetcher (:commit:`c5f7cb6f59d1f64a35dad28f386c9b1166883b81`):

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    argo_loader = ArgoDataFetcher(backend='localftp', path_ftp='/data/Argo/ftp_copy')
    argo_loader.float(6902746).to_xarray()

- Introduced global ``OPTIONS`` to set values for: cache folder, dataset (eg:`phy` or `bgc`), local ftp path, data fetcher (`erddap` or `localftp`) and user level (`standard` or `expert`). Can be used in context `with` (:commit:`83ccfb5110aa6abc6e972b92ba787a3e1228e33b`):

.. code-block:: python

    with argopy.set_options(mode='expert', datasrc='erddap'):
        ds = argopy.DataFetcher().float(3901530).to_xarray()

- Added a ``argopy.tutorial`` module to be able to load sample data for documentation and unit testing (:commit:`4af09b55a019a57fc3f1909a70e463f26f8863a1`):

.. code-block:: python

    ftproot, flist = argopy.tutorial.open_dataset('localftp')
    txtfile = argopy.tutorial.open_dataset('weekly_index_prof')

- Improved xarray *argo* accessor. Added methods for casting data types, to filter variables according to data mode, to filter variables according to quality flags. Useful methods to transform collection of points into collection of profiles, and vice versa (:commit:`14cda55f437f53cb19274324dce3e81f64bbb08f`):

.. code-block:: python

    ds = argopy.DataFetcher().float(3901530).to_xarray() # get a collection of points
    dsprof = ds.argo.point2profile() # transform to profiles
    ds = dsprof.argo.profile2point() # transform to points

- Changed License from MIT to Apache (:commit:`25f90c9cf6eab15c249c233c1677faaf5dc403c4`)

**Internal machinery**

- Add ``__all__`` to control ``from argopy import *`` (:commit:`83ccfb5110aa6abc6e972b92ba787a3e1228e33b`)

- All data fetchers inherit from class ``ArgoDataFetcherProto`` in ``proto.py`` (:commit:`44f45a5657f0ef7d06583df7142db61f82d1482e`)

- Data fetchers use default options from global OPTIONS

- In Erddap fetcher: methods to cast data type, to filter by data mode and by QC flags are now delegated to the xarray argo accessor methods.

- Data fetchers methods to filter variables according to user mode are using variable lists defined in utilities.

- ``argopy.utilities`` augmented with listing functions of: backends, standard variables and multiprofile files variables.

- Introduce custom errors in errors.py (:commit:`2563c9f0328121279a9b43220d197a622d1db12f`)

- Front-end API ArgoDataFetcher uses a more general way of auto-discovering fetcher backend and their access points. Turned of the ``deployments`` access point, waiting for the index fetcher to do that.

- Improved xarray *argo* accessor. More reliable ``point2profile`` and data type casting with ``cast_type``

**Code management**

- Add CI with github actions (:commit:`ecbf9bacded7747f27c698e90377e5ee40fc8999`)

- Contribution guideline for data fetchers (:commit:`b332495fce7f1650ae5bb8ec3148ade4c4f72702`)

- Improve unit testing (all along commits)

- Introduce code coverage (:commit:`b490ab56581d1ce0f58b44df532e35e87ecf04ff`)

- Added explicit support for python 3.6 , 3.7 and 3.8 (:commit:`58f60fe88a3aa85357754cafab8d89a4d948f35a`)


v0.1.0 (17 Mar. 2020)
---------------------

- Initial release.

- Erddap data fetcher
.. _data_viz:

Data visualisation
##################

Although **argopy** is not focus on visualisation, it provides a few functions to get you started. Plotting functions are available for both the data and index fetchers.

Trajectories
------------

.. code-block:: python

    from argopy import IndexFetcher as ArgoIndexFetcher
    idx = ArgoIndexFetcher().float([6902745, 6902746]).load()
    fig, ax = idx.plot('trajectory')
    fig, ax = idx.plot()  # Trajectory is the default plot

.. image:: _static/trajectory_sample.png

Some options are available to customise the plot, for instance:

.. code-block:: python

    from argopy import DataFetcher as ArgoDataFetcher
    idx = ArgoDataFetcher().float([6901020, 6902746, 2903359]).load()
    fig, ax = idx.plot('trajectory', style='white', palette='hls', figsize=(10,6), set_global=True)

.. image:: _static/trajectory_sample_white.png


Histograms on properties
------------------------

It is also possible to create bar plot for histograms on some data properties: 'profiler' and 'dac':

.. code-block:: python

    from argopy import IndexFetcher as ArgoIndexFetcher
    idx = ArgoIndexFetcher().region([-80,-30,20,50,'2021-01','2021-08']).load()
    fig, ax = idx.plot('dac')

.. image:: _static/bar_dac.png

.. code-block:: python

    fig, ax = idx.plot('profiler')

.. image:: _static/bar_profiler.png


Float dashboard
---------------

When working in Jupyter notebook, you can insert the EuroArgo dashboard in a cell with:

.. code-block:: python

    import argopy
    argopy.dashboard()

.. image:: _static/dashboard.png

and for a specific float, just provide its WMO:

.. code-block:: python

    import argopy
    argopy.dashboard(wmo=6902746)

.. image:: _static/dashboard_float.png
##############
Quick overview
##############

Here are some quick examples of what you can do with :py:mod:`argopy` objects.
Everything will be explained in much more details in the rest of the documentation.
.. _data_qc:

Data quality control
====================

.. contents::
   :local:

**argopy** comes with methods to help you quality control measurements. This section is probably intended for `expert` users.

Most of these methods are available through the :class:`xarray.Dataset` accessor namespace ``argo``. This means that if your dataset is `ds`, then you can use `ds.argo` to access more **argopy** functionalities.

Let's start with standard import:

.. ipython:: python
    :okwarning:

    from argopy import DataFetcher as ArgoDataFetcher


Salinity calibration
--------------------

.. currentmodule:: xarray

The Argo salinity calibration method is called OWC_, after the names of the core developers: Breck Owens, Anny Wong and Cecile Cabanes.
Historically, the OWC method has been implemented in `Matlab <https://github.com/ArgoDMQC/matlab_owc>`_ . More recently a `python version has been developed <https://github.com/euroargodev/argodmqc_owc>`_.

Preprocessing data
^^^^^^^^^^^^^^^^^^

At this point, both OWC software take as input a pre-processed version of the Argo float data to evaluate/calibrate.

**argopy** is able to perform this preprocessing and to create a *float source* data to be used by OWC software. This is made by :meth:`Dataset.argo.create_float_source`.

First, you would need to fetch the Argo float data you want to calibrate, in ``expert`` mode:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher(mode='expert').float(6902766).load().data

Then, to create the float source data, you call the method and provide a folder name to save output files:

.. ipython:: python
    :okwarning:

    ds.argo.create_float_source("float_source")

This will create the ``float_source/6902766.mat`` Matlab files to be set directly in the configuration file of the OWC software. This routine implements the same pre-processing as in the Matlab version (which is hosted on `this repo <https://github.com/euroargodev/dm_floats>`_ and ran with `this routine <https://github.com/euroargodev/dm_floats/blob/master/src/ow_source/create_float_source.m>`_). All the detailed steps of this pre-processing are given in the :meth:`Dataset.argo.create_float_source` API page.

.. note::
    If the dataset contains data from more than one float, several Matlab files are created, one for each float. This will allow you to prepare data from a collection of floats.

If you don't specify a path name, the method returns a dictionary with the float WMO as keys and pre-processed data as :class:`xarray.Dataset` as values.

.. ipython:: python
    :okwarning:

    ds_source = ds.argo.create_float_source()
    ds_source

See all options available for this method here: :meth:`Dataset.argo.create_float_source`.

The method partially relies on two others:

- :meth:`Dataset.argo.filter_scalib_pres`: to filter variables according to OWC salinity calibration software requirements. This filter modifies pressure, temperature and salinity related variables of the dataset.

- :meth:`Dataset.argo.groupby_pressure_bins`: to sub-sampled measurements by pressure bins. This is an excellent alternative to the :meth:`Dataset.argo.interp_std_levels` to avoid interpolation and preserve values of raw measurements while at the same time aligning measurements along approximately similar pressure levels (depending on the size of the bins). See more description at here: :ref:`Pressure levels: Group-by bins`.

Running the calibration
^^^^^^^^^^^^^^^^^^^^^^^

Please refer to the `OWC python software documentation <https://github.com/euroargodev/argodmqc_owc>`_.

A typical workflow would look like this:

.. code-block:: python

    import os, shutil
    from pathlib import Path

    import pyowc as owc
    import argopy
    from argopy import DataFetcher

    # Define float to calibrate:
    FLOAT_NAME = "6903010"

    # Set-up where to save OWC analysis results:
    results_folder = './analysis/%s' % FLOAT_NAME
    Path(results_folder).mkdir(parents=True, exist_ok=True)
    shutil.rmtree(results_folder) # Clean up folder content
    Path(os.path.sep.join([results_folder, 'float_source'])).mkdir(parents=True, exist_ok=True)
    Path(os.path.sep.join([results_folder, 'float_calib'])).mkdir(parents=True, exist_ok=True)
    Path(os.path.sep.join([results_folder, 'float_mapped'])).mkdir(parents=True, exist_ok=True)
    Path(os.path.sep.join([results_folder, 'float_plots'])).mkdir(parents=True, exist_ok=True)

    # fetch the default configuration and parameters
    USER_CONFIG = owc.configuration.load()

    # Fix paths to run at Ifremer:
    for k in USER_CONFIG:
        if "FLOAT" in k and "data/" in USER_CONFIG[k][0:5]:
            USER_CONFIG[k] = os.path.abspath(USER_CONFIG[k].replace("data", results_folder))
    USER_CONFIG['CONFIG_DIRECTORY'] = os.path.abspath('../data/constants')
    USER_CONFIG['HISTORICAL_DIRECTORY'] = os.path.abspath('/Volumes/OWC/CLIMATOLOGY/')  # where to find ARGO_for_DMQC_2020V03 and CTD_for_DMQC_2021V01 folders
    USER_CONFIG['HISTORICAL_ARGO_PREFIX'] = 'ARGO_for_DMQC_2020V03/argo_'
    USER_CONFIG['HISTORICAL_CTD_PREFIX'] = 'CTD_for_DMQC_2021V01/ctd_'
    print(owc.configuration.print_cfg(USER_CONFIG))

    # Create float source data with argopy:
    fetcher_for_real = DataFetcher(src='localftp', cache=True, mode='expert').float(FLOAT_NAME)
    fetcher_sample = DataFetcher(src='localftp', cache=True, mode='expert').profile(FLOAT_NAME, [1, 2])  # To reduce execution time for demo
    ds = fetcher_sample.load().data
    ds.argo.create_float_source(path=USER_CONFIG['FLOAT_SOURCE_DIRECTORY'], force='default')

    # Prepare data for calibration: map salinity on theta levels
    owc.calibration.update_salinity_mapping("", USER_CONFIG, FLOAT_NAME)

    # Set the calseries parameters for analysis and line fitting
    owc.configuration.set_calseries("", FLOAT_NAME, USER_CONFIG)

    # Calculate the fit of each break and calibrate salinities
    owc.calibration.calc_piecewisefit("", FLOAT_NAME, USER_CONFIG)

    # Results figures
    owc.plot.dashboard("", FLOAT_NAME, USER_CONFIG)

OWC references
^^^^^^^^^^^^^^

.. [OWC] See all the details about the OWC methodology in these references:

"An improved calibration method for the drift of the conductivity sensor on autonomous CTD profiling floats by θ–S climatology".
Deep-Sea Research Part I: Oceanographic Research Papers, 56(3), 450-457, 2009. https://doi.org/10.1016/j.dsr.2008.09.008

"Improvement of bias detection in Argo float conductivity sensors and its application in the North Atlantic".
Deep-Sea Research Part I: Oceanographic Research Papers, 114, 128-136, 2016. https://doi.org/10.1016/j.dsr.2016.05.007


Trajectories
------------

Topography
^^^^^^^^^^
.. currentmodule:: argopy

For some QC of trajectories, it can be useful to easily get access to the topography. This can be done with the **argopy** utility :class:`TopoFetcher`:

.. ipython:: python
    :okwarning:
    
    from argopy import TopoFetcher
    box = [-65, -55, 10, 20]
    ds = TopoFetcher(box, cache=True).to_xarray()

.. image:: _static/topography_sample.png


Combined with the fetcher property ``domain``, it now becomes easy to superimpose float trajectory with topography:

.. ipython:: python
    :okwarning:

    fetcher = ArgoDataFetcher().float(2901623)
    ds = TopoFetcher(fetcher.domain[0:4], cache=True).to_xarray()

.. code-block:: python

    fig, ax = loader.plot('trajectory', figsize=(10, 10))
    ds['elevation'].plot.contourf(levels=np.arange(-6000,0,100), ax=ax, add_colorbar=False)

.. image:: _static/trajectory_topography_sample.png


.. note::
    The :class:`TopoFetcher` can return a lower resolution topography with the ``stride`` option. See the :class:`argopy.TopoFetcher` full documentation for all the details.


Altimetry
---------
.. currentmodule:: argopy

Satellite altimeter measurements can be used to check the quality of the Argo profiling floats time series. The method compares collocated sea level anomalies from altimeter measurements and dynamic height anomalies calculated from Argo temperature and salinity profiles for each Argo float time series [Guinehut2008]_. This method is performed routinely by CLS and results are made available online.


**argopy** provides a simple access to this QC analysis with an option to the data and index fetchers :meth:`DataFetcher.plot` methods that will insert the CLS Satellite Altimeter report figure on a notebook cell.

.. code-block:: python

    fetcher = ArgoDataFetcher().float(6902745)
    fetcher.plot('qc_altimetry', embed='list')

.. image:: https://data-argo.ifremer.fr/etc/argo-ast9-item13-AltimeterComparison/figures/6902745.png

See all details about this method here: :meth:`argopy.plotters.open_sat_altim_report`


.. rubric:: References

.. [Guinehut2008] Guinehut, S., Coatanoan, C., Dhomps, A., Le Traon, P., & Larnicol, G. (2009). On the Use of Satellite Altimeter Data in Argo Quality Control, Journal of Atmospheric and Oceanic Technology, 26(2), 395-402. `10.1175/2008JTECHO648.1 <https://doi.org/10.1175/2008JTECHO648.1>`_
.. _why:

Why argopy ?
============

Surprisingly, the Argo community never provided its user base with a Python software to easily access and manipulate Argo measurements:
**argopy** aims to fill this gap.

Despite, or because, its tremendous success in data management and in developping good practices and well calibrated procedures [ADMT]_, the Argo dataset is very complex: with thousands of different variables, tens of reference tables and a `user manual <http://dx.doi.org/10.13155/29825>`_ more than 100 pages long:
**argopy** aims to help you navigate this complex realm.

For non-experts of the Argo dataset, it has become rather complicated to get access to Argo measurements.
This is mainly due to:

* Argo measurements coming from many different models of floats or sensors,
* quality control of *in situ* measurements of autonomous platforms being really a matter of ocean and data experts,
* the Argo data management workflow being distributed between more than 10 Data Assembly Centers all around the world.

Less data wrangling, more scientific analysis
---------------------------------------------

In order to ease Argo data analysis for the vast majority of **standard** users, we implemented in **argopy** different levels of verbosity and data processing to hide or simply remove variables only meaningful to **experts**.
Let **argopy** manage data wrangling, and focus on your scientific analysis.

If you don't know in which category you would place yourself, try to answer the following questions:

* [ ] what is a WMO number ?
* [ ] what is the difference between Delayed and Real Time data mode ?
* [ ] what is an adjusted parameter ?
* [ ] what a QC flag of 3 means ?

If you don't answer to more than 1 question: you probably will feel more comfortable with the *standard* user mode.

By default, all **argopy** data fetchers are set to work with a **standard** user mode, the other possible mode is **expert**.

In *standard* mode, fetched data are automatically filtered to account for their quality (only good are retained) and level of processing by the data centers (whether they looked at the data briefly or not).

Selecting user mode is further explained in the dedicated documentation section: :ref:`user-mode`.

.. [ADMT] See all the ADMT documentation here: http://www.argodatamgt.org/Documentation.. _what_is_argo:

What is Argo ?
##############

**Argo is a real-time global ocean in situ observing system**.

The ocean is a key component of the Earth climate system. It thus needs a continuous real-time monitoring to help scientists
better understand its dynamic and predict its evolution. All around the world, oceanographers have managed to join their
efforts and set up a `Global Ocean Observing System <https://www.goosocean.org>`_ among which *Argo* is a key component.

*Argo* is a global network of nearly 4000 autonomous probes measuring
pressure, temperature and salinity from the surface to 2000m depth every 10 days. The localisation of these probes is
nearly random between the 60th parallels (`see live coverage here <http://map.argo-france.fr>`_).
All probes data are collected by satellite in real-time, processed by several data centers and finally merged in a single
dataset (collecting more than 2 millions of vertical profiles data) made freely available to anyone through
a `ftp server <ftp://ftp.ifremer.fr/ifremer/argo>`_ or `monthly zip snapshots <http://dx.doi.org/10.17882/42182>`_.

The Argo international observation array was initiated in 1999 and soon revolutionized our
perspective on the large scale structure and variability of the ocean by providing seasonally and regionally unbiased
in situ temperature/salinity measurements of the ocean interior, key information that satellites can't provide
(`Riser et al, 2016 <http://dx.doi.org/10.1038/nclimate2872>`_).

The Argo array reached its full global coverage (of 1 profile per month and per 3x3 degree horizontal area) in 2007, and
continuously pursues its evolution to fulfill new scientific requirements (`Roemmich et al, 2019
<https://www.frontiersin.org/article/10.3389/fmars.2019.00439>`_). It now extents to higher latitudes and some of the
floats are able to profile down to 4000m and 6000m. New floats are also equipped with biogeochemical sensors, measuring
oxygen and chlorophyll for instance. Argo is thus providing a deluge of in situ data: more than 400 profiles per day.

Each Argo probe is an autonomous, free drifting, profiling float, i.e. a probe that can't control its trajectory but
is able to control its buoyancy and thus to move up and down the water column as it wishes. Argo floats continuously
operate the same program, or cycle, illustrated in the figure below. After 9 to 10 days of free drift at a parking
depth of about 1000m, a typical Argo float dives down to 2000m and then shoals back to the surface while measuring pressure,
temperature and salinity. Once it reaches the surface, the float sends by satellite its measurements to a data center
where they are processed in real time and made freely available on the web in less than 24h00.

*Typical 10 days program, cycle, of an Argo float*:

.. image:: _static/argofloats_cycle.png

Usage
=====

To get access to Argo data, all you need is 2 lines of codes:

.. ipython:: python
    :okwarning:

    from argopy import DataFetcher as ArgoDataFetcher
    ds = ArgoDataFetcher().region([-75, -45, 20, 30, 0, 100, '2011-01', '2011-06']).to_xarray()

In this example, we used a data fetcher to get data for a given space/time region.
We retrieved all Argo data measurements from 75W to 45W, 20N to 30N, 0db to 100db and from January to May 2011 (the max date is exclusive).
Data are returned as a collection of measurements in a :class:`xarray.Dataset`:

.. ipython:: python
    :okwarning:

    ds

.. currentmodule:: xarray

Fetched data are returned as a 1D array collection of measurements. If you prefer to work with a 2D array collection of vertical profiles, simply transform the dataset with the :class:`xarray.Dataset` accessor method :meth:`Dataset.argo.point2profile`:

.. ipython:: python
    :okwarning:

    ds = ds.argo.point2profile()
    ds

You can also fetch data for a specific float using its `WMO number <https://www.wmo.int/pages/prog/amp/mmop/wmo-number-rules.html>`_:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().float(6902746).to_xarray()

or for a float profile using the cycle number:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().profile(6902755, 12).to_xarray()
.. _metadata_fetching:

Fetching Argo meta-data
=======================

Since the Argo measurements dataset is quite complex, it comes with a collection of index files, or lookup tables with meta data. These index help you determine what you can expect before retrieving the full set of measurements. **argopy** has a specific fetcher for index:

.. ipython:: python
    :okwarning:

    from argopy import IndexFetcher as ArgoIndexFetcher

You can use the Index fetcher with the ``region`` or ``float`` access points, similarly to data fetching:

.. ipython:: python
    :suppress:

    import argopy
    ftproot = argopy.tutorial.open_dataset('localftp')[0]
    argopy.set_options(local_ftp=ftproot)

.. ipython:: python
    :okwarning:

    idx = ArgoIndexFetcher(src='localftp').float(2901623).load()
    idx.index

Alternatively, you can use :meth:`argopy.IndexFetcher.to_dataframe()`.

See :ref:`Fetching methods` for a list of all methods available for the Index fetcher... Generate API reference pages, but don't display these in tables.
.. This extra page is a work around for sphinx not having any support for
.. hiding an autosummary table.

.. autosummary::
    :toctree: generated/

    argopy

    argopy.fetchers
    argopy.fetchers.ArgoDataFetcher
    argopy.fetchers.ArgoDataFetcher.region
    argopy.fetchers.ArgoDataFetcher.float
    argopy.fetchers.ArgoDataFetcher.profile
    argopy.fetchers.ArgoDataFetcher.load
    argopy.fetchers.ArgoDataFetcher.to_xarray
    argopy.fetchers.ArgoDataFetcher.to_dataframe
    argopy.fetchers.ArgoDataFetcher.to_index
    argopy.fetchers.ArgoDataFetcher.plot
    argopy.fetchers.ArgoDataFetcher.uri
    argopy.fetchers.ArgoDataFetcher.data
    argopy.fetchers.ArgoDataFetcher.index
    argopy.fetchers.ArgoDataFetcher.dashboard
    argopy.fetchers.ArgoDataFetcher.clear_cache

    argopy.fetchers.ArgoIndexFetcher
    argopy.fetchers.ArgoIndexFetcher.region
    argopy.fetchers.ArgoIndexFetcher.float
    argopy.fetchers.ArgoIndexFetcher.profile
    argopy.fetchers.ArgoIndexFetcher.load
    argopy.fetchers.ArgoIndexFetcher.to_xarray
    argopy.fetchers.ArgoIndexFetcher.to_dataframe
    argopy.fetchers.ArgoIndexFetcher.to_csv
    argopy.fetchers.ArgoIndexFetcher.plot
    argopy.fetchers.ArgoIndexFetcher.index
    argopy.fetchers.ArgoIndexFetcher.clear_cache

    argopy.data_fetchers.erddap_data.ErddapArgoDataFetcher
    argopy.data_fetchers.erddap_data.Fetch_wmo
    argopy.data_fetchers.erddap_data.Fetch_box

    argopy.data_fetchers.localftp_data.LocalFTPArgoDataFetcher
    argopy.data_fetchers.localftp_data.Fetch_wmo
    argopy.data_fetchers.localftp_data.Fetch_box

    argopy.data_fetchers.argovis_data.ArgovisDataFetcher
    argopy.data_fetchers.argovis_data.Fetch_wmo
    argopy.data_fetchers.argovis_data.Fetch_box

    argopy.options.set_options

    argopy.tutorial.open_dataset

    argopy.utilities.monitor_status
    argopy.utilities.show_versions
    argopy.utilities.show_options
    argopy.utilities.clear_cache
    argopy.utilities.list_available_data_src
    argopy.utilities.list_available_index_src
    argopy.utilities.Chunker
    
    argopy.utilities.groupby_remap
    argopy.utilities.linear_interpolation_remap

    argopy.utilities.TopoFetcher.cname
    argopy.utilities.TopoFetcher.define_constraints
    argopy.utilities.TopoFetcher.get_url
    argopy.utilities.TopoFetcher.load
    argopy.utilities.TopoFetcher.to_xarray
    argopy.utilities.TopoFetcher.cachepath
    argopy.utilities.TopoFetcher.uri

    argopy.utilities.list_standard_variables
    argopy.utilities.list_multiprofile_file_variables
    argopy.utilities.check_localftp
    argopy.utilities.format_oneline
    argopy.utilities.is_box
    argopy.utilities.is_indexbox
    argopy.utilities.is_wmo
    argopy.utilities.check_wmo
    argopy.utilities.wmo2box

    argopy.plotters.open_dashboard
    argopy.plotters.bar_plot
    argopy.plotters.plot_trajectory
    argopy.plotters.open_sat_altim_report

    argopy.stores.filesystems.filestore
    argopy.stores.filestore.open_dataset
    argopy.stores.filestore.read_csv

    argopy.stores.filestore.open
    argopy.stores.filestore.glob
    argopy.stores.filestore.exists
    argopy.stores.filestore.store_path
    argopy.stores.filestore.register
    argopy.stores.filestore.cachepath
    argopy.stores.filestore.clear_cache
    argopy.stores.filestore.open_mfdataset

    argopy.stores.filesystems.httpstore
    argopy.stores.httpstore.open_json
    argopy.stores.httpstore.open_dataset
    argopy.stores.httpstore.read_csv
    argopy.stores.httpstore.open
    argopy.stores.httpstore.glob
    argopy.stores.httpstore.exists
    argopy.stores.httpstore.store_path
    argopy.stores.httpstore.register
    argopy.stores.httpstore.cachepath
    argopy.stores.httpstore.clear_cache
    argopy.stores.httpstore.open_mfdataset
    argopy.stores.httpstore.open_mfjson

    argopy.stores.filesystems.memorystore
    argopy.stores.memorystore.open
    argopy.stores.memorystore.glob
    argopy.stores.memorystore.exists
    argopy.stores.memorystore.store_path
    argopy.stores.memorystore.register
    argopy.stores.memorystore.cachepath
    argopy.stores.memorystore.clear_cache
    argopy.stores.memorystore.open_dataset
    argopy.stores.memorystore.open_mfdataset
    argopy.stores.memorystore.read_csv

    argopy.stores.argo_index.indexstore
    argopy.stores.argo_index.indexfilter_wmo
    argopy.stores.argo_index.indexfilter_box
    
    argopy.xarray.ArgoAccessor.point2profile
    argopy.xarray.ArgoAccessor.profile2point
    argopy.xarray.ArgoAccessor.interp_std_levels
    argopy.xarray.ArgoAccessor.groupby_pressure_bins
    argopy.xarray.ArgoAccessor.teos10
    argopy.xarray.ArgoAccessor.create_float_source
    argopy.xarray.ArgoAccessor.filter_qc
    argopy.xarray.ArgoAccessor.filter_data_mode
    argopy.xarray.ArgoAccessor.filter_scalib_pres
    argopy.xarray.ArgoAccessor.cast_types
.. _data_fetching:

Fetching Argo data
==================

To access Argo data, you need to use a data fetcher. You can import and instantiate the default argopy data fetcher
like this:

.. ipython:: python
    :okwarning:

    from argopy import DataFetcher as ArgoDataFetcher
    argo_loader = ArgoDataFetcher()
    argo_loader

Then, you can request data for a specific **space/time domain**, for a given **float** or for a given vertical **profile**.

If you fetch a lot of data, you may want to look at the :ref:`performances` section.

For a space/time domain
-----------------------

Use the fetcher access point :meth:`argopy.DataFetcher.region` to specify a domain and chain with the :meth:`argopy.DataFetcher.to_xarray` to get the data returned as :class:`xarray.Dataset`.

For instance, to retrieve data from 75W to 45W, 20N to 30N, 0db to 10db and from January to May 2011:

.. ipython:: python
    :okwarning:

    ds = argo_loader.region([-75, -45, 20, 30, 0, 10, '2011-01-01', '2011-06']).to_xarray()
    ds

Note that:

- the constraints on time is not mandatory: if not specified, the fetcher will return all data available in this region.

- the last time bound is exclusive: that's why here we specify June to retrieve data collected in May.

For one or more floats
----------------------

If you know the Argo float unique identifier number called a `WMO number <https://www.wmo.int/pages/prog/amp/mmop/wmo-number-rules.html>`_ you can use the fetcher access point :meth:`argopy.DataFetcher.float` to specify the float WMO platform number and chain with the :meth:`argopy.DataFetcher.to_xarray` to get the data returned as :class:`xarray.Dataset`.

For instance, to retrieve data for float WMO *6902746*:

.. ipython:: python
    :okwarning:

    ds = argo_loader.float(6902746).to_xarray()
    ds

To fetch data for a collection of floats, input them in a list:

.. ipython:: python
    :okwarning:

    ds = argo_loader.float([6902746, 6902755]).to_xarray()
    ds

For one or more profiles
------------------------

Use the fetcher access point :meth:`argopy.DataFetcher.profile` to specify the float WMO platform number and the profile cycle number to retrieve profiles for, then chain with the :meth:`argopy.DataFetcher.to_xarray` to get the data returned as :class:`xarray.Dataset`.

For instance, to retrieve data for the 12th profile of float WMO 6902755:

.. ipython:: python
    :okwarning:

    ds = argo_loader.profile(6902755, 12).to_xarray()
    ds

To fetch data for more than one profile, input them in a list:

.. ipython:: python
    :okwarning:

    ds = argo_loader.profile(6902755, [3, 12]).to_xarray()
    ds
.. _starting:

Getting started with <img src="https://raw.githubusercontent.com/euroargodev/argopy/master/docs/_static/argopy_logo_long.png" alt="argopy logo" width="200"/>

Import the **argopy** data fetcher:

.. ipython:: python
    :okwarning:

    from argopy import DataFetcher as ArgoDataFetcher

Then, to get access to Argo data, all you need is 1 line of code:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().region([-75, -45, 20, 30, 0, 100, '2011', '2012']).to_xarray()

In this example, we used a data fetcher to get data for a given space/time region.
We retrieved all Argo data measurements from 75W to 45W, 20N to 30N, 0db to 100db and from January to May 2011 (the max date is exclusive).
Data are returned as a collection of measurements in a `xarray.Dataset <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html>`.

.. ipython:: python
    :okwarning:

    ds

Fetched data are returned as a 1D array collection of measurements.

If you prefer to work with a 2D array collection of vertical profiles, simply transform the dataset with the `xarray.Dataset <http://xarray.pydata.org/en/stable/generated/xarray.Dataset.html>` accessor method `argo.point2profile <https://argopy.readthedocs.io/en/latest/api.html#argopy.ArgoAccessor.point2profile>`:

.. ipython:: python
    :okwarning:

    ds = ds.argo.point2profile()
    ds

You can also fetch data for a specific float using its [WMO number](<https://www.wmo.int/pages/prog/amp/mmop/wmo-number-rules.html):

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().float(6902746).to_xarray()

or for a float profile using the cycle number:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().profile(6902755, 12).to_xarray()Argo data python library
========================

**argopy** is a python library that aims to ease :ref:`Argo <what_is_argo>` data access, manipulation and visualisation
for standard users as well as Argo experts.

Documentation
-------------

**Getting Started**

* :doc:`install`
* :doc:`usage`
* :doc:`why`
* :doc:`what_is_argo`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    install
    usage
    why
    what_is_argo

**User Guide**

* :doc:`data_fetching`
* :doc:`data_sources`
* :doc:`data_manipulation`
* :doc:`visualisation`
* :doc:`user_mode`
* :doc:`metadata_fetching`
* :doc:`performances`
* :doc:`data_quality_control`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: User Guide

    data_fetching
    data_sources
    data_manipulation
    visualisation
    data_quality_control
    user_mode
    metadata_fetching
    performances

**Help & reference**

* :doc:`whats-new`
* :doc:`contributing`
* :doc:`api`

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Help & reference

    whats-new
    contributing
    api
Data sources
============

.. contents::
   :local:

Let's start with standard import:

.. ipython:: python
    :okwarning:

    import argopy
    from argopy import DataFetcher as ArgoDataFetcher

Selecting a source
------------------

**argopy** can get access to Argo data from different sources:

1. the `Ifremer erddap server <http://www.ifremer.fr/erddap>`__.

   | The erddap server database is updated daily and doesn’t require you
     to download anymore data than what you need.
   | You can select this data source with the keyword ``erddap`` and
     methods described below. The Ifremer erddap dataset is based on
     mono-profile files of the GDAC.

2. your local collection of Argo files, organised as in the `GDAC
   ftp <http://www.argodatamgt.org/Access-to-data/Argo-GDAC-ftp-and-https-servers>`__.

   | This is how you would use **argopy** with your data, as long as
     they are formatted and organised the Argo way.
   | You can select this data source with the keyword ``localftp`` and
     methods described below.

3. the `Argovis server <https://argovis.colorado.edu/>`__.

   The Argovis server database is updated daily and provides access to
   curated Argo data (QC=1 only). You can select this data source with
   the keyword ``argovis`` and methods described below.

You have several ways to specify which data source you want to use:

-  **using argopy global options**:

.. ipython:: python
    :okwarning:

    argopy.set_options(src='erddap')

-  **in a temporary context**:

.. ipython:: python
    :okwarning:

    with argopy.set_options(src='erddap'):
        loader = ArgoDataFetcher().profile(6902746, 34)

-  **with an argument in the data fetcher**:

.. ipython:: python
    :okwarning:

    loader = ArgoDataFetcher(src='erddap').profile(6902746, 34)

Setting a local copy of the GDAC ftp
------------------------------------

Data fetching with the ``localftp`` data source will require you to
specify the path toward your local copy of the GDAC ftp server with the
``local_ftp`` option.

This is not an issue for expert users, but standard users may wonder how
to set this up. The primary distribution point for Argo data, the only
one with full support from data centers and with nearly a 100% time
availability, is the GDAC ftp. Two mirror servers are available:

-  France Coriolis: ftp://ftp.ifremer.fr/ifremer/argo
-  US GODAE: ftp://usgodae.org/pub/outgoing/argo

If you want to get your own copy of the ftp server content, Ifremer
provides a nice rsync service. The rsync server “vdmzrs.ifremer.fr”
provides a synchronization service between the “dac” directory of the
GDAC and a user mirror. The “dac” index files are also available from
“argo-index”.

From the user side, the rsync service:

-  Downloads the new files
-  Downloads the updated files
-  Removes the files that have been removed from the GDAC
-  Compresses/uncompresses the files during the transfer
-  Preserves the files creation/update dates
-  Lists all the files that have been transferred (easy to use for a
   user side post-processing)

To synchronize the whole dac directory of the Argo GDAC:

.. code:: bash

   rsync -avzh --delete vdmzrs.ifremer.fr::argo/ /home/mydirectory/...

To synchronize the index:

.. code:: bash

   rsync -avzh --delete vdmzrs.ifremer.fr::argo-index/ /home/mydirectory/...

.. note::

    The first synchronisation of the whole dac directory of the Argo GDAC (365Gb) can take quite a long time (several hours).

Comparing data sources
----------------------

Features
~~~~~~~~

Each of the available data sources have their own features and
capabilities. Here is a summary:

======================= ====== ======== =======
Data source:            erddap localftp argovis
======================= ====== ======== =======
**Access Points**                       
region                  X      X        X
float                   X      X        X
profile                 X      X        X
**User mode**                           
standard                X      X        X
expert                  X      X        
**Dataset**                             
core (T/S)              X      X        X
BGC                                     
Reference data for DMQC X               
**Parallel method**                     
multi-threading         X      X        X
multi-processes                X        
Dask client                             
======================= ====== ======== =======

Fetched data and variables
~~~~~~~~~~~~~~~~~~~~~~~~~~

| You may wonder if the fetched data are different from the available
  data sources.
| This will depend on the last update of each data sources and of your
  local data.

Let's retrieve one float data from a local sample of the GDAC ftp (a sample GDAC ftp is downloaded automatically with the method :meth:`argopy.tutorial.open_dataset`):

.. ipython:: python
    :okwarning:

    # Download ftp sample and get the ftp local path:
    ftproot = argopy.tutorial.open_dataset('localftp')[0]
    
    # then fetch data:
    with argopy.set_options(src='localftp', local_ftp=ftproot):
        ds = ArgoDataFetcher().float(1900857).to_xarray()
        print(ds)

Let’s now retrieve the latest data for this float from the ``erddap`` and ``argovis`` sources:

.. ipython:: python
    :okwarning:

    with argopy.set_options(src='erddap'):
        ds = ArgoDataFetcher().float(1900857).to_xarray()
        print(ds)

.. ipython:: python
    :okwarning:

    with argopy.set_options(src='argovis'):
        ds = ArgoDataFetcher().float(1900857).to_xarray()
        print(ds)

We can see some minor differences between ``localftp``/``erddap`` vs the
``argovis`` response: this later data source does not include the
descending part of the first profile, this explains why ``argovis``
returns slightly less data.

.. _api-status:

Status of sources
-----------------

With remote, online data sources, it may happens that the data server is experiencing down time. 
With local data sources, the availability of the path is checked when it is set. But it may happens that the path points to a disk that get unmounted or unplugged after the option setting.

If you're running your analysis on a Jupyter notebook, you can use the :meth:`argopy.status` method to insert a data status monitor on a cell output. All available data sources will be monitored continuously.

.. code-block:: python

    argopy.status()

.. image:: _static/status_monitor.png
  :width: 350
  
If one of the data source become unavailable, you will see the status bar changing to something like:
  
.. image:: _static/status_monitor_down.png
  :width: 350  
  
Note that the :meth:`argopy.status` method has a ``refresh`` option to let you specify the refresh rate in seconds of the monitoring.

Last, you can check out `the following argopy status webpage that monitors all important resources to the software <https://argopy.statuspage.io>`_.
Manipulating data
=================

.. contents::
   :local:

.. currentmodule:: xarray

Once you fetched data, **argopy** comes with a handy :class:`xarray.Dataset` accessor ``argo`` to perform specific manipulation of the data. This means that if your dataset is named `ds`, then you can use `ds.argo` to access more **argopy** functions. The full list is available in the API documentation page :ref:`Dataset.argo (xarray accessor)`.

Let's start with standard import:

.. ipython:: python
    :okwarning:

    from argopy import DataFetcher as ArgoDataFetcher

Transformation
--------------

Points vs profiles
^^^^^^^^^^^^^^^^^^

By default, fetched data are returned as a 1D array collection of measurements:

.. ipython:: python
    :okwarning:

    argo_loader = ArgoDataFetcher().region([-75,-55,30.,40.,0,100., '2011-01-01', '2011-01-15'])
    ds_points = argo_loader.to_xarray()
    ds_points

If you prefer to work with a 2D array collection of vertical profiles, simply transform the dataset with :meth:`Dataset.argo.point2profile`:

.. ipython:: python
    :okwarning:

    ds_profiles = ds_points.argo.point2profile()
    ds_profiles

You can simply reverse this transformation with the :meth:`Dataset.argo.profile2point`:

.. ipython:: python
    :okwarning:

    ds = ds_profiles.argo.profile2point()
    ds

Pressure levels: Interpolation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once your dataset is a collection of vertical **profiles**, you can interpolate variables on standard pressure levels using :meth:`Dataset.argo.interp_std_levels` with your levels as input:

.. ipython:: python
    :okwarning:

    ds_interp = ds_profiles.argo.interp_std_levels([0,10,20,30,40,50])
    ds_interp

Note on the linear interpolation process : 
    - Only profiles that have a maximum pressure higher than the highest standard level are selected for interpolation.
    - Remaining profiles must have at least five data points to allow interpolation.
    - For each profile, shallowest data point is repeated to the surface to allow a 0 standard level while avoiding extrapolation.

Pressure levels: Group-by bins
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you prefer to avoid interpolation, you can opt for a pressure bins grouping reduction using :meth:`Dataset.argo.groupby_pressure_bins`. This method can be used to subsample and align an irregular dataset (pressure not being similar in all profiles) on a set of pressure bins. The output dataset could then be used to perform statistics along the N_PROF dimension because N_LEVELS will corresponds to similar pressure bins.

To illustrate this method, let's start by fetching some data from a low vertical resolution float:

.. ipython:: python
    :okwarning:

    loader = ArgoDataFetcher(src='erddap', mode='expert').float(2901623)  # Low res float
    ds = loader.load().data

Let's now sub-sample these measurements along 250db bins, selecting values from the **deepest** pressure levels for each bins:

.. ipython:: python
    :okwarning:

    bins = np.arange(0.0, np.max(ds["PRES"]), 250.0)
    ds_binned = ds.argo.groupby_pressure_bins(bins=bins, select='deep')
    ds_binned

See the new ``STD_PRES_BINS`` variable that hold the pressure bins definition.

The figure below shows the sub-sampling effect:

.. code-block:: python

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import cmocean

    fig, ax = plt.subplots(figsize=(18,6))
    ds.plot.scatter(x='CYCLE_NUMBER', y='PRES', hue='PSAL', ax=ax, cmap=cmocean.cm.haline)
    plt.plot(ds_binned['CYCLE_NUMBER'], ds_binned['PRES'], 'r+')
    plt.hlines(bins, ds['CYCLE_NUMBER'].min(), ds['CYCLE_NUMBER'].max(), color='k')
    plt.hlines(ds_binned['STD_PRES_BINS'], ds_binned['CYCLE_NUMBER'].min(), ds_binned['CYCLE_NUMBER'].max(), color='r')
    plt.title(ds.attrs['Fetched_constraints'])
    plt.gca().invert_yaxis()

.. image:: _static/groupby_pressure_bins_select_deep.png

The bin limits are shown with horizontal red lines, the original data are in the background colored scatter and the group-by pressure bins values are highlighted in red marks

The ``select`` option can take many different values, see the full documentation of :meth:`Dataset.argo.groupby_pressure_bins` , for all the details. Let's show here results from the ``random`` sampling:

.. code-block:: python

    ds_binned = ds.argo.groupby_pressure_bins(bins=bins, select='random')

.. image:: _static/groupby_pressure_bins_select_random.png


Filters
^^^^^^^

If you fetched data with the ``expert`` mode, you may want to use *filters* to help you curate the data.

- **QC flag filter**: :meth:`Dataset.argo.filter_qc`. This method allows you to filter measurements according to QC flag values. This filter modifies all variables of the dataset.
- **Data mode filter**: :meth:`Dataset.argo.filter_data_mode`. This method allows you to filter variables according to their data mode. This filter modifies the <PARAM> and <PARAM_QC> variables of the dataset.
- **OWC variables filter**: :meth:`Dataset.argo.filter_scalib_pres`. This method allows you to filter variables according to OWC salinity calibration software requirements. This filter modifies pressure, temperature and salinity related variables of the dataset.


Complementary data
------------------

TEOS-10 variables
~~~~~~~~~~~~~~~~~

You can compute additional ocean variables from `TEOS-10 <http://teos-10.org/>`_. The default list of variables is: 'SA', 'CT', 'SIG0', 'N2', 'PV', 'PTEMP' ('SOUND_SPEED', 'CNDC' are optional). `Simply raise an issue to add a new one <https://github.com/euroargodev/argopy/issues/new/choose>`_.

This can be done using the :meth:`Dataset.argo.teos10` method and indicating the list of variables you want to compute:

.. ipython:: python
    :okwarning:

    ds = ArgoDataFetcher().float(2901623).to_xarray()
    ds.argo.teos10(['SA', 'CT', 'PV'])

.. ipython:: python
    :okwarning:

    ds['SA']

Data models
-----------

By default **argopy** works with :class:`xarray.Dataset` for Argo data fetcher, and with :class:`pandas.DataFrame` for Argo index fetcher.

For your own analysis, you may prefer to switch from one to the other. This is all built in **argopy**, with the :meth:`argopy.DataFetcher.to_dataframe` and :meth:`argopy.IndexFetcher.to_xarray` methods.

.. ipython:: python
    :okwarning:

    ArgoDataFetcher().profile(6902746, 34).to_dataframe()


Saving data
===========

Once you have your Argo data as :class:`xarray.Dataset`, simply use the awesome possibilities of `xarray <http://xarray.pydata.org>`_ like :meth:`xarray.Dataset.to_netcdf` or :meth:`xarray.Dataset.to_zarr`.
Installation
============

Required dependencies
^^^^^^^^^^^^^^^^^^^^^

- xarray
- scipy
- scikit-learn
- netCDF4
- dask
- toolz
- erddapy
- fsspec
- gsw
- aiohttp
- packaging


Note that Erddapy_ is required because `erddap <https://coastwatch.pfeg.noaa.gov/erddap/information.html>`_ is the default data fetching backend.

Requirement dependencies details can be found `here <https://github.com/euroargodev/argopy/network/dependencies#requirements.txt>`_.

The **argopy** software is `continuously tested <https://github.com/euroargodev/argopy/actions?query=workflow%3Atests>`_ with under latest OS (Linux and Mac OS) and with python versions 3.6, 3.7 and 3.8.

Optional dependencies
^^^^^^^^^^^^^^^^^^^^^

For a complete **argopy** experience, the following packages are also required:

- ipython>=5.0.0
- ipywidgets>=7.5.1
- tqdm>=4.46.0
- Matplotlib>=3.0
- Cartopy>=0.17
- Seaborn>=0.9.0

Instructions
^^^^^^^^^^^^

Install the last release with conda:

.. code-block:: text

    conda install -c conda-forge argopy

or pip:

.. code-block:: text

    pip install argopy

you can also work with the latest version:

.. code-block:: text

    pip install git+http://github.com/euroargodev/argopy.git@master

.. _Erddapy: https://github.com/ioos/erddapy

.. _user-mode:

User mode: standard vs expert
=============================

**Problem**

For beginners or non-experts of the Argo dataset, it can be quite
complicated to get access to Argo measurements. Indeed, the Argo data
set is very complex, with thousands of different variables, tens of
reference tables and a `user manual <https://doi.org/10.13155/29825>`__
more than 100 pages long.

This is mainly due to:

-  Argo measurements coming from many different models of floats or
   sensors,
-  quality control of *in situ* measurements of autonomous platforms
   being really a matter of ocean and data experts,
-  the Argo data management workflow being distributed between more than
   10 Data Assembly Centers all around the world,
-  the Argo autonomous profiling floats, despite quite a simple
   principle of functioning, is a rather complex robot that needs a lot
   of data to be monitored and logged.

**Solution**

In order to ease Argo data analysis for the vast majority of standard
users, we implemented in **argopy** different levels of verbosity and
data processing to hide or simply remove variables only meaningful to
experts.

What type of user are you ?
---------------------------

If you don’t know in which user category you would place yourself, try
to answer the following questions:

-  what is a WMO number ?
-  what is the difference between Delayed and Real Time data mode ?
-  what is an adjusted parameter ?
-  what a QC flag of 3 means ?

If you answered to no more than 1 question, you probably would feel more
comfortable with the **standard** user mode. Otherwise, you can give a
try to the **expert** mode.

In **standard** mode, fetched data are automatically filtered to account
for their quality (only good are retained) and level of processing by
the data centers (whether they looked at the data briefly or not).

Setting the user mode
---------------------

Let's start with standard import:

.. ipython:: python
    :okwarning:

    import argopy
    from argopy import DataFetcher as ArgoDataFetcher


By default, all **argopy** data fetchers are set to work with a
**standard** user mode.

If you want to change the user mode, or simply makes it explicit, you
can use:

-  **argopy** global options:

.. ipython:: python
    :okwarning:

    argopy.set_options(mode='standard')

-  a temporary context:

.. ipython:: python
    :okwarning:

    with argopy.set_options(mode='standard'):
        ArgoDataFetcher().profile(6902746, 34)

-  option when instantiating the data fetcher:

.. ipython:: python
    :okwarning:

    ArgoDataFetcher(mode='standard').profile(6902746, 34)

Differences in user modes
-------------------------

To highlight that, let’s compare data fetched for one profile with each
modes.

You will note that the **standard** mode has fewer variables to let you
focus on your analysis. For **expert**, all Argo variables for you to
work with are here.

The difference is the most visible when fetching Argo data from a local
copy of the GDAC ftp, so let’s use a sample of this provided by
**argopy** tutorial datasets:

.. ipython:: python
    :okwarning:

    ftproot, flist = argopy.tutorial.open_dataset('localftp')
    argopy.set_options(local_ftp=ftproot)

In **standard** mode:

.. ipython:: python
    :okwarning:

    with argopy.set_options(mode='standard'):
        ds = ArgoDataFetcher(src='localftp').profile(6901929, 2).to_xarray()
        print(ds.data_vars)

In **expert** mode:

.. ipython:: python
    :okwarning:

    with argopy.set_options(mode='expert'):
        ds = ArgoDataFetcher(src='localftp').profile(6901929, 2).to_xarray()
        print(ds.data_vars)
Performances
============

.. contents::
   :local:

To improve **argopy** data fetching performances (in terms of time of
retrieval), 2 solutions are available:

-  :ref:`cache` fetched data, i.e. save your request locally so that you don’t have to fetch it again,
-  Use :ref:`parallel`, i.e. fetch chunks of independent data simultaneously.

These solutions are explained below.

Note that another solution from standard big data strategies would be to
fetch data lazily. But since (i) **argopy** post-processes raw Argo data
on the client side and (ii) none of the data sources are cloud/lazy
compatible, this solution is not possible (yet).

Let's start with standard import:

.. ipython:: python
    :okwarning:

    import argopy
    from argopy import DataFetcher as ArgoDataFetcher

Cache
-----

Caching data
~~~~~~~~~~~~

If you want to avoid retrieving the same data several times during a
working session, or if you fetched a large amount of data, you may want
to temporarily save data in a cache file.

You can cache fetched data with the fetchers option ``cache``.

**Argopy** cached data are persistent, meaning that they are stored
locally on files and will survive execution of your script with a new
session. **Cached data have an expiration time of one day**, since this
is the update frequency of most data sources. This will ensure you
always have the last version of Argo data.

All data and meta-data (index) fetchers have a caching system.

The argopy default cache folder is under your home directory at
``~/.cache/argopy``.

But you can specify the path you want to use in several ways:

-  with **argopy** global options:

.. code:: python

   argopy.set_options(cachedir='mycache_folder')

-  in a temporary context:

.. code:: python

   with argopy.set_options(cachedir='mycache_folder'):
       ds = ArgoDataFetcher(cache=True).profile(6902746, 34).to_xarray()

-  when instantiating the data fetcher:

.. code:: python

   ds = ArgoDataFetcher(cache=True, cachedir='mycache_folder').profile(6902746, 34).to_xarray()

.. warning::

  You really need to set the ``cache`` option to ``True``. Specifying only the ``cachedir`` won't trigger caching !

Clearing the cache
~~~~~~~~~~~~~~~~~~

If you want to manually clear your cache folder, and/or make sure your
data are newly fetched, you can do it at the fetcher level with the
``clear_cache`` method.

Start to fetch data and store them in cache:

.. ipython:: python
    :okwarning:

   fetcher = ArgoDataFetcher(cache=True, cachedir='mycache_folder').profile(6902746, 34)
   fetcher.to_xarray();

Fetched data are in the local cache folder:

.. ipython:: python
    :okwarning:

   import os
   os.listdir('mycache_folder')

where we see one hash entries for the newly fetched data and the cache
registry file ``cache``.

We can then fetch something else using the same cache folder:

.. ipython:: python
    :okwarning:

   fetcher2 = ArgoDataFetcher(cache=True, cachedir='mycache_folder').profile(1901393, 1)
   fetcher2.to_xarray();

All fetched data are cached:

.. ipython:: python
    :okwarning:

   os.listdir('mycache_folder')

Note the new hash file from *fetcher2* data.

It is important to note that we can safely clear the cache from the
first *fetcher* data, it won’t remove the *fetcher2* data:

.. ipython:: python
    :okwarning:

   fetcher.clear_cache()
   os.listdir('mycache_folder')

By using the fetcher level clear cache, you make sure that only data
fetched with it are removed, while other fetched data (with other
fetchers for instance) will stay in place.

If you want to clear the entire cache folder, whatever the fetcher used,
do it at the package level with:

.. ipython:: python
    :okwarning:

   argopy.clear_cache()

.. _parallel:

Parallel data fetching
----------------------

Sometimes you may find that your request takes a long time to fetch, or
simply does not even succeed. This is probably because you’re trying to
fetch a large amount of data.

In this case, you can try to let argopy chunks your request into smaller
pieces and have them fetched in parallel for you. This is done with the
argument ``parallel`` of the data fetcher and can be tuned using options
``chunks`` and ``chunksize``.

This goes by default like this:

.. ipython:: python
    :okwarning:

    # Define a box to load (large enough to trigger chunking):
    box = [-60, -30, 40.0, 60.0, 0.0, 100.0, "2007-01-01", "2007-04-01"]
    
    # Instantiate a parallel fetcher:
    loader_par = ArgoDataFetcher(src='erddap', parallel=True).region(box)

you can also use the option ``progress`` to display a progress bar
during fetching:

.. ipython:: python
    :okwarning:

    loader_par = ArgoDataFetcher(src='erddap', parallel=True, progress=True).region(box)
    loader_par

Then, you can fetch data as usual:

.. ipython:: python
    :okwarning:

    %%time
    ds = loader_par.to_xarray()

Number of chunks
~~~~~~~~~~~~~~~~

To see how many chunks your request has been split into, you can look at
the ``uri`` property of the fetcher, it gives you the list of paths
toward data:

.. ipython:: python
    :okwarning:

    for uri in loader_par.uri:
        print("http: ... ", "&".join(uri.split("&")[1:-2]))  # Display only the relevant part of each URLs of URI:

To control chunking, you can use the **``chunks``** option that
specifies the number of chunks in each of the *direction*:

-  ``lon``, ``lat``, ``dpt`` and ``time`` for a **region** fetching,
-  ``wmo`` for a **float** and **profile** fetching.

.. ipython:: python
    :okwarning:

    # Create a large box:
    box = [-60, 0, 0.0, 60.0, 0.0, 500.0, "2007", "2010"]
    
    # Init a parallel fetcher:
    loader_par = ArgoDataFetcher(src='erddap', 
                                 parallel=True, 
                                 chunks={'lon': 5}).region(box)
    # Check number of chunks:
    len(loader_par.uri)

This creates 195 chunks, and 5 along the longitudinale direction, as
requested.

When the ``chunks`` option is not specified for a given *direction*, it
relies on auto-chunking using pre-defined chunk maximum sizes (see
below). In the case above, auto-chunking appends also along latitude,
depth and time; this explains why we have 195 and not only 5 chunks.

To chunk the request along a single direction, set explicitly all the
other directions to ``1``:

.. ipython:: python
    :okwarning:

    # Init a parallel fetcher:
    loader_par = ArgoDataFetcher(src='erddap', 
                                 parallel=True, 
                                 chunks={'lon': 5, 'lat':1, 'dpt':1, 'time':1}).region(box)
    
    # Check number of chunks:
    len(loader_par.uri)

We now have 5 chunks along longitude, check out the URLs parameter in
the list of URIs:

.. ipython:: python
    :okwarning:

    for uri in loader_par.uri:
        print("&".join(uri.split("&")[1:-2])) # Display only the relevant URL part

.. note::
    You may notice that if you run the last command with the `argovis` fetcher, you will still have more than 5 chunks (i.e. 65). This is because `argovis` is limited to 3 months length requests. So, for this request that is 3 years long, argopy ends up with 13 chunks along time, times 5 chunks in longitude, leading to 65 chunks in total.

.. warning::
    The `localftp` fetcher and the `float` and `profile` access points of the `argovis` fetcher use a list of resources than are not chunked but fetched in parallel using a batch queue.

Size of chunks
~~~~~~~~~~~~~~

The default chunk size for each access point dimensions are:

====================== ==================
Access point dimension Maximum chunk size
====================== ==================
region / **lon**       20 deg
region / **lat**       20 deg
region / **dpt**       500 m or db
region / **time**      90 days
float / **wmo**        5
profile / **wmo**      5
====================== ==================

These default values are used to chunk data when the ``chunks``
parameter key is set to ``auto``.

But you can modify the maximum chunk size allowed in each of the
possible directions. This is done with the option
**``chunks_maxsize``**.

For instance if you want to make sure that your chunks are not larger
then 100 meters (db) in depth (pressure), you can use:

.. ipython:: python
    :okwarning:

    # Create a large box:
    box = [-60, -10, 40.0, 60.0, 0.0, 500.0, "2007", "2010"]
    
    # Init a parallel fetcher:
    loader_par = ArgoDataFetcher(src='erddap', 
                                 parallel=True, 
                                 chunks_maxsize={'dpt': 100}).region(box)
    # Check number of chunks:
    len(loader_par.uri)

Since this creates a large number of chunks, let’s do this again and
combine with the option ``chunks`` to see easily what’s going on:

.. ipython:: python
    :okwarning:

    # Init a parallel fetcher with chunking along the vertical axis alone:
    loader_par = ArgoDataFetcher(src='erddap', 
                                 parallel=True, 
                                 chunks_maxsize={'dpt': 100},
                                 chunks={'lon':1, 'lat':1, 'dpt':'auto', 'time':1}).region(box)
    
    for uri in loader_par.uri:
        print("http: ... ", "&".join(uri.split("&")[1:-2])) # Display only the relevant URL part


You can see, that the ``pres`` argument of this erddap list of URLs
define layers not thicker than the requested 100db.

With the ``profile`` and ``float`` access points, you can use the
``wmo`` keyword to control the number of WMOs in each chunks.

.. ipython:: python
    :okwarning:

    WMO_list = [6902766, 6902772, 6902914, 6902746, 6902916, 6902915, 6902757, 6902771]
    
    # Init a parallel fetcher with chunking along the list of WMOs:
    loader_par = ArgoDataFetcher(src='erddap', 
                                 parallel=True, 
                                 chunks_maxsize={'wmo': 3}).float(WMO_list)
    
    for uri in loader_par.uri:
        print("http: ... ", "&".join(uri.split("&")[1:-2])) # Display only the relevant URL part


You see here, that this request for 8 floats is split in chunks with no
more that 3 floats each.

.. note::
    At this point, there is no mechanism to chunk requests along cycle numbers for the ``profile`` access point.

Parallelization methods
~~~~~~~~~~~~~~~~~~~~~~~

They are 2 methods available to set-up your data fetching requests in
parallel:

1. `Multi-threading <https://en.wikipedia.org/wiki/Multithreading_(computer_architecture)>`__
   for all data sources,
2. `Multi-processing <https://en.wikipedia.org/wiki/Multiprocessing>`__
   for *localftp*.

Both options use a pool of
`threads <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ThreadPoolExecutor>`__
or
`processes <https://docs.python.org/3/library/concurrent.futures.html#concurrent.futures.ProcessPoolExecutor>`__
managed with the `concurrent futures
module <https://docs.python.org/3/library/concurrent.futures.html#module-concurrent.futures>`__.

The parallelization method is set with the ``parallel_method`` option of
the fetcher, which can take as values ``thread`` or ``process``.

Methods available for data sources:

=================== ====== ======== =======
**Parallel method** erddap localftp argovis
=================== ====== ======== =======
Multi-threading     X      X        X
Multi-processes            X        
=================== ====== ======== =======

Note that you can in fact pass the method directly with the ``parallel``
option, so that in practice, the following two formulations are
equivalent:

.. ipython:: python
    :okwarning:

   ArgoDataFetcher(parallel=True, parallel_method='thread')
   ArgoDataFetcher(parallel='thread')

Comparison of performances
~~~~~~~~~~~~~~~~~~~~~~~~~~

Note that to compare performances with or without the parallel option,
we need to make sure that data are not cached on the server side. To do
this, we use a very small random perturbation on the box definition,
here on the maximum latitude. This ensures that nearly the same amount of data
will be requested but not cached by the server.

.. ipython:: python
    :okwarning:

    def this_box():
        return [-60, 0, 
               20.0, 60.0 + np.random.randint(0,100,1)[0]/1000, 
               0.0, 500.0, 
               "2007", "2009"]

.. ipython:: python
    :okwarning:

    %%time
    b1 = this_box()
    f1 = ArgoDataFetcher(src='argovis', parallel=False).region(b1)
    ds1 = f1.to_xarray()

.. ipython:: python
    :okwarning:

    %%time
    b2 = this_box()
    f2 = ArgoDataFetcher(src='argovis', parallel=True).region(b2)
    ds2 = f2.to_xarray()

**This simple comparison shows that parallel request is significantly
faster than the standard one.**

Warnings
~~~~~~~~

-  Parallelizing your fetcher is useful to handle large region of data,
   but it can also add a significant overhead on *reasonable* size
   requests that may lead to degraded performances. So, we do not
   recommend for you to use the parallel option systematically.

-  You may have different dataset sizes with and without the
   ``parallel`` option. This may happen if one of the chunk data
   fetching fails. By default, data fetching of multiple resources fails
   with a warning. You can change this behaviour with the option
   ``errors`` of the ``to_xarray()`` fetcher methods, just set it to
   ``raise`` like this:

   .. code:: python

      ArgoDataFetcher(parallel=True).region(this_box()).to_xarray(errors='raise');

You can also use ``silent`` to simply hide all messages during fetching.
#############
API reference
#############

This page provides an auto-generated summary of argopy's API. For more details and examples, refer to the relevant chapters in the main part of the documentation.

.. contents::
   :local:

Top-levels functions
====================

.. currentmodule:: argopy

Fetchers
--------

.. autosummary::
    :toctree: generated/

    DataFetcher
    IndexFetcher

Fetcher access points
---------------------

.. autosummary::
   :toctree: generated/

   DataFetcher.region
   DataFetcher.float
   DataFetcher.profile

.. autosummary::
   :toctree: generated/

   IndexFetcher.region
   IndexFetcher.float
   IndexFetcher.profile

Fetcher methods
---------------

.. autosummary::
   :toctree: generated/

   DataFetcher.load
   DataFetcher.to_xarray
   DataFetcher.to_dataframe
   DataFetcher.to_index

.. autosummary::
   :toctree: generated/

   IndexFetcher.load
   IndexFetcher.to_xarray
   IndexFetcher.to_dataframe
   IndexFetcher.to_csv

Data visualisation
------------------

.. autosummary::
   :toctree: generated/

   DataFetcher.plot
   IndexFetcher.plot
   dashboard


Fetcher properties
------------------

.. autosummary::
   :toctree: generated/

   DataFetcher.uri
   DataFetcher.data
   DataFetcher.index
   IndexFetcher.index


Helpers
-------

.. autosummary::
   :toctree: generated/

   status
   TopoFetcher
   set_options
   clear_cache
   tutorial.open_dataset

Low-level functions
===================

.. currentmodule:: argopy

.. autosummary::
    :toctree: generated/

    show_versions
    utilities.list_available_data_src
    utilities.list_available_data_src
    utilities.list_available_index_src


Dataset.argo (xarray accessor)
==============================

.. currentmodule:: xarray

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor.rst

   Dataset.argo

This accessor extends :py:class:`xarray.Dataset`. Proper use of this accessor should be like:

.. code-block:: python

   >>> import xarray as xr         # first import xarray
   >>> import argopy               # import argopy (the dataset 'argo' accessor is registered)
   >>> from argopy import DataFetcher
   >>> ds = DataFetcher().float([6902766, 6902772, 6902914, 6902746]).load().data
   >>> ds.argo
   >>> ds.argo.filter_qc()


Data Transformation
-------------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Dataset.argo.point2profile
   Dataset.argo.profile2point
   Dataset.argo.interp_std_levels
   Dataset.argo.groupby_pressure_bins

Data Filters
------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Dataset.argo.filter_qc
   Dataset.argo.filter_data_mode
   Dataset.argo.filter_scalib_pres

Processing
----------

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

    Dataset.argo.teos10
    Dataset.argo.create_float_source

Misc
----

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

    Dataset.argo.uid
    Dataset.argo.cast_types

Internals
=========

.. currentmodule:: argopy

File systems
------------

.. autosummary::
    :toctree: generated/

    argopy.stores.filestore
    argopy.stores.httpstore
    argopy.stores.memorystore

.. autosummary::
    :toctree: generated/

    argopy.stores.indexstore
    argopy.stores.indexfilter_wmo
    argopy.stores.indexfilter_box

Fetcher sources
---------------

ERDDAP
^^^^^^

.. autosummary::
    :toctree: generated/

    argopy.data_fetchers.erddap_data.ErddapArgoDataFetcher
    argopy.data_fetchers.erddap_data.Fetch_wmo
    argopy.data_fetchers.erddap_data.Fetch_box

Local FTP
^^^^^^^^^

.. autosummary::
    :toctree: generated/

    argopy.data_fetchers.localftp_data.LocalFTPArgoDataFetcher
    argopy.data_fetchers.localftp_data.Fetch_wmo
    argopy.data_fetchers.localftp_data.Fetch_box

Argovis
^^^^^^^

.. autosummary::
    :toctree: generated/

    argopy.data_fetchers.argovis_data.ArgovisDataFetcher
    argopy.data_fetchers.argovis_data.Fetch_wmo
    argopy.data_fetchers.argovis_data.Fetch_box

Plotters
--------

.. autosummary::
   :toctree: generated/

    argopy.plotters.plot_trajectory
    argopy.plotters.bar_plot
    argopy.plotters.open_dashboard
    argopy.plotters.open_sat_altim_report{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessorcallable:: {{ (module.split('.')[1:] + [objname]) | join('.') }}.__call__{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessormethod:: {{ (module.split('.')[1:] + [objname]) | join('.') }}{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessorattribute:: {{ (module.split('.')[1:] + [objname]) | join('.') }}{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module.split('.')[0] }}

.. autoaccessor:: {{ (module.split('.')[1:] + [objname]) | join('.') }}