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
