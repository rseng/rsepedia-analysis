# Contribution to PyInterpolate

We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

* Reporting a bug
* Discussing the current state of the code
* Submitting a fix
* Proposing new features
* Becoming a maintainer

## Where should I start?

Here, on Github! We use github to host code, to track issues and feature requests, as well as accept pull requests. We have Discord server too and it's available here. It's the fastest way to communicate with package maintainers.

### **[Discord Server](https://discord.gg/3EMuRkj)**

---

## We Use [Github Flow](https://guides.github.com/introduction/flow/index.html), So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase (we use [Github Flow](https://guides.github.com/introduction/flow/index.html)). We actively welcome your pull requests:

1. Fork the repo and create your branch from `dev` or from `main` (preferably `dev`).
2. If you've added code that should be tested, add tests in the `test` package. We use Python's `unittest` package to perform testing.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request!

## Any contributions you make will be under the BSD 3-Clause "New" or "Revised" License
In short, when you submit code changes, your submissions are understood to be under the same [BSD 3-Clause "New" or "Revised" License] that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/szymon-datalions/pyinterpolate/issues)
We use GitHub issues to track public bugs. Report a bug by opening a new issue.

## Write bug reports with detail, background, and sample code
[This is an example](https://github.com/szymon-datalions/pyinterpolate/issues/4)

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
- Be specific!
- Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports. I'm not even kidding.

## Use a PEP8 Guidelines
[PEP8 Python Guidelines](https://www.python.org/dev/peps/pep-0008/)

## License
By contributing, you agree that your contributions will be licensed under its BSD 3-Clause "New" or "Revised" License.

## References
This document was adapted from the open-source contribution guidelines for [Facebook's Draft](https://github.com/facebook/draft-js/blob/a9316a723f9e918afde44dea68b5f9f39b7d9b00/CONTRIBUTING.md)

## Example of Contribution

1. You have an idea to speed-up computation of areal semivariance. You plan to use `multiprocessing` package for it.
2. Fork repo from `dev` branch and at the same time propose change or issue in the [project issues](https://github.com/szymon-datalions/pyinterpolate/issues). You may use two templates - one for **bug report** and other for **feature**. In this case you choose **feature**.
3. Create the new child branch from the forked `dev` branch. Name it as `dev-your-idea`. In this case `dev-areal-multiproc` is decriptive enough.
4. Code in your branch.
5. Create few unit tests in `pyinterpolate/test` directory or re-design actual tests if there is a need. For programming cases write unit tests, for mathematical and logic problems write functional tests. Use data from `sample_data` directory.
6. Multiprocessing maybe does not require new tests. But always run unittests in the `test` directory after any change in the code and check if every test has passed.
7. Run all tutorials too. Their role is not only informational. They serve as a functional test playground.
8. If everything is ok make a pull request from your forked repo.
9. And that's all! For every question use [Discord](https://discord.gg/3EMuRkj).

## Contribution by social networks

Your contribution may be other than coding itself. Questions and issues are important too. Do not be scared to write them!# Contributor Covenant Code of Conduct

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
reported by contacting the project team at simon@ml-gis-service.com. All
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
Installation
------------

Package is working with linux and mac os systems. To install it download package and open it in the terminal then type:

```
pip install pyinterpolate
```

This command runs **setup.py** file inside package and install requirements from the list provided there.

*****

_I'd like to run Jupyter Notebooks, what should I do?_

*****

There is an additional step to run this library in Jupyter Notebooks. Before using pip you have to create conda
environment and install required dependencies:

#### Step 1:

```
conda create -n [NAME OF YOUR ENV]
```

#### Step 2:

```
conda activate [NAME OF YOUR ENV]
```

#### Step 3:

```
conda install -c conda-forge python=3.7 pip notebook
```

#### Step 4:

```
pip install pyinterpolate
```

Now you are able to run library from conda notebooks.

*****

_libspatialindex_c.so dependency error_

*****

Sometimes **rtree** (and / or **GeoPandas**) which are requirements for pyinterpolate may be not installed properly
because your operating system does not have **libspatialindex_c.so** file. In this case install it from terminal:

LINUX:

```
sudo apt install libspatialindex-dev
```

MAC OS:

```
brew install spatialindex
```

*****

_How to install package with virtual environment?_

*****

Coming soon...
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02869/status.svg)](https://doi.org/10.21105/joss.02869) ![License](https://img.shields.io/github/license/DataverseLabs/pyinterpolate) ![Build Status](https://travis-ci.com/DataverseLabs/pyinterpolate.svg?branch=main) ![Documentation Status](https://readthedocs.org/projects/pyinterpolate/badge/?version=latest) ![CodeFactor](https://www.codefactor.io/repository/github/DataverseLabs/pyinterpolate/badge)

![Pyinterpolate](https://github.com/szymon-datalions/pyinterpolate/blob/main/logo.png?raw=true  "Pyinterpolate logo")

**version 0.2.5** - *Huygens Crater*
---------------------------------------

PyInterpolate is designed as the Python library for geostatistics. Its role is to provide access to spatial statistics tools used in many studies. This package helps you **interpolate spatial data** with the *Kriging* technique.

If you’re:

- GIS expert,
- geologist,
- mining engineer,
- ecologist,
- public health specialist,
- data scientist.

Then this package may be helpful for you. You could use it for:

- spatial interpolation and spatial prediction,
- alone or with machine learning libraries,
- for point and areal datasets.

Pyinterpolate allows you to perform:

1. Ordinary Kriging and Simple Kriging (spatial interpolation from points),
2. Centroid-based Kriging of Polygons (spatial interpolation from blocks and areas),
3. Area-to-area and Area-to-point Poisson Kriging of Polygons (spatial interpolation and data deconvolution from areas to points).

How it works
--------------

The package performs multiple spatial interpolation tasks. The flow of analysis is usually the same for each interpolation method:

**[1.]** Read and prepare data.

```python
from pyinterpolate.io_ops import read_point_data

point_data = read_point_data('xyz_txt_file.txt')
```

**[2.]** Analyze data, Semivariance calculation.

```python
from pyinterpolate.semivariance import calculate_semivariance

search_radius = 0.01
max_range = 0.32

experimental_semivariogram = calculate_semivariance(
	data=point_data,
	step_size=search_radius,
	max_range=max_range)
```

**[3.]** Data transformation, theoretical semivariogram.

```python
from pyinterpolate.semivariance impJezero Craterort TheoreticalSemivariogram
semivar = TheoreticalSemivariogram(points_array=point_data, empirical_semivariance=experimental_semivariogram)
number_of_ranges = 32

semivar.find_optimal_model(weighted=False, number_of_ranges=number_of_ranges)
```

**[4.]** Interpolation.

```python
from pyinterpolate.kriging import Krige

model = Krige(semivariogram_model=semivar, known_points=point_data)
unknown_point = (12.1, -5.9)

ok_pred = model.ordinary_kriging(unknown_location=unknown_point, number_of_neighbours=32)
```

**[5.]** Error and uncertainty analysis.

```python
real_val = 10  # Some real, known observation at a given point
squared_error = (real_val - ok_pred[0])**2
print(squared_error)
```

```bash
>> 48.72
```

With **pyinterpolate**, you can retrieve the point support model from areal aggregates. Example from _Tick-borne Disease Detector_ study for European Space Agency - COVID-19 population at risk mapping. It was done with the Area-to-Point Poisson Kriging technique from the package. Countries worldwide present infections as areal sums to protect the privacy of infected people. But this kind of representation introduces bias to the decision-making process. To overcome this bias, you may use Poisson Kriging. Areal aggregates of COVID-19 infection rate are transformed to new point support semivariogram created from population density blocks. As an output, we get a population at risk map:

![Covid-19 infection risk in Poland for 14th April 2020.](https://github.com/szymon-datalions/pyinterpolate/blob/main/deconvoluted_risk_areas.jpg?raw=true  "Covid-19 infection risk in Poland for 14th April 2020.")



Status
------

Beta version: the package is tested, and the main structure is preserved, but future changes are likely to occur.


Setup
-----

Setup by pip: pip install pyinterpolate / **Python 3.7** is required!

Detailed instructions on setting up the package are presented in the file [SETUP.md](https://github.com/szymon-datalions/pyinterpolate/blob/master/SETUP.md). We pointed out there most common problems related to third-party packages.

You may follow those setup steps to create a conda environment with a package for your tests:

### Recommended - conda installation

[1.] First, install system dependencies to use package (```libspatialindex_c.so```):

LINUX:

```
sudo apt install libspatialindex-dev
```

MAC OS:

```
brew install spatialindex
```

[2.] Next step is to create conda environment with Python 3.7, pip, and notebook packages and activate your environment:

```
conda create -n [YOUR NAME] -c conda-forge python=3.7 pip notebook
```

```
conda activate [YOUR NAME]
```

[3.] In the next step, install **pyinterpolate** and its dependencies with ```pip```:

```
pip install pyinterpolate
```

[4.] You are ready to use the package!

### pip installation

With **Python==3.7** and system ```libspatialindex_c.so``` dependencies, you may install the package by simple command:

```
pip install pyinterpolate
```

**Always use Virtual Environment for the installation**.

Tests and contribution
------------------------

All tests are grouped in the `test` directory. To run them, you must have installed the `unittest` package. More about test and contribution is here: [CONTRIBUTION.md](https://github.com/szymon-datalions/pyinterpolate/blob/master/CONTRIBUTION.md)


Commercial and scientific projects where the library has been used
------------------------------------------------------------------

* Tick-Borne Disease Detector (Data Lions) for the European Space Agency (2019-2020).
* B2C project related to the prediction of demand for specific flu medications,
* B2G project related to large-scale infrastructure maintenance.

Community
---------

Join our community in Discord: [Discord Server PyInterpolate](https://discord.gg/3EMuRkj)


Bibliography
------------

PyInterpolate is developed thanks to many resources, and some of them are pointed out here:

- Armstrong M., Basic Linear Geostatistics, Springer 1998,
- GIS Algorithms by Ningchuan Xiao: https://uk.sagepub.com/en-gb/eur/gis-algorithms/book241284
- Pardo-Iguzquiza E., VARFIT: a Fortran-77 program for fitting variogram models by weighted least squares, Computers & Geosciences 25, 251-261, 1999,
- Goovaerts P., Kriging and Semivariogram Deconvolution in the Presence of Irregular Geographical Units, Mathematical Geology 40(1), 101-128, 2008
- Deutsch C.V., Correcting for Negative Weights in Ordinary Kriging, Computers & Geosciences Vol.22, No.7, pp. 765-773, 1996


How to cite
-----------
Moliński, S., (2022). Pyinterpolate: Spatial interpolation in Python for point measurements and aggregated datasets. Journal of Open Source Software, 7(70), 2869, https://doi.org/10.21105/joss.02869


Requirements and dependencies (v 0.2.5)
---------------------------------------

* Python 3.7.6

* Numpy 1.18.3

* Scipy 1.4.1

* GeoPandas 0.7.0

* Fiona 1.18.13.post1 (Mac OS) / Fiona 1.8 (Linux)

* Rtree 0.9.4 (Mac OS), Rtree >= 0.8 & < 0.9 (Linux)

* Descartes 1.1.0

* Pyproj 2.6.0

* Shapely 1.7.0

* Matplotlib 3.2.1

Package structure
-----------------

High level overview:

 - [ ] pyinterpolate
    - [x] **distance** - distance calculation,
    - [x] **idw** - inverse distance weighting interpolation,
    - [x] **io_ops** - reads and prepares input spatial datasets,
    - [x] **transform** - transforms spatial datasets,
    - [x] **viz** - interpolation of smooth surfaces from points into rasters,
    - [x] **kriging** - Ordinary Kriging, Simple Kriging, Poisson Kriging: centroid based, area-to-area, area-to-point,
    - [x] **misc** - compare different kriging techniques,
    - [x] **semivariance** - calculate semivariance, fit semivariograms and regularize semivariogram,
    - [x] **tutorials** - tutorials (Basic, Intermediate and Advanced)

Functions in detail
-------------------

Pyinterpolate https://pyinterpolate.readthedocs.io/en/latest/


Development
-----------

- v 0.3 of the package that runs on the Linux, Mac, and Windows OS; updated and extended variance modeling and analysis.


Known Bugs
-----------------

- Package may crash with a huge dataset (memory issues). Operations are performed with numpy arrays, and for datasets larger than 10 000 points, there could be a memory issue ([Issue page](https://github.com/szymon-datalions/pyinterpolate/issues/64))
Data Structure
==============

Library works with geospatial datasets but they could be different from project to project. To ensure stable calculations specific data structures must be preserved by all scripts. Here we have gathered all data structures used by the library:

POINTS:
------

Points are usually described only by their coordinates and value measured at a given location:

> [coordinate x, coordinate y, value] --> [float, float, float or int]


AREAS:
------

Areas (polygons) are more complex than points and they are described by their id (optional parameter), geometry (shapely.geometry.polygon.Polygon), centroid, value:

> [area id name, [geometry ... ... ... ... ...], [centroid coordinate x, centroid coordinate y], value] --> [str or int, Polygon, list, float or int]


POINTS WITHIN AREA:
-------------------

Points within area are described by the area id and a list of all points coordinates and their values:

> [area id name, [[coordinate x, coordinate y, value], ..., [coordinate x, coordinate y, value]]]
# <Feature Title>

## Package version (main branch)
version: **Insert number here**

## Description
Short description of **feature** or **enhancement** or **issue**.

## Problem
Why this is a problem? How it affects users?

## Solution
Describe solution for the presented problem.

### Affected modules

- list of modules which are affected by the solution

### Unit tests

- list of the new unit tests for presented feature

### Package check

- [ ] All tests passed
- [ ] Documentation updated
- [ ] All tutorials are working properly

### (Optional) Additional info
Is feature related to literature? Does change require new dependencies? Any other information...


---
name: Bug report
about: Create a report to help us improve
title: "[BUG]"
labels: ''
assignees: ''

---

**Package version**
Package version which you are using.

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run function / module with parameters...

**Input Data Schema**
Describe your input data:
- types of columns,
- dataset size,
- are there missing values?
- are there duplicates?
- are there negative values?

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]

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
---
title: 'Pyinterpolate: Spatial interpolation in Python for point measurements and aggregated datasets'
tags:
  - Python
  - spatial statistics
  - spatial interpolation
  - Poisson Kriging
  - Semivariogram Deconvolution
  - public health
authors:
  - name: Szymon Moliński
    orcid: 0000-0003-3525-2104
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 20 October 2020
bibliography: paper.bib
---

# Summary

We use spatial interpolation techniques to interpolate values at unknown locations or filter and smooth existing data sources. Those methods work for point observations and areal aggregates. The basic idea behind the spatial interpolation algorithms is that every point in space can be described as a function of its neighbors’ values weighted by the relative distance from the analyzed point. It is known as Tobler's First Law of Geography, which states: *everything is related to everything else, but near things are more related than distant things* [@Tobler:1970].

The Kriging technique, originally designed for mining applications, exploits this statement formally, and nowadays, it has gained a lot of attention outside the initial area of interest. Today Kriging is a set of methods applied to problems from multiple fields: environmental science, hydrogeology, natural resources monitoring, remote sensing, epidemiology and ecology, and even computer science [@Chilès:2018]. Commonly, Kriging is used to interpolate values from point measurements or regular block units. However, the real-world datasets are often different. Especially challenging is data that represents aggregated values over polygons, for example, the administrative units [@Goovaerts:2007].

Pyinterpolate transforms areas of irregular shapes and sizes with Area-to-Area and Area-to-Point Poisson Kriging functions. Those algorithms make Pyinterpolate beneficial for social, environmental, and public health scientists because they usually deal with areal counts instead of point measurements. Moreover, the package offers basic point Kriging and Inverse Distance Weighting techniques. Those algorithms are used in every field of research where geostatistical (distance) analysis gives meaningful results. Pyinterpolate merges basic Kriging techniques with more sophisticated Area-to-Area and Area-to-Point Poisson Kriging methods.


# Statement of need

Pyinterpolate is a Python package for spatial interpolation. It performs predictions from point measurements and areal aggregates of different sizes and shapes. Pyinterpolate automates Kriging interpolation, and semivariogram regularization. The package helps with data exploration, data preprocessing, and semivariogram analysis. A researcher with geostatistical background has control over the basic modeling parameters: semivariogram models, nugget, sill, range, and the number of neighbors included in the interpolation and Kriging type. The thing that makes Pyinterpolate different from other spatial interpolation packages is the ability to perform Kriging on areas of different shapes and sizes. This type of operation is essential in social, medical and ecological sciences [@Goovaerts:2007; @Goovaerts:2008; @Kerry:2013].

## Importance of areal (block) Kriging

There are many applications where researchers need to model areal data with irregular shapes and sizes. A good example is the public health sector, where data is aggregated over administrative units for patient protection and policy-making purposes. Unfortunately, this transformation makes data analysis and modeling more complex for researchers. There are different techniques to deal with this kind of data. We can work directly with areal aggregates or fall the irregular polygons into their centroids or, finally, transform dataset into a regular grid of smaller blocks if the point-support model is available. The latter is not a way to *get back* original observations but rather a form of lossy semivariogram transformation to the point-support scale. There are reasons to do it:

1. The presence of extremely unreliable rates that typically occur for sparsely populated areas and rare events. Consider the examples with the number of leukemia cases (numerator) per population size in a given county (denominator) or the number of whales observed in a given area (numerator) per time of observation (denominator). In those cases, extreme values may be related to the fact that variance for a given area of interest is high (low number of samples) and not to the fact that the chance of the event is exceptionally high for this region.
2. The visual bias. People tend to give more importance to large blocks in contrary to the small regions.
3. The mismatch of spatial supports for aggregated data and other variables. Data for spatial modeling should have harmonized spatial scale and the same extent. The aggregated datasets are not an exception. It may lead to the trade-off where we must aggregate other variables to build a model. Unfortunately, we lost a lot of information in this case. The other problem is that administrative regions are artificial constructs and aggregation of variables may remove spatial trends from data. A downscaling of areal data into filtered population blocks may be better suited to risk estimation along with remote-sensed data or in-situ observations of correlated variables [@Goovaerts:2006].

In this context, Area-to-Area Poisson Kriging serves as the noise-filtering algorithm or areal interpolation model, and Area-to-Point Poisson Kriging interpolates and transforms values and preserves the prediction coherence (where the disaggregated estimates sum is equal to the baseline area value) [@Goovaerts:2008]. The chained-pipelines may utilize Area-to-Point Poisson Kriging, especially if scientist needs to change the support of variables. The author created a model of this type, the machine-learning pipeline with a model based on the remote-sensing data was merged with the geostatistical population-at-risk model derived from the Area-to-Point Poisson Kriging (the research outcomes are not published yet).

Alternatively to the Area-to-Area and Area-to-Point Poisson Kriging, researchers may use centroids and perform point kriging over a prepared regular point grid. However, this method has its pitfalls. Different sizes and shapes of the baseline units lead to the imbalanced number of variogram point pairs per lag. The centroid-based approach misses spatial variability of the linked variable, for example, population density over an area in the context of infection rates.



# Methodology

Pyinterpolate performs six types of spatial interpolation; inverse distance weighting and five types of Kriging:

1. **Ordinary Kriging** is a universal method for point interpolation.
2. **Simple Kriging** is a special case of point interpolation when the mean of the spatial process is known and does not vary spatially in a systematic way.
3. **Centroid-based Poisson Kriging** is used for areal interpolation and filtering. We assume that each block can collapse into its centroid. It is much faster than Area-to-Area and Area-to-Point Poisson Kriging but introduces bias related to the area's transformation into single points.
4. **Area-to-Area Poisson Kriging** is used for areal interpolation and filtering. The point-support allows the algorithm to filter unreliable rates and makes final areal representation of rates smoother.
5. **Area-to-Point Poisson Kriging** where areal support is deconvoluted in regards to the point support. Output map has a spatial resolution of the point support while coherence of analysis is preserved (sum of rates is equal to the output of Area-to-Area Poisson Kriging). It is used for point-support interpolation and data filtering.

The theory of Kriging is described in supplementary materials in the [paper repository](https://github.com/SimonMolinsky/pyinterpolate-paper/blob/main/paper/supplementary%20materials/theory_of_kriging.md) or in more detail in @Armstrong:1998. @OliverWebster:2015 point to the practical aspects of Kriging. The procedure of the interpolation with Poisson Kriging is presented in @Goovaerts:2006 and the semivariogram regularization process is described in @Goovaerts:2007.

The comparison to existing software is presented in the supplementary document [here](https://github.com/SimonMolinsky/pyinterpolate-paper/blob/main/paper/supplementary%20materials/comparison_to_gstat.md), Ordinary Kriging outcomes are compared for *gstat* and Pyinterpolate.

## Interpolation steps

The user starts with semivariogram exploration and modeling. Next, the researcher, or automatically with an algorithm, chooses the theoretical model which best fits the semivariogram. If this is done automatically, the algorithm tests linear, spherical and exponential models with different sills and ranges and the constant nugget against the experimental curve. Model performance is measured by the root mean squared error between the tested theoretical model with the experimental semivariance. 

Areal data interpolation, especially transformation from areal aggregates into point support maps, requires deconvolution of areal semivariogram. Users may do it without prior knowledge of kriging and spatial statistics because an operation is automated. The iterative procedure of the semivariogram regularization is described in detail in @Goovaerts:2007. The last step of analysis is a solution of linear Kriging equations.

Predicted data is stored as a `DataFrame` known from the *Pandas* and *GeoPandas* Python packages. Pyinterpolate allows the user to transform the point data into a regular Numpy array grid for further processing and analysis. Use case with the whole scenario is available in the [paper package repository](https://github.com/szymon-datalions/pyinterpolate-paper).

The package can automatically perform the semivariogram fitting step with a derivation of the theoretical semivariogram from the experimental curve. The semivariogram regularization is entirely automated. The process is described in @Goovaerts:2007. Users can change the derived theoretical model only by directly overwriting the derived semivariogram model parameters (nugget, sill, range, model type).

The initial field of study (epidemiology) was the reason behind the automation of the tasks related to semivariogram modeling. Pyinterpolate was initially developed for the epidemiological research, where areal aggregates of infections were transformed to point support population-at-risk maps. It is assumed that users without a broad geostatistical background may use Pyinterpolate for spatial data modeling and analysis, especially users observing processes related to the human population.

The \autoref{fig1} is an example of a full-scale process of the semivariogram regularization and Area-to-Point Poisson Kriging.

![Example use case of Pyinterpolate for the derivation of the population-at-risk map for a cancer development from the areal aggregates and the population blocks.\label{fig1}](fig1_example.png)

The repository [here](https://github.com/SimonMolinsky/pyinterpolate-paper/) presents an example of Poisson Kriging of cancer rates in North-Eastern U.S. step-by-step, with semivariogram regularization and Point Poisson Kriging functions. This repository contains three documents with additional information:

1. [*IPython* notebook with code](https://github.com/SimonMolinsky/pyinterpolate-paper/blob/main/paper-examples/example-use-case/joss%20paper%20example.ipynb).
2. [Document with detailed description of methodology](https://github.com/SimonMolinsky/pyinterpolate-paper/blob/main/paper/supplementary%20materials/example_use_case.md).
3. [Document that describes the areal data transformation process](https://github.com/SimonMolinsky/pyinterpolate-paper/blob/main/paper/supplementary%20materials/areal_data_transformation.md). This procedure follows [@Goovaerts:2007].

## Modules

Pyinterpolate has seven modules covering all operations needed to perform spatial interpolation: input/output operations, data processing, transformation, semivariogram fitting, Kriging interpolation. \autoref{fig2} shows the internal package structure.

![Structure of Pyinterpolate package.\label{fig2}](fig2_modules.png)

# Comparison to Existing Software

The main difference between Pyinterpolate and other packages is that it focuses on areal deconvolution methods and Poisson Kriging techniques useful for ecology, social science and public health studies. Potential users may choose other packages if they can perform their research with the point data interpolation.

The most similar and significant package from the Python environment is *PyKrige* [@benjamin_murphy_2020_3991907]. PyKrige is designed especially for point kriging. PyKrige supports 2D and 3D ordinary and universal Kriging. User can incorporate their own semivariogram models and use external functions (as an example from *scikit-learn* package [@scikit-learn]) to model drift in universal Kriging. The package is actively maintained.

*GRASS GIS* [@GRASS_GIS_software] is a well-established software for vector and raster data processing and analysis. GRASS contains multiple modules and a user may access them in numerous ways: GUI, command line, C API, Python APU, Jupyter Notebooks, web, QGIS or R. GRASS has three functions for spatial interpolation: 

- [`r.surf.idw`](https://grass.osgeo.org/grass78/manuals/r.surf.idw.html) and [`v.surf.idw`](https://grass.osgeo.org/grass78/manuals/v.surf.idw.html): both use Inverse Distance Weighting technique, first interpolates raster data and second vectors (points).
- [`v.surf.rst`](https://grass.osgeo.org/grass78/manuals/v.surf.rst.html) that performs surface interpolation from vector points map by splines. Spline interpolation and Kriging are compared in [@DUBRULE1984327].

*PySAL* is the next GIS / geospatial package that is used for spatial interpolation. However, PySAL is built upon the spatial graph analysis algorithms. The areal analysis is performed with the sub-module [*tobler*](https://github.com/pysal/tobler) [@eli_knaap_2020_4385980]. Moreover, the package has functions for multisource regression, where raster data is used as auxiliary information to enhance interpolation results.

The R *gstat* package is another option for spatial interpolation and spatial modeling [@PEBESMA2004683]. The package is designed for variogram modeling, simple, ordinary and universal point or block kriging (with drift), spatio-temporal kriging and sequential Gaussian (co)simulation. Gstat is a solid Kriging and spatial interpolation package and has the largest number of methods to perform spatial modeling. The main difference between gstat and Pyinterpolate is the availability of area-to-point Poisson Kriging in the latter and the difference between baseline programming languages [@Goovaerts:2007]. The functional comparison to gstat is available in the [paper repository](https://github.com/SimonMolinsky/pyinterpolate-paper).



# References

# Example Use Case – Breast Cancer Rates in Northeastern United States

This appendix to the paper presents sample pipeline of analysis. It is presented how Area-to-Point Poisson Kriging transforms areal counts of breast cancer cases in the Northeastern part of United States.

## Dataset

Breast cancer rates are taken from the *Incidence Rate Report for U.S.* counties and were clipped to the counties of the Northeastern part of U.S. [1]. Observations are age-adjusted and multiplied by 100,000 for the period 2013-2017.

Population centroids are retrieved from the *U.S. Census Blocks 2010* [2]. Breast cancer affects only females but for this example the whole population for an area was included. Raw and transformed datasets are available in a [dedicated Github repository](https://github.com/SimonMolinsky/pyinterpolate-paper/tree/main/paper-examples/example-use-case).

Presented work is Area-to-Point Poisson Kriging of Breast Cancer areal aggregates dataset and transformation of those areal aggregates into population-specific blocks (points). This process requires two main steps: **semivariogram regularization** and **Poisson Kriging**.

## 1. Read and prepare data

The initial step of analysis is data preparation. Pyinterpolate transforms passed shapefiles to numpy arrays containing geometries and values for processing. Areal data is transformed into its `id, geometry, centroid x, centroid y` and `value`; point data is transformed into `points` and their `values` within the area with specific `id`.

## 2. Analyze and test semivariance of points or areas

Experimental semivariogram model of data must be retrieved at the beginning of analysis. Semivariance of areal centroids (Figure 1) and semivariance of point support (Figure 2) should be checked to be sure that process is spatially correlated, especially at the scale of point support. User selects the maximum range of analysis - *study extent* - and step size for each lag. The package calculates experimental semivariance based on the provided input.

![Experimental semivariogram of areal centroids.\label{Figure 1}](experimental_areal.png)

*Figure 1: Experimental semivariogram of areal centroids.*

![Experimental semivariogram of point support.\label{Figure 2}](experimental_pop.png)

*Figure 2: Experimental semivariogram of point support.*

Semivariogram of areal rates shows weak spatial autocorrelation but this may be the effect of data aggregation and large differences in blocks shape and size. The semivariogram of point support presents better spatial autocorrelation pattern and it reaches sill at a distance of 100 kilometers.

## 3. Create theoretical semivariogram or regularize areal aggregated values

Semivariogram modeling is fully automated and best model is selected based on the lowest error between chosen model type from *spherical*, *linear*, *gaussian* or *exponential* models and the experimental curve.

Deconvolution of areal semivariogram is a more complex problem and it's algorithm is described in [3].  Pyinterpolate implementation divides semivariogram regularization into two parts. First part is an initial preparation of data and a development of the first optimized theoretical model. In a second step a areal semivariogram is regularized in a loop. It is a time consuming process. Computation time directly depends on the number of points of the support.

Experimental semivariogram and theoretical model of areal data along with first output of regularization may be checked before the main loop to be sure that process can be modeled with Kriging method. Figure 3 presents initial (baseline) semivariograms and Figure 4 shows those after regularization. After the procedure we are able to export model for the Poisson Kriging interpolation.

![Semivariograms after fit procedure.\label{Figure 3}](fit0.png)

*Figure 3: Semivariograms after the initial fit procedure - deviation is big and it needs to be corrected.*

![Semivariograms after regularization.\label{Figure 4}](reg.png)

*Figure 4: Semivariograms after complete regularization.*

Regularized Semivariogram has much more visible spatial component and sill is reached at 100 kilometers instead of half of this value. This model may be used for Poisson Krigin interpolation.

## 4. Build Kriging model and export output

With theoretical semivariogram we are able to model data with Kriging. Poisson Kriging model is used to estimate population at risk from areal aggregates. Area-to-Point Poisson Kriging requires us to know the semivariogram model and to assign the number of closest neighbors and maximum radius of neighborhood search.

Whole process may take a while, especially if there are many support points. Method `regularize_data()` returns *GeoDataFrame* object with `[id, geometry, estimated value, estimated prediction error, rmse]` columns. It may be plotted with *matplotlib* and as a result **population at risk** map is generated (Figure 5). Finally, point support map may be saved as a shapefile.

![Breast cancer population at risk map in Northeastern United States state.\label{Figure 5}](smoothed.png)

*Figure 5: Breast cancer population at risk map in Northeastern United States state.*

Comparison of input and output data in this example is presented in Figure 6. Output values and error of variance may be used later for reporting and / or as the elements of larger modeling infrastructure.

 ![Report output.\label{Figure 6}](fig1_example.png)
 
 *Figure 6: Sample report output.*
 
 ---
 
# Bibliography
 
**1) Cancer Data**
 
```
{,
  author = {National Cancer Institute},
  title = {Incidence Rates Table: Breast Cancer: Pennsylvania State},
  year = 2021,
  url = {https://www.statecancerprofiles.cancer.gov/incidencerates/index.php?stateFIPS=42&areatype=county&cancer=055&race=00&sex=2&age=001&stage=999&year=0&type=incd&sortVariableName=rate&sortOrder=default&output=0#results},
  urldate = {2021-02-03}
}
```

**2) Census Data**

```
{,
  author = {United States Census Bureau},
  title = {Centers of Population for the 2010 Census},
  year = 2021,
  url = {https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.2010.html},
  urldate = {2021-02-03}
}
```

**3) Goovaerts, P. Kriging and Semivariogram Deconvolution in the Presence of Irregular Geographical Units**

```
{Goovaerts:2007,
  author = {Goovaerts, Pierre},
  month = {12},
  pages = {101-128},
  title = {Kriging and Semivariogram Deconvolution in the Presence of Irregular Geographical Units},
  doi = {10.1007/s11004-007-9129-1},
  urldate = {2020-10-24},
  volume = {40},
  year = {2007},
  journal = {Mathematical Geosciences}
}
```
 # Comparison of Ordinary Kriging in `Pyinterpolate` and `gstat` package

Up to date the most popular package for Kriging is **gstat** [1] written in the R programming language. Existence of this software allows to compare **point Kriging** operations between packages. The `meuse` dataset provided with the R package **sp** [2; 3] is used in this example. This dataset contains measurements of four heavy metals concentration in the top soil in a flood plain along the river Meuse [4]. Code for this comparison is given in the dedicated notebook available in the [paper repository](https://github.com/SimonMolinsky/pyinterpolate-paper/tree/main/paper-examples/example-gstat-comparison).

From four metals concentration of zinc is used in this example. Variogram modeling and kriging are performed semi-automatically in both cases with a common setting related to the kriging itself - with the number of neighbours. Variogram modeling is different. In the case of **gstat** package variogram was derived as a Spherical model with nugget equal to the first three values of experimental semivariances, sill is equal to the mean of the last five values of experimental semivariances and range which is 1/3 of the maximum distance. Lags between steps are not equal. On the contrary, **pyinterpolate** fits a semivariogram automatically based on the lowest RMSE between theoretical model and experimental values. It is an iterative process. Function searches for the optimal range and sill pair and model type; it performs grid search. Model, sill and range with the lowest RMSE is selected as the optimal model. This process is repeated for each type of theoretical model. In this case the Spherical model is the best. Nugget is equal to zero. Interpolation grid (points) is derived from the **sp** package.

The ordinary Kriging interpolation is performed for both packages. Predictions and prediction variance error terms are calculated and compared. Scatterplot of **pyinterpolate** output versus **gstat** output is presented in the Figure 1. Calculated Pearson correlation coefficient between both outputs is extremely high and close to the 0.99 with p-value close to the 0. Small differences may be related to the slightly different ranges of both models.

 ![Correlation between predicted values from the **pyinterpolate** package and the gstat package.\label{Figure 1}](gstatxpyint.png)
 
*Figure 1: Correlation between predicted values from the pyinterpolate package and the gstat package.*

Similar pattern may be observed for error variance estimates:

 ![Correlation between the predicted error variance from the **pyinterpolate** package and the gstat package.\label{Figure 2}](gstaterrxpyinterr.png)
 
*Figure 2: Correlation between the predicted error variance from the pyinterpolate package and the gstat package.*

 
 Pattern of predicted values and variance errors are very similar in both cases. Figure 3 shows predicted output of **gstat** package and **pyinterpolate** and Figure 4 shows maps of variance errors from the both packages.
 
  ![Comparison of predicted values from gstat (on the left) and pyinterpolate (on the right) packages.\label{Figure 3}](predicted_comparison.png)
  
 *Figure 3: Comparison of predicted values from gstat (on the left) and pyinterpolate (on the right) packages.*
 
  ![Comparison of prediction errors from gstat (on the left) and pyinterpolate (on the right) packages.\label{Figure 4}](errors_comparison.png)
  
*Figure 4: Comparison of prediction errors from gstat (on the left) and pyinterpolate (on the right) packages.*
  
---
  
# Bibliography

**1) Pebesma E.J., Multivariable geostatistics in S: the gstat package**

```
{,
title = "Multivariable geostatistics in S: the gstat package",
journal = "Computers & Geosciences",
volume = "30",
number = "7",
pages = "683 - 691",
year = "2004",
issn = "0098-3004",
doi = "https://doi.org/10.1016/j.cageo.2004.03.012",
url = "http://www.sciencedirect.com/science/article/pii/S0098300404000676",
author = "Edzer J Pebesma",
keywords = "Kriging, Cokriging, Linear model of coregionalisation, Open source software, S language, Stochastic simulation"
}
```

**2) Pebesma E.J. and Bivand R.S., Classes and methods for spatial data in {R}**

```
{,
    author = {Edzer J. Pebesma and Roger S. Bivand},
    title = {Classes and methods for spatial data in {R}},
    journal = {R News},
    year = {2005},
    volume = {5},
    number = {2},
    pages = {9--13},
    month = {November},
    url = {https://CRAN.R-project.org/doc/Rnews/},
  }
```

**3) Bivand R.S., Pebesma E.J. and Gomez-Rubio V., Applied spatial data analysis with {R}, Second edition**

```
{,
    author = {Roger S. Bivand and Edzer Pebesma and Virgilio
      Gomez-Rubio},
    title = {Applied spatial data analysis with {R}, Second edition},
    year = {2013},
    publisher = {Springer, NY},
    url = {https://asdar-book.org/},
  }
```

**4) Pebesma E.J., The meuse data set: a tutorial for the gstat R package**

```
{,
    author = {Edzer J. Pebesma},
    title = {The meuse data set: a tutorial for the gstat R package},
    year = {2021},
    url = {https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf},
}
```# The Theory of Kriging

Kriging is an estimation method that gives the best unbiased linear estimates of point values or block averages [1]. It is the core method of the Pyinterpolate package.
The primary technique is the Ordinary Kriging. The value at unknown location $\hat{z}$ is estimated as a linear combination of $K$ neighbors with the observed values $z$ and weights $\lambda$ assigned to those neighbors (1).

(1) $$\hat{z} = \sum_{i=1}^{K}\lambda_{i}z_{i}$$

Weights $\lambda$ are a solution of following system of linear equations (2):

(2) $$\sum_{j=1}^{K}\lambda_{j} C(x_{i}, x_{j}) - \mu = \bar{C}(x_{i}, V); i=1, 2, ..., K$$ $$\sum_{i}^{K}\lambda_{i} = 1$$

where $C(x_{i}, x_{j})$ is a covariance between points $x_{i}$ and $x_{j}$, $\bar{C}(x_{i}, V)$ is an average covariance between point $x_{i}$ and all other points in a group ($K$ points) and $\mu$ is the Lagrange multiplier. The same system may be solved with semivariance instead of covariance (3):

(3) $$\sum_{i=1}^{K}\lambda_{j} \gamma(x_{i}, x_{j}) + \mu = \bar{\gamma}(x_{i}, V); i=1, 2, ..., K$$ $$\sum_{i}^{K}\lambda_{i} = 1$$

where $\gamma(x_{i}, x_{j})$ is a semivariance between points $x_{i}$ and $x_{j}$, $\bar{\gamma}(x_{i}, V)$ is an average semivariance between point $x_{i}$ and all other points.
Semivariance is a key concept of spatial interpolation. It is a measure of dissimilarity between observations in a function of distance. Equation (4) is an experimental semivariogram estimation formula and $\gamma_{h}$ is an experimental semivariance at lag $h$:

(4) $$\gamma_{h} = \frac{1}{2N}\sum_{i}^{N}(z_{(x_{i} + h)} - z_{x_{i}})^{2}$$

where $z_{x_{i}}$ is a value at location $x_{i}$ and $z_{(x_{i} + h)}$ is a value at a translated location in a distance $h$ from $x_{i}$.

Pyinterpolate package implements linear, spherical, exponential and Gaussian models [1]. They are fitted to the experimental curve. The model with the lowest error is used in (3) to estimate the $\gamma$ parameter.

**Simple Kriging** is another method for point interpolation in Pyinterpolate. We may use Simple Kriging when we know the process mean. This situation rarely occurs in real-world scenarios. It is observed in places where sampling density is high [1]. Simple Kriging system is defined as:

(5) $$\hat{z} = R + \mu$$

where $\mu$ is a Lagrange multiplier and $R$ is a residual at a specific location. The residual value is derived as the first element (denoted as $\boldsymbol{1}$) from:

(6) $$R = ((Z - \mu) \times \lambda)\boldsymbol{1}$$

The number of values depends on the number of neighbors in a search radius, similar to equation (1) for Ordinary Kriging. The weights $\lambda$ are the solution of the following function:

(7) $$\lambda = K^{-1}(\hat{k})$$

The $K$ denotes a semivariance matrix between each neighbor of size $NxN$. The $k$ parameter is a semivariance between unknown (interpolated) location and known points of size $Nx1$.

Users may use three types of Poisson Kriging procedure: Centroid-based Poisson Kriging, Area-to-Area Poisson Kriging and Area-to-Point Poisson Kriging. Each defines the risk over areas (or points) similarly to the equation (1). However, the algorithm estimates the weights associated with the $\lambda$ parameter with additional constraints related to the population weighting. The spatial support of each unit needs to be accounted for in both the semivariogram inference and kriging. The procedure of Poisson Kriging interpolation of areal data is presented in [2] and semivariogram deconvolution in [3].

## Bibliography

[1] Armstrong, M. (1998). Basic Linear Geostatistics. Springer. https://doi.org/10.1007/
266978-3-642-58727-6

[2] Goovaerts, P. (2006). Geostatistical analysis of disease data: Accounting for spatial support and population density in the isopleth mapping of cancer mortality risk using area-to-point poisson kriging. International Journal of Health Geographics, 5. https://doi.org/10.1186/1476-072X-5-52

[3] Goovaerts, P. (2007). Kriging and semivariogram deconvolution in the presence of irregular geographical units. Mathematical Geosciences, 40, 101–128. https://doi.org/10.1007/s11004-007-9129-1# Areal data transformation

To disaggregate areal data into the point support, one must know a regionalized variable’s point support covariance or semivariance. Then the semivariogram deconvolution is performed. In this iterative process, the experimental semivariogram of areal data is transformed to fit the semivariogram model of a linked point support variable. Journel and Huijbregts [1] presented a general approach to deconvolute regularized semivariogram:

1. Define a point-support model from inspection of the semivariogram od areal data and estimate the parameters (sill and range) using basic deconvolution rules.
2. Compute the theoretically regularized model and compare it to the experimental curve.
3. Adjust the parameters of the point-support model to bring them in line with the regularized model.

Pyinterpolate follows an extended procedure. It leads to the automatic semivariogram regularization. [2] described this process in detail. The procedure has ten steps:

1. Compute the experimental semivariogram of areal data and fit a theoretical model to it.
2. The algorithm compares a few theoretical models and calculates the error between a modeled curve and the experimental semivariogram. The algorithm selects the theoretical model with the lowest error as the initial point-support model.
3. The initial point-support model is regularized according to the procedure given in [2].
4. Quantify the deviation between the initial point-support model and the theoretically regularized model.
5. The initial point-support model, the regularized model and the associated deviation are considered optimal at this stage.
6. Iterative process begins: for each lag, the algorithm calculates the experimental values for the new point-support semivariogram. Those values are computed through a rescaling of the optimal point support model available at this stage.
7. The rescaled values are fitted to the new theoretical model in the same procedure as the second step.
8. The new theoretical model (from step 7.) is regularized.
9. Compute the difference statistic for the new regularized model (step 8.). Decide what to do next based on the value of the new difference statistic. If it is smaller than the optimal difference statistic, use the point support model (step 7.) and the associated statistic as the optimal point-support model and the optimal difference statistic. Repeat steps from 6. to 8. If the difference statistic is larger or equal to the optimal difference statistic, repeat steps 6 through 8 with a change of the rescaling weights.
10. Stop the procedure after i-th iteration whenever one of the specified criteria are met: (1) the difference statistic reaches a sufficiently small value, (2) the maximum number of iterations has been tried, (3) a small decrease in the difference statistic was recorded a given number of times.

## Bibliography

[1] Journel AG, Huijbregts CJ. Mining geostatistics. Academic Press; London: 1978. p. 600.

[2] Goovaerts, P. (2007). Kriging and semivariogram deconvolution in the presence of irregular geographical units. Mathematical Geosciences, 40, 101–128. https://doi.org/10.1007/s11004-007-9129-1# Tests - step by step tutorial

In this tutorial, we will learn about the different types of tests that could cover the scientific package. We are going to create unit and functional tests for the code developed within the [Pyinterpolate](https://pypi.org/project/pyinterpolate/) package. We see how to use Jupyter Notebooks and additional tutorial files for testing.

## Changelog

| Date | Change description | Author |
|--------|--------------------------|----------|
| 2021-09-05 | First release of tutorial | @szymon-datalions |

## Reader-friendly version

Reader-friendly version of this tutorial is available [HERE](https://ml-gis-service.com/index.php/2021/09/05/data-engineering-the-test-coverage-for-the-scientific-software/).

## Contents

- [Introduction](#introduction)
- [(Optional) Example function](#example-function)
- [Unit Tests](#unit-tests)
- [Logical and Functional Testing](#logical-and-functional-testing)
- [How to run unit tests](#how-to-run-unit-tests)
  - [PyCharm](#pycharm)
    - [Single Test](#pycharm-single-test)
    - [Multiple Tests](#pycharm-multiple-tests)
  - [Console](#console)
    - [Single Test](#console-single-test)
    - [Multiple Tests](#console-multiple-tests)
- [Testing with tutorials](#use-a-tutorial-for-testing)
- [Summary](#summary)

## Introduction

Whenever code is changed or a new feature is added we should perform tests to ensure that everything works as it should be. We have few different levels of test to perform. The list below represents all steps within the Pyinterpolate package required to be sure that everything works (hopefully) fine.

1. Write unit tests in `test` directory of the package to check if input data structure is valid or if returned values are within specific range or are of specific type. Use `unittests` package. Group tests in regards to the module within you're working. As example if you implement new Kriging technique put your tests inside `test.kriging` module. Name your test file with a prefix `test_` and follow `unittest` naming conventions of functions and classes.
2. Write a logical test where you know exactly what should be returned. This step is very important in the development of scientific software so pen & paper are equally important as keybord & monitor. Where it's possible use examples from the literature and try to achieve the same results with the same inputs. If this is not possible justify why function is usable and why your results are different than those in the literature.
3. Run all tests within the `test` package. You have two main options: use `PyCharm` or `Python` console for it.
4. Update all tutorials where the feature change or the new feature development may affect process or write the new tutorial with the new feature usage included.

To make things easier to understand we will go through the example of `calculate_seimvariance()` function from the `Pyinterpolate` and we will uncover the points above in a greater detail.

## Example function

To start with the development we must first do two things:

1. Write an equation or create block diagram algorithm of the function,
2. Prepare dataset for logic tests.

Those two steps prepare our **mental model**. In the ideal situation, we should have the mathematical equation or algorithm blocks and a sample dataset from publications. Fortunately for us, the experimental semivariogram calculation process is well described in the book **Basic Linear Geostatistics** written by *Margaret Armstrong* in 1998 (pages 47-52 for this tutorial). (If you're a geostatistican and you haven't read this book yet then do it as soon as you can. This resource is a gem among books for geostatistics and `Pyinterpolate` relies heavily on it).

Let's start from the equation for experimental semivariogram:

$$\gamma'(h) = \frac{1}{2 * N(h)} \sum_{i=1}^{N(h)} [Z(x_{i} + h) - Z(x_{i})]^2$$

where $\gamma'(h)$ is a semivariance for a given lag $h$, $x_{i}$ are the locations of samples, $Z(x_{i})$ are their values and $N(h)$ is the number of pairs $(x_{i}, x_{i + h})$. What it means in practice? We may freely translate it to the sentence that **semivariance at a given interval of distances is a halved mean squared error of all points pairs which are in a given interval**. If we understand what it means then we could go further. Based on this assumption we can preapre the first draft of our algorithm as a set of steps to do:

* (1) Read Data as an array of triplets `x, y, value`,
* (2) Calculate distances between each element from the array,
* (3) Create lags list to group semivariances (lags are separation distances $h$),
* (4a) For each lag group calculated distances to select points within a specific range (lag),
* (4b) Calculated mean squared error between all points pairs and divide it by two (each pair is taken twice a time),
* (4c) Store calculated semivariance along the lag and number of points used for calculation within array `lag, semivariance, number of points`,
* (5) Return array of lags and their semivariances.

After this step we should probably see some limits and dependencies within the algorithm output. We are ready to implement the **unit tests**.

## Unit tests

We will create a simple unit test that checks if all results are positive. (Lags can be only positive because distances are not negative. The semivariance is positive too, since it is a squared difference it must always be a positive number). First, let's create a `Python` file within the `test.semivariance` module. All files with unit tests should have a prefix `test_`. That's why we name our file `test_calculate_semivariance` and a full path to the file should be:

`pyinterpolate/test/semivariance/test_calculate_semivariance.py`

In the beginning, we must import `calculate_semivariance()` method and `unittest` module.

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance
```

To write a test we must create a `Python` class that starts with a `Test_` prefix. Usually, it is named `TestYourFunctionOrModuleName`, in our case: `TestCalculateSemivariance`. This class inherits from the `unittest.TestCase`. We can skip the explanation of what inheritance is. The key is to understand that we can use methods from `unittest.TestCase` in our class `TestCalculateSemivariance` and those methods allow us to write unit tests. Let's update our script with this new piece of information:

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):

	def test_calculate_semivariance(self):
		pass
```

The good practice with unit testing is to have data that is not dependent on external sources or processes. In other words, we use mostly static files with known datasets or artificial arrays. Those arrays may be filled with random numbers of specific distribution or hard-coded values which are simulating possible input. We are going to create one array for the sake of simplicity:

```python
INPUT = [8, 6, 4, 3, 6, 5, 7, 2, 8, 9, 5, 6, 3]
```

This array is not random. It comes from the **Basic Linear Geostatistics** and it is presented on page 48 of this book. It has an important property: it is so simple that **we can calculate semivariance by hand**... and it will be a topic of functional testing scenario. Now we consider the output test to check if all outputs are positive numbers. Function `calculate_semivariance()` returns list of triplets: `[lag, semivariance, number of point pairs]` and we must check them all. First, calculate semivariances up to lag 5 by hand.

**Lag 0:**

Calculations: _n/a_

Expected output: `[0, 0, 13]`

**Lag 1:**

Calculations: 

$$\gamma(h_{1})= \frac{1}{2*24}*2(4+4+1+9+1+4+25+36+1+16+1+9)=\frac{111}{24}=4.625$$

Expected output: `[1, 4.625, 24]`

**Lag 2:**

Calculations:

$$\gamma(h_{2})= \frac{1}{2*22}*2(16+9+4+4+1+9+1+49+9+9+4)=\frac{115}{22}=5.227$$

Expected output: `[2, 5.227, 22]`

**Lag 3:**

Calculations:

$$\gamma(h_{3})= \frac{1}{2*20}*2(25+0+1+16+16+9+4+9+4+36)=\frac{120}{20}=6.0$$

Expected output: `[3, 6.0, 20]`

**Lag 4:**

Calculations:

$$\gamma(h_{4})= \frac{1}{2*18}*2(4+1+9+1+4+16+4+16+25)=\frac{80}{18}=4.444$$

Expected output: `[4, 4.444, 18]`

**Lag 5:**

$$\gamma(h_{5})= \frac{1}{2*16}*2(9+1+4+25+9+0+1+1)=\frac{50}{16}=3.125$$

Expected output: `[5, 5.0, 16]`

We prepare a static (constant) variable with the expected output values. At this step, we import `numpy` to create an output of the same type as the output returned by the `calculate_semivariance()` function.

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):

	def test_calculate_semivariance(self):
		INPUT = [8, 6, 4, 3, 6, 5, 7, 2, 8, 9, 5, 6, 3]
		EXPECTED_OUTPUT = np.array([
			[0, 0, 13],
			[1, 4.625, 24],
			[2, 5.227, 22],
			[3, 6.0, 20],
			[4, 4.444, 18],
			[5, 3.125, 16]
		])
```

The first test is very simple: we should check if the output values from the function are positive. The thing is easy because each column from the output array MUST BE a positive number but we are going to focus only on the middle column with the experimental semivariance values. So let's begin testing! Tests in `Python` are written with the `assert` keyword. We use specific conditions with this keyword, for example, we can:

* check if values are equal,
* check if values are not equal,
* check if given value is True or False,
* check if function raises specific error,
* and more... described [here](https://docs.python.org/3/library/unittest.html#test-cases).

For this tutorial, we will perform a simple test to check if all values of the output are greater than zero. `Unittest` package has method `assertGreaterEqual(a, b)` which can be used to test this case but we rather will use `assertTrue()` to test if the every compared record in array is greater than 0, or it is `True`. How are we going to test it? `NumPy` has a great system of checking different conditions on arrays. We can test any condition with the syntax:

```python
boolean_test = my_array[my_array >= 0]
```

Variable `boolean_test` is an array of the same size as `my_array` but with `True` and `False` values at the specific positions where the condition has been met or hasn't. We can add a method `.all()` to check if every value in the boolean array is `True`:

```python
boolean_test = my_array[my_array >= 0]
boolean_output = boolean_test.all()
```

Then we can assert value to our test. In this case, let's start with the artificial array `EXPECTED_OUTPUT`. Even if we know that everything is ok we should test the output to be sure that the test is designed properly:

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):

	def test_calculate_semivariance(self):
		INPUT = [8, 6, 4, 3, 6, 5, 7, 2, 8, 9, 5, 6, 3]
		EXPECTED_OUTPUT = np.array([
			[0, 0, 13],
			[1, 4.625, 24],
			[2, 5.227, 22],
			[3, 6.0, 20],
			[4, 4.444, 18],
			[5, 3.125, 16]
		])
		expected_output_values = EXPECTED_OUTPUT[:, 1].copy()
		
		boolean_test = (expected_output_values >= 0).all()
		
		assertTrue(boolean_test, 'Test failed. Calculated values are below zero which is non-physical. Check your data.')
```

The test should definitively pass but you may want to check what'll happen if you change one value to the negative number:

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):

	def test_calculate_semivariance(self):
		INPUT = [8, 6, 4, 3, 6, 5, 7, 2, 8, 9, 5, 6, 3]
		EXPECTED_OUTPUT = np.array([
			[0, 0, 13],
			[1, 4.625, 24],
			[2, 5.227, 22],
			[3, -6.0, 20],
			[4, 4.444, 18],
			[5, 3.125, 16]
		])
		expected_output_values = EXPECTED_OUTPUT[:, 1].copy()
		
		boolean_test = (expected_output_values >= 0).all()
		
		self.assertTrue(boolean_test, 'Test failed. Calculated values are below zero which is non-physical. Check your data.')
```

Surprised!? We got the assertion error and it informs us that we haven't correctly implemented something. Let's remove the minus sign from the expected output lag's number 3 semivariance and test the baseline function, but leave the `EXPECTED_OUTPUT` variable without a test.

Function `calculate_semivariance(data, step_size, max_range)` has three parameters. From the [documentation](https://pyinterpolate.readthedocs.io/en/latest/code_documentation/semivariance.html#calculate_semivariance()) we know that:

- `data` parameter is `numpy array` of coordinates and their values,
- `step_size` parameter is a distance between lags and has a `float` type,
- `max_range` is a maximum range of analysis and has a `float` type.

From this description we see that our `INPUT` variable is not what's expected by the algorithm. We will change it and we set `step_size` and `max_range` parameters within our test case and we remove `EXPECTED_OUTPUT` variable because **for this concrete test** we do not need it. We change a name of the unit test to the `test_positive_output(self)` to make our test goal clear.

```python
import unittest
from pyinterpolate.semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):

	def test_positive_output(self):
		INPUT = np.array([
			[0, 0, 8],
			[1, 0, 6],
			[2, 0, 4], 
			[3, 0, 3], 
			[4, 0, 6], 
			[5, 0, 5], 
			[6, 0, 7], 
			[7, 0, 2], 
			[8, 0, 8], 
			[9, 0, 9], 
			[10, 0, 5], 
			[11, 0, 6], 
			[12, 0, 3]
		])
		
		# Calculate experimental semivariance
		t_step_size = 1.1
		t_max_range = 6
		
		experimental_semivariance = calculate_semivariance(INPUT, t_step_size, t_max_range)
		
		boolean_test = (experimental_semivariance >= 0).all()
		
		self.assertTrue(boolean_test, 'Test failed. Calculated values are below zero which is non-physical.')
```

Terrific! We've performed our first unit test! (By the way: it is implemented in the package).

## Logical and functional testing

The next type of testing is a functionality test. Do you remember how we've calculated semivariance manually in the previous part? We are goin to use the `EXPECTED_OUTPUT` array from the last test:

```python
EXPECTED_OUTPUT = np.array([
			[0, 0, 13],
			[1, 4.625, 24],
			[2, 5.227, 22],
			[3, 6.0, 20],
			[4, 4.444, 18],
			[5, 3.125, 16]
		])
```

Why do we need it? A **unit test is important for software development** but **scientific software development is even more rigorous! We need to prove that our function works as expected** for a wide range of scenarios. It is not an easy task. Some algorithms haven't be proven mathematically... But we have an advantage: geostatistics is a well-established discipline and we can use external sources and pen and paper for testing.

We'd created a functionality test in the previous example but we didn't use it. The idea is to create a simple example now, that can test if our program calculates semivariance correctly. We've taken an example array from the literature and we calculated semivariance values for few lags by hand. We assume that the results are correct (if you are not sure check records with the book; but I couldn't assure you that everyone is wrong but your program). 

The `EXPECTED_OUTPUT` array is a static entity against which we can compare the results of our algorithm. We can distinguish few steps which need to be done for the logic test:

1. Get the reference input and create the reference output. This could be calculated manually, derived from the literature or estimated by the other program which produces the right output.
2. Store reference output as the script values or as the non-changing file. The key is to never change values in this output throughout the testing flow.
3. Estimate the expected output from the reference input with your algorithm.
4. Compare the expected output and the reference output. If both are the same (or very close) then test is passed.

A few steps of the functionality testing can be very misleading if we consider the complexity of the whole process. There are traps in here, and sometimes it is hard to avoid them:

- the reference input covers very unusual scenario, as example it is created from the dataset with the uniform distribution but in reality data follows other types of distribution,
- the reference input and output do not exist or the reference data is not publicly available (frequent problem with publications),
- it's hard to calculate the reference output manually.

But the biggest, and the most dangerous problem which sometimes relates to the previous obstacles is the reference output assignation from the developed algorithm itself (we don't know the reference output and we create it from our algorithm itself...). Is it bad? Yes, it is. Is it sometimes the only (fast and dirty) way of development? Yes, it is. But even in this scenario we usually perform some kind of sanity check against our belief how the output should look like. Let's put aside the dirty testing and come back to our case. We have the expected output and we have the reference input. With both, we can create the functionality test within our package!

The code will be similar to the unit test, e.g.: there's the same building logic behind it. We define a function within the class that inherits from `unittest.TestCase`. Within this function, we define the reference input and the reference output. We run the algorithm on the reference input and compare results to the reference output. `Numpy` has special assertion methods and we are going to pick `asset_almost_equal()` for our test scenario. When we compare two arrays with floating-point precision it is better to use this assertion which allows a small difference between the records in the compared arrays. The code itself is presented below:

```python
import unittest
from numpy.testing import assert_almost_equal

from pyinterpolate.semivariance.semivariogram_estimation.calculate_semivariance import calculate_semivariance

class TestCalculateSemivariance(unittest.TestCase):
    """
    HERE COULD BE A THE PART FOR OTHER TESTS
    """
    def test_against_expected_value_1(self):

        REFERENCE_INPUT = np.array([
            [0, 0, 8],
            [1, 0, 6],
            [2, 0, 4],
            [3, 0, 3],
            [4, 0, 6],
            [5, 0, 5],
            [6, 0, 7],
            [7, 0, 2],
            [8, 0, 8],
            [9, 0, 9],
            [10, 0, 5],
            [11, 0, 6],
            [12, 0, 3]
        ])

        EXPECTED_OUTPUT = np.array([
			[0, 0, 13],
			[1, 4.625, 24],
			[2, 5.227, 22],
			[3, 6.0, 20],
			[4, 4.444, 18],
			[5, 3.125, 16]
		])

        # Calculate experimental semivariance
        t_step_size = 1
        t_max_range = 6

        experimental_semivariance = calculate_semivariance(REFERENCE_INPUT, t_step_size, t_max_range)

        # Get first five lags
        estimated_output = experimental_semivariance[:6, :]

        # Compare
        err_msg = 'The reference output and the estimated output are too dissimilar, check your algorithm'
        assert_almost_equal(estimated_output, EXPECTED_OUTPUT, 3, err_msg=err_msg)
```

This code is coming from the official codebase, you can look there too, to see how it is slightly rearranged within the test script.

## How to run unit tests

There are two main ways to perform a unit test: within the IDE or manually, from the terminal. We prefer to use `PyCharm` for the tests and the next example will show how those can be performed within this IDE.

### PyCharm

#### PyCharm Single test

The unit testing in `PyCharm` is very easy. You should open **Project** view on the left side. You will see a window with the project tree. Navigate to the specific file with unit tests, click the right mouse button, select `Run Unittests in test_case_`. A test is under its way. After a while, you will see the test results in the bottom window.

![Single Test Case](imgs-testing/single_test_case.jpg)

#### PyCharm Multiple tests

To run multiple tests at once do as in the previous step but right-click on the directory with tests and select `Run Unittests in tests`. All tests for `Pyinterpolate` take about a minute to go.

![Multiple Tests at Once](imgs-testing/mutiple_test_cases.jpg)

### Console

A test from the console is slightly more complex than from the IDE. It is required to be within a test directory and to have installed all required packages - in the normal scenario those are installed within the virtual environment or the conda environment. Next:

1. Activate environment with the dependencies of a package.
2. Navigate to the directory with tests.

The next steps differ in the case if we'd like to test one script or every case in the package.

#### Console Single test

* Navigate to the directory with a specific test.
* Run `unittest` script:

```
python -m unittest test_calculate_semivariance
```

#### Console Multiple tests

* Navigate to the parent directory of tests which you want to cover.
* Run `unittest` script:

```
python -m unittest
```

## Use a tutorial for testing

`Pyinterpolate` follows a specific flow of development. In the beginning, each feature is developed within the Jupyter Notebook. We build a function or class with the concrete example in mind. Then it is easy to transform a development notebook into a tutorial. About half of the tutorials were developed in this manner. The next big thing with the tutorials is that we can perform the so-called **sanity check** within them. There are bugs, especially logical, which are not covered by the unit tests or functionality testing. Tutorials cover a large part of the workflow and it makes them a great resource for program testing. Thus after the implementation of new features or changes inside old features, we always **run all affected tutorials**. Probably you saw that each tutorial has a changelog table and there are records related to the feature changes. 

The development and testing within notebooks can be summarized in a few steps:

1. Develop a feature in the notebook. (Development)
2. Write a tutorial related to this feature. (Development)
3. After any change in the feature run tutorial to check if everything works fine. (Maintenance and Testing)
4. After development of the new features test affected tutorials if there are any. (Maintenance and Testing)

We strongly encourage any developer to go within this pattern. It builds confidence about the proposed development and the understanding of the covered topics, as semivariogram modeling, kriging, and so on. **Running all tutorials before a major release of `Pyinterpolate` is mandatory, even if the changes shouldn't affect a specific tutorial workflow**. It is another kind of test and it assures that users get a very good product.

## Summary

You have learned how to:

- write the unit test,
- write and prepare the functionality test,
- run tests from the IDE and console,
- use tutorials as the playground for tests.

If you are still not sure about anything in this article then write an [Issue](https://github.com/szymon-datalions/pyinterpolate/issues) within the package repository. We will answer all your questions!PyInterpolate
=============

PyInterpolate is designed as the Python library for geostatistics. It's role is to provide access to spatial statistics tools used in a wide range of studies.

Changes by date
===============

2021-12-XX
----------

**version 0.2.5**

* neighbors selection (lags counting) has been changed,
* `TheoreticalSemivariogram` searches for optimal sill in a grid search algorithm,
* corrected error in `Krige` class; now calculation of error variance is correct.

2021-12-11
----------

**version 0.2.4**

* `self.points_values` chenged to `self.points_array` in `TheoreticalSemivariogram` class,
* `NaN` values are tested and checked in `calc_semivariance_from_pt_cloud()` function,
* new semivariogram models included in the package: **cubic**, **circular**, **power**,
* corrected calculation of the closest neighbors for kriging interpolation,
* changed `prepare_kriging_data()` function,
* the new optional parameter `check_coordinates` (**bool**) of `calc_point_to_point_distance()` function to control the coordinates uniqueness tests. This test is very resource-consuming and should be avoided in a normal work and it should be performed before data injection into the modeling pipeline.
* the new `dev/profiling/` directory to test and profile parts of a code.

2021-08-23
----------

**version 0.2.3.post1**

* the outliers removal function: you can choose side for outlier detection and remove. Default is top, available are: both, top, down,
* the outliers removal function: changed algorithm,
* new tutorial about outliers and their influence on the final model.

2021-05-13
----------

**version 0.2.3**

* more parameters to store (and access) in TheoreticalSemivariogram class,
* error weighting against the linear regression model (ax + b),
* global mean for Simple Kriging as a required parameter,
* tqdm progress bar to `RegularizedSemivariogram.transform()` and `interpolate_raster()` functions,
* refactored Semivariogram Regularization: ranges are controlled by algorithm, not an user,
* added pull request template,
* added issues templates,
* bug in spherical semivariogram model,
* experimental variogram as points (not a solid line),
* inverse distance weighting function: algorithm, tests, documentation and new tutorial,
* changed output names of regularized data (`ArealKriging.regularize_data`) from **estimated value** to **reg.est** and from **estimated prediction error** to **reg.err**,
* error related to the id column as a string removed,
* TheoreticalSemivariogram `params` attribute changed to `nugget`, `sill` and `range` attributes.

2021-03-10
----------

**version 0.2.2.post2**

* directional semivariograms methods, docs and tests added,
* check if points are within elliptical area around point of interest method, docs and tests added,
* broken dependency in `README.md` corrected.

2021-03-02
----------

**version 0.2.2.post1**

* variogram point cloud methods, tutorials, docs and tests added,
* updated tutorials and baseline datasets to show examples with spatial correlation,
* updated `README.md`: contribution, example, sample image,
* data is tested against duplicates (points with the same coordinates),
* removed bug in `interpolate_raster()` method.
*****************************
Contribution to Pyinterpolate
*****************************

We love your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

* Reporting a bug
* Discussing the current state of the code
* Submitting a fix
* Proposing new features
* Becoming a maintainer

=====================
Where should I start?
=====================

On the Github profile of pypinterpolate! We use github to host code, to track issues and feature requests, as well as accept pull requests. We have Discord server too and it's available here:

> https://discord.gg/3EMuRkj

We Use `Github Flow <https://guides.github.com/introduction/flow/index.html>`_, So All Code Changes Happen Through Pull Requests
Pull requests are the best way to propose changes to the codebase (we use `Github Flow <https://guides.github.com/introduction/flow/index.html>`_). We actively welcome your pull requests:

1. Fork the repo and create your branch from `dev` or from 'master' (preferably 'dev').
2. If you've added code that should be tested, add tests in the test package.
3. If you've changed APIs, update the documentation.
4. Ensure the test suite passes.
5. Make sure your code lints.
6. Issue that pull request!

====================================================================================
Any contributions you make will be under the BSD 3-Clause "New" or "Revised" License
====================================================================================

In short, when you submit code changes, your submissions are understood to be under the same [BSD 3-Clause "New" or "Revised" License] that covers the project. Feel free to contact the maintainers if that's a concern.

=============
Bug reporting
=============

Report bugs using Github's `issues <https://github.com/szymon-datalions/pyinterpolate/issues>`_.
We use GitHub issues to track public bugs. Report a bug by opening a new issue.

Write bug reports with detail, background, and sample code
`This is an example <https://github.com/szymon-datalions/pyinterpolate/issues/4>`_.

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
- Be specific!
- Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports. I'm not even kidding.

=====================
Use a PEP8 Guidelines
=====================

`PEP8 Python Guidelines <https://www.python.org/dev/peps/pep-0008/>`_

=======
License
=======

By contributing, you agree that your contributions will be licensed under its BSD 3-Clause "New" or "Revised" License.

==========
References
==========

This document was adapted from the open-source contribution guidelines for `Facebook's Draft <https://github.com/facebook/draft-js/blob/a9316a723f9e918afde44dea68b5f9f39b7d9b00/CONTRIBUTING.md>`_
********************
Algorithms Explained
********************

Algorithms within Pyinterpolate package can be complex. The package performs semivarioragm moideling automatically and this process may lead to the strange results. In this section algorithms from within package are described to prepare the advanced users for the in-depth analysis.

Algorithms:
-----------

.. toctree::
   :maxdepth: 1

   algorithms_documentation/Automatic Fitting of the Semivariogram Model
   algorithms_documentation/Outliers Removal
.. pyinterpolate documentation master file, created by
   sphinx-quickstart on Tue Oct 13 13:48:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

*************
Pyinterpolate
*************

.. toctree::
   :maxdepth: 1

   introduction
   setup
   tutorials.rst
   code_documentation.rst
   algorithms.rst
   contribution
   license

===================================
What can you do with Pyinterpolate?
===================================

PyInterpolate is designed as the Python library for geostatistics. It's role is to provide access to spatial statistics tools used in a wide range of studies. This package helps you **interpolate spatial data** with *Kriging* technique. In the close future you'll use more spatial interpolation tools.

=====================
Where should I start?
=====================

If you are first time here, we highly recommend to read introduction and learn about **Pyinterpolate** purpose. Otherwise it's a good idea to look into tutorials and quickly adapt code from there to your own projects.*************
Documentation
*************

Welcome into the code documentation of the Pyinterpolate package. If you are looking for every method in the package this is the best place to start.

Package modules:
----------------

.. toctree::
   :maxdepth: 2

   code_documentation/idw
   code_documentation/io_ops
   code_documentation/transform
   code_documentation/distance
   code_documentation/semivariance
   code_documentation/kriging
   code_documentation/viz
   code_documentation/misc
*********
Tutorials
*********

The easiest way to learn new package is by coding. That's why you may use tutorials to start your coding without analysis of the package and it's code.

Available tutorials:
--------------------

.. toctree::
   :maxdepth: 1

   tutorials/Variogram Point Cloud (Basic)
   tutorials/Semivariogram Estimation (Basic)
   tutorials/Ordinary and Simple Kriging (Basic)
   tutorials/How good is our Kriging model - test  it with IDW algorithm (Basic)
   tutorials/Blocks to points Ordinary Kriging interpolation (Intermediate)
   tutorials/Semivariogram Regularization (Intermediate)
   tutorials/Outliers and Their Influence on the Final Model (Intermediate)
   tutorials/Poisson Kriging - Centroid Based (Advanced)
   tutorials/Poisson Kriging - Area to Area (Advanced)
   tutorials/Poisson Kriging - Area to Point (Advanced)

.. pyinterpolate documentation master file, created by
   sphinx-quickstart on Tue Oct 13 13:48:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyinterpolate
=============

PyInterpolate is designed as the Python library for geostatistics. It's role is to provide access to spatial statistics tools used in a wide range of studies. If you're:

- GIS expert,
- geologist,
- mining engineer,
- ecologist,
- public health specialist,
- data scientist.

Then this package may be useful for you. You may use it for:

- spatial interpolation and spatial prediction,
- alone or with machine learning libraries,
- for point and areal datasets.

Package was tested in commercial and research projects:

- Tick-borne Disease Detector (prediction of areas of infection risk), research project funded by European Space Agency,
- commercial project related to the prediction of demand for specific flu medications,
- commercial project related to the large-scale infrastructure maintenance.

Pyinterpolate allows you to perform:

- Ordinary Kriging and Simple Kriging (spatial interpolation from points),
- Centroid-based Kriging of Polygons (spatial interpolation from blocks and areas),
- Area-to-area and Area-to-point Poisson Kriging of Polygons (spatial interpolation and data deconvolution from areas to points).

In the future new interpolation techniques will be introduced.


Contents:
---------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Calculate Distances
===================

.. automodule:: pyinterpolate.calculations.distances.calculate_distances
   :members:

Data Processing
===============

.. automodule:: pyinterpolate.data_processing.data_preparation.get_points_within_area
   :members:

.. automodule:: pyinterpolate.data_processing.data_preparation.prepare_areal_shapefile
   :members:

.. automodule:: pyinterpolate.data_processing.data_preparation.read_data
   :members:

Data Transformation
===================

.. automodule:: pyinterpolate.data_processing.data_transformation.get_areal_centroids
   :members:

.. automodule:: pyinterpolate.data_processing.data_transformation.prepare_kriging_data
   :members:

Data Visualization
==================

.. automodule:: pyinterpolate.data_visualization.interpolate_raster
   :members:

Kriging - Areal Kriging
=======================

.. automodule:: pyinterpolate.kriging.areal_poisson_kriging.areal_kriging
   :members:

.. automodule:: pyinterpolate.kriging.areal_poisson_kriging.centroid_based.centroid_poisson_kriging
   :members:

Kriging - Point Kriging
=======================

.. automodule:: pyinterpolate.kriging.point_kriging.kriging
   :members:

Semivariance - Semivariogram Deconvolution
==========================================

.. automodule:: pyinterpolate.semivariance.semivariogram_deconvolution.regularize_semivariogram
   :members:

Semivariance - Semivariogram Estimation
=======================================

.. automodule:: pyinterpolate.semivariance.semivariogram_estimation.calculate_semivariance
   :members:

Semivariance - Fit Semivariogram
================================

.. automodule:: pyinterpolate.semivariance.semivariogram_fit.fit_semivariance
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
