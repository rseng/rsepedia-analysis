# Kinetics Toolkit

An Open-Source Python Package to Facilitate Research in Biomechanics

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03714/status.svg)](https://doi.org/10.21105/joss.03714)

Please consult https://kineticstoolkit.uqam.ca for information.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Python installation**
Select: python.org, Anaconda, conda-forge, etc.

**Operating system**
 - OS: [e.g. macOS, Windows, Linux]
 - Version [e.g. Big Sur, 10]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
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

**Are you able and willing to help?**
Please name your ability to help. Many things are possible:
- Providing example files
- Developing the feature
- Testing the feature
- etc.
---
title: 'Kinetics Toolkit: An Open-Source Python Package to Facilitate Research in Biomechanics'

tags:
  - Python
  - biomechanics
  - kinetics
  - kinematics
  - timeseries

authors:
  - name: Félix Chénier
    orcid: 0000-0002-2085-6629
    affiliation: "1, 2"

affiliations:
 - name: Department of Physical Activity Sciences, Université du Québec à Montréal (UQAM), Montreal, Canada
   index: 1
 - name: Mobility and Adaptive Sports Research Lab, Centre for Interdisciplinary Research in Rehabilitation of Greater Montreal (CRIR), Montreal, Canada
   index: 2

date: June 7th, 2011

bibliography: JOSS.bib

---

# Summary

Kinetics Toolkit is a Python package for generic biomechanical analysis of human motion that is easily accessible by new programmers. The only prerequisite for using this toolkit is having minimal to moderate skills in Python and Numpy.

While Kinetics Toolkit provides a dedicated class for containing and manipulating data (`TimeSeries`), it loosely follows a procedural programming paradigm where processes are grouped as interrelated functions in different submodules, which is consistent with how people are generally introduced to programming. Each function has a limited and well-defined scope, making Kinetics Toolkit generic and expandable. Particular care is given to documentation, with extensive tutorials and API references. Special attention is also given to interoperability with other software programs by using Pandas Dataframes (and therefore CSV files, Excel files, etc.), JSON files or C3D files as intermediate data containers.

Kinetics Toolkit is accessible at `https://kineticstoolkit.uqam.ca` and is distributed via conda and pip.


# Statement of need

The last decade has been marked by the development of several powerful open-source software programs in biomechanics. Examples include:
OpenSim [@seth_opensimsimulatingmusculoskeletal_2018],
SimBody [@sherman_simbodymultibodydynamics_2011],
Biordb [@michaud_biorbdpythonmatlab_2021],
BiomechZoo [@dixon_biomechzooopensourcetoolbox_2017],
Pinocchio [@carpentier_pinocchiolibraryfast_2019],
FreeBody [@cleather_developmentsegmentbasedmusculoskeletal_2015],
CusToM [@muller_custommatlabtoolbox_2019],
as well as many others. However, many of these tools are rather specific (e.g., musculoskeletal modelling, neuromuscular optimization, etc.) and not especially well suited for performing generic processing of human motion data such as filtering data, segmenting cycles, changing coordinate systems, etc. Other software programs, while being open source, rely on expensive closed-source software such as Matlab (Mathworks LCC, Naticks, USA).

While Matlab has a long and successful history in biomechanical analysis, it is quickly becoming challenged by the free and open-source Python scientific ecosystem, particularly by powerful packages, including Numpy [@harris_arrayprogrammingnumpy_2020], Matplotlib [@hunter_matplotlib2dgraphics_2007], SciPy [@virtanen_scipyfundamentalalgorithms_2020] and Pandas [@mckinney_pandasfoundationalpython_2011]. Since Python is regarded as a robust introductory programming language for algorithm development [@fangohr_comparisonmatlabpython_2004], it may be an ideal tool for new programmers in biomechanics.

The Pyomeca toolbox [@martinez_pyomecaopensourceframework_2020] is a Python library for biomechanical analysis. It uses an object-oriented programming paradigm where each data class (`Angles`, `Rototrans`, `Analogs`, `Markers`) subclasses xarray [@hoyer_xarrayndlabeled_2017], and where the data processing functions are accessible as class methods. While this paradigm may be compelling from a programmer's perspective, it requires users to master xarray and object-oriented concepts such as class inheritance, which are not as straightforward to learn, especially for new programmers who may just be starting out with Python and Numpy.

With this beginner audience in mind, Kinetics Toolkit is a Python package for generic biomechanical analysis of human motion. It is a user-friendly tool for people with little experience in programming, yet elegant, fun to use and still appealing to experienced programmers. Designed with a mainly procedural programming paradigm, its data processing functions can be used directly as examples so that users can build their own scripts, functions, and even modules, and therefore make Kinetics Toolkit fit their own specific needs.


# Features

## TimeSeries

Most biomechanical data is multidimensional and vary in time. To make it easier for researchers to manipulate such data, Kinetics Toolkit provides the `TimeSeries` data class. Largely inspired by Matlab's `timeseries` and `tscollection`, this data class contains the following attributes:

- `time`: Unidimensional numpy array that contains the time;
- `data`: Dict or numpy arrays, with the arrays' first dimension corresponding to time;
- `time_info` and `data_info`: Metadata corresponding to time and data (e.g., units);
- `events`: Optional list of events.

In addition to storing data, it also provides methods to:

- manage events (e.g., `add_event`, `rename_event`);
- manage metadata (e.g., `add_data_info`, `remove_data_info`);
- split data based on time indexes, times or events (e.g., `get_ts_after_time`, `get_ts_between_events`);
- extract or combine data (e.g., `get_subset`, `merge`);
- convert from and to other formats (e.g., `from_dataframe`, `to_dataframe`)
- etc.


## Processing data

All the data processing functions are included in submodules, for example:

- `filters` to apply frequency or time-domain filters to the TimeSeries data;
- `cycles` to detect and time-normalize cycles;
- `geometry` to express points, vectors, and frames in different global coordinate systems;
- `kinematics` to work with C3D files -- thanks to the `ezc3d` library [@michaud_ezc3deasyc3d_2021] -- and to perform higher-level manipulations on markers and rigid bodies;
- etc.


## Visualizing 3D kinematics

Kinetics Toolkit provides the `Player` class, which is a simple interactive 3D visualization tool for markers, bodies and segments. The user can pan and orbit, select and follow markers, animate at different speeds and navigate in time. Since Player is based on Matplotlib, it integrates well with various setups, using either the standard Python interpreter or IPython-based environments such as Spyder or Jupyter. Being integrated with the IPython event loop, multiple Player instances can be used at the same time, without blocking the interpreter.


## Saving and loading

Kinetics Toolkit provides `save` and `load` functions to store any standard Python-type data, Numpy arrays, Pandas Series and Dataframes, TimeSeries, or lists and dictionaries that contain such data types . These data are stored into a custom `ktk.zip` (which is an archive of standard JSON files) that is easily opened in other software programs such as Matlab.


# Acknowledgements

We want to acknowledge the dedicated people involved in major software programs and packages used by Kinetics Toolkit, such as Python, Numpy, Matplotlib, Pandas, Jupyter, Pytest, Sphinx, and many others. We also wish to thank Benjamin Michaud for creating and maintaining the `ezc3d` package, which is used by Kinetics Toolkit to read and write C3D files.


# References
Some review criteria (to help the reviewers)
--------------------------------------------

### Hosting

Kinetics Toolkit source is hosted on github at this address: https://github.com/felixchenier/kineticstoolkit

### Software license

Kinetics Toolkit is distributed under the Apache 2.0 license. See LICENSE.txt in the main repository.

### Substantial scholarly effort

Kinetics Toolkit has been developed first as a Matlab toolkit between 2015 and 2019, and was then converted and expanded as a Python toolkit from 2019 up to now. It has been open source since July 2019. Initially started as a software framework for my emerging research lab, it matured in a form that I now consider helpful for lots of other researchers.

Having a formation in Engineering, I have a professor position in a Physical Activity Science department and I think I have a good understanding of the knowledge and abilities of typical students and researchers in this domain. Kinetics Toolkit is oriented directly toward this clientele, and no other package I know of have this same objective.

Kinetics Toolkit has yet to be cited in literature because I wanted to publish a paper in JOSS before publicizing the toolkit, so that it can be referred using this paper. While at our lab, we personally used and tweaked Kinetics Toolkit over the years, I now consider this toolkit mature enough to be used with confidence by other researchers.

### Documentation

All documentation, including biomechanics basics, tutorials, API, statement of need and installation instructions, community guidelines (how to contribute) are available at https://kineticstoolkit.uqam.ca.

The development website, which refers to the most up to date master branch, also documents unstable features being developed, and is available at https://kineticstoolkit.uqam.ca/master.

### Tests

Kinetics Toolkit uses several way to ensure the code integrity. For each commit on the master branch, the following items are checked:

- Types are checked using mypy;
- API documentation is tested using doctest;
- Modules and class methods are tested using unit tests;
- All tutorials are rebuilt and checked for failure using nbconvert;
- The whole website is rebuilt and checked for failure using sphinx.

Unit tests are provided into the `test` folder, and data for tests and tutorials are provided in the `data` folder.
