# Contributing to Easy SimAuto
We welcome contributions to ESA! If you find a bug, please
file an issue on [Github](https://github.com/mzy2240/ESA/issues).

The primary purpose of this document is to describe what's required to 
make a contribution to the source code.

## Fork and Pull Request
The simplest method to contribute is to first fork the repository, make
changes, and then submit a pull request.

## Starting Out
While this document attempts to be comprehensive, the best way to get 
a feel for the project is to peruse the source code. Specifically, take
a look at esa/saw.py and tests/test_saw.py.

## Conventions and Style
The following sections describe expected conventions and style for 
contributing to ESA.
### PEP-8
In general, ESA follows [PEP-8](https://www.python.org/dev/peps/pep-0008/).
Please read the PEP in full if you have not already. The good news is
you don't need to memorize everything - modern IDEs like PyCharm make 
following PEP-8 very easy.

The only time one should deviate from the PEP-8 is with regard to
function naming and function input variable naming. Functions and their
inputs should be named to match [PowerWorld's documentation](https://www.powerworld.com/WebHelp/#MainDocumentation_HTML/Simulator_Automation_Server.htm%3FTocPath%3DAutomation%2520Server%2520Add-On%2520(SimAuto)%7C_____1).

Example:
```python
def ChangeParametersMultipleElement(self, ObjectType: str, ParamList: list,
                                    ValueList: list):
...
``` 

Notice the function and input variables exactly match PowerWorld. 
However, internal variables should conform to PEP-8. E.g., following the
previous example above you may do the following:
```python
# Cast ObjectType to lower case so it matches dictionary keys. 
object_type = ObjectType.lower()
```

ESA follows the convention that attributes/methods/etc. which start 
with an underscore are private.

### Docstrings and Type Hinting
Every function should include a detailed docstring that describes what 
the function does, and also describes all parameters and return values
in detail. Additionally, the docstring should provide a direct link to
applicable PowerWorld documentation. It's also generally useful to 
document what exceptions the function/method may raise.

Docstrings should use reStructuredText (rst) format, as described in
[PEP-287](https://www.python.org/dev/peps/pep-0287/). A useful cheat
sheet can be found [here](https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html). 

Additionally, functions should utilize type hinting to help users and 
developers know explicitly what types of objects they'll be dealing 
with.

Here's an example of both type hinting and a good docstring:
```python
def OpenCase(self, FileName: Union[str, None] = None) -> None:
    """Load PowerWorld case into the automation server.

    :param FileName: Full path to the case file to be loaded. If
        None, this method will attempt to use the last FileName
        used to open a case.

    :raises TypeError: if FileName is None, and OpenCase has never
        been called before.

    `PowerWorld documentation
    <https://www.powerworld.com/WebHelp/Content/MainDocumentation_HTML/OpenCase_Function.htm>`_
    """
```
 
### Inputs and Outputs
Where applicable, the preferred return types are `pandas.DataFrame` and
`pandas.Series`. A return of `None` should be used to indicate that a 
function operated but has nothing to return.

## Unit Testing
Any and all functions, methods, classes, etc. should be tested. ESA 
uses the built-in [unittest](https://docs.python.org/3/library/unittest.html)
module, as well as [unittest.mock](https://docs.python.org/3/library/unittest.mock.html).
 
The objective is to have 100% testing coverage. Tools such as
[Coverage.py](https://coverage.readthedocs.io/en/latest/) can be used
to assess test coverage.
 
Please read the docstring in test_saw.py - it has very important 
information related to avoiding state conflicts that may arise 
between tests (in an ideal world, this wouldn't be a problem, but we 
just don't live in an ideal world :)).

Note that due to the nature of this project, many tests aren't truly
"unit" tests in that they actually call SimAuto. 

## Using the Helper Methods
The SAW class has a variety of private helper methods, prefixed with
and underscore. A developer should use these liberally, as they're 
designed to make development easy. 

### Details for the SAW class
Never call SimAuto directly - always use the `_call_simauto` helper.

Any and all DataFrames/Series that are created via output from
PowerWorld should be passed through the `clean_df_or_series` method.# Data
Use this directory to store data associated with tests.

## CandidateLines.csv
This file was provided by [Adam Birchfield](http://adambirchfield.com/)
on 2019-11-21 and is intended to be used to add lines to the Texas 
2000 bus case. ---
title: 'Easy SimAuto (ESA): A Python Package that Simplifies Interacting with PowerWorld Simulator'
tags:
  - Python
  - PowerWorld
  - Simulator
  - Automation
  - Server
  - SimAuto
  - ESA
  - Power Systems
  - Electric Grid
  - Smart Grid
  - Energy
  - Numpy
  - Pandas
authors:
  - name: Brandon L. Thayer^[The first two authors of this paper contributed equally to this software.]
    orcid: 0000-0002-6517-1295
    affiliation: "1, 2"
  - name: Zeyu Mao
    orcid: 0000-0003-0841-5123
    affiliation: 1
  - name: Yijing Liu
    orcid: 0000-0002-5104-325X
    affiliation: 1
  - name: Katherine Davis
    orcid: 0000-0002-1603-1122
    affiliation: 1
  - name: Thomas Overbye
    orcid: 0000-0002-2382-2811
    affiliation: 1
affiliations:
  - name: "Texas A&M University"
    index: 1
  - name: "Pacific Northwest National Laboratory"
    index: 2
date: 14 April 2020
bibliography: paper.bib
---

# Summary

The electric power system is an essential cornerstone of modern society,
enabling everything from the internet to refrigeration. Due to a variety
of forces including climate change, changing economics, and the digital
computer revolution, the electric grid is undergoing a period of major
change. In order to overcome current and upcoming challenges in the
electric power system, such as integrating renewable resources into a
system that was not designed for intermittent power sources,
researchers and industry practitioners must simulate the electric grid,
its component devices, and its operation.

[PowerWorld Simulator](https://www.powerworld.com/) is a commercial
power systems simulation tool that contains a suite of modeling and
simulation features including power flow simulation, contingency
analysis, transient stability simulation, and more [@powerworld]. The
Simulator Automation Server (SimAuto) add-on for PowerWorld provides an
application programming interface (API) that operates in-memory,
allowing users to rapidly configure, run, and obtain results for
simulations. PowerWorld and SimAuto are commonly used throughout the
research community as well as in industry.

SimAuto was designed to be flexible enough to work with most available
programming languages. However, the combination of this flexibility and
the in-memory nature of SimAuto communication can make using SimAuto
challenging, requiring error-checking, data type conversions, data
parsing, low-level interactions with Windows Component Object Model
(COM) objects, and more.

[Easy SimAuto (ESA)](https://github.com/mzy2240/ESA) is a Python package
that significantly simplifies interfacing with PowerWorld Simulator
[@esa]. ESA wraps all available SimAuto functions; provides high-level
helper functions to streamline workflows, and provide additional
functionality not provided by SimAuto; performs automatic error
checking, data type conversions, and data parsing; is easily installable
via Python's package installer (pip); has 100% testing coverage; and is
fully documented. Similar packages have been created in the past, but
lack functions, tests, documentation, and other useful features ESA
provides [@pypowerworld], [@matpws]. Most SimAuto users tend to write
their own one-off functions and boilerplate code for interfacing with
SimAuto. ESA eliminates this redundancy and abstracts away all the
low-level SimAuto interactions so that users can focus on performing
higher-level tasks such as automating tasks, configuring simulations,
and analyzing results.

ESA helps to meet the needs of both power system researchers and 
practitioners. As the design and operation of the electric grid becomes
more complex, researchers and developers need the ability to incorporate
their programs, algorithms, control schemes, etc. into power system
simulations. ESA enables its users to fully leverage, extend, and
automate the large depth of functionality and tools built into
PowerWorld Simulator: procedures that may have previously been
performed via a sequence of manual tasks in Simulator's graphical user
interface (GUI) can be rapidly built into Python scripts which can be
stored in version control and run with a single click. Since ESA uses
data types common to data science and scientific computing (e.g., Pandas
DataFrames and Numpy Arrays), it is well suited to both academic
research and task automation in industry. Due to ESA's use of these
common Python data types and libraries, ESA provides a much-needed
bridge between power system simulation and machine learning libraries.

ESA has already been utilized in several research projects past and
present:

- In [@gym-powerworld], [@brandon_arxiv], ESA was used to create a
standardized reinforcement learning environment for power system voltage
control. This environment was then used to carry out deep reinforcement
learning (DRL) experiments in which the algorithm attempts to learn how
to best control grid voltages under a diverse set of grid conditions 
[@drl-powerworld]. 
- In [@scenario_development], ESA was leveraged to create and simulate 
different electric grid scenarios where load, renewable generation 
levels, generation capacities, scheduled outages, and unit commitment
were all varied. The resulting scenarios were used in the
[Grid Optimization (GO) competition](https://gocompetition.energy.gov/)
hosted by the U.S. Department of Energy (DOE).
- Geomagnetic disturbances (GMDs) affect the magnetic and electric field
of the earth, inducing dc voltage sources superimposed on transmission
lines. In [@OverbyeKPEC]^[accepted, to be published after the delayed
conference takes place], a planning-based GMD mitigation strategy was
developed for large power systems. ESA is leveraged to programmatically
place GIC blocking devices in test systems per the proposed algorithm,
thus minimizing the effects of GMDs on the power grid.
- ESA is used by an ongoing research project entitled "Real Time
Monitoring Applications for the Power Grid under Geomagnetic
Disturbances (GMD)," where recently, a real-world GMD monitoring system
consisting of six magnetometers was deployed in Texas. The resulting
magnetic field measurements are coupled with ground conductivity models
to calculate real-time electric fields. These can then be fed to a grid
model of Texas using ESA to enable calculation of real-time
geomagnetically induced currents (GICs) for monitoring and
visualization.
- ESA is used by an ongoing research project entitled "Cyber Physical 
Resilient Energy Systems (CYPRES)". In this project, ESA is leveraged to
automatically map the communication system (like DNP3 outstation and 
data points) to the power system model.
- ESA is used by an ongoing research project entitled "Generalized 
Contingency Analysis Based on Graph Theory Concepts and Line Outage 
Distribution Factors (LODF)." In this project, ESA is leveraged to 
extract the topology of the power system model and obtain the LODF 
matrix.

# Acknowledgements

ESA was developed by researchers at Texas A&M University. Funding was
provided by the Texas A&M Engineering Experiment Station's Smart Grid
Center and the U.S. Department of Energy (DOE) under award DE-OE0000895.

The authors of ESA would like to also thank our fellow researchers at
Texas A&M who have provided essential feedback during ESA's creation.

# References
